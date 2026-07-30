# Title: "GO enrichment analysis of DTGs in RT data"
# Author: Aglaia Szukala
# Short description: Script to perform GO enrichment analysis of genes that are differentially targeted by sRNAs and possibly also differentially expressed

# Set project directory
# Set project root (portable; works after cloning the repository)
.root_helper <- c("scripts/00_set_project_root.R", "../scripts/00_set_project_root.R")
source(.root_helper[file.exists(.root_helper)][1])
rm(.root_helper)

# Output directory
out_dir <- "output/DT_RT/8_GOenrichment/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Needs script "4_GO_functions.r"
source("scripts/4_GO_functions.R")

# Load libraries
library(topGO)
library(GOplot)
library(xtable)
library(drake)
library(dplyr)
library(patchwork)
library(ggplot2)
pkgconfig::set_config("drake::strings_in_dots" = "literals")

###################
# Load input data #
###################

# GO terms annotation of the genome
geneID2GO <- readMappings(file = "data/GO/geneID2GO_Athaliana_genes.txt")
head(geneID2GO)

# Create GO to gene mapping file
GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)

# Read lists of DEGs 
files <- list.files(
  "data/DEGlists_RT_SzukalaEtAl2023/",
  pattern = "^DEG_.*\\.txt$",
  full.names = TRUE
)

deg_list <- lapply(files, read.table)

names(deg_list) <- sub("^DEG_|\\.txt$", "", basename(files))
names(deg_list) <- sub("\\.txt$", "", names(deg_list))
names(deg_list)

# Make a metadata to classify as within or between ecotype(s) comparison
comparisons <- sub("^DEG_", "", basename(files))
comparisons <- sub("\\.txt$", "", comparisons)

meta <- data.frame(
  comparison = comparisons,
  type = ifelse(grepl("vs", comparisons), "between_ecotype", "within_ecotype")
)

meta

# Load lists of DTRs
files <- list.files(
  "output/DT_RT/4_DT_edgeR/",
  pattern = "_DTRs_qlTest_(pos|neg)logFC\\.txt$",
  full.names = TRUE
)

# Upregulated (positive logFC)
up_files <- files[grepl("poslogFC\\.txt$", files)]

dtr_up_list <- lapply(up_files, read.table)

names(dtr_up_list) <- sub(
  "_DTRs_qlTest_poslogFC\\.txt$",
  "",
  basename(up_files)
)

# Downregulated (negative logFC)
dw_files <- files[grepl("neglogFC\\.txt$", files)]

dtr_dw_list <- lapply(dw_files, read.table)

names(dtr_dw_list) <- sub(
  "_DTRs_qlTest_neglogFC\\.txt$",
  "",
  basename(dw_files)
)

names(dtr_up_list)
names(dtr_dw_list)

# Check one
dtr_up_list["AvsA1"]

# Region-to-gene mapping
reg2gene <- read.table("data/gene2region.uniq.txt")
colnames(reg2gene) <- c("gene", "region")

####################################
# Collapse DTRs to gene-level DTGs #
####################################

# Create a mapping table to pair DEGs and DTRs lists in the correct way
deg_dtr_pairs <- data.frame(
  deg_id = c(
    "A1ud", "A3ud", "M1ud", "M3ud",
    "A1vsM1d", "A1vsM1u",
    "A3vsM3d", "A3vsM3u"
  ),
  dtr_id = c(
    "AvsA1", "AvsA3", "MvsM1", "MvsM3",
    "AvsM_1D", "AvsM_1U",
    "AvsM_3D", "AvsM_3U"
  )
)
deg_dtr_pairs

# DEG genes that are also linked to DTRs, per pair
# Keeps DEG statistics and adds DTR statistics
get_DEG_DTR_table <- function(deg_id, dtr_id, deg_list, dtr_up_list, dtr_dw_list, reg2gene) {
  
  deg_tab <- deg_list[[deg_id]]
  deg_tab$gene <- rownames(deg_tab)
  
  dtr_up <- dtr_up_list[[dtr_id]]
  dtr_up$region <- rownames(dtr_up)
  dtr_up$DTR_direction <- "DTR_up"
  
  dtr_dw <- dtr_dw_list[[dtr_id]]
  dtr_dw$region <- rownames(dtr_dw)
  dtr_dw$DTR_direction <- "DTR_down"
  
  dtr_all <- rbind(dtr_up, dtr_dw)
  
  dtr_gene <- merge(reg2gene, dtr_all, by = "region")
  
  out <- merge(deg_tab, dtr_gene, by = "gene")
  
  colnames(out) <- gsub("\\.x$", "_DEG", colnames(out))
  colnames(out) <- gsub("\\.y$", "_DTR", colnames(out))
  
  keep_cols <- c(
    "gene",
    "region",
    "logFC_DEG",
    "logCPM_DEG",
    "LR",
    "PValue_DEG",
    "padj",
    "logFC_DTR",
    "logCPM_DTR",
    "F",
    "PValue_DTR",
    "FDR",
    "DTR_direction"
  )
  
  out <- out[, keep_cols[keep_cols %in% colnames(out)]]
  
  out
}

# Run for all matched comparisons
deg_dtr_tables <- lapply(seq_len(nrow(deg_dtr_pairs)), function(i) {
  get_DEG_DTR_table(
    deg_id = deg_dtr_pairs$deg_id[i],
    dtr_id = deg_dtr_pairs$dtr_id[i],
    deg_list = deg_list,
    dtr_up_list = dtr_up_list,
    dtr_dw_list = dtr_dw_list,
    reg2gene = reg2gene
  )
})

names(deg_dtr_tables) <- deg_dtr_pairs$deg_id

## Inspect
head(deg_dtr_tables[["A1ud"]]) # no overlap!
head(deg_dtr_tables[["M1ud"]])
head(deg_dtr_tables[["A3ud"]]) # no overlap!
head(deg_dtr_tables[["M3ud"]]) # only two overlaps
head(deg_dtr_tables[["A1vsM1d"]])
head(deg_dtr_tables[["A1vsM1u"]])
head(deg_dtr_tables[["A3vsM3d"]])
head(deg_dtr_tables[["A3vsM3u"]])

## Save one file per comparison
out_dir <- "output/DT_RT/8_GOenrichment/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (nm in names(deg_dtr_tables)) {
  write.table(
    deg_dtr_tables[[nm]],
    file = file.path(out_dir, paste0("DEG_with_DTR_", nm, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

#############################################
# RUN GO TERM ENRICHMENT of genes DE and DT #
#############################################

# Output labels for filenames (match existing AvsM1d / AvsM3d naming)
go_labels <- c(
  A1ud = "A1ud",
  A3ud = "A3ud",
  M1ud = "M1ud",
  M3ud = "M3ud",
  A1vsM1d = "AvsM1d",
  A1vsM1u = "AvsM1u",
  A3vsM3d = "AvsM3d",
  A3vsM3u = "AvsM3u"
)

# Minimum number of unique DE+DT genes required to run enrichment
min_genes_go <- 10

run_go_enrichment <- function(resEdgeRtopTags, label, out_dir, geneNames, geneID2GO) {
  # One row per gene: keep DTR with lowest FDR
  resEdgeRtopTags <- resEdgeRtopTags %>%
    dplyr::arrange(gene, FDR) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  deGenes <- unique(resEdgeRtopTags$gene)
  if (length(deGenes) < min_genes_go) {
    message(
      "Skipping GO for ", label, ": only ", length(deGenes),
      " DE+DT genes (min = ", min_genes_go, ")"
    )
    return(invisible(NULL))
  }

  plan <- drake_plan(
    name_list = c("GO.ID", "Term", "Annotated", "Significant",
                  "Expected", "weight01_pval", "branch"),

    table = as.factor(geneNames) %in% deGenes,
    int_table = as.integer(table),
    int_fac_table = factor(int_table),
    fac_table = rename(table = int_fac_table, geneNames = geneNames),

    GOdata.BP = new("topGOdata", ontology = "BP",
                    allGenes = fac_table,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO),

    GOdata.MF = new("topGOdata", ontology = "MF",
                    allGenes = fac_table,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO),

    GOdata.CC = new("topGOdata", ontology = "CC",
                    allGenes = fac_table,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO),

    resultWeight01.BP = runTest(GOdata.BP, algorithm = "weight01", statistic = "fisher"),
    resultWeight01.MF = runTest(GOdata.MF, algorithm = "weight01", statistic = "fisher"),
    resultWeight01.CC = runTest(GOdata.CC, algorithm = "weight01", statistic = "fisher"),

    allRes.BP1 = GenTable(GOdata.BP, weight01_pval = resultWeight01.BP,
                          orderBy = "weight01", ranksOf = "weight01", topNodes = 100),
    allRes.BP2 = cbind(allRes.BP1, "BP"),
    allRes.BP = change_names(data = allRes.BP2, name_list = name_list),

    allRes.MF1 = GenTable(GOdata.MF, weight01_pval = resultWeight01.MF,
                          orderBy = "weight01", ranksOf = "weight01", topNodes = 100),
    allRes.MF2 = cbind(allRes.MF1, "MF"),
    allRes.MF = change_names(data = allRes.MF2, name_list = name_list),

    allRes.CC1 = GenTable(GOdata.CC, weight01_pval = resultWeight01.CC,
                          orderBy = "weight01", ranksOf = "weight01", topNodes = 100),
    allRes.CC2 = cbind(allRes.CC1, "CC"),
    allRes.CC = change_names(data = allRes.CC2, name_list = name_list),

    allRes1 = rbind(allRes.BP, allRes.MF),
    allRes = rbind(allRes1, allRes.CC),

    allGO.BP = genesInTerm(GOdata.BP),
    allGO.MF = genesInTerm(GOdata.MF),
    allGO.CC = genesInTerm(GOdata.CC),
    allGO = c(allGO.BP, allGO.MF, allGO.CC),

    SAM_ANOTATION = lapply(allGO, function(x) x[x %in% deGenes]),
    enriched_go_with_my_genes = lapply(SAM_ANOTATION[allRes[, 1]], paste0, collapse = ", "),
    enriched_go_with_my_genes.list = attach_enriched_go_genes(enriched_go_with_my_genes),

    go.dataframe = data.frame(
      Category = allRes$branch,
      ID = allRes$GO.ID,
      Term = allRes$Term,
      Genes = as.vector(enriched_go_with_my_genes.list),
      adj_pval = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE))
    ),

    ## IMPORTANT: use DEG values here
    EC.genelist = data.frame(
      ID = resEdgeRtopTags$gene,
      logFC = resEdgeRtopTags$logFC_DEG,
      adj.P.Val = resEdgeRtopTags$padj
    ),

    circ = circle_dat(go.dataframe, EC.genelist)
  )

  make(plan)
  circ <- readd(circ)

  write.table(
    circ,
    file = file.path(out_dir, paste0("Genes_DE_DT_", label, "_GO_circ.txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    col.names = TRUE
  )

  pdf(
    file = file.path(out_dir, paste0("Genes_DE_DT_", label, "_GO.pdf")),
    width = 12,
    height = 12
  )

  print(
    ggplot(
      circ[circ$adj_pval < 0.05 & circ$count > 1 &
             circ$category != "CC" & circ$category != "MF", ],
      aes(x = term, y = -log10(adj_pval))
    ) +
      geom_col(aes(fill = zscore)) +
      scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0,
        name = "Z-score"
      ) +
      geom_text(aes(label = count), hjust = -1, size = 5) +
      coord_flip() +
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size = 18),
        legend.position = "top"
      )
  )

  dev.off()
  invisible(circ)
}

# Run GO for every comparison with enough DE+DT genes
for (nm in names(deg_dtr_tables)) {
  label <- if (!is.null(go_labels[[nm]])) go_labels[[nm]] else nm
  message("GO enrichment: ", nm, " -> ", label)
  run_go_enrichment(
    resEdgeRtopTags = deg_dtr_tables[[nm]],
    label = label,
    out_dir = out_dir,
    geneNames = geneNames,
    geneID2GO = geneID2GO
  )
}
