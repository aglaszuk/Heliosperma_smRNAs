# Title: "GO enrichment analysis of DTGs in CG data"
# Author: Aglaia Szukala
# Short description: Script to perform GO enrichment analysis of genes that are differentially targeted by sRNAs and possibly also differentially expressed

# Set project directory
# Set project root (portable; works after cloning the repository)
.root_helper <- c("scripts/00_set_project_root.R", "../scripts/00_set_project_root.R")
source(.root_helper[file.exists(.root_helper)][1])
rm(.root_helper)

# Output directory
out_dir <- "output/DT_CG/8_GOenrichment/"
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
library(GOSemSim)
library(org.At.tair.db)
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

# Lists of DEGs
pairs <- c("1", "3", "4", "5")

deg_list <- lapply(pairs, function(p) {
  read.table(file = paste0("data/DEGlists_CG_SzukalaEtAl2022/AvsM_", p, "_fdr0.05.txt"))
})
names(deg_list) <- pairs

# Lists of DTRs
dtr_up_list <- lapply(pairs, function(p) {
  read.table(file = paste0("output/DT_CG/4_DT_edgeR/MvsA", p, "_DTRs_qlTest_poslogFC.txt"))
})
names(dtr_up_list) <- pairs

dtr_dw_list <- lapply(pairs, function(p) {
  read.table(file = paste0("output/DT_CG/4_DT_edgeR/MvsA", p, "_DTRs_qlTest_neglogFC.txt"))
})
names(dtr_dw_list) <- pairs

# Load region2gene mapping file
reg2gene <- read.table(file = "data/gene2region.uniq.txt")
colnames(reg2gene) <- c("gene", "region")

####################################
# Collapse DTRs to gene-level DTGs #
####################################

# Define function to match DTRs to genes
get_DEG_DTR_table <- function(pair_id, deg_list, dtr_up_list, dtr_dw_list, reg2gene) {
  
  deg_tab <- deg_list[[pair_id]]
  deg_tab$gene <- rownames(deg_tab)
  
  dtr_up <- dtr_up_list[[pair_id]]
  dtr_up$region <- rownames(dtr_up)
  dtr_up$DTR_direction <- "DTR_Mup"
  
  dtr_dw <- dtr_dw_list[[pair_id]]
  dtr_dw$region <- rownames(dtr_dw)
  dtr_dw$DTR_direction <- "DTR_Mdw"
  
  dtr_all <- rbind(dtr_up, dtr_dw)
  
  ## map DTR regions to genes
  dtr_gene <- merge(reg2gene, dtr_all, by = "region")
  
  ## merge with DEG table
  out <- merge(deg_tab, dtr_gene, by = "gene")
  
  ## cleaner column names
  colnames(out) <- gsub("\\.x$", "_DEG", colnames(out))
  colnames(out) <- gsub("\\.y$", "_DTR", colnames(out))
  
  ## reorder columns
  out <- out[, c(
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
  )]
  
  out
}

# Run for all pairs
deg_dtr_tables <- lapply(pairs, function(p) {
  get_DEG_DTR_table(
    pair_id = p,
    deg_list = deg_list,
    dtr_up_list = dtr_up_list,
    dtr_dw_list = dtr_dw_list,
    reg2gene = reg2gene
  )
})

names(deg_dtr_tables) <- pairs

# Inspect output
head(deg_dtr_tables[["1"]])
head(deg_dtr_tables[["5"]]) # no matches in pair 5

# Save one file per pair
for (p in pairs) {
  write.table(
    deg_dtr_tables[[p]],
    file = file.path(out_dir, paste0("DEG_with_DTR_pair", p, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

#############################################
# RUN GO TERM ENRICHMENT of genes DE and DT #
#############################################

# Pair 1
resEdgeRtopTags <- deg_dtr_tables[["1"]]

# Retain only the region with the lowest FDR per gene
resEdgeRtopTags <- resEdgeRtopTags %>%
  arrange(gene, padj) %>%
  distinct(gene, .keep_all = TRUE)

# Define genes for GO test
deGenes <- unique(resEdgeRtopTags$gene)

# GO analysis
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
  
  # IMPORTANT: use DEG values here, as these are biologically more meaningfull
  EC.genelist = data.frame(
    ID = resEdgeRtopTags$gene,
    logFC = resEdgeRtopTags$logFC_DEG,
    logCPM = resEdgeRtopTags$logCPM_DEG,
    P.Value = resEdgeRtopTags$PValue_DEG,
    adj.P.Val = resEdgeRtopTags$padj
  ),
  
  circ = circle_dat(go.dataframe, EC.genelist)
)

# Run GO analysis
make(plan)

# Retrieve output circ object
circ <- readd(circ)

# Save results
write.table(
  circ,
  file = file.path(out_dir, "Genes_DE_DT_pair1_GO_circ.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  col.names = TRUE
)

# Plot GO terms enrichment
pdf(
  file = file.path(out_dir, "Genes_DE_DT_pair1_GO.pdf"),
  width = 12,
  height = 12
)

ggplot(circ[circ$adj_pval < 0.05 & circ$count > 1 & circ$category != "CC" & circ$category != "MF", ],
       aes(x = term, y = -log10(adj_pval))) +
  geom_col(aes(fill = zscore)) + #category
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

dev.off()

#############################################################
# Pair 3
resEdgeRtopTags <- deg_dtr_tables[["3"]]

# Retain only the region with the lowest FDR per gene
resEdgeRtopTags <- resEdgeRtopTags %>%
  arrange(gene, padj) %>%
  distinct(gene, .keep_all = TRUE)

# Define genes for GO test
deGenes <- unique(resEdgeRtopTags$gene)

# GO analysis
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
  
  # IMPORTANT: use DEG values here, as these are biologically more meaningfull
  EC.genelist = data.frame(
    ID = resEdgeRtopTags$gene,
    logFC = resEdgeRtopTags$logFC_DEG,
    logCPM = resEdgeRtopTags$logCPM_DEG,
    P.Value = resEdgeRtopTags$PValue_DEG,
    adj.P.Val = resEdgeRtopTags$padj
  ),
  
  circ = circle_dat(go.dataframe, EC.genelist)
)

# Run GO analysis
make(plan)

# Retrieve output circ object
circ <- readd(circ)

# Save results
write.table(
  circ,
  file = file.path(out_dir, "Genes_DE_DT_pair3_GO_circ.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  col.names = TRUE
)

# Plot GO terms enrichment
pdf(
  file = file.path(out_dir, "Genes_DE_DT_pair3_GO.pdf"),
  width = 12,
  height = 12
)

ggplot(circ[circ$adj_pval < 0.05 & circ$count > 1 & circ$category != "CC" & circ$category != "MF", ],
       aes(x = term, y = -log10(adj_pval))) +
  geom_col(aes(fill = zscore)) + #category
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

dev.off()

###############################################################
# Pair 4
resEdgeRtopTags <- deg_dtr_tables[["4"]]

# Retain only the region with the lowest FDR per gene
resEdgeRtopTags <- resEdgeRtopTags %>%
  arrange(gene, padj) %>%
  distinct(gene, .keep_all = TRUE)

# Define genes for GO test
deGenes <- unique(resEdgeRtopTags$gene)

# GO analysis
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
  
  # IMPORTANT: use DEG values here, as these are biologically more meaningfull
  EC.genelist = data.frame(
    ID = resEdgeRtopTags$gene,
    logFC = resEdgeRtopTags$logFC_DEG,
    logCPM = resEdgeRtopTags$logCPM_DEG,
    P.Value = resEdgeRtopTags$PValue_DEG,
    adj.P.Val = resEdgeRtopTags$padj
  ),
  
  circ = circle_dat(go.dataframe, EC.genelist)
)

# Run GO analysis
make(plan)

# Retrieve output circ object
circ <- readd(circ)

# Save results
write.table(
  circ,
  file = file.path(out_dir, "Genes_DE_DT_pair4_GO_circ.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  col.names = TRUE
)

# Plot GO terms enrichment
pdf(
  file = file.path(out_dir, "Genes_DE_DT_pair4_GO.pdf"),
  width = 12,
  height = 12
)

ggplot(circ[circ$adj_pval < 0.05 & circ$count > 1 & circ$category != "CC" & circ$category != "MF", ],
       aes(x = term, y = -log10(adj_pval))) +
  geom_col(aes(fill = zscore)) + #category
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

dev.off()

######################################
# RUN GO TERM ENRICHMENT of all DTGs #
######################################
collapse_dtr_to_gene <- function(dtr_up_tab, dtr_dw_tab, reg2gene_tab) {

  dtr_all <- rbind(dtr_up_tab, dtr_dw_tab)
  dtr_all$region <- rownames(dtr_all)

  merged <- merge(reg2gene_tab, dtr_all, by = "region")
  merged <- merged[order(merged$gene, merged$FDR), ]

  gene_level <- merged[!duplicated(merged$gene), ]
  rownames(gene_level) <- gene_level$gene

  gene_level
}

dtg_gene_tables <- lapply(pairs, function(p) {
  collapse_dtr_to_gene(
    dtr_up_tab = dtr_up_list[[p]],
    dtr_dw_tab = dtr_dw_list[[p]],
    reg2gene_tab = reg2gene
  )
})
names(dtg_gene_tables) <- pairs

run_dtg_go_enrichment <- function(pair_id, dtg_tab, out_dir, geneNames, geneID2GO, name_list) {

  deGenes <- rownames(dtg_tab)

  plan_dtg <- drake_plan(
    name_list = name_list,

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

    ## Use DTR values from the gene-level DTG table
    EC.genelist = data.frame(
      ID = dtg_tab$gene,
      logFC = dtg_tab$logFC,
      logCPM = dtg_tab$logCPM,
      P.Value = dtg_tab$PValue,
      adj.P.Val = dtg_tab$FDR
    ),

    circ = circle_dat(go.dataframe, EC.genelist)
  )

  make(plan_dtg)
  circ <- readd(circ)

  write.table(
    dtg_tab,
    file = file.path(out_dir, paste0("DTG_pair", pair_id, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    circ,
    file = file.path(out_dir, paste0("DTG_pair", pair_id, "_GO_circ.txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    col.names = TRUE
  )

  pdf(
    file = file.path(out_dir, paste0("DTG_pair", pair_id, "_GO.pdf")),
    width = 12,
    height = 12
  )

  print(
    ggplot(
      circ[circ$adj_pval < 0.05 & circ$count > 1 & circ$category != "CC" & circ$category != "MF", ],
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
      labs(title = paste("DTG GO enrichment - pair", pair_id)) +
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

name_list <- c("GO.ID", "Term", "Annotated", "Significant",
               "Expected", "weight01_pval", "branch")

dtg_go_circ_by_pair <- list()

for (p in pairs) {
  message("Running DTG GO enrichment for pair ", p)
  dtg_go_circ_by_pair[[p]] <- run_dtg_go_enrichment(
    pair_id = p,
    dtg_tab = dtg_gene_tables[[p]],
    out_dir = out_dir,
    geneNames = geneNames,
    geneID2GO = geneID2GO,
    name_list = name_list
  )
}

#############################################################################################
# Test genetic redundancy by detecting similar GO terms enriched for different sets of DTGs #
#############################################################################################

# Note that for Figure 5 of the manuscript I redefined the overarching biological process" including a group of semantic similar GO terms.
# I also decided to retain for the figure only groups containing >= 6 semantic similar GO terms
# For this reason, the final figure differs from the output of this script

# Group related BP terms by Wang semantic similarity, then report themes enriched in >= 3 pairs.
go_sim_threshold <- 0.35   # min Wang similarity to group terms in one cluster
min_pairs_theme <- 3       # theme must be enriched in at least this many pairs

build_go_sim_matrix <- function(go_ids, go_data) {
  go_ids <- unique(go_ids)
  n <- length(go_ids)
  sim <- matrix(0, n, n, dimnames = list(go_ids, go_ids))
  diag(sim) <- 1
  if (n < 2) return(sim)

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      s <- tryCatch(
        goSim(go_ids[i], go_ids[j], semData = go_data, measure = "Wang"),
        error = function(e) 0
      )
      if (is.na(s)) s <- 0
      sim[i, j] <- s
      sim[j, i] <- s
    }
  }
  sim
}

cluster_go_by_semantic_sim <- function(sim_mat, sim_threshold = go_sim_threshold) {
  go_ids <- rownames(sim_mat)
  if (length(go_ids) == 0) {
    return(data.frame(ID = character(0), cluster = integer(0), stringsAsFactors = FALSE))
  }
  if (length(go_ids) == 1) {
    return(data.frame(ID = go_ids, cluster = 1L, stringsAsFactors = FALSE))
  }

  hc <- hclust(as.dist(1 - sim_mat), method = "average")
  clusters <- cutree(hc, h = 1 - sim_threshold)
  data.frame(
    ID = names(clusters),
    cluster = as.integer(clusters),
    stringsAsFactors = FALSE
  )
}

summarize_semantic_clusters <- function(cross_pair_annotated, min_pairs = min_pairs_theme) {
  cluster_summary <- cross_pair_annotated %>%
    group_by(cluster) %>%
    summarise(
      n_pairs = n_distinct(pair),
      pairs = paste(sort(unique(pair)), collapse = "_"),
      n_terms = n_distinct(ID),
      member_ids = paste(sort(unique(ID)), collapse = "; "),
      member_terms = paste(sort(unique(term)), collapse = " | "),
      mean_neglog10p = mean(-log10(adj_pval)),
      .groups = "drop"
    ) %>%
    filter(n_pairs >= min_pairs)

  if (nrow(cluster_summary) == 0) return(cluster_summary)

  rep_terms <- cross_pair_annotated %>%
    group_by(cluster, ID, term) %>%
    summarise(mean_p = mean(adj_pval), .groups = "drop") %>%
    group_by(cluster) %>%
    slice_min(mean_p, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::rename(theme_id = ID, theme_term = term)

  cluster_summary %>%
    left_join(rep_terms, by = "cluster") %>%
    mutate(theme_label = paste0(theme_term, " (+", pmax(n_terms - 1L, 0L), " related terms)")) %>%
    arrange(desc(n_pairs), desc(mean_neglog10p))
}

build_cluster_pair_matrix <- function(cross_pair_annotated, cluster_summary) {
  cross_pair_annotated %>%
    filter(cluster %in% cluster_summary$cluster) %>%
    group_by(cluster, pair) %>%
    summarise(
      neglog10p = max(-log10(adj_pval)),
      best_count = max(count),
      n_terms_in_pair = n_distinct(ID),
      .groups = "drop"
    ) %>%
    left_join(
      cluster_summary %>%
        dplyr::select(cluster, theme_label, theme_id, theme_term, n_pairs, member_terms),
      by = "cluster"
    )
}

plot_semantic_theme_heatmap <- function(cluster_pair_mat,
                                       title,
                                       pair_levels,
                                       pair_colors = go_pair_colors) {
  if (nrow(cluster_pair_mat) == 0) {
    message("No semantic themes to plot for: ", title)
    return(invisible(NULL))
  }

  theme_order <- cluster_pair_mat %>%
    distinct(cluster, theme_label, n_pairs) %>%
    arrange(desc(n_pairs), theme_label) %>%
    pull(theme_label)

  plot_data <- cluster_pair_mat %>%
    mutate(
      theme_label = factor(theme_label, levels = rev(theme_order)),
      pair = factor(pair, levels = pair_levels)
    )

  ggplot(plot_data, aes(x = pair, y = theme_label, fill = neglog10p)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = n_terms_in_pair), size = 3) +
    scale_fill_gradient(low = "#f7f7f7", high = "#2166ac", name = expression(-log[10](p))) +
    labs(
      title = title,
      x = "Ecotype pair",
      y = NULL,
      subtitle = "Tile labels = number of related GO terms enriched in that pair"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )
}

parse_go_gene_list <- function(gene_string) {
  if (is.na(gene_string) || !nzchar(gene_string)) return(character(0))
  unique(unlist(strsplit(gene_string, ",\\s*")))
}

build_cluster_pair_gene_tables <- function(cross_pair_annotated,
                                           cluster_summary,
                                           dtg_gene_tables,
                                           pair_levels) {
  if (nrow(cluster_summary) == 0) {
    empty_summary <- data.frame(
      cluster = integer(0),
      pair = character(0),
      theme_label = character(0),
      theme_id = character(0),
      theme_term = character(0),
      n_genes = integer(0),
      mean_logFC = numeric(0),
      median_logFC = numeric(0),
      genes = character(0),
      stringsAsFactors = FALSE
    )
    empty_detail <- data.frame(
      cluster = integer(0),
      pair = character(0),
      theme_label = character(0),
      gene = character(0),
      logFC = numeric(0),
      stringsAsFactors = FALSE
    )
    return(list(summary = empty_summary, detail = empty_detail))
  }

  theme_meta <- cluster_summary %>%
    dplyr::select(cluster, theme_label, theme_id, theme_term, n_pairs, member_terms)

  summary_rows <- list()
  detail_rows <- list()

  for (cl in theme_meta$cluster) {
    for (p in pair_levels) {
      sub <- cross_pair_annotated %>%
        filter(cluster == cl, pair == p)
      if (nrow(sub) == 0) next

      gene_vec <- unique(unlist(lapply(sub$genes, parse_go_gene_list)))
      if (length(gene_vec) == 0) next

      dtg_tab <- dtg_gene_tables[[p]]
      logfc <- dtg_tab$logFC[match(gene_vec, dtg_tab$gene)]
      names(logfc) <- gene_vec

      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        cluster = cl,
        pair = p,
        n_genes = length(gene_vec),
        mean_logFC = mean(logfc, na.rm = TRUE),
        median_logFC = stats::median(logfc, na.rm = TRUE),
        genes = paste(gene_vec, collapse = ", "),
        stringsAsFactors = FALSE
      )

      detail_rows[[length(detail_rows) + 1L]] <- data.frame(
        cluster = cl,
        pair = p,
        gene = gene_vec,
        logFC = as.numeric(logfc),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(summary_rows) == 0) {
    empty_summary <- data.frame(
      cluster = integer(0),
      pair = character(0),
      theme_label = character(0),
      theme_id = character(0),
      theme_term = character(0),
      n_genes = integer(0),
      mean_logFC = numeric(0),
      median_logFC = numeric(0),
      genes = character(0),
      stringsAsFactors = FALSE
    )
    empty_detail <- data.frame(
      cluster = integer(0),
      pair = character(0),
      theme_label = character(0),
      gene = character(0),
      logFC = numeric(0),
      stringsAsFactors = FALSE
    )
    return(list(summary = empty_summary, detail = empty_detail))
  }

  summary_tbl <- bind_rows(summary_rows) %>%
    left_join(theme_meta, by = "cluster") %>%
    dplyr::select(
      cluster, pair, theme_label, theme_id, theme_term, n_pairs, n_genes,
      mean_logFC, median_logFC, genes, member_terms
    ) %>%
    arrange(cluster, pair)

  detail_tbl <- bind_rows(detail_rows) %>%
    left_join(theme_meta %>% dplyr::select(cluster, theme_label, theme_id, theme_term), by = "cluster") %>%
    dplyr::select(cluster, pair, theme_label, theme_id, theme_term, gene, logFC) %>%
    arrange(cluster, pair, gene)

  list(summary = summary_tbl, detail = detail_tbl)
}

plot_semantic_theme_logfc_heatmap <- function(cluster_gene_summary,
                                              title,
                                              pair_levels) {
  if (nrow(cluster_gene_summary) == 0) {
    message("No semantic theme gene data to plot for: ", title)
    return(invisible(NULL))
  }

  theme_order <- cluster_gene_summary %>%
    distinct(cluster, theme_label, n_pairs) %>%
    arrange(desc(n_pairs), theme_label) %>%
    pull(theme_label)

  plot_data <- cluster_gene_summary %>%
    mutate(
      theme_label = factor(theme_label, levels = rev(theme_order)),
      pair = factor(pair, levels = pair_levels)
    )

  fill_limit <- max(abs(plot_data$mean_logFC), na.rm = TRUE)
  if (!is.finite(fill_limit) || fill_limit == 0) fill_limit <- 1

  ggplot(plot_data, aes(x = pair, y = theme_label, fill = mean_logFC)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = n_genes), size = 3) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      name = "Mean logFC"
    ) +
    labs(
      title = title,
      x = "Ecotype pair",
      y = NULL,
      subtitle = "Tile color = mean logFC of unique DTG genes; labels = gene count per cluster"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )
}

message("Building GO semantic similarity data...")
go_semdata <- godata(org.At.tair.db, keytype = "TAIR", ont = "BP")

unique_go_ids <- unique(dtg_cross_pair_go$ID)
go_sim_mat <- build_go_sim_matrix(unique_go_ids, go_semdata)
go_cluster_assignments <- cluster_go_by_semantic_sim(
  sim_mat = go_sim_mat,
  sim_threshold = go_sim_threshold
)

go_term_labels <- dtg_cross_pair_go %>%
  distinct(ID, term)

dtg_cross_pair_go_clusters <- dtg_cross_pair_go %>%
  inner_join(go_cluster_assignments, by = "ID")

dtg_semantic_clusters <- summarize_semantic_clusters(
  cross_pair_annotated = dtg_cross_pair_go_clusters,
  min_pairs = min_pairs_theme
)

print(dtg_semantic_clusters)

write.table(
  go_cluster_assignments %>%
    left_join(go_term_labels, by = "ID") %>%
    arrange(cluster, ID),
  file = file.path(out_dir, "DTG_semantic_GO_clusters.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  go_sim_mat,
  file = file.path(out_dir, "DTG_semantic_GO_sim_matrix.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  dtg_semantic_clusters,
  file = file.path(out_dir, "DTG_shared_GO_themes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(dtg_semantic_clusters) > 0) {
  dtg_theme_pair_mat <- build_cluster_pair_matrix(
    cross_pair_annotated = dtg_cross_pair_go_clusters,
    cluster_summary = dtg_semantic_clusters
  )

  write.table(
    dtg_theme_pair_mat,
    file = file.path(out_dir, "DTG_shared_GO_themes_by_pair.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  pdf(
    file.path(out_dir, "DTG_shared_GO_themes_heatmap.pdf"),
    width = 11,
    height = max(6, 0.35 * nrow(dtg_semantic_clusters) + 2)
  )
  print(plot_semantic_theme_heatmap(
    cluster_pair_mat = dtg_theme_pair_mat,
    title = paste0(
      "DTG BP themes (semantic clusters) enriched in >= ",
      min_pairs_theme, " pairs (Wang >= ", go_sim_threshold, ")"
    ),
    pair_levels = pairs,
    pair_colors = go_pair_colors
  ))
  dev.off()

  dtg_theme_gene_tables <- build_cluster_pair_gene_tables(
    cross_pair_annotated = dtg_cross_pair_go_clusters,
    cluster_summary = dtg_semantic_clusters,
    dtg_gene_tables = dtg_gene_tables,
    pair_levels = pairs
  )

  write.table(
    dtg_theme_gene_tables$summary,
    file = file.path(out_dir, "DTG_shared_GO_themes_genes_by_pair.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    dtg_theme_gene_tables$detail,
    file = file.path(out_dir, "DTG_shared_GO_themes_genes_detail.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  pdf(
    file.path(out_dir, "DTG_shared_GO_themes_heatmap_logFC.pdf"),
    width = 11,
    height = max(6, 0.35 * nrow(dtg_semantic_clusters) + 2)
  )
  print(plot_semantic_theme_logfc_heatmap(
    cluster_gene_summary = dtg_theme_gene_tables$summary,
    title = paste0(
      "DTG BP themes (semantic clusters) enriched in >= ",
      min_pairs_theme, " pairs - mean logFC of underlying genes"
    ),
    pair_levels = pairs
  ))
  dev.off()
} else {
  message(
    "No semantic GO themes found in >= ", min_pairs_theme,
    " pairs at similarity threshold ", go_sim_threshold, "."
  )
}
