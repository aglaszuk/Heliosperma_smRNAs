# Title: "Plot heatmaps of genes DE and DT by small RNAs in a RT experimental setup for pairs 1, 3"
# Author: Aglaia Szukala
# Short description: Script needs lists of DTRs and DEGs for each ecotype comparison and a table mapping genes to genomic regions.

# Load Libraries
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

# Set working directory
# Set project root (portable; works after cloning the repository)
.root_helper <- c("scripts/00_set_project_root.R", "../scripts/00_set_project_root.R")
source(.root_helper[file.exists(.root_helper)][1])
rm(.root_helper)

# Define output directory
out_dir <- "output/DT_RT/7_Heatmaps_DT_DEGs/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###################
# Load input data #
###################

# DEGs
deg_AM_1d <- read.table("data/DEGlists_RT_SzukalaEtAl2023/DEG_A1vsM1d.txt")
deg_AM_3d <- read.table("data/DEGlists_RT_SzukalaEtAl2023/DEG_A3vsM3d.txt")
deg_AM_1u <- read.table("data/DEGlists_RT_SzukalaEtAl2023/DEG_A1vsM1u.txt")
deg_AM_3u <- read.table("data/DEGlists_RT_SzukalaEtAl2023/DEG_A3vsM3u.txt")

deg_AA_1 <- read.table("data/DEGlists_RT_SzukalaEtAl2023/DEG_A1ud.txt")
deg_MM_1 <- read.table("data/DEGlists_RT_SzukalaEtAl2023/DEG_M1ud.txt")
deg_AA_3 <- read.table("data/DEGlists_RT_SzukalaEtAl2023/DEG_A3ud.txt")
deg_MM_3 <- read.table("data/DEGlists_RT_SzukalaEtAl2023/DEG_M3ud.txt")

# Region-to-gene table
reg2gene <- read.table("data/gene2region.uniq.txt")
colnames(reg2gene) <- c("gene", "region")

# DT regions: negative logFC
AvsM_1D <- read.table("output/DT_RT/4_DT_edgeR/AvsM_1D_DTRs_qlTest_neglogFC.txt")
AvsM_3D <- read.table("output/DT_RT/4_DT_edgeR/AvsM_3D_DTRs_qlTest_neglogFC.txt")
AvsM_1U <- read.table("output/DT_RT/4_DT_edgeR/AvsM_1U_DTRs_qlTest_neglogFC.txt")
AvsM_3U <- read.table("output/DT_RT/4_DT_edgeR/AvsM_3U_DTRs_qlTest_neglogFC.txt")
MvsM1   <- read.table("output/DT_RT/4_DT_edgeR/MvsM1_DTRs_qlTest_neglogFC.txt")
MvsM3   <- read.table("output/DT_RT/4_DT_edgeR/MvsM3_DTRs_qlTest_neglogFC.txt")
AvsA1   <- read.table("output/DT_RT/4_DT_edgeR/AvsA1_DTRs_qlTest_neglogFC.txt")
AvsA3   <- read.table("output/DT_RT/4_DT_edgeR/AvsA3_DTRs_qlTest_neglogFC.txt")

# DT regions: positive logFC
AvsM_1D_p <- read.table("output/DT_RT/4_DT_edgeR/AvsM_1D_DTRs_qlTest_poslogFC.txt")
AvsM_3D_p <- read.table("output/DT_RT/4_DT_edgeR/AvsM_3D_DTRs_qlTest_poslogFC.txt")
AvsM_1U_p <- read.table("output/DT_RT/4_DT_edgeR/AvsM_1U_DTRs_qlTest_poslogFC.txt")
AvsM_3U_p <- read.table("output/DT_RT/4_DT_edgeR/AvsM_3U_DTRs_qlTest_poslogFC.txt")
MvsM1_p   <- read.table("output/DT_RT/4_DT_edgeR/MvsM1_DTRs_qlTest_poslogFC.txt")
MvsM3_p   <- read.table("output/DT_RT/4_DT_edgeR/MvsM3_DTRs_qlTest_poslogFC.txt")
AvsA1_p   <- read.table("output/DT_RT/4_DT_edgeR/AvsA1_DTRs_qlTest_poslogFC.txt")
AvsA3_p   <- read.table("output/DT_RT/4_DT_edgeR/AvsA3_DTRs_qlTest_poslogFC.txt")
 
# Make contrast-specific tables

contrasts_use <- c(
  "AvsM_1D", "AvsM_3D", "AvsM_1U", "AvsM_3U",
  "AvsA1", "AvsA3", "MvsM1", "MvsM3"
)

deg_list <- list(
  AvsM_1D = deg_AM_1d,
  AvsM_3D = deg_AM_3d,
  AvsM_1U = deg_AM_1u,
  AvsM_3U = deg_AM_3u,
  AvsA1   = deg_AA_1,
  AvsA3   = deg_AA_3,
  MvsM1   = deg_MM_1,
  MvsM3   = deg_MM_3
)

dtr_up_list <- list(
  AvsM_1D = AvsM_1D_p,
  AvsM_3D = AvsM_3D_p,
  AvsM_1U = AvsM_1U_p,
  AvsM_3U = AvsM_3U_p,
  AvsA1   = AvsA1_p,
  AvsA3   = AvsA3_p,
  MvsM1   = MvsM1_p,
  MvsM3   = MvsM3_p
)

dtr_dw_list <- list(
  AvsM_1D = AvsM_1D,
  AvsM_3D = AvsM_3D,
  AvsM_1U = AvsM_1U,
  AvsM_3U = AvsM_3U,
  AvsA1   = AvsA1,
  AvsA3   = AvsA3,
  MvsM1   = MvsM1,
  MvsM3   = MvsM3
)

####################################
# Collapse DTRs to gene-level DTGs #
####################################

collapse_dtr_to_gene <- function(dtr_up_tab, dtr_dw_tab, reg2gene_tab) {
  
  dtr_all <- rbind(dtr_up_tab, dtr_dw_tab)
  dtr_all$region <- rownames(dtr_all)
  
  merged <- merge(reg2gene_tab, dtr_all, by = "region")
  merged <- merged[order(merged$gene, merged$FDR), ]
  
  gene_level <- merged[!duplicated(merged$gene), ]
  rownames(gene_level) <- gene_level$gene
  
  gene_level
}

dtg_gene_tables <- lapply(contrasts_use, function(x) {
  collapse_dtr_to_gene(
    dtr_up_tab = dtr_up_list[[x]],
    dtr_dw_tab = dtr_dw_list[[x]],
    reg2gene_tab = reg2gene
  )
})

names(dtg_gene_tables) <- contrasts_use

# Check one
head(dtg_gene_tables$AvsM_1D)
length(unique(dtg_gene_tables$AvsM_1D$gene))

#######################################
# Get anti-directional genes per pair #
# Anti-directional:                   #
# expression > 0 and targeting < 0    #
# expression < 0 and targeting > 0    #
#######################################

# Define function
get_anti_genes_one_contrast <- function(deg_tab, dtg_tab) {
  
  common_genes <- intersect(rownames(deg_tab), rownames(dtg_tab))
  
  expr_logFC <- deg_tab[common_genes, "logFC"]
  targ_logFC <- dtg_tab[common_genes, "logFC"]
  
  anti <- (expr_logFC > 0 & targ_logFC < 0) |
    (expr_logFC < 0 & targ_logFC > 0)
  
  common_genes[anti]
}

#Run
anti_genes_by_contrast <- lapply(contrasts_use, function(x) {
  get_anti_genes_one_contrast(
    deg_tab = deg_list[[x]],
    dtg_tab = dtg_gene_tables[[x]]
  )
})

names(anti_genes_by_contrast) <- contrasts_use

cat("Number of anti-directional genes per contrast:\n")
print(sapply(anti_genes_by_contrast, length))

#####################################
# Get same-direction genes per pair #
# Same-direction:                   #
# expression > 0 and targeting > 0  #
# expression < 0 and targeting < 0  #
#####################################

# Define function
get_same_genes_one_contrast <- function(deg_tab, dtg_tab) {
  
  common_genes <- intersect(rownames(deg_tab), rownames(dtg_tab))
  
  expr_logFC <- deg_tab[common_genes, "logFC"]
  targ_logFC <- dtg_tab[common_genes, "logFC"]
  
  sameg <- (expr_logFC > 0 & targ_logFC > 0) |
    (expr_logFC < 0 & targ_logFC < 0)
  
  common_genes[sameg]
}

# Run
same_genes_by_contrast <- lapply(contrasts_use, function(x) {
  get_same_genes_one_contrast(
    deg_tab = deg_list[[x]],
    dtg_tab = dtg_gene_tables[[x]]
  )
})

names(same_genes_by_contrast) <- contrasts_use

cat("Number of same-directional genes per contrast:\n")
print(sapply(same_genes_by_contrast, length))

###################################################
# Identify genes shared by at least two contrasts #
###################################################

all_anti_genes <- sort(unique(unlist(anti_genes_by_contrast)))
all_same_genes <- sort(unique(unlist(same_genes_by_contrast)))

# Per-contrast DE+DT genes (any direction), for shared-gene counting
de_dt_by_contrast <- lapply(contrasts_use, function(x) {
  unique(c(anti_genes_by_contrast[[x]], same_genes_by_contrast[[x]]))
})
names(de_dt_by_contrast) <- contrasts_use

all_de_dt_genes <- sort(unique(unlist(de_dt_by_contrast)))

gene_contrast_count <- sapply(all_de_dt_genes, function(g) {
  sum(vapply(de_dt_by_contrast, function(x) g %in% x, logical(1)))
})

shared_genes <- names(gene_contrast_count)[gene_contrast_count >= 2]

# Direction-specific contrast counts (do not mix anti and same)
anti_contrast_count <- sapply(all_anti_genes, function(g) {
  sum(vapply(anti_genes_by_contrast, function(x) g %in% x, logical(1)))
})

same_contrast_count <- sapply(all_same_genes, function(g) {
  sum(vapply(same_genes_by_contrast, function(x) g %in% x, logical(1)))
})

cat("Genes DE and DT in at least two contrasts:", length(shared_genes), "\n")

write.table(
  shared_genes,
  file = file.path(out_dir, "DE_DT_genes_shared_by_at_least_two_contrasts.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

#####################################
# Build heatmap matrix per contrats #
#####################################

make_contrast_matrix <- function(contrast_id) {
  
  genes_contrast <- sort(c(
    anti_genes_by_contrast[[contrast_id]],
    same_genes_by_contrast[[contrast_id]]
  ))
  
  genes_contrast <- unique(genes_contrast)
  
  mat <- matrix(
    NA_real_,
    nrow = length(genes_contrast),
    ncol = 2,
    dimnames = list(genes_contrast, c("Expression", "Targeting"))
  )
  
  mat[genes_contrast, "Expression"] <- deg_list[[contrast_id]][genes_contrast, "logFC"]
  mat[genes_contrast, "Targeting"] <- dtg_gene_tables[[contrast_id]][genes_contrast, "logFC"]
  
  mat
}

mat_list <- lapply(contrasts_use, make_contrast_matrix)
names(mat_list) <- contrasts_use
mat_list

# Order genes within each contrast
order_contrast_matrix <- function(mat) {
  
  if (is.null(mat) || nrow(mat) == 0) {
    return(mat)
  }
  
  expr <- mat[, "Expression"]
  targ <- mat[, "Targeting"]
  
  pattern <- case_when(
    expr > 0 & targ < 0 ~ "DEG_up_DTG_down",
    expr < 0 & targ > 0 ~ "DEG_down_DTG_up",
    expr > 0 & targ > 0 ~ "DEG_up_DTG_up",
    expr < 0 & targ < 0 ~ "DEG_down_DTG_down",
    TRUE ~ "unclassified"
  )
  
  pattern <- factor(
    pattern,
    levels = c(
      "DEG_up_DTG_down",
      "DEG_down_DTG_up",
      "DEG_up_DTG_up",
      "DEG_down_DTG_down",
      "unclassified"
    )
  )
  
  ord <- order(
    pattern,
    !(rownames(mat) %in% shared_genes),
    rownames(mat),
    na.last = TRUE
  )
  
  mat[ord, , drop = FALSE]
}

mat_list <- lapply(mat_list, order_contrast_matrix)

# Define shared color scale for logFC intensity
max_abs <- max(abs(unlist(mat_list)), na.rm = TRUE)

col_fun <- colorRamp2(
  c(-max_abs, -max_abs * 0.2, 0, max_abs * 0.2, max_abs),
  c("#2166ac", "#67a9cf", "white", "#ef8a62", "#b2182b")
)

# Define row annotations for plot
make_pattern_annotation <- function(mat) {
  
  if (is.null(mat) || nrow(mat) == 0) {
    return(NULL)
  }
  
  expr <- mat[, "Expression"]
  targ <- mat[, "Targeting"]
  
  pattern <- case_when(
    expr > 0 & targ < 0 ~ "DEG_up_DTG_down",
    expr < 0 & targ > 0 ~ "DEG_down_DTG_up",
    expr > 0 & targ > 0 ~ "DEG_up_DTG_up",
    expr < 0 & targ < 0 ~ "DEG_down_DTG_down",
    TRUE ~ "unclassified"
  )
  
  pattern <- factor(
    pattern,
    levels = c(
      "DEG_up_DTG_down",
      "DEG_down_DTG_up",
      "DEG_up_DTG_up",
      "DEG_down_DTG_down",
      "unclassified"
    )
  )
  
  rowAnnotation(
    pattern = pattern,
    shared = ifelse(rownames(mat) %in% shared_genes, "shared", "contrast_specific"),
    col = list(
      pattern = c(
        DEG_up_DTG_down = "#c1121f",
        DEG_down_DTG_up = "#15616d",
        DEG_up_DTG_up = "#57651B",
        DEG_down_DTG_down = "#6B3A46",
        unclassified = "grey90"
      ),
      shared = c(
        shared = "black",
        contrast_specific = "grey80"
      )
    ),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 8)
  )
}

ha_list <- lapply(mat_list, make_pattern_annotation)

# double-heck if matrix contains any genes
for (contrast in contrasts_use) {
  
  mat <- mat_list[[contrast]]
  
  if (is.null(mat) || nrow(mat) == 0) {
    message("Skipping ", contrast, ": no DE+DT genes")
    next
  }
  
  ha <- ha_list[[contrast]]
  
  row_gp <- gpar(
    fontsize = 6,
    fontface = ifelse(
      rownames(mat) %in% shared_genes,
      "bold",
      "plain"
    )
  )
  
  ht <- Heatmap(
    mat,
    name = "logFC",
    col = col_fun,
    left_annotation = ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = row_gp,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 9),
    column_title = paste0(contrast, " | n = ", nrow(mat)),
    na_col = "grey95",
    heatmap_legend_param = list(title = "logFC")
  )
  
  pdf(
    file.path(out_dir, paste0("heatmap_", contrast, "_DEplusDT.pdf")),
    width = 8,
    height = max(6, min(35, nrow(mat) * 0.12))
  )
  
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  
  dev.off()
  
  write.table(
    mat,
    file = file.path(out_dir, paste0("heatmap_", contrast, "_DEplusDT_matrix.txt")),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

names(which(sapply(mat_list, nrow) == 0))

###########################################
# Build and save one heatmap per contrast #
###########################################

for (contrast in contrasts_use) {
  
  mat <- mat_list[[contrast]]
  
  if (is.null(mat) || nrow(mat) == 0) {
    message("Skipping ", contrast, ": no DE+DT genes")
    next
  }
  
  ha <- ha_list[[contrast]]
  
  row_gp <- gpar(
    fontsize = 6,
    fontface = ifelse(rownames(mat) %in% shared_genes, "bold", "plain")
  )
  
  ht <- Heatmap(
    mat,
    name = "logFC",
    col = col_fun,
    left_annotation = ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = row_gp,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 9),
    column_title = paste0(contrast, " | n = ", nrow(mat)),
    na_col = "grey95",
    heatmap_legend_param = list(title = "logFC")
  )
  
  pdf(
    file.path(out_dir, paste0("heatmap_", contrast, "_DEplusDT.pdf")),
    width = 8,
    height = max(6, min(35, nrow(mat) * 0.12))
  )
  
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  
  dev.off()
  
  write.table(
    mat,
    file = file.path(out_dir, paste0("heatmap_", contrast, "_DEplusDT_matrix.txt")),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

###########################
# Export gene annotations #
###########################

gene_annotation_anti <- data.frame(
  gene = all_anti_genes,
  n_contrasts_with_anti_direction = anti_contrast_count[all_anti_genes],
  shared_across_contrasts = all_anti_genes %in% shared_genes
)

for (contrast in contrasts_use) {
  gene_annotation_anti[[paste0("in_", contrast)]] <-
    gene_annotation_anti$gene %in% anti_genes_by_contrast[[contrast]]
}

write.table(
  gene_annotation_anti,
  file = file.path(out_dir, "anti_directional_gene_annotation_all_RT_contrasts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_annotation_anti[gene_annotation_anti$n_contrasts_with_anti_direction > 1, ],
  file = file.path(out_dir, "anti_directional_genes_in_multiple_RT_contrasts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gene_annotation_same <- data.frame(
  gene = all_same_genes,
  n_contrasts_with_same_direction = same_contrast_count[all_same_genes],
  shared_across_contrasts = all_same_genes %in% shared_genes
)

for (contrast in contrasts_use) {
  gene_annotation_same[[paste0("in_", contrast)]] <-
    gene_annotation_same$gene %in% same_genes_by_contrast[[contrast]]
}

write.table(
  gene_annotation_same,
  file = file.path(out_dir, "same_directional_gene_annotation_all_RT_contrasts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_annotation_same[gene_annotation_same$n_contrasts_with_same_direction > 1, ],
  file = file.path(out_dir, "same_directional_genes_in_multiple_RT_contrasts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
