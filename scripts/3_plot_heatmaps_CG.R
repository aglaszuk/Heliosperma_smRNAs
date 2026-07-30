# Title: "Plot heatmaps of genes DE and DT by small RNAs in a common garden (CG) experimental setup for pairs 1, 3 and 4"
# Author: Aglaia Szukala
# Short description: Script needs lists of DTRs and DEGs for each ecotype comparison and a table mapping genes to genomic regions.

# Load Libraries
library(ComplexHeatmap)
library(circlize)
library(grid)

# Set working directory
# Set project root (portable; works after cloning the repository)
.root_helper <- c("scripts/00_set_project_root.R", "../scripts/00_set_project_root.R")
source(.root_helper[file.exists(.root_helper)][1])
rm(.root_helper)

# Define output directory
out_dir <- "output/DT_CG/7_Heatmaps_DT_DEGs/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

###################
# Load input data #
###################

# Table mapping regions to a gene
reg2gene <- read.table(file = "data/gene2region.uniq.txt")
colnames(reg2gene) <- c("gene", "region")

# Lists of DEGs
deg_cg_1 <- read.table(file = "data/DEGlists_CG_SzukalaEtAl2022/AvsM_1_fdr0.05.txt")
deg_cg_3 <- read.table(file = "data/DEGlists_CG_SzukalaEtAl2022/AvsM_3_fdr0.05.txt")
deg_cg_4 <- read.table(file = "data/DEGlists_CG_SzukalaEtAl2022/AvsM_4_fdr0.05.txt")

# Lists of DTRs
dtr_cg_1up <- read.table(file = "output/DT_CG/4_DT_edgeR/MvsA1_DTRs_qlTest_poslogFC.txt")
dtr_cg_1dw <- read.table(file = "output/DT_CG/4_DT_edgeR/MvsA1_DTRs_qlTest_neglogFC.txt")

dtr_cg_3up <- read.table(file = "output/DT_CG/4_DT_edgeR/MvsA3_DTRs_qlTest_poslogFC.txt")
dtr_cg_3dw <- read.table(file = "output/DT_CG/4_DT_edgeR/MvsA3_DTRs_qlTest_neglogFC.txt")

dtr_cg_4up <- read.table(file = "output/DT_CG/4_DT_edgeR/MvsA4_DTRs_qlTest_poslogFC.txt")
dtr_cg_4dw <- read.table(file = "output/DT_CG/4_DT_edgeR/MvsA4_DTRs_qlTest_neglogFC.txt")

# Make pair-specific tables
deg_list <- list(
  "1" = deg_cg_1,
  "3" = deg_cg_3,
  "4" = deg_cg_4
)

dtr_up_list <- list(
  "1" = dtr_cg_1up,
  "3" = dtr_cg_3up,
  "4" = dtr_cg_4up
)

dtr_dw_list <- list(
  "1" = dtr_cg_1dw,
  "3" = dtr_cg_3dw,
  "4" = dtr_cg_4dw
)

####################################
# Collapse DTRs to gene-level DTGs #
####################################

# Define DTR2DTG function
collapse_dtr_to_gene <- function(dtr_up_tab, dtr_dw_tab, reg2gene_tab) {
  
  dtr_all <- rbind(dtr_up_tab, dtr_dw_tab)
  dtr_all$region <- rownames(dtr_all)
  
  merged <- merge(reg2gene_tab, dtr_all, by = "region")
  merged <- merged[order(merged$gene, merged$FDR), ]
  
  gene_level <- merged[!duplicated(merged$gene), ]
  rownames(gene_level) <- gene_level$gene
  
  gene_level
}

# Define pairs to be used
pairs_use <- c("1", "3", "4")

# Retrieve DTGs tables
dtg_gene_tables <- lapply(pairs_use, function(p) {
  collapse_dtr_to_gene(
    dtr_up_tab = dtr_up_list[[p]],
    dtr_dw_tab = dtr_dw_list[[p]],
    reg2gene_tab = reg2gene
  )
})
names(dtg_gene_tables) <- pairs_use

#######################################
# Get anti-directional genes per pair #
# Anti-directional:                   #
# expression > 0 and targeting < 0    #
# expression < 0 and targeting > 0    #
#######################################

get_anti_genes_one_pair <- function(deg_tab, dtg_tab) {
  
  common_genes <- intersect(rownames(deg_tab), rownames(dtg_tab))
  
  expr_logFC <- deg_tab[common_genes, "logFC"]
  targ_logFC <- dtg_tab[common_genes, "logFC"]
  
  anti <- (expr_logFC > 0 & targ_logFC < 0) |
    (expr_logFC < 0 & targ_logFC > 0)
  
  common_genes[anti]
}

anti_genes_by_pair <- lapply(pairs_use, function(p) {
  get_anti_genes_one_pair(
    deg_tab = deg_list[[p]],
    dtg_tab = dtg_gene_tables[[p]]
  )
})
names(anti_genes_by_pair) <- pairs_use

cat("Number of anti-directional genes per pair:\n")
print(sapply(anti_genes_by_pair, length))

#####################################
# Get same-direction genes per pair #
# Same-direction:                   #
# expression > 0 and targeting > 0  #
# expression < 0 and targeting < 0  #
#####################################

get_same_genes_one_pair <- function(deg_tab, dtg_tab) {
  
  common_genes <- intersect(rownames(deg_tab), rownames(dtg_tab))
  
  expr_logFC <- deg_tab[common_genes, "logFC"]
  targ_logFC <- dtg_tab[common_genes, "logFC"]
  
  sameg <- (expr_logFC > 0 & targ_logFC > 0) |
    (expr_logFC < 0 & targ_logFC < 0)
  
  common_genes[sameg]
}

same_genes_by_pair <- lapply(pairs_use, function(p) {
  get_same_genes_one_pair(
    deg_tab = deg_list[[p]],
    dtg_tab = dtg_gene_tables[[p]]
  )
})
names(same_genes_by_pair) <- pairs_use

cat("Number of same-direction genes per pair:\n")
print(sapply(same_genes_by_pair, length))

###############################################
# Identify genes shared by at least two pairs #
###############################################

# List all antigenes and same genes
all_anti_genes <- sort(unique(unlist(anti_genes_by_pair)))
all_same_genes <- sort(unique(unlist(same_genes_by_pair)))

# Per-pair DE+DT genes (any direction), for shared-gene counting
de_dt_by_pair <- lapply(pairs_use, function(p) {
  unique(c(anti_genes_by_pair[[p]], same_genes_by_pair[[p]]))
})
names(de_dt_by_pair) <- pairs_use

all_de_dt_genes <- sort(unique(unlist(de_dt_by_pair)))

gene_pair_count <- sapply(all_de_dt_genes, function(g) {
  sum(vapply(de_dt_by_pair, function(x) g %in% x, logical(1)))
})

# Keep only genes present in at least 2 pairs
shared_genes <- names(gene_pair_count)[gene_pair_count >= 2]

cat("Genes DE and DT in at least two pairs:", length(shared_genes), "\n")

# Save as table
write.table(
  shared_genes,
  file = file.path(out_dir, "DE_DT_genes_shared_by_at_least_two_pairs.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

#################################
# Build heatmap matrix per pair #
#################################

make_pair_matrix_pair_specific <- function(pair_id) {
  
  genes_pair <- sort(c(anti_genes_by_pair[[pair_id]], same_genes_by_pair[[pair_id]]))
  
  mat <- matrix(
    NA_real_,
    nrow = length(genes_pair),
    ncol = 2,
    dimnames = list(genes_pair, c("Expression", "Targeting"))
  )
  
  mat[genes_pair, "Expression"] <- deg_list[[pair_id]][genes_pair, "logFC"]
  mat[genes_pair, "Targeting"] <- dtg_gene_tables[[pair_id]][genes_pair, "logFC"]
  
  mat
}

mat_pair1 <- make_pair_matrix_pair_specific("1")
mat_pair3 <- make_pair_matrix_pair_specific("3")
mat_pair4 <- make_pair_matrix_pair_specific("4")

########################################
# Order genes by type within each pair #
########################################

order_pair_matrix <- function(mat) {
  
  expr <- mat[, "Expression"]
  targ <- mat[, "Targeting"]
  
  pattern <- dplyr::case_when(
    expr > 0 & targ < 0 ~ "Mup_expr_Mdw_targ",
    expr < 0 & targ > 0 ~ "Mdw_expr_Mup_targ",
    expr < 0 & targ < 0 ~ "Mdw_expr_Mdw_targ",
    expr > 0 & targ > 0 ~ "Mup_expr_Mup_targ",
    TRUE ~ NA_character_
  )
  
  ord <- order(
    pattern,
    !(rownames(mat) %in% shared_genes),
    rownames(mat)
  )
  
  mat[ord, , drop = FALSE]
}

mat_pair1 <- order_pair_matrix(mat_pair1)
mat_pair3 <- order_pair_matrix(mat_pair3)
mat_pair4 <- order_pair_matrix(mat_pair4)

#################
# Plot heatmaps #
#################

# Define color scale
max_abs <- max(abs(c(mat_pair1, mat_pair3, mat_pair4)), na.rm = TRUE)

col_fun <- colorRamp2(
  c(-max_abs, -max_abs * 0.2, 0, max_abs * 0.2, max_abs),
  c("#2166ac", "#67a9cf", "white", "#ef8a62", "#b2182b")
)

# Row annotations per pair
make_pattern_annotation <- function(mat) {
  
  expr <- mat[, "Expression"]
  targ <- mat[, "Targeting"]
 
  pattern <- dplyr::case_when(
    expr > 0 & targ < 0 ~ "Mup_expr_Mdw_targ",
    expr < 0 & targ > 0 ~ "Mdw_expr_Mup_targ",
    expr < 0 & targ < 0 ~ "Mdw_expr_Mdw_targ",
    expr > 0 & targ > 0 ~ "Mup_expr_Mup_targ",
    TRUE ~ NA_character_
  )
  
  rowAnnotation(
    pattern = pattern,
    shared = ifelse(rownames(mat) %in% shared_genes, "shared", "pair_specific"),
    col = list(
      pattern = c(
        Mup_expr_Mdw_targ = "#c1121f",
        Mdw_expr_Mup_targ = "#15616d",
        Mup_expr_Mup_targ = "#57651B",
        Mdw_expr_Mdw_targ = "#6B3A46"
      ),
      shared = c(
        shared = "black",
        pair_specific = "grey80"
      )
    ),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 8)
  )
}

ha_pair1 <- make_pattern_annotation(mat_pair1)
ha_pair3 <- make_pattern_annotation(mat_pair3)
ha_pair4 <- make_pattern_annotation(mat_pair4)

# Gene name: bold shared genes
row_gp_pair1 <- gpar(
  fontsize = 6,
  fontface = ifelse(rownames(mat_pair1) %in% shared_genes, "bold", "plain")
)

row_gp_pair3 <- gpar(
  fontsize = 6,
  fontface = ifelse(rownames(mat_pair3) %in% shared_genes, "bold", "plain")
)

row_gp_pair4 <- gpar(
  fontsize = 6,
  fontface = ifelse(rownames(mat_pair4) %in% shared_genes, "bold", "plain")
)

# Build heatmaps: one heatmap per pair
ht_pair1 <- Heatmap(
  mat_pair1,
  name = "logFC",
  col = col_fun,
  left_annotation = ha_pair1,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = row_gp_pair1,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 9),
  column_title = paste0("Ecotype pair 1 | n = ", nrow(mat_pair1)),
  na_col = "grey95",
  heatmap_legend_param = list(title = "logFC")
)

ht_pair3 <- Heatmap(
  mat_pair3,
  name = "logFC",
  col = col_fun,
  left_annotation = ha_pair3,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = row_gp_pair3,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 9),
  column_title = paste0("Ecotype pair 3 | n = ", nrow(mat_pair3)),
  na_col = "grey95",
  heatmap_legend_param = list(title = "logFC")
)

ht_pair4 <- Heatmap(
  mat_pair4,
  name = "logFC",
  col = col_fun,
  left_annotation = ha_pair4,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = row_gp_pair4,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 9),
  column_title = paste0("Ecotype pair 4 | n = ", nrow(mat_pair4)),
  na_col = "grey95",
  heatmap_legend_param = list(title = "logFC")
)

#ht_all <- ht_pair1 %v% ht_pair3 %v% ht_pair4

# Save figures
pdf(
  file.path(out_dir, "heatmap_pair1.pdf"),
  width = 8,
  height = 18
)

draw(
  ht_pair1,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
  #merge_legends = TRUE
)

dev.off()

pdf(
  file.path(out_dir, "heatmap_pair3.pdf"),
  width = 8,
  height = 10
)

draw(
  ht_pair3,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
  #merge_legends = TRUE
)

dev.off()

pdf(
  file.path(out_dir, "heatmaps_pair4.pdf"),
  width = 8,
  height = 32
)

draw(
  ht_pair4,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
  #merge_legends = TRUE
)

dev.off()

# Export tables with expression and targeting logFC 
write.table(
  mat_pair1,
  file = file.path(out_dir, "heatmap_pair1_DEplusDT_pair_specific.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  mat_pair3,
  file = file.path(out_dir, "heatmap_pair3_DEplusDT_pair_specific.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  mat_pair4,
  file = file.path(out_dir, "heatmap_pair4_DEplusDT_pair_specific.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)