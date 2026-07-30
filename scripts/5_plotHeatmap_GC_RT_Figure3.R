# Title: "Plot heatmap of Figure 3 of Szukala et al. - smRNAs in Heliosperma"
# Author: Aglaia Szukala
# Short description: Heatmap of expression and sRNA-targeting logFC for genes
# present in >= 2 RT/CG comparisons, with curated annotations and consistent
# logFC direction in both Expression and Targeting.

# Load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Set project root (portable; works after cloning the repository)
.root_helper <- c("scripts/00_set_project_root.R", "../scripts/00_set_project_root.R")
source(.root_helper[file.exists(.root_helper)][1])
rm(.root_helper)

# Output directory
out_dir <- "output/Heatmap_Figure3"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Input DE+DT matrices: RT montane (1D, 3D) and CG pairs 1, 3, 4
files <- c(
  "RT-1D" = "output/DT_RT/7_Heatmaps_DT_DEGs/heatmap_AvsM_1D_DEplusDT_matrix.txt",
  "RT-3D" = "output/DT_RT/7_Heatmaps_DT_DEGs/heatmap_AvsM_3D_DEplusDT_matrix.txt",
  "CG 1"  = "output/DT_CG/7_Heatmaps_DT_DEGs/heatmap_pair1_DEplusDT_pair_specific.txt",
  "CG 3"  = "output/DT_CG/7_Heatmaps_DT_DEGs/heatmap_pair3_DEplusDT_pair_specific.txt",
  "CG 4"  = "output/DT_CG/7_Heatmaps_DT_DEGs/heatmap_pair4_DEplusDT_pair_specific.txt"
)

comparison_cols <- c(
  "RT-1D" = "#A71E34",
  "RT-3D" = "#00798C",
  "CG 1"  = "#A71E34",
  "CG 3"  = "#00798C",
  "CG 4"  = "#F59700"
)

type_cols <- c(
  "Expression" = "black",
  "Targeting"  = "grey"
)

# Read matrices into a long table
read_one <- function(file, comp_name) {
  read.delim(file, header = TRUE, row.names = 1, check.names = FALSE) %>%
    rownames_to_column("gene") %>%
    mutate(comparison = comp_name)
}

dat <- map2_dfr(files, names(files), read_one)

# Genes present in at least two comparisons
genes_keep <- dat %>%
  dplyr::distinct(gene, comparison) %>%
  dplyr::count(gene, name = "n_comparisons") %>%
  dplyr::filter(n_comparisons >= 2) %>%
  dplyr::pull(gene)

dat2 <- dat %>% dplyr::filter(gene %in% genes_keep)

# Build combined Expression / Targeting matrix
to_mat <- function(value_col) {
  dat2 %>%
    dplyr::select(gene, comparison, !!sym(value_col)) %>%
    pivot_wider(names_from = comparison, values_from = !!sym(value_col)) %>%
    column_to_rownames("gene") %>%
    as.matrix()
}

expr_mat <- to_mat("Expression")
target_mat <- to_mat("Targeting")

combined_mat <- cbind(
  expr_mat[, names(files), drop = FALSE],
  target_mat[, names(files), drop = FALSE]
)
colnames(combined_mat) <- c(
  paste0("Expr_", names(files)),
  paste0("Target_", names(files))
)
combined_mat[!is.finite(combined_mat)] <- NA

# Curated annotations (semicolon-separated)
ann <- read.csv(
  file.path(out_dir, "genes_present_in_at_least_two_comparisons_EditedAnnotation.csv"),
  sep = ";",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ann2 <- ann %>%
  dplyr::filter(
    gene %in% rownames(combined_mat),
    !is.na(`gene annotation`),
    `gene annotation` != ""
  )

# Preserve curated annotation order
combined_mat <- combined_mat[ann2$gene, , drop = FALSE]
row_labels <- paste(ann2$gene, ann2$`gene annotation`, sep = " | ")

# Keep genes with consistent logFC direction in Expression AND Targeting
is_consistent_direction <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(FALSE)
  all(x > 0) || all(x < 0)
}

expr_cols <- grep("^Expr_", colnames(combined_mat), value = TRUE)
target_cols <- grep("^Target_", colnames(combined_mat), value = TRUE)

keep <- apply(combined_mat[, expr_cols, drop = FALSE], 1, is_consistent_direction) &
  apply(combined_mat[, target_cols, drop = FALSE], 1, is_consistent_direction)

combined_mat <- combined_mat[keep, , drop = FALSE]
ann2 <- ann2[keep, , drop = FALSE]
row_labels <- row_labels[keep]

# Column annotation
ha_col <- HeatmapAnnotation(
  Type = rep(c("Expression", "Targeting"), each = length(files)),
  Comparison = rep(names(files), 2),
  col = list(
    Type = type_cols,
    Comparison = comparison_cols
  )
)

# Color scale centered on 0
lim <- max(abs(combined_mat), na.rm = TRUE)
col_fun <- colorRamp2(
  c(-lim, -3, -2, -1, 0, 1, 2, 3, lim),
  c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
    "white",
    "#FCAE91", "#FB6A4A", "#CB181D", "#67000D")
)

ht <- Heatmap(
  combined_mat,
  name = "logFC",
  col = col_fun,
  na_col = "grey90",
  top_annotation = ha_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_labels = row_labels,
  row_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  column_title = "Expression and sRNA targeting logFC",
  row_title = "Genes with consistent logFC direction"
)

pdf(file.path(out_dir, "expression_targeting_heatmap_consistent.pdf"), width = 10, height = 12)
draw(ht)
dev.off()

# Note that the final manuscript figure was manually edited to:
# 1) Reduce the set to genes showing consistent, albeit not significant, changes in direction expression
# 2) Reorder genes to mirror groups of same DE + DT direction change

write.csv(
  cbind(gene = ann2$gene, annotation = ann2$`gene annotation`, combined_mat),
  file.path(out_dir, "expression_targeting_heatmap_consistent_matrix.csv"),
  row.names = FALSE
)
