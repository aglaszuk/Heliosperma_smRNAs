# Title: "Plot histograms of DTRS, DTGs, DEGs and genes DE and DT by small RNAs in a common garden (CG) experimental setup"
# Author: Aglaia Szukala
# Short description: Script needs lists of DTRs and DEGs for each ecotype comparison and a table mapping genes to genomic regions.

# Set working directory
# Set project root (portable; works after cloning the repository)
.root_helper <- c("scripts/00_set_project_root.R", "../scripts/00_set_project_root.R")
source(.root_helper[file.exists(.root_helper)][1])
rm(.root_helper)

#Load Libraries
library(UpSetR)
library(grid)
library(ggplot2)
library(patchwork)

# Define output directory
out_dir <- "output/DT_CG/5_Histograms_DTRs_DTGs_DEGs/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Define ecotype pairs
pairs <- c("1", "3", "4", "5")

# Define pairs colors
specieColors <- c(
  "1" = "#A71E34",
  "3" = "#00798C",
  "4" = "#F59700",
  "5" = "#6E5DAC"
)

#############
# LOAD DATA #
#############
deg_list <- lapply(pairs, function(p) {
  read.table(file = paste0("data/DEGlists_CG_SzukalaEtAl2022/AvsM_", p, "_fdr0.05.txt"))
})
names(deg_list) <- pairs

dtr_up_list <- lapply(pairs, function(p) {
  read.table(file = paste0("output/DT_CG/4_DT_edgeR/MvsA", p, "_DTRs_qlTest_poslogFC.txt"))
})
names(dtr_up_list) <- pairs

dtr_dw_list <- lapply(pairs, function(p) {
  read.table(file = paste0("output/DT_CG/4_DT_edgeR/MvsA", p, "_DTRs_qlTest_neglogFC.txt"))
})
names(dtr_dw_list) <- pairs

reg2gene <- read.table(file = "data/gene2region.uniq.txt")
colnames(reg2gene) <- c("gene", "region")

# Table of trimmed counts (universe for overlap permutation tests)
trim_counts <- read.table("output/DT_CG/2_trimming/smRNA_counts_CG_trimmed.txt")
N_genes <- length(unique(reg2gene$gene[reg2gene$region %in% rownames(trim_counts)]))
N_regions <- nrow(trim_counts)
N_genes
N_regions

# Rule: for each gene, keep DTR with lowest FDR (matches heatmaps)
collapse_dtr_to_gene <- function(dtr_up_tab, dtr_dw_tab, reg2gene_tab) {
  dtr_all <- rbind(dtr_up_tab, dtr_dw_tab)
  dtr_all$region <- rownames(dtr_all)

  merged <- merge(reg2gene_tab, dtr_all, by = "region")
  merged <- merged[order(merged$gene, merged$FDR), ]

  gene_level <- merged[!duplicated(merged$gene), ]
  rownames(gene_level) <- gene_level$gene

  gene_level
}

######################################################
# Build a table summarizing amounts in each category #
######################################################

# DTG direction uses one DTR per gene (lowest FDR), matching heatmaps
get_pair_sets <- function(deg_tab, dtr_up_tab, dtr_dw_tab, reg2gene_tab) {
  deg_all <- rownames(deg_tab)
  deg_pos <- rownames(deg_tab)[deg_tab$logFC > 0]
  deg_neg <- rownames(deg_tab)[deg_tab$logFC < 0]
  
  dtr_up_regions <- rownames(dtr_up_tab)
  dtr_dw_regions <- rownames(dtr_dw_tab)
  dtr_all_regions <- unique(c(dtr_up_regions, dtr_dw_regions))
  
  dtg_tab <- collapse_dtr_to_gene(dtr_up_tab, dtr_dw_tab, reg2gene_tab)
  dtg_all <- rownames(dtg_tab)
  dtg_up <- rownames(dtg_tab)[dtg_tab$logFC > 0]
  dtg_dw <- rownames(dtg_tab)[dtg_tab$logFC < 0]
  
  list(
    total_DEG = deg_all,
    total_DTR = dtr_all_regions,
    total_DTG = dtg_all,
    total_genesDE_DT = intersect(deg_all, dtg_all),
    DEG_Mdw_DTR_Mup = intersect(deg_neg, dtg_up),
    DEG_Mup_DTR_Mdw = intersect(deg_pos, dtg_dw),
    DEG_Mdw_DTR_Mdw = intersect(deg_neg, dtg_dw),
    DEG_Mup_DTR_Mup = intersect(deg_pos, dtg_up),
    DTRs_logFC_Mup = dtr_up_regions,
    DTRs_logFC_Mdw = dtr_dw_regions,
    DEGs_logFC_Mup = deg_pos,
    DEGs_logFC_Mdw = deg_neg
  )
}

pair_sets <- lapply(pairs, function(p) {
  get_pair_sets(
    deg_tab = deg_list[[p]],
    dtr_up_tab = dtr_up_list[[p]],
    dtr_dw_tab = dtr_dw_list[[p]],
    reg2gene_tab = reg2gene
  )
})
names(pair_sets) <- pairs

# Check one
pair_sets[1][1]

# Make summary table
summary_df <- do.call(rbind, lapply(pairs, function(p) {
  data.frame(
    pair = p,
    total_DTR = length(pair_sets[[p]]$total_DTR),
    DTRs_logFC_Mup = length(pair_sets[[p]]$DTRs_logFC_Mup),
    DTRs_logFC_Mdw = length(pair_sets[[p]]$DTRs_logFC_Mdw),
    total_DTG = length(pair_sets[[p]]$total_DTG),
    total_DEG = length(pair_sets[[p]]$total_DEG),
    total_genesDE_DT = length(pair_sets[[p]]$total_genesDE_DT),
    DEGs_logFC_Mup = length(pair_sets[[p]]$DEGs_logFC_Mup),
    DEGs_logFC_Mdw = length(pair_sets[[p]]$DEGs_logFC_Mdw),
    DEG_Mup_DTR_Mdw = length(pair_sets[[p]]$DEG_Mup_DTR_Mdw),
    DEG_Mdw_DTR_Mup = length(pair_sets[[p]]$DEG_Mdw_DTR_Mup),
    DEG_Mup_DTR_Mup = length(pair_sets[[p]]$DEG_Mup_DTR_Mup),
    DEG_Mdw_DTR_Mdw = length(pair_sets[[p]]$DEG_Mdw_DTR_Mdw)
  )
}))

rownames(summary_df) <- NULL
print(summary_df)

write.table(
  summary_df,
  file = file.path(out_dir, "Summary_table_DEG_DTR_DTG.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

####################################################################
# Plot a barplot of DTRs, DTGs, DEGs and genes DE and DT by pair #
####################################################################

plot_df <- data.frame(
  pair = rep(summary_df$pair, 4),
  metric = rep(
    c("total_DTR", "total_DTG", "total_DEG", "total_genesDE_DT"), #"DEG_Mdw_DTR_Mup", "DEG_Mup_DTR_Mdw"),
    each = nrow(summary_df)
  ),
  value = c(
    summary_df$total_DTR,
    summary_df$total_DTG,
    summary_df$total_DEG,
    summary_df$total_genesDE_DT
  )
)

plot_df$pair <- factor(plot_df$pair, levels = pairs)

summary_plot <- ggplot(plot_df, aes(x = pair, y = value, fill = pair)) +
  geom_col(width = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = specieColors) +
  labs(x = "Ecotype pair", y = "Count") +
  theme_classic() +
  theme(legend.position = "none")
summary_plot

ggsave(file.path(out_dir, "Summary_barplots.pdf"), summary_plot, width = 10, height = 7)
ggsave(file.path(out_dir, "Summary_barplots.png"), summary_plot, width = 10, height = 7, dpi = 300)

#######################################################################################
# Run UpSet such that bar counts items present in ALL listed contrasts, regardless of #
# membership in other contrasts (|A∩B| >= |A∩B∩C| always holds).                      #
#######################################################################################

# Define output directory
out_dir <- "output/DT_CG/6_overlapDTRs_DTGs_DEGs/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Get all combinations of possible overlaps
get_all_combinations <- function(set_names = pairs) {
  unlist(lapply(seq_along(set_names), function(k) {
    apply(combn(set_names, k), 2, paste, collapse = "_")
  }))
}

all_combinations <- get_all_combinations()

# Define a function to separate the pairs IDs in sets of all_combinations
parse_combo_sets <- function(combo_name) {
  strsplit(combo_name, "_", fixed = TRUE)[[1]]
}

# Build a dataframe of membership of items (regions or genes) to a certain combination
build_membership_df <- function(set_list, set_names = pairs) {
  all_items <- sort(unique(unlist(set_list)))
  if (length(all_items) == 0) return(NULL)
  
  membership_cols <- setNames(
    lapply(set_names, function(nm) all_items %in% set_list[[nm]]),
    set_names
  )
  
  data.frame(item = all_items, membership_cols, check.names = FALSE)
}

# Count amount of items in intersection (inclusive, i.e. including all items also those shared with another combination)
count_inclusive_intersection <- function(membership_df, combo_sets) {
  if (length(combo_sets) == 0) return(0L)
  mat <- membership_df[, combo_sets, drop = FALSE]
  sum(apply(mat, 1, all))
}

# Retrieve regions or genes
get_inclusive_items <- function(membership_df, combo_sets) {
  if (length(combo_sets) == 0) return(character(0))
  mat <- membership_df[, combo_sets, drop = FALSE]
  membership_df$item[apply(mat, 1, all)]
}

####################################################
# Test significance of overlaps using permutations #
####################################################

# Overlap significance test parameters
# N_genes / N_regions already computed from trimmed counts above
nperm <- 100
alpha <- 0.05
set.seed(123)

# Define permutation test function
perm_test_inclusive_intersection <- function(set_sizes, combo_sets, universe_size, nperm = nperm) {
  replicate(nperm, {
    sampled_sets <- lapply(set_sizes, function(s) sample.int(universe_size, s, replace = FALSE))
    names(sampled_sets) <- names(set_sizes)
    length(Reduce(intersect, sampled_sets[combo_sets]))
  })
}

format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# Make UpSet plots for each category (DTRs, DTGs, DEGs, genes DE and DT)
plot_inclusive_upset <- function(membership_df,
                                 set_names,
                                 perm_results,
                                 pair_colors,
                                 category_name,
                                 alpha = 0.05) {
  plot_data <- perm_results[order(-perm_results$observed, perm_results$combination), , drop = FALSE]
  plot_data$id <- seq_len(nrow(plot_data))
  plot_data$id_f <- factor(plot_data$id, levels = plot_data$id)
  plot_data$bar_fill <- ifelse(plot_data$significant & plot_data$observed > 0, "#C1121F", "#404040")

  set_sizes_df <- data.frame(
    set = factor(set_names, levels = rev(set_names)),
    size = vapply(set_names, function(s) sum(membership_df[[s]]), numeric(1)),
    stringsAsFactors = FALSE
  )

  matrix_dots <- do.call(rbind, lapply(seq_len(nrow(plot_data)), function(i) {
    combo_sets <- parse_combo_sets(plot_data$combination[i])
    data.frame(
      id = i,
      set = factor(combo_sets, levels = set_names),
      stringsAsFactors = FALSE
    )
  }))

  matrix_lines <- do.call(rbind, lapply(seq_len(nrow(plot_data)), function(i) {
    combo_sets <- parse_combo_sets(plot_data$combination[i])
    y_pos <- match(combo_sets, set_names)
    if (length(y_pos) < 2) return(NULL)
    data.frame(id = i, ymin = min(y_pos), ymax = max(y_pos))
  }))

  p_sets <- ggplot(set_sizes_df, aes(x = size, y = set, fill = set)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    scale_fill_manual(values = pair_colors[set_names]) +
    labs(x = paste("Size of", category_name), y = NULL) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08)))

  p_matrix <- ggplot()
  if (!is.null(matrix_lines) && nrow(matrix_lines) > 0) {
    p_matrix <- p_matrix + geom_segment(
      data = matrix_lines,
      aes(x = id, xend = id, y = ymin, yend = ymax),
      linewidth = 0.9,
      color = "#404040"
    )
  }
  p_matrix <- p_matrix +
    geom_point(
      data = matrix_dots,
      aes(x = id, y = as.numeric(set)),
      size = 3.2,
      color = "#404040"
    ) +
    scale_y_continuous(
      breaks = seq_along(set_names),
      labels = set_names,
      limits = c(0.5, length(set_names) + 0.5)
    ) +
    scale_x_continuous(breaks = plot_data$id, labels = NULL) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )

  main_label <- if (any(plot_data$significant & plot_data$observed > 0)) {
    paste0("Inclusive overlap - ", category_name, "  (red = perm. p < ", alpha, ")")
  } else {
    paste0("Inclusive overlap - ", category_name)
  }

  p_main <- ggplot(plot_data, aes(x = id_f, y = observed, fill = bar_fill)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = observed), vjust = -0.3, size = 3) +
    scale_fill_identity() +
    labs(x = NULL, y = main_label) +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 11))

  (p_sets + p_main + plot_layout(widths = c(1, 2.8))) /
    p_matrix +
    plot_layout(heights = c(1.2, 2.5), widths = c(1, 2.8))
}

plot_upset_category <- function(set_list,
                                category_name,
                                pair_colors,
                                universe_size,
                                set_names = pairs,
                                nperm = nperm,
                                alpha = 0.05,
                                export_files = TRUE) {

  membership_df <- build_membership_df(set_list, set_names = set_names)

  if (is.null(membership_df) || nrow(membership_df) == 0) {
    message("No items found for category: ", category_name)
    return(invisible(NULL))
  }

  set_sizes <- sapply(set_list[set_names], length)

  perm_results <- lapply(all_combinations, function(combo_name) {
    combo_sets <- parse_combo_sets(combo_name)
    obs <- count_inclusive_intersection(membership_df, combo_sets)

    perm_counts <- perm_test_inclusive_intersection(
      set_sizes = set_sizes,
      combo_sets = combo_sets,
      universe_size = universe_size,
      nperm = nperm
    )

    p_val <- (sum(perm_counts >= obs) + 1) / (length(perm_counts) + 1)

    data.frame(
      combination = combo_name,
      observed = obs,
      expected = mean(perm_counts),
      p_value = p_val,
      p_label = format_p(p_val),
      significant = p_val < alpha,
      stringsAsFactors = FALSE
    )
  })

  perm_results <- do.call(rbind, perm_results)
  observed_counts <- setNames(perm_results$observed, perm_results$combination)

  if (export_files) {
    for (combo_name in all_combinations) {
      combo_sets <- parse_combo_sets(combo_name)
      items <- get_inclusive_items(membership_df, combo_sets)
      write.table(
        items,
        file = file.path(out_dir, paste0(category_name, "_", combo_name, ".txt")),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
      )
    }
  }

  write.table(
    perm_results,
    file = file.path(out_dir, paste0(category_name, "_permutation_results.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p <- plot_inclusive_upset(
    membership_df = membership_df,
    set_names = set_names,
    perm_results = perm_results,
    pair_colors = pair_colors,
    category_name = category_name,
    alpha = alpha
  )

  invisible(list(
    plot = p,
    membership_df = membership_df,
    observed_counts = observed_counts,
    permutation_results = perm_results
  ))
}

############################
# Extract UpSet categories #
############################
sets_total_DTR <- lapply(pair_sets, function(x) x$total_DTR)
sets_total_DTG <- lapply(pair_sets, function(x) x$total_DTG)
sets_total_DEG <- lapply(pair_sets, function(x) x$total_DEG)
sets_genes_DEDT <- lapply(pair_sets, function(x) x$total_genesDE_DT)
sets_DEG_Mdw_DTR_Mup <- lapply(pair_sets, function(x) x$DEG_Mdw_DTR_Mup)
sets_DEG_Mup_DTR_Mdw <- lapply(pair_sets, function(x) x$DEG_Mup_DTR_Mdw)

#######################################
# Plot and save inclusive UpSet plots #
#######################################
res_total_DEG <- plot_upset_category(
  sets_total_DEG, "total_DEG", specieColors,
  universe_size = N_genes, nperm = nperm
)
pdf(file.path(out_dir, "UpSet_total_DEG.pdf"), width = 11, height = 8, useDingbats = FALSE)
print(res_total_DEG$plot)
dev.off()

res_total_DTR <- plot_upset_category(
  sets_total_DTR, "total_DTR", specieColors,
  universe_size = N_regions, nperm = nperm
)
pdf(file.path(out_dir, "UpSet_total_DTR.pdf"), width = 11, height = 8, useDingbats = FALSE)
print(res_total_DTR$plot)
dev.off()

res_total_DTG <- plot_upset_category(
  sets_total_DTG, "total_DTG", specieColors,
  universe_size = N_genes, nperm = nperm
)
pdf(file.path(out_dir, "UpSet_total_DTG.pdf"), width = 11, height = 8, useDingbats = FALSE)
print(res_total_DTG$plot)
dev.off()

res_genes_DEDT <- plot_upset_category(
  sets_genes_DEDT, "genes_DEDT", specieColors,
  universe_size = N_genes, nperm = nperm
)
pdf(file.path(out_dir, "UpSet_genes_DEDT.pdf"), width = 11, height = 8, useDingbats = FALSE)
print(res_genes_DEDT$plot)
dev.off()

res_DEG_Mdw_DTR_Mup <- plot_upset_category(
  sets_DEG_Mdw_DTR_Mup, "DEG_Mdw_DTR_Mup", specieColors,
  universe_size = N_genes, nperm = nperm
)
pdf(file.path(out_dir, "UpSet_total_DEG_Mdw_DTR_Mup.pdf"), width = 11, height = 8, useDingbats = FALSE)
print(res_DEG_Mdw_DTR_Mup$plot)
dev.off()

res_DEG_Mup_DTR_Mdw <- plot_upset_category(
  sets_DEG_Mup_DTR_Mdw, "DEG_Mup_DTR_Mdw", specieColors,
  universe_size = N_genes, nperm = nperm
)
pdf(file.path(out_dir, "UpSet_total_DEG_Mup_DTR_Mdw.pdf"), width = 11, height = 8, useDingbats = FALSE)
print(res_DEG_Mup_DTR_Mdw$plot)
dev.off()