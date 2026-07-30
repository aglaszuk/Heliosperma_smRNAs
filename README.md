# Heliosperma smRNA differential targeting

Analysis scripts and results for **differential targeting (DT)** of genomic regions by small RNAs (smRNAs) in *Heliosperma pusillum* alpine and montane ecotypes.

Two experimental designs are analysed in parallel:

| Design | Abbreviation | Description |
|--------|--------------|-------------|
| Common garden | **CG** | Plants grown in a common environment |
| Reciprocal transplant | **RT** | Plants grown at alpine vs. montane sites |

Analyses include edgeR-based DT calling, overlaps with differentially expressed genes (DEGs), heatmaps, and GO enrichment.

## Repository structure

```
Heliosperma_smRNAs/
├── Heliosperma_smRNAs.Rproj   # Open this in RStudio
├── README.md
├── data/                      # Inputs (see data/README.md)
│   ├── DEGlists_CG_SzukalaEtAl2022/
│   ├── DEGlists_RT_SzukalaEtAl2023/
│   ├── GO/
│   ├── gene2region.uniq.txt
│   ├── smRNA_counts_CG.txt.gz
│   └── smRNA_counts_RT.txt.gz
├── scripts/                   # Analysis pipeline (run in numbered order)
└── output/
    ├── DT_CG/                 # Common-garden results
    ├── DT_RT/                 # Reciprocal-transplant results
    └── Heatmap_Figure3/       # Combined Figure 3 heatmap
```

## Requirements

- R (≥ 4.0 recommended)
- Bioconductor packages: `RUVSeq`, `edgeR`, `HTSFilter`, `topGO`, `ComplexHeatmap`
- CRAN packages: `dplyr`, `tidyverse`, `ggplot2`, `ggpubr`, `drake`, `GOplot`, `circlize`, …

Install everything once:

```bash
Rscript scripts/install_packages.R
```

## Setup

1. Clone this repository.
2. Open `Heliosperma_smRNAs.Rproj` in RStudio (sets the working directory to the project root).

Scripts auto-detect the project root; you can also run them with `Rscript` from the repository root. Count matrices are stored as `.gz` and read directly by the pipeline (see [`data/README.md`](data/README.md)).


## Pipeline (run in order)

### Common garden (CG)

| Step | Script | Output folder |
|------|--------|---------------|
| 1 | `scripts/1_DetectDT_CG.R` | `output/DT_CG/1_PCA` … `4_DT_edgeR` |
| 2 | `scripts/2_plot_histograms_and_overlaps_CG.r` | `output/DT_CG/5_Histograms_*`, `6_overlap*` |
| 3 | `scripts/3_plot_heatmaps_CG.R` | `output/DT_CG/7_Heatmaps_DT_DEGs` |
| 4 | `scripts/4_GO_enrichment_CG.R` | `output/DT_CG/8_GOenrichment` |

### Reciprocal transplant (RT)

| Step | Script | Output folder |
|------|--------|---------------|
| 1 | `scripts/1_DetectDT_RT.R` | `output/DT_RT/1_PCA` … `4_DT_edgeR` |
| 2 | `scripts/2_plot_histograms_and_overlaps_RT.r` | `output/DT_RT/5_Histograms_*`, `6_overlap*` |
| 3 | `scripts/3_plot_heatmaps_RT.R` | `output/DT_RT/7_Heatmaps_DT_DEGs` |
| 4 | `scripts/4_GO_enrichment_RT.R` | `output/DT_RT/8_GOenrichment` |

### Combined figure

| Step | Script | Output folder |
|------|--------|---------------|
| 5 | `scripts/5_plotHeatmap_GC_RT_Figure3.R` | `output/Heatmap_Figure3` |

Helper scripts (sourced by others, do not run alone):

- `scripts/1_DE_functions.R`
- `scripts/4_GO_functions.R`
- `scripts/00_set_project_root.R`

## Notes on file sizes

Count matrices are included as gzip archives (~30 MB each), which stay under GitHub’s file-size limits. Regenerable intermediates (trimmed counts, `.drake/` caches) are gitignored.

Analysis **results** under `output/` (gene lists, PDFs, GO tables) are intended to be versioned.

## Related publications

DEG lists derive from:

- **CG:** Szukala et al. (2022) — common-garden transcriptome
- **RT:** Szukala et al. (2023) — reciprocal-transplant transcriptome

Update citation details here when the smRNA manuscript is published.

## License

Add a license of your choice (e.g. MIT for code). Third-party annotations retain their original terms.
