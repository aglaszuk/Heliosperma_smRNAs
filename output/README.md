# Output directory

Results from the analysis pipeline. Subfolders mirror analysis steps:

```
output/
├── DT_CG/                 # Common garden
│   ├── 1_PCA/
│   ├── 2_trimming/        # trimmed counts gitignored (regenerated)
│   ├── 3_libSizeNormalization/
│   ├── 4_DT_edgeR/        # DTR lists (main DT results)
│   ├── 5_Histograms_DTRs_DTGs_DEGs/
│   ├── 6_overlapDTRs_DTGs_DEGs/
│   ├── 7_Heatmaps_DT_DEGs/
│   └── 8_GOenrichment/
├── DT_RT/                 # Reciprocal transplant (same step layout)
└── Heatmap_Figure3/       # Combined expression + targeting heatmap
```

Caches (`.drake/`) and large trimmed count tables are excluded via the root `.gitignore`.
