# Data directory

Inputs for the Heliosperma smRNA differential-targeting pipeline.

## Contents

| Path | Description |
|------|-------------|
| `smRNA_counts_CG.txt.gz` | smRNA count matrix (common garden); gzip-compressed |
| `smRNA_counts_RT.txt.gz` | smRNA count matrix (reciprocal transplant); gzip-compressed |
| `DEGlists_CG_SzukalaEtAl2022/` | DEG lists and related files from the common-garden study |
| `DEGlists_RT_SzukalaEtAl2023/` | DEG lists and logFC splits from the reciprocal-transplant study |
| `GO/geneID2GO_Athaliana_genes.txt` | Gene → GO term mapping (topGO) |
| `GO/geneID2description.txt` | Gene descriptions |
| `GO/helio.genome.v2.Athaliana_blast2go.gff` | Blast2GO annotation GFF |
| `gene2region.uniq.txt` | Mapping between genes and smRNA target regions |

## Notes

- Count matrices are region × sample tables of smRNA read counts. Scripts read the `.gz` files directly via `gzfile()`.
- Trimmed counts written to `output/*/2_trimming/` are gitignored; they are recreated by scripts `1_DetectDT_*.R`.
- Nested `.drake/` caches under `output/` are regenerable and ignored.
