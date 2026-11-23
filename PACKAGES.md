# Software versions

Versions installed on the server where the project was built and checked (2026-11-20). The figures and
the analysis steps were run with R 4.1.2; R 4.5.0 was used only where noted.

## R

| | Version | Used for |
|---|---|---|
| R | 4.1.2 (2021-11-01) | every script (`/usr/bin/Rscript`) |
| R | 4.5.0 (2025-04-11) | optional: S3 and S5 are pixel-identical to the proposal figures only with it (`/usr/local/bin/Rscript`; it has no ggtree) |

## R packages

Version under R 4.1.2 (the one used), and under R 4.5.0 where it differs.

| Package | R 4.1.2 | R 4.5.0 | Used by |
|---|---|---|---|
| tidyverse | 2.0.0 | | every script (loaded by `R/load.R`) |
| dplyr | 1.2.1 | | |
| tidyr | 1.3.2 | | |
| readr | 2.1.4 | 2.2.0 | |
| purrr | 1.2.2 | 1.2.1 | |
| stringr | 1.6.0 | | |
| forcats | 1.0.0 | 1.0.1 | |
| tibble | 3.2.1 | 3.3.1 | |
| ggplot2 | 4.0.3 | 4.0.2 | every figure |
| scales | 1.4.0 | | figures |
| patchwork | 1.3.2 | | figure layouts (most figures) |
| cowplot | 1.2.0 | | Figs 2, 3 |
| ggtext | 0.2.0 | | Fig 3, tree panels (`R/phylofactorization.R`) |
| ggrepel | 0.9.8 | | Fig 5 |
| ggridges | 0.5.7 | | Fig 3, S1/S2 |
| ggmagnify | 0.4.2 | not installed | Fig 3 |
| viridis | 0.6.5 | | Fig 3 |
| RColorBrewer | 1.1.3 | | S1/S2 |
| png | 0.1.9 | | `scripts/utils/compare_to_endpoint.R` |
| ggtree | 3.9.1 | not installed | tree figures: Figs 2, 5, 7, S9 (`R/phylofactorization.R`) |
| treeio | 1.18.1 | 1.34.0 | dependency of ggtree |
| tidytree | 0.4.7 | | tree tables (`R/phylofactorization.R`) |
| ape | 5.8.1 | | tree reading and statistics: Fig 5, S5/01, `R/phylofactorization.R`, `R/io_placement.R` |
| arrow | 23.0.1.2 | | parquet caches and MATOU tables: Fig 6, T1, fig6/04, S3/02, upstream 08–09 |
| data.table | 1.14.4 | 1.18.2.1 | large tables: fig3/05, fig7/02, S6/01 |
| readxl | 1.4.3 | 1.4.5 | xlsx inputs: Fig 4, S1/S2, `00_env_table.R` |
| jsonlite | 1.8.7 | 2.0.0 | jplace reading: `R/io_placement.R`, fig6/03, S4/01 |
| janitor | 2.2.1 | | Fig 4, S1/S2 |
| matrixStats | 1.5.0 | | fig6/03 |
| argparser | 0.7.3 | | upstream helpers (`eukprot_hmm_hits.R`, `reference_taxonomy.R`) |
| Biostrings | 2.62.0 | 2.78.0 | upstream helper `eukprot_hmm_hits.R` |
| Rsamtools | 2.10.0 | 2.26.0 | fig3/03 (BAMs, not in the project) |
| radEmu | 2.3.2.0 | | fastEmu dependency |
| fastEmu | 1.0.0 | | Fig 2 differential expression (`fig2/01`) |
| parallel, grid | base R | | |

## Command-line tools

| Tool | Version | Used for |
|---|---|---|
| gappa | 0.9.0 | placement steps of `00_data_generation/` (assign, extract, multiplicity, factorization, krd); `--basal-clade-name basal_branches` reproduces the old extract names |
| gappa (in nf-core/phyloplace) | 0.8.0 | graft and heat-tree outputs of the phyloplace runs |
| nextflow | 25.10.0 | upstream pipelines |
| nf-core/phyloplace | v2.0.0 (g827b85d) | MATOU placement (`00_upstream_pipelines/07`) |
| EPA-ng (in phyloplace) | 0.3.8 | placement |
| HMMER | 3.4 | hmmsearch / hmmalign / hmmbuild (phyloplace and reference trees) |
| easel (in phyloplace) | 0.48 | alignment formatting |
| seqtk (in phyloplace) | 1.4-r122 | subsets |
| beaplab/transcriptome_metaT_quantification | branch `main` | metaT quantification (BWA-mem2 + Salmon, run in its containers; `00_upstream_pipelines/01`) |
| singularity | 3.11.0 | containers of the nextflow pipelines |
| eggNOG-mapper | 2.1.12 | functional annotation (`00_upstream_pipelines/02`, `03a`, `03c`) |
| MAFFT | 7.490 | reference tree alignments |
| trimAl | 1.4.rev22 | alignment trimming |
| IQ-TREE | 2.0.7 (and 1.6.12) | reference trees |
| FastTree | 2.1.11 | trees during curation (`03a`, `03c`, `helpers/align_trim_fasttree.py`) |
| seqkit | 2.10.1 | fasta handling (reference trees, MATOU translations) |
| BLAST+ | 2.14.1 | blastn: fig3/02, S3/01; blastp: NRT2 tree (`03c`) |
| samtools | 1.13 | fig3/05 (BAM coverage) |
| Salmon (local) | 1.4.0 | local install; the quantification used the pipeline's container |
| sourmash | 4.6.1 | AAI of Fig 3, fig3/01 (conda environment `smash`) |
| Python | 3.8.19 | `helpers/align_trim_fasttree.py` |
| pdftoppm (poppler) | 22.02.0 | `compare_to_endpoint.R` |

Versions inside containers (BWA-mem2 and Salmon of the quantification pipeline) are fixed by the
pipeline revision and are not listed separately.
