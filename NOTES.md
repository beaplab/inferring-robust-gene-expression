# Working notes (formerly README.md): provenance, decisions and checks

Code and data for the paper figures. It merges two earlier projects:
- `gene_explorer` (phylogenetic placement of housekeeping genes, Flavodoxin/SIT/NRT2): `/scratch/data1/aauladell/projects_big/gene_explorer`
- `delving_into_florenciella` (transcriptome quantification of single species, *Pelagomonas* and *Florenciella*): `/scratch/data3/aauladell/project_big/delving_into_florenciella`

It also uses the EPA-ng placements of `phylotranscriptomics` (`/scratch/data3/aauladell/project_big/phylotranscriptomics`).

## Layout
```
config/paths.R            all paths (project-relative + external databases)
R/                        reusable functions; R/load.R loads libraries, paths and every R/*.R
data/raw/                 inputs: copied tables + symlinks to the large jplace folders
data/intermediate/        outputs of the long data-generation steps (copied, not rerun)
scripts/00_data_generation/  scripts that produced data/intermediate (kept for reference, long runtimes)
scripts/figN_*.R          one script per main figure
scripts/supplementary/    one script per supplementary figure / table (S1-S9, T1), with the captions in its README
scripts/utils/            compare_to_endpoint.R: pixel diff of a regenerated figure vs its endpoint PDF
results/endpoint_figures/ reference PDFs of the submitted figures
results/figures/          regenerated figures (supplementary/ for S1-S9)
results/tables/supplementary/, results/supplementary_information/   T1, T2, SI2 and the supplementary stats tables
```
Each area has a `supplementary/` subfolder for the supplementary material (`data/intermediate/supplementary/`,
`scripts/00_data_generation/supplementary/`, `results/endpoint_figures/supplementary/`).

Run any figure from the project root: `/usr/bin/Rscript scripts/figN_*.R`.
Use the system R 4.1 (`/usr/bin/Rscript`). It has tidyverse, ggtree, treeio, arrow and fastEmu. The default `Rscript` on the PATH (miniconda) has no tidyverse, and `/usr/local/bin/Rscript` (R 4.5) has no ggtree.

## Figures
| Figure | Script | Original script | Compared with the endpoint PDF |
|---|---|---|---|
| 2 | `scripts/fig2_species_vs_gene_DE.R` | GE `scripts/figures/figure1_flavodoxin.R` | 1.1% of pixels differ, only the jittered points |
| 3 | `scripts/fig3_florenciella_context.R` | DF `src/figures/combined_florenciella_figure.R` | identical |
| 4 | `scripts/fig4_hk_vs_random_predpower.R` | DF `src/eda/comparison_housekeeping_torandomgenes_pred_metaGab_carradec.R` | identical |
| 5 | `scripts/fig5_phylogenetic_tree_stats.R` | GE `scripts/figures/phylogenetic_tree_examples_matou.R` | 0.05% of pixels differ (anti-aliasing) |
| 6 | `scripts/fig6_omic_perspectives.R` | GE `scripts/figures/final_figure_hk_taxonomy_fractions.R` | identical |
| 7 | `scripts/fig7_phylofactorization_SIT_NRT2.R` | GE `scripts/figures/summary_phylofactorization.R` | differs by design (see below); against the current GE output only the jittered points differ |

The original Figs 2 and 7 jittered the points without a seed, so they cannot be reproduced point by point. The jitter now has a fixed seed, so from here on every run gives the same figure.

The scripts that made `data/intermediate/` are in `scripts/00_data_generation/`, one folder per figure. See its README for the order, the run times, and which steps were rerun and checked. Its `00_upstream_pipelines/` folder records how the data in `data/raw/` were made (metaT quantification with nextflow, reference trees, nf-core/phyloplace on MATOU, per-sample jplaces); these were not run here and would write into the original projects.

GE = gene_explorer, DF = delving_into_florenciella.

Fig 7 is reproduced with the current NRT2 nitrate data, which include the polar samples. Its NRT2 panel therefore differs from `results/endpoint_figures/figure_7_*.pdf`, which was made before the polar samples were added (gene_explorer `results/backup_nonpolar_nitrate_2026-11-18`).

## Supplementary material
Supplementary Figs S1–S9 and Tables T1–T2 follow the gene_explorer proposal
(`manuscript/supplementary_proposal/`). SI2 is a copy of a PDF. Every figure is identical to the proposal's,
apart from these differences:
- S9 differs only in the jittered points, as Fig 7 does.
- S3 and S5 are identical only when drawn with R 4.5. R 4.1 gives the same numbers with small rendering
  differences.

T1 matches the proposal's table. See `scripts/supplementary/README.md` for the scripts, the checks, the
captions and the open points.

## Data
`data/raw/quantification/`: tables of the nextflow pipeline `beaplab/transcriptome_metaT_quantification` (DF `src/02_quantification.sh`): `alignment/` (BWA-mem2) and `mapping/` (Salmon), with TPM, Numreads and gene-characteristics per transcriptome and dataset. BAMs are not included.

`data/raw/placement/`: symlinks (not copies) to:
| link | target |
|---|---|
| `hk_matou_jplaces` | GE `data/placement/matou/placements_metazoa_extracted` (HK genes × 581 MATOU metaT samples, 15 GB) |
| `hk_gene_taxonomy` | GE `data/placement/housekeeping_gene_tax` |
| `epang`, `epang_multiplicities_*` | phylotranscriptomics `data/matou_search/phyloplace/...` (Flavodoxin/SIT/NRT2) |
| `{SIT,NRT2}_jplace_raw`, `*_hkratios_jplace`, `*_metagratios_jplace` | GE `data/stats/phylofactorization_*` |
| `read_placements_metazoa_extracted` | GE `data/placement/placements_metazoa_extracted` (metaT reads placed directly, Richter et al.; Supp S4) |
| `hk_taxassign_combined` | GE `data/stats/matou/taxassign_combined` (gappa assign of the HK genes, one table per gene, 2 GB) |

`data/raw/env/tara_env.tsv.gz`: the single table of Tara sample metadata (PANGAEA variables, nutrients, MATOU sample codes, downloaded run names), made by `scripts/00_data_generation/00_env_table.R` and read with `read_tara_env()` (`R/io_env.R`).

`data/raw/quantification/statistics/sourmash_gather_output` is a symlink to DF `data/quantification/statistics/sourmash_gather_output` (sourmash gather per transcriptome × sample, 4,386 files; no script here reads it).

`data/raw/gene_data/hk_reference_alignments` is a symlink to GE `data/gene_data/27_no_sequence_ids_26` (HK reference alignments, 1.8 GB).

`data/intermediate/fig6/cache/` holds the raw caches of `fig6/04_taxonomy_compositions.R` (copied, ~75 MB), so that step reruns in minutes. `matou_taxonomy_lineage_ids.parquet` (416 MB) is a symlink to GE.
