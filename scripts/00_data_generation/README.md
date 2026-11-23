# Data generation

These scripts make the tables in `data/intermediate/` that the figure scripts read. Most of them take
hours or need data that is not in this project (BAMs, the MATOU catalogue, the 18S swarm tables). Their
outputs for the paper are already in `data/intermediate/`, so **none of them has to be run to make the
figures**. They are here to show how every number was obtained, and to rerun a step if needed.

Run them from the project root, in the order of their numbers. Each script starts with a header that
gives its inputs, outputs, run time and the original script it comes from. Shared functions are in
`R/gappa.R` (gappa wrappers, placement factorization, clade sums) and `R/ratios.R` (metaT/metaG ratios).

When a step writes jplaces (metaT/metaG or metaT/HK ratios), it writes them to `data/intermediate/`. The
next step uses those if they exist, and otherwise the copy used for the paper, symlinked in
`data/raw/placement/`. Nothing is ever written through the symlinks.

The exception is `00_upstream_pipelines/`, which is a record of how the data in `data/raw/` were made
(see below). Those scripts point to the `data/raw/` symlinks, as the originals did, so **running them would
overwrite the data of the original projects**. None was run here.

## Upstream pipelines (`00_upstream_pipelines/`)
How the quantification tables, the reference trees and the per-sample jplaces in `data/raw/` were made.
They take days, need nextflow / singularity / databases outside the project, and some steps were manual
(tree curation in iTOL), so they document the procedure rather than reproduce it byte by byte.

| Step | What it does | Output (in `data/raw/`) | Original |
|---|---|---|---|
| `01_quantification_nextflow.sh` | `beaplab/transcriptome_metaT_quantification`: metaT reads vs every transcriptome (BWA-mem2 + Salmon) | `quantification/{alignment,mapping}` | DF `src/02_quantification.sh` |
| `02_eggnog_annotation.sh` | eggNOG-mapper of the transcriptome CDS | `genomic/functional_annotation/` | DF `src/01_eggnog-annotation.sh` |
| `03a/b/c_reference_tree_<gene>.sh` | Flavodoxin, SIT, NRT2 trees: Pfam hmmsearch over EukProt v01 → curation → MAFFT/trimAl → IQ-TREE; copied as phyloplace references | `gene_data/reference_trees/<gene>/`, `placement/phyloplace_references/nonhk_genes/` | GE `generate_nonHK_gene_trees/obtain_<gene>_phylogeny.sh` |
| `04_reference_taxonomies.sh` | lineage of every reference sequence, from its genus | `placement/hk_gene_taxonomy/{SIT,NRT2,Flavodoxin}.taxonomy.txt` | GE `helper_create-nonhk-taxonomy.sh` |
| `05_matou_translations.sh` | six-frame translation of the MATOU unigenes (≥ 10 aa) | `matou_alltranslations_filt10aa.faa` (MATOU folder) | MATOU `creation_faa.sh` |
| `06_phyloplace_inputs.R` | HMMs of the non-HK alignments; phyloplace sample sheets (318 HK + 4 non-HK genes) | `placement/phyloplace_samplesheets/` | PT `helper-create-phyloplace-samplesheet_v2.R`, `creation_hmms_newgenes.sh` |
| `07_phyloplace_nextflow.sh` | `nf-core/phyloplace`: hmmsearch in MATOU → hmmalign → EPA-ng | `placement/epang/` | PT `phyloplace_search_v2.sh` |
| `08_matou_subset_tables.R` | MATOU metaT/metaG/raw-count tables restricted to the placed unigenes | `placement/matou_tables/` | PT `helper_subsetting_matou_tables.R` |
| `09_jplace_multiplicities.R` | one jplace per gene and sample, MATOU values as multiplicities | `placement/epang_multiplicities_{metaT,metaT_rawcounts,metaG}/` | PT `create_multiplicity_files.R` |
| `10_extract_metazoa.sh` | metazoan placements removed from the HK jplaces | `placement/hk_matou_jplaces/` | GE `main.nf`, process `extract_metazoans` |

`helpers/`: `eukprot_hmm_hits.R` (hmmsearch hits → fasta, replaces the three `helper_parse-<gene>-phylogeny.R`),
`align_trim_fasttree.py` (MAFFT → trimAl → FastTree), `reference_taxonomy.R`.
DF = delving_into_florenciella, PT = phylotranscriptomics, GE = gene_explorer.

Not reproducible exactly: the 20 metazoan sequences of the Flavodoxin and NRT2 trees were sampled without a
seed, the Flavodoxin and SIT selections were made by hand in iTOL (the lists are in the tree folders), and
IQ-TREE was seeded with the date. The second quantification run used a folder of hand-downloaded EukProt v3
transcriptomes that no longer exists. The exact `--no-premask` call for the 14 genes in
`hk_genes_nopremask.csv` was not recorded.

Earlier inputs, not covered here: the 318 HK reference alignments, trees, HMMs and taxonomies (Richter et
al., `external$hk_taxonomies_richter`; the taxonomy files used are identical to those), the selection and
download of genomes/transcriptomes (DF `src/00_genome-down-and-selection/`), the FAD1 tree (placed but not in
any figure), and the concatenated species tree of GE.

## Steps

### Environmental data
| Step | Output | Run time |
|---|---|---|
| `00_env_table.R` | `data/raw/env/tara_env.tsv.gz` | seconds |

The single table of Tara sample metadata (iron, nutrients, sample materials, MATOU codes and run names; read
with `read_tara_env()`, `R/io_env.R`). It is made from the originals in `external$env_root` and replaces the six
tables that were copied before:
- `envdata_Tara_All_downloaded_metaT_nisaba.tsv`
- `envdata_tara_polar_protist.tsv`
- `matou_envdata_nonpolar.tsv`
- `flavodoxin_samples_nonpolar_metadata.tsv`
- `Tara_Oceans_Pangaea_context.rds`
- `context_general.tsv`

They held the same values. The header of the script says where each column comes from.

### Fig 2: *Pelagomonas* iron DE and Flavodoxin metaT/metaG ratios
| Step | Output (`data/intermediate/fig2/`) | Run time |
|---|---|---|
| `fig2/01_pelagomonas_fastemu_iron.R` | `pelagomonas_fastemu_iron_coeffs.tsv` (copy of `fastemu_reference_comparison/all_hk/`; fits `.rds` there) | ~13 h |
| `fig2/01b`, `01c`, `01d` | reference-set comparison (paper 73-HK selection vs all 281 HK) and convergence scan, in `fastemu_reference_comparison/` | seconds; scan ~1 h |
| `fig2/02_flavodoxin_sample_iron.R` | `flavodoxin_samples_iron.tsv` (iron of every sample, polar ones added) | seconds |
| `fig2/03_flavodoxin_metaT_metaG_ratios.R` | `flavodoxin_metaT_metaG_pergene_counts.tsv.gz`, `flavodoxin_metagratio_jplace/` | ~10 min |
| `fig2/04_flavodoxin_ratio_factorization.R` | `flavodoxin_ratio_factorization/` | ~30 min |
| `fig2/05_flavodoxin_graff2023_clades.R` | `flavodoxin_graff2023_clade_ratios.tsv` | ~10 min |

### Fig 3: genetic context of *Florenciella parvula*
| Step | Output (`data/intermediate/fig3/`) | Run time |
|---|---|---|
| `fig3/01_aai_sourmash.sh` | `compar_transcriptomes.csv` | minutes |
| `fig3/02_blast_florenciella.sh` | `blast_Florenciella_allgenomes.tsv.gz` | ~1 h |
| `fig3/03_crossmapping.R` | `crossmapping_florenciella_against_all_transcriptomes.tsv.gz`, `crossmapping/` | hours (BAMs) |
| `fig3/04_expressed_genes_per_sample.R` | `allstatistics_florenciella_plus_other.tsv.gz` | ~1 h |
| `fig3/05_coverage_examples.R` | `coverage_examples.tsv` | ~1.5 min (BAMs) |

### Fig 5: tree and placement statistics
| Step | Output (`data/intermediate/fig5/`) | Run time |
|---|---|---|
| `fig5/01_tree_and_placement_stats.R` | `gene_stats_summary.tsv`, `per_sample_placements.tsv.gz`, `lwr_distribution_recalc.tsv.gz` | hours |

### Fig 6: taxonomic composition of the housekeeping genes (also feeds Fig 7)
| Step | Output | Run time |
|---|---|---|
| `fig6/01_gappa_assign_hk.sh` | `data/intermediate/hk_taxassign/` | hours |
| `fig6/02_combine_assign_summaries.sh` | `data/intermediate/hk_taxassign_combined/` (paper copy: `data/raw/placement/hk_taxassign_combined`) | minutes |
| `fig6/03_housekeeping_gene_filter.R` (needs fig5/01) | `fig6/housekeeping_gene_filter_matou.tsv` | ~6 min |
| `fig6/04_taxonomy_compositions.R` | `fig6/taxon_compositions.parquet` (raw caches in `fig6/cache/`) | ~2 min with the caches, hours without |

### Fig 7: SIT and NRT2 metaT/HK ratio factorization
| Step | Output (`data/intermediate/fig7/`) | Run time |
|---|---|---|
| `fig7/01_samples_env.R` (needs fig2/02) | `{SIT,NRT2}_samples_env.tsv`, `SIT_metadata_ironcat.tsv`, `NRT2_metadata_nitratecat.tsv` | seconds |
| `fig7/02_metaT_HK_ratios.R <SIT\|NRT2>` (needs fig6/02, fig6/03) | `<gene>_hk_normalisation/`, `<gene>_hkratio_jplace/` | ~10 min per gene |
| `fig7/03_hkratio_factorization.R <SIT\|NRT2>` | `SIT_iron_hkratio_factorization/`, `NRT2_nitrate_hkratio_factorization/` | ~30 min per gene |

`fig7/02` also takes `Flavodoxin` (Supp S7, written to `data/intermediate/supplementary/S7/`), and `fig7/03`
takes `metagratio` as a second argument (Supp S9, written to `data/intermediate/supplementary/S9/`).

Fig 4 has no data-generation step: it starts from the *Pelagomonas* TPM table in `data/raw/`.

## Supplementary (`supplementary/`)
Outputs go to `data/intermediate/supplementary/<item>/`. S1, S2, S8 and T1 have no step of their own: S1, S2
and S8 start from `data/raw/quantification/`, and T1 from fig5/01, fig6/03 and fig6/04.

| Step | Output | Run time |
|---|---|---|
| `S3/01_blast_bathycoccus_vs_matou.sh` | `data/raw/quantification/statistics/bathycoccus_prasinos_blast_92id.tsv` | ~1 h |
| `S3/02_bathycoccus_matou_tables.R` | `S3/bathycoccus_prasinos_matou_{metaT,metaG}.tsv` (paper copies in `data/raw/quantification/statistics/`) | ~1 h (streams the MATOU tables) |
| `S4/01_read_vs_seq_placement_caches.R [gene ...]` | `S4/aggregated/<gene>.rds` (7 genes) | hours |
| `S5/01_unigene_metaT_metaG_lineage.R [gene ...]` | `S5/<gene>.tsv.gz` | ~10 min per gene |
| `S6/01_metag_taxonomic_bias.R [ncores] [out]` (needs fig6/03; step 3 also fig2/03, S7/01, S9/01, fig7/02) | `S6/hk_lineage_shares.tsv.gz`, `hk_detection.tsv`, `<gene>_detection.tsv` | ~6 min with 16 cores, ~60 GB RAM |
| `S7/01_flavodoxin_hk_ratios.sh` (= `fig7/02 Flavodoxin`) | `S7/Flavodoxin_hk_normalisation/`, `S7/Flavodoxin_hkratio_jplace/` | ~2 min with 8 cores |
| `S7/02_graff2023_clade_sums.R` (needs fig2/02, S7/01) | `S7/graff2023_{metaT,metaG,hkratio}_counts.tsv`, `S7/graff2023_hkratio_pergene.tsv.gz` | ~30 min |
| `S9/01_metaT_metaG_ratios.R <SIT\|NRT2>` | `S9/<gene>_metaT_metaG_pergene_counts.tsv.gz`, `S9/<gene>_metagratio_jplace/` | ~10 min per gene |
| `fig7/03_hkratio_factorization.R <SIT\|NRT2> metagratio` (needs fig7/01, S9/01) | `S9/SIT_iron_metagratio_factorization/`, `S9/NRT2_nitrate_metagratio_factorization/` | ~30 min per gene |

Caches copied from GE and DF (not regenerated here): `S4/aggregated/`, `S5/*.tsv.gz`, the `S6` tables, the
`S7` clade tables and the per-unigene HK-ratio multiplicities (GE `graff2023_hkratio_multiplicities/`, 2,668
files, merged into `graff2023_hkratio_pergene.tsv.gz`, 10 MB, same values), and from each S9 factorization folder the
files the figure reads (the clade jplaces, 250 MB, were left out). Not copied:
- the per-gene caches of S4 (`per_name/`, `krd/`: 1.2 GB)
- the S6 `hk_assign/` cache (124 MB; step 0 remakes it in ~5 min)

New symlinks in `data/raw/placement/`:
- `{SIT,NRT2}_metagratios_jplace` → GE `phylofactorization_<gene>/ratios_jplace_scaled`
- `read_placements_metazoa_extracted` → GE `data/placement/placements_metazoa_extracted` (read placements of Richter et al., 147 GB)

## Dependencies between figures
```
upstream/01 ──► fig2/01, fig3/03-05, Fig 4          upstream/02 ──► Fig 2-4 (annotations)
upstream/03a-c ──► 04 ──┐
upstream/05 ────────────┼──► 06 ──► 07 ──► 08 ──► 09 ──► 10 ──► fig5/01, fig6/01
                        │                         └──► fig2/02-03 (Flavodoxin), fig7/02 (SIT, NRT2)
                        └──► fig7/02 (taxonomies)
fig5/01 ──► fig6/03 ──┐
fig6/01 ──► fig6/02 ──┼──► fig7/02 ──► fig7/03
fig2/02 ──► fig7/01 ──┘
fig6/02 ──► fig6/04
fig2/02 ──► fig2/03 ──► fig2/04, fig2/05
00_env_table ──► fig2/01, fig2/02, fig6/04, fig7/01

Supplementary:
fig2/02 ──► fig7/02 Flavodoxin (= S7/01) ──► S7/02
S9/01 ──► fig7/03 metagratio
fig6/03 ──► S6/01 ◄── (step 3) fig2/03, S7/01, S9/01, fig7/02
fig5/01, fig6/03, fig6/04 ──► T1
upstream/09 ──► S5/01, S9/01, S4/01 (with the read placements of Richter et al.)
```

## Checked against the paper outputs (2026-11-19)
The cheap steps, and small subsets of the long ones, were rerun and compared with the cached tables:
| Step | Check | Result |
|---|---|---|
| fig2/02 | full run | identical, except `155DCM1KKQQ13filter`, which the original also listed although it has no ratio jplace (duplicated metaT sample) |
| fig2/03 + 05 (helpers) | 3 samples (one polar) | per-gene counts, ratio jplace totals and clade sums identical |
| fig3/04 | full run | not completed (stopped after hours); not verified |
| fig3/05 | full run | the figure built from it is pixel-identical to the endpoint |
| fig5/01 | gene COQ4_mito | identical |
| fig6/01 + 02 | gene COQ4_mito (581 samples) | identical |
| fig6/03 | full run | identical |
| fig6/04 | full run from the raw caches | identical |
| fig7/01 | full run | all four tables identical |
| fig7/03 (helpers) | clades of both genes; clade sums of 2 SIT samples | identical |
| env table (2026-11-20) | `00_env_table.R`, then fig2/01 (covariates only), fig2/02, fig7/01 and the fig6/04 sample map + 18S compositions, against the six old tables | fig2/01: same 218 samples and iron categories. fig2/02: same `Sample Code`, PANGAEA id, iron and `ironcat`; `event_time` corrected (it was misparsed) and 5 `event_method` values now as in PANGAEA. fig7/01: all four tables byte-identical. fig6/04: sample map and V4/V9 compositions identical |
| upstream helpers/eukprot_hmm_hits.R | SIT hits table | same 950 sequences (names and sequences) as `SIT_eukprotv01_aa.fasta` |
| upstream 03a-c | trimmed alignments vs phyloplace MSAs | same sequences; names differ only by `\|` → `_` (now done in the scripts) |
| upstream 04 | all three genes | NRT2, Flavodoxin identical; SIT differs in 12 names that still had `\|` in 2024 |
| upstream 06 | sample sheets (no hmmbuild) | 318 HK + 14 no-premask + 4 non-HK rows, every path exists |
| upstream 09 | SIT, 2 samples (one polar), metaT | multiplicities per query identical |
| upstream 10 | ABCE, 1 sample | 5381 queries, multiplicities identical |

Supplementary steps (2026-11-20):
| Step | Check | Result |
|---|---|---|
| S3/02 | step 1 (unigene selection) | same 5284 unigenes as DF's `bathycoccus_prasinos_geneids.tsv`; but the paper's metaT/metaG tables predate that blast table: 43 of their 5106 unigenes are not selected now (MATOU extraction not rerun) |
| S4/01 | ABCE: sample matching (396 samples); extraction, assign and EDPL of 2 samples | identical to the cached rds |
| S5/01 | SIT and Flavodoxin | identical |
| S6/01 | steps 0-2, full run (into the scratchpad); step 3 for Flavodoxin | identical |
| fig7/02 (extended) | full runs of SIT and Flavodoxin | `pergene_hkratios` identical (318,114 and 981,919 rows); Flavodoxin taxonomy fix identical |
| S7/02 | 2 samples of metaT, metaG and HK ratio | clade sums and per-unigene multiplicities identical; the merged `graff2023_hkratio_pergene.tsv.gz` also identical for 2 samples (2026-11-20) |
| S9/01 | 2 SIT samples | per-gene counts and ratio jplaces identical |
| fig7/03 metagratio (helpers) | clades of both genes; clade sums of 2 SIT samples | identical |

fig2/01 (2026-11-23): Fig 2 now uses the all-HK reference set (281 genes) instead of the paper's 73-gene
selection; the all-HK run of 2026-11-20 is reused, and the flavodoxin CAH0373297, which no longer passes the
pre-filter, was tested on its own afterwards (`all_hk/run_extra.log`). The paper's table is kept in
`fastemu_reference_comparison/paper/`.

Not rerun: fig2/04, fig3/01-03, the factorization part of fig7/03 (also for S9), S3/01,
and the upstream pipelines 01, 02, 05, 07, 08 and the tree inference of 03a-c.
