# Supplementary figures and tables

One script per item, run from the project root with `/usr/bin/Rscript scripts/supplementary/<script>`.
Outputs go to `results/figures/supplementary/` and `results/tables/supplementary/`. The references (the
proposal's figures, gene_explorer `manuscript/supplementary_proposal/`) are in
`results/endpoint_figures/supplementary/`. Compare a figure against its reference with
`Rscript scripts/utils/compare_to_endpoint.R <reference.pdf> <new.pdf>`.

The slow steps behind S3–S7 and S9 are in `scripts/00_data_generation/supplementary/` (see the README there).
Their outputs are cached in `data/intermediate/supplementary/`, so every script here runs in seconds to minutes.

GE = gene_explorer, DF = delving_into_florenciella.

| Item | Script | Original | Compared with the reference |
|---|---|---|---|
| S1, S2 | `S1_S2_hk_subsets_ksweep_functional_mix.R` | DF `src/eda/comparison_housekeeping_functionalmix_and_ksweep.R` | both figures identical; stats tables identical to DF's |
| S3 | `S3_hk_vs_random_bathycoccus.R` | DF `src/eda/comparison_housekeeping_allgenes_bathy_matou.R` | identical with R 4.5; 0.26% of pixels differ with R 4.1 (same numbers) |
| S4 | `S4_read_vs_seq_placement.R` | GE `scripts/figures/visualize_read_vs_seq_placement.R` (section 1) | identical; both tables identical to GE's |
| S5 | `S5_metag_coverage_by_lineage.R` | GE `scripts/figures/metag_coverage_by_lineage.R` | identical with R 4.5; 1.2% of pixels differ with R 4.1 (legend only); tsv identical |
| S6 | `S6_metag_taxonomic_bias.R` | GE `scripts/figures/metag_taxonomic_bias.R` | identical |
| S7 | `S7_flavodoxin_hk_vs_metag.R` | GE `scripts/figures/flavodoxin_hk_vs_metag.R` | identical |
| S8a, S8b | `S8_bwa_vs_salmon.R` | DF `src/eda/helper_comparison-btw-salmon-and-bwamem.R` (`compare_samples()`) | identical (PNG) |
| S9 | `S9_SIT_NRT2_metaGratio_factorization.R` | GE `scripts/figures/summary_phylofactorization.R` | 1.5% of pixels differ, only the jittered points (as in Fig 7) |
| T1 | `T1_candidate_hk_genes.R` | new; the GE table `housekeeping_gene_summary_matou.tsv` had no script | same rows, columns and values as the proposal's T1 (rel. diff < 1e-14) |
| T2 | none: `results/tables/supplementary/T2_datasets_aligned_transcriptomes.xlsx` (copy) | | |
| SI2 | none: `results/supplementary_information/SI2_fastEmu_statistical_model.pdf` (copy) | | |

Notes:
- **R version.** S3 and S5 were drawn with R 4.5 / ggplot2 4.0.2 (`/usr/local/bin/Rscript`), which
  reproduces them pixel for pixel. The project default, R 4.1 (`/usr/bin/Rscript`), gives the same numbers
  with small rendering differences. The other scripts need R 4.1, which has ggtree.
- **S1 and S2** are one script because they share the random stream (seed 4223, S2 draws first).
- **S8.** The original picked the sample panels with `sample()` without a seed. The script fixes the
  samples to those of the paper figure.
- **S9 and Fig 7** share `manuscript_row()` / `manuscript_theme()` (`R/phylofactorization.R`). The S9 NRT2
  factorization has no polar samples (it was made before they were added), unlike Fig 7.
- **S4.** The sequence placements now hold more genes in common with the read placements than the 7 of the
  paper. The data-generation step defaults to those 7.

## Captions
From gene_explorer `manuscript/supplementary_proposal/captions.md`, with the open points noted there.

**Supplementary Figure S1.** Housekeeping vs random gene subsets across subset sizes in *Pelagomonas
calceolata*.
- (a) R² distributions between metagenomic abundance and the geometric mean of expression for 1,000 subsets
  per size k (30–100 genes), drawn from the 254 housekeeping genes (grey) or from 15,371 expressed genes
  (mean TPM 1–2000; blue). The red line marks all 254 housekeeping genes.
- (b) Median R² (± IQR/2) per k: housekeeping subsets improve with k and approach the full set, while random
  subsets stay at ~0.10.
- (c) Fraction of genes shared between two housekeeping subsets (5,000 pairs per k), with the expectation
  k/254. Part of the narrowing of the housekeeping distributions with k comes from this overlap.

*Open:* at k = 45 the median HK R² is ~0.25, against ~0.28 at k = 95, and the distribution is wider. Check
that this matches the text's "reduced overall predictive power and a more uniform distribution".

**Supplementary Figure S2.** Functional composition of housekeeping gene subsets and their predictive power
in *P. calceolata*.
- (a) Distribution of R² between metagenomic abundance and the geometric mean of expression for 5,000 random
  subsets of 95 housekeeping genes; dashed lines mark the bottom and top 10%.
- (b) Functional-category composition of the genes drawn in the low, mid and high R² subsets, compared with
  the expected composition of all 254 housekeeping genes.

*Open:* the text says "no specific gene functional enrichments". The figure shows a small but significant
difference (p = 0.001; category explains 22% of the variance).

**Supplementary Figure S3.** Housekeeping genes as an abundance proxy in *Bathycoccus prasinos*.
- (A) Adjusted R² of the linear regression of the metagenomic abundance of *B. prasinos* (median over its
  genes, per sample) on the geometric mean of the metatranscriptomic expression of 20 genes. There are 1,000
  random draws from all moderately expressed genes (between the 25th and 90th percentiles of mean
  expression, housekeeping genes excluded; blue) and 1,000 from the 65 housekeeping genes (black). Red line:
  all 65 housekeeping genes.
- (B) Example regressions: one draw from each distribution (top) and all 65 housekeeping genes (bottom).
  One point per sample; grey band: 95% confidence interval.

Housekeeping sets predict abundance better than random sets (median adjusted R² 0.49 vs 0.42; all 65 genes
0.50), but the gap is smaller than in *P. calceolata*.

**Supplementary Figure S4.** Comparison of placing metatranscriptomic reads directly with placing MATOU
assembled unigenes, for seven genes.
- (A) Total multiplicity per sample for each approach, by size fraction.
- (B) log2 ratio of sequence-based to read-based multiplicity, by size fraction and depth.

**Supplementary Figure S5.** Share of the metatranscriptomic signal of Flavodoxin, SIT and NRT2 with a
metagenomic counterpart, by lineage.
- Only metaT samples with a paired metaG sample (same station, depth and size fraction) are used: 445, 439
  and 443 of the 581 metaT samples, holding 76%, 74% and 79% of the Flavodoxin, SIT and NRT2 metaT.
- Unigenes are assigned to the lineage of the reference sequences below their best placement branch.
- Bars give the share of each lineage's normalised metaT carried by unigenes with metaG reads in the paired
  sample (blue) or without them (grey; these get a metaT/metaG ratio of 0).
- Diamonds give the same share but for the unigene × sample pairs.
- Numbers on the right give the share of each gene's metaT held by the lineage.

**Supplementary Figure S6.** The MATOU metagenome under-represents dinoflagellates.
- (A) log2 ratio between the share of each lineage among housekeeping (HK) gene metaG reads and its share
  among HK metaT reads (median and IQR over 448 paired Tara samples; Metazoa excluded).
- (B) Fraction of HK unigenes expressed in the metaT that have any metaG signal, by metaT expression
  quintile.

*Fill-in for the text:* "Xfold" ≈ ~7-fold (median log2 = −2.7 for Dinoflagellata, panel A).

**Supplementary Figure S7.** Housekeeping-gene normalisation reproduces metagenome normalisation for
flavodoxin.
- (A) log2 fold change (low/high iron) of each flavodoxin clade defined by Graff et al. (2023): clade 1
  (diel cycle), clade 2 (iron limiting), green algae + dinoflagellates, and the rest of the tree. It is
  estimated from metaT, metaG, the metaT/metaG ratio and two versions of the metaT/HK ratio. Filled points:
  Wilcoxon p.adj < 0.05 (Holm-adjusted per data type).
- (B) Per-sample agreement between the summed metaT/HK ratio and the summed metaT/metaG ratio for each
  clade. Only the genes and samples with metaG are used, so that only the denominator differs. Coloured by
  iron category; line: linear fit on log axes; samples with a summed metaT/metaG ratio of 0 are not shown.

*Open:* for "Rest of the tree", the metaT/metaG ratio is not significant while the HK ratio is.

**Supplementary Figure S8.** Gene counts from BWA-mem2 (filtered at ≥95% identity) vs Salmon for the
transcriptomes of (a) *Pelagomonas calceolata* and (b) *Florenciella*, in a subset of Tara
metatranscriptomes.
- One panel per sample; red: 1:1 line; blue: smooth fit.
- Counts are strongly correlated, with BWA giving consistently lower counts than Salmon.
- In a few samples (mostly 20–180 µm) only a handful of genes have high Salmon counts.

**Supplementary Figure S9.** Placement factorization of SIT and NRT2 on metaT/metaG ratios, for comparison
with Fig. 7 (metaT/HK ratios).
- (A) SIT over iron (low/high cutoff 0.5 nmol/L, Darwin model, 5 m). (B) NRT2 over nitrate + nitrite
  (cutoff 1 µmol/L).
- Left: reference tree, branches coloured by taxonomy. The strip next to the tips shows the clades found by
  placement factorization of the per-unigene metaT/metaG ratios: raw metaT reads over the metaG abundance of
  the same unigene in the paired metaG sample, or 0 when the unigene has no metaG.
- Right: summed metaT/metaG ratio per clade and sample (log axis) in low and high samples, with the main
  taxonomic groups of each clade below its name. Basal: placements outside every clade.
- Labels: Wilcoxon rank-sum p-value (unadjusted) and log2 fold change of the medians (high/low; NaN when a
  median is 0).

*Open:* not cited in the text yet.

**Supplementary Table T1.** Reference trees, placement quality and taxonomic composition of the 318
candidate housekeeping genes (Strassert et al. 2021) after placement of the MATOU unigenes. The column
definitions are in the proposal's `captions.md` and in the header of `T1_candidate_hk_genes.R`.

**Supplementary Table T2.** Public marine metatranscriptomic datasets to which the transcriptomes were
aligned.

**Supplementary Information SI2.** Statistical model used to estimate differential expression with fastEmu
for a single species. *Open:* section 2.3 refers to a model that is not in the file, and the text needs
editing into SI style.
