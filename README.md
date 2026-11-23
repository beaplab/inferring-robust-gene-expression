# Inferring robust gene expression from metatranscriptomes

Code and data behind the figures and supplementary material of the paper. Every figure can be redrawn
in seconds to minutes from the tables shipped here, and every table can be traced back, script by script,
to the raw sequencing data.

## What the project does

Metatranscriptomes (metaT) measure how much a gene is expressed in a sample, but a change in metaT can
come from a change in expression or from a change in the abundance of the organism that carries the
gene. The usual fix is to divide by metagenomic (metaG) abundance, which is often missing or biased. This
project tests **housekeeping (HK) genes** as the denominator instead: genes expressed at a stable level,
whose metaT then tracks the abundance of their organism.

The analysis works at two levels.

**1. Single species, from reads mapped to a transcriptome** (*Pelagomonas calceolata*, *Florenciella*,
*Bathycoccus prasinos*)
- Tara Oceans metaT reads are mapped to each transcriptome (BWA-mem2 and Salmon, nextflow pipeline
  `beaplab/transcriptome_metaT_quantification`).
- **Fig 4, S1–S3:** subsets of HK genes predict the metaG abundance of the species better than random
  gene subsets do. **S8:** BWA and Salmon counts compared.
- **Fig 3:** for *Florenciella*, its relatives (SMAGs, other strains) show which genes can be quantified
  unambiguously by mapping.
- **Fig 2 (left):** differential expression of *P. calceolata* genes between low- and high-iron samples
  (fastEmu, with all HK genes as reference).

**2. Whole community, from phylogenetic placement** (MATOU v1.5 unigene catalogue, Tara Oceans)
- The MATOU unigenes of 318 candidate HK genes and of three target genes (Flavodoxin, SIT, NRT2) are
  placed on reference trees (nf-core/phyloplace: hmmsearch, hmmalign, EPA-ng). Each sample becomes one
  jplace, with the unigene metaT (or metaG) values as multiplicities (gappa).
- **Fig 5, T1:** how good the reference trees and placements of the HK genes are.
- **Fig 6, S6:** what taxonomic composition the HK genes see, against 18S amplicons and the whole MATOU
  metaG/metaT catalogues.
- **Figs 2 (right), 7, S5, S7, S9:** metaT of the target genes normalised either by metaG (metaT/metaG
  ratio) or by the HK genes of the same lineage (metaT/HK ratio). Placement factorization (gappa) then
  finds the clades of the tree whose expression changes with iron or nitrate.
- **S4:** placing reads directly vs placing assembled unigenes.

## Quick start

```bash
# from the project root
/usr/bin/Rscript scripts/fig4_hk_vs_random_predpower.R          # one main figure
/usr/bin/Rscript scripts/supplementary/S6_metag_taxonomic_bias.R  # one supplementary figure
```

- Run every script **from the project root**. Each script starts with `source('R/load.R')`, which loads
  the libraries, the paths (`config/paths.R`) and all the shared functions in `R/`.
- Use the R and packages listed in [`PACKAGES.md`](PACKAGES.md). The figures were made with R 4.1.2
  (`/usr/bin/Rscript` on the original server). S3 and S5 are pixel-identical to the paper only with R 4.5.
- Outputs go to `results/figures/` (pdf + png) and `results/tables/`.
- To compare a new figure with the paper version:
  `Rscript scripts/utils/compare_to_endpoint.R results/endpoint_figures/<paper>.pdf results/figures/<new>.pdf`.
  It counts the pixels that differ and writes a `_diff.png`.

## Layout

```
README.md                  this file
PACKAGES.md                R, R packages and command-line tools, with versions
NOTES.md                   working notes: provenance, decisions and checks made while building the project
config/paths.R             every path, in one place (p('key', ...) builds a path)
R/                         shared functions, loaded by R/load.R
scripts/
  figN_*.R                 one script per main figure (2-7)
  supplementary/           one script per supplementary figure / table (S1-S9, T1); README with captions
  00_data_generation/      the slow steps that made data/intermediate/ (README with order and run times)
    00_env_table.R         the Tara sample metadata table
    00_upstream_pipelines/ how data/raw/ was made (quantification, trees, placement); documentation only
    fig2/ ... fig7/, supplementary/S3 ... S9/
  utils/compare_to_endpoint.R
data/
  raw/                     inputs (tables copied from the original projects; large jplace folders symlinked)
  intermediate/            outputs of 00_data_generation/, one folder per figure (+ supplementary/)
results/
  figures/, tables/        what the scripts write
  endpoint_figures/        the figures as submitted (references)
  supplementary_information/
```

## Main figures

| Figure | Script | Shows |
|---|---|---|
| 2 | `scripts/fig2_species_vs_gene_DE.R` | iron-related differential expression in *P. calceolata* (fastEmu) and of the Flavodoxin clades on the reference tree (metaT/metaG ratio) |
| 3 | `scripts/fig3_florenciella_context.R` | genetic context of *Florenciella parvula* (AAI, cross-mapping) and which of its genes can be quantified |
| 4 | `scripts/fig4_hk_vs_random_predpower.R` | HK gene subsets vs random gene subsets as predictors of *P. calceolata* abundance |
| 5 | `scripts/fig5_phylogenetic_tree_stats.R` | reference trees and placement quality of the HK genes and the three target genes |
| 6 | `scripts/fig6_omic_perspectives.R` | taxonomic composition seen by the HK genes vs 18S-V4/V9 and MATOU metaG/metaT |
| 7 | `scripts/fig7_phylofactorization_SIT_NRT2.R` | placement factorization of SIT (iron) and NRT2 (nitrate) on metaT/HK ratios |

Supplementary figures S1–S9 and tables T1–T2 are described, with their captions, in
[`scripts/supplementary/README.md`](scripts/supplementary/README.md). SI2 (the fastEmu model) is in
`results/supplementary_information/`.

Differences from the submitted figures:
- **Jitter.** Figs 2 and 7 and S9 jittered their points without a seed. The scripts now set one, so the
  points move slightly against the paper PDFs but stay fixed from run to run.
- **Fig 7.** It uses the NRT2 nitrate data with the polar samples, added after the submitted version.

## How the data were made

Three levels, each documented by the scripts that made it:

1. **`data/raw/`: inputs.**
   - metaT quantification tables (`quantification/`)
   - functional annotations (`genomic/`)
   - gene lists and reference trees (`gene_data/`)
   - per-sample jplaces (`placement/`)
   - Tara sample metadata (`env/tara_env.tsv.gz`)

   How they were made is in `scripts/00_data_generation/00_upstream_pipelines/`. Those steps take days,
   need nextflow, singularity and external databases, and included manual curation (trees in iTOL), so
   they document the procedure rather than being meant to be rerun.
2. **`data/intermediate/`: the slow analysis steps.**
   - per-sample ratios and gappa placement factorizations
   - fastEmu fits
   - taxonomic assignments and composition caches

   Their scripts are in `scripts/00_data_generation/<figure>/`, numbered in run order. Their outputs are
   shipped, so they don't need to be rerun to make the figures. The README there gives the dependencies,
   the run times and which steps were rerun and checked against the shipped outputs.
3. **`results/`: the figures and tables**, made by `scripts/figN_*.R` and `scripts/supplementary/`.

Main external resources: the Tara Oceans metaT/metaG data, the MATOU v1.5 unigene catalogue (Carradec et
al. 2018), EukProt v1, the 318 candidate HK genes of Strassert et al. (2021) with the reference
alignments of Richter et al., the Tara PANGAEA environmental context, and the Flavodoxin clades of Graff
et al. (2023).

## Moving the project

The GitHub repository holds only the code and documentation: `data/` and `results/` are not
versioned (`data/intermediate/` alone is >100 GB) and will be deposited separately.

All relative paths are in `config/paths.R` (`paths`). The absolute ones, needed only by
`00_data_generation/`, are in `external`. Large placement folders in `data/raw/placement/` and a few
others are symlinks into the original projects on the server where this was built; `NOTES.md` lists
them. To use the project elsewhere, copy the targets of the symlinks you need, or just the figure scripts
with `data/intermediate/`, which is enough to redraw every figure.

## Software

R 4.1.2 with tidyverse 2.0.0, ggtree 3.9.1, arrow 23.0.1, fastEmu 1.0.0 / radEmu 2.3.2, and gappa 0.9.0
for the placement steps. The full list, with versions and with the tools of the upstream pipelines, is in
[`PACKAGES.md`](PACKAGES.md).
