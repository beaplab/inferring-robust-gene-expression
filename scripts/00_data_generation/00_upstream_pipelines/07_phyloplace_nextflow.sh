#!/usr/bin/env bash
# Upstream 07: phylogenetic placement of the MATOU unigenes on the reference trees of the 318 HK genes and
# of SIT, NRT2 and Flavodoxin (+ FAD1), with the nextflow pipeline nf-core/phyloplace.
#
# Per gene (row of the sample sheet): hmmsearch of its HMM against the six-frame MATOU translations
# (05) -> hits aligned to the reference alignment (hmmalign) -> placed on the reference tree with
# EPA-ng -> gappa. Every unigene is placed once (no abundances yet; those are added in 09).
#
# Input:  data/raw/placement/phyloplace_samplesheets/*.csv (06)
#         matou_alltranslations_filt10aa.faa (05, external$matou_translations)
# Output: data/raw/placement/{epang,hmmer,...}; this project uses
#         data/raw/placement/epang/<gene>.epa_result.jplace.gz (symlink to phylotranscriptomics
#         data/matou_search/phyloplace/epang)
# Run time: days. Needs nextflow and singularity.
# Original: phylotranscriptomics scripts/matou_phyloplace/phyloplace_search_v2.sh
#
# NOT RUN IN THIS PROJECT. If run, it writes into phylotranscriptomics through the epang symlink.
#
# Usage (from the project root): bash scripts/00_data_generation/00_upstream_pipelines/07_phyloplace_nextflow.sh

set -euo pipefail

SHEETS=data/raw/placement/phyloplace_samplesheets
MATOU_FAA=/home/aauladell/datasets/common_databases/MATOU/matou_alltranslations_filt10aa.faa
OUT=data/raw/placement

# HK genes
nextflow run nf-core/phyloplace \
    --phylosearch_input ${SHEETS}/hk_genes.csv \
    --search_fasta ${MATOU_FAA} \
    -profile singularity --outdir ${OUT}

# Some genes (hk_genes_nopremask.csv: e.g. HSP70, HSP90, with unspecific sequences) had to be placed with
# the EPA-ng option --no-premask. The exact call was not recorded in phylotranscriptomics.

# Non-HK genes: FAD1, Flavodoxin, NRT2, SIT
nextflow run nf-core/phyloplace \
    --phylosearch_input ${SHEETS}/nonhk_genes.csv \
    --search_fasta ${MATOU_FAA} \
    -profile singularity --outdir ${OUT}
