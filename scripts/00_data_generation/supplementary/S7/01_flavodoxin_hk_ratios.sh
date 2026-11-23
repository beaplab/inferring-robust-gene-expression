#!/usr/bin/env bash
# Supp S7, step 01: metaT/HK ratios of the Flavodoxin unigenes, with the same script as SIT and NRT2 in
# Fig 7 (scripts/00_data_generation/fig7/02_metaT_HK_ratios.R; for Flavodoxin it fills in the 25 reference
# sequences with NA taxonomy first).
#
# Output: data/intermediate/supplementary/S7/Flavodoxin_hk_normalisation/ (assign/per_query.tsv is also
#         read by S6/01 step 3) and Flavodoxin_hkratio_jplace/
#         (the copy used for the paper: data/raw/placement/Flavodoxin_hkratios_jplace/Flavodoxin)
# Run time: ~10 min with 16 cores. Needs gappa.
# Original: gene_explorer scripts/analysis/generate_ratios_metaT_HK.R Flavodoxin
#
# Usage (from the project root): bash scripts/00_data_generation/supplementary/S7/01_flavodoxin_hk_ratios.sh [ncores]
set -euo pipefail
Rscript scripts/00_data_generation/fig7/02_metaT_HK_ratios.R Flavodoxin "${1:-16}"
