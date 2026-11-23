#!/usr/bin/env bash
# Upstream 05: protein database searched by phyloplace (07): the MATOU v1.5 unigenes translated in the
# six frames, keeping translations of >= 10 aa (seqkit translate -f 6 -m 10 --trim).
#
# Input:  MATOU-v1.5.fna.gz (external$matou in config/paths.R)
# Output: matou_alltranslations_filt10aa.faa in the same folder (external$matou_translations, 70 GB)
# Run time: hours.
# Original: MATOU folder, creation_faa.sh
#
# NOT RUN IN THIS PROJECT.
#
# Usage (from the project root): bash scripts/00_data_generation/00_upstream_pipelines/05_matou_translations.sh

set -euo pipefail

MATOU=/home/aauladell/datasets/common_databases/MATOU

pigz -dk ${MATOU}/MATOU-v1.5.fna.gz
# all 6 frames, >= 10 aa, stop symbols trimmed
seqkit translate -m 10 -f 6 --trim ${MATOU}/MATOU-v1.5.fna > ${MATOU}/matou_alltranslations_filt10aa.faa
