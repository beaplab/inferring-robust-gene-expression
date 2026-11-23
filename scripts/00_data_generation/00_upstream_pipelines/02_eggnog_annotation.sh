#!/usr/bin/env bash
# Upstream 02: functional annotation (eggNOG-mapper) of the CDS of every transcriptome / genome.
#
# Input:  data/raw/genomic/transcriptomes_nucleotide/*.fna.gz
#         eggNOG database (external$eggnog_db in config/paths.R)
# Output: data/raw/genomic/functional_annotation/<transcriptome>.emapper.annotations (+ .hits,
#         .seed_orthologs). Used by Fig 2 (P. calceolata DE, flavodoxins), Fig 3 (fig3/03) and Fig 4.
# Run time: hours. Needs eggnog-mapper (conda environment "eggnog" in the original run).
# Original: delving_into_florenciella src/01_eggnog-annotation.sh (it also annotated the hand-downloaded
#           EukProt v3 transcriptomes, a folder that no longer exists)
#
# NOT RUN IN THIS PROJECT. If run, it overwrites the annotations in data/raw/genomic/functional_annotation.
#
# Usage (from the project root): bash scripts/00_data_generation/00_upstream_pipelines/02_eggnog_annotation.sh

set -euo pipefail

IN=data/raw/genomic/transcriptomes_nucleotide
OUT=data/raw/genomic/functional_annotation
DB=/scratch/data1/aauladell/databases/eggnog
mkdir -p ${OUT}

for transcriptome in ${IN}/*.fna.gz; do
    name=$(basename ${transcriptome} .fna.gz)
    # CDS input: translated and searched (--genepred search)
    emapper.py --itype CDS --genepred search --cpu 12 \
        -i ${transcriptome} --data_dir ${DB} -o ${OUT}/${name}
done
