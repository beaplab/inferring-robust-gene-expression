#!/usr/bin/env bash
# Supp S3, step 01: MATOU v1.5 unigenes of Bathycoccus prasinos, by blastn of its transcriptome (EukProt
# EP00893) against the MATOU nucleotide catalogue, >= 92% identity.
#
# Input:  data/raw/genomic/transcriptomes_nucleotide/EP00893_Bathycoccus_prasinos.fna.gz
#         MATOU blast database (external$matou/MATOU-v1.5blast, made with
#         makeblastdb -in MATOU-v1.5.fna -dbtype nucl -out MATOU-v1.5blast)
# Output: data/intermediate/supplementary/S3/bathycoccus_prasinos_blast_92id.tsv
#         (the copy used for the paper: data/raw/quantification/statistics/bathycoccus_prasinos_blast_92id.tsv)
# Run time: hours (15 threads).
# Original: delving_into_florenciella src/analysis/comparing_genomes_w_MATOUv1_5.sh, run over five genomes;
#   the script says IDENTITY=90, but the table used is the _92id run (lowest identity in it: 92.0), so
#   the identity was 92 for that run. The coverage filter of the original (awk qcovs >= 80 -> filtered_ids)
#   is applied in step 02 instead.
#
# Usage (from the project root): bash scripts/00_data_generation/supplementary/S3/01_blast_bathycoccus_vs_matou.sh

set -euo pipefail

OUT=data/intermediate/supplementary/S3
DB=/home/aauladell/datasets/common_databases/MATOU/MATOU-v1.5blast
IDENTITY=92
mkdir -p ${OUT}

query=$(mktemp --suffix .fna)
zcat data/raw/genomic/transcriptomes_nucleotide/EP00893_Bathycoccus_prasinos.fna.gz > ${query}

blastn -query ${query} -db ${DB} \
    -out ${OUT}/bathycoccus_prasinos_blast_${IDENTITY}id.tsv \
    -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore" \
    -num_threads 15 -perc_identity ${IDENTITY} -evalue 1e-5
rm ${query}
