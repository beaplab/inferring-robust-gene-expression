#!/usr/bin/env bash
# Fig 3A: amino-acid identity (AAI) between the Florenciella transcriptomes, the three Florenciella
# SMAGs and Pelagomonas calceolata (outgroup), estimated by sourmash from protein k-mer sketches.
#
# Input:  predicted proteins of each transcriptome / genome (delving_into_florenciella data/genomic_data,
#         not copied into this project): EukProt transcriptomes (*.fasta), SMAG gene models
#         (*.gmove.pep.faa), NCBI genomes (*/*.faa.gz) and the Pelagomonas proteome
# Output: data/intermediate/fig3/compar_transcriptomes.csv (pairwise AAI matrix; it also has the
#         TOPAZ genomes, which the figure drops)
# Run time: minutes. Needs sourmash (conda environment "smash" in the original run).
# Original: delving_into_florenciella src/00_genome-down-and-selection/01_comparison-ani-genomes.sh
#           (which also compared the genomes at the nucleotide level, not used in the figure)
#
# Usage (from the project root): bash scripts/00_data_generation/fig3/01_aai_sourmash.sh

set -euo pipefail

GENOMIC=/scratch/data3/aauladell/project_big/delving_into_florenciella/data/genomic_data
SIGS=$(mktemp -d)
OUT=data/intermediate/fig3
mkdir -p ${OUT}

# protein sketches of every proteome
sourmash sketch protein \
    ${GENOMIC}/genomes/*/*.faa.gz \
    ${GENOMIC}/genomes/*.gmove.pep.faa \
    ${GENOMIC}/transcriptomes/*.fasta \
    ${GENOMIC}/genomes/Pelagomonas_outgroup/pelagomonas_aa.faa \
    --output-dir ${SIGS}

# pairwise AAI (--ani on protein sketches)
sourmash compare --ani --csv ${OUT}/compar_transcriptomes.csv ${SIGS}/*.sig

rm -r ${SIGS}
