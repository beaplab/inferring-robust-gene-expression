#!/usr/bin/env bash
# Upstream 01: quantification of the metatranscriptomes against each transcriptome / genome, with the
# nextflow pipeline beaplab/transcriptome_metaT_quantification (branch main).
#
# For each transcriptome (nucleotide CDS) and each metatranscriptome of the fastq sheet: reads aligned
# with BWA-mem2 (alignment/) and quantified with Salmon (mapping/); per transcriptome and dataset, tables
# of TPM, read counts (Numreads) and gene characteristics, plus the filtered BAMs.
#
# Input:  data/raw/genomic/transcriptomes_nucleotide/*.fna.gz   (symlink to delving_into_florenciella)
#         data/raw/quantification/fastq_sheets/*.csv             (group, name, fastq_r1, fastq_r2, single_end,
#                                                                 sig; fastqs under /scratch/datasets_symbolic_links)
# Output: data/raw/quantification/{alignment,mapping}/
#         The project keeps only the tables (TPM_*, Numreads_*, gene-characteristics_*); the BAMs used by
#         fig3/03 and fig3/05 stay in delving_into_florenciella data/quantification/alignment/<transcriptome>/.
# Used by: Fig 2 (fig2/01), Fig 3 (fig3/03-05), Fig 4.
# Run time: days. Needs nextflow and the pipeline's containers.
# Original: delving_into_florenciella src/02_quantification.sh (+ src/02quant.config ->
#           config/nextflow/quantification.config)
#
# NOT RUN IN THIS PROJECT. If run, it writes (BAMs included, hundreds of GB) into data/raw/quantification.
#
# Usage (from the project root): bash scripts/00_data_generation/00_upstream_pipelines/01_quantification_nextflow.sh

set -euo pipefail

OUT=data/raw/quantification
SHEETS=data/raw/quantification/fastq_sheets

# 1. Florenciella transcriptomes, SMAGs and outgroups against every metatranscriptome dataset
nextflow run beaplab/transcriptome_metaT_quantification \
    --fastq_sheet ${SHEETS}/metatranscriptomes_datasets.csv \
    --transcriptome "data/raw/genomic/transcriptomes_nucleotide/*.fna.gz" \
    --outdir ${OUT} \
    -r main \
    -resume

# 2. Other EukProt v3 transcriptomes downloaded by hand (e.g. Rhizochromulina sp. CCMP1243, Phaeocystis
#    cordata). NOTE: this folder no longer exists; the call is kept as a record.
nextflow run beaplab/transcriptome_metaT_quantification \
    --fastq_sheet ${SHEETS}/metatranscriptomes_datasets.csv \
    --transcriptome "/home/aauladell/datasets/genomes_assemblies_and_genes/2022_eukprot_v3/some_downloaded_by-hand_datsets/*.fna.gz" \
    --outdir ${OUT} \
    -r main \
    -resume

# 3. Pelagomonas calceolata against the Tara samples of Carradec et al. (paired-end;
#    TPM_Pelagomonas_calceolata_quant-over_2012_carradec_tara_paired-end.tsv, used by fig2/01), 15 cpus
nextflow run beaplab/transcriptome_metaT_quantification \
    --fastq_sheet ${SHEETS}/only_new_carradec.csv \
    --transcriptome "data/raw/genomic/transcriptomes_nucleotide/Pelagomonas_calceolata.fna.gz" \
    --outdir ${OUT} \
    -c config/nextflow/quantification.config \
    -r main
