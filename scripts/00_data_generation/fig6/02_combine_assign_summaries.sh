#!/usr/bin/env bash
# Figs 6 and 7: one table per housekeeping gene with the per-sample summaries of 01_gappa_assign_hk.sh,
# so later steps read one file per gene instead of ~580.
#
# Input:  data/intermediate/hk_taxassign/<gene>/<sample>summary.tsv
# Output: data/intermediate/hk_taxassign_combined/<gene>_summary.tsv
#         (gene, sample, taxopath, n_queries, expression, expression_lwr) and <gene>.nsamples
#         (the copy used for the paper is data/raw/placement/hk_taxassign_combined)
# Run time: minutes. Genes already combined from the same number of samples are skipped.
# Original: gene_explorer scripts/analysis/combine_taxassign_summaries_matou.sh
#
# Usage (from the project root): bash scripts/00_data_generation/fig6/02_combine_assign_summaries.sh [n_parallel_jobs]

set -euo pipefail

NJOBS=${1:-16}
INDIR=data/intermediate/hk_taxassign
OUTDIR=data/intermediate/hk_taxassign_combined
mkdir -p ${OUTDIR}

combine_gene() {
    gene=$1; indir=$2; outdir=$3
    files=( ${indir}/${gene}/*summary.tsv )
    [[ -e ${files[0]} ]] || return 0
    nfiles=${#files[@]}

    # already combined from the same number of samples
    if [[ -s ${outdir}/${gene}_summary.tsv && -f ${outdir}/${gene}.nsamples ]] \
        && [[ $(cat ${outdir}/${gene}.nsamples) -eq ${nfiles} ]]; then
        return 0
    fi

    # header of the first file only
    awk 'FNR > 1 || NR == 1' "${files[@]}" > ${outdir}/${gene}_summary.tmp \
    && mv ${outdir}/${gene}_summary.tmp ${outdir}/${gene}_summary.tsv \
    && echo ${nfiles} > ${outdir}/${gene}.nsamples
}
export -f combine_gene

ls ${INDIR} | xargs -P ${NJOBS} -I{} bash -c 'combine_gene "$@"' _ {} ${INDIR} ${OUTDIR}
