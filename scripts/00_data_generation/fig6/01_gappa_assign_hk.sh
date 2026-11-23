#!/usr/bin/env bash
# Figs 6 and 7: taxonomic assignment of the MATOU metaT placements of every housekeeping (HK) gene, per
# sample, and the metaT expression summed per taxonomic path.
#
# For each gene and sample:
#   gappa examine assign --best-hit   taxonomic path of the edge with the highest LWR of each query
#   gappa edit multiplicity           metaT expression (multiplicity) of each query
#   summary: per taxonomic path, number of queries (n_queries), summed expression (expression) and
#            summed expression x best-hit LWR (expression_lwr)
#
# Input:  data/raw/placement/hk_matou_jplaces/<gene>/<sample>filter_basal_branches.jplace.gz
#         data/raw/placement/hk_gene_taxonomy/<gene>.taxonomy.txt (EukProt-style path of each reference)
# Output: data/intermediate/hk_taxassign/<gene>/<sample>{per_query.tsv,profile.tsv,multiplicities.csv,summary.tsv}
# Run time: hours (~190,000 gappa runs) with 16 parallel jobs. Samples already done are skipped.
# Original: gene_explorer scripts/analysis/calculate_gappa-assign_matou.sh
#
# Usage (from the project root):
#   bash scripts/00_data_generation/fig6/01_gappa_assign_hk.sh [n_parallel_jobs] [gene ...]   (default: all genes)

set -euo pipefail

NJOBS=${1:-16}
shift || true
GENES=("$@")

TAXDIR=data/raw/placement/hk_gene_taxonomy
PLACEDIR=data/raw/placement/hk_matou_jplaces
OUTDIR=data/intermediate/hk_taxassign

assign_sample() {
    gene=$1; jplace=$2; taxfile=$3; outdir=$4
    sample=$(basename ${jplace} filter_basal_branches.jplace.gz)

    # already done
    [[ -s ${outdir}/${sample}summary.tsv ]] && return 0

    {
        gappa examine assign --jplace-path ${jplace} --taxon-file ${taxfile} \
            --out-dir ${outdir} --file-prefix ${sample} \
            --per-query-results --best-hit --threads 1 --allow-file-overwriting &&
        gappa edit multiplicity --jplace-path ${jplace} --write-multiplicity-file \
            --out-dir ${outdir} --file-prefix ${sample} --threads 1 --allow-file-overwriting
    } > ${outdir}/${sample}.log 2>&1 || { echo "FAILED: ${gene} ${sample}" >&2; return 0; }

    # per_query.tsv (name, LWR, ..., taxopath in column 6) + multiplicities.csv (sample, name, multiplicity)
    # -> one row per taxopath
    awk -F'\t' -v OFS='\t' -v gene=${gene} -v sample=${sample} '
        NR == FNR { mult[$2] = $3; next }
        FNR == 1  { next }
        {
            if (!($1 in mult)) { print "No multiplicity for " $1 > "/dev/stderr"; next }
            n[$6]++; expr[$6] += mult[$1]; exprlwr[$6] += $2
        }
        END {
            print "gene", "sample", "taxopath", "n_queries", "expression", "expression_lwr"
            for (t in n) print gene, sample, t, n[t], expr[t], exprlwr[t]
        }' \
        ${outdir}/${sample}multiplicities.csv ${outdir}/${sample}per_query.tsv > ${outdir}/${sample}summary.tmp \
    && mv ${outdir}/${sample}summary.tmp ${outdir}/${sample}summary.tsv
}
export -f assign_sample

[[ ${#GENES[@]} -eq 0 ]] && GENES=( $(ls ${PLACEDIR}) )

for gene in "${GENES[@]}"; do
    taxfile=${TAXDIR}/${gene}.taxonomy.txt
    if [[ ! -f ${taxfile} ]]; then echo "No taxonomy file for ${gene}, skipping" >&2; continue; fi
    mkdir -p ${OUTDIR}/${gene}
    for jplace in ${PLACEDIR}/${gene}/*.jplace.gz; do
        echo "${gene} ${jplace} ${taxfile} ${OUTDIR}/${gene}"
    done
done | xargs -P ${NJOBS} -n 4 bash -c 'assign_sample "$@"' _
