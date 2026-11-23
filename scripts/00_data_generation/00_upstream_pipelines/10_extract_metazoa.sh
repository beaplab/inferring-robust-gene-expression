#!/usr/bin/env bash
# Upstream 10: metazoan placements removed from the per-sample metaT jplaces of every HK gene.
#
# Per gene, the references whose lineage contains Metazoa form one clade (their lineage cut before
# ";Metazoa"); gappa edit extract splits every sample jplace into the queries placed in that clade and
# the rest ("basal_branches"; queries with < 95% of their mass in either go to "uncertain" and are dropped).
# The basal_branches jplaces, i.e. the non-metazoan placements, are the HK jplaces used downstream.
#
# Input:  data/raw/placement/epang_multiplicities_metaT/<gene>/<sample>filter.jplace.gz (09)
#         data/raw/placement/hk_gene_taxonomy/<gene>.taxonomy.txt
# Output: data/raw/placement/hk_matou_jplaces/<gene>/<sample>filter_basal_branches.jplace.gz
#         (symlink to gene_explorer data/placement/matou/placements_metazoa_extracted)
# Used by: fig5/01, fig6/01 and fig6/03.
# Run time: hours. Needs gappa.
# Original: gene_explorer main.nf, process extract_metazoans (scripts/modules/gappa_jplace_modifications.nf),
#   run on data/placement/matou_version (= epang_multiplicities_metaT). main.nf now holds the parameters of
#   a later metaG run, and the metaT outputs have no publishDir in the repository (copied by hand); this
#   script runs the process command directly. The original called `gappa prepare extract` (old gappa,
#   where the basal clade was named "basal_branches"); in gappa >= 0.8 it is `gappa edit extract` and the
#   name is set with --basal-clade-name.
# Checked (2026-11-20): gene ABCE, one sample: same queries and multiplicities as the jplace used.
#
# NOT RUN IN THIS PROJECT. If run, it writes into gene_explorer through the hk_matou_jplaces symlink.
#
# Usage (from the project root):
#   bash scripts/00_data_generation/00_upstream_pipelines/10_extract_metazoa.sh [n_parallel_jobs] [gene ...]

set -euo pipefail

NJOBS=${1:-10}
shift || true
GENES=("$@")

IN=data/raw/placement/epang_multiplicities_metaT
TAX=data/raw/placement/hk_gene_taxonomy
OUT=data/raw/placement/hk_matou_jplaces

extract_gene() {
    gene=$1; in=$2; tax=$3; out=$4
    mkdir -p ${out}/${gene}
    clades=$(mktemp)
    # reference <tab> lineage up to (not including) Metazoa: one clade with all metazoan references
    grep 'Metazoa' ${tax}/${gene}.taxonomy.txt | sed 's/;Metazoa.*//g' > ${clades}
    for jplace in ${in}/${gene}/*.jplace.gz; do
        gappa edit extract --jplace-path ${jplace} --clade-list-file ${clades} \
            --samples-file-prefix $(basename ${jplace} .jplace.gz)_ --samples-out-dir ${out}/${gene} \
            --basal-clade-name basal_branches --threads 1 --allow-file-overwriting > /dev/null
    done
    # only the non-metazoan part is kept
    find ${out}/${gene} -name '*.jplace' ! -name '*_basal_branches.jplace' -delete
    gzip -f ${out}/${gene}/*_basal_branches.jplace
    rm ${clades}
}
export -f extract_gene

[[ ${#GENES[@]} -eq 0 ]] && GENES=( $(ls ${IN}) )
printf '%s\n' "${GENES[@]}" | xargs -P ${NJOBS} -I{} bash -c 'extract_gene "$@"' _ {} ${IN} ${TAX} ${OUT}
