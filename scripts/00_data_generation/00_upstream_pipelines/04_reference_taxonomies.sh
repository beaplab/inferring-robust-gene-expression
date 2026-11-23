#!/usr/bin/env bash
# Upstream 04: taxonomy files (reference name <tab> lineage) of the SIT, NRT2 and Flavodoxin references,
# used by gappa examine assign (fig7/02) and by the metazoan extraction (10).
#
# Input:  data/raw/placement/phyloplace_references/nonhk_genes/<gene>_msa.fasta (03a/b/c)
#         data/raw/gene_data/taxonomy_by_species.txt
# Output: data/raw/placement/hk_gene_taxonomy/<gene>.taxonomy.txt
#         The 318 HK gene taxonomies in the same folder are those of Richter et al. (external
#         $hk_taxonomies_richter, identical files), not made here.
# Run time: seconds.
# Original: gene_explorer scripts/data_generation/generate_nonHK_gene_trees/helper_create-nonhk-taxonomy.sh
#           + scripts/data_generation/helper_create-housekeeping-taxonomy.R (helpers/reference_taxonomy.R)
# Checked (2026-11-20): NRT2 and Flavodoxin byte-identical to the files used; SIT differs only in 12 names
#           that still had "|" in 2024 (handled in fig7/02).
#
# NOT RUN IN THIS PROJECT. If run, it overwrites the files in gene_explorer through the hk_gene_taxonomy symlink.
#
# Usage (from the project root): bash scripts/00_data_generation/00_upstream_pipelines/04_reference_taxonomies.sh

set -euo pipefail

REFS=data/raw/placement/phyloplace_references/nonhk_genes
OUT=data/raw/placement/hk_gene_taxonomy

for gene in SIT NRT2 Flavodoxin; do
    Rscript scripts/00_data_generation/00_upstream_pipelines/helpers/reference_taxonomy.R \
        --fasta ${REFS}/${gene}_msa.fasta --out ${OUT}/${gene}.taxonomy.txt
done
