# Supplementary Table T1: reference trees, placement quality and taxonomic composition of the 318
# candidate housekeeping (HK) genes (Strassert et al. 2021) after placement of the MATOU unigenes.
#
# Per gene, the HK gene filter table (fig6/03: placement metrics of fig5/01, coverage of the key marine
# groups, fragment references, filter flags), the total branch length of the reference tree (fig5/01)
# and three composition columns, computed from the taxonomic composition of each gene in every sample
# (fig6/04):
#   other_frac  mean fraction of "Other" (Fornicata, Metamonada, Breviatea) over all samples
#   dist_metaT  Hellinger distance between the gene's mean composition and that of the whole MATOU metaT
#               catalogue (0: identical, 1: no overlap); as in gene_explorer
#               scripts/figures/pca_taxonomy_housekeeping_filter_matou.R ("All samples" panel)
#   dominant    most abundant lineage of the mean composition, with its fraction
# reasons: the failed criteria, in this order: LBA, taxgap, fragment, manual (these define the core
# subset, keep_core) and the flags lowLWR* (lwr_median < 0.6) and Other>0.4* (other_frac > 0.4).
# Rows: kept genes first, then removed ones, each by functional category and gene name.
# Mean composition: every sample with the same weight, taxa absent from a sample counted as 0.
#
# Input:  data/intermediate/fig6/housekeeping_gene_filter_matou.tsv, data/intermediate/fig6/taxon_compositions.parquet,
#         data/intermediate/fig5/gene_stats_summary.tsv
# Original: gene_explorer results/tables/housekeeping_gene_summary_matou.tsv (written interactively; no
#   script in the repository) -> manuscript/supplementary_proposal/T1_candidate_HK_genes_summary.tsv
#   (columns trimmed).
# Run time: seconds.
# Output: results/tables/supplementary/T1_candidate_HK_genes_summary.tsv

source('R/load.R')
library(arrow, warn.conflicts = FALSE)

gene_filter <- read_tsv(p('intermediate', 'fig6', 'housekeeping_gene_filter_matou.tsv'), show_col_types = FALSE)
compositions <- read_parquet(p('intermediate', 'fig6', 'taxon_compositions.parquet'))
tree_stats <- read_tsv(p('intermediate', 'fig5', 'gene_stats_summary.tsv'), show_col_types = FALSE) %>%
    select(houseg, branch_length)
OTHER_MAX <- 0.4

# mean composition of each source (gene or reference catalogue) over its samples; sources x taxa
mean_composition <- compositions %>%
    complete(nesting(source, samplename), taxon, fill = list(rel_value = 0)) %>%
    group_by(source, taxon) %>%
    summarize(rel_value = mean(rel_value), .groups = 'drop')

hellinger <- mean_composition %>%
    mutate(rel_value = sqrt(rel_value)) %>%
    pivot_wider(names_from = taxon, values_from = rel_value, values_fill = 0)
mat <- as.matrix(select(hellinger, -source))
rownames(mat) <- hellinger$source
dist_metaT <- sqrt(rowSums(sweep(mat, 2, mat['MATOU metaT (all genes)', ])^2)) / sqrt(2)

composition_columns <- mean_composition %>%
    group_by(source) %>%
    summarize(other_frac = sum(rel_value[taxon == 'Other']),
              dominant = sprintf('%s (%.2f)', taxon[which.max(rel_value)], max(rel_value)),
              .groups = 'drop') %>%
    mutate(dist_metaT = dist_metaT[source])

reason_flags <- c(flag_lba = 'LBA', flag_taxgap = 'taxgap', flag_fragment = 'fragment', flag_manual = 'manual',
                  flag_lowlwr = 'lowLWR*', flag_other = str_c('Other>', OTHER_MAX, '*'))

t1 <- gene_filter %>%
    left_join(tree_stats, by = 'houseg') %>%
    left_join(composition_columns, by = c('houseg' = 'source')) %>%
    mutate(flag_other = other_frac > OTHER_MAX) %>%
    rowwise() %>%
    mutate(reasons = str_c(reason_flags[unlist(c_across(all_of(names(reason_flags))))], collapse = ', ')) %>%
    ungroup() %>%
    arrange(desc(keep_core), label, houseg) %>%
    select(houseg, label, longbranch_frac, wmed_pend, top2_frac, min_cov, worst_group, frac_fragment,
           lwr_median, median_multiplicity_per_sample, ntips, branch_length, other_frac, dist_metaT,
           dominant, reasons)
stopifnot(nrow(t1) == 318, !anyNA(t1$dist_metaT))

dir.create(p('tables', 'supplementary'), showWarnings = FALSE)
write_tsv(t1, p('tables', 'supplementary', 'T1_candidate_HK_genes_summary.tsv'), na = '')
