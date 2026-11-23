# Fig 5 (and the HK gene filter of Fig 6): per-gene statistics of the reference trees and of the
# placements of the MATOU metaT unigenes on them, for every housekeeping (HK) gene and SIT, NRT2 and
# Flavodoxin.
#
# Per gene, from its per-sample jplaces (metazoan placements already removed):
#   tree:       number of tips and nodes, total branch length, treeness (share of the branch length on
#               internal branches), branch length per tip; from the tree stored in the jplace
#   placements: per sample, number of pqueries and summed multiplicity (metaT)
#   LWR:        for each pquery, the like-weight ratio (LWR) of its best placement; histogram (bins of
#               0.05), mean, multiplicity-weighted mean, median, mode, share > 0.5 and > 0.9; mean number
#               of alternative placements; median pendant length of the best placement
#   lwr_type:   low (mode <= 0.10), intermediate (<= 0.5) or high LWR
#
# Input:  data/raw/placement/hk_matou_jplaces/<gene>/<sample>filter_basal_branches.jplace.gz (15 GB)
#         data/raw/gene_data/gene_uniprot-functional-summary.tsv (functional label of each HK gene)
# Output: data/intermediate/fig5/tree_stats_cache/<gene>.rds  one per gene (a rerun only does the
#                                                              missing genes)
#         data/intermediate/fig5/gene_stats_summary.tsv       one row per gene
#         data/intermediate/fig5/per_sample_placements.tsv.gz
#         data/intermediate/fig5/lwr_distribution_recalc.tsv.gz
#         The three tables are only written once every gene has its cache file.
# Run time: hours (reads every jplace; ~1 min per gene with 4 cores).
# Original: gene_explorer scripts/statistics/stat_phylogenetic_tree_summary_matou.R
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/fig5/01_tree_and_placement_stats.R [gene ...]   (default: all genes)

source('R/load.R')
library(parallel)

jplace_dir <- p('placement', 'hk_matou_jplaces')
out_dir <- p('intermediate', 'fig5')
cache_dir <- file.path(out_dir, 'tree_stats_cache')
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

target_genes <- c('SIT', 'NRT2', 'Flavodoxin')
lwr_breaks <- round(seq(0, 1, by = 0.05), 2)

all_genes <- list.files(jplace_dir)
genes <- commandArgs(trailingOnly = TRUE)
if (length(genes) == 0) genes <- all_genes

# Statistics of one gene ------------------------------------------------------------------------------
gene_statistics <- function(gene) {
    cache_file <- file.path(cache_dir, str_c(gene, '.rds'))
    if (file.exists(cache_file)) return(readRDS(cache_file))

    files <- list.files(file.path(jplace_dir, gene), pattern = 'jplace.gz$', full.names = TRUE)
    samples <- str_remove(basename(files), '_basal_branches.jplace.gz')
    parsed <- map(files, read_pqueries)
    pq <- map2_dfr(parsed, samples, ~mutate(.x$pq, sample = .y))

    per_sample <- pq %>%
        group_by(sample) %>%
        summarize(n_pqueries = n(), total_multiplicity = sum(multiplicity)) %>%
        # samples without placements are kept, with 0
        right_join(tibble(sample = samples), by = 'sample') %>%
        replace_na(list(n_pqueries = 0, total_multiplicity = 0)) %>%
        mutate(houseg = gene, .before = 1)

    lwr_hist <- pq %>%
        mutate(bin_start = lwr_breaks[findInterval(lwr1, lwr_breaks, rightmost.closed = TRUE)]) %>%
        group_by(bin_start) %>%
        summarize(count = n(), weighted_mass = sum(multiplicity)) %>%
        right_join(tibble(bin_start = head(lwr_breaks, -1)), by = 'bin_start') %>%
        replace_na(list(count = 0, weighted_mass = 0)) %>%
        arrange(bin_start) %>%
        mutate(frac = count / sum(count),
               weighted_frac = weighted_mass / sum(weighted_mass),
               houseg = gene, .before = 1)

    lwr_summary <- pq %>%
        summarize(lwr_mean = mean(lwr1),
                  lwr_weighted_mean = weighted.mean(lwr1, multiplicity),
                  lwr_median = median(lwr1),
                  frac_lwr_gt0.5 = mean(lwr1 > 0.5),
                  frac_lwr_gt0.9 = mean(lwr1 > 0.9),
                  mean_n_alt_placements = mean(n_alt),
                  median_pendant_length = median(pendant1)) %>%
        mutate(lwr_mode = lwr_hist$bin_start[which.max(lwr_hist$count)],
               houseg = gene, .before = 1)

    res <- list(tree = tree_stats(parsed[[1]]$tree) %>% mutate(houseg = gene, .before = 1),
                per_sample = per_sample,
                lwr_hist = lwr_hist,
                lwr_summary = lwr_summary)
    saveRDS(res, cache_file)
    res
}

done <- mclapply(genes, function(g) try(gene_statistics(g)), mc.cores = 4, mc.preschedule = FALSE)
failed <- genes[map_lgl(done, ~inherits(.x, 'try-error'))]
if (length(failed) > 0) warning('Failed genes: ', str_c(failed, collapse = ', '))

# Summary tables, only with every gene done ---------------------------------------------------------------
missing <- all_genes[!file.exists(file.path(cache_dir, str_c(all_genes, '.rds')))]
if (length(missing) > 0) {
    message(length(missing), ' genes without statistics yet; summary tables not written')
    quit(save = 'no')
}

results <- map(file.path(cache_dir, str_c(all_genes, '.rds')), readRDS)
per_sample <- map_dfr(results, 'per_sample')
write_tsv(per_sample, file.path(out_dir, 'per_sample_placements.tsv.gz'))
write_tsv(map_dfr(results, 'lwr_hist'), file.path(out_dir, 'lwr_distribution_recalc.tsv.gz'))

placements <- per_sample %>%
    group_by(houseg) %>%
    summarize(total_placements = sum(n_pqueries),
              median_placements_per_sample = median(n_pqueries),
              n_samples_with_placements = sum(n_pqueries > 0),
              # before the sum, which overwrites the column
              median_multiplicity_per_sample = median(total_multiplicity),
              total_multiplicity = sum(total_multiplicity))

functional_labels <- read_tsv(p('gene_data', 'gene_uniprot-functional-summary.tsv'), show_col_types = FALSE) %>%
    select(houseg, label)

map_dfr(results, 'tree') %>%
    left_join(placements, by = 'houseg') %>%
    left_join(map_dfr(results, 'lwr_summary'), by = 'houseg') %>%
    mutate(lwr_type = case_when(lwr_mode <= 0.10 ~ 'low LWR',
                                lwr_mode > 0.10 & lwr_mode <= 0.5 ~ 'intermediate LWR',
                                lwr_mode > 0.5 ~ 'high LWR',
                                TRUE ~ 'other'),
           gene_group = if_else(houseg %in% target_genes, 'target', 'housekeeping')) %>%
    left_join(functional_labels, by = 'houseg') %>%
    write_tsv(file.path(out_dir, 'gene_stats_summary.tsv'))
