# Fig 6 (and the HK denominators of Fig 7): filter of the MATOU housekeeping (HK) genes.
#
# Removal criteria (keep_core = FALSE if any):
#   flag_lba       long-branch attraction: > 50% of the placement mass on the 1% longest edges of the tree
#                  and the reads far from the references (multiplicity-weighted median pendant length > 2)
#   flag_taxgap    taxonomic gap: < 20% of the reference species of at least one of Chlorophyta, diatoms,
#                  dinoflagellates, haptophytes and ciliates present among the references
#   flag_fragment  fragmentary references: > 15% of the references shorter than half the median length
#   flag_manual    erf3b (most of its abundance on a 99 aa reference fragment)
# Only flagged, not removed: flag_lowlwr, median LWR < 0.6 (imprecise placement, not wrong abundance).
# The placement metrics are the median over four samples (both depths, several size fractions).
#
# Input:  data/intermediate/fig5/gene_stats_summary.tsv            (fig5/01_tree_and_placement_stats.R)
#         data/raw/placement/hk_matou_jplaces/<gene>/              (placements)
#         data/raw/gene_data/hk_reference_alignments/<gene>.fasta  (reference sequences)
#         data/raw/gene_data/taxonomy_by_species.txt               (lineage of each reference species)
# Output: data/intermediate/fig6/housekeeping_gene_filter_matou.tsv (one row per HK gene: flags,
#         keep_core and the metrics)
# Run time: ~5 min with 16 cores.
# Original: gene_explorer scripts/analysis/housekeeping_gene_filter_matou.R
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig6/03_housekeeping_gene_filter.R [ncores]

source('R/load.R')
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
ncores <- if (length(args) >= 1) as.integer(args[1]) else 16

jplace_dir <- p('placement', 'hk_matou_jplaces')
fasta_dir <- p('gene_data', 'hk_reference_alignments')

# samples for the placement metrics: both depths and several size fractions
samples <- c('128DCM1GGZZ14filter', '100SUR1SSUU14filter', '110SUR1QQSS14filter', '122SUR1MMQQ14filter')

# lineages whose reference coverage is computed (patterns on the reference taxonomy); the taxonomic gap
# uses the main marine ones (key_groups)
groups <- c(Chlorophyta = 'Chlorophyta', Diatoms = 'Bacillariophyta|Diatom', Dinoflagellata = 'Dinoflagellata',
            Ciliophora = 'Ciliophora', Haptophyta = 'Haptophyta', Cryptophyceae = 'Cryptophyceae',
            Bigyra = 'Bigyra', Rhizaria = 'Rhizaria', Pelagophyceae = 'Pelagophyceae')
key_groups <- c('Chlorophyta', 'Diatoms', 'Dinoflagellata', 'Haptophyta', 'Ciliophora')
manual_removal <- c('erf3b')

gene_stats <- read_tsv(p('intermediate', 'fig5', 'gene_stats_summary.tsv'), show_col_types = FALSE)
tax <- read_tsv(p('gene_data', 'taxonomy_by_species.txt'), col_names = c('sp', 'tax'), show_col_types = FALSE)
species_per_group <- map_int(groups, ~sum(grepl(.x, tax$tax)))

# Reference sequences: share of the species of each group present, and fragments
reference_metrics <- function(gene) {
    lines <- readLines(file.path(fasta_dir, str_c(gene, '.fasta')))
    header_i <- grep('^>', lines)
    ids <- sub(' .*', '', sub('^>', '', lines[header_i]))
    seqs <- map_chr(seq_along(header_i), ~paste(lines[(header_i[.x] + 1):(if (.x < length(header_i)) header_i[.x + 1] - 1 else length(lines))],
                                                collapse = ''))
    len <- nchar(gsub('[-.]', '', seqs))   # alignment gaps removed
    lineages <- tax$tax[match(ids, tax$sp)]
    coverage <- map_dbl(groups, ~sum(grepl(.x, lineages))) / species_per_group
    bind_cols(as_tibble(as.list(coverage)),
              tibble(frac_fragment = mean(len < 0.5 * median(len)),
                     med_len = median(len)))
}

# Placements of one sample: best placement of each pquery, weighted by multiplicity
placement_metrics <- function(gene, sample) {
    f <- file.path(jplace_dir, gene, str_c(sample, '_basal_branches.jplace.gz'))
    if (!file.exists(f)) return(NULL)
    j <- jsonlite::fromJSON(f, simplifyVector = FALSE)
    if (length(j$placements) == 0) return(NULL)
    fields <- unlist(j$fields)
    m <- str_match_all(j$tree, ':([0-9.eE-]+)\\{(\\d+)\\}')[[1]]
    edges <- tibble(edge = as.integer(m[, 3]), len = as.numeric(m[, 2]))
    long_edge <- quantile(edges$len, 0.99)
    pq <- map_dfr(j$placements, function(p) {
        lwr <- map_dbl(p$p, ~.x[[which(fields == 'like_weight_ratio')]])
        best <- p$p[[which.max(lwr)]]
        tibble(edge = as.integer(best[[which(fields == 'edge_num')]]),
               pend = best[[which(fields == 'pendant_length')]],
               mult = if (is.null(p$nm)) length(p$n) else sum(map_dbl(p$nm, 2)))
    }) %>%
        left_join(edges, by = 'edge')
    tibble(top2_frac = pq %>% count(edge, wt = mult) %>% slice_max(n, n = 2, with_ties = FALSE) %>%
               pull(n) %>% sum() / sum(pq$mult),
           longbranch_frac = sum(pq$mult[pq$len >= long_edge]) / sum(pq$mult),
           wmed_pend = matrixStats::weightedMedian(pq$pend, pq$mult))
}

gene_metrics <- function(gene) {
    placements <- map_dfr(samples, ~placement_metrics(gene, .x)) %>%
        summarize(across(everything(), median))
    bind_cols(tibble(houseg = gene), reference_metrics(gene), placements)
}

hk_genes <- gene_stats %>%
    filter(gene_group == 'housekeeping', file.exists(file.path(fasta_dir, str_c(houseg, '.fasta')))) %>%
    pull(houseg)
metrics <- bind_rows(mclapply(hk_genes, gene_metrics, mc.cores = ncores))

gene_filter <- gene_stats %>%
    select(houseg, label, ntips, median_multiplicity_per_sample, lwr_median) %>%
    inner_join(metrics, by = 'houseg') %>%
    mutate(min_cov = pmin(!!!syms(key_groups)),
           worst_group = key_groups[apply(select(., all_of(key_groups)), 1, which.min)],
           flag_lba = longbranch_frac > 0.5 & wmed_pend > 2,
           flag_taxgap = min_cov < 0.2,
           flag_fragment = frac_fragment > 0.15,
           flag_manual = houseg %in% manual_removal,
           flag_lowlwr = lwr_median < 0.6,
           keep_core = !(flag_lba | flag_taxgap | flag_fragment | flag_manual)) %>%
    mutate(across(where(is.double), ~signif(.x, 3))) %>%
    select(houseg, label, keep_core, flag_lba, flag_taxgap, flag_fragment, flag_manual, flag_lowlwr,
           longbranch_frac, wmed_pend, top2_frac, min_cov, worst_group, frac_fragment,
           lwr_median, median_multiplicity_per_sample, ntips) %>%
    arrange(desc(keep_core), houseg)

dir.create(p('intermediate', 'fig6'), showWarnings = FALSE, recursive = TRUE)
write_tsv(gene_filter, p('intermediate', 'fig6', 'housekeeping_gene_filter_matou.tsv'))

gene_filter %>%
    summarize(genes = n(), kept = sum(keep_core), long_branch = sum(flag_lba), taxonomic_gap = sum(flag_taxgap),
              fragment = sum(flag_fragment), manual = sum(flag_manual), low_lwr_flag = sum(flag_lowlwr)) %>%
    print()
