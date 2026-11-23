# Fig 2, check of step 01: how much the fastEmu results of P. calceolata depend on the reference set.
#
# Compares the runs of 01_pelagomonas_fastemu_iron.R in data/intermediate/fig2/fastemu_reference_comparison/
# (made with the earlier version of 01, which took the reference set as an argument):
#   selection/  reference set = the 73 HK genes of selection_houseg_denoms_bycorr-07.tsv found in Y (the paper)
#   all_hk/     reference set = every HK gene
# each with the robust test against fastEmu's data-driven reference (as in the paper; coeffs.tsv) and
# against the same HK reference set (coeffs_samereference.tsv), plus the cached paper table
# (paper/pelagomonas_fastemu_iron_coeffs.tsv; samples not aligned, see 01). all_hk/ (data-driven) is the
# table now used by Fig 2; its flavodoxin CAH0373297 was tested afterwards (always_test in 01).
#
# Significant: p <= 0.01, as in Fig 2. Direction as in Fig 2: log2 fold change high/low iron
# (estimate is the low_iron coefficient, natural log, so log2FC = -estimate / log(2)).
#
# Output: results/tables/fig2_fastemu_reference_comparison.tsv (one row per run: genes tested,
#         significant, up / down, reference set) and fig2_fastemu_reference_pairs.tsv (pairs of runs:
#         overlap of the significant genes, correlation and sign agreement of the estimates); the
#         flavodoxins of Fig 2 in each run are printed.
# Run time: seconds.
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig2/01b_compare_fastemu_references.R

source('R/load.R')

pval_cutoff <- 0.01
cmp_dir <- p('intermediate', 'fig2', 'fastemu_reference_comparison')
flavodoxins <- c('CAH0373297', 'CAH0367870', 'CAH0367856')

runs <- tribble(
    ~run,                              ~file,
    'paper (cached)',                  file.path(cmp_dir, 'paper', 'pelagomonas_fastemu_iron_coeffs.tsv'),
    'selection, robust data-driven',   file.path(cmp_dir, 'selection', 'pelagomonas_fastemu_iron_coeffs.tsv'),
    'selection, robust same ref.',     file.path(cmp_dir, 'selection', 'pelagomonas_fastemu_iron_coeffs_samereference.tsv'),
    'all HK, robust data-driven',      file.path(cmp_dir, 'all_hk', 'pelagomonas_fastemu_iron_coeffs.tsv'),
    'all HK, robust same ref.',        file.path(cmp_dir, 'all_hk', 'pelagomonas_fastemu_iron_coeffs_samereference.tsv')
) %>%
    filter(file.exists(file))

coefs <- runs %>%
    mutate(data = map(file, ~read_tsv(.x, show_col_types = FALSE))) %>%
    select(-file) %>%
    unnest(data) %>%
    mutate(log2fc = -estimate / log(2),
           significant = !is.na(pval) & pval <= pval_cutoff,
           run = factor(run, levels = runs$run))

# The first fits (before the robust test): estimates against each HK reference set
fits <- tibble(reference = c('selection', 'all_hk')) %>%
    mutate(file = file.path(cmp_dir, reference, 'pelagomonas_fastemu_iron_fit.rds')) %>%
    filter(file.exists(file))
if (nrow(fits) == 2) {
    fits <- fits %>%
        mutate(coef = map(file, ~as_tibble(read_rds(.x)$coef))) %>%
        select(-file) %>%
        unnest(coef) %>%
        mutate(to_test = (lower > 0 | upper < 0) & abs(estimate) >= 0.4)
    wide <- fits %>%
        select(reference, category, estimate, to_test) %>%
        pivot_wider(names_from = reference, values_from = c(estimate, to_test))
    message('First fits: ', nrow(wide), ' genes; estimate correlation ',
            round(cor(wide$estimate_selection, wide$estimate_all_hk), 4),
            ', median shift (all_hk - selection) ', round(median(wide$estimate_all_hk - wide$estimate_selection), 4),
            ', sd of the difference ', round(sd(wide$estimate_all_hk - wide$estimate_selection), 4))
    message('Genes passed to the robust test: selection ', sum(wide$to_test_selection), ', all_hk ',
            sum(wide$to_test_all_hk), ', both ', sum(wide$to_test_selection & wide$to_test_all_hk))
}

# Per run ---------------------------------------------------------------------------------------
summary_runs <- coefs %>%
    group_by(run) %>%
    summarize(genes = n(),
              tested = sum(!is.na(pval)),
              up_high_iron = sum(significant & log2fc > 0),
              down_high_iron = sum(significant & log2fc < 0),
              reference_genes = sum(reference_category),
              median_log2fc_significant = median(log2fc[significant]),
              n_significant = sum(significant),
              .groups = 'drop') %>%
    relocate(n_significant, .after = tested)
print(summary_runs, width = Inf)

# Pairs of runs ---------------------------------------------------------------------------------
pair_stats <- function(a, b) {
    x <- filter(coefs, run == a); y <- filter(coefs, run == b)
    sx <- x$category[x$significant]; sy <- y$category[y$significant]
    both <- inner_join(select(x, category, log2fc, sig_x = significant), select(y, category, log2fc, sig_y = significant),
                       by = 'category')
    both_sig <- filter(both, sig_x & sig_y)
    tibble(run_a = a, run_b = b,
           significant_a = length(sx), significant_b = length(sy),
           significant_both = length(intersect(sx, sy)),
           jaccard = length(intersect(sx, sy)) / length(union(sx, sy)),
           cor_log2fc_all = cor(both$log2fc.x, both$log2fc.y),
           median_shift_log2fc = median(both$log2fc.y - both$log2fc.x),
           sign_agreement_both_significant = mean(sign(both_sig$log2fc.x) == sign(both_sig$log2fc.y)))
}
pairs <- (if (nrow(runs) > 1) combn(as.character(runs$run), 2, simplify = FALSE) else list()) %>%
    map(~pair_stats(.x[1], .x[2])) %>%
    bind_rows()
print(pairs, width = Inf, n = Inf)

# Flavodoxins of Fig 2 ------------------------------------------------------------------------------
coefs %>%
    filter(str_detect(category, str_c(flavodoxins, collapse = '|'))) %>%
    transmute(run, gene = str_extract(category, 'CAH\\d+'), log2fc = round(log2fc, 2), pval = signif(pval, 2), significant) %>%
    arrange(gene, run) %>%
    print(n = Inf)

dir.create(p('tables'), showWarnings = FALSE)
write_tsv(summary_runs, p('tables', 'fig2_fastemu_reference_comparison.tsv'))
write_tsv(pairs, p('tables', 'fig2_fastemu_reference_pairs.tsv'))
