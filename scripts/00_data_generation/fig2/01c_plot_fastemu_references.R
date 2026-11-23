# Fig 2 (left panel), drawn for each fastEmu reference set, to compare them by eye.
#
# mode = robust (default): one column per robust score test, genes with p <= 0.01
#   Paper                     the cached table (selection first fit, data-driven reference)
#   Selection, data-driven    73 selected HK genes in the first fit, data-driven reference in the test
#   Selection, same ref.      73 selected HK genes in the first fit and in the test
#   All HK, data-driven       all 281 HK genes in the first fit, data-driven reference in the test
#   All HK, same ref.         all 281 HK genes in the first fit and in the test
# mode = first_fit: paper, and the two first fits (no p-values: a gene counts when its 95% interval
#   excludes 0)
# Top: every such gene (estimate and interval, violin per direction, counts). Bottom: the specific
# cases of Fig 2 (flavodoxins, ferredoxins, ISIP2A): every one that counts in at least one column, plus
# the non-significant flavodoxin of Fig 2, with the same rows in each column (faded where it does not
# count). Log2 fold change high/low iron, as in Fig 2.
#
# Input:  data/intermediate/fig2/fastemu_reference_comparison/paper/pelagomonas_fastemu_iron_coeffs.tsv (paper)
#         data/intermediate/fig2/fastemu_reference_comparison/{selection,all_hk}/
#           pelagomonas_fastemu_iron_fit_coeffs.tsv, pelagomonas_fastemu_iron_coeffs{,_samereference}.tsv
#         (earlier version of 01_pelagomonas_fastemu_iron.R: <selection|all_hk> both ...; all_hk data-driven
#         is the table of Fig 2)
# Output: results/figures/fig2_fastemu_reference_comparison_{robust,first_fit}.{pdf,png}
# Run time: seconds.
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig2/01c_plot_fastemu_references.R [robust|first_fit]

source('R/load.R')
library(patchwork)

pval_cutoff <- 0.01
cmp_dir <- p('intermediate', 'fig2', 'fastemu_reference_comparison')
isip2a <- 'CAH0367801.1|CAH0371101.1|CAH0367865.1|CAH0376816.1|CAH0367864.1'
flavodoxin_ns <- '9618'
# the flavodoxins highlighted in Fig 2
highlight <- c(CAH0373297 = '#3887d7', CAH0367870 = '#e31c1c', CAH0367856 = '#2eb019')
beige <- '#efe9df'

annotation <- read_tsv(p('annotation', 'Pelagomonas_calceolata.emapper.annotations'),
                       comment = '##', show_col_types = FALSE) %>%
    select(category = `#query`, PFAMs)

# fastEmu gives the low_iron coefficient (natural log): log2 fold change high/low = -x / log(2), which
# swaps the bounds of the interval. ci_excludes_0 is computed before.
to_log2 <- function(coef) {
    coef %>%
        mutate(ci_excludes_0 = lower > 0 | upper < 0,
               estimate = -estimate / log(2),
               new_lower = -upper / log(2), upper = -lower / log(2), lower = new_lower,
               direction = ifelse(estimate > 0, 'up', 'down')) %>%
        select(-new_lower)
}

mode <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(mode)) mode <- 'robust'
stopifnot(mode %in% c('robust', 'first_fit'))

read_robust <- function(file, column) {
    read_tsv(file, show_col_types = FALSE) %>%
        to_log2() %>% mutate(column = column, counts = !is.na(pval) & pval <= pval_cutoff)
}
read_first_fit <- function(reference) {
    read_tsv(file.path(cmp_dir, reference, 'pelagomonas_fastemu_iron_fit_coeffs.tsv'), show_col_types = FALSE) %>%
        to_log2() %>% mutate(column = reference, counts = ci_excludes_0)
}
paper <- read_robust(file.path(cmp_dir, 'paper', 'pelagomonas_fastemu_iron_coeffs.tsv'), 'paper')

if (mode == 'robust') {
    columns <- c(paper = 'Paper\n(selection; data-driven ref. in the test)',
                 sel_dd = 'Selection (73 HK)\ndata-driven ref. in the test',
                 sel_same = 'Selection (73 HK)\nsame ref. in the test',
                 all_dd = 'All HK (281)\ndata-driven ref. in the test',
                 all_same = 'All HK (281)\nsame ref. in the test')
    coefs <- bind_rows(
        paper,
        read_robust(file.path(cmp_dir, 'selection', 'pelagomonas_fastemu_iron_coeffs.tsv'), 'sel_dd'),
        read_robust(file.path(cmp_dir, 'selection', 'pelagomonas_fastemu_iron_coeffs_samereference.tsv'), 'sel_same'),
        read_robust(file.path(cmp_dir, 'all_hk', 'pelagomonas_fastemu_iron_coeffs.tsv'), 'all_dd'),
        read_robust(file.path(cmp_dir, 'all_hk', 'pelagomonas_fastemu_iron_coeffs_samereference.tsv'), 'all_same'))
    title_mode <- 'robust score tests, p ≤ 0.01'
} else {
    columns <- c(paper = 'Paper: robust test, p ≤ 0.01\n(data-driven reference)',
                 selection = 'First fit: 73 selected HK genes\n(95% CI excludes 0)',
                 all_hk = 'First fit: all 281 HK genes\n(95% CI excludes 0)')
    coefs <- bind_rows(paper, read_first_fit('selection'), read_first_fit('all_hk'))
    title_mode <- 'first fits (95% CI excludes 0) next to the paper'
}
coefs <- coefs %>%
    select(column, category, estimate, lower, upper, direction, counts) %>%
    mutate(column = factor(columns[column], levels = columns))

theme_panel <- function() {
    theme_minimal(base_size = 9) +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
              panel.grid.minor = element_blank(), strip.text = element_text(size = 8.5, face = 'bold'))
}

# Top: all genes that count ---------------------------------------------------------------------
counting <- filter(coefs, counts)
n_labels <- counting %>%
    count(column, direction) %>%
    mutate(x = ifelse(direction == 'down', -2.6, 2.6), y = 1.45, label = str_c(n, ' ', direction))

plot_all <- ggplot(counting) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_pointrange(aes(y = 1, x = estimate, xmin = lower, xmax = upper),
                    position = position_jitter(height = 0.1, seed = 1),
                    alpha = 0.04, size = 0.15, linewidth = 0.25) +
    geom_violin(aes(y = 1, x = estimate, group = direction), fill = beige, linewidth = 0.3,
                quantile.linetype = 3, position = 'identity') +
    geom_text(data = n_labels, aes(x = x, y = y, label = label), size = 3.6, color = '#444444') +
    facet_wrap(~column, nrow = 1) +
    coord_cartesian(xlim = c(-3.5, 3.5)) +
    labs(x = 'Log2 fold change over iron presence (high / low)') +
    theme_panel()

# Bottom: specific cases ------------------------------------------------------------------------
cases <- coefs %>%
    left_join(annotation, by = 'category') %>%
    mutate(label = case_when(str_detect(PFAMs, 'Flavodoxin_1') ~ 'Flavodoxin',
                             str_detect(PFAMs, 'Fer2') ~ 'Ferredoxin',
                             str_detect(category, isip2a) ~ 'ISIP2A')) %>%
    filter(!is.na(label)) %>%
    group_by(category) %>%
    filter(any(counts) | str_detect(category, flavodoxin_ns)) %>%
    ungroup() %>%
    mutate(gene = str_extract(category, 'CAH\\d+'),
           color = coalesce(highlight[gene], 'black'))
rows <- cases %>%
    distinct(category, label, gene) %>%
    arrange(label, gene) %>%
    mutate(row = row_number(), row_label = str_c(label, ' ', gene))
cases <- left_join(cases, select(rows, category, row, row_label), by = 'category')

plot_cases <- ggplot(cases, aes(y = row, color = I(color))) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_pointrange(aes(x = estimate, xmin = lower, xmax = upper, alpha = counts), size = 0.2, linewidth = 0.5) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = 'none') +
    scale_y_continuous(breaks = rows$row, labels = rows$row_label, trans = 'reverse') +
    facet_wrap(~column, nrow = 1) +
    coord_cartesian(xlim = c(-3.5, 3.5)) +
    labs(x = 'Log2 fold change over iron presence (high / low)',
         title = 'Specific cases (faded: not counted in that column; coloured: the flavodoxins highlighted in Fig 2)') +
    theme_panel() +
    theme(axis.text.y = element_text(size = 7), strip.text = element_blank(),
          plot.title = element_text(size = 8.5))

fig <- (plot_all / plot_cases) +
    plot_layout(heights = c(0.55, 0.45)) +
    plot_annotation(title = str_c('Pelagomonas calceolata regulation with iron: fastEmu reference sets compared (', title_mode, ')'),
                    subtitle = 'Same genes and samples. First fit: which HK genes anchor zero and pick the genes to test (interval excludes 0, |estimate| ≥ 0.4; the paper tested all with the interval excluding 0). Test: the reference of the score test.',
                    theme = theme(plot.title = element_text(size = 12, face = 'bold'),
                                  plot.subtitle = element_text(size = 9, color = 'grey30')))

save_figure(fig, str_c('fig2_fastemu_reference_comparison_', mode), width = if (mode == 'robust') 18 else 12, height = 8.5)
