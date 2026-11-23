# Supplementary Figure S3: housekeeping (HK) genes as a proxy of the abundance of Bathycoccus prasinos
# (the Fig 4 analysis, on MATOU instead of a transcriptome quantification).
#
# Expression: MATOU metaT of the B. prasinos unigenes (blast hits of its transcriptome), summarised per
# transcript (median over its unigenes) and sample. Abundance (Y): median MATOU metaG over its unigenes
# per sample, in samples where >= 100 unigenes are expressed.
#   A: adjusted R² of Y ~ geometric mean of the expression of 20 genes, for 1000 draws from the
#      moderately expressed non-HK genes (25th-90th percentile of mean expression) and from the HK genes;
#      red line: all HK genes
#   B: example regressions (draw 5 of each set) and all HK genes; values in units of 1e-6
#
# Inputs (data/raw/quantification/statistics, made by scripts/00_data_generation/supplementary/S3/):
#   bathycoccus_prasinos_blast_92id.tsv, bathycoccus_prasinos_matou_{metaT,metaG}.tsv
# Original: delving_into_florenciella src/eda/comparison_housekeeping_allgenes_bathy_matou.R (redrawn
#   2026-11-18: adjusted R² in the histograms, legend). Its exhaustive / top-N search settings were unused.
# Run time: ~1 min.
# The MATOU tables come from an earlier blast run than bathycoccus_prasinos_blast_92id.tsv (43 of their 5106
# unigenes are not hits in it now; see scripts/00_data_generation/supplementary/S3/02_bathycoccus_matou_tables.R).
# R: the paper figure was drawn with R 4.5 / ggplot2 4.0.2 (/usr/local/bin/Rscript), which reproduces it
#    pixel for pixel; R 4.1 (/usr/bin/Rscript) gives the same numbers with small rendering differences.
# Output: results/figures/supplementary/S3_randsampling_housekeeping_bathy.{pdf,png}

source('R/load.R')
library(patchwork)

K <- 20               # genes per draw
N_DRAWS <- 1000
MIN_EXPRESSED <- 100  # unigenes expressed in a sample for it to be used
EXPR_QUANTILES <- c(0.25, 0.90)   # moderately expressed genes (values are relative abundances, not TPM)

fill_colors <- c(random = '#7EBDC2', housekeeping = '#231f20')
allhk_color <- '#bb4430'

set.seed(4223)

# --- data --------------------------------------------------------------------------------------------
blast <- read_tsv(p('quant_stats', 'bathycoccus_prasinos_blast_92id.tsv'), show_col_types = FALSE,
                  col_names = c('qseqid', 'sseqid', 'pident', 'length', 'qlen', 'qcovs', 'evalue', 'bitscore')) %>%
    mutate(sseqid = as.numeric(str_remove(sseqid, 'MATOU-v1.5.')))

# sample codes without the year (last two characters)
read_matou_bathy <- function(file) {
    read_tsv(p('quant_stats', file), show_col_types = FALSE) %>%
        mutate(samplename = str_sub(samplename, 0, -3))
}
metaT <- read_matou_bathy('bathycoccus_prasinos_matou_metaT.tsv')
metaG <- read_matou_bathy('bathycoccus_prasinos_matou_metaG.tsv')

# transcript x sample expression: median over the MATOU unigenes hit by the transcript
expr_table <- metaT %>%
    left_join(blast, by = c('geneid' = 'sseqid'), relationship = 'many-to-many') %>%
    filter(!is.na(qseqid)) %>%
    group_by(qseqid, samplename) %>%
    summarize(value = median(value), .groups = 'drop') %>%
    pivot_wider(names_from = samplename, values_from = value, values_fill = 0) %>%
    rename(Gene = qseqid)

annotation <- read_annotation_housekeeping(p('annotation', 'EP00893_Bathycoccus_prasinos.emapper.annotations'))
housekeeping_vars <- housekeeping_ids(annotation)

expr_no_houseg <- expr_table %>% filter(!Gene %in% housekeeping_vars)
mean_expr <- rowMeans(expr_no_houseg[, -1])
bounds <- quantile(mean_expr, probs = EXPR_QUANTILES)
genes_good_ab <- expr_no_houseg$Gene[mean_expr >= bounds[1] & mean_expr <= bounds[2]]

random_mat <- expression_to_sample_matrix(expr_table, genes_good_ab)

# Y: median metaG of the unigenes per sample
samples_expressed <- metaT %>%
    group_by(samplename) %>%
    summarize(n_expressed = sum(value > 0), .groups = 'drop') %>%
    filter(n_expressed >= MIN_EXPRESSED)
metaG_median <- metaG %>%
    filter(samplename %in% samples_expressed$samplename) %>%
    group_by(samplename) %>%
    summarize(median_abundance = median(value), .groups = 'drop')

samples_common <- intersect(rownames(random_mat), metaG_median$samplename)
random_mat <- random_mat[samples_common, ]
Y <- metaG_median$median_abundance[match(samples_common, metaG_median$samplename)]

hk_mat <- expression_to_sample_matrix(expr_table, housekeeping_vars)[samples_common, ]
n_hk <- ncol(hk_mat)
k_hk <- min(K, n_hk - 1)

# --- random draws ----------------------------------------------------------------------------------
atrandom     <- sample_geomean_subsets(random_mat, Y, k = K, n_draws = N_DRAWS)
housekeeping <- sample_geomean_subsets(hk_mat, Y, k = k_hk, n_draws = N_DRAWS)
all_hk <- fit_geomean_model(colnames(hk_mat), hk_mat, Y)

draws <- bind_rows(random = atrandom, housekeeping = housekeeping, .id = 'comparison') %>%
    mutate(comparison = factor(comparison, levels = names(fill_colors)))
draws %>%
    group_by(comparison) %>%
    summarize(n = n(), median_adj_r2 = median(adj_r2)) %>%
    mutate(all_hk_adj_r2 = all_hk$adj_r2) %>%
    print()

# --- figure ----------------------------------------------------------------------------------------
fill_labels <- c(random = str_glue('{N_DRAWS} draws of {K} genes\nfrom all genes'),
                 housekeeping = str_glue('{N_DRAWS} draws of {k_hk} genes\nfrom the {n_hk} housekeeping genes'))

plot_comparison <- draws %>%
    ggplot(aes(adj_r2)) +
    geom_histogram(aes(fill = comparison), bins = 50, alpha = 0.75, position = 'identity') +
    geom_vline(aes(xintercept = all_hk$adj_r2, colour = 'allhk'), linewidth = 1.2) +
    scale_colour_manual(values = c(allhk = allhk_color),
                        labels = c(allhk = str_glue('All {n_hk}\nhousekeeping genes')), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_manual(values = fill_colors, labels = fill_labels, name = NULL) +
    labs(y = 'Count', x = expression(paste('Adjusted ', R^{2}))) +
    guides(colour = guide_legend(order = 2), fill = guide_legend(order = 1)) +
    theme_bw() +
    theme(legend.position = 'inside', legend.position.inside = c(0.02, 0.98),
          legend.justification = c(0, 1), legend.background = element_rect(fill = alpha('white', 0.8), colour = NA),
          legend.key.spacing.y = unit(4, 'pt'), legend.spacing.y = unit(2, 'pt'),
          panel.grid.minor = element_blank(), plot.margin = margin(5, 15, 5, 5))

# metaG vs geometric mean of a gene set, in units of 1e-6
plot_relationship <- function(genes, gene_mat, colorline, title) {
    gm <- geom_mean_rows(gene_mat[, genes, drop = FALSE])
    adj <- fit_geomean_model(genes, gene_mat, Y)$adj_r2
    tibble(gm = gm * 1e6, metag = Y * 1e6) %>%
        ggplot(aes(gm, metag)) +
        geom_point(shape = 21, size = 1.2) +
        geom_smooth(method = 'lm', formula = y ~ x, color = colorline) +
        scale_x_continuous(n.breaks = 4) +
        scale_y_continuous(n.breaks = 4) +
        labs(title = title, subtitle = sprintf('adj. R² = %.2f', adj),
             y = expression(paste('metaG median abundance (', 10^{-6}, ')')),
             x = expression(paste('metaT geometric mean (', 10^{-6}, ')'))) +
        theme_classic() +
        theme(plot.title = element_text(size = 9.5, colour = colorline, face = 'bold'),
              plot.subtitle = element_text(size = 8.5), axis.text = element_text(size = 7.5))
}
example_genes <- function(draws, i = 5) str_split_1(draws$genes[i], ', ')

examples <- ((plot_relationship(example_genes(atrandom), random_mat, fill_colors[['random']], str_glue('{K} random genes')) +
                  theme(axis.title.x = element_blank()) +
              plot_relationship(example_genes(housekeeping), hk_mat, fill_colors[['housekeeping']],
                                str_glue('{k_hk} housekeeping genes')) +
                  theme(axis.title = element_blank())) /
             plot_relationship(colnames(hk_mat), hk_mat, allhk_color, str_glue('All {n_hk} housekeeping genes'))) +
    plot_layout(heights = c(0.4, 0.6), axis_titles = 'collect')

final <- (plot_comparison | examples) +
    plot_layout(widths = c(0.55, 0.45)) +
    plot_annotation(tag_levels = list(c('A', 'B', '', '')))

dir.create(p('figures', 'supplementary'), showWarnings = FALSE)
save_figure(final, 'supplementary/S3_randsampling_housekeeping_bathy', width = 12, height = 6.5)
