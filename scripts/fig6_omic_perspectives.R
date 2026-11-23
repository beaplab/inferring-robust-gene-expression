# Figure 6: taxonomic composition of the housekeeping (HK) genes compared with the other omic views
# of the same Tara samples (18S-V4 and V9 amplicons, MATOU metaG and metaT whole catalogues), per size
# fraction. Metazoa removed everywhere.
#   A: mean composition per size fraction (samples present in every dataset)
#   B: PCA (Hellinger) of the composition of every HK gene, fitted on the HK genes only; the other
#      datasets and the HK median are projected on it. Arrows: the 6 groups with the longest loadings.
#
# Inputs (data/intermediate/fig6, made by scripts/00_data_generation/fig6/):
#   taxon_compositions.parquet              03_taxonomy_compositions.R
#   housekeeping_gene_filter_matou.tsv      02_housekeeping_gene_filter.R (the genes of panel B)
# Original: gene_explorer scripts/figures/final_figure_hk_taxonomy_fractions.R (version "allfractions")
# Output: results/figures/fig6_omic_perspectives.{pdf,png}

source('R/load.R')
suppressPackageStartupMessages({
    library(arrow, warn.conflicts = FALSE)
    library(patchwork)
})

fractions <- c(GGMM = '0.8-5 µm', MMQQ = '5-20 µm', QQSS = '20-180 µm', SSUU = '180-2000 µm')

summary_sources <- c('18S-V4', '18S-V9', 'MATOU metaG (all genes)', 'MATOU metaT (all genes)')
# panel A uses the samples present in every summary dataset and in n_random_genes random HK genes
n_random_genes <- 10

fig_width <- 8.27
fig_height <- 7
ink <- '#0b0b0b'; ink2 <- '#52514e'; grid_col <- '#e4e3df'

theme_final <- function() {
    theme_minimal(base_size = 9) +
        theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2, size = 8),
              axis.title = element_text(size = 9),
              panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
              plot.tag = element_text(face = 'bold', size = 13), plot.tag.position = 'topleft',
              plot.title = element_text(face = 'bold', size = 10, lineheight = 1), plot.title.position = 'plot',
              strip.text = element_text(face = 'bold', size = 8.5, hjust = 0),
              legend.title = element_text(size = 8.5), legend.text = element_text(size = 8),
              legend.key.size = unit(3.5, 'mm'), legend.margin = margin(0, 0, 0, 0),
              legend.position = 'bottom', plot.margin = margin(12, 14, 4, 4))
}

# Data ----------------------------------------------------------------------------------------
# One row per source (dataset or HK gene), sample and taxon group: rel_value = relative abundance
compositions <- read_parquet(p('intermediate', 'fig6', 'taxon_compositions.parquet')) %>%
    mutate(code = str_extract(samplename, '[A-Z]{4}(?=\\d{2}$)')) %>%   # size fraction code
    filter(code %in% names(fractions)) %>%
    mutate(fraction = factor(fractions[code], levels = fractions))
gene_filter <- read_tsv(p('intermediate', 'fig6', 'housekeeping_gene_filter_matou.tsv'), show_col_types = FALSE)
hk_label <- str_subset(unique(compositions$source), '^Housekeeping genes')

# Samples present in every dataset (the summary datasets and 10 random HK genes)
hk_genes <- setdiff(unique(compositions$source), c(summary_sources, hk_label))
random_genes <- select_random_genes(sort(hk_genes), n_random_genes)
datasets <- c(summary_sources, random_genes)
complete_samples <- compositions %>%
    filter(source %in% datasets) %>%
    distinct(samplename, source) %>%
    count(samplename) %>%
    filter(n == length(datasets)) %>%
    pull(samplename)

ref_labels <- c('18S-V4' = '18S-V4', '18S-V9' = '18S-V9', 'MATOU metaG (all genes)' = 'metaG',
                'MATOU metaT (all genes)' = 'metaT')
ref_shapes <- c('18S-V4' = 24, '18S-V9' = 25, 'metaG' = 22, 'metaT' = 23, 'HK median' = 21)

# Mean composition per source, fraction and taxon (absent taxa count as 0)
mean_composition <- function(df) {
    df %>%
        select(source, fraction, samplename, taxon, rel_value) %>%
        complete(nesting(source, fraction, samplename), taxon, fill = list(rel_value = 0)) %>%
        group_by(source, fraction, taxon) %>%
        summarize(rel_value = mean(rel_value), .groups = 'drop')
}

# A. Mean composition per fraction --------------------------------------------------------------
pooled <- compositions %>%
    filter(samplename %in% complete_samples, source %in% c(summary_sources, hk_label)) %>%
    mean_composition() %>%
    mutate(taxon = factor(taxon, levels = c(taxon_groups, 'Other')),
           source = if_else(source == hk_label, 'HK median', ref_labels[source]),
           source = factor(source, levels = rev(c(ref_labels, 'HK median'))))

panel_a <- ggplot(pooled, aes(rel_value, source, fill = taxon)) +
    geom_col(colour = 'white', linewidth = 0.2, width = 0.8) +
    scale_x_continuous(labels = scales::percent, expand = c(0, 0), breaks = c(0, 0.5, 1)) +
    scale_fill_manual(values = taxa_colors, drop = FALSE, name = NULL,
                      guide = guide_legend(nrow = 2, byrow = TRUE)) +
    facet_wrap(~fraction, nrow = 1) +
    labs(tag = 'A', title = 'Mean composition per size fraction', x = 'Mean relative abundance', y = NULL) +
    theme_final() +
    theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(2.2, 'lines'))

# B. PCA of the HK genes per fraction -----------------------------------------------------------
comp <- compositions %>%
    filter(source %in% c(gene_filter$houseg, summary_sources, hk_label)) %>%
    mean_composition() %>%
    pivot_wider(names_from = taxon, values_from = rel_value, values_fill = 0)

mat <- sqrt(as.matrix(select(comp, -source, -fraction)))  # Hellinger
is_gene <- comp$source %in% gene_filter$houseg
pca <- prcomp(mat[is_gene, ], center = TRUE, scale. = FALSE)
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

scores <- bind_cols(select(comp, source, fraction), as_tibble(predict(pca, mat)[, 1:2])) %>%
    mutate(dominant = colnames(mat)[max.col(mat, ties.method = 'first')],
           dominant = factor(dominant, levels = names(taxa_colors)))
genes <- filter(scores, source %in% gene_filter$houseg)
refs <- filter(scores, !source %in% gene_filter$houseg) %>%
    mutate(source = if_else(source == hk_label, 'HK median', ref_labels[source]))

# loadings scaled to the spread of the genes
k <- 0.9 * max(abs(c(genes$PC1, genes$PC2))) / max(abs(pca$rotation[, 1:2]))
arrows <- as_tibble(pca$rotation[, 1:2], rownames = 'taxon') %>%
    mutate(len = sqrt(PC1^2 + PC2^2)) %>%
    slice_max(len, n = 6) %>%
    mutate(xend = k * PC1, yend = k * PC2, xlab = 1.12 * xend, ylab = 1.12 * yend)

panel_b <- ggplot(genes, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = 'grey90') +
    geom_vline(xintercept = 0, colour = 'grey90') +
    geom_segment(data = arrows, aes(x = 0, y = 0, xend = xend, yend = yend),
                 colour = 'grey55', linewidth = 0.3, arrow = arrow(length = unit(1.5, 'mm'))) +
    geom_point(aes(fill = dominant), shape = 21, size = 1.4, colour = 'white', stroke = 0.2, alpha = 0.9) +
    geom_label(data = arrows, aes(x = xlab, y = ylab, label = taxon),
               colour = 'grey30', fill = alpha('white', 0.8), label.size = 0,
               size = 2.6, fontface = 'italic', label.padding = unit(0.1, 'lines')) +
    geom_point(data = refs, aes(shape = source), fill = 'black', colour = 'white', size = 2.4, stroke = 0.4) +
    scale_fill_manual(values = taxa_colors, drop = TRUE, guide = 'none') +  # colours of A
    scale_shape_manual(values = ref_shapes, breaks = names(ref_shapes), name = NULL,
                       guide = guide_legend(nrow = 1, override.aes = list(fill = 'black', colour = 'black', size = 2.4))) +
    facet_wrap(~fraction, nrow = 1) +
    coord_cartesian(clip = 'off') +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.1))) +
    labs(tag = 'B', title = 'Taxonomic composition of the housekeeping genes (PCA), coloured by dominant group',
         x = str_c('PC1 (', var_expl[1], '%)'), y = str_c('PC2 (', var_expl[2], '%)')) +
    theme_final() +
    theme(panel.grid = element_blank(), panel.border = element_rect(colour = grid_col, fill = NA, linewidth = 0.4),
          panel.spacing.x = unit(0.6, 'lines'))

# Figure ---------------------------------------------------------------------------------------
figure <- panel_a / panel_b + plot_layout(heights = c(1, 1.15))
save_figure(figure, 'fig6_omic_perspectives', width = fig_width, height = fig_height, pdf_device = cairo_pdf)
