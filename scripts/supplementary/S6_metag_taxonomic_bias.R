# Supplementary Figure S6: the MATOU metagenome under-represents dinoflagellates.
#   A. log2 ratio between the share of each lineage in housekeeping (HK) gene metaG and its share in HK
#      metaT: median (dot) and IQR (bar) over the paired Tara samples where the lineage is > 1% of the HK
#      metaT (448 pairs; Metazoa excluded). 0: same share; -2: 4x under-represented in metaG.
#   B. fraction of the HK unigene x sample pairs with metaT that have any metaG, by metaT quintile;
#      Dinoflagellata, diatoms, Haptophyta and Pelagophyceae highlighted, other lineages in grey.
#
# Input:  data/intermediate/supplementary/S6/{hk_lineage_shares.tsv.gz, hk_detection.tsv}
#         (scripts/00_data_generation/supplementary/S6/01_metag_taxonomic_bias.R)
# Original: gene_explorer scripts/figures/metag_taxonomic_bias.R. There, a per size fraction heatmap and a
#   Flavodoxin / SIT / NRT2 detection panel were also built but left out of the figure; they are not
#   drawn here (their numbers are in the S6 tables).
# Run time: seconds.
# Output: results/figures/supplementary/S6_metag_taxonomic_bias.{pdf,png}

source('R/load.R')
library(patchwork)

indir <- p('intermediate', 'supplementary', 'S6')

# Palette: diverging red (under-represented in metaG) <-> blue (over), grey midpoint; fixed colours for
# the highlighted lineages; grey for context
ink <- '#0b0b0b'; ink2 <- '#52514e'; grid_col <- '#e4e3df'; context <- '#b9b8b2'
div_low <- '#d03b3b'; div_mid <- '#f0efec'; div_high <- '#2a78d6'
highlight <- c(Dinoflagellata = '#eb6834', Diatoms = '#2a78d6', Haptophyta = '#1baf7a', Pelagophyceae = '#6250d6')

theme_bias <- function() {
    theme_minimal(base_size = 9) +
        theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
              panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
              plot.title = element_text(face = 'bold', size = 10), plot.subtitle = element_text(colour = ink2, size = 8),
              plot.title.position = 'plot', strip.text = element_text(face = 'bold', hjust = 0),
              legend.title = element_text(size = 8), legend.text = element_text(size = 7))
}

log2_ratio <- function(g, t) log2((g + 1e-6) / t)

# A. lineage bias --------------------------------------------------------------------------------------
overall <- read_tsv(file.path(indir, 'hk_lineage_shares.tsv.gz'), show_col_types = FALSE) %>%
    filter(share_metaT > 0.01) %>%
    mutate(l2 = log2_ratio(share_metaG, share_metaT)) %>%
    group_by(lineage) %>%
    summarise(n = n(), med = median(l2), q25 = quantile(l2, .25), q75 = quantile(l2, .75)) %>%
    mutate(lineage = fct_reorder(lineage, med))

pA <- ggplot(overall, aes(med, lineage)) +
    geom_vline(xintercept = 0, colour = ink2, linewidth = 0.4) +
    geom_segment(aes(x = q25, xend = q75, yend = lineage), colour = context, linewidth = 2, lineend = 'round') +
    geom_point(aes(fill = med), shape = 21, size = 3, colour = 'white', stroke = 0.6) +
    geom_text(aes(x = q75 + 0.2, label = str_c('n=', n)), hjust = 0, size = 2.3, colour = ink2) +
    scale_fill_gradient2(low = div_low, mid = div_mid, high = div_high, midpoint = 0, limits = c(-3.3, 3.3),
                         oob = scales::squish, breaks = -3:3, name = 'log2 metaG share /\nmetaT share') +
    scale_x_continuous(breaks = -4:2, expand = expansion(mult = c(0.03, 0.12))) +
    labs(title = 'A. Lineage share in HK metaG vs HK metaT',
         subtitle = 'Median (dot) and IQR (bar) over paired samples where the lineage is > 1% of HK metaT\n0 = same share in metaG as in metaT; −2 = 4× under-represented in metaG',
         x = 'log2(share in metaG / share in metaT)', y = NULL) +
    theme_bias() + theme(legend.position = 'none', panel.grid.major.y = element_blank())

# B. HK detection by metaT quintile ---------------------------------------------------------------------
detq <- read_tsv(file.path(indir, 'hk_detection.tsv'), show_col_types = FALSE) %>%
    select(lineage, Q1:Q5) %>%
    pivot_longer(Q1:Q5, names_to = 'quintile', values_to = 'detected') %>%
    mutate(q = as.integer(str_remove(quintile, 'Q')),
           group = if_else(lineage %in% names(highlight), lineage, 'Other lineages'))
labels_B <- detq %>% filter(q == 5, group != 'Other lineages')

pB <- ggplot(detq, aes(q, detected, group = lineage)) +
    geom_line(data = ~filter(.x, group == 'Other lineages'), colour = context, linewidth = 0.5) +
    geom_line(data = ~filter(.x, group != 'Other lineages'), aes(colour = group), linewidth = 1) +
    geom_point(data = ~filter(.x, group != 'Other lineages'), aes(colour = group), size = 2.2) +
    geom_text(data = labels_B, aes(label = lineage), hjust = -0.15, size = 2.6, colour = ink) +
    scale_colour_manual(values = highlight, name = NULL) +
    scale_x_continuous(breaks = 1:5, labels = c('Q1\nlowest', 'Q2', 'Q3', 'Q4', 'Q5\nhighest'),
                       expand = expansion(add = c(0.15, 1.6))) +
    coord_cartesian(clip = 'off') +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(title = 'B. HK unigenes with any metaG, by metaT level',
         subtitle = 'Unigene × sample pairs with metaT; quintiles of metaT over all lineages\ngrey: the other lineages',
         x = 'metaT quintile', y = 'with metaG') +
    theme_bias() + theme(legend.position = 'bottom', panel.grid.major.x = element_blank())

fig <- pA + pB + plot_layout(widths = c(1.15, 1)) +
    plot_annotation(title = 'The MATOU metaG under-represents dinoflagellates',
                    subtitle = 'Housekeeping (HK) genes in 448 paired metaT/metaG Tara samples; Metazoa excluded',
                    theme = theme(plot.title = element_text(face = 'bold', size = 13, colour = ink),
                                  plot.subtitle = element_text(size = 9, colour = ink2)))

dir.create(p('figures', 'supplementary'), showWarnings = FALSE)
save_figure(fig, 'supplementary/S6_metag_taxonomic_bias', width = 12, height = 6, dpi = 200, pdf_device = cairo_pdf)
