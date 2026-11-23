# Figure 2: regulation of iron-related genes in Pelagomonas calceolata (left) and of the Flavodoxin
# variants across the reference tree (right), linked by a callout.
#
# Left: genes of P. calceolata differentially expressed between low and high iron (fastEmu: first fit
#   anchored on all HK genes, robust score test against fastEmu's data-driven reference set), with the
#   flavodoxins, ferredoxins and ISIP2A as specific cases.
# Right: Flavodoxin reference tree with the Graff et al. 2023 clades and per-clade violins of the
#   summed metaT/metaG ratio in low/high iron samples; the P. calceolata flavodoxins of the left panel
#   are marked on the tree in the same colours, and the first factor of the placement factorization of
#   the ratios is outlined over the clade strip.
#
# Inputs (data/intermediate/fig2, made by scripts/00_data_generation/fig2/):
#   pelagomonas_fastemu_iron_coeffs.tsv          01_pelagomonas_fastemu_iron.R
#   flavodoxin_ratio_factorization/              04_flavodoxin_ratio_factorization.R
#   flavodoxin_graff2023_clade_ratios.tsv        05_flavodoxin_graff2023_clades.R
# Original: gene_explorer scripts/figures/figure1_flavodoxin.R
# Sized for A4 landscape: saved at its print size (270 x 153 mm), so font sizes are in print points.
# Output: results/figures/fig2_species_vs_gene_DE.{pdf,png}

source('R/load.R')
suppressPackageStartupMessages({
    library(ggtree)
    library(patchwork)
    library(cowplot)
})

# Settings --------------------------------------------------------------------

# Pelagomonas flavodoxins highlighted in both panels
highlight <- tibble(gene = c('CAH0373297', 'CAH0367870', 'CAH0367856'),
                    color = c('#3887d7', '#e31c1c', '#2eb019'))

# Genes with a fastEmu p-value <= this are shown in the left panel
pval_cutoff <- 0.01

# ISIP2A genes of P. calceolata
isip2a <- 'CAH0367801.1|CAH0371101.1|CAH0367865.1|CAH0376816.1|CAH0367864.1'
# Flavodoxin highlighted on the tree, shown as a specific case even if not significant, with its p-value
flavodoxin_ns <- '9618'

# First factors of the placement factorization (metaT/metaG ratio, iron category) outlined on the
# tree over the clade strip; 0 for none. Only the first: F2 is a 2-tip factor and F3 changes with
# the pseudo-count
n_factors <- 1
# dim outline colours, one per factor, distinct from the clade colours of the strip
factor_colors <- c('#8a6f4d', '#4d7f8a', '#b0605a')

beige <- '#efe9df'
frame_color <- '#c5c5c5'
arrow_color <- '#3d311e'

# Text sizes in pt, at print size
sizes <- list(heading = 14, title = 11, base = 9, small = 7.5, counts = 13, callout = 8.5)

# Positions on the canvas (0-1, from the bottom left): x, y, width, height
layout <- list(
    left = c(0.005, 0.01, 0.4515, 0.87),
    right = c(0.4615, 0.01, 0.5335, 0.87),
    # over the empty right end of the left panel, straddling the gap
    callout = c(0.3885, 0.24, 0.12, 0.24),
    heading_y = 0.94
)

canvas <- list(width = 270, height = 153)  # mm

fig2_dir <- p('intermediate', 'fig2')
factorization_dir <- file.path(fig2_dir, 'flavodoxin_ratio_factorization')

# Flavodoxin panel: Graff et al. 2023 clades drawn on the tree of the ratio factorization
graff_labels <- c(clade1_dielcycle = 'Clade 1 (diel cycle)',
                  clade2_ironlimiting = 'Clade 2 (iron limiting)',
                  greenalgae_dinos = 'Green algae + dinoflagellates')
flavodoxin <- list(
    gene = 'Flavodoxin', variable = 'iron', ylab = 'metaT/metaG ratio',
    tree = file.path(factorization_dir, 'flavo-categorical-factors_tree.newick'),
    clades = p('gene_data', 'graff2023_flavodoxin_clades.tsv'),
    counts = file.path(fig2_dir, 'flavodoxin_graff2023_clade_ratios.tsv'),
    strip_title = 'Graff et al.\n2023 clades',
    facet_labels = graff_labels,
    clade_colors = c(clade1_dielcycle = '#ffe599', clade2_ironlimiting = '#b4a7d6',
                     greenalgae_dinos = '#b6d7a8'),
    exclude_samples = character(),
    # short clade names on the strip (the facets carry the full ones), just enough room for them
    strip_labels = function(clade) clade_short[clade],
    strip_vertical = TRUE,
    strip_room = 0.16)

clade_short <- c(clade1_dielcycle = 'Clade 1', clade2_ironlimiting = 'Clade 2',
                 greenalgae_dinos = 'Green algae\n+ dinos', basal = 'Basal')

# Pelagomonas data ------------------------------------------------------------

# fastEmu coefficients of the iron category per gene, with their PFAM annotation. fastEmu estimates
# low vs high iron in natural log, so the sign is flipped (positive = more expressed with high iron)
# and the values turned into log2.
annotation <- read_tsv(p('annotation', 'Pelagomonas_calceolata.emapper.annotations'),
                       comment = '##', show_col_types = FALSE) %>%
    select(category = `#query`, PFAMs)

diffexp <- read_tsv(file.path(fig2_dir, 'pelagomonas_fastemu_iron_coeffs.tsv'), show_col_types = FALSE) %>%
    left_join(annotation, by = 'category') %>%
    mutate(across(c(estimate, lower, upper), ~ log2(exp(-.x))),
           relationship_nutrients = ifelse(estimate > 0, 'upregulated', 'downregulated')) %>%
    select(category, estimate, lower, upper, pval, PFAMs, relationship_nutrients)
write_tsv(diffexp, p('tables', 'fig2_pelagomonas_iron_diffexp.tsv'))

significant <- filter(diffexp, !is.na(pval), pval <= pval_cutoff)

# Left panel: Pelagomonas -----------------------------------------------------

theme_left <- function() {
    theme_minimal(base_size = sizes$base) +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
              axis.text.x = element_text(size = sizes$small),
              plot.title = element_text(size = sizes$base))
}

# All significant genes: their estimate and interval, and a violin per direction
counts <- count(significant, relationship_nutrients)
plot_updown <- significant %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_pointrange(aes(y = 1, x = estimate, xmin = lower, xmax = upper),
                    position = position_jitter(height = 0.1, seed = 1),
                    alpha = 0.05, size = 0.2, linewidth = 0.3) +
    geom_violin(aes(y = 1, x = estimate, group = relationship_nutrients),
                fill = beige, linewidth = 0.3, quantile.linetype = 3, position = 'identity') +
    annotate('text', y = c(1.5, 1.4), x = c(-0.5, 1), size = sizes$counts / .pt, color = '#444444',
             label = counts$n[match(c('downregulated', 'upregulated'), counts$relationship_nutrients)]) +
    coord_cartesian(xlim = c(-3.5, 3.5)) +
    xlab('Log2 fold change over iron presence') +
    theme_left()

# Clade of each highlighted gene in the Flavodoxin tree (basal: outside the Graff et al. clades),
# shown in both panels with the clade colour of the tree strip and the violins
clade_colors <- clade_palette(read_tip_clades(flavodoxin)$clade, flavodoxin)
highlight <- highlight %>%
    left_join(read_tip_clades(flavodoxin) %>% transmute(gene = str_extract(label, 'CAH\\d+'), clade),
              by = 'gene') %>%
    mutate(clade = replace_na(clade, 'basal'), clade_color = clade_colors[clade])

# Specific cases: flavodoxins, ferredoxins and ISIP2A, plus the non-significant flavodoxin
selected <- significant %>%
    bind_rows(filter(diffexp, str_detect(category, flavodoxin_ns), !category %in% significant$category)) %>%
    mutate(label = case_when(str_detect(PFAMs, 'Flavodoxin_1') ~ 'Flavodoxin',
                             str_detect(PFAMs, 'Fer2') ~ 'Ferredoxin',
                             str_detect(category, isip2a) ~ 'ISIP2A')) %>%
    filter(!is.na(label)) %>%
    arrange(label) %>%
    mutate(label = ifelse(str_detect(category, flavodoxin_ns) & pval > pval_cutoff,
                          str_c('Flavodoxin (p=', signif(pval, 2), ')'), label),
           row = row_number(),
           gene = str_extract(category, 'CAH\\d+')) %>%
    left_join(highlight, by = 'gene') %>%
    mutate(color = replace_na(color, 'black'),
           # highlighted genes also name their tree clade
           label = ifelse(is.na(clade), label, str_c(label, ' · ', clade_short[clade])))

plot_cases <- selected %>%
    ggplot(aes(y = row, color = color)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_pointrange(aes(x = estimate, xmin = lower, xmax = upper), size = 0.25, linewidth = 0.5) +
    geom_text(data = filter(selected, is.na(clade)), aes(label = label, x = estimate),
              nudge_y = 0.45, size = sizes$small / .pt, color = 'black') +
    # highlighted genes: label on a chip of their clade colour
    geom_label(data = filter(selected, !is.na(clade)), aes(label = label, x = estimate, fill = I(clade_color)),
               nudge_y = 0.45, size = sizes$small / .pt, color = 'black', border.colour = NA,
               label.padding = unit(0.5, 'mm'), label.r = unit(0.5, 'mm')) +
    scale_color_identity() +
    coord_cartesian(xlim = c(-3.5, 3.5)) +
    xlab('Log2 fold change over iron presence') +
    ggtitle('Specific cases') +
    theme_left()

panel_left <- (plot_updown / plot_cases) +
    plot_layout(axes = 'collect', heights = c(0.68, 0.42)) +
    plot_annotation(title = bquote(italic('Pelagomonas calceolata') ~ 'regulation with iron concentration'),
                    theme = theme(plot.title = element_text(size = sizes$title, hjust = 0.5)))

# Right panel: Flavodoxin -----------------------------------------------------

tree <- plot_tree(flavodoxin, text_size = sizes$small)
# Branches in a single colour, without taxonomy, so the clade colours link the tree to the violins
branch_color <- 'grey30'
branch_colors <- setNames(rep(branch_color, length(group_colors) + 1), c(names(group_colors), 'mixed'))
tree_xmax <- max(tree$data$x)
n_tips <- sum(tree$data$isTip)
# Pelagomonas genes of the left panel, which are tips of the reference tree
genes_on_tree <- tree$data %>%
    filter(isTip) %>%
    mutate(gene = str_extract(label, 'CAH\\d+')) %>%
    inner_join(select(highlight, gene, color), by = 'gene')
# Each factor's clade (the side away from the root), as one run of consecutive tips
factor_boxes <- read_tsv(file.path(factorization_dir, 'flavo-categorical-factor_taxa.csv'), show_col_types = FALSE) %>%
    filter(Factor <= n_factors, RootSide == 0) %>%
    inner_join(filter(tree$data, isTip), by = c(Taxon = 'label')) %>%
    group_by(Factor) %>%
    summarize(ymin = min(y) - 0.5, ymax = max(y) + 0.5) %>%
    # fitted to the strip (tree_xmax * 1.02-1.06), the first factor just outside the later ones
    mutate(xmin = tree_xmax * ifelse(Factor == 1, 1.012, 1.018),
           xmax = tree_xmax * ifelse(Factor == 1, 1.068, 1.062),
           # the first label spelled out, so the F2, F3... below read as factors too
           label = ifelse(Factor == 1, 'Factor 1', str_c('F', Factor)),
           color = factor_colors[Factor])
tree <- tree +
    geom_rect(data = factor_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, color = factor_boxes$color, linewidth = 0.45) +
    geom_label(data = factor_boxes, aes(x = xmin, y = ymax, label = label),
               inherit.aes = FALSE, hjust = 1, vjust = 1, size = (sizes$small - 0.5) / .pt,
               fontface = 'bold', color = factor_boxes$color, fill = 'white',
               border.colour = factor_boxes$color, linewidth = 0.3,
               label.padding = unit(0.4, 'mm'), label.r = unit(0, 'mm')) +
    geom_point(data = genes_on_tree, aes(x = x, y = y, fill = I(color)), inherit.aes = FALSE,
               shape = 21, size = 2.2, stroke = 0.3, color = 'white') +
    # replaces the taxonomy scale; the clade colours stay for the leaders to the strip
    scale_color_manual(values = c(branch_colors, clade_colors), guide = 'none') +
    guides(color = 'none') +
    # key for the dots, in the empty space at the top left of the tree (left of the shallow tips)
    annotate('point', x = tree_xmax * c(0.02, 0.07, 0.12), y = n_tips * 0.975, shape = 21, size = 2.2,
             stroke = 0.3, color = 'white', fill = highlight$color) +
    annotate('text', x = 0, y = n_tips * 0.955, hjust = 0, vjust = 1,
             size = sizes$small / .pt, color = 'grey20', fontface = 'italic', label = 'P. calceolata') +
    annotate('text', x = 0, y = n_tips * 0.955 - n_tips * 0.03, hjust = 0, vjust = 1,
             lineheight = 0.95, size = sizes$small / .pt, color = 'grey20',
             label = 'flavodoxins\n(specific cases)')

# Facet titles: clade name, and a summary of the taxonomy only where it adds to the name
facet_titles <- c(clade1_dielcycle = 'Clade 1 (diel cycle)\nDiatomeae',
                  clade2_ironlimiting = 'Clade 2 (iron limiting)\nStramenopiles,\nAlveolata',
                  greenalgae_dinos = 'Green algae +\ndinoflagellates',
                  basal = 'Basal')

# log2 fold change high/low iron, the same direction as the left panel
violins <- plot_violins(flavodoxin, text_size = sizes$small - 0.5, point_size = 0.15, log2_fc = TRUE) +
    facet_wrap(~clade, nrow = 2, scales = 'free_y', labeller = as_labeller(facet_titles)) +
    # log tick labels a bit larger, so their superscript exponents stay readable
    theme(axis.text.x = element_text(size = sizes$small),
          axis.text.y = element_text(size = sizes$small + 1),
          axis.title = element_text(size = sizes$base),
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(2, 'mm'),
          # long facet titles may extend past their panel instead of being clipped
          strip.clip = 'off',
          plot.margin = margin(2, 2, 2, 0))

# free(): the tree takes the full height instead of aligning to the violin panels below their titles
panel_right <- (free(tree) | violins) +
    plot_layout(widths = c(1.15, 1)) +
    plot_annotation(title = 'Flavodoxin gene variants regulation with iron',
                    theme = theme(plot.title = element_text(size = sizes$title, hjust = 0.5)))

# Callout ---------------------------------------------------------------------

# Two questions linked by a loop of arrows
callout <- ggplot() +
    geom_curve(aes(x = 0.78, y = 0.64, xend = 0.22, yend = 0.56), curvature = 0.7,
               arrow = arrow(length = unit(1.5, 'mm'), type = 'closed'),
               linewidth = 0.8, color = arrow_color) +
    geom_curve(aes(x = 0.22, y = 0.36, xend = 0.78, yend = 0.44), curvature = 0.7,
               arrow = arrow(length = unit(1.5, 'mm'), type = 'closed'),
               linewidth = 0.8, color = arrow_color) +
    geom_label(aes(x = 0.5, y = c(0.92, 0.1),
                   label = c('Distribution\nother genes?',
                             'Flavodoxins present\nmultiple patterns.\nResponse in a\nlarger context?')),
               size = sizes$callout / .pt, fill = beige, border.colour = NA, lineheight = 0.95,
               label.padding = unit(1.2, 'mm'), label.r = unit(1, 'mm')) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = 'off') +
    theme_void()

# Figure ----------------------------------------------------------------------

draw_panel <- function(figure, plot, pos) {
    figure +
        draw_grob(grid::rectGrob(gp = grid::gpar(col = frame_color, fill = 'white', lwd = 0.8)),
                  pos[1], pos[2], pos[3], pos[4]) +
        draw_plot(plot, pos[1] + 0.005, pos[2] + 0.008, pos[3] - 0.01, pos[4] - 0.016)
}

figure <- ggdraw() %>%
    draw_panel(panel_left, layout$left) %>%
    draw_panel(panel_right, layout$right) +
    draw_label('Gene expression regulation in a species', fontface = 'bold', size = sizes$heading,
               x = layout$left[1] + layout$left[3] / 2, y = layout$heading_y) +
    draw_label('Functional gene analysis', fontface = 'bold', size = sizes$heading,
               x = layout$right[1] + layout$right[3] / 2, y = layout$heading_y) +
    draw_plot(callout, layout$callout[1], layout$callout[2], layout$callout[3], layout$callout[4]) +
    theme(plot.background = element_rect(fill = 'white', color = NA))

save_figure(figure, 'fig2_species_vs_gene_DE', width = canvas$width, height = canvas$height, units = 'mm')
