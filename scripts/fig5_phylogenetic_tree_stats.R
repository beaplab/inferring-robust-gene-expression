# Figure 5: statistics of the reference trees of the housekeeping (HK) genes and of the MATOU metaT
# placements on them, with the three target genes (Flavodoxin, NRT2, SIT) for comparison.
#   A: number of tips, treeness, median multiplicity per sample and median LWR of every gene
#      (violin = HK genes coloured by functional category; red = target genes; labels = 3 extremes)
#   B: toy diagram explaining the likelihood weight ratio (LWR) and the multiplicity
#   C: placements of one sample (128DCM1GGZZ14) on the trees of rpl5 and RICTOR
#      (bubbles: best placement of each pquery summed per edge; size = multiplicity, colour = LWR)
#
# Inputs:
#   data/intermediate/fig5/gene_stats_summary.tsv   scripts/00_data_generation/fig5/01_tree_and_placement_stats.R
#   data/raw/placement/hk_matou_jplaces/<gene>/128DCM1GGZZ14filter_basal_branches.jplace.gz
# Original: gene_explorer scripts/figures/phylogenetic_tree_examples_matou.R (figure 10)
# Output: results/figures/fig5_phylogenetic_tree_stats.{pdf,png}

source('R/load.R')
suppressPackageStartupMessages({
    library(ggtree)
    library(patchwork)
    library(ggrepel)
})

example_genes <- c('rpl5', 'RICTOR')
example_sample <- '128DCM1GGZZ14filter'

gene_stats <- read_tsv(p('intermediate', 'fig5', 'gene_stats_summary.tsv'), show_col_types = FALSE)

# LWR is a magnitude between 0 and 1: single hue, light to dark
lwr_scale <- scale_colour_gradient(low = '#c6dbef', high = '#08306b', limits = c(0, 1), name = 'LWR')

a4_theme <- theme_minimal(base_size = 8) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = 'bold', size = 9),
          plot.tag = element_text(face = 'bold', size = 11))

# A: gene statistics -----------------------------------------------------------

target_col <- '#e34948'
gene_stats_cat <- gene_stats %>%
    mutate(label = case_when(gene_group == 'target' ~ 'Target gene',
                             is.na(label) ~ 'Other',
                             TRUE ~ label),
           # top functional categories keep a hue, the rest go to Other
           label_short = fct_lump_n(label, n = 6, other_level = 'Other') %>%
               as.character() %>%
               if_else(gene_group == 'target', 'Target gene', .) %>%
               fct_infreq() %>%
               fct_relevel('Other', 'Target gene', after = Inf))

top_categories <- setdiff(levels(gene_stats_cat$label_short), c('Other', 'Target gene'))
category_colors <- c('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#6250d6')[seq_along(top_categories)] %>%
    set_names(top_categories) %>%
    c('Other' = 'grey70', 'Target gene' = target_col)

stat_labels <- c(ntips = 'Number of tips',
                 treeness = 'Treeness',
                 median_multiplicity_per_sample = 'Median multiplicity',
                 lwr_median = 'Median LWR')
log_stats <- c('median_multiplicity_per_sample')

targets_df <- filter(gene_stats_cat, gene_group == 'target')
# fixed jitter so that points and labels of the extreme genes match
set.seed(42)
housekeeping_df <- filter(gene_stats_cat, gene_group == 'housekeeping') %>%
    mutate(x_jit = 1 + runif(n(), -0.15, 0.15))

draw_violins <- function(selcol) {
    extremes_df <- housekeeping_df %>%
        filter(!is.na(.data[[selcol]])) %>%
        mutate(rank = rank(.data[[selcol]], ties.method = 'first')) %>%
        filter(rank <= 3 | rank > n() - 3)
    p <- ggplot(housekeeping_df, aes(y = .data[[selcol]])) +
        geom_violin(aes(x = 1), colour = 'grey60', fill = 'grey95', width = 0.8, linewidth = 0.3) +
        geom_point(aes(x = x_jit, colour = label_short), size = 0.5, alpha = 0.7) +
        # median as a short horizontal bar
        stat_summary(aes(x = 0.9, xend = 1.1, yend = after_stat(y)), fun = median,
                     geom = 'segment', linewidth = 0.6, colour = 'black') +
        geom_text_repel(data = extremes_df, aes(x = x_jit, label = houseg), colour = 'grey30',
                        size = 1.9, nudge_x = -0.35, direction = 'y', hjust = 1,
                        box.padding = 0.3, force = 3, segment.size = 0.2,
                        min.segment.length = 0, max.overlaps = Inf, seed = 42) +
        # target genes in their own column, close to the housekeeping violin
        geom_point(data = targets_df, aes(x = 1.85, colour = label_short), size = 1.3) +
        geom_text_repel(data = targets_df, aes(x = 1.85, label = houseg), colour = target_col,
                        size = 2.1, nudge_x = 0.25, direction = 'y', hjust = 0,
                        segment.size = 0.2, min.segment.length = 0, seed = 42) +
        # targets keep their colour but are not listed in the legend
        scale_colour_manual(values = category_colors, breaks = c(top_categories, 'Other'),
                            name = 'Housekeeping genes functional category') +
        scale_x_continuous(breaks = NULL, limits = c(0.35, 2.75)) +
        labs(title = stat_labels[[selcol]], x = NULL, y = NULL) +
        a4_theme +
        theme(panel.grid.major.x = element_blank(),
              plot.title = element_text(face = 'plain', size = 9))
    # extra room on top so the labels of the maximum genes do not overlap
    y_expand <- expansion(mult = c(0.03, 0.12))
    p + if (selcol %in% log_stats) scale_y_log10(expand = y_expand) else scale_y_continuous(expand = y_expand)
}

violin_plots <- map(names(stat_labels), draw_violins)
violin_plots[[1]] <- violin_plots[[1]] + labs(tag = 'A')
violins <- wrap_plots(violin_plots, nrow = 1) +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom',
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7.5, face = 'bold'),
          legend.title.position = 'top',
          legend.key.size = unit(3, 'mm')) &
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2)))
violins <- violins +
    plot_annotation(title = 'Statistics of multiple phylogenetic trees',
                    theme = theme(plot.title = element_text(face = 'bold', size = 10, hjust = 0.5)))

# B: LWR and multiplicity diagram -------------------------------------------------
# Two sequence variants placed on a toy tree: colour = sequence, bubble size = multiplicity of the
# sequence, opacity and label = LWR of each candidate edge (the LWRs of one sequence sum to 1)

toy_tree <- ape::read.tree(text = '((((A:1,B:1):1,C:2):1,(D:1.5,E:1.5):1.5):1,((F:1,G:1):1.5,H:2.5):1.5);')
toy_p <- ggtree(toy_tree, colour = 'grey40', linewidth = 0.6)
toy_edges <- toy_p$data %>%
    filter(node != parent) %>%
    mutate(x_mid = x - branch.length / 2) %>%
    select(node, x_mid, y)
clade_node <- ape::getMRCA(toy_tree, c('A', 'B', 'C'))
node_of <- function(tips) if (length(tips) == 1) which(toy_tree$tip.label == tips) else ape::getMRCA(toy_tree, tips)

toy_placements <- bind_rows(
    # spread over the edges of a clade
    tibble(seq_id = 'seq1',
           node = c(node_of('A'), node_of('B'), node_of(c('A', 'B')), node_of('C')),
           lwr = c(0.45, 0.30, 0.15, 0.10),
           multiplicity = 3),
    # almost all on a single edge
    tibble(seq_id = 'seq2',
           node = c(node_of('F'), node_of('G')),
           lwr = c(0.95, 0.05),
           multiplicity = 9)) %>%
    left_join(toy_edges, by = 'node')

# hues outside the blue LWR scale of panel C, one per sequence, so they are not read as tree placements
box_colors <- c(seq1 = '#e66101', seq2 = '#5e3c99')
seq_notes <- toy_placements %>%
    group_by(seq_id) %>%
    summarize(y = mean(y)) %>%
    mutate(label = c(seq1 = 'Sequence 1\nmultiplicity 3\nlow LWR',
                     seq2 = 'Sequence 2\nmultiplicity 9\nhigh LWR')[seq_id])

explanation_box <- toy_p +
    geom_hilight(node = clade_node, fill = 'grey90', extend = 0.5) +
    geom_tree(colour = 'grey40', linewidth = 0.4) +
    geom_tiplab(size = 2.4, offset = 0.08) +
    geom_point(data = toy_placements, aes(x = x_mid, y = y, size = multiplicity,
                                          colour = seq_id, alpha = lwr)) +
    # LWR to the left of each bubble, further away for the bigger bubbles
    geom_label(data = toy_placements, aes(x = x_mid - 0.12 - 0.022 * multiplicity, y = y, label = lwr),
               hjust = 1, size = 2.1, colour = 'grey20', fill = 'grey98',
               label.size = 0, label.padding = unit(0.08, 'lines')) +
    geom_text(data = seq_notes, aes(x = 4.6, y = y, label = label, colour = seq_id),
              hjust = 0, size = 2.1, lineheight = 0.95, fontface = 'bold') +
    scale_colour_manual(values = box_colors) +
    scale_alpha_continuous(range = c(0.35, 1), limits = c(0, 1)) +
    scale_size_area(max_size = 7, limits = c(0, 10)) +
    coord_cartesian(xlim = c(0, 6.6), clip = 'off') +
    labs(title = 'Explanation of likelihood\nweight ratio (LWR)\nand multiplicity', tag = 'B') +
    a4_theme +
    theme(legend.position = 'none',
          plot.title = element_text(face = 'bold', size = 8),
          panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
          plot.background = element_rect(fill = 'grey98', colour = 'grey60', linewidth = 0.4),
          panel.background = element_rect(fill = 'grey98', colour = NA),
          plot.margin = margin(5, 5, 5, 5))

# C: placements of one sample on two trees ------------------------------------------

examples <- map(example_genes, ~read_placements(
    p('placement', 'hk_matou_jplaces', .x, str_c(example_sample, '_basal_branches.jplace.gz')))) %>%
    set_names(example_genes)

# best placement per pquery, aggregated by edge: bubble size = multiplicity,
# colour = multiplicity-weighted mean LWR of the pqueries placed there
edge_summary <- map(examples, ~.x$pq %>%
                        filter(best) %>%
                        distinct(pquery, .keep_all = TRUE) %>%
                        group_by(edge_num) %>%
                        summarize(lwr = weighted.mean(lwr, multiplicity),
                                  multiplicity = sum(multiplicity),
                                  n_pqueries = n()))
mult_limits <- range(map(edge_summary, 'multiplicity'))

# Midpoint of every edge in the ggtree layout, keyed by jplace edge number
edge_positions <- function(p) {
    p$data %>%
        mutate(edge_num = as.integer(str_extract(label, '(?<=_E)\\d+$')),
               x_mid = x - branch.length / 2) %>%
        filter(!is.na(edge_num), !is.na(x_mid)) %>%
        select(edge_num, x_mid, y)
}

example_tree <- function(gene) {
    p <- ggtree(examples[[gene]]$tree, layout = 'circular', colour = 'grey60', linewidth = 0.15)
    bubbles <- left_join(edge_summary[[gene]], edge_positions(p), by = 'edge_num') %>%
        arrange(multiplicity)
    p +
        geom_point(data = bubbles, aes(x = x_mid, y = y, size = multiplicity, colour = lwr), alpha = 0.85) +
        lwr_scale +
        scale_size(range = c(0.3, 4), limits = mult_limits) +
        # no radial padding, the tree reaches the edge of the panel
        scale_x_continuous(expand = c(0, 0)) +
        theme_void() +
        theme(legend.position = 'none')
}

# polar coordinates leave an empty ring around the tree that margins cannot remove,
# so the tree is drawn as an inset larger than its cell (crop) under the gene title
tree_cell <- function(gene, tag = NULL, crop = 0.14, lwr_legend = FALSE) {
    # invisible point so the cell can carry the LWR colourbar of the trees
    ggplot(tibble(x = 0, y = 0, lwr = 0.5), aes(x, y, colour = lwr)) +
        geom_point(alpha = 0) +
        lwr_scale +
        guides(colour = guide_colourbar(barwidth = unit(2, 'mm'), barheight = unit(15, 'mm'))) +
        labs(title = gene, tag = tag, x = NULL, y = NULL) +
        a4_theme +
        theme(panel.grid = element_blank(), axis.text = element_blank(),
              plot.title = element_text(face = 'plain', size = 9, hjust = 0.5),
              # tag pushed inwards so it does not touch the explanation box
              plot.tag.position = c(0.06, 0.99),
              legend.position = if (lwr_legend) 'inside' else 'none',
              legend.position.inside = c(1.02, 0.15),
              legend.title = element_text(size = 7),
              legend.text = element_text(size = 6),
              plot.margin = margin(0, 0, 0, 0)) +
        inset_element(example_tree(gene), left = -crop, right = 1 + crop,
                      bottom = -crop, top = 1 + crop * 0.6, clip = FALSE)
}

# Figure (A4) ----------------------------------------------------------------------

# both example trees form a single panel (wrapped so the tree panels are not aligned to the panel of
# the explanation box, which shrank them)
trees <- wrap_elements(full = (tree_cell('rpl5', 'C') | tree_cell('RICTOR', lwr_legend = TRUE)) +
                           plot_annotation(title = 'Phylogenetic placements over one sample',
                                           theme = theme(plot.title = element_text(face = 'bold', size = 10, hjust = 0.5))))
bottom_row <- explanation_box + trees + plot_layout(widths = c(0.9, 2.5))

figure <- wrap_elements(violins) / wrap_elements(bottom_row) + plot_layout(heights = c(1.6, 1.05))

save_figure(figure, 'fig5_phylogenetic_tree_stats', width = 7.5, height = 9)
