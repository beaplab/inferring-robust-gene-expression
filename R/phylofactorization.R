# Plots of the placement analyses of Flavodoxin (Fig 2), SIT and NRT2 (Fig 7): the reference tree
# with its branches coloured by taxonomy and the clades as a strip next to the tips, and per-clade
# violins of the summed counts/ratios between low and high values of an environmental variable.
#
# Original: gene_explorer scripts/figures/phylofactorization_plots.R + phylofactorization_panels.R
#
# A "panel" is a list describing one gene and set of clades:
#   gene, variable  gene name; environmental variable (a name of `env_variables`)
#   tree            newick of the reference tree
#   clades          two-column tsv (tip label, clade), no header
#   counts          tsv with columns clade, count, sample and the environmental variable
#   ylab            y axis title of the violins
#   strip_title     header of the clade strip next to the tree
#   strip_labels    function(clade) -> text on the strip
#   facet_labels    (optional) named vector clade -> name shown in the facets
#   clade_colors    (optional) named vector clade -> colour, instead of the default slots
#   strip_room      (optional) room right of the tree for the strip labels, as a fraction of its width
#   strip_vertical  (optional) TRUE to write the strip labels vertically, along the strip
#   strip_title_face (optional) font face of strip_title, e.g. "italic"; default plain
#   colored_labels  (optional) TRUE to write the strip labels (bold) and the facet titles in the clade colours
#   strip_guide     (optional) TRUE to draw a thin line behind the strip, from its title down to the bottom
#   y_break_step    (optional) powers of ten between y breaks of the violins (log1p axis); default 2
#   basal_composition (optional) FALSE to leave the taxonomic composition out of the basal facet title
#   xlab            (optional) x axis title of the violins, instead of "<Variable> category"
#   exclude_samples (optional) samples left out of the violins

# Environmental variables and their low/high cutoffs.
# iron_darwin_5m is in µmol/L (PANGAEA "Iron, total [µmol/l]"), so 5e-04 = 0.5 nmol/L
env_variables <- list(
    iron     = list(column = 'iron_darwin_5m', threshold = 5e-04),
    silicate = list(column = 'silicate_q2', threshold = 2),
    nitrate  = list(column = 'nitrate_nitrite_q2', threshold = 1)
)

# Panel of the clades found by a placement factorization (gappa analyze placement-factorization),
# as written by scripts/00_data_generation/fig7/03_ratio_placement_factorization.R:
# <prefix>factors_tree.newick, taxaall.tsv (clades) and the count table
factorization_panel <- function(dir, prefix, gene, variable, counts, ylab) {
    list(gene = gene, variable = variable, ylab = ylab,
         tree = file.path(dir, str_c(prefix, 'factors_tree.newick')),
         clades = file.path(dir, 'taxaall.tsv'),
         counts = file.path(dir, counts),
         strip_title = 'Placement\nfactorization\nclades',
         # clade_N -> "N" on the strip, facets keep the clade name
         strip_labels = function(clade) str_remove(clade, 'clade_'),
         exclude_samples = character())
}

# Palettes --------------------------------------------------------------------
# Categorical palettes checked for colour-vision-deficiency separation

# Clades, shared by the tree strip and the violins, assigned in order; basal = reads outside any clade
clade_slots <- c('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4')
# A panel can bring its own clade colours (panel$clade_colors) instead of the slots
clade_palette <- function(clades, panel) {
    clades <- sort(setdiff(clades, 'basal'))
    colors <- if (is.null(panel$clade_colors)) setNames(clade_slots[seq_along(clades)], clades)
              else panel$clade_colors[clades]
    c(colors, basal = '#a3a3a3')
}

# Taxonomy of the branches
group_colors <- c('Diatomeae' = '#1f7a3a', 'Other Stramenopiles' = '#2b9fc4', 'Alveolata' = '#b5532c',
                  'Archaeplastida' = '#7db546', 'Rhizaria' = '#6250d6', 'Haptophyta' = '#c48f00',
                  'Opisthokonta' = '#a63d6e', 'Other' = '#b8b8b8')
# Branches leading to more than one group
mixed_color <- 'grey45'
# Points of the violins
point_color <- 'grey20'

theme_phylofactorization <- function() {
    theme_bw() +
        theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
              strip.text = element_text(size = 8, lineheight = 0.95))
}

# Taxonomy --------------------------------------------------------------------

# EukProt v1 data sets table (supergroup and taxogroup of every EukProt ID), read once
eukprot_table <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            cache <<- read_tsv(p('gene_data', 'EukProt_included_data_sets.v01.2019_05_08.txt'),
                               show_col_types = FALSE) %>%
                select(EukProt_ID, Supergroup_UniEuk, Taxogroup_UniEuk)
        }
        cache
    }
})

# Tip taxonomy: supergroups, with diatoms split from the other stramenopiles and the
# plant and opisthokont lineages pooled
tip_taxonomy <- function(tree) {
    tibble(label = tree$tip.label, EukProt_ID = str_extract(label, '^EP\\d+')) %>%
        left_join(eukprot_table(), by = 'EukProt_ID') %>%
        mutate(group = case_when(
            Taxogroup_UniEuk == 'Diatomeae' ~ 'Diatomeae',
            Supergroup_UniEuk == 'Stramenopiles' ~ 'Other Stramenopiles',
            Supergroup_UniEuk %in% c('Chloroplastida', 'Rhodophyta', 'Glaucophyta') ~ 'Archaeplastida',
            Supergroup_UniEuk %in% c('Metazoa', 'Choanoflagellata', 'Filasterea', 'Ichthyosporea',
                                     'Pluriformea', 'Fungi', 'Rotosphaerida') ~ 'Opisthokonta',
            Supergroup_UniEuk %in% names(group_colors) ~ Supergroup_UniEuk,
            TRUE ~ 'Other'),
            group = factor(group, levels = names(group_colors)))
}

read_tip_clades <- function(panel) {
    read_tsv(panel$clades, col_names = c('label', 'clade'), show_col_types = FALSE)
}

# Taxonomic composition of each clade's reference tips (basal: tips outside every clade): the
# groups making up at least 10% of them, largest first, at most three
clade_composition <- function(panel) {
    tree <- ape::read.tree(panel$tree)
    tip_taxonomy(tree) %>%
        left_join(read_tip_clades(panel), by = 'label') %>%
        mutate(clade = replace_na(clade, 'basal')) %>%
        count(clade, group) %>%
        group_by(clade) %>%
        mutate(share = n / sum(n)) %>%
        filter(share >= 0.1) %>%
        slice_max(share, n = 3, with_ties = FALSE) %>%
        # non-breaking spaces keep each group name on one line when wrapped
        summarize(composition = str_c(str_replace_all(group, ' ', ' '), collapse = ', ')) %>%
        deframe()
}

# size of the clade name relative to the composition line in coloured facet titles
facet_name_scale <- 1.3

# Facet titles: clade name, and its taxonomic composition below in smaller text
# text_size: size in pt of the facet titles, used to size the composition line of coloured titles
clade_facet_labels <- function(panel, clades, text_size = NULL) {
    names_shown <- coalesce(unname(panel$facet_labels[clades]), clades)
    composition <- coalesce(clade_composition(panel)[clades], '')
    if (isFALSE(panel$basal_composition)) composition[clades == 'basal'] <- ''
    if (isTRUE(panel$colored_labels)) {
        # markdown for ggtext::element_markdown(): clade name bold in its colour, composition below;
        # narrow facets, so group names may break across lines here
        colors <- clade_palette(clades, panel)[clades]
        return(setNames(str_glue("<span style='color:{colors};font-size:{(text_size %||% 8) * facet_name_scale}pt'>",
                                 '**{names_shown}**</span><br>',
                                 "{str_replace_all(str_wrap(str_replace_all(composition, ' ', ' '), 15), '\n', '<br>')}"),
                        clades))
    }
    setNames(str_c(names_shown, '\n', str_wrap(composition, 30)), clades)
}

# Tree ------------------------------------------------------------------------

# Reference tree, branches coloured by taxonomy, with the clades as a strip next to the tips.
# text_size: strip label size in pt; NULL keeps the default (about 8.5 pt)
plot_tree <- function(panel, text_size = NULL) {
    label_size <- if (is.null(text_size)) 3 else text_size / .pt
    tree <- ape::read.tree(panel$tree)
    tip_clades <- read_tip_clades(panel)
    clade_colors <- clade_palette(tip_clades$clade, panel)

    # A node takes a taxonomic group when all its descendant tips belong to it
    tree_df <- as_tibble(tree) %>%
        left_join(select(tip_taxonomy(tree), label, group), by = 'label')
    node_group <- map_chr(tree_df$node, function(n) {
        tips <- as.character(tidytree::offspring(tree_df, n, tiponly = TRUE, self_include = TRUE)$group)
        if (n_distinct(tips) == 1) tips[1] else 'mixed'
    })
    tree_data <- tree_df %>%
        mutate(group = factor(node_group, levels = c(names(group_colors), 'mixed'))) %>%
        tidytree::as.treedata()

    p <- ggtree::ggtree(tree_data, aes(color = group), size = 0.3)

    # Clades as a strip next to the tips
    xmax <- max(p$data$x)
    strip <- p$data %>%
        filter(isTip) %>%
        left_join(tip_clades, by = 'label') %>%
        arrange(y)
    # One rectangle per run of consecutive tips of the same clade
    runs <- strip %>%
        mutate(run = consecutive_id(clade)) %>%
        filter(!is.na(clade)) %>%
        group_by(clade, run) %>%
        summarize(ymin = min(y) - 0.5, ymax = max(y) + 0.5, n = n(), .groups = 'drop')
    # Labels on each clade's longest run, pushed apart when small clades are close together
    min_gap <- nrow(strip) * 0.04 * if (isTRUE(panel$colored_labels)) label_size * 1.4 / 3 else 1
    labels <- runs %>%
        group_by(clade) %>%
        slice_max(n, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(y = (ymin + ymax) / 2) %>%
        arrange(y) %>%
        mutate(y0 = y, y = accumulate(y, ~max(.y, .x + min_gap)))
    # Coloured labels: centre each group of pushed-apart labels on its clades instead of pushing it
    # only upwards, so the labels stay next to their strips
    if (isTRUE(panel$colored_labels)) {
        labels <- labels %>%
            mutate(group = cumsum(c(TRUE, diff(y) > min_gap + 1e-9))) %>%
            group_by(group) %>%
            mutate(y = y - (mean(y) - mean(y0))) %>%
            ungroup()
    }

    p <- p +
        # Dotted leaders from the tips of each clade out to the strip, in the clade colour
        geom_segment(data = filter(strip, !is.na(clade)),
                     aes(x = x, xend = xmax * 1.02, y = y, yend = y, color = clade),
                     inherit.aes = FALSE, linetype = 'dotted', linewidth = 0.25, alpha = 0.6,
                     show.legend = FALSE) +
        # Invisible segments with every group draw the taxonomy legend (the tree layer would only
        # draw the groups it has), so one tree's legend covers both
        geom_segment(data = tibble(group = factor(names(group_colors), levels = names(group_colors))),
                     aes(x = 0, xend = 0, y = 1, yend = 1, color = group), inherit.aes = FALSE,
                     linewidth = 0) +
        # Thin line behind the strip, joining it to its title
        (if (isTRUE(panel$strip_guide))
             annotate('segment', x = xmax * 1.04, xend = xmax * 1.04, y = 0.5, yend = nrow(strip) * 1.01,
                      color = 'grey80', linewidth = 0.3)
         else NULL) +
        geom_rect(data = runs, aes(xmin = xmax * 1.02, xmax = xmax * 1.06, ymin = ymin, ymax = ymax,
                                   fill = clade), inherit.aes = FALSE) +
        geom_text(data = labels, aes(x = xmax * 1.08, y = y, label = panel$strip_labels(clade)),
                  inherit.aes = FALSE, lineheight = 0.9,
                  size = if (isTRUE(panel$colored_labels)) label_size * 1.4 else label_size,
                  color = if (isTRUE(panel$colored_labels)) clade_colors[labels$clade] else 'grey20',
                  fontface = if (isTRUE(panel$colored_labels)) 'bold' else 'plain',
                  angle = if (isTRUE(panel$strip_vertical)) 90 else 0,
                  hjust = if (isTRUE(panel$strip_vertical)) 0.5 else 0,
                  vjust = if (isTRUE(panel$strip_vertical)) 1 else 0.5) +
        # Header naming what the strip shows
        annotate('text', x = xmax * 1.04, y = nrow(strip) + nrow(strip) * 0.015, label = panel$strip_title,
                 vjust = 0, hjust = 0.5, size = label_size * 0.93, lineheight = 0.9, color = 'grey20',
                 fontface = panel$strip_title_face %||% 'plain') +
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.08))) +
        # The colour scale also holds the clade colours for the leaders; only taxonomy is in the legend
        scale_color_manual(values = c(group_colors, mixed = mixed_color, clade_colors),
                           breaks = names(group_colors),
                           name = 'Taxonomy') +
        scale_fill_manual(values = clade_colors) +
        # Room on the right for the strip labels
        scale_x_continuous(expand = expansion(mult = c(0.01, panel$strip_room %||%
                                                            if (is.null(panel$facet_labels)) 0.04 else 0.25))) +
        # Clades are identified by the strip numbers and the violin facets, so no fill legend;
        # taxonomy legend inside the tree panel (empty top-left corner)
        guides(fill = 'none', color = guide_legend(override.aes = list(linewidth = 2))) +
        theme(legend.position = 'inside', legend.position.inside = c(0.02, 0.98),
              legend.justification.inside = c(0, 1), legend.background = element_blank())
    p$layers[[1]]$show.legend <- FALSE
    p
}

# Clade violins ---------------------------------------------------------------

# y axis for heavy-tailed values: log10(1 + y), as counts and ratios can be 0 or close to it
# top: room above the data, as a fraction of its range, for the p/FC labels
# step: powers of ten between breaks
clade_y_scale <- function(top = 0.25, step = 2) {
    scale_y_continuous(trans = 'log1p', breaks = 10^seq(0, 12, by = step),
                       labels = scales::label_log(), expand = expansion(mult = c(0.02, top)))
}

# basal last; clades with no reads in a sample don't have a row, so the factor is built from the data
order_clades <- function(counts) {
    counts %>%
        mutate(clade = factor(clade),
               clade = if ('basal' %in% levels(clade)) fct_relevel(clade, 'basal', after = Inf) else clade)
}

# Per-clade violins between the low/high categories of the environmental variable.
# counts: clade, count and envcat (two levels, low first).
# text_size: size of the p/FC labels and the facet titles in pt; NULL keeps the defaults (8 pt)
# log2_fc: show log2(high/low) of the medians (positive = higher with high values of the variable)
# instead of the default low/high ratio
plot_clade_violins <- function(counts, panel, text_size = NULL, point_size = 0.4, log2_fc = FALSE) {
    label_size <- if (is.null(text_size)) 2.8 else text_size / .pt
    counts <- counts %>%
        filter(!is.na(envcat)) %>%
        mutate(envcat = factor(envcat, labels = c('low', 'high'))) %>%
        order_clades()
    clade_colors <- clade_palette(levels(counts$clade), panel)

    # Wilcoxon rank-sum test (unadjusted) and median fold change, per clade
    stats <- counts %>%
        group_by(clade) %>%
        summarize(p = wilcox.test(count ~ envcat)$p.value,
                  fc = median(count[envcat == 'low']) / median(count[envcat == 'high'])) %>%
        mutate(label = if (log2_fc) str_c(scales::label_pvalue(add_p = TRUE)(p), '\nlog2 FC ', round(-log2(fc), 2))
                       else str_c(scales::label_pvalue(add_p = TRUE)(p), '\nFC ', round(fc, 1)))

    counts %>%
        ggplot(aes(x = envcat, y = count)) +
        geom_violin(aes(fill = clade, alpha = envcat), color = 'grey20', linewidth = 0.3,
                    quantile.linetype = 1, show.legend = FALSE) +
        # jitter with a fixed seed so the figure is the same in every run (the original had none)
        geom_point(position = position_jitter(width = 0.1, seed = 1), size = point_size, alpha = 0.4,
                   color = point_color) +
        geom_text(data = stats, aes(label = label, x = I(0.5), y = I(if (isTRUE(panel$colored_labels)) 0.985 else 0.93)),
                  vjust = if (isTRUE(panel$colored_labels)) 1 else 0.5, size = label_size, lineheight = 0.9) +
        facet_wrap(~clade, nrow = 2, scales = 'free_y',
                   labeller = as_labeller(clade_facet_labels(panel, levels(counts$clade), text_size))) +
        scale_fill_manual(values = clade_colors) +
        scale_alpha_manual(values = c(low = 0.45, high = 0.9)) +
        clade_y_scale(top = if (isTRUE(panel$colored_labels)) 0.45 else 0.25,
                      step = panel$y_break_step %||% 2) +
        xlab(panel$xlab %||% str_glue('{str_to_sentence(panel$variable)} category')) +
        ylab(panel$ylab) +
        theme_phylofactorization() +
        theme(legend.position = 'bottom') +
        (if (!is.null(text_size)) theme(strip.text = element_text(size = text_size)) else theme()) +
        (if (isTRUE(panel$colored_labels))
             theme(strip.text = ggtext::element_markdown(size = text_size, lineheight = 1.1, halign = 0.5))
         else theme())
}

# Reads the panel's count table, adds the low/high category and draws the violins
# text_size, point_size: NULL keeps the defaults of plot_clade_violins()
plot_violins <- function(panel, text_size = NULL, point_size = NULL, log2_fc = FALSE) {
    config <- env_variables[[panel$variable]]
    read_tsv(panel$counts, show_col_types = FALSE) %>%
        filter(!sample %in% panel$exclude_samples) %>%
        mutate(envcat = factor(.data[[config$column]] >= config$threshold, levels = c(FALSE, TRUE))) %>%
        plot_clade_violins(panel, text_size = text_size, point_size = point_size %||% 0.4, log2_fc = log2_fc) +
        theme(legend.position = 'none')
}

# One row: tree on the left, its clade violins on the right
# taxonomy_legend: "inside" the tree panel or "none"
panel_row <- function(panel, taxonomy_legend = 'inside', text_size = NULL, log2_fc = FALSE) {
    tree <- plot_tree(panel, text_size = text_size)
    if (taxonomy_legend == 'none') tree <- tree + guides(color = 'none')
    tree | plot_violins(panel, text_size = text_size, log2_fc = log2_fc)
}

# Manuscript styling of the factorization figures (Fig 7, Supplementary Fig S9) -------------------------
# Text sizes in pt at the saved size (14 in wide), about half of it once fitted to an A4 page width

cutoff_labels <- c(iron = 'Iron (cutoff: 0.5 nmol/L, Darwin model)',
                   silicate = 'Silicate (cutoff: 2 µmol/L)',
                   nitrate = 'Nitrate + nitrite (cutoff: 1 µmol/L)')
manuscript_text <- list(base = 20, strip = 13, title = 22)

manuscript_panel <- function(panel) {
    clades <- str_c('clade_', 1:5)
    modifyList(panel, list(
        facet_labels = c(setNames(str_c('Clade ', str_remove(clades, 'clade_')), clades), basal = 'Basal'),
        colored_labels = TRUE,
        strip_title_face = 'italic',
        strip_guide = TRUE,
        # the strip only holds single digits
        strip_room = 0.1,
        basal_composition = FALSE,
        xlab = cutoff_labels[[panel$variable]]))
}

# One row (tree + violins) of a manuscript figure; tag: panel letter, in bold before the title
manuscript_row <- function(panel, tag, title, ...) {
    row <- panel_row(manuscript_panel(panel), text_size = manuscript_text$strip, log2_fc = TRUE, ...) +
        patchwork::plot_layout(widths = c(1.4, 1))
    row[[1]] <- row[[1]] +
        ggtitle(str_glue('**{tag}** {title}')) +
        # room above the tips for the strip title
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.16)))
    row
}

manuscript_theme <- function() {
    theme(text = element_text(size = manuscript_text$base),
          axis.text = element_text(size = manuscript_text$base * 0.8),
          axis.title = element_text(size = manuscript_text$base - 2),
          legend.text = element_text(size = manuscript_text$base * 0.75),
          legend.title = element_text(size = manuscript_text$base * 0.85),
          plot.title = ggtext::element_markdown(size = manuscript_text$title),
          plot.title.position = 'plot')
}
