# Figure 3: genetic context of Florenciella parvula, and what it means for quantifying its expression
# by mapping metatranscriptomic reads to its transcriptome.
#   A: amino-acid identity (AAI, sourmash) between the Florenciella transcriptomes, three Florenciella
#      SMAGs and Pelagomonas calceolata (outgroup); dendrogram = average-linkage clustering on 1 - AAI,
#      branches coloured by the gene classification of panel D that each species drives
#   B: crossmapping of the 400 most expressed F. parvula genes with every other species
#   C: coverage of three example genes (Calmodulin, RPL26, RplB), F. parvula vs crossmapping species
#   D: classification of the F. parvula genes by their crossmapping
#   E: number of expressed genes per sample and species (Tara, Carradec et al. 2022); red: samples with
#      >= 1000 expressed genes
#   F: number of expressed genes of F. parvula vs Florenciella sp. RCC1587 and RCC1007 per sample
#
# Inputs (data/intermediate/fig3, made by scripts/00_data_generation/fig3/):
#   compar_transcriptomes.csv                                 01_aai_sourmash.sh
#   blast_Florenciella_allgenomes.tsv.gz                      02_blast_florenciella.sh
#   crossmapping_florenciella_against_all_transcriptomes.tsv.gz  03_crossmapping.R
#   allstatistics_florenciella_plus_other.tsv.gz              04_expressed_genes_per_sample.R
#   coverage_examples.tsv                                     05_coverage_examples.R
#   + data/raw/quantification/alignment/TPM_EP00618_Florenciella_parvula_quant-over_2012_carradec_tara.tsv
# Original: delving_into_florenciella src/figures/combined_florenciella_figure.R
# Output: results/figures/fig3_florenciella_context.{pdf,png}

source('R/load.R')
suppressPackageStartupMessages({
    library(patchwork)
    library(viridis)
    library(ggtext)
    library(ggridges)
    library(ggmagnify)
})

fig3_dir <- p('intermediate', 'fig3')

# Names -------------------------------------------------------------------------

# file / transcriptome names -> short names, the keys used to join the tables
shortnames <- c(
    EP00618_Florenciella_parvula = 'Florenciella parvula',
    EP00619_Florenciella_sp_RCC1007 = 'RCC1007',
    EP00620_Florenciella_sp_RCC1587 = 'RCC1587',
    EP00621_Florenciella_sp_RCC1693 = 'RCC1693',
    pelagomonas_aa = 'Pelagomonas',
    Pelagomonas_calceolata = 'Pelagomonas',
    TARA_ARC_108_MAG_00262 = 'TARAARC108',
    TARA_MED_95_MAG_00475 = 'TARAMED95',
    TARA_SOC_28_MAG_00069 = 'TARASOC28',
    TARA_ARC_108_MAG_00262_CDS = 'TARAARC108',
    TARA_MED_95_MAG_00475_CDS = 'TARAMED95',
    TARA_SOC_28_MAG_00069_CDS = 'TARASOC28'
)

# display labels (markdown, rendered with ggtext::element_markdown)
display_names <- c(
    'Florenciella parvula' = '*Florenciella parvula*',
    RCC1007 = '*Florenciella* sp. RCC1007',
    RCC1587 = '*Florenciella* sp. RCC1587',
    RCC1693 = '*Florenciella* sp. RCC1693',
    Pelagomonas = '*Pelagomonas calceolata*',
    TARAARC108 = 'TARAARC108',
    TARAMED95 = 'TARAMED95',
    TARASOC28 = 'TARASOC28'
)

# Gene classification of panel D, also used for the dendrogram of panel A
status_levels <- c('Specific of\nspecies level', 'Crossmapping with \nsister species',
                   'Established at\ngenus level', 'To be removed')
status_colors <- setNames(inferno(5)[2:5], status_levels)
species_status <- c(
    'Florenciella parvula' = status_levels[1],
    RCC1587 = status_levels[2],
    RCC1693 = status_levels[2],
    RCC1007 = status_levels[3],
    TARAARC108 = status_levels[3],
    TARAMED95 = status_levels[3],
    TARASOC28 = status_levels[3],
    Pelagomonas = status_levels[4]
)

# Crossmapping categories (03_crossmapping.R), from the reads of a F. parvula gene that also map to the
# other species, divided by those that do not (summed over the samples). 'only reads crossmapped!'
# (every read crossmaps) is pooled with the large category.
crossmap_levels <- c('no crossmapping', 'small amount of crossmapping (<10%)',
                     'intermediate situation (10-50%)', 'large crossmapping (>50%)')
pool_only_crossmapped <- function(category) {
    ifelse(category == 'only reads crossmapped!', 'large crossmapping (>50%)', category)
}

# Left indent (pt) for the titles of panels C and E, so they start at the same x as panel A's title,
# which sits right of the dendrogram column
panel_title_indent <- 36

# Data --------------------------------------------------------------------------

anidf <- read_csv(file.path(fig3_dir, 'compar_transcriptomes.csv'), show_col_types = FALSE) %>%
    mutate(species = colnames(.)) %>%
    filter(!str_detect(species, 'TOPAZ')) %>%
    select(!contains('TOPAZ')) %>%
    select(-species)
colnames(anidf) <- basename(colnames(anidf)) %>% str_remove('.faa.gz|.faa|.gmove.pep.faa|.fasta')
anidf <- rename_with(anidf, ~ shortnames[.x])

ani_mat <- as.matrix(anidf)
rownames(ani_mat) <- colnames(ani_mat)

# Species order shared by panels A and B: average-linkage clustering of 1 - AAI, so the dendrogram
# next to panel A is a crossing-free picture of the row order
species_hc <- hclust(as.dist(1 - ani_mat), method = 'average')
taxorder <- species_hc$labels[species_hc$order]

florenciella_tpm <- read_tsv(p('quant_alignment', 'TPM_EP00618_Florenciella_parvula_quant-over_2012_carradec_tara.tsv'),
                             show_col_types = FALSE)

# F. parvula genes ranked by mean TPM over the samples
rank_ab_gene <- florenciella_tpm %>%
    pivot_longer(names_to = 'gene_name', values_to = 'tpm', cols = -sample) %>%
    group_by(gene_name) %>%
    summarize(mean.tpm = mean(tpm)) %>%
    arrange(-mean.tpm) %>%
    mutate(rank = 1:nrow(.))

crossmapping <- read_tsv(file.path(fig3_dir, 'crossmapping_florenciella_against_all_transcriptomes.tsv.gz'),
                         show_col_types = FALSE)

# F. parvula genes with a blastn hit (>= 80% identity over >= 100 nt) in each other species
selection_present <- read_tsv(file.path(fig3_dir, 'blast_Florenciella_allgenomes.tsv.gz'), show_col_types = FALSE) %>%
    mutate(subject = recode(subject, Pelagomonas_calceolata = 'Pelagomonas', !!!as.list(shortnames))) %>%
    distinct(subject, qseqid) %>%
    mutate(presence = 'yes')

# Expressed genes (TPM > 0) per sample, species and dataset
allspecies <- read_tsv(file.path(fig3_dir, 'allstatistics_florenciella_plus_other.tsv.gz'), show_col_types = FALSE) %>%
    filter(str_detect(taxon, 'TARA|Florenciella|Pelagomonas')) %>%
    mutate(taxon = factor(shortnames[taxon], levels = taxorder))

coverage <- read_tsv(file.path(fig3_dir, 'coverage_examples.tsv'), show_col_types = FALSE)

# Panel A: AAI heatmap + dendrogram ---------------------------------------------

ani_long <- anidf %>%
    mutate(genomes = colnames(.)) %>%
    pivot_longer(names_to = 'genomes2', values_to = 'ani', cols = -genomes) %>%
    mutate(genomes = factor(genomes, levels = taxorder),
           genomes2 = factor(genomes2, levels = taxorder),
           textcolor = ifelse(ani > 0.75, 'white', 'black'))

ani_heatmap_gg <- ggplot(ani_long, aes(x = genomes2, y = genomes)) +
    geom_tile(aes(fill = ani), color = 'white') +
    geom_text(aes(label = sprintf('%.2f', ani), color = textcolor), size = 3) +
    scale_fill_gradient(low = '#E8C4C4', high = 'firebrick', name = 'AAI',
                        guide = guide_colorbar(barwidth = unit(6, 'cm'), barheight = unit(0.4, 'cm'))) +
    scale_color_identity() +
    scale_y_discrete(position = 'right', labels = display_names) +
    labs(caption = 'TARAARC108, TARAMED95 and TARASOC28 are partial SMAG genomes.',
         # wrapped: patchwork clips text to its own viewport
         title = str_wrap('What is the overall genomic relatedness of the species?', 40)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y.right = element_markdown(hjust = 1),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'bottom',
          plot.margin = margin(r = 6, t = 6),
          plot.caption = element_text(hjust = 0, size = 8, color = 'grey30', margin = margin(t = 2)))

# Segments of a dendrogram (hclust) with the leaves at leaf_x. A branch keeps the group colour of its
# leaves while they all share it, and is black (NA) above where groups merge.
build_dendro_segments <- function(hc, leaf_x, leaf_group) {
    n <- length(hc$order)
    node_x <- numeric(n - 1)
    node_group <- rep(NA_character_, n - 1)

    x_of <- function(node) if (node < 0) leaf_x[[hc$labels[-node]]] else node_x[node]
    y_of <- function(node) if (node < 0) 0 else hc$height[node]
    group_of <- function(node) if (node < 0) leaf_group[[hc$labels[-node]]] else node_group[node]

    segments <- vector('list', n - 1)
    for (i in seq_len(n - 1)) {
        left <- hc$merge[i, 1]
        right <- hc$merge[i, 2]
        x_left <- x_of(left)
        x_right <- x_of(right)
        node_x[i] <- (x_left + x_right) / 2
        g_left <- group_of(left)
        g_right <- group_of(right)
        node_group[i] <- if (identical(g_left, g_right)) g_left else NA_character_

        segments[[i]] <- tibble(
            x = c(x_left, x_right, x_left),
            y = c(y_of(left), y_of(right), hc$height[i]),
            xend = c(x_left, x_right, x_right),
            yend = c(hc$height[i], hc$height[i], hc$height[i]),
            group = c(g_left, g_right, node_group[i]))
    }
    bind_rows(segments)
}

n_species <- length(taxorder)
dendro_segments <- build_dendro_segments(species_hc, setNames(seq_along(taxorder), taxorder), species_status)

dendro_gg <- ggplot(dendro_segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = group),
                 linewidth = 0.9, lineend = 'square') +
    # dots at the leaf tips: the within-Florenciella branches are very short
    geom_point(data = tibble(x = seq_along(taxorder), y = 0, group = species_status[taxorder]),
               aes(x = x, y = y, color = group), size = 3) +
    scale_color_manual(values = status_colors, na.value = 'black', guide = 'none') +
    # 0.4/0.6 matches the default discrete-axis expansion, so leaves line up with the heatmap rows
    scale_x_continuous(limits = c(0.4, n_species + 0.6), expand = c(0, 0)) +
    scale_y_reverse(expand = expansion(mult = c(0.02, 0.05))) +
    coord_flip(clip = 'off') +
    labs(tag = 'A') +
    theme_void() +
    theme(plot.margin = margin(r = 2, t = 6, b = 6))

ani_row <- (dendro_gg | ani_heatmap_gg) + plot_layout(widths = c(0.4, 5))

# legend stacked below the caption (ggplot draws the caption below the legend)
ani_legend <- cowplot::get_legend(ani_heatmap_gg)
ani_gg <- (ani_row & theme(legend.position = 'none')) /
    wrap_elements(full = ani_legend) +
    plot_layout(heights = c(5, 2))

# Panel B: crossmapping of the most expressed genes ----------------------------

# Crossmapping category of the top-ranked F. parvula genes with one species
rank_crossmapping <- function(species) {
    rank_ab_gene %>%
        filter(rank < 2000) %>%
        left_join(filter(crossmapping, comparison == species), by = 'gene_name') %>%
        mutate(category_cross = factor(pool_only_crossmapped(category_cross), levels = crossmap_levels))
}

comparisons <- intersect(taxorder, unique(crossmapping$comparison))
names(comparisons) <- comparisons
rank_cross <- map(comparisons, rank_crossmapping)

# 'gene absent' (no BLAST hit in the other species) is white, the four crossmapping categories a ramp
category_levels <- c('gene absent // no correspondence in genome', crossmap_levels)
category_labels <- c('gene absent\n(no correspondence\nin genome)', str_wrap(crossmap_levels, 15))
category_colors <- setNames(c('white', colorRampPalette(c('#E8C4C4', 'firebrick'))(4)), category_labels)

# The legend is split in two: 'gene absent' on its own (white key with a grey border) and the four
# crossmapping categories under a 'gene present' title
absent_legend <- cowplot::get_legend(
    ggplot(tibble(x = 1, y = 1, f = category_labels[1]), aes(x, y, fill = f)) +
        geom_tile(color = 'grey50', linewidth = 0.5) +
        scale_fill_manual(values = category_colors, name = NULL) +
        theme_minimal() +
        theme(legend.position = 'bottom') +
        guides(fill = guide_legend(keywidth = unit(1.4, 'cm'))))

plot_crossmap_heatmap <- function(rank_cross, rank_pos = 400) {
    all <- bind_rows(rank_cross, .id = 'comparison') %>%
        left_join(selection_present, by = c('gene_name' = 'qseqid', 'comparison' = 'subject')) %>%
        mutate(comparison = factor(comparison, levels = taxorder),
               category_cross = ifelse(is.na(presence), category_levels[1], as.character(category_cross)),
               category_cross = factor(category_cross, levels = category_levels, labels = category_labels)) %>%
        filter(rank < rank_pos)

    # F. parvula is the reference, never a comparison: its row is painted solid
    reference_row_y <- as.numeric(factor('Florenciella parvula', levels = taxorder))
    reference_row <- tibble(ymin = reference_row_y - 0.5, ymax = reference_row_y + 0.5)

    ggplot(all, aes(y = comparison, x = rank)) +
        geom_rect(data = reference_row, inherit.aes = FALSE, fill = 'firebrick',
                  aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf)) +
        geom_tile(aes(fill = category_cross)) +
        scale_fill_manual(values = category_colors, breaks = category_labels[-1],
                          name = 'Gene present', drop = FALSE) +
        scale_x_continuous(expand = expansion(mult = c(0, 0))) +
        scale_y_discrete(drop = FALSE) +
        theme_minimal() +
        xlab('*Florenciella parvula* genes ranked by expression') +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_markdown(),
              axis.ticks.y = element_blank(),
              plot.title = element_markdown(),
              legend.position = 'bottom',
              plot.margin = margin(l = 2)) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                                   keywidth = unit(1.4, 'cm'),
                                   theme = theme(legend.key.spacing.x = unit(0.5, 'cm'),
                                                 legend.key.spacing.y = unit(0.25, 'cm')),
                                   title.position = 'top'))
}

crossmap_gg <- plot_crossmap_heatmap(rank_cross, rank_pos = 400) +
    labs(tag = 'B', title = 'How does this relatedness impacts crossmapping between *Florenciella parvula* and<br>other species?')

# legend in its own row, with the same height ratio as panel A, so both heatmaps stay top-aligned
crossmap_legend <- cowplot::plot_grid(cowplot::get_legend(crossmap_gg), absent_legend,
                                      nrow = 1, rel_widths = c(2, 1))
crossmap_gg <- (crossmap_gg + theme(legend.position = 'none')) /
    wrap_elements(full = crossmap_legend) +
    plot_layout(heights = c(5, 2))

# Panel C: coverage of three example genes ----------------------------------------

parvula_fill <- '#86bbd8'
crossmap_fill <- '#2f4858'

# Coverage of a F. parvula gene (light) over that of its crossmapping counterpart (dark), if any
plot_coverage <- function(example_id, title, fill, textcolor) {
    cov <- filter(coverage, example == example_id)
    parvula <- filter(cov, species == 'EP00618_Florenciella_parvula')
    other <- filter(cov, species != 'EP00618_Florenciella_parvula')
    ggplot(parvula) +
        (if (nrow(other) > 0)
             geom_rect(data = other, fill = crossmap_fill, aes(xmin = from, xmax = to, ymin = 0, ymax = coverage))
         else NULL) +
        geom_rect(fill = parvula_fill, aes(xmin = from, xmax = to, ymin = 0, ymax = coverage)) +
        theme_minimal() +
        theme(plot.title = element_textbox_simple(size = 9, fill = fill, color = textcolor,
                                                   padding = margin(4, 4, 4, 4),
                                                   margin = margin(b = 6)),
              plot.margin = margin(t = 14, r = 6, b = 6, l = 6)) +
        labs(x = 'Gene position', y = 'Coverage') +
        ggtitle(title)
}

rank11 <- plot_coverage('rank11', '**Large crossmapping, outer group**<br>Calmodulin; ranked #11.',
                        fill = '#fcffa4ff', textcolor = 'black') +
    annotate(geom = 'label', y = 13000, x = 500, label = 'Crossmapping species',
             size = 4, fill = crossmap_fill, color = 'white') +
    annotate(geom = 'label', y = 10000, x = 500, label = 'Florenciella parvula',
             size = 4, fontface = 'italic', fill = parvula_fill, color = 'white')

# the RCC1007 coverage is drawn again on top: in this gene it lies above the F. parvula one
rank238 <- plot_coverage('rank238', '**Large crossmapping, within *Florenciella* genus**<br>RPL26; ranked #238.',
                         fill = '#F98C0AFF', textcolor = 'black') +
    geom_rect(data = filter(coverage, example == 'rank238', species != 'EP00618_Florenciella_parvula'),
              fill = crossmap_fill, aes(xmin = from, xmax = to, ymin = 0, ymax = coverage))

rank22 <- plot_coverage('rank22', '**No crossmapping**<br>RplB; ranked #22.',
                        fill = '#56106eff', textcolor = 'white')

coverage_row <- (rank11 | rank238 | rank22) + plot_layout(axis_titles = 'collect')

# header plot carrying the tag and title above the wrapped row
panel_C_header <- ggplot() +
    theme_void() +
    labs(title = 'Examples of this crossmapping at the gene level', tag = 'C') +
    theme(plot.title.position = 'plot',
          plot.title = element_text(margin = margin(b = 5.5, l = panel_title_indent)))
panel_C <- (panel_C_header / wrap_elements(full = coverage_row)) + plot_layout(heights = c(0.4, 15))

# Panel D: gene classification ---------------------------------------------------
# Large crossmapping with the outer group (Pelagomonas): remove the gene; large crossmapping with all
# the other Florenciella transcriptomes: expression established at the genus level; with RCC1587
# only: crossmapping with the sister species; anything else: specific of F. parvula.

gene_flags <- crossmapping %>%
    mutate(category_cross = pool_only_crossmapped(category_cross)) %>%
    select(comparison, gene_name, category_cross) %>%
    # 24 genes have several functional annotations (one row each)
    distinct(.keep_all = TRUE) %>%
    pivot_wider(names_from = comparison, values_from = category_cross, id_cols = gene_name) %>%
    group_by(gene_name) %>%
    summarize(classification = case_when(
        Pelagomonas == 'large crossmapping (>50%)' ~ status_levels[4],
        RCC1587 == 'large crossmapping (>50%)' &
            RCC1693 == 'large crossmapping (>50%)' &
            RCC1007 == 'large crossmapping (>50%)' ~ status_levels[3],
        RCC1587 == 'large crossmapping (>50%)' ~ status_levels[2],
        TRUE ~ status_levels[1])) %>%
    mutate(classification = factor(classification, levels = status_levels))
write_tsv(mutate(gene_flags, classification = str_replace_all(classification, '\\s*\n', ' ')),
          p('tables', 'fig3_florenciella_gene_classification.tsv'))

panel_D <- gene_flags %>%
    count(classification) %>%
    ggplot(aes(x = classification, y = n, fill = classification)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.4, size = 3) +
    scale_fill_manual(values = status_colors) +
    labs(x = NULL, y = 'Number of genes', tag = 'D', title = 'Impact in the gene classifcation?') +
    theme_minimal() +
    theme(legend.position = 'none', axis.text.x = element_text(size = 8)) +
    scale_y_log10()

# Panel E: expressed genes per sample -------------------------------------------------

totals <- allspecies %>%
    group_by(taxon) %>%
    summarize(total_samples_above_1000 = sum(occurrence_genes >= 1000))

panel_E <- allspecies %>%
    ggplot(aes(x = occurrence_genes, y = taxon, height = after_stat(density))) +
    geom_vline(xintercept = 1000, linetype = 3, color = 'red') +
    geom_density_ridges(stat = 'binline', fill = 'grey', color = 'transparent', bins = 30,
                        scale = 0.95, draw_baseline = FALSE) +
    geom_text(data = totals,
              mapping = aes(x = 1100, y = taxon, label = str_c(total_samples_above_1000, ' samples')),
              hjust = 'inward', nudge_y = -0.3, color = 'red', size = 3.5, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(100, 20000)) +
    scale_x_log10(breaks = c(100, 1000, 2000, 5000, 10000)) +
    scale_y_discrete(labels = display_names) +
    theme_minimal() +
    theme(axis.text.y.left = element_markdown(hjust = 1),
          plot.title.position = 'plot',
          plot.title = element_text(margin = margin(b = 5.5, l = panel_title_indent))) +
    labs(x = 'Number of expressed genes',
         y = 'Number of samples',
         tag = 'E',
         title = str_wrap('When an species expresses enough genes to establish robust patterns?', 50))

# Panel F: expressed genes of F. parvula vs two sister species -------------------------

specific_compar <- allspecies %>%
    filter(taxon %in% c('Florenciella parvula', 'Pelagomonas', 'RCC1007', 'RCC1587')) %>%
    select(-max_TPM, -median_TPM_over_zero) %>%
    pivot_wider(names_from = taxon, values_from = occurrence_genes, id_cols = c(dataset, sample))

# zoom on the samples with few expressed genes (from), drawn enlarged at (to)
magnify_from <- c(xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)

plot_species_vs_parvula <- function(species, label, magnify_to) {
    specific_compar %>%
        ggplot(aes(.data[[species]], `Florenciella parvula`)) +
        geom_abline(slope = 1, intercept = 0, linetype = 3, color = 'grey50') +
        geom_point(alpha = 0.5, na.rm = TRUE, color = 'grey30') +
        theme_minimal() +
        geom_magnify(from = magnify_from, to = magnify_to) +
        theme(axis.title.x = element_markdown(), axis.title.y = element_markdown()) +
        ylab('*Florenciella parvula*<br>number of expressed genes') +
        xlab(str_c(label, '<br>number of expressed genes'))
}

panel_F_plot <- (plot_species_vs_parvula('RCC1587', '*Florenciella* sp. RCC1587', c(0, 5000, 13000, 20000)) +
                     plot_species_vs_parvula('RCC1007', '*Florenciella* sp. RCC1007', c(10000, 15000, 10000, 20000)) +
                     plot_layout(axes = 'collect'))
panel_F_header <- ggplot() + theme_void() +
    labs(title = 'Can we distinguish expression patterns between the species?', tag = 'F')
panel_F <- (panel_F_header / panel_F_plot) +
    plot_layout(heights = c(0.4, 15)) &
    theme(plot.tag = element_text(size = 14, face = 'bold'))

# Assemble ----------------------------------------------------------------------------
# Top row wrapped so that its panels only align with each other (A with B), not with the long y axis
# labels of panel E below

top_row <- wrap_elements(full = (ani_gg | crossmap_gg) +
                             plot_layout(widths = c(1.2, 1.9)) &
                             theme(plot.tag = element_text(size = 14, face = 'bold')))
bottom_row <- (panel_C | panel_D) + plot_layout(widths = c(3, 1.3))
breadth_row <- (panel_E | panel_F) + plot_layout(widths = c(1, 2))

figure <- (top_row / bottom_row / breadth_row) +
    plot_layout(heights = c(2, 1.3, 1.6)) &
    theme(plot.tag = element_text(size = 14, face = 'bold'))

save_figure(figure, 'fig3_florenciella_context', width = 16, height = 13)
