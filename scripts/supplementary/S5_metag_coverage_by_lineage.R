# Supplementary Figure S5: share of the metaT of Flavodoxin, SIT and NRT2 with a metagenomic counterpart,
# by lineage.
#
# Only metaT samples with a paired metaG sample (same station, depth and size fraction) are used; the
# facet title gives how many there are and how much of the gene's metaT they hold. Per lineage of the
# placement branch, the bars split the normalised metaT into unigenes with metaG reads in the paired
# sample (blue) and without them (grey: these get a metaT/metaG ratio of 0); the diamond is the same
# share counted over unigene x sample pairs. Right: share of the gene's metaT held by the lineage.
#
# Input:  data/intermediate/supplementary/S5/<gene>.tsv.gz
#         (scripts/00_data_generation/supplementary/S5/01_unigene_metaT_metaG_lineage.R)
# Original: gene_explorer scripts/figures/metag_coverage_by_lineage.R (plot of section 4 of
#           scripts/analysis/phyloplacement_ratio_factorization_comparison.qmd)
# Run time: ~30 s.
# R: the paper figure was drawn with R 4.5 / ggplot2 4.0.2 (/usr/local/bin/Rscript), which reproduces it
#    pixel for pixel; R 4.1 (/usr/bin/Rscript) gives the same numbers with small rendering differences.
# Output: results/figures/supplementary/S5_metag_coverage_by_lineage.{pdf,png,tsv} (tsv: plotted numbers)

source('R/load.R')

genes <- c('Flavodoxin', 'SIT', 'NRT2')
covered_colors <- c('metaT with metaG in the paired sample' = '#2a78d6',
                    'metaT without metaG (ratio = 0)' = '#d4d4d4')
mixed_label <- 'Branches shared by several groups'
outfile <- 'supplementary/S5_metag_coverage_by_lineage'

data <- map(set_names(genes), ~read_tsv(p('intermediate', 'supplementary', 'S5', str_c(.x, '.tsv.gz')),
                                         show_col_types = FALSE))

# Summaries ---------------------------------------------------------------------------------------------

group_levels <- c('All lineages', names(group_colors), mixed_label)

# per gene and lineage: unigene x sample pairs, share of pairs and of metaT with metaG, share of the gene
coverage <- imap_dfr(data, function(d, gene) {
    paired <- filter(d, metag_sample)
    bind_rows(paired, mutate(paired, group = 'All lineages')) %>%
        mutate(group = recode(group, mixed = mixed_label)) %>%
        group_by(group) %>%
        summarize(pairs = n(),
                  pairs_covered = mean(has_metag),
                  mass_covered = sum(metaT[has_metag]) / sum(metaT),
                  share_of_gene = sum(metaT) / sum(paired$metaT),
                  .groups = 'drop') %>%
        mutate(gene = gene)
}) %>%
    mutate(gene = factor(gene, levels = genes),
           group = factor(group, levels = group_levels))

sample_info <- imap_dfr(data, function(d, gene) {
    tibble(gene = gene,
           samples = n_distinct(d$sample_code),
           paired_samples = n_distinct(d$sample_code[d$metag_sample]),
           paired_mass = sum(d$metaT[d$metag_sample]) / sum(d$metaT))
}) %>%
    mutate(gene = factor(gene, levels = genes),
           facet = str_glue('{gene}\n{paired_samples}/{samples} metaT samples with paired metaG\n',
                            '(holding {scales::percent(paired_mass, 1)} of the {gene} metaT)'))

dir.create(p('figures', 'supplementary'), showWarnings = FALSE)
write_tsv(left_join(coverage, sample_info, by = 'gene') %>% select(-facet), str_c(p('figures', outfile), '.tsv'))

# Figure --------------------------------------------------------------------------------------------------

# every lineage on the same row in all genes, so a single y axis
plot_df <- coverage %>% mutate(row = fct_rev(group))
bars <- plot_df %>%
    transmute(gene, row, group,
              `metaT with metaG in the paired sample` = mass_covered,
              `metaT without metaG (ratio = 0)` = 1 - mass_covered) %>%
    pivot_longer(c(-gene, -row, -group), names_to = 'part', values_to = 'value') %>%
    mutate(part = factor(part, levels = rev(names(covered_colors))))
inside <- 0.2   # covered share above which its label goes inside the bar (white)

fig <- ggplot(plot_df, aes(y = row)) +
    geom_col(data = bars, aes(x = value, fill = part), width = 0.72) +
    geom_point(aes(x = pairs_covered, shape = 'Unigene × sample pairs with metaG (unweighted)'),
               size = 2, fill = 'black', colour = 'white', stroke = 0.4) +
    geom_text(aes(x = mass_covered, label = scales::percent(mass_covered, 1)),
              hjust = 1.15, size = 2.6, colour = 'white', fontface = 'bold',
              data = ~filter(.x, mass_covered >= inside)) +
    geom_text(aes(x = mass_covered, label = scales::percent(mass_covered, 1)),
              hjust = -0.2, size = 2.6, colour = 'grey15', fontface = 'bold',
              data = ~filter(.x, mass_covered < inside)) +
    geom_text(aes(x = 1.03, label = scales::percent(share_of_gene, 0.1)), hjust = 0, size = 2.5, colour = 'grey30') +
    geom_text(data = sample_info, aes(x = 1.03, y = Inf, label = 'Share of\ngene metaT'),
              hjust = 0, vjust = -0.25, size = 2.3, colour = 'grey30', lineheight = 0.9, inherit.aes = FALSE) +
    facet_wrap(~gene, nrow = 1) +
    scale_fill_manual(values = covered_colors, breaks = names(covered_colors), name = NULL) +
    scale_shape_manual(values = c('Unigene × sample pairs with metaG (unweighted)' = 23), name = NULL) +
    scale_x_continuous(labels = scales::label_percent(), breaks = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1.22), expand = expansion(add = c(0, 0))) +
    coord_cartesian(clip = 'off') +
    labs(x = 'Share of the lineage metaT (normalised expression), in samples with a paired metaG sample',
         y = 'Lineage of the placement branch',
         title = 'How much of the metaT has a metagenomic counterpart') +
    theme_minimal(base_size = 9) +
    theme(legend.position = 'bottom',
          strip.text = element_text(hjust = 0, size = 10, face = 'bold', margin = margin(b = 14)),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.x = unit(1.6, 'lines'),
          plot.title = element_text(face = 'bold', size = 12),
          axis.text.y = element_text(colour = 'grey15'),
          plot.title.position = 'plot', plot.margin = margin(5, 30, 5, 5))

save_figure(fig, outfile, width = 11, height = 4.2, pdf_device = cairo_pdf)
