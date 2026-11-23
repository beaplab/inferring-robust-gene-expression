# Supplementary Figure S4: placing metaT reads directly vs placing the MATOU assembled unigenes, for 7
# housekeeping genes (same reference trees, metazoan placements removed in both).
#   A. total multiplicity (read pairs) per sample for each approach, by size fraction (dashed: 1:1)
#   B. log2 ratio of sequence-based to read-based multiplicity, by size fraction and depth
#
# Input:  data/intermediate/supplementary/S4/aggregated/<gene>.rds
#         (scripts/00_data_generation/supplementary/S4/01_read_vs_seq_placement_caches.R)
# Original: gene_explorer scripts/figures/visualize_read_vs_seq_placement.R (section 1; its sections 2-3,
#           placement resolution and taxonomic assignment, are not in the paper)
# Run time: ~20 s.
# Output: results/figures/supplementary/S4_read_vs_seq_placement_recovery.{pdf,png}
#         results/tables/supplementary/read_vs_seq_placement_{sample_summary,recovery_per_gene}.tsv

source('R/load.R')
library(patchwork)

fractions <- c('0o8-5', '5-20', '20-180', '180-2000', '0o8-inf')
depth_colors <- c(SUR = '#9ecae1', DCM = '#3182bd')
tables_dir <- p('tables', 'supplementary')

# Data ------------------------------------------------------------------------------------------------
cache_files <- list.files(p('intermediate', 'supplementary', 'S4', 'aggregated'), '\\.rds$', full.names = TRUE) %>%
    set_names(basename(.) %>% str_remove('\\.rds$'))
cache <- map(cache_files, readRDS)
bind_cache <- function(element) map_dfr(cache, element)

samples <- bind_cache('samples') %>%
    mutate(fraction = factor(fraction, levels = fractions))
sample_info <- distinct(samples, gene, sample, depth, fraction)

taxonomy <- bind_cache('taxonomy') %>%
    mutate(method = factor(method, levels = c('reads', 'sequences'))) %>%
    left_join(sample_info, by = c('gene', 'sample'))

# Amount of data recovered per sample and method ------------------------------------------------------
sample_summary <- taxonomy %>%
    group_by(gene, sample, depth, fraction, method) %>%
    summarize(total_multiplicity = sum(multiplicity), n_units = sum(n_units), .groups = 'drop') %>%
    complete(nesting(gene, sample, depth, fraction), method, fill = list(total_multiplicity = 0, n_units = 0)) %>%
    left_join(samples %>% select(gene, sample, seq_total_before_extraction), by = c('gene', 'sample')) %>%
    mutate(seq_total_before_extraction = if_else(method == 'sequences', seq_total_before_extraction, NA_real_))
write_tsv(sample_summary, file.path(tables_dir, 'read_vs_seq_placement_sample_summary.tsv'))

recovery_wide <- sample_summary %>%
    select(gene, sample, depth, fraction, method, total_multiplicity, n_units) %>%
    pivot_wider(names_from = method, values_from = c(total_multiplicity, n_units)) %>%
    mutate(log2_ratio_seq_reads = log2((total_multiplicity_sequences + 1) / (total_multiplicity_reads + 1)))

recovery_gene <- recovery_wide %>%
    group_by(gene) %>%
    summarize(n_samples = n(),
              total_reads = sum(total_multiplicity_reads),
              total_sequences = sum(total_multiplicity_sequences),
              ratio_seq_reads_pooled = total_sequences / total_reads,
              median_log2_ratio = median(log2_ratio_seq_reads),
              spearman_rho = cor(total_multiplicity_reads, total_multiplicity_sequences, method = 'spearman'),
              units_reads = sum(n_units_reads),
              units_sequences = sum(n_units_sequences),
              .groups = 'drop') %>%
    left_join(sample_summary %>%
                  filter(method == 'sequences') %>%
                  group_by(gene) %>%
                  summarize(seq_metazoa_removed_frac = 1 - sum(total_multiplicity) / sum(seq_total_before_extraction)),
              by = 'gene')
write_tsv(recovery_gene, file.path(tables_dir, 'read_vs_seq_placement_recovery_per_gene.tsv'))
print(recovery_gene)

# Figure ------------------------------------------------------------------------------------------------
recovery_scatter <- recovery_wide %>%
    ggplot(aes(total_multiplicity_reads + 1, total_multiplicity_sequences + 1, color = fraction)) +
    geom_abline(linetype = 'dashed', color = 'grey50') +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_x_continuous(breaks = scales::breaks_extended(n = 3), labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    facet_wrap(~gene, nrow = 1, scales = 'free') +
    labs(x = 'Total multiplicity, read placement (+1)',
         y = 'Total multiplicity,\nsequence placement (+1)',
         color = 'Size fraction (µm)') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

recovery_ratio <- recovery_wide %>%
    ggplot(aes(fraction, log2_ratio_seq_reads, fill = depth)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
    geom_boxplot(outlier.size = 0.6) +
    facet_wrap(~gene, nrow = 1) +
    scale_fill_manual(values = depth_colors) +
    labs(x = 'Size fraction (µm)', y = 'log2(sequences / reads)\ntotal multiplicity', fill = 'Depth') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

recovery_plot <- recovery_scatter / recovery_ratio +
    plot_annotation(title = 'Amount of data recovered: read vs sequence placement', tag_levels = 'A')

dir.create(p('figures', 'supplementary'), showWarnings = FALSE)
save_figure(recovery_plot, 'supplementary/S4_read_vs_seq_placement_recovery', width = 11, height = 7.5)
