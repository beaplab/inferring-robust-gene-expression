# Supplementary Figure S7: Flavodoxin, housekeeping-gene (HK) vs metaG normalisation, per Graff et al.
# 2023 clade.
#   A. log2 fold change of the median low/high iron (iron_darwin_5m, cutoff 0.5 nmol/L) of the clade sums,
#      for metaT, metaG, the metaT/metaG ratio and two versions of the metaT/HK ratio (all samples and
#      unigenes; or only the samples and unigenes with metaG, as the metaT/metaG ratio). Filled: Wilcoxon
#      p.adj < 0.05 (Holm, per data type).
#   B. per-sample agreement of the summed metaT/HK ratio (same samples, metaG unigenes) with the summed
#      metaT/metaG ratio, per clade (clade sums of 0 left out: log axes).
# In "Rest of the tree" the metaT/metaG ratio is not significant, but the HK ratio is.
#
# Inputs: data/intermediate/supplementary/S7/graff2023_{metaT,metaG,hkratio}_counts.tsv and
#         graff2023_hkratio_pergene.tsv.gz (scripts/00_data_generation/supplementary/S7/)
#         data/intermediate/fig2/flavodoxin_graff2023_clade_ratios.tsv, flavodoxin_metaT_metaG_pergene_counts.tsv.gz
# Original: gene_explorer scripts/figures/flavodoxin_hk_vs_metag.R (chunks fc-plot and sample-agreement of
#           scripts/analysis/phyloplacement_flavodoxin_hkratio_analysis.qmd)
# Run time: ~1 min.
# Output: results/figures/supplementary/S7_flavodoxin_HK_vs_metaG.{pdf,png}

source('R/load.R')
library(patchwork)

s7_dir <- p('intermediate', 'supplementary', 'S7')
iron <- env_variables$iron   # iron_darwin_5m, threshold 5e-4
catlevels <- c('low_iron', 'high_iron')
category_colors <- c(low_iron = '#d95f02', high_iron = '#1b9e77')
category_labels <- c(low_iron = 'Low iron', high_iron = 'High iron')
fit_color <- '#2a78d6'

graff_labels <- c(clade1_dielcycle = 'Clade 1 (diel cycle)',
                  clade2_ironlimiting = 'Clade 2 (iron limiting)',
                  greenalgae_dinos = 'Green algae + dinoflagellates',
                  basal = 'Rest of the tree')
clade_labeller <- as_labeller(\(x) str_wrap(graff_labels[x], 16))

# Data ------------------------------------------------------------------------------------------------

load_counts <- function(file) {
    read_tsv(file, show_col_types = FALSE) %>%
        mutate(envcat = case_when(.data[[iron$column]] < iron$threshold ~ catlevels[1],
                                  .data[[iron$column]] >= iron$threshold ~ catlevels[2]),
               envcat = factor(envcat, levels = catlevels))
}

clade_metaT <- load_counts(file.path(s7_dir, 'graff2023_metaT_counts.tsv'))
clade_metaG <- load_counts(file.path(s7_dir, 'graff2023_metaG_counts.tsv'))
clade_metagratio <- load_counts(p('intermediate', 'fig2', 'flavodoxin_graff2023_clade_ratios.tsv'))
clade_hkratio <- load_counts(file.path(s7_dir, 'graff2023_hkratio_counts.tsv'))
clade_hkratio_shared <- filter(clade_hkratio, sample %in% clade_metagratio$sample)

pergene_metag <- read_tsv(p('intermediate', 'fig2', 'flavodoxin_metaT_metaG_pergene_counts.tsv.gz'), show_col_types = FALSE) %>%
    mutate(has_metag = is.finite(count_metat / count_metag))

# HK ratio summed over the unigenes with metaG only (same unigenes and samples as the metaT/metaG ratio)
hk_pergene_clade <- read_tsv(file.path(s7_dir, 'graff2023_hkratio_pergene.tsv.gz'), show_col_types = FALSE)
clade_hkratio_metaggenes <- hk_pergene_clade %>%
    semi_join(pergene_metag %>% filter(has_metag) %>% transmute(sample = str_remove(sample, 'filter$'), geneid),
              by = c('sample', 'geneid')) %>%
    group_by(sample, clade) %>%
    summarize(count_mg = sum(count), .groups = 'drop') %>%
    right_join(clade_hkratio_shared, by = c('sample', 'clade')) %>%
    mutate(count = replace_na(count_mg, 0)) %>%
    select(-count_mg)

# A. fold changes -------------------------------------------------------------------------------------

compare_categories <- function(counts) {
    counts %>%
        filter(!is.na(envcat)) %>%
        group_by(clade) %>%
        summarize(fc = median(count[envcat == catlevels[1]]) / median(count[envcat == catlevels[2]]),
                  p = wilcox.test(count ~ envcat)$p.value,
                  .groups = 'drop') %>%
        mutate(p.adj = p.adjust(p, method = 'holm'))
}

stats_all <- list(`metaT` = clade_metaT,
                  `metaG` = clade_metaG,
                  `metaT/metaG ratio` = clade_metagratio,
                  `metaT/HK ratio: same samples, metaG genes` = clade_hkratio_metaggenes,
                  `metaT/HK ratio: all samples, all genes` = clade_hkratio) %>%
    map(compare_categories) %>%
    bind_rows(.id = 'data') %>%
    mutate(data = fct_rev(fct_inorder(data)),
           clade = factor(clade, levels = names(graff_labels)))

pA <- stats_all %>%
    mutate(signif = p.adj < 0.05) %>%
    ggplot(aes(log2(fc), data)) +
    geom_vline(xintercept = 0, color = 'grey60') +
    geom_segment(aes(x = 0, xend = log2(fc), yend = data), color = 'grey70') +
    geom_point(aes(shape = signif), size = 2.5, color = fit_color) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = 'No', `TRUE` = 'Yes'),
                       name = 'Wilcoxon p.adj < 0.05') +
    facet_wrap(~clade, nrow = 1, labeller = clade_labeller) +
    labs(x = 'log2 fold change low/high iron (> 0: higher at low iron)', y = NULL) +
    theme_phylofactorization() +
    theme(legend.position = 'bottom')

# B. agreement per sample -----------------------------------------------------------------------------

paired <- clade_metagratio %>%
    select(clade, sample, ratio_metag = count, envcat) %>%
    inner_join(clade_hkratio_metaggenes %>% select(clade, sample, ratio_hk_mg = count), by = c('clade', 'sample'))

pB <- paired %>%
    filter(ratio_metag > 0) %>%
    mutate(clade = factor(clade, levels = names(graff_labels))) %>%
    ggplot(aes(ratio_metag, ratio_hk_mg)) +
    geom_point(aes(color = envcat), size = 0.8, alpha = 0.6) +
    geom_smooth(method = 'lm', formula = y ~ x, color = 'grey30', linewidth = 0.5, se = FALSE) +
    facet_wrap(~clade, nrow = 1, scales = 'free', labeller = clade_labeller) +
    scale_x_log10(labels = scales::label_log()) +
    scale_y_log10(labels = scales::label_log()) +
    scale_color_manual(values = category_colors, labels = category_labels, na.value = 'grey70', name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = 'Summed metaT/metaG ratio', y = 'Summed metaT/HK ratio\n(same samples, metaG genes)') +
    theme_phylofactorization() +
    theme(legend.position = 'bottom')

# free(): panel B is not aligned to the long y labels of panel A
fig <- free(pA) / pB + plot_annotation(tag_levels = 'A')

dir.create(p('figures', 'supplementary'), showWarnings = FALSE)
save_figure(fig, 'supplementary/S7_flavodoxin_HK_vs_metaG', width = 10, height = 8)
