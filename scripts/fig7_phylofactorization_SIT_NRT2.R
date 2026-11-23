# Figure 7: placement factorization of SIT over iron and of NRT2 over nitrate, on metaT/HK ratios.
#
# Per gene: reference tree (branches coloured by taxonomy, strip next to the tips = clades found by the
# placement factorization of the per-unigene metaT/HK ratios against the low/high category of the
# environmental variable) and, next to it, the summed metaT/HK ratio per clade and sample in low and
# high samples. Labels: Wilcoxon rank-sum p-value (unadjusted) and log2 fold change of the medians
# (high/low).
#
# Inputs (data/intermediate/fig7, made by scripts/00_data_generation/fig7/):
#   SIT_iron_hkratio_factorization/, NRT2_nitrate_hkratio_factorization/
#       03_ratio_placement_factorization.R (tree, clades, ratio_df_w_metadata.tsv)
# NRT2 nitrate includes the polar samples (nitrate matched by station and layer, see
# scripts/00_data_generation/fig7/01_sample_environment.R), so panel B differs from the endpoint figure,
# which was made without them.
# Original: gene_explorer scripts/figures/summary_phylofactorization.R (+ phylofactorization_panels.R)
# Output: results/figures/fig7_phylofactorization_SIT_NRT2.{pdf,png}

source('R/load.R')
suppressPackageStartupMessages({
    library(ggtree)
    library(patchwork)
})

fig7_dir <- p('intermediate', 'fig7')
hk_ylab <- 'Summed metaT/HK ratio'

panels <- list(
    SIT = factorization_panel(file.path(fig7_dir, 'SIT_iron_hkratio_factorization'), 'sit-categorical-',
                              gene = 'SIT', variable = 'iron', counts = 'ratio_df_w_metadata.tsv', ylab = hk_ylab),
    NRT2 = factorization_panel(file.path(fig7_dir, 'NRT2_nitrate_hkratio_factorization'), 'nrt2-categorical-',
                               gene = 'NRT2', variable = 'nitrate', counts = 'ratio_df_w_metadata.tsv', ylab = hk_ylab)
) %>%
    map(~modifyList(.x, list(y_break_step = 1)))   # one power of ten between the violin y breaks

# manuscript_row() and manuscript_theme(): R/phylofactorization.R (shared with Supplementary Fig S9)

# Figure ----------------------------------------------------------------------

# Both trees share the taxonomy palette, so its legend is only drawn in the first one
figure <- (manuscript_row(panels$SIT, 'A', 'SIT placement factorization over iron concentration') /
           manuscript_row(panels$NRT2, 'B', 'NRT2 placement factorization over nitrate concentration',
                          taxonomy_legend = 'none')) &
    manuscript_theme()

save_figure(figure, 'fig7_phylofactorization_SIT_NRT2', width = 14, height = 13, pdf_device = cairo_pdf)
