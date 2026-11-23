# Supplementary Figure S9: placement factorization of SIT and NRT2 on metaT/metaG ratios, for comparison
# with Fig 7 (metaT/HK ratios); same layout.
#   A: SIT over iron (low/high cutoff 0.5 nmol/L, Darwin model, 5 m)
#   B: NRT2 over nitrate + nitrite (cutoff 1 µmol/L)
# Left: reference tree, branches coloured by taxonomy; the strip next to the tips shows the clades found
# by the placement factorization of the per-unigene metaT/metaG ratios (raw metaT reads over the metaG
# abundance of the same unigene in the paired metaG sample; 0 without metaG). Right: summed ratio per
# clade and sample in low and high samples. Labels: Wilcoxon rank-sum p-value (unadjusted) and log2 fold
# change of the medians (high/low; NaN when a median is 0).
#
# Inputs (data/intermediate/supplementary/S9): SIT_iron_metagratio_factorization/,
#   NRT2_nitrate_metagratio_factorization/ (tree, clades, ratio_df_w_metadata.tsv), made by
#   scripts/00_data_generation/supplementary/S9/01_metaT_metaG_ratios.R and
#   scripts/00_data_generation/fig7/03_hkratio_factorization.R <gene> metagratio
# NRT2 here has no polar samples (factorized before they were added), unlike Fig 7.
# Original: gene_explorer scripts/figures/summary_phylofactorization.R
#           ("summary_phylofactorization_SIT-iron-metaGratio_NRT2-nitrate-metaGratio")
# Run time: ~1 min.
# Output: results/figures/supplementary/S9_SIT-NRT2_metaGratio_factorization.{pdf,png}

source('R/load.R')
suppressPackageStartupMessages({
    library(ggtree)
    library(patchwork)
})

s9_dir <- p('intermediate', 'supplementary', 'S9')
metag_ylab <- 'Summed metaT/metaG ratio'

panels <- list(
    SIT = factorization_panel(file.path(s9_dir, 'SIT_iron_metagratio_factorization'), 'sit-categorical-',
                              gene = 'SIT', variable = 'iron', counts = 'ratio_df_w_metadata.tsv', ylab = metag_ylab),
    NRT2 = factorization_panel(file.path(s9_dir, 'NRT2_nitrate_metagratio_factorization'), 'nrt2-categorical-',
                               gene = 'NRT2', variable = 'nitrate', counts = 'ratio_df_w_metadata.tsv', ylab = metag_ylab)
)

# Both trees share the taxonomy palette, so its legend is only drawn in the first one
figure <- (manuscript_row(panels$SIT, 'A', 'SIT placement factorization over iron concentration') /
           manuscript_row(panels$NRT2, 'B', 'NRT2 placement factorization over nitrate concentration',
                          taxonomy_legend = 'none')) &
    manuscript_theme()

dir.create(p('figures', 'supplementary'), showWarnings = FALSE)
save_figure(figure, 'supplementary/S9_SIT-NRT2_metaGratio_factorization', width = 14, height = 13, pdf_device = cairo_pdf)
