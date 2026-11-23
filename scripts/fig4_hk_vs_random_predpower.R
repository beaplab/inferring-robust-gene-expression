# Figure 4: random subsets of housekeeping (HK) genes predict the abundance of Pelagomonas calceolata
# better than random subsets of all genes.
#
# For subsets of K_MAIN genes, drawn either from the expressed non-HK genes (mean TPM 1-2000) or from
# the HK genes, fit metaG relative abundance ~ geometric mean of the metaT TPM of the subset
# (Carradec et al. 2022 Tara samples) and compare the R² distributions with the R² of all HK genes.
#   A: histograms of the R² of 1000 subsets of each kind; red line = adjusted R² of all HK genes
#   B: example regressions (one random subset, one HK subset, all HK genes)
#
# Original: delving_into_florenciella src/eda/comparison_housekeeping_torandomgenes_pred_metaGab_carradec.R
# Run time: ~1 min. Output: results/figures/fig4_hk_vs_random_predpower_k60.{pdf,png}

source('R/load.R')
library(patchwork)

K_MAIN   <- 60     # genes per subset (95 in the first version of the figure)
N_DRAWS  <- 1000   # subsets per gene pool
MIN_TPM  <- 1      # expressed non-HK genes: mean TPM in [MIN_TPM, MAX_TPM]
MAX_TPM  <- 2000

col_random <- '#7EBDC2'
col_hk     <- '#231f20'
col_allhk  <- '#bb4430'

set.seed(4223)

# --- data ---------------------------------------------------------------------------------------
tpm_table <- read_tsv(p('quant_stats', 'Pelago_GeneExpression_TaraOceans_TPM.tab'), show_col_types = FALSE)

# Carradec et al. 2022 (Commun Biol), Supplementary Data 4, sheet 10: metaG relative abundance per sample
envdata <- readxl::read_xlsx(p('quant_stats', '42003_2022_3939_MOESM4_ESM.xlsx'), sheet = 10, skip = 1) %>%
    janitor::clean_names()

annotation <- read_annotation_housekeeping(p('annotation', 'Pelagomonas_calceolata_RCC100.emapper.annotations'))
housekeeping_vars <- housekeeping_ids(annotation)

# expressed non-HK genes
tpm_no_houseg <- tpm_table %>% filter(!Gene %in% housekeeping_vars)
mean_tpm <- rowMeans(tpm_no_houseg[, -1])
genes_good_ab <- tpm_no_houseg$Gene[mean_tpm >= MIN_TPM & mean_tpm <= MAX_TPM]

random_mat <- expression_to_sample_matrix(tpm_table, genes_good_ab, trim_sample_suffix = 2)
hk_mat     <- expression_to_sample_matrix(tpm_table, housekeeping_vars, trim_sample_suffix = 2)

envdata_filt <- envdata %>%
    select(sample_name_2, relative_abundance_percent) %>%
    filter(sample_name_2 %in% rownames(random_mat))

random_mat <- random_mat[envdata_filt$sample_name_2, ]
hk_mat     <- hk_mat[envdata_filt$sample_name_2, ]
Y <- envdata_filt$relative_abundance_percent

# --- random subsets ---------------------------------------------------------------------------------
atrandom     <- sample_geomean_subsets(random_mat, Y, k = K_MAIN, n_draws = N_DRAWS)
housekeeping <- sample_geomean_subsets(hk_mat, Y, k = K_MAIN, n_draws = N_DRAWS)
results_allhouseg <- fit_geomean_model(colnames(hk_mat), hk_mat, Y)

allrandcompar <- bind_rows(random = atrandom, housekeeping = housekeeping, .id = 'comparison')

write_tsv(allrandcompar, p('tables', str_c('fig4_subsets_k', K_MAIN, '.tsv')))

# --- panel A: R² distributions ---------------------------------------------------------------------
# NOTE: as in the original figure, the histograms show r2 while the axis says adjusted R² and the
# red line is the adjusted R² of all HK genes (R² - adjusted R² = 0.004-0.006 here).
plot_comparison <- allrandcompar %>%
    ggplot(aes(r2)) +
    geom_histogram(aes(fill = comparison), bins = 30, show.legend = FALSE) +
    annotate(geom = 'text', x = if (K_MAIN == 95) 0.15 else 0.17, y = if (K_MAIN == 95) 150 else 80,
             label = str_wrap('Subsampling all genes', width = 10), color = col_random) +
    annotate(geom = 'text', x = if (K_MAIN == 95) 0.22 else 0.25, y = if (K_MAIN == 95) 130 else 165,
             label = str_wrap('Subsampling housekeeping genes', width = 20), color = col_hk) +
    annotate(geom = 'text', x = 0.33, y = 170,
             label = str_wrap(sprintf('All %d housekeeping genes', ncol(hk_mat)), width = 7),
             color = col_allhk) +
    geom_vline(xintercept = results_allhouseg$adj_r2, color = col_allhk, linewidth = 2) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 10))) +
    scale_fill_manual(values = c(housekeeping = col_hk, random = col_random)) +
    ylab('Count') +
    xlab(expression(paste('Adjusted ', R^{2}))) +
    ggtitle(str_wrap('A random housekeeping gene sampling presents higher predictive power'))

# --- panel B: example regressions (5th subset of each kind) ---------------------------------------
subset_genes <- function(df, i) str_split(df$genes[i], ', ')[[1]]

section2 <- ((plot_geomean_relationship(subset_genes(atrandom, 5), random_mat, Y, col_random) +
                  theme(axis.title.x = element_blank()) +
                  plot_geomean_relationship(subset_genes(housekeeping, 5), hk_mat, Y, col_hk) +
                  theme(axis.title.x = element_blank()) +
                  plot_layout(axis_titles = 'collect')) /
                 plot_geomean_relationship(colnames(hk_mat), hk_mat, Y, col_allhk)) +
    plot_layout(heights = c(0.3, 0.7), axis_titles = 'collect')

final <- (plot_comparison + section2) + plot_layout(widths = c(0.65, 0.35))

save_figure(final, str_c('fig4_hk_vs_random_predpower_k', K_MAIN), width = 12, height = 7)
