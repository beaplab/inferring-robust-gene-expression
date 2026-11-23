# Supplementary Figures S1 and S2: housekeeping (HK) gene subsets as a proxy of the abundance of
# Pelagomonas calceolata (same data and model as Fig 4: metaG relative abundance ~ geometric mean of the
# metaT TPM of a gene subset, Carradec et al. 2022 Tara samples).
#
#   S2 (part A, run first): is the spread of the R² of HK subsets driven by functional categories?
#       5000 subsets of 95 HK genes; (a) R² histogram with the bottom / top 10%, (b) functional
#       composition of the genes drawn in low / mid / high subsets vs all HK genes. Permutation tests
#       and the variance of subset R² explained by category composition vs by gene identity.
#   S1 (part B): HK vs random subsets for k = 30..100 genes (1000 subsets per k and pool); (a) R²
#       distributions, (b) median R² ± IQR/2, (c) genes shared between two HK subsets (redundancy).
#
# Both parts are in one script because they share the random number stream (seed 4223, part A first):
# run separately, the subsets of S1 would not be those of the paper.
# R² = cor(Y, geomean)^2, identical to the R² of lm(Y ~ geomean) and much faster.
#
# Original: delving_into_florenciella src/eda/comparison_housekeeping_functionalmix_and_ksweep.R
# Run time: ~5 min.
# Output: results/figures/supplementary/S1_housekeeping_vs_random_ksweep.{pdf,png},
#         results/figures/supplementary/S2_housekeeping_functional_mix.{pdf,png},
#         results/tables/supplementary/housekeeping_{functional_mix_stats,functional_mix_pergene,vs_random_ksweep_stats}.tsv

source('R/load.R')
library(patchwork)
library(ggridges)

PSEUDOCOUNT <- 1e-6
MIN_TPM <- 1          # expressed non-HK genes: mean TPM in [MIN_TPM, MAX_TPM]
MAX_TPM <- 2000

K_A <- 95             # part A (S2): genes per subset, number of subsets
N_A <- 5000
K_VALUES <- seq(30, 100, by = 5)   # part B (S1): subset sizes, subsets per size and pool,
N_B      <- 1000                   # pairs of HK subsets for the redundancy
N_PAIRS  <- 5000

cols <- c(housekeeping = '#231f20', random = '#7EBDC2', all = '#bb4430')
fig_dir <- 'supplementary'
table_dir <- p('tables', 'supplementary')
dir.create(p('figures', fig_dir), showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(4223)

# --- data (as Fig 4) --------------------------------------------------------------------------------
tpm_table <- read_tsv(p('quant_stats', 'Pelago_GeneExpression_TaraOceans_TPM.tab'), show_col_types = FALSE)
envdata <- readxl::read_xlsx(p('quant_stats', '42003_2022_3939_MOESM4_ESM.xlsx'), sheet = 10, skip = 1) %>%
    janitor::clean_names()
annotation <- read_annotation_housekeeping(p('annotation', 'Pelagomonas_calceolata_RCC100.emapper.annotations'))
housekeeping_vars <- housekeeping_ids(annotation)

tpm_no_houseg <- tpm_table %>% filter(!Gene %in% housekeeping_vars)
mean_tpm <- rowMeans(tpm_no_houseg[, -1])
genes_good_ab <- tpm_no_houseg$Gene[mean_tpm >= MIN_TPM & mean_tpm <= MAX_TPM]

random_mat <- expression_to_sample_matrix(tpm_table, genes_good_ab, trim_sample_suffix = 2)
envdata_filt <- envdata %>%
    select(sample_name_2, relative_abundance_percent) %>%
    filter(sample_name_2 %in% rownames(random_mat))
Y <- envdata_filt$relative_abundance_percent

hk_mat <- expression_to_sample_matrix(tpm_table, housekeeping_vars, trim_sample_suffix = 2)
log_random <- log(as.matrix(random_mat[envdata_filt$sample_name_2, ]) + PSEUDOCOUNT)
log_hk     <- log(as.matrix(hk_mat[envdata_filt$sample_name_2, ]) + PSEUDOCOUNT)
n_hk <- ncol(log_hk)
lbl_expected <- sprintf('Expected (all %d)', n_hk)

# R² of Y ~ geometric mean of the columns idx of a log matrix
r2_subset <- function(log_mat, idx) cor(exp(rowMeans(log_mat[, idx, drop = FALSE])), Y)^2

# n random subsets of k columns: their indices and R²
sample_subsets <- function(log_mat, k, n) {
    idx <- replicate(n, sample.int(ncol(log_mat), k), simplify = FALSE)
    list(k = k, idx = idx, r2 = vapply(idx, function(i) r2_subset(log_mat, i), numeric(1)))
}

r2_all_hk <- r2_subset(log_hk, seq_len(n_hk))

# =====================================================================================================
# Part A (S2): functional composition of HK subsets
# =====================================================================================================

# functional category of each HK gene (categories with < 5 genes lumped into "Other")
gene_labels <- annotation %>%
    select(gene = 1, label) %>%
    filter(gene %in% colnames(log_hk)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    mutate(label = fct_lump_min(label, min = 5, other_level = 'Other'))
gene_labels <- gene_labels[match(colnames(log_hk), gene_labels$gene), ]
stopifnot(identical(gene_labels$gene, colnames(log_hk)))

hk_A <- sample_subsets(log_hk, K_A, N_A)

# subset x gene incidence (1: gene drawn in the subset)
incidence <- matrix(0L, nrow = N_A, ncol = n_hk, dimnames = list(NULL, colnames(log_hk)))
for (i in seq_len(N_A)) incidence[i, hk_A$idx[[i]]] <- 1L

# per gene: mean R² of the subsets with it minus without it, and its R² alone
r2_in  <- colSums(incidence * hk_A$r2) / colSums(incidence)
r2_out <- colSums((1 - incidence) * hk_A$r2) / colSums(1 - incidence)
single_r2 <- vapply(seq_len(n_hk), function(j) r2_subset(log_hk, j), numeric(1))
gene_effects <- gene_labels %>%
    mutate(delta_r2 = r2_in - r2_out, single_r2 = single_r2, n_in = colSums(incidence))

# category effect on the per-gene values: Kruskal-Wallis, permutation of the labels, eta²
between_var <- function(x, g) {
    m <- tapply(x, g, mean); n <- tapply(x, g, length)
    sum(n * (m - mean(x))^2)
}
perm_test_labels <- function(x, g, n_perm = 9999) {
    obs <- between_var(x, g)
    null <- replicate(n_perm, between_var(x, sample(g)))
    (sum(null >= obs) + 1) / (n_perm + 1)
}
kw_eta2 <- function(kw, n, k) max(0, (kw$statistic - k + 1) / (n - k))

n_cat <- nlevels(gene_effects$label)
kw_delta  <- kruskal.test(delta_r2 ~ label, data = gene_effects)
kw_single <- kruskal.test(single_r2 ~ label, data = gene_effects)
perm_delta  <- perm_test_labels(gene_effects$delta_r2, gene_effects$label)
perm_single <- perm_test_labels(gene_effects$single_r2, gene_effects$label)
eta_delta  <- kw_eta2(kw_delta, nrow(gene_effects), n_cat)
eta_single <- kw_eta2(kw_single, nrow(gene_effects), n_cat)

# composition of the bottom vs top 10% subsets: chi² against 999 permutations of the R²
cat_onehot <- model.matrix(~ label - 1, data = gene_effects)
colnames(cat_onehot) <- levels(gene_effects$label)
cat_counts <- incidence %*% cat_onehot   # subsets x categories
tail_chisq <- function(r2) {
    q <- quantile(r2, c(0.1, 0.9))
    low  <- colSums(cat_counts[r2 <= q[1], , drop = FALSE])
    high <- colSums(cat_counts[r2 >= q[2], , drop = FALSE])
    suppressWarnings(chisq.test(rbind(low, high))$statistic)
}
obs_tail  <- tail_chisq(hk_A$r2)
null_tail <- replicate(999, tail_chisq(sample(hk_A$r2)))
perm_tail <- (sum(null_tail >= obs_tail) + 1) / (999 + 1)

# share of the variance of subset R² explained by category composition vs by gene identity
var_expl_category <- summary(lm(hk_A$r2 ~ cat_counts))$adj.r.squared
var_expl_gene     <- summary(lm(hk_A$r2 ~ incidence))$adj.r.squared

q_A <- quantile(hk_A$r2, c(0.1, 0.9))
bin_A <- cut(hk_A$r2, breaks = c(-Inf, q_A, Inf), labels = c('Low (bottom 10%)', 'Mid', 'High (top 10%)'))
expected_prop <- table(gene_effects$label) / nrow(gene_effects)

composition <- map_dfr(levels(bin_A), function(b) {
    cnt <- colSums(cat_counts[bin_A == b, , drop = FALSE])
    tibble(bin = b, label = names(cnt), n = cnt, prop = cnt / sum(cnt))
}) %>%
    bind_rows(tibble(bin = lbl_expected, label = names(expected_prop),
                     n = as.numeric(table(gene_effects$label)), prop = as.numeric(expected_prop))) %>%
    mutate(bin = factor(bin, levels = c(lbl_expected, levels(bin_A))),
           label = factor(label, levels = levels(gene_effects$label)))

tail_oe <- composition %>%
    filter(bin != 'Mid') %>%
    select(bin, label, prop) %>%
    pivot_wider(names_from = bin, values_from = prop) %>%
    transmute(label,
              log2OE_low  = log2(`Low (bottom 10%)` / .data[[lbl_expected]]),
              log2OE_high = log2(`High (top 10%)` / .data[[lbl_expected]]))

gene_effects %>%
    group_by(label) %>%
    summarise(n_genes = n(), median_delta_r2 = median(delta_r2), median_single_r2 = median(single_r2),
              .groups = 'drop') %>%
    left_join(tail_oe, by = 'label') %>%
    mutate(kw_p_delta = kw_delta$p.value, perm_p_delta = perm_delta, eta2_delta = eta_delta,
           kw_p_single = kw_single$p.value, perm_p_single = perm_single, eta2_single = eta_single,
           tail_chisq = obs_tail, perm_p_tail = perm_tail,
           var_expl_category = var_expl_category, var_expl_gene = var_expl_gene) %>%
    write_tsv(file.path(table_dir, 'housekeeping_functional_mix_stats.tsv'))
write_tsv(gene_effects, file.path(table_dir, 'housekeeping_functional_mix_pergene.tsv'))

message(sprintf('S2: category ~ delta R2: KW p = %.3g, perm p = %.3g; top vs bottom 10%%: perm p = %.3g; ',
                kw_delta$p.value, perm_delta, perm_tail),
        sprintf('variance explained by category %.1f%%, by gene %.1f%%', 100 * var_expl_category, 100 * var_expl_gene))

cat_pal <- setNames(colorRampPalette(RColorBrewer::brewer.pal(12, 'Paired'))(n_cat), levels(gene_effects$label))

pA_hist <- tibble(r2 = hk_A$r2, bin = bin_A) %>%
    ggplot(aes(r2, fill = bin)) +
    geom_histogram(bins = 50, boundary = 0) +
    geom_vline(xintercept = q_A, linetype = 2) +
    scale_fill_manual(values = c('#7EBDC2', 'grey70', '#231f20'), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = expression(R^2), y = 'Count', title = sprintf('%d random subsets of %d housekeeping genes', N_A, K_A)) +
    theme_bw() + theme(legend.position = 'bottom')

pA_comp <- composition %>%
    ggplot(aes(bin, prop, fill = label)) +
    geom_col(color = 'white', linewidth = 0.2) +
    scale_fill_manual(values = cat_pal, name = 'Category') +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = NULL, y = 'Gene draws per category', subtitle = sprintf('Top vs bottom 10%%: perm. p = %.3g', perm_tail)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig_S2 <- (pA_hist / pA_comp + plot_layout(heights = c(0.45, 0.55))) +
    plot_annotation(tag_levels = 'a',
                    title = 'Functional composition of housekeeping subsets vs predictive power (R²)',
                    caption = sprintf(paste0('Variance of subset R² explained by functional-category composition: %.1f%%; ',
                                             'by individual gene identity: %.1f%%\n(adjusted R² of linear models on %d subsets).'),
                                      100 * var_expl_category, 100 * var_expl_gene, N_A))
save_figure(fig_S2, file.path(fig_dir, 'S2_housekeeping_functional_mix'), width = 9, height = 9, pdf_device = cairo_pdf)

# =====================================================================================================
# Part B (S1): HK vs random subsets across subset sizes
# =====================================================================================================

# fraction of genes shared by two subsets drawn independently (expected: k / pool size)
subset_redundancy <- function(idx, k, n_pairs = N_PAIRS) {
    a <- sample.int(length(idx), n_pairs, replace = TRUE)
    b <- sample.int(length(idx), n_pairs, replace = TRUE)
    keep <- a != b
    map2_dbl(a[keep], b[keep], ~length(intersect(idx[[.x]], idx[[.y]])) / k)
}

sweep_list <- map(K_VALUES, function(k) {
    hk  <- sample_subsets(log_hk, k, N_B)
    rnd <- sample_subsets(log_random, k, N_B)
    list(r2 = bind_rows(tibble(k = k, group = 'housekeeping', r2 = hk$r2),
                        tibble(k = k, group = 'random', r2 = rnd$r2)),
         redundancy = tibble(k = k, shared = subset_redundancy(hk$idx, k)))
})
sweep      <- map_dfr(sweep_list, 'r2')
redundancy <- map_dfr(sweep_list, 'redundancy')

stats_redundancy <- redundancy %>%
    group_by(k) %>%
    summarise(median_shared = median(shared), q25_shared = quantile(shared, 0.25),
              q75_shared = quantile(shared, 0.75), mean_shared = mean(shared),
              expected_shared = unique(k) / n_hk, .groups = 'drop')

stats_B <- sweep %>%
    group_by(k) %>%
    group_modify(function(d, key) {
        h <- d$r2[d$group == 'housekeeping']
        r <- d$r2[d$group == 'random']
        w  <- wilcox.test(h, r)
        ks <- suppressWarnings(ks.test(h, r))
        auc <- unname(w$statistic) / (length(h) * length(r))
        tibble(median_hk = median(h), median_random = median(r), median_diff = median(h) - median(r),
               sd_hk = sd(h), sd_random = sd(r), iqr_hk = IQR(h), iqr_random = IQR(r),
               auc = auc, cliffs_delta = 2 * auc - 1, ks_D = unname(ks$statistic), wilcox_p = w$p.value)
    }) %>%
    ungroup() %>%
    mutate(wilcox_p_adj = p.adjust(wilcox_p, method = 'BH')) %>%
    left_join(stats_redundancy, by = 'k')
write_tsv(stats_B, file.path(table_dir, 'housekeeping_vs_random_ksweep_stats.tsv'))

pB_ridges <- sweep %>%
    ggplot(aes(r2, factor(k), fill = group)) +
    geom_density_ridges(alpha = 0.6, scale = 1.6, color = 'white', linewidth = 0.2,
                        quantile_lines = TRUE, quantiles = 2) +
    geom_vline(xintercept = r2_all_hk, color = cols['all'], linewidth = 1) +
    annotate('text', x = r2_all_hk, y = Inf, vjust = 1.2, hjust = -0.05, size = 3,
             label = sprintf('All %d housekeeping', n_hk), color = cols['all']) +
    scale_fill_manual(values = cols[c('housekeeping', 'random')], name = NULL) +
    labs(x = expression(R^2), y = 'Genes per subset (k)') +
    theme_bw() + theme(legend.position = 'bottom')

pB_redundancy <- redundancy %>%
    ggplot(aes(k, shared)) +
    geom_boxplot(aes(group = k), fill = 'grey85', color = cols['housekeeping'],
                 outlier.size = 0.3, outlier.alpha = 0.3, width = 3) +
    geom_line(data = stats_redundancy, aes(k, expected_shared), color = cols['all'], linetype = 2) +
    annotate('text', x = min(K_VALUES), y = max(redundancy$shared), hjust = 0, vjust = 1, size = 3,
             color = cols['all'], label = sprintf('Expected: k / %d', n_hk)) +
    scale_x_continuous(breaks = K_VALUES) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = 'Genes per subset (k)', y = 'Genes shared between\ntwo housekeeping subsets',
         subtitle = sprintf('Housekeeping subset redundancy (%s pairs per k)', format(N_PAIRS, big.mark = ','))) +
    theme_bw()

pB_median <- stats_B %>%
    select(k, housekeeping = median_hk, random = median_random, iqr_hk, iqr_random) %>%
    pivot_longer(c(housekeeping, random), names_to = 'group', values_to = 'median') %>%
    mutate(iqr = ifelse(group == 'housekeeping', iqr_hk, iqr_random)) %>%
    ggplot(aes(k, median, color = group, fill = group)) +
    geom_ribbon(aes(ymin = median - iqr / 2, ymax = median + iqr / 2), alpha = 0.2, color = NA) +
    geom_line() + geom_point() +
    geom_hline(yintercept = r2_all_hk, color = cols['all'], linetype = 2) +
    scale_color_manual(values = cols[c('housekeeping', 'random')], name = NULL) +
    scale_fill_manual(values = cols[c('housekeeping', 'random')], name = NULL) +
    scale_x_continuous(breaks = K_VALUES) +
    labs(x = 'Genes per subset (k)', y = expression(paste('Median ', R^2, ' (± IQR/2)'))) +
    theme_bw() + theme(legend.position = 'none')

fig_S1 <- (pB_ridges | (pB_median / pB_redundancy)) +
    plot_layout(widths = c(0.55, 0.45)) +
    plot_annotation(tag_levels = 'a',
                    title = 'Housekeeping vs random gene subsets across subset sizes',
                    caption = str_wrap(paste0(
                        'Each k: ', N_B, ' subsets per group. Housekeeping subsets are drawn from ', n_hk, ' genes, ',
                        'random subsets from ', ncol(log_random), ' genes (mean TPM 1-2000); ',
                        'housekeeping variance partly shrinks with k because subsets overlap more.'), 150))
save_figure(fig_S1, file.path(fig_dir, 'S1_housekeeping_vs_random_ksweep'), width = 13, height = 8, pdf_device = cairo_pdf)
