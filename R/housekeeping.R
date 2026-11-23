# Housekeeping (HK) genes in single-species transcriptomes, and the geometric mean of their
# expression as a proxy of the species abundance.

# eggNOG annotation of a transcriptome joined to the HK gene list.
# HK genes are the rows with a non-NA `label` (functional category); the link goes through the
# eggNOG Preferred_name (data/raw/gene_data/linkage_eggnogdb/correspondence_annotations.tsv).
read_annotation_housekeeping <- function(emapper_file) {
    corr_genes_to_eggnog <- read_tsv(p('gene_data', 'linkage_eggnogdb', 'correspondence_annotations.tsv'),
                                     show_col_types = FALSE) %>%
        distinct(Preferred_name, .keep_all = TRUE)

    read_tsv(emapper_file, comment = '##', show_col_types = FALSE) %>%
        left_join(corr_genes_to_eggnog, by = 'Preferred_name')
}

# Gene ids (first column of the annotation) of the HK genes
housekeeping_ids <- function(annotation) {
    annotation %>% filter(!is.na(label)) %>% pull(1)
}

# Row-wise geometric mean for a matrix of values (samples x genes)
geom_mean_rows <- function(mat, pseudocount = 1e-6) {
    mat <- mat + pseudocount
    exp(rowMeans(log(mat)))
}

# Fit lm(Y ~ geometric mean of `genes`) and return adjusted R², R² and AIC
fit_geomean_model <- function(genes, gene_mat, Y) {
    gm  <- geom_mean_rows(gene_mat[, genes, drop = FALSE])
    fit <- lm(Y ~ gm)
    s   <- summary(fit)
    list(adj_r2 = s$adj.r.squared,
         r2     = s$r.squared,
         aic    = AIC(fit))
}

# Draw `n_draws` random subsets of `k` genes (columns of gene_mat), fit the geometric-mean model
# for each and return one row per subset (k, genes, adj_r2, r2, aic).
# Subsets are drawn with replicate(sample()), so with the same seed the draws are the same as in the
# original scripts (delving_into_florenciella src/eda/comparison_housekeeping_torandomgenes_*.R).
sample_geomean_subsets <- function(gene_mat, Y, k, n_draws = 1000) {
    gene_names <- colnames(gene_mat)
    n_draws <- min(n_draws, choose(length(gene_names), k))
    combos <- replicate(n_draws, sample(gene_names, k, replace = FALSE), simplify = FALSE)
    combos <- unique(lapply(combos, sort))

    map_dfr(combos, function(genes) {
        scores <- fit_geomean_model(genes, gene_mat, Y)
        tibble(k = k,
               genes = paste(sort(genes), collapse = ', '),
               adj_r2 = scores$adj_r2,
               r2 = scores$r2,
               aic = scores$aic)
    })
}

# Scatter of the metaG abundance against the geometric mean of a set of genes, with a linear fit
plot_geomean_relationship <- function(genes, gene_mat, Y, colorline) {
    gm <- geom_mean_rows(gene_mat[, genes, drop = FALSE])

    tibble(gm, sample = names(gm), metag = Y) %>%
        ggplot(aes(gm, metag)) +
        geom_point(shape = 21) +
        geom_smooth(method = 'lm', formula = y ~ x, color = colorline) +
        ylab('metaG relative abundance (%)') +
        xlab('metaT geometric mean') +
        theme_classic()
}
