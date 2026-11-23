# Fig 2 (left): genes of Pelagomonas calceolata differentially expressed between low and high iron
# Tara samples, with fastEmu (Clausen & Willis; radEmu with a reference set of genes).
#
# 1. Samples: Tara Oceans metaT (dataset 2012_carradec_tara, paired-end quantification) where P. calceolata is well covered: > 20% of its genes
#    expressed (TPM > 0)
# 2. Iron category: Darwin model iron at 5 m, high_iron if > 5e-04, low_iron otherwise
# 3. Genes: total TPM > 10, plus every housekeeping (HK) gene, ferredoxin and flavodoxin
# 4. Reference set of the first fit (the genes assumed not to change, whose robust median sets the zero of
#    the estimates): every HK gene of P. calceolata (281)
# 5. fastEmu fit (estimates and confidence intervals for every gene); then the robust score test for the
#    genes whose interval excludes 0 and |estimate| >= 0.4, plus the three flavodoxins of Fig 2
#    (always_test), reusing the fitted B
#
# Input:  data/raw/quantification/alignment/TPM_Pelagomonas_calceolata_quant-over_2012_carradec_tara_paired-end.tsv
#         data/raw/env/tara_env.tsv.gz (iron per downloaded run, name_downloaded)
#         data/raw/genomic/functional_annotation/Pelagomonas_calceolata.emapper.annotations
#         data/raw/gene_data/linkage_eggnogdb/correspondence_annotations.tsv (eggNOG name -> HK gene)
# Output: <out_dir>/pelagomonas_fastemu_iron_fit.rds        first fit (all genes)
#         <out_dir>/pelagomonas_fastemu_iron_fit_coeffs.tsv its coefficients
#         <out_dir>/pelagomonas_fastemu_iron_robust.rds     robust test
#         <out_dir>/pelagomonas_fastemu_iron_coeffs.tsv     its coefficients (read by Fig 2)
# Run time: ~13 h (first fit ~1 h, robust test with the data-driven reference ~1 h, with the HK reference
#           ~11 h); ~28 GB of RAM. Needs fastEmu and radEmu.
# Original: delving_into_florenciella src/analysis/differential_expression.R (section "Categorical data
#           with most genes" and "Robust test"; the continuous-iron fit there is not used in the paper)
#
# NOTES
# - The submitted paper anchored the first fit on the 73 HK genes of selection_houseg_denoms_bycorr-07.tsv
#   (data/intermediate/fig2/fastemu_reference_comparison/paper/). That subset is no longer used: every HK gene
#   is. The two only shift all estimates by a constant (-0.07 log2), but the all-HK intervals are narrower, so
#   more genes pass to the robust test. 01b-01d compare the two (run with the earlier version of this
#   script, which took the reference set as an argument).
# - The table read by Fig 2 (data/intermediate/fig2/pelagomonas_fastemu_iron_coeffs.tsv) is a copy of
#   fastemu_reference_comparison/all_hk/pelagomonas_fastemu_iron_coeffs.tsv. That run (7,458 genes, FAU
#   included) predates always_test; the flavodoxin CAH0373297 was tested afterwards with the reuse below.
# - The model is fitted on TPM, not on read counts, as in the original.
# - The robust test is run without reference_set, as in the original, so fastEmu uses its default: a
#   data-driven reference set of the 30 genes (of all genes, here none of them HK) closest to 0 in the fitted
#   B, i.e. the genes that change least relative to the HK genes. robust_reference = same (or both) also runs
#   it against the HK reference set, written as pelagomonas_fastemu_iron_{robust,coeffs}_samereference.*
# - Reuse: the fit is reused if pelagomonas_fastemu_iron_fit.rds is already in out_dir (delete it to refit).
#   If a coefficients table of the robust test is already there, only the genes to test that have no p-value
#   in it are tested, with the same call (the same B gives the same data-driven reference set), and their
#   score_stat and pval are written into that table; this run's rds is not saved.
# - The original passed the TPM matrix (219 samples) and the environmental table (218 samples, in
#   another order) to fastEmuFit without aligning them. Here both are restricted to the shared samples
#   and put in the same order. The rerun (2026-11-20) gives the same estimates as the cached table
#   (largest difference 0.007 log2), so the order did not matter.
# - The paper's table tested every gene whose interval excludes 0, without the |estimate| >= 0.4 filter
#   written in the original script (and here).
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/fig2/01_pelagomonas_fastemu_iron.R [data_driven|same|both] [out_dir] [exclude.tsv]
#   (default: data_driven data/intermediate/fig2; exclude.tsv: genes (column category) left out of the
#   model, e.g. fastemu_reference_comparison/convergence_scan/excluded_genes.tsv)

source('R/load.R')
suppressPackageStartupMessages({
    library(radEmu)
    library(fastEmu)
})

min_breadth <- 0.20        # share of the genes with TPM > 0 for a sample to be kept
min_total_tpm <- 10        # genes with a lower total TPM over the samples are left out
iron_threshold <- 5e-04    # iron_darwin_5m above this: high_iron
robust_min_estimate <- 0.4 # minimum |estimate| for the robust test
# tested whatever their first fit: the flavodoxins highlighted in Fig 2
always_test <- c('CAH0373297', 'CAH0367870', 'CAH0367856')

args <- commandArgs(trailingOnly = TRUE)
robust_reference <- if (length(args) >= 1) args[1] else 'data_driven'
out_dir <- if (length(args) >= 2) args[2] else p('intermediate', 'fig2')
exclude <- if (length(args) >= 3) read_tsv(args[3], show_col_types = FALSE)$category else character(0)
stopifnot(robust_reference %in% c('data_driven', 'same', 'both'))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Samples ----------------------------------------------------------------------------------
tpm <- read_tsv(p('quant_alignment', 'TPM_Pelagomonas_calceolata_quant-over_2012_carradec_tara_paired-end.tsv'),
                show_col_types = FALSE)
breadth <- rowSums(tpm[, -1] > 0) / (ncol(tpm) - 1)
tpm <- tpm[breadth > min_breadth, ]

# 2. Iron category ----------------------------------------------------------------------------
env <- read_tara_env() %>%
    filter(name_downloaded %in% tpm$sample) %>%
    distinct(name_downloaded, .keep_all = TRUE) %>%
    mutate(iron_darwin_5m = ifelse(iron_darwin_5m > iron_threshold, 'high_iron', 'low_iron')) %>%
    column_to_rownames('name_downloaded')

# 3. Genes ------------------------------------------------------------------------------------
annotation <- read_annotation_housekeeping(p('annotation', 'Pelagomonas_calceolata.emapper.annotations')) %>%
    rename(gene_name = `#query`)
hk_genes <- annotation %>% filter(!is.na(houseg)) %>% pull(gene_name) %>% intersect(colnames(tpm))
ferredoxins_flavodoxins <- annotation %>%
    filter(str_detect(Description, '^Ferredoxin$|^Flavodoxin$')) %>%
    pull(gene_name)
abundant <- names(which(colSums(tpm[, -1]) > min_total_tpm))

Y <- tpm[, c('sample', unique(c(abundant, hk_genes, ferredoxins_flavodoxins)))] %>%
    column_to_rownames('sample')
# without the excluded genes, and without samples or genes with only zeros
Y <- Y[, !colnames(Y) %in% exclude]
Y <- Y[rowSums(Y) != 0, colSums(Y) != 0]
if (length(exclude) > 0) message(length(exclude), ' genes excluded')

# same samples, same order in Y and in the covariates
shared <- intersect(rownames(Y), rownames(env))
Y <- Y[shared, ]
env <- env[shared, , drop = FALSE]

# 4. Reference set: every HK gene -------------------------------------------------------------
reference_set <- intersect(hk_genes, colnames(Y))
message('Reference set: ', length(reference_set), ' HK genes; Y: ', nrow(Y), ' samples x ', ncol(Y), ' genes')

# 5. Fits -------------------------------------------------------------------------------------
fit_file <- file.path(out_dir, 'pelagomonas_fastemu_iron_fit.rds')
if (file.exists(fit_file)) {
    message('Reusing ', fit_file)
    fit <- read_rds(fit_file)
} else {
    fit <- fastEmuFit(reference_set = reference_set, formula = ~iron_darwin_5m, data = env,
                      Y = Y, run_score_tests = FALSE, verbose = TRUE)
    write_rds(fit, fit_file)
    write_tsv(fit$coef, file.path(out_dir, 'pelagomonas_fastemu_iron_fit_coeffs.tsv'))
}

to_test <- fit$coef %>%
    as_tibble() %>%
    filter(lower > 0 | upper < 0, abs(estimate) >= robust_min_estimate) %>%
    pull(category) %>%
    union(str_subset(colnames(Y), str_c(always_test, collapse = '|')))

# column 2 of the design matrix (after the intercept) is the iron category
robust_test <- function(which_reference) {
    suffix <- if (which_reference == 'same') '_samereference' else ''
    coef_file <- file.path(out_dir, str_c('pelagomonas_fastemu_iron_coeffs', suffix, '.tsv'))
    previous <- if (file.exists(coef_file)) read_tsv(coef_file, show_col_types = FALSE) else NULL
    genes <- if (is.null(previous)) to_test else intersect(to_test, previous$category[is.na(previous$pval)])
    message(which_reference, ': ', length(genes), ' genes to test',
            if (!is.null(previous)) str_c(' (the others already in ', coef_file, ')') else '')
    if (length(genes) == 0) return(invisible())
    robust <- fastEmuFit(reference_set = if (which_reference == 'same') reference_set else 'data_driven',
                         formula = ~iron_darwin_5m, data = env, Y = Y, B = fit$B,
                         test_kj = data.frame(k = 2, j = which(colnames(Y) %in% genes)))
    if (is.null(previous)) {
        write_rds(robust, file.path(out_dir, str_c('pelagomonas_fastemu_iron_robust', suffix, '.rds')))
        write_tsv(robust$coef, coef_file)
    } else {
        new <- as_tibble(robust$coef)
        # same genes and reference set as the table being completed
        stopifnot(identical(new$category, previous$category),
                  identical(new$reference_category, previous$reference_category))
        message('largest estimate difference to the previous table: ', signif(max(abs(new$estimate - previous$estimate)), 2))
        previous %>%
            rows_update(new %>% filter(category %in% genes) %>% select(category, score_stat, pval), by = 'category') %>%
            write_tsv(coef_file)
    }
}
robust_references <- if (robust_reference == 'both') c('data_driven', 'same') else robust_reference
invisible(parallel::mclapply(robust_references, robust_test, mc.cores = length(robust_references)))
