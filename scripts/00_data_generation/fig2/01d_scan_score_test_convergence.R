# Fig 2 (left): genes of Pelagomonas calceolata differentially expressed between low and high iron
# Tara samples, with fastEmu (Clausen & Willis; radEmu with a reference set of genes).
#
# 1. Samples: Tara Oceans metaT (dataset 2012_carradec_tara, paired-end quantification) where P. calceolata is well covered: > 20% of its genes
#    expressed (TPM > 0)
# 2. Iron category: Darwin model iron at 5 m, high_iron if > 5e-04, low_iron otherwise
# 3. Genes: total TPM > 10, plus every housekeeping (HK) gene, ferredoxin and flavodoxin
# 4. Reference set (the genes assumed not to change, against which the others are compared): the HK
#    genes of selection_houseg_denoms_bycorr-07.tsv (reference = selection, the paper) or every HK gene
#    (reference = all_hk)
# 5. fastEmu fit (estimates and confidence intervals for every gene); then the robust score test only
#    for the genes whose interval excludes 0 and |estimate| >= 0.4, reusing the fitted B
#
# Input:  data/raw/quantification/alignment/TPM_Pelagomonas_calceolata_quant-over_2012_carradec_tara_paired-end.tsv
#         data/raw/env/tara_env.tsv.gz (iron per downloaded run, name_downloaded)
#         data/raw/genomic/functional_annotation/Pelagomonas_calceolata.emapper.annotations
#         data/raw/gene_data/linkage_eggnogdb/correspondence_annotations.tsv (eggNOG name -> HK gene)
#         data/raw/gene_data/selection_houseg_denoms_bycorr-07.tsv (reference set)
# Output: data/intermediate/fig2/pelagomonas_fastemu_iron_fit.rds        first fit (all genes)
#         data/intermediate/fig2/pelagomonas_fastemu_iron_robust.rds     robust test
#         data/intermediate/fig2/pelagomonas_fastemu_iron_coeffs.tsv     its coefficients (read by Fig 2)
# Run time: hours. Needs fastEmu and radEmu.
# Original: delving_into_florenciella src/analysis/differential_expression.R (section "Categorical data
#           with most genes" and "Robust test"; the continuous-iron fit there is not used in the paper)
#
# NOTES
# - The model is fitted on TPM, not on read counts, as in the original.
# - The robust test is run without reference_set, as in the original, so fastEmu uses its default: a
#   data-driven reference set of 30 genes (none of them HK genes in the paper run). The HK reference set
#   thus only selects the genes to test, and the coefficients and p-values of the paper are relative to
#   the data-driven set. robust_reference = same (or both) also runs it against the HK reference set,
#   written as pelagomonas_fastemu_iron_{robust,coeffs}_samereference.*
# - The fit is reused if pelagomonas_fastemu_iron_fit.rds is already in out_dir (delete it to refit).
# - The original passed the TPM matrix (219 samples) and the environmental table (218 samples, in
#   another order) to fastEmuFit without aligning them. Here both are restricted to the shared samples
#   and put in the same order. The rerun (2026-11-20) gives the same estimates as the cached table
#   (largest difference 0.007 log2), so the order did not matter.
# - The cached table tested every gene whose interval excludes 0 (2,302 genes), without the
#   |estimate| >= 0.4 filter written in the original script (and here), which keeps 2,035 of them.
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/fig2/01_pelagomonas_fastemu_iron.R [selection|all_hk] [data_driven|same|both] [out_dir]
#   (default: selection data_driven data/intermediate/fig2, the paper run)

source('R/load.R')
suppressPackageStartupMessages({
    library(radEmu)
    library(fastEmu)
})

min_breadth <- 0.20        # share of the genes with TPM > 0 for a sample to be kept
min_total_tpm <- 10        # genes with a lower total TPM over the samples are left out
iron_threshold <- 5e-04    # iron_darwin_5m above this: high_iron
robust_min_estimate <- 0.4 # minimum |estimate| for the robust test

args <- commandArgs(trailingOnly = TRUE)
reference <- if (length(args) >= 1) args[1] else 'selection'
robust_reference <- if (length(args) >= 2) args[2] else 'data_driven'
out_dir <- if (length(args) >= 3) args[3] else p('intermediate', 'fig2')
stopifnot(reference %in% c('selection', 'all_hk'), robust_reference %in% c('data_driven', 'same', 'both'))
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
# without samples or genes with only zeros
Y <- Y[rowSums(Y) != 0, colSums(Y) != 0]

# same samples, same order in Y and in the covariates
shared <- intersect(rownames(Y), rownames(env))
Y <- Y[shared, ]
env <- env[shared, , drop = FALSE]

# 4. Reference set ----------------------------------------------------------------------------
reference_hk <- if (reference == 'selection') {
    read_tsv(p('gene_data', 'selection_houseg_denoms_bycorr-07.tsv'), show_col_types = FALSE)$houseg
} else {
    unique(na.omit(annotation$houseg))
}
reference_set <- annotation %>%
    filter(houseg %in% reference_hk) %>%
    pull(gene_name) %>%
    intersect(colnames(Y))
message('Reference set (', reference, '): ', length(reference_set), ' genes; Y: ', nrow(Y), ' samples x ',
        ncol(Y), ' genes')

# ---- scan: replay the robust score tests gene by gene, recording convergence messages -------------
library(parallel)
variant <- args[4]                       # data_driven | same
ncores <- as.integer(args[5])
scan_dir <- out_dir
cmp_dir <- p('intermediate', 'fig2', 'fastemu_reference_comparison', reference)
fit <- read_rds(file.path(cmp_dir, 'pelagomonas_fastemu_iron_fit.rds'))
message('fit class: ', paste(class(fit), collapse = ' '))
to_test <- fit$coef %>% as_tibble() %>% filter(lower > 0 | upper < 0, abs(estimate) >= robust_min_estimate) %>% pull(category)
js <- which(colnames(Y) %in% to_test)
if (variant == 'data_driven') {
    # the reference set the robust run chose (reference_category of its table), fitted once
    dd <- read_tsv(file.path(cmp_dir, 'pelagomonas_fastemu_iron_coeffs.tsv'), show_col_types = FALSE)
    ref_dd <- which(colnames(Y) %in% dd$category[dd$reference_category])
    message('data-driven reference: ', length(ref_dd), ' genes')
    fit <- fastEmuFit(reference_set = ref_dd, formula = ~iron_darwin_5m, data = env, Y = Y, B = fit$B,
                      run_score_tests = FALSE)
}
if (nzchar(Sys.getenv('SCAN_N'))) js <- head(js, as.integer(Sys.getenv('SCAN_N')))
message(length(js), ' genes to scan with ', ncores, ' cores')
chunks <- split(js, cut(seq_along(js), ncores * 4, labels = FALSE))
scan_chunk <- function(chunk) {
    log <- list(); current <- NA_integer_; t0 <- proc.time()[['elapsed']]; times <- c()
    note <- function(type, text) log[[length(log) + 1]] <<- tibble(j = current, type = type, text = text)
    res <- try(withCallingHandlers(
        fastEmuFit(formula = ~iron_darwin_5m, data = env, Y = Y, fitted_model = fit, fastEmu_refit = FALSE,
                   test_kj = data.frame(k = 2, j = chunk), verbose = TRUE),
        message = function(m) {
            msg <- conditionMessage(m)
            i <- str_match(msg, '^Running score test (\\d+) of')[, 2]
            if (!is.na(i)) {
                now <- proc.time()[['elapsed']]
                if (!is.na(current)) times[as.character(current)] <<- now - t0
                current <<- chunk[as.integer(i)]; t0 <<- now
            } else note('message', str_squish(msg))
            invokeRestart('muffleMessage')
        },
        warning = function(w) { note('warning', str_squish(conditionMessage(w))); invokeRestart('muffleWarning') }),
        silent = TRUE)
    if (!is.na(current)) times[as.character(current)] <- proc.time()[['elapsed']] - t0
    pv <- if (inherits(res, 'try-error')) tibble(j = chunk, pval = NA_real_, error = as.character(res)) else
        as_tibble(res$coef) %>% mutate(j = row_number()) %>% filter(j %in% chunk) %>% transmute(j, pval, error = NA_character_)
    list(log = bind_rows(log), times = tibble(j = as.integer(names(times)), seconds = unname(times)), pval = pv)
}
out <- mclapply(chunks, scan_chunk, mc.cores = ncores, mc.preschedule = FALSE)
logs <- bind_rows(map(out, 'log')); times <- bind_rows(map(out, 'times')); pvals <- bind_rows(map(out, 'pval'))
per_gene <- tibble(j = js, category = colnames(Y)[js]) %>%
    left_join(times, by = 'j') %>% left_join(pvals, by = 'j') %>%
    left_join(logs %>% group_by(j) %>% summarize(n_retry = sum(str_detect(text, 'retrying with smaller')),
                                                 n_sandwich = sum(str_detect(text, 'sandwich')),
                                                 failed = any(str_detect(text, 'failed to converge across')),
                                                 other = str_c(unique(text[!str_detect(text, 'retrying with smaller|sandwich|failed to converge across')]), collapse = ' | '),
                                                 .groups = 'drop'), by = 'j') %>%
    mutate(across(c(n_retry, n_sandwich), ~replace_na(.x, 0L)), failed = replace_na(failed, FALSE))
write_tsv(per_gene, file.path(scan_dir, str_glue('scan_{reference}_{variant}.tsv')))
write_tsv(logs, file.path(scan_dir, str_glue('scan_{reference}_{variant}_messages.tsv')))
message('done: ', sum(per_gene$n_retry > 0), ' genes with retries, ', sum(per_gene$failed), ' failed, ',
        sum(!is.na(per_gene$error)), ' errors; median ', round(median(per_gene$seconds, na.rm = TRUE), 1), ' s per test')
