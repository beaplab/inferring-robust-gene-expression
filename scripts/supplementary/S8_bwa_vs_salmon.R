# Supplementary Figure S8: gene read counts from BWA-mem2 (alignment, filtered at >= 95% identity) vs
# Salmon (mapping) for the same transcriptome and metatranscriptomes, both from the nextflow
# quantification (scripts/00_data_generation/00_upstream_pipelines/01_quantification_nextflow.sh).
#   S8a: Pelagomonas calceolata, Tara Oceans samples of Carradec et al.
#   S8b: Florenciella parvula, Tara polar samples
# One panel per sample; only genes counted by both tools. Red: 1:1 line; blue: smooth fit (ggplot2
# default, gam). The samples were drawn at random among those where BWA detects more than a given
# fraction of the genes (25% for P. calceolata, 10% for F. parvula): about 200 and 76 samples qualify.
#
# Original: delving_into_florenciella src/eda/helper_comparison-btw-salmon-and-bwamem.R (compare_samples();
#   the Phaeocystis, crossmapping and "weird genes" parts of that script are not in the supplement).
#   The original drew the samples without a seed (20 for P. calceolata in the version of the paper,
#   9 in the current script), so the samples shown in the paper figures are listed here explicitly; the
#   script checks that they all pass the detection filter.
# Run time: ~4 min (reading the wide count tables).
# Output: results/figures/supplementary/S8a_bwa_vs_salmon_pelagomonas.{pdf,png},
#         results/figures/supplementary/S8b_bwa_vs_salmon_florenciella.{pdf,png}

source('R/load.R')

# samples shown in the paper figures
samples_pelagomonas <- c('038_0o8-5_DCM', '038_0o8-5_SUR', '046_0o8-5_SUR', '046_20-180_SUR', '052_0o8-5_SUR',
                         '065_0o8-5_SUR', '102_0o8-5_SUR', '102_0o8-inf_SUR', '110_0o8-5_DCM', '124_0o8-5_SUR',
                         '124_0o8-inf_SUR', '125_0o8-5_SUR', '128_5-20_SUR', '131_0o8-5_DCM', '137_20-180_DCM',
                         '137_20-180_SUR', '142_0o8-5_SUR', '146_0o8-5_SUR', '146_20-180_SUR', '149_20-180_SUR')
samples_florenciella <- c('ERR6349574', 'ERR6349592', 'ERR6349594', 'ERR6349666', 'ERR6349696', 'ERR6349699',
                          'ERR6349711', 'ERR6349734', 'ERR6349767')

# Numreads tables (sample x gene) of both tools for one transcriptome and dataset
read_counts <- function(dir, file, value) {
    read_tsv(file.path(dir, file), show_col_types = FALSE) %>%
        pivot_longer(-sample, names_to = 'gene', values_to = value)
}

compare_samples <- function(file, min_genes_detected, shown) {
    comparison <- read_counts(paths$quant_mapping, file, 'count_salmon') %>%
        left_join(read_counts(paths$quant_alignment, file, 'count_bwa'), by = c('sample', 'gene')) %>%
        filter(count_salmon != 0 & count_bwa != 0)

    # samples where BWA counts more than min_genes_detected of the genes
    n_genes <- n_distinct(comparison$gene)
    good_samples <- comparison %>%
        group_by(sample) %>%
        filter(sum(count_bwa > 0) > n_genes * min_genes_detected) %>%
        pull(sample) %>%
        unique()
    stopifnot(all(shown %in% good_samples))

    comparison %>%
        filter(sample %in% shown) %>%
        ggplot(aes(count_salmon, count_bwa)) +
        geom_abline(slope = 1, intercept = 0, color = 'red') +
        geom_point() +
        geom_smooth() +   # ggplot2 default, as the original: gam (y ~ s(x, bs = 'cs')) for >= 1000 points
        facet_wrap(~sample, scales = 'free') +
        xlab('Salmon counts') +
        ylab('BWA counts')
}

pelagomonas <- compare_samples('Numreads_Pelagomonas_calceolata_quant-over_2012_carradec_tara.tsv', 0.25, samples_pelagomonas)
florenciella <- compare_samples('Numreads_EP00618_Florenciella_parvula_quant-over_2021_tara_polar.tsv', 0.10, samples_florenciella)

dir.create(p('figures', 'supplementary'), showWarnings = FALSE)
save_figure(pelagomonas, 'supplementary/S8a_bwa_vs_salmon_pelagomonas', width = 11, height = 9)
save_figure(florenciella, 'supplementary/S8b_bwa_vs_salmon_florenciella', width = 7, height = 7)
