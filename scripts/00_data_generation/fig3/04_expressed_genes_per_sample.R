# Fig 3E-F: number of expressed genes (TPM > 0) of each transcriptome / genome in each sample.
#
# Input:  data/raw/quantification/alignment/TPM_<taxon>_quant-over_<dataset>.tsv
#         (samples x genes, quantification pipeline)
# Output: data/intermediate/fig3/allstatistics_florenciella_plus_other.tsv.gz
#         (taxon, dataset, sample, occurrence_genes, max_TPM, median_TPM_over_zero)
#   NOTE: median_TPM_over_zero is median(tpm > 0), i.e. whether more than half of the genes are
#   expressed, not a median TPM (kept as in the original; the figure does not use it).
# Run time: ~5 min (reads every TPM table).
# Original: delving_into_florenciella src/figures/general_statistics_transcriptomes.R
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig3/04_expressed_genes_per_sample.R

source('R/load.R')

tpm_files <- list.files(p('quant_alignment'), pattern = '^TPM_.*\\.tsv$', full.names = TRUE)

sample_statistics <- function(path) {
    read_tsv(path, show_col_types = FALSE) %>%
        pivot_longer(names_to = 'gene', values_to = 'tpm', cols = -sample) %>%
        group_by(sample) %>%
        summarize(occurrence_genes = sum(tpm > 0),
                  max_TPM = max(tpm),
                  median_TPM_over_zero = median(tpm > 0, na.rm = TRUE))
}

tibble(path = tpm_files) %>%
    mutate(taxon = str_extract(basename(path), 'TPM_(.*)_quant-over', group = 1),
           dataset = str_extract(basename(path), 'quant-over_(.*).tsv', group = 1),
           statistics = map(path, sample_statistics)) %>%
    select(-path) %>%
    unnest(statistics) %>%
    write_tsv(p('intermediate', 'fig3', 'allstatistics_florenciella_plus_other.tsv.gz'))
