# Environmental data of the Tara samples: data/raw/env/tara_env.tsv.gz, one row per PANGAEA sample and
# downloaded metaT run (scripts/00_data_generation/00_env_table.R). Key columns:
#   sample_id_pangaea   PANGAEA sample id
#   matou_sample_code   MATOU v1.5 sample code (e.g. 4SUR1GGMM14filter), for the 441 non-polar MATOU samples
#   name_downloaded     name of the downloaded metaT run (Carradec et al. 2018 and Tara polar), NA otherwise
#   in_context_general  sample of the Tara 18S context table (the samples pooled per MATOU sample in Fig 6)
# Every sample has one row per run, so use distinct(sample_id_pangaea, .keep_all = TRUE) for one per sample.
# The types are guessed from the whole table, so every script parses it the same way.
read_tara_env <- function() {
    read_tsv(p('env', 'tara_env.tsv.gz'), guess_max = Inf, show_col_types = FALSE)
}
