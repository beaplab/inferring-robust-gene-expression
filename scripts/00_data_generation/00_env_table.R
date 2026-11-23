# Environmental data: the single table of Tara sample metadata read by every step (R/io_env.R,
# read_tara_env()).
#
# One row per PANGAEA sample and downloaded metaT run (36,454 rows, 34,776 samples), with
#   - the PANGAEA water-column variables of the sample (Tara_Oceans_Pangaea_context.rds, $watercolumn:
#     iron_darwin_5m, station_label, env_feature, sample_material, ...) and three nutrients of $nut
#     (phosphate_q2, nitrate_nitrite_q2, silicate_q2), as joined by 02_formatting_Tara_forall_metaT_protist_samples.R
#     into envdata_Tara_All_downloaded_metaT_nisaba.tsv, which is taken as it is (text, so no value changes)
#     except event_time: there it was misparsed (2009-09-07T15:30 became 2020-09-09 07:15:30), so it is
#     taken back from $watercolumn
#   - name_downloaded: the downloaded run names of that table (Carradec et al. 2018 Tara runs, Tara polar
#     run accessions)
#   - matou_sample_code: MATOU v1.5 sample code of the 441 non-polar MATOU samples, from the MATOU sample
#     list (codes_matou_samples.xlsx, Tara Oceans ID -> Sample Code), with the "filter" of the jplace names.
#     Polar MATOU samples have no PANGAEA id (fig2/02 builds their rows from the code).
#   - in_context_general: the sample is in the Tara 18S context table (context_general.tsv, 1,341 samples);
#     Fig 6 (fig6/04) pools only these per MATOU sample
#
# It replaces six tables that were copied into data/raw/env/ and held the same values:
#   envdata_Tara_All_downloaded_metaT_nisaba.tsv, envdata_tara_polar_protist.tsv (iron of the polar
#   stations), matou_envdata_nonpolar.tsv and flavodoxin_samples_nonpolar_metadata.tsv (the MATOU samples),
#   Tara_Oceans_Pangaea_context.rds ($nut nitrate), context_general.tsv (sample_material)
#
# Input:  external$env_root (~/datasets/environmental_data_from_datasets):
#           tara_Alldata_including_polar_protists/Tara_Oceans_Pangaea_context.rds (PANGAEA context, from
#             N. Henry, https://gitlab.com/tara-and-friends-euk-metab/tara-oceans-metab-context)
#           tara_Alldata_including_polar_protists/envdata_Tara_All_downloaded_metaT_nisaba.tsv
#           tara_Alldata_including_polar_protists/codes_matou_samples.xlsx
#           tara_general_environmental_data/context_general.tsv
# Output: data/raw/env/tara_env.tsv.gz
# Run time: seconds.
# Original: tara_Alldata_including_polar_protists/02_formatting_Tara_forall_metaT_protist_samples.R and
#           03_formattint_Tara_for_matou1@5.R (in env_root)
#
# Usage (from the project root): Rscript scripts/00_data_generation/00_env_table.R [out_file]

source('R/load.R')

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1) args[1] else p('env', 'tara_env.tsv.gz')
tara_dir <- file.path(external$env_root, 'tara_Alldata_including_polar_protists')

# PANGAEA metadata + nutrients + run names, as text (4 rows have an embedded null in event_method, as
# in the original: readr warns and keeps the text before it)
env <- read_tsv(file.path(tara_dir, 'envdata_Tara_All_downloaded_metaT_nisaba.tsv'),
                col_types = cols(.default = col_character()))

# event_time as in PANGAEA
pangaea_time <- readRDS(file.path(tara_dir, 'Tara_Oceans_Pangaea_context.rds'))$watercolumn$values %>%
    as_tibble() %>%
    distinct(sample_id_pangaea, event_time)
stopifnot(!anyDuplicated(pangaea_time$sample_id_pangaea))

# MATOU sample codes
matou <- readxl::read_xlsx(file.path(tara_dir, 'codes_matou_samples.xlsx'), skip = 1) %>%
    transmute(sample_id_pangaea = `Tara Oceans ID`, matou_sample_code = str_c(`Sample Code`, 'filter'))
stopifnot(!anyDuplicated(matou$sample_id_pangaea), all(matou$sample_id_pangaea %in% env$sample_id_pangaea))

context_ids <- read_tsv(file.path(external$env_root, 'tara_general_environmental_data', 'context_general.tsv'),
                        col_select = sample_id_pangaea, show_col_types = FALSE)$sample_id_pangaea

env <- env %>%
    mutate(event_time = pangaea_time$event_time[match(sample_id_pangaea, pangaea_time$sample_id_pangaea)],
           in_context_general = sample_id_pangaea %in% context_ids) %>%
    left_join(matou, by = 'sample_id_pangaea') %>%
    relocate(sample_id_pangaea, matou_sample_code, name_downloaded, in_context_general)

write_tsv(env, out_file, na = '')
message(nrow(env), ' rows, ', n_distinct(env$sample_id_pangaea), ' PANGAEA samples, ',
        n_distinct(na.omit(env$matou_sample_code)), ' MATOU samples, ',
        n_distinct(env$sample_id_pangaea[env$in_context_general]), ' in context_general -> ', out_file)
