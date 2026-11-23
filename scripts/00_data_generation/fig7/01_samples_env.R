# Fig 7: environmental data and low/high category of the MATOU samples used in the SIT (iron) and NRT2
# (nitrate) analyses. The sample set is that of Fig 2 (563 samples, polar ones included).
#
# SIT:  iron (Darwin model, 5 m) for every sample (02_flavodoxin_sample_iron.R of Fig 2), and silicate
#       joined by PANGAEA sample id (not available for the polar samples).
#       Category: low_iron if iron_darwin_5m < 5e-04, high_iron otherwise; only samples with a raw SIT
#       metaT jplace.
# NRT2: nitrate + nitrite (nitrate_nitrite_q2) of the MATOU sample when there is one ("exact"); otherwise
#       (all polar samples have no PANGAEA id) the median of the PANGAEA nutrient values of the same
#       station and depth layer ("station_layer"). 95 polar and 3 non-polar samples get it that way, and
#       72 have none and are left out.
#       Category: low_nitrate if < 1 µmol/L, high_nitrate otherwise.
#
# Input:  data/intermediate/fig2/flavodoxin_samples_iron.tsv
#         data/raw/env/tara_env.tsv.gz   (silicate; nitrate per MATOU sample and per station and layer)
#         data/raw/placement/SIT_jplace_raw/SIT                       (samples with SIT metaT)
# Output: data/intermediate/fig7/SIT_samples_env.tsv, SIT_metadata_ironcat.tsv,
#         data/intermediate/fig7/NRT2_samples_env.tsv, NRT2_metadata_nitratecat.tsv
# Run time: seconds.
# Original: gene_explorer scripts/analysis/phyloplacement_sit_analysis.qmd (chunk category-metadata),
#           scripts/analysis/build_nrt2_nitrate_env.R,
#           scripts/analysis/phyloplacement_nrt2_analysis.qmd (chunk category-metadata)
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig7/01_samples_env.R [out_dir] [fig2_dir]

source('R/load.R')

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else p('intermediate', 'fig7')
fig2_dir <- if (length(args) >= 2) args[2] else p('intermediate', 'fig2')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Tara Polar Circle stations are numbered from 155 on
is_polar <- function(sample_code) as.integer(str_extract(sample_code, '^\\d+')) >= 155

# MATOU layer codes -> PANGAEA env_feature codes
layer_codes <- c(SUR = 'SRF', DCM = 'DCM', MES = 'MES', MXL = 'MIX', IZZ = 'MIX', ZZZ = 'ZZZ', OTH = 'OTH')

# SIT: iron and silicate ------------------------------------------------------------------------
# one row per PANGAEA sample
env <- read_tara_env() %>%
    distinct(sample_id_pangaea, .keep_all = TRUE)

silicate <- select(env, sample_id_pangaea, silicate_q2)

sit_env <- read_tsv(file.path(fig2_dir, 'flavodoxin_samples_iron.tsv'), show_col_types = FALSE, guess_max = 1e4) %>%
    select(`Sample Code`, sample_id_pangaea, iron_darwin_5m) %>%
    left_join(silicate, by = 'sample_id_pangaea') %>%
    mutate(polar = is_polar(`Sample Code`)) %>%
    select(-sample_id_pangaea)
write_tsv(sit_env, file.path(out_dir, 'SIT_samples_env.tsv'))

sit_env %>%
    mutate(ironcat = env_category(iron_darwin_5m, 'iron', env_variables$iron$threshold)) %>%
    filter(!is.na(ironcat),
           file.exists(p('placement', 'SIT_jplace_raw', 'SIT', str_c(`Sample Code`, '.jplace.gz')))) %>%
    select(`Sample Code`, ironcat) %>%
    write_tsv(file.path(out_dir, 'SIT_metadata_ironcat.tsv'))

# NRT2: nitrate -----------------------------------------------------------------------------------
samples <- sit_env %>%
    select(`Sample Code`, iron_darwin_5m, polar) %>%
    mutate(sample = str_remove(`Sample Code`, 'filter.*'),
           station = as.integer(str_extract(sample, '^[0-9]+')),
           layer = layer_codes[str_extract(sample, '^[0-9]+([A-Z]{3})', group = 1)])

exact <- env %>%
    filter(!is.na(matou_sample_code)) %>%
    transmute(sample = str_remove(matou_sample_code, 'filter.*'), exact = nitrate_nitrite_q2) %>%
    filter(!is.na(exact)) %>%
    distinct(sample, .keep_all = TRUE)

# PANGAEA nutrient values of every sample of the same station and layer
station_layer <- env %>%
    transmute(station = as.integer(str_extract(station_label, '[0-9]+')),
              layer = str_extract(env_feature, '^\\[([A-Z]+)\\]', group = 1),
              nitrate_nitrite_q2) %>%
    filter(!is.na(nitrate_nitrite_q2)) %>%
    group_by(station, layer) %>%
    summarize(station_layer = median(nitrate_nitrite_q2),
              n_values = n_distinct(nitrate_nitrite_q2), .groups = 'drop')

nrt2_env <- samples %>%
    left_join(exact, by = 'sample') %>%
    left_join(station_layer, by = c('station', 'layer')) %>%
    mutate(nitrate_source = case_when(!is.na(exact) ~ 'exact',
                                      !is.na(station_layer) ~ 'station_layer'),
           nitrate_nitrite_q2 = coalesce(exact, station_layer),
           nitrate_n_values = if_else(nitrate_source == 'station_layer', n_values, 1L)) %>%
    select(`Sample Code`, iron_darwin_5m, nitrate_nitrite_q2, nitrate_source, nitrate_n_values, polar)
write_tsv(nrt2_env, file.path(out_dir, 'NRT2_samples_env.tsv'))

nrt2_env %>%
    select(`Sample Code`, nitrate_nitrite_q2) %>%
    distinct(`Sample Code`, .keep_all = TRUE) %>%
    mutate(nitratecat = env_category(nitrate_nitrite_q2, 'nitrate', env_variables$nitrate$threshold)) %>%
    filter(!is.na(nitratecat)) %>%
    write_tsv(file.path(out_dir, 'NRT2_metadata_nitratecat.tsv'))

count(nrt2_env, polar, nitrate_source) %>% print()
