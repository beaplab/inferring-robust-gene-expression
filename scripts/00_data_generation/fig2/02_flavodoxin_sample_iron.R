# Fig 2 (right) and Fig 7: iron of every MATOU sample used in the Flavodoxin, SIT and NRT2 analyses,
# polar ones included.
#
# Non-polar samples (stations 4-153) have a PANGAEA id (matou_sample_code of the env table), so they come
# with their PANGAEA metadata, iron included. Polar Circle samples (stations 155-210) have no MATOU code
# -> PANGAEA id, so they are added here: iron_darwin_5m is a model value at 5 m, one per station, so it is
# taken from the env table by station. Their size fraction and layer are parsed from the
# MATOU code: <station><layer><replicate><lower><lower><upper><upper><year>, e.g. 173 DCM 1 KK QQ 13.
# The polar samples added are those with both a metaT and a metaG Flavodoxin jplace (the only ones that
# can get a metaT/metaG ratio).
#
# Input:  data/raw/env/tara_env.tsv.gz
#         data/raw/placement/epang_multiplicities_metaT_rawcounts/Flavodoxin, epang_multiplicities_metaG/Flavodoxin
# Output: data/intermediate/fig2/flavodoxin_samples_iron.tsv (`Sample Code`, iron_darwin_5m, ironcat, ...)
# Run time: seconds.
# Original: gene_explorer scripts/figures/flavodoxin_iron_ratios.R (it took the polar samples from the
#           clade counts of the ratio jplaces; same 121 samples, plus 155DCM1KKQQ13, which has no ratio
#           jplace)
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig2/02_flavodoxin_sample_iron.R [out_dir]

source('R/load.R')

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else p('intermediate', 'fig2')

iron_threshold <- 5e-04

# size-fraction letters (each repeated twice in the code; Z = no upper limit)
size_letters <- c(G = 0.8, K = 3, M = 5, Q = 20, S = 180, U = 2000)
# layer -> PANGAEA env_feature; IZZ (station 168) is the integrated 0-100 m sample, a MIX one in the
# polar table
layer_features <- c(SUR = '[SRF] surface water layer (ENVO:00010504)',
                    DCM = '[DCM] deep chlorophyll maximum layer (ENVO:01000326)',
                    MES = '[MES] marine water layer (ENVO:01000295) within the mesopelagic zone (ENVO:00000213)',
                    MXL = '[MIX] marine epipelagic wind mixed layer (ENVO:01000061)',
                    IZZ = '[MIX] marine epipelagic wind mixed layer (ENVO:01000061)')

env <- read_tara_env()

nonpolar <- env %>%
    filter(!is.na(matou_sample_code)) %>%
    distinct(sample_id_pangaea, .keep_all = TRUE) %>%
    select(-in_context_general) %>%
    rename(`Sample Code` = matou_sample_code) %>%
    relocate(`Sample Code`, iron_darwin_5m) %>%
    mutate(sample = str_remove(`Sample Code`, 'filter'))

# one value per station
polar_iron <- env %>%
    distinct(station_label, iron_darwin_5m) %>%
    filter(!is.na(iron_darwin_5m))

pairs <- paired_jplaces(p('placement', 'epang_multiplicities_metaT_rawcounts', 'Flavodoxin'),
                        p('placement', 'epang_multiplicities_metaG', 'Flavodoxin'))

polar <- tibble(sample = str_remove(basename(pairs$metat_jplace), 'filter\\.jplace\\.gz$')) %>%
    anti_join(nonpolar, by = 'sample') %>%
    extract(sample, into = c('station', 'layer', 'lower', 'upper'),
            regex = '^(\\d+)([A-Z]{3})\\d([A-Z])\\3([A-Z])\\4\\d{2}$', remove = FALSE) %>%
    mutate(`Sample Code` = str_c(sample, 'filter'),
           station_label = str_c('TARA_', str_pad(station, 3, pad = '0')),
           env_feature = layer_features[layer],
           size_fraction_lower_threshold = sprintf('%.2f', size_letters[lower]),
           size_fraction_upper_threshold = ifelse(upper == 'Z',
                                                  str_c('>', size_fraction_lower_threshold),
                                                  sprintf('%.2f', size_letters[upper]))) %>%
    select(-station, -layer, -lower, -upper) %>%
    left_join(polar_iron, by = 'station_label', relationship = 'many-to-one')

samples <- bind_rows(nonpolar, polar) %>%
    mutate(ironcat = ifelse(iron_darwin_5m < iron_threshold, 'low_iron', 'high_iron'))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(samples, file.path(out_dir, 'flavodoxin_samples_iron.tsv'))
message(nrow(nonpolar), ' non-polar + ', nrow(polar), ' polar samples')
