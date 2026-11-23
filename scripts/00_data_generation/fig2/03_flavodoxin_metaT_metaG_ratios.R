# Fig 2 (right): per-unigene metaT/metaG ratio of the Flavodoxin unigenes (MATOU), written as new
# jplaces with the ratios as multiplicities, so the placements stay and only the weights change.
#
# For every sample with both a metaT and a metaG Flavodoxin jplace (R/ratios.R):
#   ratio = metaT reads / metaG value / 1e6   (metaT reads per metaG RPM; 0 without metaG)
# Samples whose ratios sum to <= 1 are dropped, and only the samples of 02_flavodoxin_sample_iron.R are
# kept, so that the placement factorization (04) gets exactly one metadata row per jplace.
#
# Input:  data/raw/placement/epang_multiplicities_metaT_rawcounts/Flavodoxin  metaT read counts
#         data/raw/placement/epang_multiplicities_metaG/Flavodoxin           metaG values
#         data/intermediate/fig2/flavodoxin_samples_iron.tsv                 (02)
# Output: data/intermediate/fig2/flavodoxin_metaT_metaG_pergene_counts.tsv.gz
#           (sample, geneid, count_metat, count_metag)
#         data/intermediate/fig2/flavodoxin_metagratio_jplace/<sample code>.jplace.gz
#           (the copy used for the paper is data/raw/placement/Flavodoxin_metagratios_jplace/Flavodoxin)
# Run time: ~10 min with 16 cores. Needs gappa.
# Original: gene_explorer scripts/analysis/generate_ratios_metaT_metaG_phylogenetictrees.R (Flavodoxin section)
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/fig2/03_flavodoxin_metaT_metaG_ratios.R [ncores]

source('R/load.R')
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
ncores <- if (length(args) >= 1) as.integer(args[1]) else 16

metat_dir <- p('placement', 'epang_multiplicities_metaT_rawcounts', 'Flavodoxin')
metag_dir <- p('placement', 'epang_multiplicities_metaG', 'Flavodoxin')
out_jplace_dir <- p('intermediate', 'fig2', 'flavodoxin_metagratio_jplace')

# 1. metaT and metaG per unigene --------------------------------------------------------------
pairs <- paired_jplaces(metat_dir, metag_dir)
counts <- mcmapply(metat_metag_counts, pairs$metat_jplace, pairs$metag_jplace,
                   SIMPLIFY = FALSE, mc.cores = ncores) %>%
    bind_rows()
write_tsv(counts, p('intermediate', 'fig2', 'flavodoxin_metaT_metaG_pergene_counts.tsv.gz'))

# 2. Ratios, only for the samples with iron -------------------------------------------------------
iron_samples <- read_tsv(p('intermediate', 'fig2', 'flavodoxin_samples_iron.tsv'), show_col_types = FALSE)$`Sample Code`
ratios <- metat_metag_ratios(counts) %>%
    filter(sample %in% iron_samples)

# 3. Ratio jplaces: the metaT jplace of the sample with the ratios as multiplicities ----------------
unlink(out_jplace_dir, recursive = TRUE)   # no jplace of an earlier run is left behind
statuses <- split(ratios, ratios$sample) %>%
    mclapply(function(df) write_jplace_multiplicities(file.path(metat_dir, str_c(df$sample[1], '.jplace.gz')),
                                                      select(df, geneid, ratio), out_jplace_dir),
             mc.cores = ncores)
failed <- keep(statuses, ~inherits(.x, 'try-error'))
message(length(statuses) - length(failed), ' ratio jplaces written to ', out_jplace_dir,
        if (length(failed) > 0) str_c(', ', length(failed), ' failed') else '')
