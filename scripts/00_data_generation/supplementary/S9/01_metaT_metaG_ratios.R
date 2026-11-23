# Supp S9, step 01: per-unigene metaT/metaG ratio of the SIT or NRT2 unigenes (MATOU), written as new
# jplaces with the ratios as multiplicities (as fig2/03 for Flavodoxin).
#
# For every sample with both a metaT and a metaG jplace of the gene (R/ratios.R):
#   ratio = metaT reads / metaG value / 1e6   (metaT reads per metaG RPM; 0 without metaG)
# Samples whose ratios sum to <= 1 are dropped (gappa's placement factorization refuses them).
#
# Input:  data/raw/placement/<gene>_jplace_raw/<gene>             metaT read counts (the
#         epang_multiplicities_metaT_rawcounts jplaces of the samples with metadata)
#         data/raw/placement/epang_multiplicities_metaG/<gene>     metaG values
# Output: data/intermediate/supplementary/S9/<gene>_metaT_metaG_pergene_counts.tsv.gz
#           (sample, geneid, count_metat, count_metag; also used by S6/01 step 3)
#         data/intermediate/supplementary/S9/<gene>_metagratio_jplace/<sample code>.jplace.gz
#           (the copy used for the paper: data/raw/placement/<gene>_metagratios_jplace/<gene>)
# Next:   Rscript scripts/00_data_generation/fig7/03_hkratio_factorization.R <gene> metagratio
# Run time: ~10 min per gene with 16 cores. Needs gappa.
# Original: gene_explorer scripts/analysis/generate_ratios_metaT_metaG_phylogenetictrees.R (SIT and NRT2
#           sections, output ratios_jplace_scaled)
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/supplementary/S9/01_metaT_metaG_ratios.R <SIT|NRT2> [ncores] [out_dir]

source('R/load.R')
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
gene <- if (length(args) >= 1) args[1] else 'SIT'
ncores <- if (length(args) >= 2) as.integer(args[2]) else 16
out_dir <- if (length(args) >= 3) args[3] else p('intermediate', 'supplementary', 'S9')
stopifnot(gene %in% c('SIT', 'NRT2'))

metat_dir <- p('placement', str_c(gene, '_jplace_raw'), gene)
metag_dir <- p('placement', 'epang_multiplicities_metaG', gene)
out_jplace_dir <- file.path(out_dir, str_c(gene, '_metagratio_jplace'))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. metaT and metaG per unigene --------------------------------------------------------------------
pairs <- paired_jplaces(metat_dir, metag_dir)
counts <- mcmapply(metat_metag_counts, pairs$metat_jplace, pairs$metag_jplace,
                   SIMPLIFY = FALSE, mc.cores = ncores) %>%
    bind_rows()
write_tsv(counts, file.path(out_dir, str_c(gene, '_metaT_metaG_pergene_counts.tsv.gz')))

# 2. Ratio jplaces: the metaT jplace of the sample with the ratios as multiplicities ---------------------
ratios <- metat_metag_ratios(counts)
unlink(out_jplace_dir, recursive = TRUE)   # no jplace of an earlier run is left behind
statuses <- split(ratios, ratios$sample) %>%
    mclapply(function(df) write_jplace_multiplicities(file.path(metat_dir, str_c(df$sample[1], '.jplace.gz')),
                                                      select(df, geneid, ratio), out_jplace_dir),
             mc.cores = ncores)
failed <- keep(statuses, ~inherits(.x, 'try-error'))
message(length(statuses) - length(failed), ' ratio jplaces written to ', out_jplace_dir,
        if (length(failed) > 0) str_c(', ', length(failed), ' failed') else '')
