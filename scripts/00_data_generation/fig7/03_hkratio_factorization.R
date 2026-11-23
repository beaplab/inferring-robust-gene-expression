# Fig 7: placement factorization of the metaT/HK ratio jplaces of SIT (against low/high iron) and NRT2
# (against low/high nitrate), and the summed ratio per factorization clade and sample (the violins).
#
# 1. Placement factorization of the ratio jplaces (ratio_placement_factorization() in R/gappa.R: samples
#    with a category and a ratio jplace; pseudo-count scaled to the size of the ratios)
# 2. Clades: each reference sequence goes to the highest-numbered factor whose non-root side contains it
# 3. gappa edit extract of every ratio jplace into one jplace per clade ("basal": outside the clades;
#    "uncertain": spanning several, dropped), and the ratios summed per clade and sample
#
# Input:  data/intermediate/fig7/<gene>_hkratio_jplace (02; if missing, the copy used for the paper,
#         data/raw/placement/<gene>_hkratios_jplace/<gene>)
#         data/raw/placement/<gene>_jplace_raw/<gene> (metaT read counts, pseudo-count only)
#         data/intermediate/fig7/SIT_metadata_ironcat.tsv / NRT2_metadata_nitratecat.tsv and
#         <gene>_samples_env.tsv (01)
# Output: data/intermediate/fig7/SIT_iron_hkratio_factorization/, NRT2_nitrate_hkratio_factorization/
#         <gene>-categorical-* (gappa), taxaall.tsv (clades), pseudocount.txt, metadata_factorization.tsv,
#         ratio_df_w_metadata.tsv (clade, count, sample + environmental data; read by Fig 7)
#         (the per-clade jplaces and multiplicities go to clade_ratio_jplace/ and clade_ratio_multiplicities/)
# Run time: ~30 min per gene with 8 threads. Needs gappa.
# Original: gene_explorer scripts/analysis/ratio_placement_factorization.R (run as "SIT iron hkratio" and
#           "NRT2 nitrate hkratio"; and "... metagratio" for Supplementary Fig S9)
#
# With "metagratio" as second argument, the same on the metaT/metaG ratio jplaces (Supplementary Fig S9):
#   input  data/intermediate/supplementary/S9/<gene>_metagratio_jplace (S9/01; if missing, the paper copy
#          data/raw/placement/<gene>_metagratios_jplace/<gene>)
#   output data/intermediate/supplementary/S9/<gene>_<variable>_metagratio_factorization/
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/fig7/03_hkratio_factorization.R <SIT|NRT2> [hkratio|metagratio]

source('R/load.R')

args <- commandArgs(trailingOnly = TRUE)
gene <- if (length(args) >= 1) args[1] else 'SIT'
ratio <- if (length(args) >= 2) args[2] else 'hkratio'
stopifnot(ratio %in% c('hkratio', 'metagratio'))
analyses <- list(SIT = list(variable = 'iron', metadata = 'SIT_metadata_ironcat.tsv'),
                 NRT2 = list(variable = 'nitrate', metadata = 'NRT2_metadata_nitratecat.tsv'))
stopifnot(gene %in% names(analyses))
variable <- analyses[[gene]]$variable
catcol <- str_c(variable, 'cat')

work_dir <- if (ratio == 'hkratio') p('intermediate', 'fig7') else p('intermediate', 'supplementary', 'S9')
ratio_dir <- file.path(work_dir, str_glue('{gene}_{ratio}_jplace'))
if (!dir.exists(ratio_dir)) ratio_dir <- p('placement', str_glue('{gene}_{ratio}s_jplace'), gene)
out_dir <- file.path(work_dir, str_glue('{gene}_{variable}_{ratio}_factorization'))
prefix <- str_glue('{str_to_lower(gene)}-categorical-')

# 1. Placement factorization -----------------------------------------------------------------------
samples <- read_tsv(p('intermediate', 'fig7', analyses[[gene]]$metadata), show_col_types = FALSE) %>%
    filter(!is.na(.data[[catcol]]))
used <- ratio_placement_factorization(samples, catcol,
                                      ratio_dir = ratio_dir,
                                      metat_dir = p('placement', str_c(gene, '_jplace_raw'), gene),
                                      out_dir = out_dir, prefix = prefix)

# 2. Clades ----------------------------------------------------------------------------------------
clade_file <- file.path(out_dir, 'taxaall.tsv')
factorization_clades(file.path(out_dir, str_c(prefix, 'factor_taxa.csv')), clade_file)

# 3. Summed ratio per clade and sample, with the environmental data ------------------------------------
mult_dir <- gappa_clade_multiplicities(jplaces = used$ratio_jplace, clade_file = clade_file,
                                       jplace_dir = file.path(out_dir, 'clade_ratio_jplace'),
                                       mult_dir = file.path(out_dir, 'clade_ratio_multiplicities'))

env <- read_tsv(p('intermediate', 'fig7', str_c(gene, '_samples_env.tsv')), show_col_types = FALSE) %>%
    mutate(sample = str_remove(`Sample Code`, 'filter.*')) %>%
    select(sample, any_of(unique(c(env_variables[[variable]]$column, 'iron_darwin_5m', 'silicate_q2',
                                   'nitrate_nitrite_q2', 'polar')))) %>%
    distinct(sample, .keep_all = TRUE)

counts <- sum_clade_multiplicities(mult_dir) %>%
    inner_join(env, by = 'sample')
write_tsv(counts, file.path(out_dir, 'ratio_df_w_metadata.tsv'))
message('Count table: ', n_distinct(counts$sample), ' samples, clades ', str_flatten_comma(sort(unique(counts$clade))))
