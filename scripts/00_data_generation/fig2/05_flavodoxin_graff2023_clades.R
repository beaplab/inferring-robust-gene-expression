# Fig 2 (right): summed metaT/metaG ratio of the Flavodoxin unigenes per Graff et al. 2023 clade and
# sample (the violins).
#
# gappa edit extract splits every ratio jplace into one jplace per clade of the reference tree
# (graff2023_flavodoxin_clades.tsv: reference sequence -> clade1_dielcycle, clade2_ironlimiting,
# greenalgae_dinos); queries placed outside the clades go to "basal", those spanning several clades to
# "uncertain" (dropped). The multiplicities (ratios) of each clade jplace are then summed.
#
# Input:  data/intermediate/fig2/flavodoxin_metagratio_jplace (03; if missing, the copy used for the paper,
#         data/raw/placement/Flavodoxin_metagratios_jplace/Flavodoxin)
#         data/raw/gene_data/graff2023_flavodoxin_clades.tsv
#         data/intermediate/fig2/flavodoxin_samples_iron.tsv (02)
# Output: data/intermediate/fig2/flavodoxin_graff2023_clade_ratios.tsv (clade, count, sample, iron_darwin_5m)
# Run time: ~10 min. Needs gappa.
# Original: gene_explorer scripts/analysis/phyloplacement_flavodoxin_analysis.qmd (steps 4-5 of the
#           pipeline section, ratio jplaces; the same was done there for the metaT and metaG jplaces)
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig2/05_flavodoxin_graff2023_clades.R

source('R/load.R')

ratio_dir <- p('intermediate', 'fig2', 'flavodoxin_metagratio_jplace')
if (!dir.exists(ratio_dir)) ratio_dir <- p('placement', 'Flavodoxin_metagratios_jplace', 'Flavodoxin')
work_dir <- p('intermediate', 'fig2', 'flavodoxin_graff2023_extraction')

samples <- read_tsv(p('intermediate', 'fig2', 'flavodoxin_samples_iron.tsv'), show_col_types = FALSE) %>%
    mutate(sample = str_remove(`Sample Code`, 'filter.*')) %>%
    select(sample, iron_darwin_5m) %>%
    distinct(sample, .keep_all = TRUE)

mult_dir <- gappa_clade_multiplicities(jplaces = list.files(ratio_dir, pattern = '\\.jplace\\.gz$', full.names = TRUE),
                                       clade_file = p('gene_data', 'graff2023_flavodoxin_clades.tsv'),
                                       jplace_dir = file.path(work_dir, 'clade_jplace'),
                                       mult_dir = file.path(work_dir, 'clade_multiplicities'))

sum_clade_multiplicities(mult_dir) %>%
    inner_join(samples, by = 'sample') %>%
    write_tsv(p('intermediate', 'fig2', 'flavodoxin_graff2023_clade_ratios.tsv'))
