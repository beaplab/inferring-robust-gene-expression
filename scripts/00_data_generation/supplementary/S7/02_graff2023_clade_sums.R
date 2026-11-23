# Supp S7, step 02: Flavodoxin metaT, metaG and metaT/HK ratio summed per Graff et al. 2023 clade and sample
# (the metaT/metaG ratio sums are those of Fig 2: data/intermediate/fig2/flavodoxin_graff2023_clade_ratios.tsv).
#
# As fig2/05: gappa edit extract splits each per-sample jplace into one jplace per clade of
# graff2023_flavodoxin_clades.tsv ("basal": outside the clades; "uncertain": spanning several, dropped)
# and the multiplicities of each clade jplace are summed. Joined with the iron of the sample; metaG sample
# codes differ from the metaT ones in the two digits before "filter", so they are matched without them.
#
# Input:  data/raw/placement/epang_multiplicities_metaT_rawcounts/Flavodoxin   metaT read counts
#         data/raw/placement/epang_multiplicities_metaG/Flavodoxin             metaG values
#         data/intermediate/supplementary/S7/Flavodoxin_hkratio_jplace (01; if missing, the copy used for
#         the paper, data/raw/placement/Flavodoxin_hkratios_jplace/Flavodoxin)
#         data/raw/gene_data/graff2023_flavodoxin_clades.tsv, data/intermediate/fig2/flavodoxin_samples_iron.tsv
# Output: data/intermediate/supplementary/S7/graff2023_{metaT,metaG,hkratio}_counts.tsv
#           (clade, count, sample, iron_darwin_5m)
#         data/intermediate/supplementary/S7/graff2023_hkratio_pergene.tsv.gz (sample, clade, geneid, count:
#           the HK-ratio multiplicities per unigene, read by the figure for the "same samples, metaG genes"
#           version of the HK ratio; one table instead of gappa's one file per sample and clade)
#         the clade jplaces and the per-file multiplicities go to graff2023_extraction/ (not kept)
# Run time: ~30 min. Needs gappa.
# Original: gene_explorer scripts/analysis/phyloplacement_flavodoxin_analysis.qmd (steps 4-5) and
#           phyloplacement_flavodoxin_hkratio_analysis.qmd (extraction of the HK-ratio jplaces)
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/supplementary/S7/02_graff2023_clade_sums.R [out_dir]

source('R/load.R')

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else p('intermediate', 'supplementary', 'S7')
work_dir <- file.path(out_dir, 'graff2023_extraction')
clade_file <- p('gene_data', 'graff2023_flavodoxin_clades.tsv')

hkratio_dir <- p('intermediate', 'supplementary', 'S7', 'Flavodoxin_hkratio_jplace')
if (!dir.exists(hkratio_dir)) hkratio_dir <- p('placement', 'Flavodoxin_hkratios_jplace', 'Flavodoxin')

iron <- read_tsv(p('intermediate', 'fig2', 'flavodoxin_samples_iron.tsv'), show_col_types = FALSE)

# jplace folder -> clade sums with iron; sample_regex: what to strip from the jplace sample code
clade_sums <- function(jplace_dir, name, sample_regex = 'filter.*', mult_dir = file.path(work_dir, str_c(name, '_multiplicities'))) {
    message('Clade sums: ', name)
    samples <- iron %>%
        mutate(sample = str_remove(`Sample Code`, sample_regex)) %>%
        select(sample, iron_darwin_5m) %>%
        distinct(sample, .keep_all = TRUE)
    gappa_clade_multiplicities(jplaces = list.files(jplace_dir, pattern = '\\.jplace\\.gz$', full.names = TRUE),
                               clade_file = clade_file,
                               jplace_dir = file.path(work_dir, str_c(name, '_jplace')),
                               mult_dir = mult_dir)
    sum_clade_multiplicities(mult_dir, sample_regex) %>%
        inner_join(samples, by = 'sample') %>%
        write_tsv(file.path(out_dir, str_glue('graff2023_{name}_counts.tsv')))
}

clade_sums(p('placement', 'epang_multiplicities_metaT_rawcounts', 'Flavodoxin'), 'metaT')
clade_sums(p('placement', 'epang_multiplicities_metaG', 'Flavodoxin'), 'metaG', sample_regex = '[0-9][0-9]filter.*')
clade_sums(hkratio_dir, 'hkratio')
read_clade_multiplicities(file.path(work_dir, 'hkratio_multiplicities')) %>%
    write_tsv(file.path(out_dir, 'graff2023_hkratio_pergene.tsv.gz'))
unlink(work_dir, recursive = TRUE)
