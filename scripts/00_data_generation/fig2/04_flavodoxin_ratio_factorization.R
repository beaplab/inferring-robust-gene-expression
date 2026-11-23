# Fig 2 (right): placement factorization of the Flavodoxin metaT/metaG ratio jplaces against the iron
# category (low/high). Fig 2 uses its tree, to outline the first factors over the Graff et al. 2023
# clades.
#
# gappa analyze placement-factorization finds, one at a time, the edges of the reference tree whose
# split of the placement mass (balance) best separates low- from high-iron samples. The pseudo-count is
# scaled to the size of the ratios (factorization_pseudocount() in R/gappa.R).
#
# Input:  data/intermediate/fig2/flavodoxin_metagratio_jplace (03; if missing, the copy used for the paper,
#         data/raw/placement/Flavodoxin_metagratios_jplace/Flavodoxin)
#         data/raw/placement/epang_multiplicities_metaT_rawcounts/Flavodoxin (pseudo-count only)
#         data/intermediate/fig2/flavodoxin_samples_iron.tsv (02)
# Output: data/intermediate/fig2/flavodoxin_ratio_factorization/
#         flavo-categorical-* (factors_tree.newick, factor_taxa.csv, ...), pseudocount.txt,
#         metadata_factorization.tsv, jplace_list.txt, placement-factorization.log
# Run time: ~30 min with 8 threads. Needs gappa.
# Original: gene_explorer scripts/analysis/phyloplacement_flavodoxin_analysis.qmd (steps 3-4 of the
#           pipeline section)
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig2/04_flavodoxin_ratio_factorization.R

source('R/load.R')

ratio_dir <- p('intermediate', 'fig2', 'flavodoxin_metagratio_jplace')
if (!dir.exists(ratio_dir)) ratio_dir <- p('placement', 'Flavodoxin_metagratios_jplace', 'Flavodoxin')

samples <- read_tsv(p('intermediate', 'fig2', 'flavodoxin_samples_iron.tsv'), show_col_types = FALSE) %>%
    filter(!is.na(ironcat)) %>%
    distinct(`Sample Code`, .keep_all = TRUE)

ratio_placement_factorization(samples, catcol = 'ironcat',
                              ratio_dir = ratio_dir,
                              metat_dir = p('placement', 'epang_multiplicities_metaT_rawcounts', 'Flavodoxin'),
                              out_dir = p('intermediate', 'fig2', 'flavodoxin_ratio_factorization'),
                              prefix = 'flavo-categorical-')
