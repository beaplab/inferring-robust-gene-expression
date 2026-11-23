# Helper of 04: taxonomy file of a reference alignment, for gappa (reference name <tab> lineage).
#
# The lineage of each reference comes from its genus: the second field of the EukProt-style name
# (EP00003_Mantamonas_plastica_... -> Mantamonas), matched to the genus of the species in
# data/raw/gene_data/taxonomy_by_species.txt (first match kept). References whose genus is not there
# (e.g. sequences added from the literature) get NA.
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/00_upstream_pipelines/helpers/reference_taxonomy.R \
#       --fasta <alignment> --out <gene>.taxonomy.txt
# Original: gene_explorer scripts/data_generation/helper_create-housekeeping-taxonomy.R
# Checked (2026-11-20): from the current phyloplace MSAs it gives the NRT2 and Flavodoxin taxonomy files
# byte-identical; SIT differs only in 12 names with "|" (the SIT taxonomy was made in 2024 from an
# alignment whose names still had "|", later replaced by "_"; fig7/02 handles it).

suppressPackageStartupMessages({
    library(tidyverse)
    library(argparser)
})
source('config/paths.R')

args <- arg_parser('Taxonomy file of a reference alignment') %>%
    add_argument('--fasta', help = 'reference alignment (fasta)') %>%
    add_argument('--out', help = 'output taxonomy file') %>%
    parse_args()

species_lineages <- read_tsv(p('gene_data', 'taxonomy_by_species.txt'), col_names = c('species', 'taxrank'),
                             show_col_types = FALSE) %>%
    mutate(genus = str_remove(species, '_.*'))

# reference names: the fasta headers
tibble(names = str_remove(grep('^>', readLines(args$fasta), value = TRUE), '^>')) %>%
    mutate(genus = map_chr(names, ~str_split_1(.x, pattern = '_')[2])) %>%
    left_join(species_lineages, by = 'genus', relationship = 'many-to-many') %>%
    select(names, taxrank) %>%
    distinct(names, .keep_all = TRUE) %>%
    write_tsv(args$out, col_names = FALSE)
