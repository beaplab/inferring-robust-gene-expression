# Upstream 08: MATOU abundance tables restricted to the unigenes placed by phyloplace (07), so that 09 does
# not stream the full tables (25 GB metaT, 9 GB metaG) once per gene.
#
# The query names of the jplaces are <MATOU-v1.5.><geneid>_frame=<n>; the geneid is matched against:
#   MATOU-v1.5.metaT.occurrences.gz      metaT, normalised (geneid, samplename, value)
#   MATOU-v1.5.metaG.occurrences.gz      metaG, normalised (geneid, samplename, value)
#   MATOU-v1.5_rawcounts_metaT.tsv.gz    metaT read counts (geneID, geneLength, sampleName, readCount),
#                                        renamed here to geneid, geneLength, samplename, value
# HK and non-HK genes go to separate tables, as in the original (the non-HK ones were placed later).
#
# Input:  data/raw/placement/epang/<gene>.epa_result.jplace.gz (07); MATOU tables (external$matou)
# Output: data/raw/placement/matou_tables/ (symlinks to phylotranscriptomics data/matou_search/):
#   matou_metaT_table_seqs_hk.csv.gz, matou_metaG_table_seqs_hk.csv.gz,
#   matou_metaT_table_seqs_nonhk.tsv.gz, matou_metaG_table_seqs_nonhk.tsv.gz, matou_metaT_rawcounts.tsv.gz
#   (all CSV, whatever the extension: the names are those of the files used)
# Run time: hours (streams the MATOU tables with arrow).
# Original: phylotranscriptomics scripts/matou_phyloplace/helper_subsetting_matou_tables.R (the non-HK
#           part was commented out after its run; the raw-count columns were renamed outside the script)
#
# NOT RUN IN THIS PROJECT. If run, it overwrites the tables in phylotranscriptomics through the symlinks.
#
# Usage (from the project root): Rscript scripts/00_data_generation/00_upstream_pipelines/08_matou_subset_tables.R

source('R/load.R')
library(arrow, warn.conflicts = FALSE)

nonhk_genes <- c('SIT', 'NRT2', 'FAD1', 'Flavodoxin')
out_dir <- p('placement', 'matou_tables')

# MATOU geneids of the queries of a jplace
jplace_geneids <- function(jplace) {
    read_jplace(jplace)$placements %>%
        map(~unlist(.x$n %||% map(.x$nm, 1))) %>%
        unlist() %>%
        str_remove('MATOU-v1.5.') %>%
        str_remove('_frame=.*') %>%
        as.double() %>%
        unique()
}

jplaces <- list.files(p('placement', 'epang'), pattern = '\\.jplace\\.gz$', full.names = TRUE)
is_nonhk <- str_remove(basename(jplaces), '\\.epa_result\\.jplace\\.gz$') %in% nonhk_genes
geneids <- list(hk = unique(unlist(map(jplaces[!is_nonhk], jplace_geneids))),
                nonhk = unique(unlist(map(jplaces[is_nonhk], jplace_geneids))))

occurrences <- function(file) open_dataset(file.path(external$matou, file), format = 'tsv')

subset_table <- function(file, ids, out_file) {
    occurrences(file) %>%
        filter(geneid %in% ids) %>%
        write_csv_arrow(file.path(out_dir, out_file))
}

subset_table('MATOU-v1.5.metaT.occurrences.gz', geneids$hk, 'matou_metaT_table_seqs_hk.csv.gz')
subset_table('MATOU-v1.5.metaG.occurrences.gz', geneids$hk, 'matou_metaG_table_seqs_hk.csv.gz')
subset_table('MATOU-v1.5.metaT.occurrences.gz', geneids$nonhk, 'matou_metaT_table_seqs_nonhk.tsv.gz')
subset_table('MATOU-v1.5.metaG.occurrences.gz', geneids$nonhk, 'matou_metaG_table_seqs_nonhk.tsv.gz')

# read counts, all genes
occurrences('MATOU-v1.5_rawcounts_metaT.tsv.gz') %>%
    filter(geneID %in% c(geneids$hk, geneids$nonhk)) %>%
    select(geneid = geneID, geneLength, samplename = sampleName, value = readCount) %>%
    write_csv_arrow(file.path(out_dir, 'matou_metaT_rawcounts.tsv.gz'))
