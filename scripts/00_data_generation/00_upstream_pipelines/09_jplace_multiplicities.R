# Upstream 09: one jplace per gene and MATOU sample, with the abundance of each unigene in that sample as
# its multiplicity.
#
# For each gene (epang/<gene>.epa_result.jplace.gz, every placed unigene) and each sample of a MATOU table:
#   1. gappa edit filter --keep-names: only the queries of the unigenes present in the sample
#      (all six frames of a unigene that were placed get its value)
#   2. gappa edit multiplicity: the table value as multiplicity
# Output files are named <sample>filter.jplace.gz (the "filter" comes from gappa edit filter and is part of
# the sample codes everywhere downstream).
#
# Three versions, one per table (08):
#   epang_multiplicities_metaT          metaT, normalised MATOU values (HK and non-HK tables)
#                                       -> fig5/01, fig6/01 (via 10), fig7/02
#   epang_multiplicities_metaT_rawcounts metaT read counts (all genes) -> fig2/02-03 (Flavodoxin),
#                                       and the SIT/NRT2 jplace_raw of fig7 (copies, pruned to the samples
#                                       with metadata, gene_explorer scripts/analysis/phyloplacement_sit.sh)
#   epang_multiplicities_metaG          metaG, normalised (non-HK genes; HK genes commented out in the
#                                       original) -> fig2/02-03
#
# Input:  data/raw/placement/epang/ (07), data/raw/placement/matou_tables/ (08)
# Output: data/raw/placement/epang_multiplicities_{metaT,metaT_rawcounts,metaG}/<gene>/<sample>filter.jplace.gz
#         (symlinks to phylotranscriptomics data/matou_search/phyloplace/)
# Run time: days for all genes (~580 samples x 322 genes x 2 gappa calls).
# Original: phylotranscriptomics scripts/matou_phyloplace/create_multiplicity_files.R. There, each version
#   was run by commenting the others out, the metaG block looped over an undefined `sitgenes` (meant:
#   the non-HK genes), and every call shared /tmp/multiplicities_sel.txt (unsafe in parallel). Here the
#   three versions are explicit and each call uses its own temporary files.
# Checked (2026-11-20): gene SIT, 2 samples, metaT: same multiplicities per query as the jplaces used.
#
# NOT RUN IN THIS PROJECT. If run, it writes into phylotranscriptomics through the symlinks.
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/00_upstream_pipelines/09_jplace_multiplicities.R [gene ...]  (default: all)

source('R/load.R')
library(arrow, warn.conflicts = FALSE)

nonhk_genes <- c('SIT', 'NRT2', 'FAD1', 'Flavodoxin')
tables <- p('placement', 'matou_tables')

# Query names of a jplace and the MATOU geneid of each (several frames per geneid are possible)
jplace_queries <- function(jplace) {
    tibble(namesjplace = unlist(map(read_jplace(jplace)$placements, ~unlist(.x$n %||% map(.x$nm, 1))))) %>%
        mutate(geneid = as.double(str_remove(str_remove(namesjplace, '_frame=.*'), 'MATOU-v1.5.'))) %>%
        distinct()
}

# Abundance of the given geneids in a MATOU table, one tibble per sample
table_values <- function(table_file, geneids) {
    open_dataset(table_file, format = 'csv') %>%
        filter(geneid %in% geneids) %>%
        select(geneid, samplename, value) %>%
        collect() %>%
        split(.$samplename)
}

# jplace of one sample: the gene jplace filtered to the sample's unigenes, with their values as multiplicities
sample_jplace <- function(gene_jplace, queries, values, out_dir) {
    sample <- values$samplename[1]
    mult <- values %>%
        inner_join(queries, by = 'geneid', relationship = 'many-to-many') %>%
        select(namesjplace, value)
    tmp <- tempfile('sample_jplace')
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE))
    writeLines(mult$namesjplace, file.path(tmp, 'keep.txt'))
    run_command(str_glue('gappa edit filter --jplace-path {gene_jplace} --out-dir {tmp} --file-prefix {sample} ',
                         '--keep-names {tmp}/keep.txt --threads 1 --allow-file-overwriting > /dev/null 2>&1'))
    write_jplace_multiplicities(file.path(tmp, str_c(sample, 'filter.jplace')), mult, out_dir)
}

gene_multiplicities <- function(gene, table_file, out_name) {
    gene_jplace <- p('placement', 'epang', str_c(gene, '.epa_result.jplace.gz'))
    queries <- jplace_queries(gene_jplace)
    out_dir <- p('placement', out_name, gene)
    walk(table_values(table_file, queries$geneid), ~sample_jplace(gene_jplace, queries, .x, out_dir))
}

genes <- commandArgs(trailingOnly = TRUE)
if (length(genes) == 0) genes <- str_remove(list.files(p('placement', 'epang'), pattern = '\\.jplace\\.gz$'),
                                            '\\.epa_result\\.jplace\\.gz$')

for (gene in genes) {
    message(gene)
    nonhk <- gene %in% nonhk_genes
    gene_multiplicities(gene, file.path(tables, if (nonhk) 'matou_metaT_table_seqs_nonhk.tsv.gz' else 'matou_metaT_table_seqs_hk.csv.gz'),
                        'epang_multiplicities_metaT')
    gene_multiplicities(gene, file.path(tables, 'matou_metaT_rawcounts.tsv.gz'), 'epang_multiplicities_metaT_rawcounts')
    if (nonhk) gene_multiplicities(gene, file.path(tables, 'matou_metaG_table_seqs_nonhk.tsv.gz'), 'epang_multiplicities_metaG')
}
