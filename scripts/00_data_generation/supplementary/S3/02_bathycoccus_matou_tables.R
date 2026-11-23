# Supp S3, step 02: MATOU metaT and metaG abundance of the Bathycoccus prasinos unigenes.
#
# 1. Unigenes: best blast hit(s) of each B. prasinos transcript (step 01) with query coverage >= 80%,
#    identity >= 92% and alignment > 200 nt; only those with a Viridiplantae MATOU taxonomy are kept
# 2. Their values in the MATOU metaT and metaG occurrence tables (normalised abundance per sample)
#
# Input:  data/intermediate/supplementary/S3/bathycoccus_prasinos_blast_92id.tsv (01), or the paper copy
#         in data/raw/quantification/statistics/
#         external$matou: MATOU-v1.5.taxonomy.tsv.gz, MATOU-v1.5.meta{T,G}.occurrences.gz
# Output: data/intermediate/supplementary/S3/bathycoccus_prasinos_{geneids,matou_metaT,matou_metaG}.tsv
#         (the copies used for the paper: data/raw/quantification/statistics/)
# Run time: ~1 h (streams the MATOU tables, 25 and 9 GB).
# Original: delving_into_florenciella src/eda/bathycoccus_metaG-MATOU_vs_metaT.R (the extraction is
#   commented out there after its run)
# Check (2026-11-20): step 1 gives the same 5284 unigenes as DF's bathycoccus_prasinos_geneids.tsv. The
#   paper's metaT / metaG tables (March 2026) were extracted with an earlier blast run than the surviving
#   _92id table (September 2026): 43 of their 5106 unigenes are not in the current selection, so a rerun
#   changes S3 slightly.
#
# Usage (from the project root): Rscript scripts/00_data_generation/supplementary/S3/02_bathycoccus_matou_tables.R

source('R/load.R')
library(arrow, warn.conflicts = FALSE)

out_dir <- p('intermediate', 'supplementary', 'S3')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
blast_file <- file.path(out_dir, 'bathycoccus_prasinos_blast_92id.tsv')
if (!file.exists(blast_file)) blast_file <- p('quant_stats', 'bathycoccus_prasinos_blast_92id.tsv')

# 1. Unigenes ---------------------------------------------------------------------------------------
geneids <- read_tsv(blast_file, show_col_types = FALSE,
                    col_names = c('qseqid', 'sseqid', 'pident', 'length', 'qlen', 'qcovs', 'evalue', 'bitscore')) %>%
    filter(qcovs >= 80, pident >= 92, length > 200) %>%
    group_by(qseqid) %>%
    top_n(n = 1, wt = pident) %>%          # ties kept
    pull(sseqid) %>%
    str_remove('MATOU-v1.5.')

bathy_geneids <- read_tsv(file.path(external$matou, 'MATOU-v1.5.taxonomy.tsv.gz'), show_col_types = FALSE) %>%
    filter(geneID %in% geneids, str_detect(taxLineage, 'Viridiplantae')) %>%
    pull(geneID)
write_lines(str_c('MATOU-v1.5.', bathy_geneids), file.path(out_dir, 'bathycoccus_prasinos_geneids.tsv'))

# 2. MATOU abundance ---------------------------------------------------------------------------------
matou_values <- function(file) {
    open_dataset(file.path(external$matou, file), format = 'tsv') %>%
        filter(geneid %in% as.numeric(bathy_geneids)) %>%
        collect()
}
write_tsv(matou_values('MATOU-v1.5.metaT.occurrences.gz'), file.path(out_dir, 'bathycoccus_prasinos_matou_metaT.tsv'))
write_tsv(matou_values('MATOU-v1.5.metaG.occurrences.gz'), file.path(out_dir, 'bathycoccus_prasinos_matou_metaG.tsv'))
