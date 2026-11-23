# Supp S6, step 01: taxonomic bias of the MATOU metaG relative to metaT, measured on the housekeeping (HK)
# genes.
#
# HK expression per gene copy should vary little between lineages, so the share of a lineage in HK metaG
# should follow its share in HK metaT. A lineage much rarer in metaG than in metaT is under-sampled by the
# metaG (e.g. a huge genome dilutes the reads over non-coding DNA), and its metaT/metaG ratios rest on few
# or missing metaG values.
#
#   0. Taxonomy of every HK unigene: gappa examine assign --best-hit on the full jplace of each HK gene
#      (placement does not depend on the sample: one run per gene), coarse lineage from the taxopath
#   1. Share of each lineage in HK metaG vs in HK metaT, per paired sample (and size fraction)
#   2. metaG detection of the HK unigenes expressed in metaT, by lineage and global metaT quintile
#   3. The same detection, and the metaT/metaG level by lineage, for Flavodoxin, SIT and NRT2 (in the
#      tables, not in the figure); a gene is skipped if its inputs are missing
# Metazoa are dropped from the HK tables (as in the HK placements, metazoa extracted). metaT and metaG
# samples are paired when their codes differ only in the last two digits; codes with more than one
# sample are skipped (448 pairs).
#
# Input:  data/raw/placement/epang/<HK gene>.epa_result.jplace.gz, data/raw/placement/hk_gene_taxonomy/
#         data/raw/placement/matou_tables/matou_meta{T,G}_table_seqs_hk.csv.gz
#         data/intermediate/fig6/housekeeping_gene_filter_matou.tsv (the 318 HK genes)
#         step 3: per-gene metaT/metaG counts (Flavodoxin: intermediate/fig2; SIT, NRT2: S9/01) and the
#         gappa assign of the gene's unigenes (Flavodoxin: S7/01; SIT, NRT2: fig7/02)
# Output: data/intermediate/supplementary/S6/hk_lineage_shares.tsv.gz, hk_detection.tsv, <gene>_detection.tsv
#         hk_assign/<HK gene>/per_query.tsv.gz (cache of step 0, 124 MB, not copied into the project)
# Run time: ~15 min with 16 cores; ~60 GB RAM (HK metaT table).
# Original: gene_explorer scripts/analysis/metag_taxonomic_bias.R (printed summary of the first run:
#           gene_explorer results/metag_taxonomic_bias.txt)
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/supplementary/S6/01_metag_taxonomic_bias.R [ncores] [out_dir]

source('R/load.R')
suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
ncores <- if (length(args) >= 1) as.integer(args[1]) else 16
outdir <- if (length(args) >= 2) args[2] else p('intermediate', 'supplementary', 'S6')
setDTthreads(ncores)
options(width = 200)

hk_jplace_dir <- p('placement', 'epang')
hk_taxonomy_dir <- p('placement', 'hk_gene_taxonomy')
matou_tables <- p('placement', 'matou_tables')
hk_filter_file <- p('intermediate', 'fig6', 'housekeeping_gene_filter_matou.tsv')
assign_dir <- file.path(outdir, 'hk_assign')
dir.create(assign_dir, showWarnings = FALSE, recursive = TRUE)

# step 3 inputs per target gene: per-gene counts and gappa assign of its unigenes
target_inputs <- list(
    Flavodoxin = c(counts = p('intermediate', 'fig2', 'flavodoxin_metaT_metaG_pergene_counts.tsv.gz'),
                   assign = p('intermediate', 'supplementary', 'S7', 'Flavodoxin_hk_normalisation', 'assign', 'per_query.tsv')),
    SIT = c(counts = p('intermediate', 'supplementary', 'S9', 'SIT_metaT_metaG_pergene_counts.tsv.gz'),
            assign = p('intermediate', 'fig7', 'SIT_hk_normalisation', 'assign', 'per_query.tsv')),
    NRT2 = c(counts = p('intermediate', 'supplementary', 'S9', 'NRT2_metaT_metaG_pergene_counts.tsv.gz'),
             assign = p('intermediate', 'fig7', 'NRT2_hk_normalisation', 'assign', 'per_query.tsv')))

# Helpers ------------------------------------------------------------------------------------------

# Coarse lineages from the EukProt-style taxopaths (first match wins); coarser than R/taxonomy_groups.R
lineage_of <- function(tp) {
    fcase(grepl('Metazoa', tp), 'Metazoa',
          grepl('core syndiniales', tp), 'Syndiniales',
          grepl('Dinoflagellata', tp), 'Dinoflagellata',
          grepl('Ciliophora', tp), 'Ciliophora',
          grepl('Alveolata', tp), 'Other Alveolata',
          grepl('Diatomeae', tp), 'Diatoms',
          grepl('Pelagophyceae', tp), 'Pelagophyceae',
          grepl('Ochrophyta', tp), 'Other Ochrophyta',
          grepl('Stramenopiles', tp), 'Other Stramenopiles',
          grepl('Rhizaria', tp), 'Rhizaria',
          grepl('Haptophyta', tp), 'Haptophyta',
          grepl('Chlorophyta', tp), 'Chlorophyta',
          grepl('Cryptophyceae|Cryptista', tp), 'Cryptophyta',
          grepl('Choanoflagellata', tp), 'Choanoflagellata',
          grepl('Nucletmycea', tp), 'Fungi',
          lengths(strsplit(tp, ';', fixed = TRUE)) <= 3, 'Unresolved (<= Sar level)',
          default = 'Other')
}
# "MATOU-v1.5.28796_frame=-1" -> 28796, the geneid of the MATOU tables
geneid_of <- function(name) as.integer(sub('^MATOU-v1\\.5\\.(\\d+)_frame=.*$', '\\1', name))
# "100DCM1GGMM14" -> "100DCM1GGMM": shared by the metaT and metaG sample of a filter
pair_code <- function(s) sub('\\d\\d(filter)?$', '', s)
# size fraction: last four letters of the pair code (GGMM 0.8-5, MMQQ 5-20, QQSS 20-180, ...)
fraction_of <- function(code) substr(code, nchar(code) - 3, nchar(code))
log2_ratio <- function(g, t) log2((g + 1e-6) / t)

# 0. HK unigene taxonomy -----------------------------------------------------------------------------
# "|" is replaced in the taxonomy as in the tree (leaf labels cannot have it)
message('Step 0: taxonomy of the HK unigenes')
hk_genes <- fread(hk_filter_file)$houseg

assign_gene <- function(g) {
    o <- file.path(assign_dir, g)
    if (file.exists(file.path(o, 'per_query.tsv.gz'))) return(0)
    dir.create(o, showWarnings = FALSE)
    tax <- file.path(o, 'tax.txt')
    writeLines(gsub('|', '_', readLines(file.path(hk_taxonomy_dir, paste0(g, '.taxonomy.txt'))), fixed = TRUE), tax)
    status <- system(sprintf(paste('gappa examine assign --jplace-path %s --taxon-file %s --out-dir %s',
                                   '--per-query-results --best-hit --threads 1 --allow-file-overwriting > %s 2>&1'),
                             file.path(hk_jplace_dir, paste0(g, '.epa_result.jplace.gz')), tax, o, file.path(o, 'assign.log')))
    if (status == 0) system(sprintf('gzip -f %s', file.path(o, 'per_query.tsv')))
    unlink(c(tax, file.path(o, c('profile.tsv', 'labelled_tree.newick'))))
    status
}
status <- unlist(mclapply(hk_genes, assign_gene, mc.cores = ncores))
if (any(status != 0)) message('  gappa assign failed for: ', paste(hk_genes[status != 0], collapse = ', '))

hkmap <- rbindlist(mclapply(hk_genes[status == 0], function(g)
    fread(file.path(assign_dir, g, 'per_query.tsv.gz'), select = c('name', 'taxopath'), na.strings = NULL)[, hk := g],
    mc.cores = ncores))
hkmap[, geneid := geneid_of(name)][, lineage := lineage_of(taxopath)]
hkmap <- unique(hkmap[, .(geneid, hk, lineage)], by = 'geneid')
message('  ', nrow(hkmap), ' HK unigenes with taxonomy')

# HK metaT and metaG, paired samples ---------------------------------------------------------------
message('Reading the HK metaT and metaG tables')
read_matou <- function(f) {
    d <- fread(cmd = paste('zcat', f), colClasses = c('integer', 'character', 'double'))
    d[hkmap, on = 'geneid', nomatch = NULL]
}
mt <- read_matou(file.path(matou_tables, 'matou_metaT_table_seqs_hk.csv.gz'))
mg <- read_matou(file.path(matou_tables, 'matou_metaG_table_seqs_hk.csv.gz'))

# samples whose pair code is unique (one metaT and one metaG sample per filter)
unambiguous <- function(d) {
    u <- unique(d[, .(samplename)])[, code := pair_code(samplename)]
    u[, if (.N == 1) .SD, by = code]
}
paired <- merge(unambiguous(mt), unambiguous(mg), by = 'code', suffixes = c('_t', '_g'))
message('  ', nrow(paired), ' paired metaT/metaG samples')
mt <- mt[paired[, .(samplename = samplename_t, code)], on = 'samplename', nomatch = NULL][lineage != 'Metazoa']
mg <- mg[paired[, .(samplename = samplename_g, code)], on = 'samplename', nomatch = NULL][lineage != 'Metazoa']

# 1. Share of each lineage in HK metaG vs metaT -------------------------------------------------------
message('Step 1: lineage shares')
shares <- merge(mt[, .(metaT = sum(value)), by = .(code, lineage)],
                mg[, .(metaG = sum(value)), by = .(code, lineage)], by = c('code', 'lineage'), all = TRUE)
shares[is.na(metaT), metaT := 0][is.na(metaG), metaG := 0]
shares[, `:=`(share_metaT = metaT / sum(metaT), share_metaG = metaG / sum(metaG)), by = code]
shares[, fraction := fraction_of(code)]
fwrite(shares, file.path(outdir, 'hk_lineage_shares.tsv.gz'), sep = '\t')

cat('\n== 1. HK genes: log2(share in metaG / share in metaT) per lineage (samples where the lineage is > 1% of HK metaT)\n')
print(shares[share_metaT > 0.01,
       .(n_samples = .N, med_share_metaT = round(median(share_metaT), 3), med_share_metaG = round(median(share_metaG), 3),
         med_log2 = round(median(log2_ratio(share_metaG, share_metaT)), 2),
         q25 = round(quantile(log2_ratio(share_metaG, share_metaT), .25), 2),
         q75 = round(quantile(log2_ratio(share_metaG, share_metaT), .75), 2)),
       by = lineage][order(med_log2)])

# 2. metaG detection of the HK unigenes ----------------------------------------------------------------
# Per unigene and sample with metaT: does the paired metaG sample have any value for it? Stratified by
# global metaT quintile, so lineages are compared at the same expression level.
message('Step 2: HK detection')
mt[mg, on = .(geneid, code), metaG := i.value]
mt[, detected := !is.na(metaG)]
mt[, metaT_bin := cut(log10(value), breaks = quantile(log10(value), 0:5 / 5), include.lowest = TRUE, labels = paste0('Q', 1:5))]
hk_detection <- merge(
    mt[, .(n = .N, detected = round(mean(detected), 3),
           med_log10_T_over_G = round(median(log10(value[detected] / metaG[detected])), 2)), by = lineage],
    dcast(mt[, .(d = round(mean(detected), 3)), by = .(lineage, metaT_bin)], lineage ~ metaT_bin, value.var = 'd'),
    by = 'lineage')
hk_allmed <- mt[detected == TRUE, median(log10(value / metaG))]
fwrite(hk_detection, file.path(outdir, 'hk_detection.tsv'), sep = '\t')

cat('\n== 2. HK unigenes x samples with metaT: fraction with metaG, overall and by metaT quintile (Q1 lowest)\n')
print(hk_detection[order(detected)])
rm(mt, mg); invisible(gc())

# 3. Target genes ---------------------------------------------------------------------------------------
# count_metat are raw metaT reads and count_metag MATOU metaG values (other units than the HK tables), so
# the metaT/metaG level is compared as the offset of each lineage from the median of its own gene
# (offset_gene), next to the same offset in the HK genes (offset_HK).
message('Step 3: Flavodoxin, SIT, NRT2')
for (g in names(target_inputs)) {
    files <- target_inputs[[g]]
    if (!all(file.exists(files))) {
        message('  ', g, ': inputs missing, skipped (', paste(files[!file.exists(files)], collapse = ', '), ')')
        next
    }
    d <- fread(files[['counts']])
    a <- fread(files[['assign']], select = c('name', 'taxopath'), na.strings = NULL)
    d <- d[a, on = c(geneid = 'name'), nomatch = NULL][, lineage := lineage_of(taxopath)][count_metat > 0]
    d[, detected := !is.na(count_metag) & count_metag > 0]
    d[, metaT_bin := paste0('Q', as.integer(cut(log10(count_metat), breaks = unique(quantile(log10(count_metat), 0:5 / 5)),
                                                 include.lowest = TRUE)))]
    gene_med <- d[detected == TRUE, median(log10(count_metat / count_metag))]
    res <- d[, .(n = .N, metaT_share = round(sum(count_metat) / d[, sum(count_metat)], 3), detected = round(mean(detected), 3),
                 offset_gene = round(median(log10(count_metat[detected] / count_metag[detected])) - gene_med, 2)), by = lineage]
    res <- merge(res, dcast(d[, .(v = round(mean(detected), 3)), by = .(lineage, metaT_bin)], lineage ~ metaT_bin, value.var = 'v'),
                 by = 'lineage')
    res <- merge(res, hk_detection[, .(lineage, offset_HK = round(med_log10_T_over_G - hk_allmed, 2))], by = 'lineage', all.x = TRUE)
    fwrite(res, file.path(outdir, paste0(g, '_detection.tsv')), sep = '\t')
}
