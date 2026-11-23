# Fig 7 (and Supplementary Fig S7 for Flavodoxin): metaT of each SIT / NRT2 / Flavodoxin unigene divided by the metaT of the housekeeping (HK) genes of the
# same lineage in the same sample, written as new jplaces with the ratios as multiplicities.
#
#   ratio(query q, sample s) = metaT(q, s) / GM_g[ HKexpr(g, s, T) ]
#
#   T       the deepest taxon in the taxonomic path of q (best hit of gappa examine assign) with a reliable
#           HK denominator in sample s (step 5); if none is reliable, its parent, up to Eukaryota
#   HKexpr  summed metaT of the queries of HK gene g whose best hit falls within T
#   GM_g    geometric mean over the HK genes that have references within T (so they can detect it) and
#           are expressed in T in s
#
# The SIT / NRT2 tree and the HK trees have different topologies, so their placements are linked through
# the taxonomy: both are assigned with gappa examine assign against references labelled with the same
# (EukProt-style) taxonomy.
#
# Steps
#   1. Reference taxonomy of the gene ("|" -> "_", as in the tree; NA references of Flavodoxin filled in)
#   2. gappa examine assign of all queries of the gene (best hit), once, on the jplace without multiplicities
#   3. HK expression per HK gene, sample and taxon (cumulative over the descendants of the taxon)
#   4. Which taxa each HK gene can detect (has references in)
#   5. Denominator per sample and taxon: geometric mean, valid if >= MIN_GENES HK genes expressed and
#      >= MIN_FRAC of those able to detect the taxon
#   6. metaT of the gene per query and sample (numerator)
#   7. Ratio, with the deepest valid denominator per query and sample
#   8. Ratio jplaces (gappa edit multiplicity)
#
# Units: numerator and denominator are the normalised MATOU metaT values (library size and gene length),
# so the library size cancels in the ratio.
#
# Input:  data/raw/placement/epang/<gene>.epa_result.jplace.gz                all queries of the gene
#         data/raw/placement/epang_multiplicities_metaT/<gene>/                per-sample metaT jplaces
#         data/raw/placement/hk_gene_taxonomy/<gene or HK gene>.taxonomy.txt   reference taxonomies
#         HK gappa assign summaries, one table per HK gene: data/intermediate/hk_taxassign_combined
#           (fig6/02_combine_assign_summaries.sh; if missing, the copy used for the paper,
#           data/raw/placement/hk_taxassign_combined)
#         data/intermediate/fig6/housekeeping_gene_filter_matou.tsv            HK gene list (all genes,
#           its keep_core filter is not applied)
# Output: <out>/<gene>_hk_normalisation/
#           <gene>.taxonomy_fixed.txt, assign/ (gappa), hk_taxon_expression.tsv.gz,
#           hk_denominators.tsv.gz, pergene_hkratios.tsv.gz (one row per query and sample)
#         <out>/<gene>_hkratio_jplace/<sample code>.jplace.gz
#           (the copy used for the paper is data/raw/placement/<gene>_hkratios_jplace/<gene>)
#         <out>: data/intermediate/fig7 for SIT and NRT2, data/intermediate/supplementary/S7 for Flavodoxin
# Run time: ~10 min with 16 cores. Needs gappa.
# Original: gene_explorer scripts/analysis/generate_ratios_metaT_HK.R
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/fig7/02_metaT_HK_ratios.R <SIT|NRT2|Flavodoxin> [ncores] [out]

source('R/load.R')
suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
gene <- if (length(args) >= 1) args[1] else 'SIT'
ncores <- if (length(args) >= 2) as.integer(args[2]) else 16
stopifnot(gene %in% c('SIT', 'NRT2', 'Flavodoxin'))
out_base <- if (length(args) >= 3) args[3] else if (gene == 'Flavodoxin') p('intermediate', 'supplementary', 'S7') else
    p('intermediate', 'fig7')

# Reliability of a denominator (step 5)
MIN_GENES <- 20
MIN_FRAC <- 0.5

# Paths ------------------------------------------------------------------------------------------
full_jplace <- p('placement', 'epang', str_c(gene, '.epa_result.jplace.gz'))
metat_dir <- p('placement', 'epang_multiplicities_metaT', gene)
taxonomy_dir <- p('placement', 'hk_gene_taxonomy')
hk_summary_dir <- p('intermediate', 'hk_taxassign_combined')
if (!dir.exists(hk_summary_dir)) hk_summary_dir <- p('placement', 'hk_taxassign_combined')
hk_filter_file <- p('intermediate', 'fig6', 'housekeeping_gene_filter_matou.tsv')

out_dir <- file.path(out_base, str_c(gene, '_hk_normalisation'))
jplace_out_dir <- file.path(out_base, str_c(gene, '_hkratio_jplace'))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
unlink(jplace_out_dir, recursive = TRUE)   # no jplace of an earlier run is left behind

# All ancestors of each taxonomic path, itself included, with their depth. Taxa are kept as full
# paths, so equal names in different lineages stay apart.
#   "Eukaryota;Sar;Stramenopiles" -> Eukaryota (1), Eukaryota;Sar (2), Eukaryota;Sar;Stramenopiles (3)
taxon_prefixes <- function(taxopaths) {
    taxopaths <- unique(taxopaths)
    parts <- strsplit(taxopaths, ';', fixed = TRUE)
    data.table(taxopath = rep(taxopaths, lengths(parts)),
               taxon = unlist(lapply(parts, function(x) vapply(seq_along(x), function(i)
                   paste(x[seq_len(i)], collapse = ';'), ''))),
               depth = unlist(lapply(parts, seq_along)))
}

# 1. Reference taxonomy ----------------------------------------------------------------------------
# Tree leaf labels cannot contain "|", so it became "_" in the tree; 12 SIT references keep it in the
# taxonomy file, and gappa assign stops on unlabelled leaves.
# Flavodoxin: 25 references have NA as taxonomy, the diatom flavodoxins added from Graff et al. 2023
# (Thalassiosira, Phaeodactylum, Fragilariopsis, Seminavis) and the Pelagomonas calceolata genome CDS
# (lcl_CAKKNE...). With NA, queries placed next to them get no lineage (precisely in the diatom clades of
# interest) and gappa examine assign crashes on the full jplace. Their taxopaths are taken from the same
# genera in the other references (Seminavis, absent there, as a raphid pennate like Phaeodactylum).
message('Step 1: reference taxonomy of ', gene)
diatoms <- 'Eukaryota;Diaphoretickes;Sar;Stramenopiles;Gyrista;Ochrophyta;BD clade;Diatomeae;Bacillariophytina'
na_fixes <- tribble(
    ~pattern,          ~taxopath_fixed,
    '^Thalassiosira',  str_c(diatoms, ';Mediophyceae;Thalassiosirophycidae;Thalassiosira'),
    '^Phaeodactylum',  str_c(diatoms, ';Bacillariophyceae;Bacillariophycidae;Phaeodactylum'),
    '^Fragilariopsis', str_c(diatoms, ';Bacillariophyceae;Bacillariophycidae;Fragilariopsis'),
    '^Seminavis',      str_c(diatoms, ';Bacillariophyceae;Bacillariophycidae;Seminavis'),
    '^lcl_CAKKNE',     'Eukaryota;Diaphoretickes;Sar;Stramenopiles;Gyrista;Ochrophyta;Pelagophyceae;Pelagomonadales;Pelagomonas')

ref_tax <- read_tsv(file.path(taxonomy_dir, str_c(gene, '.taxonomy.txt')), col_names = c('ref', 'taxopath'),
                    col_types = 'cc', na = character()) %>%
    mutate(ref = str_replace_all(ref, fixed('|'), '_'))
# SIT and NRT2 have no NA: the patterns only match Flavodoxin references
for (i in seq_len(nrow(na_fixes))) {
    fix <- ref_tax$taxopath == 'NA' & str_detect(ref_tax$ref, na_fixes$pattern[i])
    ref_tax$taxopath[fix] <- na_fixes$taxopath_fixed[i]
}
stopifnot(!anyDuplicated(ref_tax$ref), !any(ref_tax$taxopath %in% c('NA', '')))
taxonomy_fixed <- file.path(out_dir, str_c(gene, '.taxonomy_fixed.txt'))
write_tsv(ref_tax, taxonomy_fixed, col_names = FALSE)

# 2. Taxonomic assignment of the queries of the gene --------------------------------------------------
# The placement of a query does not depend on the sample (the per-sample jplaces are subsets of the
# full one), so one run on the full jplace assigns them all. --best-hit: taxopath of the edge with the
# highest LWR, as for the HK genes.
message('Step 2: gappa examine assign of the ', gene, ' queries')
assign_dir <- file.path(out_dir, 'assign')
dir.create(assign_dir, showWarnings = FALSE)
run_command(str_glue('gappa examine assign --jplace-path {full_jplace} --taxon-file {taxonomy_fixed} ',
                     '--out-dir {assign_dir} --per-query-results --best-hit --threads 1 ',
                     '--allow-file-overwriting > {assign_dir}/assign.log 2>&1'))

# queries on the root edge get an empty taxopath: their taxon is the root, Eukaryota
query_assign <- fread(file.path(assign_dir, 'per_query.tsv'), select = c('name', 'taxopath'), na.strings = NULL)
query_assign[taxopath == '', taxopath := 'Eukaryota']
stopifnot(!anyDuplicated(query_assign$name))

# candidate denominators: every taxon on the path of some query
query_prefixes <- taxon_prefixes(query_assign$taxopath)
needed_taxa <- unique(query_prefixes$taxon)

# 3. HK expression per HK gene, sample and taxon --------------------------------------------------------
# Expression of a taxon T = summed metaT of the HK queries whose best hit is T or a descendant. HK queries
# assigned above T (e.g. only to Sar for a diatom T) do not count for T.
message('Step 3: HK expression per taxon')
hk_genes <- read_tsv(hk_filter_file, show_col_types = FALSE)$houseg
hk_genes <- hk_genes[file.exists(file.path(hk_summary_dir, str_c(hk_genes, '_summary.tsv'))) &
                     file.exists(file.path(taxonomy_dir, str_c(hk_genes, '.taxonomy.txt')))]
message('  ', length(hk_genes), ' HK genes')

hk_taxon_expression <- function(hk_gene) {
    summ <- fread(file.path(hk_summary_dir, str_c(hk_gene, '_summary.tsv')),
                  select = c('sample', 'taxopath', 'expression'), na.strings = NULL)
    summ[taxopath == '', taxopath := 'Eukaryota']
    prefixes <- taxon_prefixes(summ$taxopath)[taxon %in% needed_taxa]
    summ[prefixes, on = 'taxopath', allow.cartesian = TRUE, nomatch = NULL][
        , .(expression = sum(expression)), by = .(sample, taxon)][
        , gene := hk_gene][]
}
hk_expr <- rbindlist(mclapply(hk_genes, hk_taxon_expression, mc.cores = ncores))
fwrite(hk_expr, file.path(out_dir, 'hk_taxon_expression.tsv.gz'), sep = '\t')

# 4. Which taxa each HK gene can detect -------------------------------------------------------------
# An HK gene without references in a lineage cannot place reads there: its zero is structural, not
# biological, so only the genes with references within T enter the denominator of T
message('Step 4: HK gene representation per taxon')
hk_representation <- rbindlist(map(hk_genes, function(hk_gene) {
    ref_paths <- fread(file.path(taxonomy_dir, str_c(hk_gene, '.taxonomy.txt')), header = FALSE,
                       col.names = c('ref', 'taxopath'), na.strings = NULL)$taxopath
    data.table(gene = hk_gene, taxon = intersect(taxon_prefixes(ref_paths)$taxon, needed_taxa))
}))
n_represented <- hk_representation[, .(n_rep = .N), by = taxon]

# 5. Denominator per sample and taxon ---------------------------------------------------------------
# Geometric mean, so every HK gene weighs the same whatever its expression level. Zeros are left out
# (log(0)); MIN_GENES and MIN_FRAC keep it from resting on the few genes that happen to be detected.
message('Step 5: denominators')
denominators <- hk_expr[hk_representation, on = .(gene, taxon), nomatch = NULL][
    expression > 0,
    .(n_expr = .N, gm_hk = exp(mean(log(expression)))), by = .(sample, taxon)][
    n_represented, on = 'taxon', nomatch = NULL][
    , `:=`(depth = lengths(strsplit(taxon, ';', fixed = TRUE)),
           valid = n_expr >= MIN_GENES & n_expr / n_rep >= MIN_FRAC)][]
fwrite(denominators, file.path(out_dir, 'hk_denominators.tsv.gz'), sep = '\t')

# 6. metaT of the gene per query and sample ------------------------------------------------------------
message('Step 6: ', gene, ' metaT')
gene_metat <- list.files(metat_dir, pattern = 'jplace.gz$', full.names = TRUE) %>%
    mclapply(gappa_multiplicities, mc.cores = ncores) %>%
    bind_rows() %>%
    rename(metaT = count) %>%
    as.data.table()
# jplace samples are "<code>filter", HK samples "<code>"
gene_metat[, sample := str_remove(jplace_sample, 'filter$')]

# 7. Ratio ----------------------------------------------------------------------------------------
# per query and sample, the deepest taxon of its path with a valid denominator in that sample: where a
# lineage is rare and the HK genes miss it, the query falls back to its parent
message('Step 7: ratios')
gene_metat <- query_assign[gene_metat, on = 'name']
stopifnot(!anyNA(gene_metat$taxopath))

chosen <- unique(gene_metat[, .(sample, taxopath)])[
    query_prefixes, on = 'taxopath', allow.cartesian = TRUE, nomatch = NULL][
    denominators[valid == TRUE, .(sample, taxon, gm_hk, n_expr, n_rep)], on = .(sample, taxon), nomatch = NULL][
    order(-depth), .SD[1], by = .(sample, taxopath)]
setnames(chosen, c('taxon', 'depth'), c('denominator_taxon', 'denominator_depth'))

pergene <- chosen[gene_metat, on = .(sample, taxopath)][
    , `:=`(query_depth = lengths(strsplit(taxopath, ';', fixed = TRUE)),
           ratio = metaT / gm_hk)][]

# samples without any valid denominator (no HK data, or too few HK genes even at Eukaryota)
no_denominator <- pergene[is.na(ratio), .(n_queries = .N), by = sample]
if (nrow(no_denominator) > 0) message('  ', sum(no_denominator$n_queries), ' queries in ', nrow(no_denominator),
                                      ' samples without a valid denominator, removed')
pergene <- pergene[!is.na(ratio)]
pergene[, .(jplace_sample, sample, name, taxopath, query_depth, denominator_taxon, denominator_depth,
            n_rep, n_expr, metaT, gm_hk, ratio)] %>%
    fwrite(file.path(out_dir, 'pergene_hkratios.tsv.gz'), sep = '\t')

# 8. Ratio jplaces: the metaT jplace of each sample with the ratios as multiplicities -------------------
message('Step 8: ratio jplaces')
statuses <- split(pergene, by = 'jplace_sample') %>%
    mclapply(function(df) write_jplace_multiplicities(file.path(metat_dir, str_c(df$jplace_sample[1], '.jplace.gz')),
                                                      df[, .(name, ratio)], jplace_out_dir),
             mc.cores = ncores)
failed <- keep(statuses, ~inherits(.x, 'try-error'))
message('  ', length(statuses) - length(failed), ' jplaces written to ', jplace_out_dir,
        if (length(failed) > 0) str_c(', ', length(failed), ' failed') else '')

# Share of the metaT normalised at each depth of the denominator taxon
pergene[, .(metaT = sum(metaT)), by = denominator_depth][
    , share := round(metaT / sum(metaT), 3)][order(denominator_depth)] %>%
    print()
