# Supp S4, step 01: read placement vs sequence placement, per-gene caches.
#
# For 7 HK genes in both datasets, compares the placement on the same reference tree of
#   - reads:     metaT reads placed directly (Richter et al.; metazoa extracted, keeping the basal branches)
#   - sequences: MATOU unigenes, with the metaT raw counts of each unigene as multiplicity
# The slow part (gappa) runs here and leaves, per gene, an aggregated cache read by the figure script.
#
#   1. Samples: MATOU sample codes (<station><depth><replicate><fraction><year>filter) matched to the read
#      sample names (<station 3 digits>_<fraction>_<depth>); with two replicates / years, the first replicate
#   2. gappa edit extract: the metazoan placements are removed from the sequence jplaces, as for the reads
#      (only the basal branches kept; same command as 00_upstream_pipelines/10_extract_metazoa.sh)
#   3. gappa examine assign (per query) + gappa examine edpl, per sample and method: one row per read /
#      unigene with multiplicity, EDPL, best-hit LWR and taxopath, and the deepest taxopath with an
#      accumulated LWR >= CONF_THRESHOLD. A unigene placed in several frames keeps only the frame with the
#      highest best-hit LWR (otherwise its raw counts are counted twice).
#   4. gappa analyze krd: KR distance between all the jplaces of a gene (both methods)
# Read multiplicities are 0.5 for each mate of a pair, so both methods are in read pairs (assuming the
# MATOU raw counts are pairs too).
#
# Input:  data/raw/placement/read_placements_metazoa_extracted/<gene>/   read placements (Richter et al.)
#         data/raw/placement/epang_multiplicities_metaT_rawcounts/<gene>/ sequence placements
#         data/raw/placement/hk_gene_taxonomy/<gene>.taxonomy.txt
# Output: data/intermediate/supplementary/S4/aggregated/<gene>.rds, a list of tibbles with gene, sample
#           and method: samples (matched samples, and the sequence multiplicity before removing the
#           metazoa), taxonomy (multiplicity and n_units per best_taxopath x conf_depth), lwr and edpl
#           (multiplicity per bin, lower bound), krd (pairs of jplaces), settings
#         Caches (each level only recomputed when missing; not copied into the project):
#           S4/seq_metazoa_extracted/<gene>/, S4/per_name/<method>/<gene>/, S4/krd/
#         Changing CONF_THRESHOLD needs the per_name and aggregated caches deleted.
# Run time: hours for the 7 genes with 16 cores (per_name: ~900 MB, krd: ~300 MB). Needs gappa.
# Original: gene_explorer scripts/analysis/compare_read_vs_sequence_placement.R (there also the resolution
#           and taxonomy figures of scripts/figures/visualize_read_vs_seq_placement.R use these caches)
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/supplementary/S4/01_read_vs_seq_placement_caches.R [gene ...]
#   (default: the 7 genes of the paper)

source('R/load.R')
library(parallel)

N_CORES <- 16
CONF_THRESHOLD <- 0.9   # accumulated LWR needed to call a rank resolved
LWR_BIN <- 0.01         # bin widths of the LWR / EDPL distributions in the cache
EDPL_BIN <- 0.01

seq_dir <- p('placement', 'epang_multiplicities_metaT_rawcounts')
read_dir <- p('placement', 'read_placements_metazoa_extracted')
tax_dir <- p('placement', 'hk_gene_taxonomy')
s4_dir <- p('intermediate', 'supplementary', 'S4')
seq_extracted_dir <- file.path(s4_dir, 'seq_metazoa_extracted')
aggregated_dir <- file.path(s4_dir, 'aggregated')

# MATOU fraction codes -> size fraction of the read sample names
fraction_codes <- c(GGMM = '0o8-5', MMQQ = '5-20', QQSS = '20-180', SSUU = '180-2000', GGZZ = '0o8-inf')

# 1. Samples --------------------------------------------------------------------------------------------
match_samples <- function(gene) {
    seq_samples <- tibble(seq_path = list.files(file.path(seq_dir, gene), '\\.jplace\\.gz$', full.names = TRUE)) %>%
        mutate(seq_sample = basename(seq_path) %>% str_remove('filter\\.jplace\\.gz$')) %>%
        extract(seq_sample, c('station', 'depth', 'replicate', 'fraction_code', 'year'),
                '^(\\d+)(SUR|DCM)(\\d)([A-Z]{4})(\\d{2})$', remove = FALSE, convert = TRUE) %>%
        filter(fraction_code %in% names(fraction_codes)) %>%
        mutate(fraction = fraction_codes[fraction_code],
               sample = sprintf('%03d_%s_%s', station, fraction, depth)) %>%
        arrange(sample, replicate, desc(year)) %>%
        distinct(sample, .keep_all = TRUE)

    read_samples <- tibble(read_path = list.files(file.path(read_dir, gene), '\\.jplace\\.gz$', full.names = TRUE)) %>%
        mutate(sample = basename(read_path) %>% str_remove('_basal_branches\\.jplace\\.gz$'))

    inner_join(seq_samples, read_samples, by = 'sample') %>%
        transmute(gene, sample, seq_sample, station, depth, fraction, seq_path, read_path)
}

# total multiplicity of a jplace (names without multiplicity count 1)
jplace_multiplicity <- function(path) {
    json <- jsonlite::read_json(path)
    map_dbl(json$placements, function(pq) {
        if (!is.null(pq$nm)) sum(map_dbl(pq$nm, ~as.double(.x[[2]]))) else length(pq$n)
    }) %>% sum()
}

# 2. Metazoa removal of the sequence placements ------------------------------------------------------------
extract_metazoa <- function(gene, seq_sample, seq_path) {
    outdir <- file.path(seq_extracted_dir, gene)
    out <- file.path(outdir, str_glue('{seq_sample}filter_basal_branches.jplace.gz'))
    if (file.exists(out)) return(out)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    tmpdir <- tempfile()
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))
    taxfile <- file.path(tax_dir, str_glue('{gene}.taxonomy.txt'))

    # clade list: the metazoan references (no metazoa: the jplace is used as it is)
    system(str_glue("grep 'Metazoa' {taxfile} | sed 's/;Metazoa.*//g' > {tmpdir}/metazoan_vars.txt || true"))
    if (file.size(file.path(tmpdir, 'metazoan_vars.txt')) == 0) {
        file.copy(seq_path, out)
        return(out)
    }
    status <- system2('gappa', c('edit', 'extract', '--jplace-path', seq_path,
                                 '--clade-list-file', file.path(tmpdir, 'metazoan_vars.txt'),
                                 '--basal-clade-name', 'basal_branches',
                                 '--samples-file-prefix', str_glue('{seq_sample}filter_'),
                                 '--samples-out-dir', tmpdir, '--threads', '1'),
                      stdout = FALSE, stderr = FALSE)
    basal <- file.path(tmpdir, str_glue('{seq_sample}filter_basal_branches.jplace'))
    if (status != 0 || !file.exists(basal)) {
        warning('gappa edit extract failed for ', seq_path)
        return(NA_character_)
    }
    system2('gzip', basal)
    file.copy(str_c(basal, '.gz'), out)
    out
}

# 3. Taxonomic assignment and EDPL per read / unigene ----------------------------------------------------
assign_placements <- function(jplace, gene, method, sample) {
    outdir <- file.path(s4_dir, 'per_name', method, gene)
    out <- file.path(outdir, str_glue('{sample}.tsv.gz'))
    if (file.exists(out)) return(out)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    tmpdir <- tempfile()
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    status_assign <- system2('gappa', c('examine', 'assign', '--jplace-path', jplace,
                                        '--taxon-file', file.path(tax_dir, str_glue('{gene}.taxonomy.txt')),
                                        '--out-dir', tmpdir, '--per-query-results', '--threads', '1'),
                             stdout = FALSE, stderr = FALSE)
    status_edpl <- system2('gappa', c('examine', 'edpl', '--jplace-path', jplace, '--out-dir', tmpdir, '--threads', '1'),
                           stdout = FALSE, stderr = FALSE)
    if (status_assign != 0 || status_edpl != 0) {
        warning('gappa failed for ', jplace)
        return(NA_character_)
    }

    per_query <- read_tsv(file.path(tmpdir, 'per_query.tsv'), col_types = 'cddddc') %>%
        mutate(depth = str_count(taxopath, ';') + 1)
    best <- per_query %>%
        arrange(name, desc(fract), desc(depth)) %>%
        distinct(name, .keep_all = TRUE) %>%
        select(name, best_lwr = fract, best_taxopath = taxopath)
    conf <- per_query %>%
        filter(afract >= CONF_THRESHOLD) %>%
        arrange(name, desc(depth)) %>%
        distinct(name, .keep_all = TRUE) %>%
        select(name, conf_taxopath = taxopath, conf_depth = depth)
    # names of a pquery with several names are joined by ';' in per_query, one row per name in the EDPL list
    edpl <- read_csv(file.path(tmpdir, 'edpl_list.csv'), col_types = 'ccdd') %>%
        select(name = Pquery, multiplicity = Multiplicity, edpl = EDPL)

    best %>%
        left_join(conf, by = 'name') %>%
        separate_rows(name, sep = ';') %>%
        inner_join(edpl, by = 'name') %>%
        write_tsv(out)
    out
}

read_per_name <- function(path, method) {
    df <- read_tsv(path, col_types = cols(name = 'c', best_lwr = 'd', best_taxopath = 'c', conf_taxopath = 'c',
                                          conf_depth = 'i', multiplicity = 'd', edpl = 'd'))
    if (method == 'sequences') {
        # one frame per unigene: the one with the highest best-hit LWR
        df <- df %>%
            mutate(unit = str_remove(name, '_frame=.*')) %>%
            arrange(unit, desc(best_lwr)) %>%
            distinct(unit, .keep_all = TRUE) %>%
            select(-unit)
    }
    df
}

# 4. KR distance between all the jplaces of a gene, each pair once ------------------------------------------
run_krd <- function(gene, jplaces) {
    out <- file.path(s4_dir, 'krd', str_glue('{gene}_krd-list.csv'))
    if (file.exists(out)) return(out)
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    status <- system2('gappa', c('analyze', 'krd', '--jplace-path', jplaces, '--matrix-format', 'list', '--normalize',
                                 '--out-dir', dirname(out), '--file-prefix', str_glue('{gene}_'), '--threads', N_CORES),
                      stdout = FALSE, stderr = FALSE)
    if (status != 0) {
        warning('gappa analyze krd failed for ', gene)
        return(NA_character_)
    }
    file.rename(file.path(dirname(out), str_glue('{gene}_krd_matrix.csv')), out)
    out
}

# Per-gene aggregation -------------------------------------------------------------------------------------
aggregate_gene <- function(gene) {
    out <- file.path(aggregated_dir, str_glue('{gene}.rds'))
    if (file.exists(out)) return(out)
    message('== ', gene)

    samples <- match_samples(gene)
    message(nrow(samples), ' matched samples')
    samples$seq_extracted_path <- unlist(mcmapply(extract_metazoa, samples$gene, samples$seq_sample, samples$seq_path,
                                                  mc.cores = N_CORES))
    samples$seq_total_before_extraction <- unlist(mclapply(samples$seq_path, jplace_multiplicity, mc.cores = N_CORES))

    jobs <- bind_rows(samples %>% transmute(gene, sample, method = 'reads', jplace = read_path),
                      samples %>% transmute(gene, sample, method = 'sequences', jplace = seq_extracted_path)) %>%
        filter(!is.na(jplace))
    message('gappa examine assign / edpl on ', nrow(jobs), ' jplaces')
    jobs$per_name_path <- unlist(mcmapply(assign_placements, jobs$jplace, jobs$gene, jobs$method, jobs$sample,
                                          mc.cores = N_CORES))
    jobs <- filter(jobs, !is.na(per_name_path))

    per_name <- jobs %>%
        mutate(data = mcmapply(read_per_name, per_name_path, method, SIMPLIFY = FALSE, mc.cores = N_CORES)) %>%
        select(gene, sample, method, data) %>%
        unnest(data)

    taxonomy <- per_name %>%
        group_by(gene, sample, method, best_taxopath, conf_depth) %>%
        summarize(multiplicity = sum(multiplicity), n_units = n(), .groups = 'drop')
    lwr <- per_name %>%
        mutate(best_lwr = floor(round(best_lwr / LWR_BIN, 6)) * LWR_BIN) %>%
        group_by(gene, sample, method, best_lwr) %>%
        summarize(multiplicity = sum(multiplicity), .groups = 'drop')
    edpl <- per_name %>%
        filter(!is.na(edpl)) %>%
        mutate(edpl = floor(round(edpl / EDPL_BIN, 6)) * EDPL_BIN) %>%
        group_by(gene, sample, method, edpl) %>%
        summarize(multiplicity = sum(multiplicity), .groups = 'drop')

    message('gappa analyze krd')
    krd_path <- run_krd(gene, jobs$jplace)
    krd <- if (is.na(krd_path)) tibble() else {
        jplace_names <- jobs %>% transmute(sample, method, jplace_name = basename(jplace) %>% str_remove('\\.jplace\\.gz$'))
        read_tsv(krd_path, col_names = c('first', 'second', 'krd'), col_types = 'ccd') %>%
            inner_join(rename_with(jplace_names, ~str_c(.x, '_1'), -jplace_name), by = c('first' = 'jplace_name')) %>%
            inner_join(rename_with(jplace_names, ~str_c(.x, '_2'), -jplace_name), by = c('second' = 'jplace_name')) %>%
            filter(str_c(sample_1, method_1) < str_c(sample_2, method_2)) %>%
            transmute(gene = !!gene, sample_1, method_1, sample_2, method_2, krd)
    }

    saveRDS(list(samples = select(samples, -seq_path, -read_path, -seq_extracted_path),
                 taxonomy = taxonomy, lwr = lwr, edpl = edpl, krd = krd,
                 settings = list(conf_threshold = CONF_THRESHOLD, lwr_bin = LWR_BIN, edpl_bin = EDPL_BIN)),
            out)
    out
}

# Run ------------------------------------------------------------------------------------------------------
# genes: the arguments, otherwise the 7 genes of the paper figure (those in both datasets when it was made;
# the sequence placements now hold more genes in common with the reads, e.g. ALG11, ap1m1, ARP2, ...)
PAPER_GENES <- c('ABCE', 'ABT1', 'AGX', 'FAM', 'psma6', 'rps5', 'rps9')
dir.create(aggregated_dir, recursive = TRUE, showWarnings = FALSE)
args_genes <- commandArgs(trailingOnly = TRUE)
genes <- if (length(args_genes) > 0) args_genes else PAPER_GENES
available <- intersect(list.files(seq_dir), list.files(read_dir)) %>%
    keep(~file.exists(file.path(tax_dir, str_glue('{.x}.taxonomy.txt'))))
if (length(setdiff(genes, available)) > 0) stop('Not in both datasets: ', str_flatten_comma(setdiff(genes, available)))
message('Genes: ', str_flatten_comma(genes))
walk(genes, aggregate_gene)
