# Fig 6: taxonomic composition of every Tara sample with housekeeping (HK) gene assignments, as seen
# by each omic dataset, in the shared taxonomic groups of R/taxonomy_groups.R (Metazoa removed):
#   18S-V4 and 18S-V9 amplicons (Tara swarm tables)
#   MATOU metaG and metaT, whole gene catalogue (taxonomy of every unigene)
#   each HK gene (gappa assign of its metaT placements, 01-02) and the median of the HK genes
#
# Steps
#   1. Raw caches (slow, built only if missing; delete a file to recompute it), in data/intermediate/fig6/cache/:
#        tara_<marker>_sample_lineage.parquet      18S reads per PANGAEA sample and lineage
#        matou_taxonomy_lineage_ids.parquet         unigene -> lineage id
#        matou_lineages.parquet                     lineage id -> lineage
#        matou_<metaG|metaT>_sample_lineage.parquet value per MATOU sample and lineage (the 25 GB metaT
#                                                   and 9 GB metaG occurrence tables, streamed with arrow)
#        hk_gene_sample_lineage.parquet             HK metaT per gene, sample and lineage (depth
#                                                   median_tax_level)
#   2. HK genes: relative expression per gene and sample; their median composition per sample
#   3. 18S: MATOU sample -> PANGAEA samples (sample_material); replicates pooled
#   4. MATOU metaG / metaT: metaT codes match exactly; metaG codes differ in the replicate and year
#      digits, so they are matched on station, depth and size fraction, pooling the samples that share it
#   5. One table with every composition
#
# Input:  Tara 18S-V4 / V9 swarm tables and taxonomies (external, see below)
#         MATOU v1.5 taxonomy and metaG / metaT occurrences (external$matou)
#         HK gappa assign summaries: data/intermediate/hk_taxassign_combined (02; if missing, the copy used
#           for the paper, data/raw/placement/hk_taxassign_combined)
#         data/raw/env/tara_env.tsv.gz (sample correspondences: MATOU code, sample_material, PANGAEA id)
# Output: data/intermediate/fig6/taxon_compositions.parquet (samplename, source, taxon, rel_value)
# Run time: hours the first time (MATOU occurrence tables), minutes once the raw caches exist.
# Original: gene_explorer scripts/analysis/build_taxonomy_caches_18S_matou.R (without the random unigene
#           pool and the per-gene relative-expression cache, not used in the paper)
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig6/04_taxonomy_compositions.R

source('R/load.R')
library(arrow, warn.conflicts = FALSE)
set_cpu_count(16)

# Paths ------------------------------------------------------------------------------------------
swarm_files <- list(
    V4 = c(table = '/scratch/data4/aauladell/tara_v4_swarm/TARA-Oceans_18S-V4_Swarm-Mumu_table.tsv.gz',
           taxo = '/scratch/data4/aauladell/tara_v4_swarm/TARA-Oceans_18S-V4_Swarm-Mumu_taxo.tsv.gz'),
    V9 = c(table = '/home/aauladell/datasets/common_databases/tara_v9_swarm/TARA-Oceans_18S-V9_Swarm-Mumu_table.tsv.gz',
           taxo = '/home/aauladell/datasets/common_databases/tara_v9_swarm/TARA-Oceans_18S-V9_Swarm_taxo.tsv.gz'))
matou_taxonomy_file <- file.path(external$matou, 'MATOU-v1.5.taxonomy.tsv.gz')
matou_occurrences_file <- function(dataset) file.path(external$matou, str_c('MATOU-v1.5.', dataset, '.occurrences.gz'))
hk_summary_dir <- p('intermediate', 'hk_taxassign_combined')
if (!dir.exists(hk_summary_dir)) hk_summary_dir <- p('placement', 'hk_taxassign_combined')

cache_dir <- p('intermediate', 'fig6', 'cache')
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
tara_caches <- c(V4 = file.path(cache_dir, 'tara_v4_sample_lineage.parquet'),
                 V9 = file.path(cache_dir, 'tara_v9_sample_lineage.parquet'))
matou_tax_cache <- file.path(cache_dir, 'matou_taxonomy_lineage_ids.parquet')
matou_lineage_cache <- file.path(cache_dir, 'matou_lineages.parquet')
matou_caches <- c(metaG = file.path(cache_dir, 'matou_metaG_sample_lineage.parquet'),
                  metaT = file.path(cache_dir, 'matou_metaT_sample_lineage.parquet'))
hk_cache <- file.path(cache_dir, 'hk_gene_sample_lineage.parquet')

# 1. Raw caches ----------------------------------------------------------------------------------

# Tara 18S: reads per PANGAEA sample and lineage
for (marker in names(tara_caches)) {
    if (file.exists(tara_caches[[marker]])) next
    message('Summarising the ', marker, ' swarm table per sample and lineage')
    taxo <- read_tsv_arrow(swarm_files[[marker]][['taxo']], col_select = c(amplicon, closest_hits_lca))
    counts <- read_tsv_arrow(swarm_files[[marker]][['table']], as_data_frame = FALSE)
    sample_cols <- setdiff(names(counts), c('amplicon', 'maxrelabund', 'total', 'spread', 'sequence'))

    lineage <- tibble(amplicon = as.vector(counts$amplicon)) %>%
        left_join(taxo, by = 'amplicon') %>%
        pull(closest_hits_lca) %>%
        replace_na('unassigned')
    count_mat <- counts %>% select(all_of(sample_cols)) %>% as.data.frame() %>% as.matrix()

    # amplicon x sample -> lineage x sample
    rowsum(count_mat, lineage) %>%
        as_tibble(rownames = 'lineage') %>%
        pivot_longer(-lineage, names_to = 'sample_id_pangaea', values_to = 'reads') %>%
        filter(reads > 0) %>%
        write_parquet(tara_caches[[marker]])
    rm(counts, count_mat)
}

# MATOU taxonomy, lineage strings replaced by integer ids
if (!file.exists(matou_tax_cache) || !file.exists(matou_lineage_cache)) {
    message('Building the MATOU taxonomy cache')
    taxonomy <- read_tsv_arrow(matou_taxonomy_file, col_select = c(geneID, taxLineage), as_data_frame = FALSE)
    lineages <- taxonomy %>%
        distinct(taxLineage) %>%
        collect() %>%
        transmute(lineage_id = row_number(), lineage = taxLineage)
    write_parquet(lineages, matou_lineage_cache)
    taxonomy %>%
        left_join(arrow_table(lineages), by = c('taxLineage' = 'lineage')) %>%
        transmute(geneid = cast(geneID, int64()), lineage_id = cast(lineage_id, int32())) %>%
        write_parquet(matou_tax_cache)
    rm(taxonomy)
}

# MATOU metaG / metaT: value per sample and lineage
for (dataset in names(matou_caches)) {
    if (file.exists(matou_caches[[dataset]])) next
    message('Summarising the ', dataset, ' occurrences per sample and lineage (slow)')
    open_dataset(matou_occurrences_file(dataset), format = 'tsv', skip = 1,
                 schema = schema(geneid = int64(), samplename = utf8(), value = float64())) %>%
        inner_join(open_dataset(matou_tax_cache), by = 'geneid') %>%
        group_by(samplename, lineage_id) %>%
        summarise(value = sum(value), n_genes = n()) %>%
        collect() %>%
        write_parquet(matou_caches[[dataset]])
}

# HK genes: metaT per gene, sample and lineage (cut at depth median_tax_level)
if (!file.exists(hk_cache)) {
    message('Building the HK gene x sample x lineage cache')
    gene_lineages <- open_dataset(list.files(hk_summary_dir, pattern = '_summary.tsv$', full.names = TRUE),
                                  format = 'tsv',
                                  col_types = schema(gene = utf8(), sample = utf8(), taxopath = utf8(),
                                                     n_queries = float64(), expression = float64(),
                                                     expression_lwr = float64())) %>%
        group_by(gene, samplename = sample, taxopath) %>%
        summarise(expression = sum(expression)) %>%
        collect()
    taxopaths <- tibble(taxopath = unique(gene_lineages$taxopath)) %>%
        filter(!is_metazoa(taxopath)) %>%
        mutate(lineage = fine_lineage(taxopath))
    gene_lineages %>%
        inner_join(taxopaths, by = 'taxopath') %>%
        group_by(gene, samplename, lineage) %>%
        summarize(expression = sum(expression), .groups = 'drop') %>%
        write_parquet(hk_cache)
}

# 2. HK genes ------------------------------------------------------------------------------------
hk_lineages <- read_parquet(hk_cache)
hk_taxa <- tibble(lineage = unique(hk_lineages$lineage)) %>%
    mutate(taxon = assign_taxon_group(lineage)) %>%
    filter(!is.na(taxon))

# relative expression per gene and sample (SIT, NRT2 and Flavodoxin are not HK genes)
gene_rel <- hk_lineages %>%
    inner_join(hk_taxa, by = 'lineage') %>%
    filter(!gene %in% non_housekeeping_genes) %>%
    group_by(gene, samplename) %>%
    mutate(rel_value = expression / sum(expression)) %>%
    ungroup() %>%
    select(gene, samplename, lineage, taxon, expression, rel_value)

candidate_samples <- unique(gene_rel$samplename)

genes <- gene_rel %>%
    group_by(samplename, taxon, source = gene) %>%
    summarize(rel_value = sum(rel_value), .groups = 'drop')
hk_median <- gene_rel %>%
    median_composition(samplename) %>%
    mutate(source = str_c('Housekeeping genes (median of ', n_distinct(gene_rel$gene), ')'))

# 3. 18S compositions ------------------------------------------------------------------------------
# MATOU sample code -> sample_material (TARA_<station>_<depth>_<size fraction>) -> PANGAEA samples
fraction_bounds <- c(GG = '0.8', KK = '3', MM = '5', QQ = '20', SS = '180', UU = '2000', ZZ = '>')
depth_codes <- c(SUR = 'SRF', DCM = 'DCM', MES = 'MES', MXL = 'MIX')
matou_code_to_material <- function(code) {
    parts <- str_match(code, '^(\\d+)([A-Z]{3})\\d+([A-Z]{2})([A-Z]{2})\\d{2}$')
    upper <- fraction_bounds[parts[, 5]]  # '>' gives e.g. '0.8->', as in the context table
    str_c('TARA_', str_pad(parts[, 2], 3, pad = '0'), '_', depth_codes[parts[, 3]], '_',
          fraction_bounds[parts[, 4]], '-', upper)
}

# the MATOU metadata first (it has the special cases, e.g. OTH samples), the parsed code otherwise
tara_env <- read_tara_env()
matou_materials <- tara_env %>%
    filter(!is.na(matou_sample_code)) %>%
    transmute(samplename = str_remove(matou_sample_code, 'filter$'), sample_material) %>%
    distinct()
# only the samples of the 18S context table (other PANGAEA samples share some of their sample_material)
context_samples <- tara_env %>%
    filter(in_context_general) %>%
    distinct(sample_id_pangaea, sample_material)
sample_map <- tibble(samplename = candidate_samples) %>%
    left_join(matou_materials, by = 'samplename') %>%
    mutate(sample_material = coalesce(sample_material, matou_code_to_material(samplename))) %>%
    left_join(context_samples, by = 'sample_material', relationship = 'many-to-many')

composition_18S <- function(marker) {
    swarm <- read_parquet(tara_caches[[marker]])
    lineages <- tibble(lineage = unique(swarm$lineage)) %>%
        filter(!is_metazoa(lineage)) %>%
        mutate(taxon = assign_taxon_group(lineage)) %>%
        filter(!is.na(taxon))
    swarm %>%
        inner_join(lineages, by = 'lineage') %>%
        inner_join(sample_map, by = 'sample_id_pangaea', relationship = 'many-to-many') %>%
        group_by(samplename, taxon) %>%
        summarize(reads = sum(reads), .groups = 'drop_last') %>%
        mutate(rel_value = reads / sum(reads)) %>%
        ungroup() %>%
        transmute(samplename, taxon, rel_value, source = str_c('18S-', marker))
}

# 4. MATOU whole-catalogue compositions -------------------------------------------------------------
# unigenes resolved only to 'Eukaryota;Opisthokonta' are mostly unclassified metazoans (they dominate the
# large size fractions), so they go with Metazoa
matou_taxa <- read_parquet(matou_lineage_cache) %>%
    filter(!is_metazoa(lineage), !str_detect(lineage, ';Opisthokonta;$')) %>%
    mutate(taxon = assign_taxon_group(lineage)) %>%
    filter(!is.na(taxon)) %>%
    select(lineage_id, taxon)

# station + depth + size fraction, without replicate and year (e.g. 194SUR1QQSS13 -> 194SURQQSS)
sample_key <- function(code) str_replace(code, '^(\\d+[A-Z]{3})\\d+([A-Z]{4})\\d{2}$', '\\1\\2')

composition_matou <- function(dataset) {
    sample_taxa <- read_parquet(matou_caches[[dataset]]) %>%
        inner_join(matou_taxa, by = 'lineage_id')
    if (dataset == 'metaT') {
        sample_taxa <- filter(sample_taxa, samplename %in% candidate_samples)
    } else {
        sample_taxa <- sample_taxa %>%
            mutate(key = sample_key(samplename)) %>%
            inner_join(tibble(ref_sample = candidate_samples, key = sample_key(candidate_samples)),
                       by = 'key', relationship = 'many-to-many') %>%
            mutate(samplename = ref_sample)
    }
    sample_taxa %>%
        group_by(samplename, taxon) %>%
        summarize(value = sum(value), .groups = 'drop_last') %>%
        mutate(rel_value = value / sum(value)) %>%
        ungroup() %>%
        transmute(samplename, taxon, rel_value, source = str_c('MATOU ', dataset, ' (all genes)'))
}

# 5. All compositions ------------------------------------------------------------------------------
compositions <- bind_rows(composition_18S('V4'), composition_18S('V9'),
                          composition_matou('metaG'), composition_matou('metaT'),
                          hk_median, genes) %>%
    select(samplename, source, taxon, rel_value)
write_parquet(compositions, p('intermediate', 'fig6', 'taxon_compositions.parquet'))

compositions %>%
    distinct(source, samplename) %>%
    filter(!source %in% gene_rel$gene) %>%
    count(source) %>%
    print()
