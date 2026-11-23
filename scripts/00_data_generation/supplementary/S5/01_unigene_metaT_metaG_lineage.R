# Supp S5, step 01: for Flavodoxin, SIT and NRT2, one row per unigene x metaT sample with its normalised
# metaT, the metaG of the same unigene in the paired metaG sample and the lineage of its placement.
#
#   - paired metaG sample: same station, depth and size fraction (sample codes equal except the last two
#     digits before "filter"); a metaT sample is usable only when the pair has exactly one metaT and one
#     metaG jplace (as in the metaT/metaG ratios)
#   - has_metag: the unigene has metaG > 0 in the paired sample
#   - lineage (group): taxonomy group (tip_taxonomy(), R/phylofactorization.R) of the reference tips below
#     the best EPA-ng placement edge of the unigene; "mixed" when they belong to more than one group
#
# Input:  data/raw/placement/epang/<gene>.epa_result.jplace.gz (all unigenes, placements and tree)
#         data/raw/placement/epang_multiplicities_metaT/<gene>/   normalised MATOU metaT per sample
#         data/raw/placement/epang_multiplicities_metaG/<gene>/   normalised MATOU metaG per sample
# Output: data/intermediate/supplementary/S5/<gene>.tsv.gz
#         (sample_code, pair, geneid, metaT, metag_sample, metaG, has_metag, group)
# Run time: ~10 min (8 cores).
# Original: gene_explorer scripts/figures/metag_coverage_by_lineage.R (build_gene(), its cache in
#           data/stats/metag_coverage_by_lineage/)
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/supplementary/S5/01_unigene_metaT_metaG_lineage.R [gene ...] [--out <dir>]

source('R/load.R')
library(ape)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- p('intermediate', 'supplementary', 'S5')
if ('--out' %in% args) {
    out_dir <- args[which(args == '--out') + 1]
    args <- args[-c(which(args == '--out'), which(args == '--out') + 1)]
}
genes <- if (length(args) > 0) args else c('Flavodoxin', 'SIT', 'NRT2')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Placement lineage ---------------------------------------------------------------------------------

# taxonomy group of the reference tips below every edge ("mixed" when more than one); the tree comes
# from read_jplace_tree() (R/io_placement.R), with the jplace edge numbers in the node labels
edge_groups <- function(tree) {
    n_tip <- Ntip(tree)
    labels <- c(tree$tip.label, tree$node.label)
    edge <- as.integer(str_extract(labels, '(?<=_E)\\d+$'))
    tree$tip.label <- str_remove(tree$tip.label, '_E\\d+$')
    tip_group <- as.character(tip_taxonomy(tree)$group)
    groups <- vector('list', n_tip + Nnode(tree))
    groups[seq_len(n_tip)] <- as.list(tip_group)
    postorder <- reorder(tree, 'postorder')$edge
    for (i in seq_len(nrow(postorder))) {
        groups[[postorder[i, 1]]] <- union(groups[[postorder[i, 1]]], groups[[postorder[i, 2]]])
    }
    tibble(edge = edge, group = map_chr(groups, ~if (length(.x) == 1) .x else 'mixed')) %>%
        filter(!is.na(edge))
}

# best placement edge (highest LWR) of every unigene of the gene, and its lineage
best_placements <- function(gene) {
    full <- read_jplace(gzfile(p('placement', 'epang', str_c(gene, '.epa_result.jplace.gz'))))
    fields <- unlist(full$fields)
    edge_i <- which(fields == 'edge_num')
    lwr_i <- which(fields == 'like_weight_ratio')
    best <- vapply(full$placements, function(pq) {
        lwr <- vapply(pq$p, function(p) p[[lwr_i]], 0)
        pq$p[[which.max(lwr)]][[edge_i]]
    }, numeric(1))
    tibble(geneid = vapply(full$placements, function(pq) pq$n[[1]], ''), edge = as.integer(best)) %>%
        left_join(edge_groups(read_jplace_tree(full$tree)), by = 'edge')
}

# Multiplicities --------------------------------------------------------------------------------------

# value of every unigene in one sample jplace
read_multiplicities <- function(path) {
    pq <- read_jplace(gzfile(path))$placements
    tibble(sample_code = str_remove(basename(path), '\\.jplace\\.gz$'),
           geneid = vapply(pq, function(p) if (is.null(p$nm)) p$n[[1]] else p$nm[[1]][[1]], ''),
           value = vapply(pq, function(p) if (is.null(p$nm)) NA_real_ else as.numeric(p$nm[[1]][[2]]), 0))
}

read_multiplicity_dir <- function(dir) {
    files <- list.files(dir, pattern = 'jplace.gz$', full.names = TRUE)
    bind_rows(parallel::mclapply(files, read_multiplicities, mc.cores = 8)) %>%
        # pair: station, depth and fraction, shared by the metaT (..14filter) and metaG (..1Xfilter) jplaces
        mutate(pair = str_remove(sample_code, '[0-9]{2}filter$'))
}

# One gene ---------------------------------------------------------------------------------------------

for (gene in genes) {
    message(gene)
    metaT <- read_multiplicity_dir(p('placement', 'epang_multiplicities_metaT', gene)) %>%
        mutate(value = replace_na(value, 1))   # gappa counts a pquery without multiplicity as 1
    metaG <- read_multiplicity_dir(p('placement', 'epang_multiplicities_metaG', gene))

    good_pairs <- inner_join(count(distinct(metaT, sample_code, pair), pair, name = 'n_t'),
                             count(distinct(metaG, sample_code, pair), pair, name = 'n_g'), by = 'pair') %>%
        filter(n_t == 1, n_g == 1) %>%
        pull(pair)

    metaT %>%
        transmute(sample_code, pair, geneid, metaT = value, metag_sample = pair %in% good_pairs) %>%
        left_join(metaG %>% filter(pair %in% good_pairs) %>% select(pair, geneid, metaG = value),
                  by = c('pair', 'geneid')) %>%
        mutate(has_metag = metag_sample & !is.na(metaG) & metaG > 0) %>%
        left_join(select(best_placements(gene), geneid, group), by = 'geneid') %>%
        write_tsv(file.path(out_dir, str_c(gene, '.tsv.gz')))
}
