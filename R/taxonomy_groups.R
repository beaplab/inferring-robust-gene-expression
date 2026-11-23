# Shared taxonomic groups for the barplots, so that the gappa assign
# (housekeeping genes), MATOU (NCBI taxonomy) and 18S-V4/V9 (EukRibo taxonomy)
# compositions use the same groups and colours.
# The groups are the main marine protist supergroups present in all the
# datasets. Haptista is used instead of Haptophyta because NCBI resolves many
# MATOU genes only to 'Eukaryota;Haptista' (it also takes the few
# centrohelids). Anaerobic groups (Fornicata, Metamonada, Breviatea) are not
# groups of their own and go to 'Other': they are absent from 18S and MATOU and
# only appear in some housekeeping genes, as misplacements onto their
# long-branch references.
#
# Main functions:
#   assign_taxon_group(lineage)  # lineage with ';' or '|' separators
#   median_composition(df, ...)  # housekeeping-gene median, see below
# Original: gene_explorer scripts/figures/taxonomy_groups.R

taxon_groups <- c('Alveolata', 'Stramenopiles', 'Rhizaria', 'Opisthokonta',
                  'Haptista', 'Chloroplastida', 'Cryptista', 'Telonema', 'Eukaryota')

taxa_colors <- set_names(c('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300',
                           '#4a3aa7', '#0d366b', '#8c564b', '#b5b5b5'),
                         c(taxon_groups, 'Other'))

# lineage names (any taxonomy) for each group, checked in this order so the
# more specific ones win
group_patterns <- c(
    'Telonema'       = 'Telonemia|Telonemida|Telonemea|Telonema',
    'Alveolata'      = 'Alveolata',
    'Stramenopiles'  = 'Stramenopiles',
    'Rhizaria'       = 'Rhizaria',
    'Haptista'       = 'Haptista|Haptophyta',
    'Cryptista'      = 'Cryptista|Cryptophyta|Cryptophyceae',
    'Chloroplastida' = 'Chloroplastida|Viridiplantae',
    'Opisthokonta'   = 'Opisthokonta'
)

# ranks of a lineage, without empty ones, the NCBI 'root;cellular organisms'
# and the '*' placeholders of the 18S taxonomies
lineage_ranks <- function(lineage) {
    str_split(lineage, '[;|]') %>%
        map(~.x[!.x %in% c('', 'root', 'cellular organisms', '*')])
}

assign_taxon_group <- function(lineage) {
    ranks <- lineage_ranks(lineage)

    group <- rep(NA_character_, length(lineage))
    for (g in names(group_patterns)) {
        hit <- is.na(group) &
            map_lgl(ranks, ~any(str_detect(.x, str_c('^(', group_patterns[[g]], ')$'))))
        group[hit] <- g
    }

    is_euk <- map_lgl(ranks, ~'Eukaryota' %in% .x)
    only_euk <- map_lgl(ranks, ~length(.x) > 0 && last(.x) == 'Eukaryota')
    case_when(!is.na(group) ~ group,
              only_euk ~ 'Eukaryota',   # unresolved below Eukaryota
              is_euk ~ 'Other',
              TRUE ~ NA_character_)     # not eukaryotic / unassigned
}

# Housekeeping-gene median ----------------------------------------------------

# depth of the lineages on which the median is taken (e.g.
# Eukaryota;Diaphoretickes;Sar;Alveolata;Dinoflagellata)
median_tax_level <- 5

fine_lineage <- function(lineage, level = median_tax_level) {
    map_chr(lineage_ranks(lineage), ~str_c(head(.x, level), collapse = ';'))
}

# Composition across genes: the median over genes of the relative expression of
# each fine lineage (0 in the genes where it is absent), summed into the plot
# groups and rescaled to 1. Taking the median per fine lineage, not per group,
# makes it independent of how broad the groups are: lineages present in only
# some genes (misplacements) get a median of 0 even when pooled into 'Other'.
#   df:  gene, lineage (fine), taxon, rel_value and the columns in ...
#   ...: columns identifying a composition (e.g. samplename)
median_composition <- function(df, ...) {
    id_cols <- names(select(df, ...))
    lineage_taxa.df <- distinct(df, lineage, taxon)

    df %>%
        group_by(across(all_of(c('gene', id_cols, 'lineage')))) %>%
        summarize(rel_value = sum(rel_value), .groups = 'drop') %>%
        complete(nesting(!!!syms(c('gene', id_cols))), lineage, fill = list(rel_value = 0)) %>%
        group_by(across(all_of(c(id_cols, 'lineage')))) %>%
        summarize(rel_value = median(rel_value), .groups = 'drop') %>%
        inner_join(lineage_taxa.df, by = 'lineage') %>%
        group_by(across(all_of(c(id_cols, 'taxon')))) %>%
        summarize(rel_value = sum(rel_value), .groups = 'drop_last') %>%
        mutate(rel_value = rel_value / sum(rel_value)) %>%
        ungroup()
}

# the gappa placements had the metazoan reads removed, do the same elsewhere
is_metazoa <- function(lineage) {
    str_detect(lineage, '(^|[;|])Metazoa([;|]|$)')
}

# Housekeeping gene selection --------------------------------------------------

# genes in the catalogue that are functional markers, not housekeeping genes
non_housekeeping_genes <- c('SIT', 'Flavodoxin', 'NRT2')

is_ribosomal <- function(gene) str_detect(gene, regex('^(rp[ls]\\d|rplp|rpsa)', ignore_case = TRUE))
is_proteasome <- function(gene) str_detect(gene, regex('^(psm[abcd]|psd\\d)', ignore_case = TRUE))

# n random housekeeping genes, always including at least one ribosomal and
# one proteasome gene
select_random_genes <- function(genes, n = 10, seed = 42) {
    genes <- setdiff(genes, non_housekeeping_genes)
    set.seed(seed)
    ribosomal <- sample(genes[is_ribosomal(genes)], 1)
    proteasome <- sample(genes[is_proteasome(genes)], 1)
    others <- sample(setdiff(genes, c(ribosomal, proteasome)), n - 2)
    c(ribosomal, proteasome, others)
}
