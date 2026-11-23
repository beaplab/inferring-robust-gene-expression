# metaT/metaG ratios per unigene (MATOU), used by the data-generation scripts of Fig 2.
# The per-sample jplaces of a gene carry, as multiplicities, the metaT read counts
# (epang_multiplicities_metaT_rawcounts) or the metaG values (epang_multiplicities_metaG) of each unigene.
# MATOU metaT and metaG sample codes differ only in the two digits before "filter"
# (e.g. 4SUR1GGMM14filter / 4SUR1GGMM11filter), so samples are paired on the code without them.

# Sample code without the two digits before "filter": the key shared by metaT and metaG
matou_pair_key <- function(sample_code) str_remove(sample_code, '[0-9][0-9]filter.*$')

# jplaces of a folder, named by their pairing key. Keys with more than one jplace (155DCM1KKQQ13 and 14
# in metaT, 7DCM1MMQQ11 and 12 in metaG) cannot be paired unambiguously and are left out.
unique_sample_jplaces <- function(dir) {
    files <- list.files(dir, pattern = '\\.jplace\\.gz$', full.names = TRUE)
    keys <- matou_pair_key(basename(files))
    duplicated_keys <- unique(keys[duplicated(keys)])
    if (length(duplicated_keys) > 0) message('Skipping duplicated samples: ', str_c(duplicated_keys, collapse = ', '))
    keep <- !keys %in% duplicated_keys
    set_names(files[keep], keys[keep])
}

# metaT and metaG jplaces of the samples present in both: tibble(key, metat_jplace, metag_jplace)
paired_jplaces <- function(metat_dir, metag_dir) {
    metat <- unique_sample_jplaces(metat_dir)
    metag <- unique_sample_jplaces(metag_dir)
    keys <- intersect(names(metat), names(metag))
    tibble(key = keys, metat_jplace = unname(metat[keys]), metag_jplace = unname(metag[keys]))
}

# metaT read count and metaG value of each unigene of one sample:
# tibble(sample = metaT sample code, geneid, count_metat, count_metag); count_metag is NA for unigenes
# without metaG
metat_metag_counts <- function(metat_jplace, metag_jplace) {
    metag <- gappa_multiplicities(metag_jplace) %>% select(geneid = name, count_metag = count)
    counts <- gappa_multiplicities(metat_jplace) %>%
        select(sample = jplace_sample, geneid = name, count_metat = count) %>%
        left_join(metag, by = 'geneid')
    # one row per unigene: a duplicated one would get two multiplicities in the ratio jplace
    stopifnot(!anyDuplicated(counts$geneid))
    counts
}

# Ratio per unigene, as metaT reads per metaG RPM:
#   ratio = count_metat / count_metag / 1e6
# The MATOU metaG value is reads / total reads of the library (a proportion), so the division by 1e6
# puts the "per million" back. Not rounded: gappa takes decimal multiplicities.
# - Unigenes without metaG (NA) or with a metaG value of 0 (Inf) get a ratio of 0, not dropped:
#   gappa edit multiplicity --multiplicity-file only overwrites the queries listed, so a dropped unigene
#   would keep its metaT read count in the ratio jplace.
# - Samples whose ratios sum to <= 1 are dropped: gappa analyze placement-factorization refuses them
#   ("Cannot calculate balance data on trees with normalized masses").
metat_metag_ratios <- function(counts) {
    ratios <- counts %>%
        mutate(ratio = count_metat / count_metag / 1e6,
               ratio = if_else(is.finite(ratio), ratio, 0)) %>%
        group_by(sample) %>%
        mutate(total = sum(ratio)) %>%
        ungroup()

    too_small <- distinct(filter(ratios, total <= 1), sample, total)
    if (nrow(too_small) > 0) message('Dropping samples with a summed ratio <= 1: ',
                                     str_c(too_small$sample, ' (', signif(too_small$total, 2), ')', collapse = ', '))
    ratios %>%
        filter(total > 1) %>%
        select(sample, geneid, ratio)
}
