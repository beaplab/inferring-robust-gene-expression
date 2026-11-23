# Fig 3B-D: crossmapping between F. parvula and every other transcriptome / genome.
#
# The metatranscriptomic reads of each sample were mapped separately against each transcriptome
# (quantification pipeline, one BAM per sample and transcriptome). A read of a F. parvula gene
# "crossmaps" when the same read (same name) also maps to the other transcriptome. Per F. parvula gene,
# summed over the samples shared by both:
#   ratio = crossmapped reads / not crossmapped reads
#   category: 0 'no crossmapping', <= 0.1 'small amount of crossmapping (<10%)',
#             >= 0.5 'large crossmapping (>50%)', otherwise 'intermediate situation (10-50%)',
#             no read without crossmapping 'only reads crossmapped!'
#   (the labels say %, but the ratio is crossmapped / not crossmapped, not a share of the reads)
#
# Input:  BAMs of the quantification pipeline (delving_into_florenciella data/quantification/alignment/
#         <transcriptome>/<sample>.filtered.bam*, not copied into this project)
#         data/raw/genomic/functional_annotation/EP00618_Florenciella_parvula_nucnames.tsv (eggNOG)
# Output: data/intermediate/fig3/crossmapping/<comparison>.rds   per sample, gene and fate, with the
#                                                                genes of the other species that took
#                                                                the reads (used to pick the examples of
#                                                                05_coverage_examples.R)
#         data/intermediate/fig3/crossmapping_florenciella_against_all_transcriptomes.tsv.gz
#         (comparison, gene_name, reads_nocrossmap, reads_crossmap, ratio, category_cross + eggNOG)
# Run time: hours (reads every BAM). Needs Rsamtools.
# Original: delving_into_florenciella src/analysis/calculate_crossmapping-transcriptomes.R (one run
#           per comparison, crossmapping_calculation.sh) + src/analysis/helper_crossmapping_Florenciella.R
#
# Usage (from the project root): Rscript scripts/00_data_generation/fig3/03_crossmapping.R

source('R/load.R')
suppressPackageStartupMessages(library(Rsamtools))

bam_dir <- file.path(external$delving_florenciella, 'data/quantification/alignment')
reference <- 'EP00618_Florenciella_parvula'

# comparison name -> BAM folder of the other transcriptome / genome
comparisons <- c(RCC1587 = 'EP00620_Florenciella_sp_RCC1587',
                 RCC1693 = 'EP00621_Florenciella_sp_RCC1693',
                 RCC1007 = 'EP00619_Florenciella_sp_RCC1007',
                 TARAARC108 = 'TARA_ARC_108_MAG_00262_CDS',
                 TARAMED95 = 'TARA_MED_95_MAG_00475_CDS',
                 TARASOC28 = 'TARA_SOC_28_MAG_00069_CDS',
                 Pelagomonas = 'Pelagomonas_calceolata',
                 Bathycoccus = 'EP00893_Bathycoccus_prasinos')

out_dir <- p('intermediate', 'fig3', 'crossmapping')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# read name -> gene it maps to
bam_reads <- function(path) {
    bam <- scanBam(path, param = ScanBamParam(what = c('qname', 'rname')))[[1]]
    tibble(qname = bam$qname, refname = as.character(bam$rname))
}

# Per F. parvula gene of one sample: reads that also map to the other transcriptome (crossmapping =
# TRUE) or not, and the genes of the other transcriptome that took them
sample_crossmapping <- function(bam_ref, bam_query) {
    bam_reads(bam_ref) %>%
        left_join(bam_reads(bam_query), by = 'qname', suffix = c('', '.other')) %>%
        mutate(crossmapping = !is.na(refname.other)) %>%
        group_by(gene_name = refname, crossmapping) %>%
        summarize(reads = n(), genes_cross_set = list(unique(refname.other)), .groups = 'drop')
}

crossmapping_per_sample <- function(other) {
    bams <- function(dir) list.files(file.path(bam_dir, dir), pattern = '\\.bam', recursive = TRUE)
    # .bam / .bam.sorted files, without their indexes
    shared <- intersect(bams(reference), bams(other)) %>% str_subset('\\.(bai|csi)$', negate = TRUE)
    message(other, ': ', length(shared), ' shared samples')
    tibble(sample = shared) %>%
        mutate(stats = map(sample, ~sample_crossmapping(file.path(bam_dir, reference, .x),
                                                       file.path(bam_dir, other, .x)))) %>%
        unnest(stats)
}

classify_crossmapping <- function(per_sample) {
    per_sample %>%
        group_by(gene_name, crossmapping) %>%
        summarize(reads = sum(reads), .groups = 'drop') %>%
        pivot_wider(names_from = crossmapping, values_from = reads) %>%
        rename(reads_nocrossmap = `FALSE`, reads_crossmap = `TRUE`) %>%
        mutate(ratio = ifelse(is.na(reads_crossmap), 0, reads_crossmap / reads_nocrossmap),
               category_cross = case_when(
                   ratio == 0 ~ 'no crossmapping',
                   ratio <= 0.1 ~ 'small amount of crossmapping (<10%)',
                   ratio >= 0.5 ~ 'large crossmapping (>50%)',
                   is.na(ratio) ~ 'only reads crossmapped!',
                   TRUE ~ 'intermediate situation (10-50%)'))
}

eggnog <- read_tsv(p('annotation', 'EP00618_Florenciella_parvula_nucnames.tsv'), show_col_types = FALSE)

crossmapping <- imap(comparisons, function(other, name) {
    per_sample <- crossmapping_per_sample(other)
    write_rds(per_sample, file.path(out_dir, str_c(name, '.rds')))
    classify_crossmapping(per_sample)
}) %>%
    bind_rows(.id = 'comparison') %>%
    left_join(eggnog, by = 'gene_name')

write_tsv(crossmapping, p('intermediate', 'fig3', 'crossmapping_florenciella_against_all_transcriptomes.tsv.gz'))
crossmapping %>% count(comparison, category_cross) %>% print(n = Inf)
