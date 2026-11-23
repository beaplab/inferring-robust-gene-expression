# Fig 3C: per-base read coverage of the three example genes of Florenciella parvula and of the genes
# of other species that take reads from them (crossmapping).
#
# Examples (genes ranked by mean TPM of F. parvula in the Tara samples, Carradec et al. 2022):
#   rank 11  Calmodulin   F. parvula in 082_0o8-5_SUR, with its Pelagomonas calceolata counterpart
#            (lcl|CAKKNE010000005.1_cds_CAH0375635.1_11956) in 004_0o8-5_SUR
#   rank 238 RPL26        F. parvula in 082_0o8-5_SUR, with its Florenciella sp. RCC1007 counterpart
#            (CAMNT_0053288775) in 004_0o8-5_SUR
#   rank 22  RplB         F. parvula in 084_0o8-5_SUR (no crossmapping)
# The counterparts are the genes that took most of the crossmapped reads (genes_cross_set of
# data/intermediate/fig3/crossmapping/<comparison>.rds, written by 03_crossmapping.R).
# NOTE: as in the original figure, the counterpart coverage comes from sample 004_0o8-5_SUR and the
# F. parvula coverage from 082/084_0o8-5_SUR.
#
# Input: BAMs of the quantification pipeline (not copied into this project), read with
#        samtools view <region> | bedtools genomecov -bga
# Output: data/intermediate/fig3/coverage_examples.tsv
#         (example, species, genename, from, to, coverage)
# Run time: ~1.5 min. Needs samtools and bedtools.
# Original: delving_into_florenciella src/figures/combined_florenciella_figure.R (panel C section)

source('R/load.R')

bam_dir <- file.path(external$delving_florenciella, 'data/quantification/alignment')

# F. parvula genes ranked by mean TPM
rank_genes <- read_tsv(p('quant_alignment', 'TPM_EP00618_Florenciella_parvula_quant-over_2012_carradec_tara.tsv'),
                       show_col_types = FALSE) %>%
    pivot_longer(names_to = 'gene_name', values_to = 'tpm', cols = -sample) %>%
    group_by(gene_name) %>%
    summarize(mean.tpm = mean(tpm)) %>%
    arrange(-mean.tpm) %>%
    mutate(rank = row_number())
gene_at <- function(r) rank_genes$gene_name[rank_genes$rank == r]

regions <- tribble(
    ~example, ~species,                           ~bam,                                                         ~gene,
    'rank11',  'EP00618_Florenciella_parvula',     'EP00618_Florenciella_parvula/082_0o8-5_SUR.filtered.bam.sorted', gene_at(11),
    'rank11',  'Pelagomonas_calceolata',           'Pelagomonas_calceolata/004_0o8-5_SUR.filtered.bam',              'lcl|CAKKNE010000005.1_cds_CAH0375635.1_11956',
    'rank238', 'EP00618_Florenciella_parvula',     'EP00618_Florenciella_parvula/082_0o8-5_SUR.filtered.bam.sorted', gene_at(238),
    'rank238', 'EP00619_Florenciella_sp_RCC1007',  'EP00619_Florenciella_sp_RCC1007/004_0o8-5_SUR.filtered.bam.sorted', 'CAMNT_0053288775',
    'rank22',  'EP00618_Florenciella_parvula',     'EP00618_Florenciella_parvula/084_0o8-5_SUR.filtered.bam.sorted', gene_at(22)
)

# Coverage of one gene (bedgraph: one row per run of bases with the same coverage, zeros included)
gene_coverage <- function(bam, gene) {
    cmd <- str_glue("samtools view -b {shQuote(file.path(bam_dir, bam))} {shQuote(gene)} | bedtools genomecov -bga -ibam -")
    data.table::fread(cmd = cmd, col.names = c('genename', 'from', 'to', 'coverage')) %>%
        as_tibble() %>%
        filter(genename == gene)
}

coverage <- regions %>%
    mutate(data = map2(bam, gene, gene_coverage)) %>%
    select(example, species, data) %>%
    unnest(data)

dir.create(p('intermediate', 'fig3'), showWarnings = FALSE, recursive = TRUE)
write_tsv(coverage, p('intermediate', 'fig3', 'coverage_examples.tsv'))
coverage %>% count(example, species) %>% print()
