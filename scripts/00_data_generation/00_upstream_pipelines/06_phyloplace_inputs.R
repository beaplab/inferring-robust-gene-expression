# Upstream 06: inputs of nf-core/phyloplace (07): HMM profiles of the non-HK genes and the sample sheets.
#
# Each row of a sample sheet is one gene: target, hmm (profile used to find its unigenes in MATOU),
# refseqfile (reference alignment), refphylogeny (reference tree), model (LG+R8), taxonomy.
#   HK genes (318): alignments, trees (RAxML) and HMMs of Richter et al. (phyloplace_references/hk_genes)
#   non-HK genes (FAD1, Flavodoxin, NRT2, SIT): trees and alignments of 03a/b/c (FAD1 not in the paper);
#   their HMMs are built here with hmmbuild
# NOTE: in the original non-HK sheet the taxonomy column of every gene points to the ABCE taxonomy (a
#   placeholder); phyloplace's own taxonomic assignment of the non-HK genes is therefore not meaningful and
#   is not used (fig7/02 runs gappa examine assign with the right taxonomy). Kept as in the original.
# The 14 genes of the "no-premask" sheet were rerun with the EPA-ng --no-premask option (see 07).
#
# Input:  data/raw/placement/phyloplace_references/{hk_genes/{fastas,hmm,phylogenies}, nonhk_genes}
#         data/raw/placement/hk_gene_taxonomy/<gene>.taxonomy.txt
# Output: data/raw/placement/phyloplace_references/nonhk_genes/<gene>.hmm
#         data/raw/placement/phyloplace_samplesheets/{hk_genes,nonhk_genes,hk_genes_nopremask}.csv
# Run time: minutes. Needs hmmbuild (HMMER).
# Original: phylotranscriptomics scripts/helper-create-phyloplace-samplesheet_v2.R (HK sheet),
#           scripts/matou_phyloplace/creation_hmms_newgenes.sh (hmmbuild); the non-HK and no-premask sheets
#           (data/matou_search/nonhk_samplesheets.csv, with_incorrect_seqs_sheet.csv) were written by hand
#
# NOT RUN IN THIS PROJECT (the HMMs would overwrite the copies in phyloplace_references/nonhk_genes).
#
# Usage (from the project root): Rscript scripts/00_data_generation/00_upstream_pipelines/06_phyloplace_inputs.R

source('R/load.R')

refs <- p('placement', 'phyloplace_references')
sheet_dir <- p('placement', 'phyloplace_samplesheets')
dir.create(sheet_dir, showWarnings = FALSE)

# HMM profiles of the non-HK reference alignments -----------------------------------------------------
nonhk_genes <- c('FAD1', 'Flavodoxin', 'NRT2', 'SIT')
for (gene in nonhk_genes) {
    run_command(str_glue('hmmbuild -n {gene} {refs}/nonhk_genes/{gene}.hmm {refs}/nonhk_genes/{gene}_msa.fasta'))
}

# Sample sheets --------------------------------------------------------------------------------------------
hk_sheet <- tibble(target = str_remove(list.files(file.path(refs, 'hk_genes', 'fastas')), '\\.fasta$')) %>%
    mutate(hmm = file.path(refs, 'hk_genes', 'hmm', str_c(target, '.hmm')),
           refseqfile = file.path(refs, 'hk_genes', 'fastas', str_c(target, '.fasta')),
           refphylogeny = file.path(refs, 'hk_genes', 'phylogenies', str_c(target, '.fasta.raxml.bestTree')),
           model = 'LG+R8',
           taxonomy = p('placement', 'hk_gene_taxonomy', str_c(target, '.taxonomy.txt')))
write_csv(hk_sheet, file.path(sheet_dir, 'hk_genes.csv'))

# genes whose placement was rerun without premasking (07)
nopremask <- c('ABCE', 'HSP70C', 'MAT1A', 'RBM19', 'RFC2', 'RPS21', 'calm', 'gnb2L1', 'hsp90', 'psmc6',
               'rpl11', 'rpl3', 'rpl35a', 'rpl36a')
write_csv(filter(hk_sheet, target %in% nopremask), file.path(sheet_dir, 'hk_genes_nopremask.csv'))

tibble(target = nonhk_genes) %>%
    mutate(hmm = file.path(refs, 'nonhk_genes', str_c(target, '.hmm')),
           refseqfile = file.path(refs, 'nonhk_genes', str_c(target, '_msa.fasta')),
           refphylogeny = file.path(refs, 'nonhk_genes', str_c(target, '_iqtree.treefile')),
           model = 'LG+R8',
           taxonomy = p('placement', 'hk_gene_taxonomy', 'ABCE.taxonomy.txt')) %>%   # placeholder, see NOTE
    write_csv(file.path(sheet_dir, 'nonhk_genes.csv'))
