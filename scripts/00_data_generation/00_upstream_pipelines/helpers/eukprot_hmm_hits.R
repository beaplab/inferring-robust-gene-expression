# Helper of 03a/b/c: protein sequences of the EukProt v01 hits of a Pfam HMM.
#
# Reads the hmmsearch tables of all proteomes, already concatenated in one file with the proteome name
# as first column (`grep --with-filename` in 03a/b/c): proteome <tab> hmmsearch --tblout line.
# Keeps hits with full-sequence E-value <= 1e-10. With --metazoa N, keeps every non-metazoan hit and N
# random metazoan ones (Supergroup_UniEuk of the EukProt table), so animals do not swamp the tree.
# Writes them as fasta, named <proteome>_<protein id>.
#
# Usage (from the project root):
#   Rscript scripts/00_data_generation/00_upstream_pipelines/helpers/eukprot_hmm_hits.R \
#       --hits <eukprot_all_hits.txt> --out <fasta> [--metazoa 20] [--seed N]
# Original: gene_explorer scripts/data_generation/generate_nonHK_gene_trees/helper_parse-<gene>-phylogeny.R
#   (one copy per gene: Flavodoxin and NRT2 subsampled 20 metazoan hits, SIT kept all). The originals
#   subsampled without a seed, so the reference sets cannot be regenerated exactly; --seed is new.

suppressPackageStartupMessages({
    library(tidyverse)
    library(Biostrings)
    library(argparser)
})
source('config/paths.R')

args <- arg_parser('EukProt sequences of the hits of an hmmsearch') %>%
    add_argument('--hits', help = 'proteome <tab> hmmsearch --tblout line') %>%
    add_argument('--out', help = 'output fasta') %>%
    add_argument('--metazoa', help = 'number of metazoan hits kept at random (0 = all hits, no subsampling)',
                 default = 0) %>%
    add_argument('--seed', help = 'seed for the metazoan subsample (none in the original runs)', default = NA) %>%
    parse_args()

max_evalue <- 1e-10

hits <- read_tsv(args$hits, col_names = c('proteome', 'hits'), show_col_types = FALSE) %>%
    separate(col = hits, sep = ' +', extra = 'drop',
             into = c('target', 'accession', 'query_name', 'pfam', 'fullseq_eval', 'fullseq_score', 'fullseq_bias')) %>%
    mutate(fullseq_eval = as.numeric(fullseq_eval)) %>%
    filter(fullseq_eval <= max_evalue)

if (args$metazoa > 0) {
    eukprot <- read_tsv(p('gene_data', 'EukProt_included_data_sets.v01.2019_05_08.txt'), show_col_types = FALSE)
    hits <- hits %>%
        mutate(EukProt_ID = str_remove(proteome, '_.*')) %>%
        left_join(eukprot, by = 'EukProt_ID')
    if (!is.na(args$seed)) set.seed(as.integer(args$seed))
    hits <- bind_rows(filter(hits, Supergroup_UniEuk != 'Metazoa'),
                      hits %>%
                          filter(Supergroup_UniEuk == 'Metazoa') %>%
                          distinct(target, .keep_all = TRUE) %>%
                          slice_sample(n = args$metazoa))
}

# Sequences of the hits of one proteome, named <proteome>_<protein id>
proteome_hits <- function(proteome, targets) {
    seqs <- readAAStringSet(file.path(external$eukprot_proteins, str_c(proteome, '.fasta')))
    names(seqs) <- str_remove(names(seqs), ' .*')
    seqs <- seqs[names(seqs) %in% targets]
    names(seqs) <- str_c(proteome, '_', names(seqs))
    seqs
}

sequences <- hits %>%
    group_by(proteome) %>%
    summarize(targets = list(target)) %>%
    { map2(.$proteome, .$targets, proteome_hits) } %>%
    do.call(c, .)
writeXStringSet(sequences, args$out, width = 80)
message(length(sequences), ' sequences from ', n_distinct(hits$proteome), ' proteomes written to ', args$out)
