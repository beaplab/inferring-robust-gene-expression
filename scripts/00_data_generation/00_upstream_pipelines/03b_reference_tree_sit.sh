#!/usr/bin/env bash
# Upstream 03b: SIT (silicon transporter) reference tree.
#
# 1. Candidates: hmmsearch of the SIT Pfam domain PF03842 (E <= 1e-10) over the EukProt v01 proteomes;
#    all hits (no metazoan subsampling); >= 100 aa
# 2. Sequences that broke the tree in a first look in iTOL removed (SIT_incorrect_seqs.txt, by hand)
# 3. Alignment MAFFT L-INS-i (--localpair --maxiterate 100), trimAl -gt 0.3
# 4. IQ-TREE (ModelFinder, 100 standard bootstraps) -> tree/SIT_iqtree_cleaned_nogappyout
# 5. Tree and trimmed alignment copied as the phyloplace references
#
# Input:  EukProt v01 proteomes (external$eukprot_proteins)
#         data/raw/gene_data/reference_trees/SIT/SIT_incorrect_seqs.txt (made by hand, only exists there)
# Output: data/raw/gene_data/reference_trees/SIT/ (tree/SIT_iqtree_cleaned_nogappyout.treefile, ...)
#         data/raw/placement/phyloplace_references/nonhk_genes/SIT_{iqtree.treefile,msa.fasta}
# Run time: hours. Needs hmmsearch, seqkit, mafft, trimal, iqtree.
# Original: gene_explorer scripts/data_generation/generate_nonHK_gene_trees/obtain_SIT_phylogeny.sh
#           + helper_parse-SIT-phylogeny.R (helpers/eukprot_hmm_hits.R here; checked: same 950 sequences)
#   Left out: a first IQ-TREE on the uncleaned alignment (default MAFFT), replaced by this one; a commented
#   DNA tree from the Durkin et al. 2016 sequences; the EukProt download (figshare 12417881).
#
# NOT RUN IN THIS PROJECT (IQ-TREE was seeded with the date, so the tree cannot be regenerated exactly).
# If run, it writes into gene_explorer through the reference_trees symlink.
# Checked (2026-11-20): the phyloplace alignment (<gene>_msa.fasta) is this trimmed alignment with "|" in
#       the sequence names replaced by "_" (tree leaf labels cannot hold "|"); sequences identical.
#
# Usage (from the project root): bash scripts/00_data_generation/00_upstream_pipelines/03b_reference_tree_sit.sh

set -euo pipefail

GENE=SIT
PFAM=PF03842
OUT=data/raw/gene_data/reference_trees/${GENE}
EUKPROT=/scratch/data4/aauladell/eukprot_v01/proteins
HELPERS=scripts/00_data_generation/00_upstream_pipelines/helpers
mkdir -p ${OUT}/eukprot_hmm_search ${OUT}/tree

# 1. Candidates -----------------------------------------------------------------------------------------
wget "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/${PFAM}?annotation=hmm" -O ${OUT}/${PFAM}.hmm.gz
gzip -df ${OUT}/${PFAM}.hmm.gz

for proteome in ${EUKPROT}/*; do
    hmmsearch -E 1e-10 --tblout ${OUT}/eukprot_hmm_search/$(basename ${proteome} .fasta).txt \
        ${OUT}/${PFAM}.hmm ${proteome}
done
grep --with-filename -v '#' ${OUT}/eukprot_hmm_search/*.txt |
    sed "s@${OUT}/eukprot_hmm_search/@@" | sed 's/.txt:/\t/' > ${OUT}/eukprot_all_hits.txt

Rscript ${HELPERS}/eukprot_hmm_hits.R --hits ${OUT}/eukprot_all_hits.txt --out ${OUT}/${GENE}_eukprotv01_aa.fasta
seqkit seq -m 100 ${OUT}/${GENE}_eukprotv01_aa.fasta > ${OUT}/${GENE}_eukprotv01_100aamin.fasta

# 2-3. Cleaning and alignment --------------------------------------------------------------------------------
FASTA=${OUT}/${GENE}_eukprotv01_100aamin.fasta
CLEAN=${FASTA/.fasta/_cleaned.fasta}
ALIGNED=${CLEAN/.fasta/_align.fasta}
ALIGNMENT=${ALIGNED/.fasta/_trimed.fasta}

seqkit grep -n -v -f ${OUT}/SIT_incorrect_seqs.txt ${FASTA} > ${CLEAN}
mafft --localpair --maxiterate 100 ${CLEAN} > ${ALIGNED}
trimal -in ${ALIGNED} -out ${ALIGNMENT} -gt 0.3

# 4. Tree ----------------------------------------------------------------------------------------------
iqtree -s ${ALIGNMENT} -st AA -pre ${OUT}/tree/${GENE}_iqtree_cleaned_nogappyout \
    -b 100 -seed $(date +%s) -mem 10GB -nt 15 -wbtl

# 5. phyloplace references ------------------------------------------------------------------------------
REFS=data/raw/placement/phyloplace_references/nonhk_genes
cp ${OUT}/tree/${GENE}_iqtree_cleaned_nogappyout.treefile ${REFS}/${GENE}_iqtree.treefile
sed '/^>/ s/|/_/g' ${ALIGNMENT} > ${REFS}/${GENE}_msa.fasta
