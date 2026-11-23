#!/usr/bin/env bash
# Upstream 03c: NRT2 (high-affinity nitrate transporter) reference tree.
#
# 1. Candidates: hmmsearch of PF07690 (Major Facilitator Superfamily, the domain of NRT2; E <= 1e-10) over
#    the EukProt v01 proteomes; all non-metazoan hits + 20 random metazoan ones; >= 100 aa
# 2. MFS is a huge family, so only the hits > 70% identical (blastp) to a NRT2 of Ocaña-Pallarès et al.
#    (fig5_NRT2_completeseq.fasta; https://pubmed.ncbi.nlm.nih.gov/30789903/) are kept
# 3. eggNOG annotation of the reference NRT2s; preliminary FastTrees (MAFFT -> trimAl -gt 0.3 -> FastTree)
#    with and without the reference NRT2s
# 4. IQ-TREE (ModelFinder, 100 standard bootstraps) on the EukProt sequences only -> tree/NRT2_iqtree_eukprotonly
#    (the tree used; an IQ-TREE with the reference NRT2s was also made, not used)
# 5. Tree and trimmed alignment copied as the phyloplace references
#
# Input:  EukProt v01 proteomes (external$eukprot_proteins), eggNOG db (external$eggnog_db)
#         data/raw/gene_data/reference_trees/NRT2/fig5_NRT2_completeseq.fasta (from the paper, only exists there)
# Output: data/raw/gene_data/reference_trees/NRT2/ (tree/NRT2_iqtree_eukprotonly.treefile, ...)
#         data/raw/placement/phyloplace_references/nonhk_genes/NRT2_{iqtree.treefile,msa.fasta}
# Run time: hours. Needs hmmsearch, seqkit, blastp, emapper.py, mafft, trimal, FastTree, iqtree.
# Original: gene_explorer scripts/data_generation/generate_nonHK_gene_trees/obtain_NRT2_phylogeny.sh
#           + helper_parse-NRT2-phylogeny.R (helpers/eukprot_hmm_hits.R here)
#   The folder also has PF00909.hmm and selection_itol_from_eukprot_w_ocañaseqs.txt, which the original
#   script does not use.
#
# Checked (2026-11-20): the phyloplace alignment (NRT2_msa.fasta) is this trimmed alignment with "|" in
#       the sequence names replaced by "_" (tree leaf labels cannot hold "|"); sequences identical.
# NOT RUN IN THIS PROJECT, and not reproducible as such (metazoan subsample without a seed, IQ-TREE seeded
# with the date). If run, it writes into gene_explorer through the reference_trees symlink.
#
# Usage (from the project root): bash scripts/00_data_generation/00_upstream_pipelines/03c_reference_tree_nrt2.sh

set -euo pipefail

GENE=NRT2
PFAM=PF07690
OUT=data/raw/gene_data/reference_trees/${GENE}
EUKPROT=/scratch/data4/aauladell/eukprot_v01/proteins
EGGNOG_DB=/scratch/data1/aauladell/databases/eggnog
HELPERS=scripts/00_data_generation/00_upstream_pipelines/helpers
mkdir -p ${OUT}/eukprot_hmm_search ${OUT}/tree

# 1. Candidates -----------------------------------------------------------------------------------------
wget "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/${PFAM}?annotation=hmm" -O ${OUT}/${PFAM}.hmm.gz
gzip -df ${OUT}/${PFAM}.hmm.gz

for proteome in ${EUKPROT}/E*; do
    hmmsearch -E 1e-10 --tblout ${OUT}/eukprot_hmm_search/$(basename ${proteome} .fasta).txt \
        ${OUT}/${PFAM}.hmm ${proteome}
done
grep --with-filename -v '#' ${OUT}/eukprot_hmm_search/*.txt |
    sed "s@${OUT}/eukprot_hmm_search/@@" | sed 's/.txt:/\t/' > ${OUT}/eukprot_all_hits.txt

Rscript ${HELPERS}/eukprot_hmm_hits.R --hits ${OUT}/eukprot_all_hits.txt \
    --out ${OUT}/${GENE}_eukprotv01_aa.fasta --metazoa 20
seqkit seq -m 100 ${OUT}/${GENE}_eukprotv01_aa.fasta > ${OUT}/${GENE}_eukprotv01_100aamin.fasta

# 2. Only MFS hits close to a known NRT2 ------------------------------------------------------------------
blastp -query ${OUT}/fig5_NRT2_completeseq.fasta -subject ${OUT}/${GENE}_eukprotv01_100aamin.fasta \
    -outfmt 6 -out ${OUT}/blast_comparison_${PFAM}_vs_edu.txt
awk '$3 > 70' ${OUT}/blast_comparison_${PFAM}_vs_edu.txt | cut -f2 > ${OUT}/filtered_version_blast.txt
seqkit grep -f ${OUT}/filtered_version_blast.txt ${OUT}/${GENE}_eukprotv01_100aamin.fasta > ${OUT}/${GENE}_eukprotv01_100min_sub.fasta

# 3. Annotation and preliminary trees -----------------------------------------------------------------------
emapper.py --cpu 12 -i ${OUT}/fig5_NRT2_completeseq.fasta --data_dir ${EGGNOG_DB} -o ${OUT}/annotation_ocaña.tsv
python3 ${HELPERS}/align_trim_fasttree.py -i ${OUT}/${GENE}_eukprotv01_100min_sub.fasta ${OUT}/fig5_NRT2_completeseq.fasta \
    -o ${OUT}/preliminar_results -n eukprot_w_ocañaseqs
python3 ${HELPERS}/align_trim_fasttree.py -i ${OUT}/${GENE}_eukprotv01_100min_sub.fasta \
    -o ${OUT}/preliminar_results -n eukprot_only
ALIGNMENT=${OUT}/preliminar_results/eukprot_only_trimed.fasta

# 4. Final tree (EukProt sequences only) ---------------------------------------------------------------
iqtree -s ${ALIGNMENT} -st AA -pre ${OUT}/tree/${GENE}_iqtree_eukprotonly \
    -b 100 -seed $(date +%s) -mem 10GB -nt 20 -wbtl

# 5. phyloplace references ------------------------------------------------------------------------------
REFS=data/raw/placement/phyloplace_references/nonhk_genes
cp ${OUT}/tree/${GENE}_iqtree_eukprotonly.treefile ${REFS}/${GENE}_iqtree.treefile
sed '/^>/ s/|/_/g' ${ALIGNMENT} > ${REFS}/${GENE}_msa.fasta   # "|" not allowed in leaf labels
