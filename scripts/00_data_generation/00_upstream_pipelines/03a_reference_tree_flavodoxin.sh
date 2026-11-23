#!/usr/bin/env bash
# Upstream 03a: Flavodoxin reference tree (the tree on which the MATOU unigenes are placed, 07).
#
# 1. Candidates: hmmsearch of the flavodoxin Pfam domain PF00258 (E <= 1e-10) over the EukProt v01
#    proteomes; all non-metazoan hits + 20 random metazoan ones; >= 100 aa
# 2. eggNOG annotation of the candidates (guides the manual selection)
# 3. Literature sequences added: diatom flavodoxins of Graff et al. 2022 (eLife;
#    flavodoxins_graff_elife_2022.fasta) and the P. calceolata flavodoxins (adding_pelagomonas_seqs.tsv,
#    translated from its CDS)
# 4. Preliminary tree (MAFFT -> trimAl -gt 0.3 -> FastTree), inspected in iTOL with the eggNOG
#    annotation: the true flavodoxins were listed by hand in secondsel_w_eggnog_annot.tsv
# 5. Selected EukProt sequences + literature sequences -> MAFFT -> trimAl -> FastTree (preliminar_results_v2)
# 6. Final tree: IQ-TREE on that trimmed alignment (ModelFinder, 100 standard bootstraps)
# 7. Tree and trimmed alignment copied as the phyloplace references (06 builds the HMM from the alignment)
#
# Input:  EukProt v01 proteomes (external$eukprot_proteins), eggNOG db (external$eggnog_db)
#         data/raw/genomic/transcriptomes_nucleotide/Pelagomonas_calceolata.fna.gz
#         curated files in data/raw/gene_data/reference_trees/Flavodoxin/ (made by hand, only exist there):
#           flavodoxins_graff_elife_2022.fasta, adding_pelagomonas_seqs.tsv, secondsel_w_eggnog_annot.tsv
# Output: data/raw/gene_data/reference_trees/Flavodoxin/ (tree/Flavodoxin_iqtree.treefile, ...)
#         data/raw/placement/phyloplace_references/nonhk_genes/Flavodoxin_{iqtree.treefile,msa.fasta}
# Run time: hours. Needs hmmsearch, seqkit, emapper.py, mafft, trimal, FastTree, iqtree.
# Original: gene_explorer scripts/data_generation/generate_nonHK_gene_trees/obtain_Flavodoxin_phylogeny.sh
#           + helper_parse-Flavodoxin-phylogeny.R (helpers/eukprot_hmm_hits.R here)
#
# NOT RUN IN THIS PROJECT, and not reproducible as such: the metazoan subsample had no seed, step 4 is a
# manual selection, and IQ-TREE was seeded with the date. If run, it writes into gene_explorer through the
# reference_trees symlink.
# Checked (2026-11-20): the phyloplace alignment (<gene>_msa.fasta) is this trimmed alignment with "|" in
#       the sequence names replaced by "_" (tree leaf labels cannot hold "|"); sequences identical.
#
# Usage (from the project root): bash scripts/00_data_generation/00_upstream_pipelines/03a_reference_tree_flavodoxin.sh

set -euo pipefail

GENE=Flavodoxin
PFAM=PF00258
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
# one table: proteome <tab> hmmsearch line
grep --with-filename -v '#' ${OUT}/eukprot_hmm_search/*.txt |
    sed "s@${OUT}/eukprot_hmm_search/@@" | sed 's/.txt:/\t/' > ${OUT}/eukprot_all_hits.txt

Rscript ${HELPERS}/eukprot_hmm_hits.R --hits ${OUT}/eukprot_all_hits.txt \
    --out ${OUT}/${GENE}_eukprotv01_aa.fasta --metazoa 20
seqkit seq -m 100 ${OUT}/${GENE}_eukprotv01_aa.fasta > ${OUT}/${GENE}_eukprotv01_100aamin.fasta

# 2. eggNOG annotation ------------------------------------------------------------------------------------
emapper.py --cpu 12 -i ${OUT}/${GENE}_eukprotv01_100aamin.fasta --data_dir ${EGGNOG_DB} \
    -o ${OUT}/annotation_eukprotv01.tsv

# 3. Literature sequences ---------------------------------------------------------------------------------
seqkit grep -f ${OUT}/adding_pelagomonas_seqs.tsv \
    data/raw/genomic/transcriptomes_nucleotide/Pelagomonas_calceolata.fna.gz |
    seqkit translate > ${OUT}/pelago_flavo.fasta

# 4. Preliminary tree (inspected by hand -> secondsel_w_eggnog_annot.tsv) ------------------------------------
python3 ${HELPERS}/align_trim_fasttree.py \
    -i ${OUT}/${GENE}_eukprotv01_100aamin.fasta ${OUT}/flavodoxins_graff_elife_2022.fasta ${OUT}/pelago_flavo.fasta \
    -o ${OUT}/preliminar_results -n eukprot_w_graff-elife-2022_pelago

# 5. Selected sequences ----------------------------------------------------------------------------------
seqkit grep -f ${OUT}/secondsel_w_eggnog_annot.tsv ${OUT}/${GENE}_eukprotv01_100aamin.fasta > ${OUT}/flavodoxin_selected.fasta
python3 ${HELPERS}/align_trim_fasttree.py \
    -i ${OUT}/flavodoxin_selected.fasta ${OUT}/flavodoxins_graff_elife_2022.fasta ${OUT}/pelago_flavo.fasta \
    -o ${OUT}/preliminar_results_v2 -n subseteuk_wadditional-seqs-biblio
ALIGNMENT=${OUT}/preliminar_results_v2/subseteuk_wadditional-seqs-biblio_trimed.fasta

# 6. Final tree -----------------------------------------------------------------------------------------
iqtree -s ${ALIGNMENT} -st AA -pre ${OUT}/tree/${GENE}_iqtree \
    -b 100 -seed $(date +%s) -mem 10GB -nt 20 -wbtl

# 7. phyloplace references ------------------------------------------------------------------------------
REFS=data/raw/placement/phyloplace_references/nonhk_genes
cp ${OUT}/tree/${GENE}_iqtree.treefile ${REFS}/${GENE}_iqtree.treefile
sed '/^>/ s/|/_/g' ${ALIGNMENT} > ${REFS}/${GENE}_msa.fasta   # "|" not allowed in leaf labels
