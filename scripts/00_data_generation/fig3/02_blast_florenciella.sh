#!/usr/bin/env bash
# Fig 3B: which F. parvula genes have a counterpart in each other transcriptome / genome
# (white "gene absent" cells of the crossmapping heatmap).
#
# blastn of every F. parvula transcript against the nucleotide version of each other transcriptome,
# keeping hits with >= 80% identity over >= 100 nt.
#
# Input:  delving_into_florenciella data/genomic_data/transcriptomes/nucleotide_version/*.fna.gz
#         (not copied into this project)
# Output: data/intermediate/fig3/blast_Florenciella_allgenomes.tsv.gz
#         (query, subject, then the 12 standard blastn -outfmt 6 columns)
# Run time: ~1 h. Needs blastn (BLAST+).
# Original: delving_into_florenciella src/00_genome-down-and-selection/02_comparison-blast-genomes.sh
#           + src/analysis/helper_clean-blast-florenciella-vs-allgenomes.R
#   The original read the 12-column output with 17 column names (a leftover of a custom -outfmt), so
#   only qseqid, sseqid, pident and length had the right name; those were the only ones used.
#
# Usage (from the project root): bash scripts/00_data_generation/fig3/02_blast_florenciella.sh

set -euo pipefail

NUCL=/scratch/data3/aauladell/project_big/delving_into_florenciella/data/genomic_data/transcriptomes/nucleotide_version
QUERY=EP00618_Florenciella_parvula
OUT=data/intermediate/fig3/blast_Florenciella_allgenomes.tsv.gz
TMP=$(mktemp -d)

for subject_path in ${NUCL}/*.fna.gz; do
    subject=$(basename ${subject_path} .fna.gz)
    # F. parvula against itself is not a comparison
    [[ ${subject} == ${QUERY} ]] && continue
    blastn -query <(gzip -dc ${NUCL}/${QUERY}.fna.gz) \
           -subject <(gzip -dc ${subject_path}) \
           -outfmt 6 \
           -out ${TMP}/${subject}.tsv
done

# keep hits with >= 80% identity over >= 100 nt; prepend query and subject
{
    printf 'query\tsubject\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n'
    for f in ${TMP}/*.tsv; do
        awk -F'\t' -v OFS='\t' -v q=${QUERY} -v s=$(basename $f .tsv) \
            '$3 >= 80 && $4 >= 100 { print q, s, $0 }' $f
    done
} | gzip > ${OUT}

rm -r ${TMP}
