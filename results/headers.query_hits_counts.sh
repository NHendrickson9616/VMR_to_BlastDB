#!/usr/bin/env bash
#
# create header-line TSV file
#
echo "q_accession q_gbk_name hits best_hit_status genus qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore" | sed 's/ /\t/g' > headers.query_hits_counts.tsv

