#!/usr/bin/env bash
#
# show blast hits for a query accessions
#
if [ -z "$1" ]; then
    echo "SYNTAX: $0 QUERY_ACCESSION [hit_species_name]"
    exit 1
fi

# parse args
QACC=$1
shift

if [ ! -z "$1" ]; then 
    OR_HIT_SPECIES_NAME="|$1"
    shift
fi

# scan
(
    cat results/headers.query_hits_counts.tsv
    egrep -n "($QACC$OR_HIT_SPECIES_NAME)" results/*/a/summary.bitscore.query_hits_counts.tsv
) \
    | cut -f 1,3-7,8,9,17 \
    | column -t \
    | egrep -B 1 --color=yes "($QACC|q_accession|blastn[0-9]+$OR_HIT_SPECIES_NAME)"
