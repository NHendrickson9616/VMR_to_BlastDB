#!/usr/bin/env bash
#
# summarize results for each blast setting
#
mkdir -p logs

SORTS="raw bitscore evalue"
echo SORTS=\"$SORTS\"

MODES=$(ls -d */a | cut -d / -f 1)
echo MODES=\"$MODES\"

for SORT in $SORTS; do
    echo "# -------------------"
    echo "# SORT=$SORT"
    echo "# -------------------"
    for MODE in $MODES; do
	echo "  # MODE=$MODE"
	echo "  "sbatch --job-name ICTV_BLAST_extract_results-${MODE}-${SORT} ./extract_results.sh -m $MODE -s $SORT
	sbatch --job-name ICTV_BLAST_extract_results-${MODE}-${SORT} ./extract_results.sh -m $MODE -s $SORT
    done
done
