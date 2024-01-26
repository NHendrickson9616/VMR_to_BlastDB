#!/usr/bin/env bash
#
# join mismatch-by-genus tables across all runs
#
SORTS="bitscore evalue raw"
for SORT in $SORTS; do

    echo "# -------------------------"
    echo "# SORT: $SORT"
    echo "# -------------------------"

    FINAL=mismatches_by_genus_type.$SORT.all.tsv

    # first two
    L1=blastn10
    L2=blastn11
    IN1=$L1/a/summary.$SORT.mismatch_by_genus.tsv
    IN2=$L2/a/summary.$SORT.mismatch_by_genus.tsv
    TMP=.merge_by_genus.$L2.tsv

    echo -e "KEY\t$L1\t$L2" > $TMP
    join -t "$(echo -e "\t")" -e 0 -a 1 -a 2 -o "0 1.4 2.4" $IN1 $IN2  >> $TMP

    head $TMP

    #
    # iterate
    #
    LABELS="blastn13 blastn15 blastn16 blastn17 blastn18  blastn20  blastn22  blastn24 blastn28"

    for L2 in $LABELS; do
	echo "# ------------- $L ----------------"
	IN1=$TMP
	IN2=$L2/a/summary.$SORT.mismatch_by_genus.tsv
	TMP=.merge_by_genus.$L2.tsv

	join -t "$(echo -e "\t")" -e 0 -a 1 -a 2 \
	    <(sort -k1,1 $IN1) \
	    <((echo -e "KEY\t$L2";cut -f 1,4 $IN2)|sort -k1,1) \
	    > $TMP

	(grep KEY $TMP;head -2 $TMP)|column -t
    done

    #
    # finalize
    #
    cp $TMP $FINAL
    wc -l $FINAL
done
