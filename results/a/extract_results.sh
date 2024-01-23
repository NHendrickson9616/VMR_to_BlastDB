(grep "hits found" */*.csv \
    | sed 's/\//\t/;s/.csv:# /\t/;s/ hits found.*//;' \
    > test_out.hits.tsv) &


(fgrep -A 1 "hits found" */*.csv  \
    | fgrep ".csv-" \
    | sed 's/\//\t/;s/.csv-/\t/;s/,/\t/g;' \
    | awk '{if(NF==7){print "no_hits\t"$1"\t"$2;next;};if($3==$4){print "match\t"$0} else {print "mismatch\t"$0}}' \
    > test_out.match.tsv) & 
