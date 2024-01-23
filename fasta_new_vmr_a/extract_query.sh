awk '(FNR==1){print FILENAME,$0}' */*.fa  | sed 's/\//\t/;s/.fa >/\t/;s/ /\t/;s/ /\t/;' > query.tsv
