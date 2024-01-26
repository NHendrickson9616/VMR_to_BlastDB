#!/usr/bin/env bash
#
# 
#

echo "# check counts"
wc -l test_out.hits.tsv test_out.match.tsv  ../../fasta_new_vmr_a/query.tsv 


#
# problem in test_out.match.tsv
# when no hits, we're sucking in the comment after the hit line.x

# 
TAB=$(echo -e "\t")
echo "TEST TAB: a${TAB}b"
echo "#"
echo "# test_out.hits.tsv"
sort -k2,2 test_out.hits.tsv | sed 's/\t/[_]/g' | head -5
echo "#"
echo "# test_out.match.tsv "
sort -k3,3 test_out.match.tsv |  sed 's/\t/[_]/g' | head -5
echo "#"
echo "# join"
join   -t "$TAB" -1 2 -2 3  <(sort -k2,2 test_out.hits.tsv) <(sort -k3,3 test_out.match.tsv)| sort -k1,1  >  test_out.merged.tsv
head test_out.merged.tsv|  sed 's/\t/[_]/g' #|column -t

echo "#"
echo "# query"
sort -k2,2 ../../fasta_new_vmr_a/query.tsv | head -5

echo "#"
echo "# merge query and results"
join -t "$TAB" -1 2 -2 1  <(sort -k2,2 ../../fasta_new_vmr_a/query.tsv) <(cat test_out.merged.tsv)| sort -k1,1  >  test_out.query.merged.tsv 
head test_out.query.merged.tsv
