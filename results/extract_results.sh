#!/usr/bin/env bash
#
# Analyzed Blast results (-fmt 7) from A sequences vs E-seq-BASTdb
#
# as generated by ../test_e_accessions.sh
#
#SBATCH --job-name=ICTV_BLAST_extract_results
#SBATCH --output=logs/log.%J.%x.out
#SBATCH --error=logs/log.%J.%x.out
#
# Number of tasks needed for this job. Generally, used with MPI jobs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=500  # in Megabytes
# 
# partition and time
# Time format = HH:MM:SS, DD-HH:MM:SS
#SBATCH --partition=amd-hdr100 --time=00-00:20:00
##SBATCH --partition=amd-hdr100 --time=06-06:00:00
##SBATCH --partition=medium --time=40:00:00
#
TAB=$(echo -ne "\t")
BLAST_ABBREV="default"
MODE="bitscore"  # bitscore, raw



# parse args
if [[ "$1" == "-m" && ! -z "$1"  ]]; then BLAST_ABBREV=$2; echo BLAST_ABBREV=$BLAST_ABBREV; shift 2; fi
if [[ "$1" == "-s" && ! -z "$1"  ]]; then MODE=$2; echo MODE=$MODE; shift 2; fi

# build DIR paths
RES_DIR=./$BLAST_ABBREV/a
echo "RES_DIR=$RES_DIR"

# files
QUERY_INFO=../fasta_new_vmr_a/query.tsv

if [ $MODE == "bitscore" ]; then HITS_SUFFIX="bitscore.csv"; fi
if [ $MODE == "evalue" ];   then HITS_SUFFIX="evalue.csv"; fi
if [ $MODE == "raw" ]; then      HITS_SUFFIX="raw.txt"; fi
echo MODE=$MODE
echo HITS_SUFFIX=$HITS_SUFFIX

HIT_COUNTS=$RES_DIR/summary.$MODE.hits.tsv
NO_HITS=$RES_DIR/summary.$MODE.no_hits.tsv
TOP_HITS=$RES_DIR/summary.$MODE.match.tsv
MISMATCH_SUMMARY=$RES_DIR/summary.$MODE.mismatch_by_genus.tsv
MISMATCH_ONELINE=$RES_DIR/summary.$MODE.mismatch_totals.tsv
MERGE_HITS_COUNTS=$RES_DIR/summary.$MODE.hits_counts.tsv
MERGE_QUERY_HITS_COUNTS=$RES_DIR/summary.$MODE.query_hits_counts.tsv


echo "# --------------------------------------------------"
echo "# grep for hit counts"
#
# output lines from grep look like
# ./default/a/Agatevirus/KJ010548.bitscore.csv:# 38 hits found
#
(grep "hits found" $RES_DIR/*/*.$HITS_SUFFIX \
    | cut -d / --output-delimiter "$TAB" -f 4- \
    | sed "s/.$HITS_SUFFIX:# /\t/"';s/ hits found.*//;' \
    > $HIT_COUNTS) &

#
# output lines from grep look like
# ./default/a/Agatevirus/KJ010548.bitscore.csv:# 38 hits found
# ./default/a/Agatevirus/KJ010548.bitscore.csv-Agatevirus_Bp8pC,Nitunavirus_grass,76.041,2162,484,28,47809,49956,89704,87563,0.0,1092
# --
#
echo "# grep for top hit "
(fgrep -A 1 "hits found" $RES_DIR/*/*.$HITS_SUFFIX  \
    | fgrep ".${HITS_SUFFIX}-" \
    | cut -d / --output-delimiter "$TAB" -f 4- \
    | cut -d , --output-delimiter "$TAB" -f 1- \
    | sed "s/.${HITS_SUFFIX}-/\t/;" \
    | awk -f extract_results.awk \
    > $TOP_HITS) & 

#
echo "# wait, while those run in parallel"
#
wait
wc -l $HIT_COUNTS
wc -l $TOP_HITS

#
# summarize queries with no hits
#
echo "# --------------------------------------------------"
echo "# no_hit summary"
echo "#"
wc -l $HIT_COUNTS
awk 'BEGIN{GENUS=1;ACCESSION=2;HITS=3}($HITS==0){print $GENUS}' $HIT_COUNTS \
    | sort \
    | uniq -c \
    > $NO_HITS
echo -n "No hits: "
awk 'BEGIN{sum=0}{sum = sum + $1}END{print sum " queries in " NR " genera"}' $NO_HITS

echo "# --------------------------------------------------"
echo "# summarize mismatches"
echo "#"
echo "# $TOP_HITS -> $MISMATCH_SUMMARY"
cut -f 1,2 $TOP_HITS | grep -wv match |  sort -k2,2 -k1,1 | uniq -c | awk 'BEGIN{OFS="\t"}{print $3"."$2,$3,$2,$1}'| sort -k1,1 > $MISMATCH_SUMMARY
cat $MISMATCH_SUMMARY | awk 'BEGIN{KEY=1;GENUS=2;ERR=3;CT=4}{sum=sum+$CT;cts[$ERR]=cts[$ERR]+$CT}END{printf sum " total mismatches: ";for(i in cts){printf cts[i]" "i","};print ""}' > $MISMATCH_ONELINE
echo "Mismatches: " $(cat $MISMATCH_ONELINE)


echo "# --------------------------------------------------"
echo "# head $HIT_COUNTS"
echo "#"
sort -k2,2 $HIT_COUNTS | sed 's/\t/[_]/g' | head -5
echo "#"
echo "# head $TOP_HITS"
sort -k3,3 $TOP_HITS |  sed 's/\t/[_]/g' | head -5
echo "#"
echo "# join hit_counts and top_match: $MERGE_HITS_COUNTS"
echo "#"
join   -t "$TAB" -1 1 -2 3 <(sort -k2,2 $HIT_COUNTS| cut -f 2,3) <(sort -k3,3 $TOP_HITS)| sort -k1,1  >  $MERGE_HITS_COUNTS
head $MERGE_HITS_COUNTS |  sed 's/\t/[_]/g' #|column -t

echo "#"
echo "# query"
sort -k2,2 $QUERY_INFO | head -5

echo "#"
echo "# merge query and results"
join -t "$TAB" -1 1 -2 1  <(sort -k2,2 $QUERY_INFO| cut -f 2,5) <(cat $MERGE_HITS_COUNTS )| sort -k1,1  >  $MERGE_QUERY_HITS_COUNTS
head $MERGE_QUERY_HITS_COUNTS
