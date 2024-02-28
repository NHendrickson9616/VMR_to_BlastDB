#!/usr/bin/env bash
#
#
# try different blast approachs for a problem genome
#
#SBATCH --job-name=ICTV_BLAST_test_Kayvirus_JX08032_vs_db
#SBATCH --output=logs/log.%J.%x.out
#SBATCH --error=logs/log.%J.%x.err
#
# Number of tasks needed for this job. Generally, used with MPI jobs
# Time format = HH:MM:SS, DD-HH:MM:SS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000  # in Megabytes
# 
# partition and time
#SBATCH --partition=amd-hdr100 --time=00-12:00:00
##SBATCH --partition=amd-hdr100 --time=06-06:00:00
##SBATCH --partition=medium --time=40:00:00

OUT_DIR=results/blastn10_test/a/Keyvirus
QUERY_FA=./fasta_new_vmr_a/Kayvirus/JX080302.fa

if [ -z "$(which blastn 2>/dev/null)" ]; then
    echo module load Anaconda3
    module load Anaconda3
    echo conda activate conda/vmr_openpyxl3
    conda activate conda/vmr_openpyxl3
fi

mkdir -p $OUT_DIR

# F="results/blastn10_test/a/Keyvirus/JX080302.5.maxtargetseqs10.maxhsps10.xml"; if [ ! -s "$F" ]; then echo "empty/non-exist"; elif [ ! -z "$(grep '</BlastOutput>' $F)" ]; then echo "full"; else echo "partial"; fi

for FMT in 0 1 2 3 4 5 6 7 "7 qseqid qlen sseqid slen bitscore length" "7 qseqid qlen sseqid slen bitscore length qcovs" 8 9 10 11 ; do
    FMT_SAFE=$(echo -n $FMT | sed 's/ /_/g')
    echo "#"
    echo "# FMT=$FMT"
    echo "# FMT_SAFE==$FMT_SAFE"
    echo "#"
    LOG=logs/test_Kayvirus_JX080302_${FMT_SAFE}.txt
    OUT=$OUT_DIR/JX080302.${FMT_SAFE}.txt

    if [ -s "$OUT" ]; then
	echo "FMT=$FMT  SKIP: exists & non-empty ($OUT)"
    else
	echo blastn -task blastn -word_size 10 -db ./blast/ICTV_VMR_e -query $QUERY_FA -out $OUT -outfmt "$FMT" 2\>\&1 \| tee $LOG
             blastn -task blastn -word_size 10 -db ./blast/ICTV_VMR_e -query $QUERY_FA -out $OUT -outfmt "$FMT" 2>&1    | tee $LOG
    fi
done

for FMT in 0 7 "7 qseqid qlen sseqid slen length"; do
    FMT_SAFE=maxhsps1.$(echo -n $FMT | sed 's/ /_/g')
    echo "#"
    echo "# FMT=$FMT"
    echo "# FMT_SAFE==$FMT_SAFE"
    echo "#"
    LOG=logs/test_Kayvirus_JX080302_${FMT_SAFE}.txt
    OUT=$OUT_DIR/JX080302.${FMT_SAFE}.txt

    if [ -s "$OUT" ]; then
	echo "FMT=$FMT  SKIP: exists & non-empty ($OUT)"
    else
	echo blastn -task blastn -word_size 10 -db ./blast/ICTV_VMR_e -query $QUERY_FA -out $OUT -max_hsps 1 -outfmt "$FMT" 2\>\&1 \| tee $LOG
             blastn -task blastn -word_size 10 -db ./blast/ICTV_VMR_e -query $QUERY_FA -out $OUT -max_hsps 1 -outfmt "$FMT" 2>&1    | tee $LOG
    fi
done

FMT_ARGS="-outfmt 5  -max_target_seqs 10 -max_hsps 10"
for FMT_SAFE in 5.maxtargetseqs10.maxhsps10; do
    echo "#"
    echo "# FMT=$FMT"
    echo "# FMT_SAFE==$FMT_SAFE"
    echo "#"
    LOG=logs/test_Kayvirus_JX080302_${FMT_SAFE}.txt
    OUT=$OUT_DIR/JX080302.${FMT_SAFE}.xml

    if [ -s "$OUT" ]; then
	echo "FMT=$FMT  SKIP: exists & non-empty ($OUT)"
    else
	echo blastn -task blastn -word_size 10 -db ./blast/ICTV_VMR_e -query $QUERY_FA -out $OUT $FMT_ARGS 2\>\&1 \| tee $LOG
             blastn -task blastn -word_size 10 -db ./blast/ICTV_VMR_e -query $QUERY_FA -out $OUT $FMT_ARGS 2>&1    | tee $LOG
    fi
done

FMT_ARGS="-outfmt 5  -max_target_seqs 10"
for FMT_SAFE in 5.hit10; do
    echo "#"
    echo "# FMT=$FMT"
    echo "# FMT_SAFE==$FMT_SAFE"
    echo "#"
    LOG=logs/test_Kayvirus_JX080302_${FMT_SAFE}.txt
    OUT=$OUT_DIR/JX080302.${FMT_SAFE}.xml

    if [ -s "$OUT" ]; then
	echo "FMT=$FMT  SKIP: exists & non-empty ($OUT)"
    else
	echo blastn -task blastn -word_size 10 -db ./blast/ICTV_VMR_e -query $QUERY_FA -out $OUT $FMT_ARGS 2\>\&1 \| tee $LOG
             blastn -task blastn -word_size 10 -db ./blast/ICTV_VMR_e -query $QUERY_FA -out $OUT $FMT_ARGS 2>&1    | tee $LOG
    fi
done

