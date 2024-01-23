#!/usr/bin/env bash
#
# query each "A" (additional) sequence against the blastdb of "E" sequences. 
# tabulate matches and mis-matches
#
# (This would be much faster if we merged all the query fastas)
#
#
# 20240122 runtime 
#
#SBATCH --job-name=ICTV_BLAST_A_fastas_vs_E_db
#SBATCH --output=log.%J.%x.out
#SBATCH --error=log.%J.%x.out
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
#

QUERY_DIR=./fasta_new_vmr_a
echo "QUERY_DIR:   $QUERY_DIR"
RESULTS_DIR=./results/a
echo "RESULTS_DIR: $RESULTS_DIR"
mkdir -p $RESULTS_DIR
BLAST_DB=./blast/ICTV_VMR_e
echo "BLAST_DB:    $BLAST_DB"

if [ -z "$(which blastn 2>/dev/null)" ]; then
    echo "#module load BLAST"
    module load BLAST
fi


echo "#"
echo "# enumerate A seqs"
echo "# "
echo '# find fasta_new_vmr_a -name "*.fa"'
A_FASTAS=$(find fasta_new_vmr_a -name "*.fa")
echo '# found ' $(find fasta_new_vmr_a -name "*.fa"|wc -l) " A fastas"


echo "#"
echo "# run blast"
echo "#" 
for QUERY in $A_FASTAS; do
    
    GENUS=$(basename $(dirname $QUERY))
    ACCESS=$(basename $QUERY .fa)
    RESULT_DIR=$RESULTS_DIR/$GENUS
    RESULTS_FILE="$RESULT_DIR/$ACCESS.csv"

    if [[ -e $RESULTS_FILE && $RESULTS_FILE -nt $QUERY ]]; then
	echo "${GENUS}/${ACCESS}: SKIP"
    else
	echo -n "${GENUS}/${ACCESS}: blast " 
	mkdir -p $RESULT_DIR
	blastn -db $BLAST_DB -query $QUERY -out ${RESULTS_FILE} -outfmt '7 delim=,'    
	grep "hits found" $RESULTS_FILE
    fi
done
