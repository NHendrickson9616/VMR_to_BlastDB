#!/usr/bin/env bash
#
# given a new MSL
#  download from box
#  parse
#  pull fasta from NCBI & reformat
#  build E record blastdb
#  blast A records against E-blastdb
#  tabulate results

# better step logging if sbatch'ed
SRUN=""
if [ ! -z "$SLURM_JOB_ID" ]; then 
    SRUN=srun
    echo "SLURM_JOB_ID=$SLURM_JOB_ID"
    echo "seff $SLURM_JOB_ID"
    echo "sacctx -j $SLURM_JOB_ID"
fi
BLAST_MODES="blastn10 blastn11 blastn13"
BLAST_MODES=$BLAST_MODES

echo "#"
echo "# download VMRs from box into ./VMRs/"
echo "#"
echo  ./pull_VMRs_from_box.sh
$SRUN ./pull_VMRs_from_box.sh

echo "#"
echo "# get most recent MSL.xlsx file"
echo "#"

MSL_XLSX=$(ls -rt VMRs/VMR_MSL*.xlsx | tail -1)
if [ ! -z "$1" ]; then 
    MSL_XLSX=$1
fi
echo MSL_XLSX=$MSL_XLSX

echo "#"
echo "# parse VMR "
echo "#"
echo  ./VMR_to_fasta.py -mode VMR -VMR_file_name $VMR_XLSX
$SRUN ./VMR_to_fasta.py -mode VMR -VMR_file_name $VMR_XLSX

echo "#"
echo "# pull and format fastas"
echo "#"
echo "# remove formated fastas"
echo 'find . -name "*.fa" -exec rm {} +'
$SRUN find . -name "*.fa" -exec rm {} +

echo "# downlaod and format E fastas - in serial"
$SRUN ./download_fasta_files $VRM_XLSX

echo "# downlaod and format A fastas - in serial"
$SRUN ./download_fasta_files_a  $VRM_XLSX

echo "#"
echo "# build blast db: 3-4 minutes""
echo "#"
echo "makedatabase.sh"
$SRUN makedatabase.sh

echo "#"
echo "# test blast A records vs E_blastdb: 30-40 minutes, in parallel"
echo "#"
echo "sbatch_modes_by_genus_test_e.sh $BLAST_MODES"
$SRUN sbatch_modes_by_genus_test_e.sh $BLAST_MODES

echo "# wait for parallel blast jobs to complete"
echo "squeue_email ICTV_BLAST_A_fastas_vs_E_db_ done for $0"
$SRUN squeue_email ICTV_BLAST_A_fastas_vs_E_db_ done for $0

echo "#"
echo "# summarize results"
echo "#"
echo "cd results/"
cd results/

for BLAST_MODE in $BLAST_MODES; do
    echo "./extract_results.sh -m $BLAST_MODE -s bitscore"
    $SRUN ./extract_results.sh -m $BLAST_MODE -s bitscore
done

echo "#" 
echo "# announce we're done"
echo "#"
for BLAST_MODE in $BLAST_MODES; do
    echo -e "$BLAST_MODE\t$(sed 's/[:,]/\t/g' ./$BLAST_MODE/a/summary.bitscore.mismatch_totals.tsv)"
done


