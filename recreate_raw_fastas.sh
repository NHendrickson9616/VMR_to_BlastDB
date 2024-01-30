#!/usr/bin/env bash
#
# un-process FASTAS to make them back into 
# raw NCBI download fastas 
# (remove species from front of fasta header)
#
# The .[0-9] part of the original accession is lost, so we replace it with ".?"
#
#SBATCH --job-name=ICTV_VMR_rebuild_raw_fastas_from_processed_fastas
#SBATCH --output=log.%J.%x.out
#SBATCH --error=log.%J.%x.out
#
# Number of tasks needed for this job. Generally, used with MPI jobs
# Time format = HH:MM:SS, DD-HH:MM:SS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=amd-hdr100 --time=00-12:00:00
##SBATCH --partition=amd-hdr100 --time=06-06:00:00
##SBATCH --partition=medium --time=40:00:00
#
# Number of CPUs allocated to each task. 
#
# Mimimum memory required per allocated  CPU  in  MegaBytes. 
#  last run was 402M
#SBATCH --mem-per-cpu=1000
#

FASTA_DIRS="fasta_new_vmr fasta_new_vmr_a"
if [ ! -z "$*" ]; then FASTA_DIRS="$*"; fi
echo FASTA_DIRS=$FASTA_DIRS

for FASTA_DIR in $FASTA_DIRS; do
    
    echo "# ----------------------------------------------------------------------"
    echo "#"
    echo "# FASTA_DIR=$FASTA_DIR"
    echo "#"
    echo "# ----------------------------------------------------------------------"


    echo "Total .fa's: " $(find $FASTA_DIR -name "*.fa" | wc -l)
    
    for FA in $(find $FASTA_DIR -name "*.fa"); do
	RAW=$(dirname $FA)/$(basename $FA .fa).raw
	echo "# ................................." 
	echo "# $FA" 

	echo "perl -pe 's/^>[^ ]+ ([^ ]+) />$1.? /' $FA > $RAW"
	perl -pe 's/^>[^ ]+ ([^ ]+) />$1.? /' $FA > $RAW

	head -1 $FA $RAW
    done
done

