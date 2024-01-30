#!/usr/bin/env bash
#
# query each "A" (additional) sequence against the blastdb of "E" sequences. 
#
# split by genus
#
#
DIRS="fasta_new_vmr fasta_new_vmr_a"
for DIR in $DIRS; do
    echo DIR=$DIR

    GENERA=$(cd $DIR;ls -ls|grep drwx|sed 's/.* //g')
    for GENUS in $GENERA; do
	echo -n "   $DIR $GENUS "
	sbatch --job-name ICTV_VMR_rebuild_raw_fasta_${DIR}_${GENUS} ./recreate_raw_fastas.sh $DIR/$GENUS
    done
done
