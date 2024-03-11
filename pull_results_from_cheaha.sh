#!/usr/bin/env bash
#
rsync --progress -hav curtish@cheaha.rc.uab.edu:/scratch/curtish/ictv/VMR/VMR_to_BlastDB/processed_accessions_a.tsv .
rsync --progress -hav curtish@cheaha.rc.uab.edu:/scratch/curtish/ictv/VMR/VMR_to_BlastDB/processed_accessions_e.tsv .
rsync --progress -hav curtish@cheaha.rc.uab.edu:/scratch/curtish/ictv/VMR/VMR_to_BlastDB/fasta_new_vmr/vmr_e.fa fasta_new_vmr
rsync --progress -hav curtish@cheaha.rc.uab.edu:/scratch/curtish/ictv/VMR/VMR_to_BlastDB/fasta_new_vmr_a ./
rsync --progress -hav curtish@cheaha.rc.uab.edu:/scratch/curtish/ictv/VMR/VMR_to_BlastDB/results/blastn10 results
rsync --progress -hav curtish@cheaha.rc.uab.edu:/scratch/curtish/ictv/VMR/VMR_to_BlastDB/results/blastn7  results
