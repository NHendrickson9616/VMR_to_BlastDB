#!/bin/bash

echo "# concatenate individual fasta's into all.fa"
# if any fastas are updated, rebuild master fasta
if [ "$(find fasta_new_vmr -newer fasta_new_vmr/vmr_e.fa|wc -l)" -gt 0 ]; then
    echo 'cat ./fasta_new_vmr/*/*.fa > ./fasta_new_vmr/vmr_e.fa'
    cat ./fasta_new_vmr/*/*.fa > ./fasta_new_vmr/vmr_e.fa
else
    echo "SKIP: ./fasta_new_vmr/vmr_e.fa is up-to-date."
fi


echo "# Make the BLAST database"
if [ "$(which makeblastdb 2>/dev/null)" == "" ]; then 
    echo "module load BLAST"
    module load BLAST
fi

echo 'makeblastdb -in ./fasta_new_vmr/vmr_e.fa -input_type "fasta" -title "ICTV VMR refseqs" -out "./blast/ICTV_VMR_e" -dbtype "nucl"'
makeblastdb -in ./fasta_new_vmr/vmr_e.fa -input_type "fasta" -title "ICTV VMR refseqs" -out "./blast/ICTV_VMR_e" -dbtype "nucl"

echo "# Example usage:"
echo "# blastn -db ./blast/ICTV_VMR_e -query ./fasta_new_vmr/Eponavirus/MG711462.fa  -out ./results/e/Eponavirus/MG711462.out"
