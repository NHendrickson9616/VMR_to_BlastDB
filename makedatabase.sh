#!/bin/bash
# Make the BLAST database
makeblastdb -in ./fasta/all.fa -input_type "fasta" -title "ICTV refseqs" -out "./blast/ICTV_refseqs" -dbtype "nucl"

# Example usage:
# blastn -db ./blast/ICTV_refseqs -query ./fasta/U64512.fasta -out ./output/results.out