#!/bin/bash
query=$1
blastn -db ./blast/ICTV_refseqs -query $query -out ./output/results.out -outfmt '6 delim=,'

