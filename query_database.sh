#!/bin/bash
query=$1
results_file=$2
results_dir=$(basename "$results_file")
mkdir -p "${results_dir}"

if [ -z "$(which blastn 2>/dev/null)" ]; then
    echo "#module load BLAST"
    module load BLAST
fi

echo "blastn -db ./blast/ICTV_VMR_e -query $query -out ${results_file} -outfmt '7 delim=,'"
      blastn -db ./blast/ICTV_VMR_e -query $query -out ${results_file} -outfmt '7 delim=,'

