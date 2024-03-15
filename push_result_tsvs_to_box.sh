#!/usr/bin/env bash
#
VMR="VMR_MSL38_v3_FINAL.2024-03-14"
BLAST=blastn10
#
if [ -z "$(which rclone 2>/dev/null)" ]; then
    echo module load rclone
    module load rclone
fi

#
# VMR parsing results
# 
SRC_URL="."
DEST_URL="box:/Virus Knowledgebase/VMR-blast/reports/$VMR/"

for FILE in vmr.tsv processed_accessions_e.tsv processed_accessions_a.tsv tabulation.${BLAST}.tsv; do 
     echo rclone --verbose copy "$SRC_URL" "$DEST_URL" --filter "- conda/" --filter "- results/" --filter "+ ${FILE}" --filter "- **"
          rclone --verbose copy "$SRC_URL" "$DEST_URL" --filter "- conda/" --filter "- results/" --filter "+ ${FILE}" --filter "- **"
done

#
# blast results
# 
SRC_URL="./results/${BLAST}/a"
DEST_URL="box:/Virus Knowledgebase/VMR-blast/reports/${VMR}/${BLAST}/"

echo rclone --verbose copy "$SRC_URL" "$DEST_URL" --filter "- summary.*" --filter "+ *.tsv" 
     rclone --verbose copy "$SRC_URL" "$DEST_URL" --filter "- summary.*" --filter "+ *.tsv" 

