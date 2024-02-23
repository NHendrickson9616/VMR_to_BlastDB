#!/usr/bin/env zsh
#
# pull data/results files from LTS
#
# these files are big and not stored in git
#

BUCKET=lts:/lefkowitz-ictv-vmr

if [ -z "$(which rclone 2>/dev/null)" ]; then 
    echo module load rclone
    module load rclone
fi
cat <<EOF 
# ----------------------------------------------------------------------
#
# fasta downloads: E records
#
# ----------------------------------------------------------------------
EOF
DIR=fasta_new_vmr
echo rclone copy $* ${BUCKET}/${DIR} ${DIR}
     rclone copy $* ${BUCKET}/${DIR} ${DIR}
cat <<EOF 
# ----------------------------------------------------------------------
#
# fasta downloads: A records
#
# ----------------------------------------------------------------------
EOF
DIR=fasta_new_vmr_a
echo rclone copy $* ${BUCKET}/${DIR} ${DIR}
     rclone copy $* ${BUCKET}/${DIR} ${DIR}
cat <<EOF 
# ----------------------------------------------------------------------
# 
# blastdb
#
# ----------------------------------------------------------------------
EOF
DIR=blast
echo rclone copy $* ${BUCKET}/${DIR} ${DIR}
     rclone copy $* ${BUCKET}/${DIR} ${DIR}
cat <<EOF 
# ----------------------------------------------------------------------
#
# blast results
#
# ----------------------------------------------------------------------
EOF
DIR=results
echo rclone copy $* ${BUCKET}/${DIR} ${DIR}
     rclone copy $* ${BUCKET}/${DIR} ${DIR}
cat <<EOF
# ----------------------------------------------------------------------
# DONE
# ----------------------------------------------------------------------
EOF
