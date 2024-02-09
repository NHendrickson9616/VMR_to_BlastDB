#!/usr/bin/bash
#
# copy all data and temp files to LTS:/lefkowtiz-ictv-vmr
#
# IE, most everything that's NOT in git
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
echo rclone copy $* fasta_new_vmr   ${BUCKET}/fasta_new_vmr
rclone copy $* fasta_new_vmr   ${BUCKET}/fasta_new_vmr

cat <<EOF 
# ----------------------------------------------------------------------
#
# fasta downloads: A records
#
# ----------------------------------------------------------------------
EOF
echo rclone copy $* fasta_new_vmr_a ${BUCKET}/fasta_new_vmr_a
rclone copy $* fasta_new_vmr_a ${BUCKET}/fasta_new_vmr_a

cat <<EOF 
# ----------------------------------------------------------------------
# 
# blastdb
#
# ----------------------------------------------------------------------
EOF
echo rclone copy $* blast ${BUCKET}/blast
rclone copy $* blast ${BUCKET}/blast

cat <<EOF 
# ----------------------------------------------------------------------
#
# blast results
#
# ----------------------------------------------------------------------
EOF
echo rclone copy $* results ${BUCKET}/results
rclone copy $* results ${BUCKET}/results

cat <<EOF
# ----------------------------------------------------------------------
#
# this script itself
# (which is in git)
#
# ----------------------------------------------------------------------
EOF
echo rclone copy $* $0 ${BUCKET}/$0
rclone copy $* $0 ${BUCKET}/$0

