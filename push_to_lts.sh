#!/usr/bin/bash
#
# copy all data and temp files to LTS:/lefkowtiz-ictv-vmr
#
# IE, most everything that's NOT in git
#

DEST=lts:/lefkowitz-ictv-vmr

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
echo rclone copy $* fasta_new_vmr   ${DEST}/fasta_new_vmr
rclone copy $* fasta_new_vmr   ${DEST}/fasta_new_vmr

cat <<EOF 
# ----------------------------------------------------------------------
#
# fasta downloads: A records
#
# ----------------------------------------------------------------------
EOF
echo rclone copy $* fasta_new_vmr_a ${DEST}/fasta_new_vmr_a
rclone copy $* fasta_new_vmr_a ${DEST}/fasta_new_vmr_a

cat <<EOF 
# ----------------------------------------------------------------------
# 
# blastdb
#
# ----------------------------------------------------------------------
EOF
echo rclone copy $* blast ${DEST}/blast
rclone copy $* blast ${DEST}/blast

cat <<EOF 
# ----------------------------------------------------------------------
#
# blast results
#
# ----------------------------------------------------------------------
EOF
echo rclone copy $* results ${DEST}/results
rclone copy $* results ${DEST}/results

cat <<EOF
# ----------------------------------------------------------------------
#
# this script itself
# (which is in git)
#
# ----------------------------------------------------------------------
EOF
echo rclone copy $* $0 ${DEST}/$0
rclone copy $* $0 ${DEST}/$0

