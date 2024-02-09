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
DIR=fasta_new_vmr
echo rclone copy $* ${DIR} ${BUCKET}/${DIR}
     rclone copy $* ${DIR} ${BUCKET}/${DIR}

cat <<EOF 
# ----------------------------------------------------------------------
#
# fasta downloads: A records
#
# ----------------------------------------------------------------------
EOF
DIR=fasta_new_vmr_a
echo rclone copy $* ${DIR} ${BUCKET}/${DIR}
     rclone copy $* ${DIR} ${BUCKET}/${DIR}

cat <<EOF 
# ----------------------------------------------------------------------
# 
# blastdb
#
# ----------------------------------------------------------------------
EOF
DIR=blast
echo rclone copy $* ${DIR} ${BUCKET}/${DIR}
     rclone copy $* ${DIR} ${BUCKET}/${DIR}

cat <<EOF 
# ----------------------------------------------------------------------
#
# blast results
#
# ----------------------------------------------------------------------
EOF
DIR=results
echo rclone copy $* ${DIR} ${BUCKET}/${DIR}
     rclone copy $* ${DIR} ${BUCKET}/${DIR}

cat <<EOF
# ----------------------------------------------------------------------
#
# this script itself
# (which is in git)
#
# ----------------------------------------------------------------------
EOF
DIR=$0
echo rclone copy $* ${DIR} ${BUCKET}/${DIR}
     rclone copy $* ${DIR} ${BUCKET}/${DIR}


cat <<EOF
# ----------------------------------------------------------------------
# DONE
# ----------------------------------------------------------------------
EOF
