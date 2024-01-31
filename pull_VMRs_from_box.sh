#!/usr/bin/bash
#
#
module load rclone

rclone --verbose copy "box:/Virus Knowledgebase/VMR-blast/VMRs" ./VMRs/

