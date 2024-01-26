#!/usr/bin/env bash
#
# Create conda env 
#
conda create \
	-p ./conda/vmr_openpyxl3_blast215 \
	-c bioconda -c conda-forge \
	pandas openpyxl=3 numpy biopython blast=2.15 xlrd

