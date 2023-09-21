#!/usr/bin/env bash
#
# Create conda env 
#
conda create \
	-p ./conda/vmr \
	-c bioconda -c conda-forge \
	pandas openpyxl=2.5.7 numpy biopython blast xlrd

