import pandas
import subprocess
from Bio import Entrez
import time
from urllib import error
import argparse
def fetch_fasta():
    bad_accessions = []
    Accessions = pandas.read_excel('processed_accessions.xlsx')
    Accessions = Accessions['Accession_IDs'].values[1:3]
    all_reads = []
    # Fetchs FASTA data for every accession number
    for accession_ID in Accessions:
        Entrez.email = 'halko@uab.edu'
        try:
            handle = Entrez.efetch(db="nuccore", id=accession_ID, rettype="xml", retmode="text")
            read = handle.read()
            read = read.split('taxname')[1]
            read = read.split(',')[0]
            
            file = open('return_xml.txt','w')
            file.writelines(str(read))
            file.close()
        # prints out accession that got though cleaning
        except(error.HTTPError):
            print(accession_ID)
            bad_accessions = bad_accessions+[accession_ID]
        
        

        
fetch_fasta()