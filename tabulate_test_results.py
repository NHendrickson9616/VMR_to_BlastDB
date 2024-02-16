#!/usr/bin/env python3
#
# tabulate QC Blast results of
#  A records vs BlastDB of E records
#
#
# future directions:
#   https://plumbum.readthedocs.io/en/latest/
#
import argparse
import pandas as pd
import re
import os
from pprint import pprint
from Bio import SeqIO

parser = argparse.ArgumentParser(description="")
parser.add_argument('-verbose',help="printout details during run",action=argparse.BooleanOptionalAction)
parser.add_argument('-blast',help="blast parameter abbreviation",default="blastn10")
args = parser.parse_args()

# hard coded things to move into ARGS
A_ACC_FILE="processed_accessions_a.tsv"
E_FASTA_FILE="fasta_new_vmr/vmr_e.fa"
E_FASTA_FILE="fasta_new_vmr/vmr_e.test4.fa"
BLAST_OUT_DIR="results"

#
# load list of A accessions
#
queryDf = pd.read_csv(A_ACC_FILE, sep="\t",
                           usecols=["Species","Accession_IDs","segment","Genus"],
                           dtype={'Species':'string','Accession_IDs':'string', 'segment':'string', 'Genus':'string'}
#                           ,index_col=["Accession_IDs"]
                           )
queryDf["acc_count"] = 1

print("loaded",A_ACC_FILE,": {0} rows, {1} columns.".format(*queryDf.shape))
print("genera: ",queryDf["Genus"].nunique())
print("species: ",queryDf["Species"].nunique())
print(queryDf.columns)
#print(queryDf.index)
#genera = queryDf.groupby('Genus').sum()
#print("Genera:")
##genera.shape
#print(genera.head(10))


if args.verbose:
#    for qacc in queryDf["Accession_IDs"].head(10):
#        print("Query:",qacc, queryDf.loc[qacc, "Species"])
    for qidx in range(0,len(queryDf.index)):
        print("Query:",qidx, queryDf["Accession_IDs"][qidx], queryDf["Species"][qidx],"[",queryDf["segment"][qidx],"]")

#
# load list of E accessions in the db
#
blastdb_headers = []
with open(E_FASTA_FILE, "r") as handle:
    # Iterate over each record in the file
    for record in SeqIO.parse(handle, "fasta"):
        # parse FASTA header
        # split VMR vs NCBI
        vmr_fields = re.split('[>#]',record.id)
        vmr_acc = vmr_fields[0]
        vmr_species = vmr_fields[1]
        vmr_seg = vmr_fields[2] if len(vmr_fields) > 2 else str(None)
        # split VMR vs NCBI
        (ncbi_acc,ncbi_header) = record.description.split(" ",1)[1].split(" ",1)
        if args.verbose and False: 
            #
            # check
            #
            print("ncbi_acc:\t",ncbi_acc)
            print("ncbi_header:\t",ncbi_header)
            print("vmr_acc:\t",vmr_acc)
            print("vmr_species:\t",vmr_species)
            print("vmr_seg:\t",vmr_seg)
        #
        # stash
        #
        blastdb_headers.append([vmr_acc, vmr_species, vmr_seg, ncbi_header])

#
# convert blastdb headers to panada DF
#
blastdb_header_columns=['vmr_acc', 'vmr_species', 'vmr_seg', 'ncbi_header']
blastdb_df = pd.DataFrame(blastdb_headers, columns=blastdb_header_columns, dtype=str)
blastdb_df.set_index(['vmr_species','vmr_seg'], inplace=True)
print(blastdb_df)


#
# iterate over query (A) records, parse blast output
#
def parse_blast_fmt7(filename):
    # Placeholder for hits data
    hits_data = []
    num_hits = 0
    no_hits = False
    blast_complete = False
    status = None

    # Open and read the BLAST results file
    if not os.path.exists(filename):
        status = "missing"
    else:
        with open(filename, 'r') as file:
            for line in file:
                # Skip comment lines but check for hits information
                if line.startswith('#'):
                    match = re.search(r"# (\d+) hits found", line)
                    if match:
                        num_hits = int(match.group(1))
                        if num_hits == 0:
                            no_hits = True
                    if '# BLAST processed 1 queries' in line:
                        blast_complete = True
                    continue

                # Parse tab-delimited hit lines
                hits_data.append(line.strip().split(','))
    
    # Convert hits data to DataFrame if there are hits
    if hits_data:
        df = pd.DataFrame(hits_data)
        print("\tloaded hits:", df.shape)
        df.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        num_hits = len(df)
    else:
        df = pd.DataFrame()  # Empty DataFrame if no hits

    # top level status
    #print("\tpre-status; nun_hits:", num_hits, " blast_complete:", blast_complete, ", status:", status)
    if status is None:
        #print("\tstatus is None")
        if not blast_complete:
            #print("\tstatus is trunc")
            status = "trunc"
        elif num_hits == 0:
            #print("\tstatus is no_hits")
            status = "no_hits"
        else:
            #print("\tstatus is ok")
            status = "OK"

    return status, num_hits, df

def analyze_blast_result(queryDf, qidx, blastdb_df, status, num_hits, hit_df):
    #
    # if blast didn't work, just note taht
    #
    if not status == "OK":
        queryDf["blast_status"][qidx] = status
        if args.verbose: print("\tblast_status = ", status)
        return

    #
    # figure out if correct result was in the db
    #
    qspecies=queryDf["Species"][qidx]
    qseg=queryDf["segment"][qidx]
    if args.verbose: print("\tdb_check: species=",qspecies, ", seg=", qseg)
    pprint(blastdb_df.loc[(qspecies,qseg)])
    

for qidx in range(0,len(queryDf.index)):
    qacc = queryDf["Accession_IDs"][qidx]
    qgenus = queryDf["Genus"][qidx]
    print("Query:",qidx, qgenus, qacc, queryDf["Species"][qidx],"[",queryDf["segment"][qidx],"]")
    
    #
    # read and parse blast results (may be empty!)
    #
    query_results_file = "./results/"+args.blast+"/a/"+qgenus+"/"+qacc+".raw.txt"
    print("\tReading "+query_results_file)
    status, num_hits, df_hits = parse_blast_fmt7(query_results_file)
    print("\tStatus: ", status," hits:", num_hits)

    #
    # compare the results
    #
    analyze_blast_result(queryDf, qidx, blastdb_df, status, num_hits, df_hits)
    exit(1)
