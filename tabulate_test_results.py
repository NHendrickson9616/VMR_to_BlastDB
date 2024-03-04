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

# blast output parsing & HSP collapsing
from Bio.Blast import NCBIXML
import bioframe as bf
import bioframe.vis

import itertools

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#
# parse ARGS
#
parser = argparse.ArgumentParser(description="")
parser.add_argument('-verbose',help="printout details during run",action=argparse.BooleanOptionalAction)
parser.add_argument('-blast',help="blast parameter abbreviation",default="blastn10")
args = parser.parse_args()

# hard coded things to move into ARGS
A_ACC_FILE="processed_accessions_a.tsv"
E_FASTA_FILE="fasta_new_vmr/vmr_e.fa"
#E_FASTA_FILE="fasta_new_vmr/vmr_e.test4.fa"
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
#queryDf = queryDf.head(4)
print("loaded",A_ACC_FILE,": {0} rows, {1} columns.".format(*queryDf.shape))
print("genera: ",queryDf["Genus"].nunique())
print("species: ",queryDf["Species"].nunique())
print("column names: ",queryDf.columns)
#print(queryDf.index)
#genera = queryDf.groupby('Genus').sum()
#print("Genera:")
##genera.shape
#print(genera.head(10))

# replace NA segments with ""
queryDf.fillna({"segment": ""}, inplace=True)


#
# print list of query (A) sequences
#
#if args.verbose:
#    for qacc in queryDf["Accession_IDs"].head(10):
#        print("Query:",qacc, queryDf.loc[qacc, "Species"])
#    for qidx in queryDf.index:
#        print("Query:",qidx, queryDf.at[qidx,"Accession_IDs"], queryDf.at[qidx,"Species"],"[",queryDf.at[qidx,"segment"],"]")

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
        vmr_seg = vmr_fields[2] if len(vmr_fields) > 2 else ""
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
# stats
#
print("loaded:",E_FASTA_FILE,": {0} rows.".format(len(blastdb_headers)))
#
# convert blastdb headers to panada DF
#
blastdb_header_columns=['vmr_acc', 'vmr_species', 'vmr_seg', 'ncbi_header']
blastdb_df = pd.DataFrame(blastdb_headers, columns=blastdb_header_columns, dtype=str)
blastdb_df.set_index(['vmr_species','vmr_seg'], inplace=True)

print(blastdb_df)

#
# get query sequence length
#
def get_fasta_length(filename):
    with open(filename, "r") as handle:
        # Iterate over each record in the file
        for record in SeqIO.parse(handle, "fasta"):
            sequence_length=len(record.seq)
            return sequence_length
    return None

#
# parse XML blast results
#
# build per-query-subject coverage score by overlapping all HSPs
#

def blast_compute_hsp_coverage(blast_query):
    #
    # get query info
    #
    query_length = blast_query.query_length
    query_def    = blast_query.query
    # parse our headers
    query_parts  = query_def.split(" ")
    query_name   = query_parts[0]

    #
    # convert HSPs intervals into BioFrame
    #
    query_frames = []
    subject_frames = []
    for subject in blast_query.alignments:
        subject_def = subject.hit_def
        # parse our fasta headers
        subject_parts = subject_def.split(" ")
        subject_name = subject_parts[0]
        for hsp in subject.hsps:
            query_frames.append([subject_name, hsp.query_start, hsp.query_end])
            #subject_frames.append([subject_name, hsp.subject_start, hsp.subject_end])
    query_hsp_df = pd.DataFrame(query_frames,columns=['chrom','start','end'])
    #subject_hsp_df = pd.DataFrame(subject_frames,columns=['chrom','start','end'])


    print("QUERY_FRAMES:", query_hsp_df)#
    # compute query coverage per subject
    #
    subject_stats = pd.DataFrame(columns=['pcov'])

    scale = 100/query_length

    subject_defs = set(query_hsp_df["chrom"])
    
    for subject_def in subject_defs:
        if args.verbose: 
            print("# Hit Def:",subject_def, "######################################################################")

        # get HSPs for this subject
        subject_intervals = query_hsp_df[query_hsp_df["chrom"]==subject_def]

        # scale query coverage from 0-1
        ref_df = pd.DataFrame([[subject_def,1,query_length]],columns=['chrom','start','end'])
        scaled_hits = subject_intervals.copy()
        scaled_hits.start *= scale
        scaled_hits.end *= scale

        print("# pre-merge ----------------------------------------------------------------------")
        bf.vis.plot_intervals(scaled_hits, show_coords=True, xlim=(0,100))
        plt.title('PRE-MERGE: ['+query_name+'] vs ['+subject_def+']');
        fname ='figures/'+query_name+'.vs.'+subject_def+'.pre.pdf'
        plt.savefig(fname)
        print("wrote: ", fname)
        #plt.show()

        print("# post-merge ----------------------------------------------------------------------")
        hits_merged =  bf.merge(subject_intervals, min_dist=0)

        scaled_hits_merged = hits_merged.copy()
        scaled_hits_merged.start *= scale
        scaled_hits_merged.end *= scale

        bf.vis.plot_intervals(scaled_hits_merged, show_coords=True, xlim=(0,100))
        plt.title('POST-MERGE: ['+query_name+'] vs ['+subject_def+']');
        fname ='figures/'+query_name+'.vs.'+subject_def+'.post.pdf'
        plt.savefig(fname)
        print("wrote: ", fname)
        #plt.show()

        print("# post-merge coverage---------------------------------------------------------------")
        subject_cov_df = bf.coverage(ref_df,hits_merged)
        print("# subject_cov_df: ", subject_cov_df)
        subject_cov = subject_cov_df.coverage[0]
        print("# subject_cov: ", subject_cov)
        subject_pcov = subject_cov/query_length
        print("# subject_pcov: ", subject_pcov)

        subject_stats.loc[subject_def] = subject_pcov

        
    return subject_stats.sort_values(by=['pcov'],ascending=False)


    
def parse_blast_xml_collapse_hsps(filename):
    #
    # load XML
    #
    blast_records = NCBIXML.parse(open(filename, 'r'))
    if args.verbose:
        print("#########################################################################")
        print("## Loaded: ",filename)

    #
    # iterate over queries (should be 1)
    #
    for blast_query in blast_records:
        query_length=blast_query.query_length
        query_def = blast_query.query
        # parse our headers
        query_parts = query_def.split(" ")
        query_name = query_parts[0]
        if args.verbose:
            print("## query-len=",blast_record.query_length,
                  " query-ID:", blast_record.query_id,
                  " query-def:", blast_record.query)

        subject_stats = blast_compute_hsp_coverage(blast_query)

        print(subject_stats)
        
    return subject_stats
    
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
        col_names=    ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        hits_df = pd.DataFrame(hits_data,columns=col_names)#,dtype={col: int for col in range(2,len(col_names))})
        # convert to numeric
        for col in col_names[2:]:
            hits_df[col] = pd.to_numeric(hits_df[col])
        print("\tloaded hits:", hits_df.shape)
        #pprint(hits_df.info())
        num_hits = len(hits_df)
    else:
        hits_df = pd.DataFrame()  # Empty DataFrame if no hits

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

    #
    # unpack sseqid into species & seg
    #
    if len(hits_df.index)>0:
        split_values = hits_df['sseqid'].str.split('#', expand=True)
        hits_df['saccession'] = split_values[0]
        hits_df['sspecies'] = split_values[1]
        # if there are no sseqid's with a segment, then fill that column with empty
        if split_values.shape[1] > 2:
            hits_df['sseg'] = split_values[2].fillna("")
        else:
            hits_df['sseg'] = ""

    return status, num_hits, hits_df

def analyze_blast_result(queryDf, qidx, blastdb_df, status, num_hits, hits_df):
    if args.verbose: print("# analyze_blast_result(",queryDf.at[qidx,"Accession_IDs"],", status=", status, ", num_hits=", num_hits,"/",len(hits_df.index),")")
    if args.verbose: print("\tQuery Len: {qlen}".format(qlen=queryDf.at[qidx,"qlen"]))

    #
    # if blast didn't work, just note taht
    #
    if not status == "OK":
        queryDf.at[qidx,"blast_status"] = status
        if args.verbose: print("\tblast_status = ", status)
        return

    #
    # figure out if correct result was in the db
    #
    qspecies=queryDf.at[qidx,"Species"]
    qspecies_=qspecies.replace(" ","_")
    qseg=queryDf.at[qidx,"segment"]
    if args.verbose: print("\tdb_check: species='{qspecies}', seg='{qseg}'".format(qspecies=qspecies_, qseg=qseg))
    species_matches_df = blastdb_df.query("vmr_species==@qspecies_")
    if args.verbose: print("\tE species records: ", len(species_matches_df.index))
    # check this species exists in the db
    if len(species_matches_df.index) < 1: 
        status = "E_NO_SPECIES"
    else:
        species_seg_matches_df = species_matches_df.query("vmr_seg==@qseg")
        if args.verbose: print("\tE species/seg records: ", len(species_seg_matches_df.index))
        # check this segment exists in the db
        if len(species_seg_matches_df.index) < 1:
            status = "E_SPECIES_NO_SEG"
        elif len(hits_df.index) < 1:
            status = "BLAST_NO_HITS"
        else:
            print("\tchecking {0} blast results...".format(len(hits_df.index)))

            #
            # filter to 60% of query length
            #
            min_length = queryDf.at[qidx,"qlen"] * 0.60
            filt_hits_df = hits_df.query("length > @min_length")
            print("\tfilter for length > {min_length} leaves {hits} blast results...".format(min_length=min_length,hits=len(filt_hits_df.index)))
            
            
            #
            # ok, actually analyze the blast results
            #
            # columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
            #         ['saccession', 'sspecies', 'sseg'
            match_species_df = filt_hits_df.query("sspecies==@qspecies_")
            print("\tBLAST species hits: ", len(match_species_df.index))
            match_species_seg_df = match_species_df.query("sseg==@qseg")
            print("\tBLAST species-seg hits: ", len(match_species_seg_df.index))
            
    
    return status

for qidx in queryDf.index[0:10]:
    qacc = queryDf.at[qidx,"Accession_IDs"]
    qgenus = queryDf.at[qidx,"Genus"]
    print("Query:",qidx, qgenus, qacc, queryDf.at[qidx,"Species"],"[",queryDf.at[qidx,"segment"],"]")
    
    #
    # get query sequence length
    #
    query_fasta_file = "./fasta_new_vmr_a/{genus}/{accession}.fa".format(genus=qgenus, accession=qacc)
    qlen = get_fasta_length(query_fasta_file)
    queryDf.at[qidx,"qlen"] = int(qlen)
    if args.verbose: print("\tQuery Len: {qlen} from {qfasta}".format(qlen=queryDf.at[qidx,"qlen"], qfasta=query_fasta_file))
    
    #
    # read and parse blast results (may be empty!)
    #
    query_results_file = "./results/"+args.blast+"/a/"+qgenus+"/"+qacc+"..hit20.xml"
    #query_results_file = "./results/"+args.blast+"/a/"+qgenus+"/"+qacc+".raw.txt"
    print("\tReading "+query_results_file)
    hits_df = parse_blast_xml_collapse_hsps(query_results_file)
    print("\thits:", len(hits_df.index))

    #
    # compare the results
    #
    status = analyze_blast_result(queryDf,qidx, blastdb_df, status, num_hits, hits_df)

    print("\tHIT: ", status)
    #exit(1)
