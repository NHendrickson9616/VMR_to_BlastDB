#!/usr/bin/env python3
#
#
# ChatGPT prompt: please write python code that will load NCBI blasts
# output in XML format, using BioPython, and integrate through the
# HSPs for each hit, adding the query start to end intervals defined
# by each HSP to a hit-specific BioFrame, using the hit_def as the
# chromosome name, then use BioFrame::merge to compute that actual
# coverage from all merged HSPs for each hit.
#

from pprint import pprint

from Bio.Blast import NCBIXML

#from Bio import SearchIO
#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
#from Bio.SeqFeature import FeatureLocation
#from Bio.Graphics import GenomeDiagram
#from Bio.Graphics.GenomeDiagram import FeatureSet, GraphSet, Diagram

import bioframe as bf
import bioframe.vis

import itertools

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd


def load_blast_xml(xml_file):
    """
    Load NCBI BLAST XML output.
    """
    blast_records = NCBIXML.parse(open(xml_file, 'r'))
    return blast_records

def integrate_hsp_to_frame(blast_record):
    """
    Integrate HSPs for each hit and add query start to end intervals defined by each HSP to a hit-specific BioFrame.
    """
    print("integrate_hsp_to_frame()")
    frames = []
    for hit in blast_record.alignments:
        hit_def = hit.hit_def
        hit_parts = hit_def.split(" ")
        hit_name = hit_parts[0]
        for hsp in hit.hsps:
            frames.append([hit_name, hsp.query_start, hsp.query_end])

    #print(frames)

    return pd.DataFrame(frames,columns=['chrom','start','end'])

def compute_coverage(query_name, query_length, frames):
    print("# compute_coverage(frames)")
    """
    Compute actual coverage from all merged HSPs for each hit using BioFrame::merge.
    """
    hit_defs = set(frames["chrom"])
    #print("hit_defs:", hit_defs)
    scale = 100/query_length

    hit_stats = pd.DataFrame(columns=['pcov'])

    
    print("shape(frames):", frames.shape)
    for hit_def in hit_defs:
        print("# Hit Def:",hit_def, "######################################################################")

        ref_df = pd.DataFrame([[hit_def,1,query_length]],columns=['chrom','start','end'])
        #print("# ref_df = ",ref_df)
        
        hit_intervals = frames[frames["chrom"]==hit_def]
        print("# pre-merge ----------------------------------------------------------------------")
        #print(hit_intervals )
    
        scaled_hits = hit_intervals.copy()
        scaled_hits.start *= scale
        scaled_hits.end *= scale
        
        #print("scaled:", scaled_hits)
        bf.vis.plot_intervals(scaled_hits, show_coords=True, xlim=(0,100))
        plt.title('PRE-MERGE: ['+query_name+'] vs ['+hit_def+']');
        fname ='figures/'+hit_def+'.pre.pdf'
        plt.savefig(fname)
        print("wrote: ", fname)
        #plt.show()
    
        print("# post-merge ----------------------------------------------------------------------")
        hits_merged =  bf.merge(hit_intervals, min_dist=0)

        scaled_hits_merged = hits_merged.copy()
        scaled_hits_merged.start *= scale
        scaled_hits_merged.end *= scale

        bf.vis.plot_intervals(scaled_hits_merged, show_coords=True, xlim=(0,100))
        plt.title('POST-MERGE: ['+query_name+'] vs ['+hit_def+']');
        fname ='figures/'+hit_def+'.post.pdf'
        plt.savefig(fname)
        print("wrote: ", fname)
        #plt.show()

        print("# post-merge coverage---------------------------------------------------------------")
        hit_cov_df = bf.coverage(ref_df,hits_merged)
        print("# hit_cov_df: ", hit_cov_df)
        hit_cov = hit_cov_df.coverage[0]
        print("# hit_cov: ", hit_cov)
        hit_pcov = hit_cov/query_length
        print("# hit_pcov: ", hit_pcov)

        hit_stats.loc[hit_def] = hit_pcov
    print("hit_stats: ")
    pprint(hit_stats.sort_values(by='pcov',ascending=False))
# test
if 1==0:
    df1 = pd.DataFrame([
        ['chr1', 10000, 50000],
        ['chr1', 30000, 80000],
        ['chr1', 80000, 100000],
        ['chr1', 120000, 140000]],
                       columns=['chrom', 'start', 'end']
                       )

    scale = 100/max(df1.end)
    sdf1 = df1.copy()
    sdf1.start *= scale
    sdf1.end *= scale
    bf.vis.plot_intervals(sdf1, show_coords=True, xlim=(0,100))
    plt.title('bedFrame1 intervals');
    plt.show()

    df1m = bf.merge(df1, min_dist=0)
    print(df1m)
    sdf1m =df1m.copy()
    sdf1m.start *= scale
    sdf1m.end *= scale
    bf.vis.plot_intervals(sdf1m, show_coords=True, xlim=(0,100))
    plt.title('bedFrame1 MERGED intervals');
    plt.show()

# Example usage
#blast_xml_file = "your_blast_output.xml"
blast_xml_file = 'results/blastn10_test/a/Keyvirus/JX080302.5.xml'
#blast_xml_file = 'results/blastn10_test/a/Keyvirus/JX080302.5.maxtargetseqs10.maxhsps10.xml'
blast_xml_file = 'results/blastn10_test/a/Keyvirus/JX080302.5.hit10.xml'
blast_records = load_blast_xml(blast_xml_file)
for blast_record in blast_records:
    query_length=blast_record.query_length
    query_def = blast_record.query
    query_parts = query_def.split(" ")
    query_name = query_parts[0]
    
    print("#########################################################################")
    print("### RECORD: query-len=",blast_record.query_length, " query-ID:", blast_record.query_id," query-def:", blast_record.query)
    print("###")
    frames = integrate_hsp_to_frame(blast_record)
    #print("# FRAMES ----------------------------------------------------------------------")
    #pprint(frames)
    merged_frames = compute_coverage(query_name, query_length,frames)

