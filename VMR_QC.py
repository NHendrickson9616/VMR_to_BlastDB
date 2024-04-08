#!/usr/bin/env python3
import pandas
import subprocess
import time
from urllib import error
import argparse
import numpy as np
import sys
import os
import pathlib # for stem=basename(.txt)

# Class needed to load args from files. 
class LoadFromFile (argparse.Action):
    def __call__ (self, parser, namespace, values, option_string = None):
        with values as f:
            # parse arguments in the file and store them in the target namespace
            parser.parse_args(f.read().split(), namespace)
parser = argparse.ArgumentParser(description="")

#setting arguments.
parser.add_argument('-verbose',help="printout details during run",action=argparse.BooleanOptionalAction) #sets args.verbose to true
parser.add_argument('-file',help="optional argument. Name of the file to get arguments from.",type=open, action=LoadFromFile)
parser.add_argument("-email",help="email for Entrez to use when fetching Fasta files")
parser.add_argument("-mode",help="what function to do. Options: VMR,fasta,db")
parser.add_argument("-ea",help="Fetch E or A records (Exemplars or AdditionalIsolates)", default="E")
parser.add_argument("-VMR_file_name",help="name of the VMR file to load.",default="VMR_E_data.xlsx")
parser.add_argument("-query",help="Name of the fasta file to query the database")
args = parser.parse_args()
if args.mode != 'fasta' and args.mode != "VMR" and args.mode != "db":
    print("Valid mode not selected. Options: VMR,fasta,db",file=sys.stderr)
#Takes forever to import so only imports if it's going to be needed
if args.mode == 'fasta':
    from Bio import Entrez
#Catching error
if args.mode == "db":
    if args.query == None:
        print("Database Query mode is selected but no fasta file was specified! Please set the '-fasta_file_name' or change mode.",file=sys.stderr)

VMR_file_name_tsv = './vmr.tsv'
VMR_hack_file_name = "./fixed_vmr_"+args.ea.lower()+".tsv"
processed_accession_file_name ="./processed_accessions_"+args.ea.lower()+".tsv"


###############################################################################################################
# Loads excel from https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/ and puts it into a DataFrame
# NOTE: URL is incorrect. 
############################################################################################################### 
# DataFrame['column name'] = provides entire column
# DataFrame['column name'][0,1,2,3,4,5 etc] provides row for that column
# 
#
def load_VMR_data():
    if args.verbose: print("load_VMR_data()")
    if args.verbose: print("  opening", args.VMR_file_name)

    # Importing excel sheet as a DataFrame. Requires xlrd and openpyxl package
    try:
        raw_vmr_data = pandas.read_excel(args.VMR_file_name,engine='openpyxl')
        if args.verbose: print("VMR data loaded: {0} rows, {1} columns.".format(*raw_vmr_data.shape))

        # list of the columns to extract from raw_vmr_data
        vmr_cols_needed = ['Exemplar or additional isolate','Sort','Isolate Sort','Virus GENBANK accession','Species','Genome coverage','Genus']

        for col_name in list(raw_vmr_data.columns):
            if col_name in vmr_cols_needed:
                print("    "+col_name+" [NEEDED]")
            else:
                print("    "+col_name)
        for col_name in vmr_cols_needed:
            if not col_name in list(raw_vmr_data.columns):
                print("    "+col_name+" [!MISSING!]")

    except(FileNotFoundError):
        print("The VMR file specified does not exist! Make sure the path set by '-VMR_file_name' is correct.",file=sys.stderr)
    

    # save As TSV for diff'ing
    if os.path.exists(VMR_file_name_tsv) and os.path.getmtime(VMR_file_name_tsv) > os.path.getmtime(args.VMR_file_name):
        if args.verbose: print("  SKIP writing", VMR_file_name_tsv)
    else:
        if args.verbose: print("  writing", VMR_file_name_tsv)
        raw_vmr_data.to_csv(VMR_file_name_tsv,sep='\t', index=False)

    # compiling new dataframe from vmr_cols_needed
    #truncated_vmr_data = raw_vmr_data[vmr_cols_needed]

    # DataFrame.loc is helpful for indexing by row. Allows expression as an argument. Here, 
    # it finds every row where 'E' is in column 'Exemplar or additional isolate' and returns 
    # only the columns specified. 
    #vmr_data = truncated_vmr_data.loc[truncated_vmr_data['Exemplar or additional isolate']==args.ea.upper(),['Sort','Isolate Sort','Species','Virus GENBANK accession',"Genome coverage","Genus"]]
    vmr_data = raw_vmr_data.loc[raw_vmr_data['Exemplar or additional isolate']==args.ea.upper(),vmr_cols_needed[1:]]

    # only works when I reload the vmr_data, probably not necessary. have to look into why it's doing this. 
    if args.verbose: print("Writing"+VMR_hack_file_name,": workaround - filters VMR down to "+args.ea.upper()+" records only")
    if args.verbose: print("\tcolumns: ",vmr_data.columns)
    
    vmr_data.to_csv(VMR_hack_file_name, sep='\t')
    if args.verbose: print("Loading",VMR_hack_file_name)
    narrow_vmr_data = pandas.read_csv(VMR_hack_file_name,sep='\t')
    if args.verbose: print("   Read {0} rows, {1} columns.".format(*narrow_vmr_data.shape))
    if args.verbose: print("   columns:", list(narrow_vmr_data.columns))

    # Removing Genome Coverage column from the returned value. 
    truncated_vmr_data = narrow_vmr_data[vmr_cols_needed[1:]]
    
    #truncated_vmr_data = truncated_vmr_data.drop(columns=['Exemplar or additional isolate'])
    if args.verbose: print("   Truncated: {0} rows, {1} columns.".format(*truncated_vmr_data.shape))
    
    return truncated_vmr_data


def test_accession_IDs(df):
    if args.verbose: print("test_accession_IDs()")
    if args.verbose: print("\tcolumns: ",df.columns)
##############################################################################################################
# Cleans Accession numbers assuming the following about the accession numbers:
# 1. Each Accession Number is 6-8 characters long
# 2. Each Accession Number contains at least 3 numbers
# 3. Each Accession Number contains at most 3 letters
# 4. Accession Numbers in the same block are seperated by a ; or a , or a :
##############################################################################################################
    # defining new DataFrame before hand
    processed_accession_IDs = pandas.DataFrame(columns=['Species','Accession_IDs','segment',"Genus","Sort","Isolate Sort","Original_Accession_String","Errors"])
    # for loop for every entry in given processed_accessionIDs
    for entry_count in range(0,len(df.index)):
        # Find current species in this row
        Species=df['Species'][entry_count]
        # Find current Genus in this row
        Genus=df['Genus'][entry_count]
        #initial processing of raw accession numbers.
        orig_accession_str = df['Virus GENBANK accession'][entry_count]
        #removing whitespace.
        accession_ID = str(orig_accession_str).replace(" ","")
        # instead of trying to split by commas and semicolons, I just replace the commas with semicolons. 
        accession_ID.replace(",",";")
        accession_ID = accession_ID.split(';')
        
        #split by colons too
        accession_ID = [accession_part.split(':') for accession_part in accession_ID]
    
        # for loop for every ";" split done
        for accession_seg in accession_ID:
            segment = None
            errors= ""
            #flag long strings as being suspect. Most aren't 10 characters long. 
            if len(accession_seg) > 10:
                errors = 'suspiciously long accession number or segment name. Please verify its correct:'+str(accession_seg)
                print(errors,file=sys.stderr)

                # for loop for every ":" split done
            for accession_part in accession_seg:
                
                number_count = 0
                letter_count = 0
                # counting letters
                for char in accession_part:
                    if char in 'qwertyuiopasdfghjklzxcvbnm':
                        letter_count = letter_count+1
                # counting numbers
                    elif char in '1234567890':
                        number_count = number_count+1
                #checks if current selection fits what an accession number should be
                if len(str(accession_part)) == 8 or 6 and letter_count<3 and number_count>3:
                    processed_accession_ID = accession_part
                    processed_accession_IDs.loc[len(processed_accession_IDs.index)] = [Species,processed_accession_ID,segment,Genus,
                                                                                       df['Sort'][entry_count],
                                                                                       df['Isolate Sort'][entry_count],
                                                                                       df['Virus GENBANK accession'][entry_count],
                                                                                       errors ]
                    #print("'"+processed_accession_ID+"'"+' has been cleaned.')
                else:
                    segment = accession_part
    return processed_accession_IDs               


#######################################################################################################################################
# Utilizes Biopython's Entrez API to fetch FASTA data from Accession numbers. 
# Prints Accession Numbers that failed to 'clean' correctly
# 
# this should use epost to work in batches
#######################################################################################################################################  
def fetch_fasta(processed_accession_file_name):
    if args.verbose: print("fetch_fasta(",processed_accession_file_name,")")

    # make sure the output directory exists
    fasta_dir = "./fasta_new_vmr"
    if args.ea == "a": 
        fasta_dir = fasta_dir+"_a"
    if not os.path.exists(fasta_dir):
        # Create the directory if it doesn't exist
        os.makedirs(fasta_dir)
        if args.verbose: print(f"Directory '{fasta_dir}' created successfully.")

    bad_accessions_fname="./bad_accessions_"+args.ea.lower()+".tsv"

    #Check to see if fasta data exists and, if it does, loads the accessions numbers from it into an np array.
    if args.verbose: print("  loading:", processed_accession_file_name)
    Accessions = pandas.read_csv(processed_accession_file_name,sep='\t')

    all_reads = []
    bad_accessions = pandas.DataFrame(columns=Accessions.columns)

    # NCBI Entrez Session setup
    entrez_sleep = 0.34 # 3 requests per second with email authN
    Entrez.email = args.email
    if "NCBI_API_KEY" in os.environ:
        # use API_KEY  authN (10 queries per second)
        entrez_sleep = 0.1 # 10 requrests per second with API_KEY
        Entrez.api_key = os.environ["NCBI_API_KEY"]
        if args.verbose: print("NCBI Entrez 10/second with NCBI_API_KEY")
    else: 
        # use email authN
        if args.verbose: print("NCBI Entrez 3/second with email=",args.email)

    # Fetches FASTA data for every accession number
    count = 0
    for accession_ID in Accessions['Accession_IDs']:
            row = Accessions.loc[count]
            Species = row.iloc[1]
            accession_ID = row.iloc[2]
            segment = row.iloc[3]
            genus_name = row.iloc[4]

            # emtpy cell becomes float:NaN!
            if segment != segment: 
                segment = ""
            if genus_name != genus_name: 
                genus_name = ""
                
            # fasta_file_name
            genus_dir = fasta_dir+"/"+genus_name
            if genus_name == "":
                genus_dir = fasta_dir+"/"+"no_genus"
            accession_raw_file_name = genus_dir+"/"+str(accession_ID)+".raw"
            accession_fa_file_name = genus_dir+"/"+str(accession_ID)+".fa"
            
            # make sure dir exists
            if not os.path.exists(genus_dir):
                # Create the directory if it doesn't exist
                os.makedirs(genus_dir)
                if args.verbose: print(f"Directory '{genus_dir}' created successfully.")
    
            # check if the raw file exists
            if os.path.exists(accession_raw_file_name):
                if args.verbose: print("[FETCH]  SKIP NCBI fetch for {accession_raw_file_name}".format(**locals()))
            else:
                raw_file = open(accession_raw_file_name,'w')
                if args.verbose: print("[FETCH]  EXEC NCBI fetch for {accession_raw_file_name}".format(**locals()))
                try:
                    # fetch FASTA from NCBI
                    handle = Entrez.efetch(db="nuccore", id=accession_ID, rettype="fasta", retmode="text")

                    # limit requests: 3/second with email, 10/second with API_KEY
                    time.sleep(entrez_sleep)

                    # prints out accession that got though cleaning
                    if args.verbose: print('    fasta for '+accession_ID+ ' obtained.')

                    # prints out accession that got though cleaning
                    raw_fa = handle.read()
                    raw_file.write(raw_fa);
                    raw_file.close()
                    if args.verbose: print('    wrote: '+accession_raw_file_name)

                except:
                    print("    [ERR] Accession ID"+"'"+str(accession_ID)+"'"+"did not get properly cleaned. Accession Cleaning Heuristic needs editing.",file=sys.stderr)
                    bad_accessions.append(row)

            # check if processed fasta is out of date
            if os.path.getsize(accession_raw_file_name) == 0:
                if args.verbose: print("[FORMAT] SKIP/ERROR raw files is empty for {accession_fa_file_name}".format(**locals()))
            elif os.path.exists(accession_fa_file_name) and os.path.getmtime(accession_fa_file_name) > os.path.getmtime(accession_raw_file_name):
                if args.verbose: print("[FORMAT] SKIP reformat header for {accession_fa_file_name}".format(**locals()))
            else:
                if args.verbose: print("[FORMAT] EXEC reformat header for {accession_fa_file_name}".format(**locals()))
                
                # open local raw genbank fasta
                raw_file = open(accession_raw_file_name,'r')
                raw_fa = raw_file.read()
                raw_file.close()

                # open local (header modified) version
                fa_file  = open(accession_fa_file_name,'w')

                # parse out header and seq
                fa_desc = raw_fa.split("\n")[0].replace(">","")
                ncbi_accession = fa_desc.split(" ",1)[0]
                fa_seq =    raw_fa.split("\n",1)[1]

                # build ICTV-modified header
                #  ACCESSION#VMR_SPECIES[#VMR_SEG] NCBI_HEADER
                if str(segment).lower() == "":
                    desc_line = ">"+str(ncbi_accession)+"#"+str(Species.replace(" ","_"))                 +" "+fa_desc
                else:
                    desc_line = ">"+str(ncbi_accession)+"#"+str(Species.replace(" ","_"))+"#"+str(segment)+" "+fa_desc
                if args.verbose: print("    ", desc_line)

                # write ICTV formated header to fasta
                fa_file.write(desc_line+"\n"+fa_seq)
                fa_file.close()
                if args.verbose: print('    wrote: '+accession_fa_file_name)
                
            count=count+1

    # wrap up and report errors
    print("Bad_Accession count:", len(bad_accessions.index))
    pandas.DataFrame.to_csv(bad_accessions,bad_accessions_fname,sep='\t')
    print("Wrote to ", bad_accessions_fname)
    
#######################################################################################################################################
# Calls makedatabase.sh. Uses 'all.fa'
#######################################################################################################################################
def make_database():
    if args.verbose: print("make_database()")
    p1 = subprocess.run("makedatabase.sh")
    
#######################################################################################################################################
# BLAST searches a given FASTA file and returns DataFrame rows with from the accession numbers. Returns in order of significance.  
#######################################################################################################################################
# How closely related
# Run 'A' viruses
# Compare to members of same species
# BLAST score -- 
# Top hit
# check to see if many seg return same virus

def query_database(path_to_query):
    results_dir="./results/e"
    results_file=results_dir+"/"+pathlib.Path(path_to_query).stem+".csv"

    if args.verbose: print("query_database("+path_to_query+")")
    if args.verbose: print("   run(query_database.sh "+path_to_query+" "+results_file)
    p1 = subprocess.run(["bash","query_database.sh",path_to_query,results_file])
    """
    if args.verbose: print("   reading: "+results_file")
    results = open(results_file,"r")
    result_text = results.readlines()
    results.close()
    print(res)
    """
    # set count to 20 since thats where result summary starts.
    """
    count = 20
    hits = []
    while True:
        current_line = result_text[count]
        count = count +1
        #checks to see if a "." is in the line and assumes from 20, it's an accession number. 
        if "." in current_line and ">" not in current_line:
            Accession_Number = current_line.split(" ")[0]
            Accession_Number = Accession_Number.split(".")[0]
            if args.verbose: print("  reading: "+processed_accesion_file_name)
            Isolates = pandas.read_csv(processed_accession_file_name,sep='\t')
           
            hits = hits+[Isolates.loc[Isolates["Accession_IDs"] == Accession_Number]]
        elif ">" in current_line:
            break
    return hits 
    """

def main():
    if args.verbose: print("main()")

    if args.mode == "VMR" or None:
        print("# load VMR")
        vmr_data = load_VMR_data()

        if args.verbose: print("# testing accession IDs")
        tested_accessions_ids = test_accession_IDs(vmr_data)
        
        if args.verbose: print("Writing", processed_accession_file_name)
        if args.verbose: print("\tColumn: ", tested_accessions_ids.columns)
        pandas.DataFrame.to_csv(tested_accessions_ids, processed_accession_file_name, sep='\t')

    if args.mode == "fasta" or None:
        print("# pull FASTAs from NCBI")
        fetch_fasta(processed_accession_file_name)

    if args.mode == "db" or None:
        print("# Query local VMR-E BLASTdb")
        query_database(args.query)

main()

if args.verbose: print("Done.")





    






