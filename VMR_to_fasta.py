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
parser.add_argument('-verbose',help="printout details during run",action=argparse.BooleanOptionalAction)
parser.add_argument('-file',help="optional argument. Name of the file to get arguments from.",type=open, action=LoadFromFile)
parser.add_argument("-email",help="email for Entrez to use when fetching Fasta files")
parser.add_argument("-mode",help="what function to do. Options: VMR,fasta,db")
parser.add_argument("-ea",help="Fetch E or A records (Exemplars or AdditionalIsolates)", default="E")
parser.add_argument("-fasta_file_name",help="Name of the fasta file to output",default="./fasta/all_alt.fa")
parser.add_argument("-VMR_file_name",help="name of the VMR file to load.",default="VMR_E_data.xlsx")
parser.add_argument("-query",help="Name of the fasta file to query the database")
args = parser.parse_args()
mode = args.mode
if mode != 'fasta' and mode != "VMR" and mode != "db":
    print("Valid mode not selected. Options: VMR,fasta,db",file=sys.stderr)
#Takes forever to import so only imports if it's going to be needed
if mode == 'fasta':
    from Bio import Entrez
query = args.query
fasta_file_name = args.fasta_file_name
VMR_file_name = args.VMR_file_name
#Catching error
if mode == "db":
    if args.query == None:
        print("Database Query mode is selected but no fasta file was specified! Please set the '-fasta_file_name' or change mode.",file=sys.stderr)
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
    if args.verbose: print("  opening", VMR_file_name)

    # Importing excel sheet as a DataFrame. Requires xlrd and openpyxl package
    try:
        raw_vmr_data = pandas.read_excel(VMR_file_name,engine='openpyxl')

    # list of the columns to extract from raw_vmr_data
        vmr_cols_needed = ['Virus GENBANK accession','Species','Exemplar or additional isolate','Genome coverage','Genus']

        if args.verbose: print("VMR data loaded.")
    except(FileNotFoundError):
        print("The VMR file specified does not exist! Make sure the path set by '-VMR_file_name' is correct.",file=sys.stderr)
    

    # compiling new dataframe from vmr_cols_needed
    truncated_vmr_data = pandas.DataFrame(raw_vmr_data[[col_name for col_name in vmr_cols_needed]])
    # DataFrame.loc is helpful for indexing by row. Allows expression as an argument. Here, 
    # it finds every row where 'E' is in column 'Exemplar or additional isolate' and returns 
    # only the columns specified. 
    vmr_data = truncated_vmr_data.loc[truncated_vmr_data['Exemplar or additional isolate']==args.ea.upper(),['Species','Virus GENBANK accession',"Genome coverage","Genus"]]

    # only works when I reload the vmr_data, probably not necessary. have to look into why it's doing this. 
    VMR_hack_file_name = "fixed_vmr_"+args.ea+".xlsx"
    if args.verbose: print("Writing"+VMR_hack_file_name,": workaround - filters VMR down to "+args.ea.upper()+" records only")
    vmr_data.to_excel(VMR_hack_file_name)
    if args.verbose: print("Loading",VMR_hack_file_name)
    raw_vmr_data = pandas.read_excel(VMR_hack_file_name,engine='openpyxl')
    # Removing Genome Coverage column from the returned value. 
    vmr_cols_needed = ['Virus GENBANK accession','Species',"Genus"]
    truncated_vmr_data = pandas.DataFrame(raw_vmr_data[[col_name for col_name in vmr_cols_needed]])

    
    #truncated_vmr_data = truncated_vmr_data.drop(columns=['Exemplar or additional isolate'])
    
    return truncated_vmr_data


def test_accession_IDs(df):
    if args.verbose: print("test_accession_IDs()")
##############################################################################################################
# Cleans Accession numbers assuming the following about the accession numbers:
# 1. Each Accession Number is 6-8 characters long
# 2. Each Accession Number contains at least 3 numbers
# 3. Each Accession Number contains at most 3 letters
# 4. Accession Numbers in the same block are seperated by a ; or a , or a :
##############################################################################################################
    # defining new DataFrame before hand
    processed_accession_IDs = pandas.DataFrame(columns=['Species','Accession_IDs','segment',"Genus"])
    # for loop for every entry in given processed_accession_IDs
    for entry_count in range(0,len(df.index)):
        # Find current species in this row
        Species=df['Species'][entry_count]
        # Find current Genus in this row
        Genus=df['Genus'][entry_count]
        #initial processing of raw accession numbers.
        accession_ID = df['Virus GENBANK accession'][entry_count]
        #removing whitespace.
        accession_ID = str(accession_ID).replace(" ","")
        # instead of trying to split by commas and semicolons, I just replace the commas with semicolons. 
        accession_ID.replace(",",";")
        accession_ID = accession_ID.split(';')
        
        
        #split by colons too

        accession_ID = [accession_part.split(':') for accession_part in accession_ID]
    
        for accession_part in accession_ID:
            #flag long strings as being suspect. Most aren't 10 characters long. 
            if len(accession_part) > 10:
                print('suspiciously long accession number or segment name. Please verify its correct:'+str(accession_part),file=sys.stderr)


        # for loop for every ";" split done
        for accession_seg in accession_ID:
            segment = None
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
                    processed_accession_IDs.loc[len(processed_accession_IDs.index)] = [Species,processed_accession_ID,segment,Genus]
                    #print("'"+processed_accession_ID+"'"+' has been cleaned.')
                else:
                    segment = accession_part
    return processed_accession_IDs               
#test_accession_IDs(truncated_vmr_data).to_excel('processed_accessions.xlsx')
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
    
    #Check to see if fasta data exists and, if it does, loads the accessions numbers from it into an np array.
    bad_accessions = []

    if args.verbose: print("  loading:", processed_accession_file_name)
    Accessions = pandas.read_excel(processed_accession_file_name,engine='openpyxl')

    all_reads = []

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
            accession_file_name = genus_dir+"/"+str(accession_ID)+".fa"
            
            # make sure dir exists
            if not os.path.exists(genus_dir):
                # Create the directory if it doesn't exist
                os.makedirs(genus_dir)
                if args.verbose: print(f"Directory '{genus_dir}' created successfully.")
    
            try:
                fa_file = open(accession_file_name,'r')
                fa_file.close()
                print("file found for "+genus_name+"/"+str(accession_ID),"[SKIP]")
            except(FileNotFoundError):
                fa_file = open(accession_file_name,'w')
                print("fetching fasta file for "+genus_name+"/"+str(accession_ID))
                try:
                    # fetch FASTA from NCBI
                    handle = Entrez.efetch(db="nuccore", id=accession_ID, rettype="fasta", retmode="text")

                    # limit requests: 3/second with email, 10/second with API_KEY
                    time.sleep(entrez_sleep)

                    print('fasta for '+accession_ID+ ' obtained.')

                    # prints out accession that got though cleaning
                
                    read = handle.read()
                    desc_line = read.split("\n")[0].split(" ",1)[1]
                    if str(segment).lower() == "":
                        desc_line = ">"+""+str(Species.replace(" ","_"))                 +" "+str(accession_ID)+" "+desc_line.replace(">","")
                    else:
                        desc_line = ">"+""+str(Species.replace(" ","_"))+"#"+str(segment)+" "+str(accession_ID)+" "+desc_line.replace(">","")
                    print(desc_line)
                    fa_file.write(desc_line+"\n"+read.split("\n",1)[1])
                
                    fa_file.close()
                except:
                    print("Accession ID"+"'"+str(accession_ID)+"'"+"did not get properly cleaned. Accession Cleaning Heuristic needs editing.",file=sys.stderr)
                    bad_accessions = bad_accessions+[accession_ID]
            count=count+1


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

    p1 = subprocess.run(["bash","query_database.sh",path_to_query,results_file])
    """
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
            Isolates = pandas.read_excel('processed_accessions.xlsx')
           
            hits = hits+[Isolates.loc[Isolates["Accession_IDs"] == Accession_Number]]
        elif ">" in current_line:
            break
    return hits 
    """

def main():
    if args.verbose: print("main()")

    processed_accession_file_name ="processed_accessions.xlsx"
    if args.ea == "a": 
        processed_accession_file_name ="processed_accessions_a.xlsx"        

    if mode == "VMR" or None:
        print("# load VMR")
        vmr_data = load_VMR_data()

        if args.verbose: print("# testing accession IDs")
        tested_accessions_ids = test_accession_IDs(vmr_data)
        
        if args.verbose: print("Writing", processed_accession_file_name)
        pandas.DataFrame.to_excel(tested_accessions_ids,processed_accession_file_name)

    if mode == "fasta" or None:
        print("# pull FASTAs from NCBI")
        fetch_fasta(processed_accession_file_name)

    if mode == "db" or None:
        print("# Query local VMR-E BLASTdb")
        query_database(args.query)

main()

if args.verbose: print("Done.")





    







