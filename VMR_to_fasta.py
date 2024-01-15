#!/usr/bin/env python3
import pandas
import subprocess
import time
from urllib import error
import argparse
import numpy as np
import sys
# Class needed to load args from files. 
class LoadFromFile (argparse.Action):
    def __call__ (self, parser, namespace, values, option_string = None):
        with values as f:
            # parse arguments in the file and store them in the target namespace
            parser.parse_args(f.read().split(), namespace)
parser = argparse.ArgumentParser(description="")
#setting arguments.
parser.add_argument('-file',help="optional argument. Name of the file to get arguments from.",type=open, action=LoadFromFile)
parser.add_argument("-email",help="email for Entrez to use when fetching Fasta files")
parser.add_argument("-mode",help="what function to do. Options: VMR,fasta,db")
parser.add_argument("-fasta_file_name",help="Name of the fasta file to output",default="/fasta/all_alt.fa")
parser.add_argument("-VMR_file_name",help="name of the VMR file to load.",default="VMR_E_data.xlsx")
parser.add_argument("-query",help="Name of the fasta file to query the database")
args = parser.parse_args()
email = args.email
mode = args.mode
#if mode != 'fasta' or "VMR" or "db":
#    print("Valid mode not selected. Options: VMR,fasta,db",file=sys.stderr)
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
############################################################################################################### 
# DataFrame['column name'] = provides entire column
# DataFrame['column name'][0,1,2,3,4,5 etc] provides row for that column
# 
#
def load_VMR_data():
    # Importing excel sheet as a DataFrame. Requires xlrd and openpyxl package
    try:
        raw_vmr_data = pandas.read_excel(VMR_file_name,engine='openpyxl')

    # list of the columns to extract from raw_vmr_data
        vmr_cols_needed = ['Virus GENBANK accession','Species','Exemplar or additional isolate','Genome coverage','Genus']

        print("VMR data loaded.")
    except(FileNotFoundError):
        print("The VMR file specified does not exist! Make sure the path set by '-VMR_file_name' is correct.",file=sys.stderr)
    

    # compiling new dataframe from vmr_cols_needed
    truncated_vmr_data = pandas.DataFrame(raw_vmr_data[[col_name for col_name in vmr_cols_needed]])
    # DataFrame.loc is helpful for indexing by row. Allows expression as an argument. Here, 
    # it finds every row where 'E' is in column 'Exemplar or additional isolate' and returns 
    # only the columns specified. 
    vmr_data = truncated_vmr_data.loc[truncated_vmr_data['Exemplar or additional isolate']=='E',['Species','Virus GENBANK accession',"Genome coverage","Genus"]]
    # only works when I reload the vmr_data, probably not necessary. have to look into why it's doing this. 
    vmr_data.to_excel("fixed_vmr.xlsx")
    raw_vmr_data = pandas.read_excel('fixed_vmr.xlsx',engine='openpyxl')
    # Removing Genome Coverage column from the returned value. 
    vmr_cols_needed = ['Virus GENBANK accession','Species',"Genus"]
    truncated_vmr_data = pandas.DataFrame(raw_vmr_data[[col_name for col_name in vmr_cols_needed]])

    
    #truncated_vmr_data = truncated_vmr_data.drop(columns=['Exemplar or additional isolate'])
    
    return truncated_vmr_data


def test_accession_IDs(df):
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
# Utilizes Biopython's Entrez API to fetch FASTA data from Accession numbers. Prints Accession Numbers that failed to 'clean' correctly
#######################################################################################################################################  
def parse_taxname(raw_xml):
    # Fetchs FASTA data for every accession number
    read = raw_xml
    read = read.split('taxname')[1]
    read = read.split(',')[0]
    return read
def fetch_fasta():
    #Check to see if fasta data exists and, if it does, loads the accessions numbers from it into an np array.
    bad_accessions = []
    Accessions = pandas.read_excel('processed_accessions.xlsx',engine='openpyxl')
    all_reads = []
    # Fetches FASTA data for every accession number
    count = 0
    for accession_ID in Accessions['Accession_IDs']:
            row = Accessions.loc[count]
            Species = row[1]
            accession_ID = row[2]
            segment = row[3]
            print(segment)
            if str(segment).lower() == "nan":
                segment = ""
            #needed to limit requests to 3 per second. can be upped to 10 with auth token. 
            time.sleep(.34)
            # loading fasta_file into memory
            try:

                fa_file = open("fasta_new_vmr/"+accession_ID+".fa",'r')
                fa_file.close()
                print("file found for "+str(accession_ID))

            except(FileNotFoundError):
                fa_file = open("fasta_new_vmr/"+str(accession_ID)+".fa",'w')
                print("creating fasta file for "+str(accession_ID))
                Entrez.email = email
                try:
                
                    handle = Entrez.efetch(db="nuccore", id=accession_ID, rettype="fasta", retmode="text")

                    print('fasta for '+accession_ID+ ' obtained.')

            # prints out accession that got though cleaning
                
                    read = handle.read()
                    desc_line = read.split("\n")[0].split(" ",1)[1]
                    if str(segment).lower() == "":
                        desc_line = ">"+""+str(Species.replace(" ","_"))+" "+str(accession_ID)+" "+desc_line.replace(">","")
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
    p1 = subprocess.run(["bash","query_database.sh",path_to_query])
    """
    results = open("output/results.out","r")
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

    if mode == "VMR" or None:
        pandas.DataFrame.to_excel(test_accession_IDs(load_VMR_data()),"processed_accessions.xlsx")

    if mode == "fasta" or None:
        fetch_fasta()
    if mode == "db" or None:
        query_database("query/"+query)

main()






    







