import pandas
import subprocess
from Bio import Entrez
import time
from urllib import error
import argparse
# variables
parser = argparse.ArgumentParser()
parser.add_argument('-email',help="email for entrez to use")
parser.add_argument('-VMR_file_name',help='path to VMR file to use.')
parser.add_argument('-FASTA_file_name',help='path to output FASTA')
args = parser.parse_args()
email = "halko@uab.edu"
VMR_file_name = 'VMR_E_data.xlsx'
FASTA_file_name = 'all.fa'
query = 'MN876842.fasta'

###############################################################################################################
# Loads excel from https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/ and puts it into a DataFrame
############################################################################################################### 
def load_VMR_data():
    # Importing excel sheet as a DataFrame. Requires xlrd and openpyxl package
    raw_vmr_data = pandas.read_excel(VMR_file_name)
    

    # list of the columns to extract from raw_vmr_data
    vmr_cols_needed = ['Virus GENBANK accession','Species']
    # compiling new dataframe from vmr_cols_needed
    truncated_vmr_data = pandas.DataFrame(raw_vmr_data[[col_name for col_name in vmr_cols_needed]])
    return truncated_vmr_data

##############################################################################################################
# Cleans Accession numbers assuming the following about the accession numbers:
# 1. Each Accession Number is 6-8 characters long
# 2. Each Accession Number contains at least 3 numbers
# 3. Each Accession Number contains at most 3 letters
# 4. Accession Numbers in the same block are seperated by a ; or a :
##############################################################################################################

def test_accession_IDs(df):
    # defining new DataFrame before hand
    processed_accession_IDs = pandas.DataFrame(columns=['Species','Accession_IDs','segment'])
    # for loop for every entry in given dataframe
    for entry_count in range(0,len(df.index)):
        # Find current species in this row
        Species=df['Species'][entry_count]
        #initial processing of raw accession numbers.
        accession_ID = df['Virus GENBANK accession'][entry_count]
        accession_ID = str(accession_ID).replace(" ","")
        accession_ID = accession_ID.split(';')
        accession_ID = [accession_part.split(':') for accession_part in accession_ID]
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
                #checks if current selection fits schema of a accession number
                if len(str(accession_part)) == 8 or 6 and letter_count<3 and number_count>3:
                    processed_accession_ID = accession_part
                    if processed_accession_ID in processed_accession_IDs['Accession_IDs'].values:
                        print(processed_accession_ID+" is already in array! Duplicate found.")
                    processed_accession_IDs.loc[len(processed_accession_IDs.index)] = [Species,processed_accession_ID,segment]
                else:
                    segment = accession_part
    return processed_accession_IDs               
#test_accession_IDs(truncated_vmr_data).to_excel('processed_accessions.xlsx')
#######################################################################################################################################
# Utilizes Biopython's Entrez API to fetch FASTA data from Accession numbers. Prints Accession Numbers that failed to 'clean' correctly
#######################################################################################################################################  
def fetch_fasta():
    bad_accessions = []
    Accessions = pandas.read_excel('processed_accessions.xlsx')
    all_reads = []
    # Fetchs FASTA data for every accession number
    for accession_ID in Accessions['Accession_IDs']:
        time.sleep(1)
        fa_file = open('fasta/all_alt.fa','a')
        Entrez.email = email
        try:
            handle = Entrez.efetch(db="nuccore", id=accession_ID, rettype="fasta", retmode="text")
        # prints out accession that got though cleaning
        except(error.HTTPError):
            print(accession_ID)
            bad_accessions = bad_accessions+[accession_ID]
        read = handle.read()
        fa_file.write(read)
        fa_file.close()
#######################################################################################################################################
# Calls makedatabase.sh. Uses 'all.fa'
#######################################################################################################################################
def make_database():
    p1 = subprocess.run("makedatabase.sh")
    
#######################################################################################################################################
# BLAST searches a given FASTA file and returns DataFrame rows with from the accession numbers. Returns in order of significance.  
#######################################################################################################################################
def query_database(path_to_query):
    p1 = subprocess.run(["bash","query_database.sh",path_to_query])

    results = open("output/results.out","r")
    result_text = results.readlines()
    results.close()
    # set count to 20 since thats where result summary starts. 
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

def main():
    try:
        pandas.read_excel('processed_accessions.xlsx')
    except:
        test_accession_IDs(load_VMR_data())
    try:
        open("fasta/all_alt.fa")
    except:
        fetch_fasta()
    query_database("fasta/"+query)

main()






    







