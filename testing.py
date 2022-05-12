import pandas
import numpy
import os
master_list = pandas.read_excel('Master_Species.xlsx',sheet_name="ICTV2020 Master Species List#36")
col = master_list.columns
Master_Species = master_list["Species"].to_numpy()
vmr_source = pandas.read_excel("VMR_Source_1.xlsx")
VMr_list = pandas.read_excel("fixed_vmr.xlsx")
VMR_Species = VMr_list['Species'].to_numpy()
Genome_status = VMr_list['Genome coverage'].to_numpy()
print(len(VMR_Species))
print(len(Master_Species))
arr=numpy.setdiff1d(Master_Species,VMR_Species)
print(VMr_list.loc[VMr_list['Species']=="Acidianus filamentous virus 3",["Genome coverage"]]["Genome coverage"].values[0])
print(arr) 

def validate():
    raw_master_species = pandas.read_excel("Master_Species.xlsx",sheet_name="ICTV2020 Master Species List#36")
    col = raw_master_species.columns
    Master_Species = raw_master_species["Species"].to_numpy()
    fasta_dirs=os.listdir("fasta/")
    virus_names = []
    for fasta_dir in fasta_dirs:
        
        if ".fa" in fasta_dir:
            fasta_file = open("fasta/"+fasta_dir,"r")
            header = fasta_file.readline()
            virus_name = header.split(" ")[0]
            virus_name = virus_name.split("#")[0]
            virus_name = virus_name.replace("_"," ")
            virus_name = virus_name.replace(">","")
            #virus_name = virus_name.replace("-","")
            if virus_name not in Master_Species:
                if virus_name == "":
                    print("The fasta file "+str(fasta_dir) +" is empty!")
                else:
                    print("The virus "+virus_name+" is present but not in the Master Species file. file name: "+str(fasta_dir))
            virus_names = virus_names+[virus_name.lower()]
            fasta_file.close()
    for master_specie in Master_Species:
        if master_specie.lower() not in virus_names:
        
            if VMr_list.loc[VMr_list['Species']==master_specie,["Genome coverage"]]["Genome coverage"].values[0] != "No entry in Genbank":
                
                print("The virus "+master_specie+" does not have a fasta file but should according to the VMR")
            

            
            
validate()




    
        
