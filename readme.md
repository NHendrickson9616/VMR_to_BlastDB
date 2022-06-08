<h1>VMR_to_BlastDB</h1>

This is a tool to extract data from the Virial Metadata Resource(VMR) published by the ICTV. 

<h2>Requirements</h2>

1. Python 3.*

3. Pandas, including openpyxl. 

5. Numpy

7. Biopython

9. A VMR excel file placed in the directory.

NOTE: For Blast capabilities, a verison of BLASTDB also needs to be installed. 

<h2>Usage</h3>
  My main script is VMR_to_fasta.py. The process is broken into steps via a argument dubbed "mode". 
  This isn't ideal, and eventually the script will infer mode based on the arguments given. 
  
  To parse the VMR and extract Accession numbers(mode:VMR):
  
    VMR_to_fasta.py -mode VMR -VMR_file_name [PATH_TO_VMR]
  
  To download fasta data(mode:fasta):
  
    python VMR_to_fasta.py -mode fasta -fasta_file_name [PATH_TO_EXCEL]
    
  TODO: Fix db mode. 

