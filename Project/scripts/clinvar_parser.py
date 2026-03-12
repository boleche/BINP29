#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
clinvar_parser.py

Description: Parsing ClinVar variant data from "variant_summary.txt" found in ClinVar database and converting it to a more usable format.

User-defined functions:  
    1. filtered_rows

Non-standard modules: 
    Path: imported from pathlib to define working directory.
    sys: imported to safely exit the program when necessary.
    argparse: imported to define command line arguments using flags.
    os: imported to split output file names by the extension.
    
Procedure:
    1. Parse the input file and filter rows based on specified criteria.
    2. Save the entries of interest for each filtered row in a list and append that list to a list of filtered rows.
    3. Check that the headers of interest are in the expected location in the input file
    4. Write the clean headers and the filtered rows to an output file in a tab delimited format.
    5. Print a final message to the user with the input and output file names.

Input: variant_summary.txt (downloaded from ClinVar database and unzipped)
Output: ClinVar_parsed.tsv (or specified output file name)

Usage: python3 clinvar_parser.py -i variant_summary.txt -o ClinVar_parsed.tsv -d output_directory
Version: 1.00
Date: 2026-03-10
Name: Emma Bolech
"""

#%%
###############################################################################
# Importing necessary non-standard modules.
###############################################################################

from pathlib import Path

import sys

import argparse

import os

#%%
###############################################################################
# Defining functions.
###############################################################################

'''
filtered_rows function: takes the input file and filters the rows based on the following criteria:
    1. Type is "single nucleotide variant"
    2. ClinicalSignificance is "Pathogenic" or "Likely pathogenic"
    3. RS# (dbSNP) is not "-"
    4. Assembly is "GRCh37"
The function saves the entries of interest for each row in a list and appends that list to a list of filtered rows. The function returns the list of filtered rows.
'''

def filtered_rows(input_file:str) -> list:
# save each filtered row as a list within a list of lists (filtered_rows)
    with open(input_file, "r") as file:
        next(file) #skip the first line with the headers
        filtered_rows = list() #create an empty list to save the filtered rows
        for line in file:
            columns = line.strip().split("\t") #split the line into a list by tabs
            entries = list() #create an empty list to save the entries of interest for each line
            if (columns[1] == "single nucleotide variant" and
                columns[6] == "Pathogenic" or columns[6] == "Likely pathogenic" and
                columns[9] != "-" and
                columns[16] == "GRCh37"):
                entries.append(columns[1]) #append the entry to the entries list
                entries.append(columns[4])
                entries.append(columns[6])
                entries.append(columns[9])
                entries.append(columns[13])   
                entries.append(columns[16])
                entries.append(columns[18])
                entries.append(columns[19])
                entries.append(columns[21])
                entries.append(columns[22])
                entries.append(columns[24])
                entries.append(columns[30])
                entries.append(columns[31])
                entries.append(columns[32])
                entries.append(columns[33])
                filtered_rows.append(entries) #append the entries list to the filtered rows list
        return filtered_rows #return the list of filtered rows



#%%
############################################################################### 
# Establish working directory with Path.
###############################################################################

#establish working directory and save into workingdirectory variable
workingdirectory = Path.cwd()   
print(f"This is your working directory: {workingdirectory}")    


#%%
###############################################################################
# Setting up command line arguments and flags using argparse.
###############################################################################

'''
argparse is used here to create a parser for this specific script. The second argument
(after the script name) is defined to be the input_file and should be flagged with -f in 
the command line. This argument is required. The help section is written to provide guidance to the user.
The third argument is optional and is the output_file. By default it is ClinVar_parsed.txt. 

'''

parser = argparse.ArgumentParser(prog = "clinvar_parser.py", description = "Parsing ClinVar variant data from 'variant_summary.txt' found in ClinVar database and converting it to a more usable format.")
parser.add_argument("-f","--input_file", required=True, help = "Please provide an input ClinVar variant_summary.txt file. This can be fount at: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz. Please unzip the file before using it as an input.")
parser.add_argument("-o","--output_file", nargs="?", default="ClinVar_parsed.tsv", help = "Specify an output file. Default is ClinVar_parsed.tsv")
parser.add_argument("-d", "--output-dir", required=False, default='.', help='Output directory')

args = parser.parse_args()  #setting the above variables as args

'''
Checks whether the input files exist in the working directory. Exits the program if not.
'''

input_path = Path(args.input_file)
try:
    if not input_path.exists():
        raise Exception ("- Error: The provided input ClinVar datafile does not exist or is missing from the working directory.\nExiting...")
except Exception as e:
    print(e)
    sys.exit()
    

#%%
###############################################################################
# Inspecting input files.
###############################################################################

# checking for empty input file
with open(args.input_file, "r") as file:
    total_lines = list()
    for line in file:
        total_lines.append(line) #append each line in a list
    
    length_input = len(total_lines) #save the number of total lines
    try:
        if length_input == 0: #if 0 assume empty and throw an error
            raise Exception("- Error: The input file is empty.\nExiting...")
    except Exception as e:
        print(e)
        sys.exit() #safely exit

# checking that first head line contains 43 tab deliminated columns
with open(args.input_file) as file:
    for line in file:
        columns = line.strip().split("\t") #split the first line into a list by tabs
        if len(columns) != 43: #if the length of the list is not 43 throw an error
            print("- Error: The first line of the input file does not contain 43 tab delimited columns.\nExiting...")
            sys.exit() #safely exit
        else:
            print("- Check: The first line of the input file contains 43 tab delimited columns.")
            break #break out of the for loop if this is satisfied

#%%
###############################################################################
# Saving the header and writing out filtered entry rows to output file.
###############################################################################

# check whether the headers of interest are in the right location and if yes append to clean headers list
with open(args.input_file, "r") as file:
    for line in file:
        headers = line.strip().split("\t") #split the first line into a list by tabs
        clean_headers = list() #create an empty list to save the clean headers
        if (headers[1] == "Type" and 
        headers[4] == "GeneSymbol" and
        headers[6] == "ClinicalSignificance" and
        headers[9] == "RS# (dbSNP)" and
        headers[13] == "PhenotypeList" and
        headers[16] == "Assembly" and
        headers[18] == "Chromosome" and
        headers[19] == "Start" and
        headers[21] == "ReferenceAllele" and
        headers[22] == "AlternateAllele" and
        headers[24] == "ReviewStatus" and
        headers[30] == "VariationID" and
        headers[31] == "PositionVCF" and
        headers[32] == "ReferenceAlleleVCF" and
        headers[33] == "AlternateAlleleVCF"):
            clean_headers.append(headers[1]) #append the header to the clean headers list
            clean_headers.append(headers[4])
            clean_headers.append(headers[6])
            clean_headers.append(headers[9])
            clean_headers.append(headers[13])
            clean_headers.append(headers[16])
            clean_headers.append(headers[18])
            clean_headers.append(headers[19])
            clean_headers.append(headers[21])
            clean_headers.append(headers[22])
            clean_headers.append(headers[24])
            clean_headers.append(headers[30])
            clean_headers.append(headers[31])
            clean_headers.append(headers[32])
            clean_headers.append(headers[33])
            break #break out of the for loop if this is satisfied
        else:
            print("- Error: The headers of interest are not in the expected location in the input file.\nExiting...")
            sys.exit() #safely exit 
    print("- Check: The headers of interest are in the expected location in the input file.")
   
###############################################################################

filtered_rows_list = filtered_rows(args.input_file) #save the list of filtered rows from the function   

###############################################################################

# join the output directory and output file name to create the full output path
output_path = os.path.join(args.output_dir, args.output_file)

with open(output_path, "w") as file:
    file.write("\t".join(clean_headers) + "\n") #write the clean headers to the output file with tabs in between and a newline at the end
    for row in filtered_rows_list:
        file.write("\t".join(row) + "\n") #write each row of entries to the output file with tabs in between and a newline at the end   



#%%
###############################################################################
# Final message.
###############################################################################

#final print message that tells user what the input and output files are
print(f"Your input file is: {args.input_file}\nYour output file is: {args.output_file}.")

print("Success!")


