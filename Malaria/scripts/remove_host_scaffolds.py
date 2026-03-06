#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
remove_host_scaffolds.py

Description: Removing specified host scaffolds from the input genome fasta file. The host identified scaffold names are provided in a separate text file. The output is a new fasta file with the host identified scaffolds removed.

User-defined functions:  
    1. host_id_scaffolds
    2. genome_filtering

Non-standard modules: 
    Path: imported from pathlib to define working directory.
    sys: imported to safely exit the program when necessary.
    re: imported to use regex modules to check for matching characters with ex. re.findall.
    math: imported to use natrual log equation.
    pandas: imported to store data in DataFrames and output files using meaningful formatting
    argparse: imported to define command line arguments using flags.
    os: imported to split output file names by the extension.
    
Procedure:
    1. Define functions.
    2. Establish working directory with Path.
    3. Set up command line arguments and flags using argparse.
    4. Inspect input files.
    5. Filter genome.
    6. Final message.

Input: Ht_genome_filtered.genome combined_scaffolds.txt
Output: Ht_malaria_only.genome (or specified output file name)

Usage: python3 remove_host_scaffolds.py Ht_genome_filtered.genome combined_scaffolds.txt Ht_malaria_only.genome
*** last argument is optional 
*** -o if last argument is not given "_.genome" is used by default

Version: 1.00
Date: 2026-03-04
Name: Emma Bolech
"""

#%%
###############################################################################
# Importing necessary non-standard modules.
###############################################################################

from pathlib import Path

import sys

import re

import argparse

import os

#%%
###############################################################################
# Defining functions.
###############################################################################

'''
host_id_scaffolds takes the input txt file and creates a list of the host identified scaffold names. 
This list is used in the genome_filtering function to filter out the host identified scaffolds from the input fasta file.
'''

def host_id_scaffolds(textfile: str) -> list: 
    with open(textfile, "r") as file:
        scaffolds_list = list() #create an empty list to append the scaffold names to
        for line in file:
            line = line.strip() #strip whitespace
            scaffolds_list.append(line) #append the scaffold names to the list

    return scaffolds_list #return the list of scaffold names


'''
genome_filtering takes the input fasta file and the list of host identified scaffold names 
and creates a dictionary of the contig names and sequences for the contigs that are not in the list 
of host identified scaffolds.
'''

def genome_filtering(fastafile: str, scaffolds_list: list) -> dict:
    with open(fastafile, "r") as file:
        genome_dict = dict() #create an empty dictionary to store the contig names and sequences
        for line in file:
            line = line.strip() #strip whitespace
            if line.startswith(">"): #if the line starts with a > it is a header line
                name_match = re.search(r"^>?(contig\d+)\b", line)
                if name_match:
                    name = name_match.group(1) #save the contig name
                    if name not in scaffolds_list: #if the contig name is not in the list of host identified scaffolds
                        genome_dict[line] = "" #add the contig name as a key to the dictionary with an empty string as the value
                        sequence = next(file).strip() #save the next line as the sequence
                        genome_dict[line] = sequence #add the sequence as the value for the contig name
                    else: #if the contig name is in the list of host identified scaffolds
                        continue #skip the line and the next line (the sequence)
            else: #if the line does not start with a > it is a sequence line
                continue #skip the line

        return genome_dict #return the dictionary of contig names and sequences
    


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
The third argument is defined to be the second input_file and should be flagged with -i in the command line.
This is required. The fourth argument is optional and is the output_file. By default it is _.genome. 

'''

parser = argparse.ArgumentParser(prog = "remove_host_scaffolds.py", description = "Removing host scaffolds from the input file.")
parser.add_argument("-f","--input_file", required=True, help = "Please provide an input genome fasta file. This should be a .genome file with the same format as the Ht_genome_filtered.genome file provided in the repository. The contig names should match those found in provided combined_scaffolds.txt file.")
parser.add_argument("-i", "--input_file02", required=True, help = "Please provide an input file with the same format as the combined_scaffolds.txt file provided in the repository. This file should contain the host identified contig names in a single column. The contig names should match those found in provided Ht_genome_filtered.genome file.")
parser.add_argument("-o","--output_file", nargs="?", default="no_host.genome", help = "Specify an output file. Default is no_host.genome")

args = parser.parse_args()  #setting the above variables as args

'''
Checks whether the input files exist in the working directory. Exits the program if not.
'''

input_path = Path(args.input_file)
try:
    if not input_path.exists():
        raise Exception ("- Error: The provided input fasta file does not exist or is missing from the working directory.\nExiting...")
except Exception as e:
    print(e)
    sys.exit()

input_path = Path(args.input_file02)
try:
    if not input_path.exists():
        raise Exception ("- Error: The provided input .txt file does not exist or is missing from the working directory.\nExiting...")
except Exception as e:
    print(e)
    sys.exit()
    

#%%
###############################################################################
# Inspecting input files.
###############################################################################

#checking for empty input file
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

#checking that the first line is a header in the FASTA file
with open(args.input_file) as file:
    for line in file:
        #if the first line starts with a > print the following check
        if line.strip().startswith(">"):
            print("- Check: The first line of the input file is a valid FASTA header line.")
            #break out of the for loop if this is satisfied
            break
        else:
            #if the first line does not start with an > print an error and exit the program safely
            print("- Error: The first line of the input file is not a valid FASTA header line.\nExiting...")
            sys.exit()

###############################################################################

#checking for 1 column in input txt file
with open (args.input_file02, "r") as file:
    
    for line in file:
        line_list = line.strip().split("\t")     #splits columns into a list by tabs
        try:
            if len(line_list) > 1: #if the length of columns is greater than 1 throw an error
                raise Exception("- Error: The input file is more than one column. Please inspect.\nExiting...")
        except Exception as e:
            print(e)
            sys.exit() #safely exit
            
    print("- Check: Input text file is correctly limited to 1 column.") #print in the end if passed
        

#%%
###############################################################################
# Filtering genome.
###############################################################################

'''
Applying the functions to the input fasta file and outputing a final filtered genome file.
'''

host_names = host_id_scaffolds(args.input_file02) #save the list of host identified scaffold names from the input txt file


genome = genome_filtering(args.input_file, host_names)

with open(args.output_file, "w") as file:
    for key, value in genome.items():
        file.write(f"{key}\n{value}\n") #write the contig name and sequence to the output file with a newline in between    
    

###############################################################################

# check whether output file contains any of the omitted contigs

with open(args.output_file, "r") as file:
    output_lines = list()
    for line in file:
        output_lines.append(line.strip()) #append each line in a list
    
    for scaffold in host_names:
        if scaffold in output_lines:
            print(f"- Error: The output file contains the following omitted contig: {scaffold}\nExiting...")
            sys.exit() #safely exit



#%%
###############################################################################
# Final message.
###############################################################################

#final print message that tells user what the input and output files are
print(f"Your input files are: {args.input_file} and {args.input_file02}\nYour output file is: {args.output_file}.")

print("Success!")


