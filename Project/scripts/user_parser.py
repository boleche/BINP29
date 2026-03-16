#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
user_parser.py

Description: This script is designed to parse SNP files from multiple consumer genomics providers 
(23andMe, AncestryDNA, and generic formats) into a standardized format. 
The script checks the format of the input file, extracts the relevant information (rsid, chromosome, position, allele1, allele2, and zygosity), 
and saves the standardized data to an output file in a tab-delimited format.

User-defined functions:  
    1. check_23_ancestry
    2. parse_23_andme
    3. parse_ancestry
    4. parse_format1
    5. parse_format2
    6. parse_format3
    7. parse_format4
    8. parse_generic


Non-standard modules: 
    1. pandas
    2. argparse
    3. os
    4. sys
    5. pathlib
    
Procedure:
    1. First, the script imports the necessary non-standard modules (pandas, argparse, os, sys, pathlib).
    2. Next, it defines several functions to check the format of the input file and parse it accordingly.
    3. The script then sets up command line arguments using argparse to specify the input file, output file name, and output directory.
    4. It checks if the input file exists and constructs the output path.
    5. The script determines the file type (23andMe, AncestryDNA, or generic) and calls the appropriate parser function to get a standardized dataframe.
    6. Finally, it saves the standardized dataframe to a tab-delimited file and prints a success message.

Input: A SNP file from 23andMe, AncestryDNA, or a generic format containing columns for rsid, chromosome, position, and genotype information.
Output: A standardized tab-delimited file containing columns for rsid, chromosome, position, allele1, allele2, and zygosity.

Usage: python3 user_parser.py -i input.snp -o output.tsv
Version: 1.00
Date: 2026-03-12
Name: Emma Bolech
"""

#%%
###############################################################################
# Importing necessary non-standard modules.
###############################################################################

import pandas as pd
import argparse
import os
import sys
from pathlib import Path


#%%
###############################################################################
# Defining functions.
###############################################################################

# function to check if the input file is from 23andMe or AncestryDNA or a generic format
def check_23_ancestry(user_file:str) -> str:
    with open(user_file, 'r') as f:
        for line in f:
            if "23andme" in line.lower():
                return "23andMe"
                break
            elif "ancestry" in line.lower():
                return "AncestryDNA"
                break
            else:
                return "Generic"
                break
                
# function to parse 23andMe files and return a dataframe with the relevant columns and zygosity information
def parse_23_andme(user_file:str) -> pd.DataFrame:
    df = pd.read_csv(user_file, sep="\t", comment="#", header=None, names=["rsid", "chromosome", "position", "genotype"])
    df["allele1"] = df["genotype"].str[0]
    df["allele2"] = df["genotype"].str[1]
    df["zygosity"] = df.apply(lambda row: "Homozygous" if row["allele1"] == row["allele2"] else "Heterozygous", axis=1)
    return df[["rsid", "chromosome", "position", "allele1", "allele2", "zygosity"]]

# function to parse AncestryDNA files and return a dataframe with the relevant columns and zygosity information
def parse_ancestry(user_file:str) -> pd.DataFrame:
    df = pd.read_csv(user_file, sep="\t", comment="#")
    df.columns = ["rsid", "chromosome", "position", "allele1", "allele2"]
    df["zygosity"] = df.apply(lambda row: "Homozygous" if row["allele1"] == row["allele2"] else "Heterozygous", axis=1)
    return df[["rsid", "chromosome", "position", "allele1", "allele2", "zygosity"]] 


# parser for """RSID, CHROMOSOME, POSITION, RESULT (e.g. AA, TT)""" format
# ex. test user 3 (comma separated):
def parse_format1(user_file:str) -> pd.DataFrame:
    df = pd.read_csv(user_file, sep=",", comment="#", header=0, names=["RSID", "CHROMOSOME", "POSITION", "RESULT"])
    df = df.rename(columns={"RSID": "rsid", "CHROMOSOME": "chromosome", "POSITION": "position", "RESULT": "genotype"})  
    df["allele1"] = df["genotype"].str[0]
    df["allele2"] = df["genotype"].str[1]
    df["zygosity"] = df.apply(lambda row: "Homozygous" if row["allele1"] == row["allele2"] else "Heterozygous", axis=1)
    return df[["rsid", "chromosome", "position", "allele1", "allele2", "zygosity"]]


# parser for """# name, chromosome, position, allele1, allele2""" format
# ex. test user 2 (comma separated):
def parse_format2(user_file:str) -> pd.DataFrame:
    df = pd.read_csv(user_file, sep=",", comment="#", header=0, names=["name", "chromosome", "position", "allele1", "allele2"])
    df = df.rename(columns={"name": "rsid"})
    df["zygosity"] = df.apply(lambda row: "Homozygous" if row["allele1"] == row["allele2"] else "Heterozygous", axis=1)
    return df[["rsid", "chromosome", "position", "allele1", "allele2", "zygosity"]]

# parser for """# rsid, chromosome, position, genotype (e.g. TT, AA)""" format
# ex. test user 1 (tab separated):
def parse_format3(user_file:str) -> pd.DataFrame:
    df = pd.read_csv(user_file, sep="\t", comment="#", header=0, names=["rsid", "chromosome", "position", "genotype"])
    df["allele1"] = df["genotype"].str[0]
    df["allele2"] = df["genotype"].str[1]
    df["zygosity"] = df.apply(lambda row: "Homozygous" if row["allele1"] == row["allele2"] else "Heterozygous", axis=1)
    return df[["rsid", "chromosome", "position", "allele1", "allele2", "zygosity"]]

# parser for """rsid, chromosome, position, allele1, allele2 (no comment header)""" format
# test user 4 and 5 (tab separated) no header:
def parse_format4(user_file:str) -> pd.DataFrame:
    df = pd.read_csv(user_file, sep="\t", header=0, names=["rsid", "chromosome", "position", "allele1", "allele2"])
    df["zygosity"] = df.apply(lambda row: "Homozygous" if row["allele1"] == row["allele2"] else "Heterozygous", axis=1)
    return df[["rsid", "chromosome", "position", "allele1", "allele2", "zygosity"]]
    

# function to check format of generic file and call the appropriate parser
def parse_generic(user_file:str) -> pd.DataFrame:
    with open(user_file, "r") as f:
        first_line = f.readline()
        if "RESULT" in first_line:
            return parse_format1(user_file)
        elif "name" in first_line:
            return parse_format2(user_file)
        elif "genotype" in first_line:
            return parse_format3(user_file)
        else:
            return parse_format4(user_file)

# main function to determine file type and call the appropriate parser, returning a standardized dataframe
def parse_user_file(user_file:str) -> pd.DataFrame:
    file_type = check_23_ancestry(user_file)
    if file_type == "23andMe":
        df = parse_23_andme(user_file)
    elif file_type == "AncestryDNA":
        df = parse_ancestry(user_file)
    else:
        df = parse_generic(user_file)
    
    # remove non rs ids
    df = df[df["rsid"].str.startswith("rs")]
    # keep only rows where both alleles are standard nucleotides
    valid_nucs = {"A", "T", "C", "G"}
    df = df[df["allele1"].isin(valid_nucs) & df["allele2"].isin(valid_nucs)]

    return df

#%%
###############################################################################
# Setting up command line arguments and flags using argparse.
###############################################################################

# adding wrapper so i can call this in the app.py streamlit app without it trying to parse arguments when the app is run
if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="user_parser.py", description="Parse SNP files from multiple consumer genomics providers into a standardized format.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to input SNP file (23andMe, AncestryDNA, or generic format).")
    parser.add_argument("-o", "--output_file", nargs="?", default="standardized_snps.tsv", help="Output file name. Default is standardized_snps.tsv")
    parser.add_argument("-d", "--output_dir", required=False, default=".", help="Output directory. Default is current directory.")
    args = parser.parse_args()

    # check input file exists
    try:
        if not Path(args.input_file).exists():
            raise Exception(f"- Error: Input file does not exist: {args.input_file}\nExiting...")
    except Exception as e:
        print(e)
        sys.exit()

    output_path = os.path.join(args.output_dir, args.output_file)

    #%%
    ###############################################################################
    # Determine file type and output parsed df.
    ###############################################################################

    # check if the input file is from 23andMe, AncestryDNA, or a generic format and call the appropriate parser
    df = parse_user_file(args.input_file)


    #%%
    ###############################################################################
    # Output parsed df.
    ###############################################################################

    # save the standardized dataframe to a tab-delimited file

    df.to_csv(output_path, sep="\t", index=False)
    print(f"Successfully parsed {args.input_file} and saved standardized data to {output_path}")



