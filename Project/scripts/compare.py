#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compare.py

Description: Compares the outputs of AADR_parsing.py and clinvar_user_match.py 
to check for consistency in disease state determinations.

User-defined functions:  
    1. 

Non-standard modules: 
    1. pandas
    2. argparse
    3. os
    4. sys
    5. pathlib
    
Procedure:
    1. Load AADR_parsing output and clinvar_user_match output.
    2. Merge the two datasets on rsid and individual_id.
    3. Compare disease state determinations between the two datasets.
    4. Output a TSV file with the comparison results.

Input: Outputs of AADR_parsing.py and clinvar_user_match.py
Output: TSV file with comparison results

Usage: python3 compare.py -i <clinvar_user_match_output> -a <AADR_parsing_output> [-o <output_file>] [-d <output_dir>]
Name: Emma Bolech
Date: 2026-03-13

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
# Setting up command line arguments and flags using argparse.
###############################################################################

'''
argparse is used here to create a parser for this specific script.
    -i / --input_file:      path to parsed user SNP file (output of user_parser.py)
    -a / --AADR_file:       path to AADR_parsed.tsv file (output of AADR_parsing.py)
    -o / --output_file:     optional, output TSV file name (default: comparison_matches.tsv)
    -d / --output_dir:      optional, output directory (default: current directory)
'''

parser = argparse.ArgumentParser(prog="compare.py", description="Compare outputs of AADR_parsing.py and clinvar_user_match.py.")
parser.add_argument("-i", "--input_file", required=True, help="Path to parsed user SNP file (output of clinvar_user_match.py).")
parser.add_argument("-a", "--aadr_file", required=True, help="Path to AADR_parsed.tsv file (output of AADR_parsing.py).")
parser.add_argument("-o", "--output_file", nargs="?", default="comparison_matches.tsv", help="Output file name. Default is comparison_matches.tsv")
parser.add_argument("-d", "--output_dir", required=False, default=".", help="Output directory. Default is current directory.")
args = parser.parse_args()

for path, name in [(args.input_file, "Input file"), (args.aadr_file, "AADR file")]:
    try:
        if not Path(path).exists():
            raise Exception(f"- Error: {name} does not exist: {path}\nExiting...")
    except Exception as e:
        print(e)
        sys.exit()

output_path = os.path.join(args.output_dir, args.output_file)


#%%
###############################################################################
# Loading in user SNPs and AADR SNPs.
###############################################################################


user = pd.read_csv(args.input_file, sep="\t", dtype=str)
aadr = pd.read_csv(args.aadr_file, sep="\t", dtype=str)

# validate required columns
aadr_required = {"individual_id", "rsid", "disease_state", "GeneSymbol", "ClinicalSignificance", "PhenotypeList"}
user_required = {"rsid", "disease_state"}

# check if required columns are present in the dataframes
if not aadr_required.issubset(set(aadr.columns)):
    print(f"- Error: AADR file is missing required columns: {aadr_required - set(aadr.columns)}\nExiting...")
    sys.exit()

if not user_required.issubset(set(user.columns)):
    print(f"- Error: User file is missing required columns: {user_required - set(user.columns)}\nExiting...")
    sys.exit()  

#%%
###############################################################################
# Matching SNPs based on rsid.
###############################################################################

# slim down to just the columns we need for matching and rename for clarity
user_slim = user[["rsid", "disease_state"]].rename(columns={"disease_state": "user_disease_state"})

# merge the AADR and user dataframes on rsid to find matches
merged = aadr.merge(user_slim, on="rsid", how="inner")

# check how many matches we got and how many unique rsIDs are shared between the two datasets
print(f"--Check: {merged['rsid'].nunique()} shared rsID(s), {len(merged)} total row(s) after join")

# rearranging output columns and renaming for clarity
output = merged[[
    "individual_id",
    "rsid",
    "GeneSymbol",
    "ClinicalSignificance",
    "PhenotypeList",
    "disease_state",
    "user_disease_state"
]].rename(columns={"disease_state": "aadr_disease_state"})


#%%
###############################################################################
# Saving output.
###############################################################################

# save the output to a TSV file
output.to_csv(output_path, sep="\t", index=False)

print(f"--Done: Saved {len(output)} row(s) to {output_path}")






