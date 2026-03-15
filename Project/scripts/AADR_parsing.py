#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AADR_parsing.py

Description: Parsing ClinVar data and matching it to the AADR data to extract genotype information for the matched SNPs. Filtering out modern samples and 
adding carrier status information for recessive conditions.

User-defined functions:  
    1. load_clinvar
    2. load_plink_files
    3. find_zygosity
    4. determine_disease_state

Non-standard modules: 
    1. sqlite3
    2. pandas
    3. numpy
    4. pandas_plink
    5. argparse
    6. os
    7. sys
    8. pathlib
    
Procedure:
    1. Load the AADR data using the load_plink_files function and the ClinVar data using the load_clinvar function.
    2. Replace the rs in the rsid column of the bim dataframe with an empty 
        string to get just the numeric part of the rs number, and create new columns in the bim dataframe.
    3. Load the bim and clinvar dataframes into an in-memory SQLite database and query the database to 
        find matches between the bim and clinvar tables based on rs ID and allele info, storing the results in a new dataframe called sql_match.
    4. Extract the rows from the G matrix (SNPs) that match what is specified in the sql_match dataframe 
        and store those rows in a new variable called G_matched.
    5. Create a new dataframe called G_df from the G_matched variable and add SNP metadata columns from the sql_match dataframe.
    6. Melt the G_df dataframe from wide to long format to create a new dataframe
        called sql_match_genotypes that contains one row per SNP x individual, and drop any rows with missing genotypes.
    7. Filter out modern samples from the sql_match_genotypes dataframe using the second column of the modern_samples.txt file, and print the number of rows remaining after filtering.
    8. Apply the find_zygosity and determine_disease_state functions to the sql_match_genotypes
        dataframe to add new columns for genotype and disease state information.
    9. Save the sql_match_genotypes dataframe to a new file called 'ClinVar_AADR_genotypes.tsv' in a 
        tab delimited format and also save it to the SQLite database as a new table called 'clinvar_aadr_genotypes', 
        and print the number of rows saved to the file.
    

Input: PLINK files, ClinVar data, modern samples file
Output: ClinVar_AADR_genotypes.tsv

Usage: python3 AADR_parsing.py -p <plink_prefix> -c <clinvar_file> -m <modern_samples> [-o <output_file>] [-d <output_dir>]
Date: 2026-03-10
Name: Emma Bolech
"""
#%%
###############################################################################
# Importing necessary non-standard modules.
###############################################################################

import sqlite3
import pandas as pd
import numpy as np
from pandas_plink import read_plink
import argparse
import os
import sys
from pathlib import Path

#%%
###############################################################################
# Defining functions.
###############################################################################

# load ClinVar parsed data
def load_clinvar(clinvar_path):
    clinvar = pd.read_csv(clinvar_path, sep='\t', dtype={
        'Type': str,
        'GeneSymbol': str,
        'ClinicalSignificance': str,
        'RS# (dbSNP)': str,
        'PhenotypeList': str,
        'Assembly': str,
        'Chromosome': str,
        'Start': int,
        'ReferenceAllele': str,
        'AlternateAllele': str,
        'ReviewStatus': str,
        'VariationID': str,
        'PositionVCF': str,
        'ReferenceAlleleVCF': str,
        'AlternateAlleleVCF': str,
        'Dominance': str
    })
    return clinvar

# load PLINK files
def load_plink_files(prefix):
    # read the PLINK files using the pandas_plink library and rename the columns in the bim and fam dataframes to more descriptive names
    bim, fam, G = read_plink(prefix)
    bim_df = bim.rename(columns={
        "chrom": "chromosome",
        "snp":   "rsid",
        "cm":    "cM",
        "pos":   "position",
        "a0":    "alt_allele",
        "a1":    "ref_allele"
    })
    fam_df = fam.rename(columns={
        "fid":    "family_id",
        "iid":    "individual_id",
        "father": "paternal_id",
        "mother": "maternal_id",
        "gender": "sex",
        "trait":  "phenotype"
    })
    return G, bim_df, fam_df


# building the genotype rows (for recessive carrier info)
# create a new column in the sql_match dataframe called 'genotype' that contains the genotype information for each matched SNP based on the values in the G_matched matrix

def find_zygosity(G_col_value: int) -> str:
    if G_col_value == 0:
        return 'homozygous reference'
    elif G_col_value == 1:
        return 'heterozygous'
    elif G_col_value == 2:
        return 'homozygous alternate'
    else:
        return 'unknown'
    
# function which adds to a column called 'disease_state' in the sql_match dataframe which indicates whether the individual is affected by a recessive condition
def determine_disease_state(genotype: str, dominance: str) -> str:
    if dominance == "recessive":
        if genotype == "heterozygous":
            return "carrier"
        elif genotype == "homozygous alternate":
            return "affected"
        else:
            return "unaffected"
    elif dominance == "dominant":
        if genotype == "heterozygous" or genotype == "homozygous alternate":
            return "affected"
        else:
            return "unaffected"
    elif dominance == "unknown":
        if genotype == "homozygous alternate":
            return "affected (unknown dominance)"
        elif genotype == "heterozygous":
            return "potential carrier or affected (unknown dominance)"
        else:
            return "unaffected (unknown dominance)"
    else:
        return "NA"
    
###############################################################################
# Setting up command line arguments and flags using argparse.
###############################################################################

'''
argparse is used here to create a parser for this specific script.
    -p / --plink_prefix:    required, path to PLINK files (without extension)
    -c / --clinvar_file:    required, path to ClinVar_parsed.tsv
    -m / --modern_samples:  required, path to modern_samples.txt
    -o / --output_file:     optional, output TSV file name (default: ClinVar_AADR_genotypes.tsv)
    -d / --output_dir:      optional, output directory (default: current directory)
'''

parser = argparse.ArgumentParser(prog='AADR_parsing.py', description='Merging ClinVar pathogenic variant data with ancient DNA genotype data from the AADR dataset.')
parser.add_argument('-p', '--plink_prefix', required=True, help='Path to PLINK files without extension (e.g. ../resources/AADR/v54.1_1240K_public)')
parser.add_argument('-c', '--clinvar_file', required=True, help='Path to ClinVar_parsed.tsv file.')
parser.add_argument('-m', '--modern_samples', required=True, help='Path to modern_samples.txt file.')
parser.add_argument('-o', '--output_file', nargs='?', default='ClinVar_AADR_genotypes.tsv', help='Output file name. Default is ClinVar_AADR_genotypes.tsv')
parser.add_argument('-d', '--output_dir', required=False, default='.', help='Output directory. Default is current directory.')
args = parser.parse_args()

# check that all input files exist
for path, name in [(args.clinvar_file,'ClinVar file'),(args.modern_samples,  'Modern samples file'),]:
    input_path = Path(path)
    try:
        if not input_path.exists():
            raise Exception(f"- Error: The provided {name} does not exist or is missing: {path}\nExiting...")
    except Exception as e:
        print(e)
        sys.exit()

# check plink files exist (.bed, .bim, .fam)
for ext in ['.bed', '.bim', '.fam']:
    plink_path = Path(args.plink_prefix + ext)
    try:
        if not plink_path.exists():
            raise Exception(f"- Error: PLINK file missing: {plink_path}\nExiting...")
    except Exception as e:
        print(e)
        sys.exit()

output_path = os.path.join(args.output_dir, args.output_file)

#%%
###############################################################################
# Loading in AADR data and parsed clinvar.
###############################################################################

# load the AADR data
G, bim, fam = load_plink_files(args.plink_prefix)

# load clinvar data
clinvar = load_clinvar(args.clinvar_file)

# keeping only rows where the reference and alternate alleles are different (i.e. true variants) to speed up matching and avoid false matches
clinvar = clinvar[clinvar["ReferenceAlleleVCF"] != clinvar["AlternateAlleleVCF"]]


#%%
###############################################################################
# Prepare bim file for merge.
###############################################################################

# replace the rs in the rsid column with an empty string to get just the numeric part of the rs number
# creates a new column in the bim dataframe called rs_num that contains just the numeric part of the rs number
bim["rs_num"] = bim["rsid"].str.replace("rs", "", regex=False)

# replace allele info with actual A C T or G 
allele_map = {"1": "A", "2": "C", "3": "G", "4": "T"}
bim["ref_allele_nuc"] = bim["ref_allele"].map(allele_map)
bim["alt_allele_nuc"] = bim["alt_allele"].map(allele_map)

# forcing all to strings
bim["rs_num"] = bim["rs_num"].astype(str)
bim["ref_allele_nuc"] = bim["ref_allele_nuc"].astype(str)
bim["alt_allele_nuc"] = bim["alt_allele_nuc"].astype(str)

# this will ignore any na, N, R, and Y allele values in the alt and ref allele columns, as these are not valid alleles and would not match with the ancient DNA data in the bim dataframe

#%%
###############################################################################
# Load into SQLite db and match SNPs.
###############################################################################

# load into SQlite database
# create an in-memory SQLite database and load the bim and clinvar dataframes into it as tables named 'bim' and 'clinvar', respectively
connection = sqlite3.connect(':memory:')
bim.to_sql('bim', connection, index=False, if_exists='replace') # index true because my bim file does not have indexes
clinvar.to_sql('clinvar', connection, index=False, if_exists='replace')

# query the database to find matches between the bim and clinvar tables based on rs ID and allele info
sql_match = pd.read_sql_query("""
    SELECT
        b.i,
        b.rsid,
        b.chromosome,
        b.position,
        b.ref_allele_nuc,
        b.alt_allele_nuc,
        c.GeneSymbol,
        c.ClinicalSignificance,
        c."RS# (dbSNP)",
        c.PhenotypeList,
        c.Dominance
    FROM bim b
    INNER JOIN clinvar c
        ON b.rs_num = c."RS# (dbSNP)" 
        AND b.alt_allele_nuc = c.AlternateAlleleVCF
        AND b.ref_allele_nuc = c.ReferenceAlleleVCF
""", connection)

print(f"--Check: Matched {len(sql_match)} SNPs")


#%%
###############################################################################
# Extract matched genotypes.
###############################################################################

# extract the rows from the G matrix (SNPs) that match what is specified in the sql_match dataframe (which contains the matched SNPs based on rs ID and allele info) and store those rows in a new variable called G_matched
G_matched = G[sql_match['i'].values, :].compute()


#%%
###############################################################################
# Build the output dataframe (vectorized for speed).
###############################################################################

# create a new dataframe called G_df from the G_matched variable and add SNP metadata columns from the sql_match dataframe
G_df = pd.DataFrame(G_matched, columns=fam["individual_id"].values)

# add SNP metadata columns
for col in ['i', 'rsid', 'chromosome', 'position', 'ref_allele_nuc',
            'alt_allele_nuc', 'GeneSymbol', 'ClinicalSignificance',
            'PhenotypeList', 'Dominance']:
    G_df[col] = sql_match[col].values

# melt from wide to long (one row per SNP x individual)
snp_cols = ['i', 'rsid', 'chromosome', 'position', 'ref_allele_nuc',
            'alt_allele_nuc', 'GeneSymbol', 'ClinicalSignificance',
            'PhenotypeList', 'Dominance']

sql_match_genotypes = G_df.melt(
    id_vars=snp_cols,
    var_name='individual_id',
    value_name='dosage'
)

# drop missing genotypes
sql_match_genotypes = sql_match_genotypes.dropna(subset=['dosage'])
sql_match_genotypes['dosage'] = sql_match_genotypes['dosage'].astype(int)

# filter out modern samples using second column of modern_samples.txt
modern_samples = pd.read_csv(args.modern_samples, sep='\t', header=None)
modern_sample_ids = set(modern_samples.iloc[:, 1].values)
sql_match_genotypes = sql_match_genotypes[~sql_match_genotypes['individual_id'].isin(modern_sample_ids)]
print(f"--Check: Rows after filtering modern samples: {len(sql_match_genotypes)}") 

# apply zygosity and carrier status vectorized using existing functions
sql_match_genotypes['genotype'] = sql_match_genotypes['dosage'].map({
    0: 'homozygous reference',
    1: 'heterozygous',
    2: 'homozygous alternate'
})
sql_match_genotypes['disease_state'] = sql_match_genotypes.apply(
    lambda r: determine_disease_state(r['genotype'], r['Dominance']), axis=1)


#%%
###############################################################################
# Saving output.
###############################################################################

sql_match_genotypes.to_csv(output_path, sep='\t', index=False)
sql_match_genotypes.to_sql('clinvar_aadr_genotypes', connection, index=False, if_exists='replace')

print(f"--Done: Saved {len(sql_match_genotypes)} rows to {output_path}")


