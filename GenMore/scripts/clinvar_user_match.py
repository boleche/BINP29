#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
clinvar_user_match.py

Description: Matches parsed user SNP data against ClinVar pathogenic variants and determines disease state.

User-defined functions:  
    1. load_clinvar
    2. determine_disease_state

Non-standard modules: 
    1. sqlite3
    2. pandas
    3. argparse
    4. os
    5. sys
    6. pathlib
    
Procedure:
    1. Load ClinVar data.
    2. Load user SNP data.
    3. Match SNPs against ClinVar variants.
    4. Determine disease state for each matched variant.

Input: Parsed user SNP data, ClinVar pathogenic variants
Output: TSV file with disease state information

Usage: python3 clinvar_user_match.py -i <user_snp_file> -c <clinvar_file> [-o <output_file>] [-d <output_dir>]
Name: Emma Bolech
Date: 2026-03-13

"""
#%%
###############################################################################
# Importing necessary non-standard modules.
###############################################################################

import sqlite3
import pandas as pd
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

# function which adds to a column called 'disease_state' in the sql_match dataframe which indicates whether the individual is affected by a recessive condition
def determine_disease_state(zygosity: str, dominance: str, alt_allele_match: bool) -> str:
    zygosity = str(zygosity).strip().lower()
    dominance = str(dominance).strip().lower()
    
    if not alt_allele_match:
        return "unaffected (no alt allele match)"
    if dominance == "recessive":
        if zygosity == "heterozygous":
            return "carrier"
        elif zygosity == "homozygous":
            return "affected"
        else:
            return "unaffected"
    elif dominance == "dominant":
        if zygosity in ("heterozygous", "homozygous"):
            return "affected"
        else:
            return "unaffected"
    elif dominance == "unknown":
        if zygosity == "homozygous":
            return "affected (unknown dominance)"
        elif zygosity == "heterozygous":
            return "potential carrier or affected (unknown dominance)"
        else:
            return "unaffected (unknown dominance)"
    else:
        return "NA"
    

# function to match user SNPs to ClinVar variants using an in-memory SQLite database for efficiency
def match_user_to_clinvar(user_snps: pd.DataFrame, clinvar: pd.DataFrame) -> pd.DataFrame:
    # strip rs prefix for matching
    user_snps["rs_num"] = user_snps["rsid"].str.replace("rs", "", regex=False).astype(str)
    clinvar["RS# (dbSNP)"] = clinvar["RS# (dbSNP)"].astype(str)

    # load into in-memory SQLite and match
    connection = sqlite3.connect(":memory:")
    user_snps.to_sql("user_snps", connection, index=False, if_exists="replace")
    clinvar.to_sql("clinvar", connection, index=False, if_exists="replace")

    sql_match = pd.read_sql_query("""
        SELECT
            u.rsid,
            u.chromosome,
            u.position,
            u.allele1,
            u.allele2,
            u.zygosity,
            c.GeneSymbol,
            c.ClinicalSignificance,
            c."RS# (dbSNP)",
            c.PhenotypeList,
            c.Dominance,
            c.ReferenceAlleleVCF,
            c.AlternateAlleleVCF
        FROM user_snps u
        INNER JOIN clinvar c
            ON u.rs_num = c."RS# (dbSNP)"
            AND (u.allele1 = c.AlternateAlleleVCF OR
                 u.allele2 = c.AlternateAlleleVCF)
    """, connection)

    if sql_match.empty:
        return sql_match

    sql_match["disease_state"] = sql_match.apply(
        lambda r: determine_disease_state(r["zygosity"], r["Dominance"], True), axis=1
    )
    return sql_match


    
###############################################################################
# Setting up command line arguments and flags using argparse.
###############################################################################

'''
argparse is used here to create a parser for this specific script.
    -i / --input_file:      required, path to parsed user SNP file (output of user_parser.py)
    -c / --clinvar_file:    required, path to ClinVar_parsed.tsv
    -o / --output_file:     optional, output TSV file name (default: ClinVar_user_genotypes.tsv)
    -d / --output_dir:      optional, output directory (default: current directory)
'''


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="clinvar_user_match.py", description="Match parsed user SNP data against ClinVar pathogenic variants.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to parsed user SNP file (output of user_parser.py).")
    parser.add_argument("-c", "--clinvar_file", required=True, help="Path to ClinVar_parsed.tsv file.")
    parser.add_argument("-o", "--output_file", nargs="?", default="ClinVar_user_genotypes.tsv", help="Output file name. Default is ClinVar_user_genotypes.tsv")
    parser.add_argument("-d", "--output_dir", required=False, default=".", help="Output directory. Default is current directory.")
    args = parser.parse_args()

    for path, name in [(args.input_file, "Input file"), (args.clinvar_file, "ClinVar file")]:
        try:
            if not Path(path).exists():
                raise Exception(f"- Error: {name} does not exist: {path}\nExiting...")
        except Exception as e:
            print(e)
            sys.exit()

    output_path = os.path.join(args.output_dir, args.output_file)

    #%%
    ###############################################################################
    # Loading in User data and parsed clinvar.
    ###############################################################################

    user_snps = pd.read_csv(args.input_file, sep="\t", dtype=str)
    clinvar = load_clinvar(args.clinvar_file)

    # keeping only rows where the reference and alternate alleles are different (i.e. true variants) to speed up matching and avoid false matches
    clinvar = clinvar[clinvar["ReferenceAlleleVCF"] != clinvar["AlternateAlleleVCF"]]

    # validate required user columns
    required_cols = {"rsid", "chromosome", "position", "allele1", "allele2", "zygosity"}
    missing = required_cols - set(user_snps.columns)
    if missing:
        print(f"Error: User SNP file is missing required columns: {missing}\nExiting...")
        sys.exit(1)


    #%%
    ###############################################################################
    # Load into SQLite db and match SNPs.
    ###############################################################################

    # match user SNPs to ClinVar variants using an in-memory SQLite database for efficiency
    sql_match = match_user_to_clinvar(user_snps, clinvar)

    # if no matches were found, print a message and exit without writing an output file
    if sql_match.empty:
        print("No matches found. Output file will not be written.")
        sys.exit(0)



    #%%
    ###############################################################################
    # Saving output.
    ###############################################################################

    # save the output to a TSV file
    sql_match.to_csv(output_path, sep='\t', index=False)
    print(f"--Done: Saved {len(sql_match)} row(s) to {output_path}")

