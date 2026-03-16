#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
app.py

Description: Streamlit app script. Run this to open the app locally.

User-defined functions:  
    1. 

Non-standard modules: 
    1. 
    
Procedure:
    1. 

Input: 
Output:

Usage: streamlit run app.py
Name: Emma Bolech
Date: 2026-03-14

"""
#%%
###############################################################################
# Importing necessary non-standard modules.
###############################################################################

import sys
import streamlit as st
import pandas as pd
import numpy as np
import pathlib
from pathlib import Path



#%%
###############################################################################
# App setup.
###############################################################################

st.set_page_config(
    page_title="GenMore",
    page_icon="🧬",
    layout="centered",
    initial_sidebar_state="expanded"
)


st.markdown("""
    <h1 style="font-family: 'Georgia', serif; 
               font-size: 2.5em; 
               font-weight: bold; 
               letter-spacing: 2px;
               color: #2E7D8C;">
        🧬 GenMore
    </h1>
    <h2 style="font-family: 'Georgia', serif;
               font-weight: normal;
               color: #444;">
        Genetic Risk Assessment and Ancient DNA Disease Susceptibility Explorer
    </h2>
    <p style="font-family: 'Georgia', serif; font-size: 1.1em; color: #666;">
        Upload your genetic data to explore your disease-associated variants against ClinVar variants. 
        Compare your results with ancient DNA samples from the AADR dataset to see how your genetic 
        risk compares with ancient populations.
    </p>
    <hr>
""", unsafe_allow_html=True)


#%%
###############################################################################
# Loading in reference data.
###############################################################################

# define constants for data loading
SCRIPTS_DIR = Path(__file__).parent / "scripts" 
DATA_DIR    = Path(__file__).parent / "results"
# add scripts directory to sys.path for importing functions from scripts
sys.path.insert(0, str(SCRIPTS_DIR))


# importing functions from the following scripts
from user_parser import parse_user_file
from clinvar_user_match import load_clinvar, match_user_to_clinvar
from compare import compare_aadr_to_user


# caching ref data loading for performance
@st.cache_resource(show_spinner="Loading ClinVar and AADR reference data...")

def load_data():
    # loading parsed clinvar tsv and parsed AADR tsv
    clinvar_path = DATA_DIR / "ClinVar_parsed.tsv"
    aadr_path = DATA_DIR / "AADR_clinvar_matches.tsv"

    # check for existence of data files
    for path, name in [(clinvar_path, "ClinVar data"), (aadr_path, "AADR data")]:
        try:
            if not path.exists():
                raise Exception(f"- Error: {name} file does not exist: {path}\nExiting...")
        except Exception as e:
            st.error(str(e))
            sys.exit()


    clinvar = load_clinvar(clinvar_path)
    # adding in filtration step for where ref allele = alt allele
    clinvar = clinvar[clinvar["ReferenceAlleleVCF"] != clinvar["AlternateAlleleVCF"]]

    aadr = pd.read_csv(aadr_path, sep="\t", dtype=str)

    return clinvar, aadr

# call on the function to load in the two reference datasets (clinvar and aadr)
clinvar_data, aadr_data = load_data()

#%%
###############################################################################
# Loading in user data.
###############################################################################

# adding in an information sidebar for types of user data accepted
with st.sidebar:
    st.header("Accepted Data Formats")
    st.markdown(
        "- **23andMe**: raw data file (txt format) downloaded from your 23andMe account.\n"
        "- **AncestryDNA**: raw data file (txt format) downloaded from your AncestryDNA account.\n"
        "- **txt/csv/tsv**: file containing your genetic variants.\n"
        "\n\n"
        "Please ensure that your data file is in one of the accepted formats and contains the necessary information for analysis (e.g., rsIDs, genotypes)."
    )   


# file input section for user to upload their genetic data file
st.header("Step 1: Upload Your Genetic Data")

# adding in a file uploader for the user to upload their genetic data file
uploaded_file = st.file_uploader("Choose a file", type=["txt", "csv", "tsv"], help="Accepted formats: 23andMe raw data, AncestryDNA raw data, .txt, .csv files.")

# if no file is uploaded, display a message and stop the app from running further
if not uploaded_file:
    st.info("Waiting for file upload...")
    st.stop()


#%%
###############################################################################
# Parsing the user data file.
###############################################################################

st.header("Step 2: Standardizing Your Genetic Data")

@st.cache_data(show_spinner="Parsing your genetic data...")

# function to run the existing parser functions on the uploaded file and return a standardized dataframe
def run_parser(uploaded_file):
    # read the uploaded file into a temporary location so it can be parsed by the existing parser functions
    temp_path = SCRIPTS_DIR / "temp_user_file.txt"
    with open(temp_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    # call the parser function to parse the user file and return a standardized dataframe
    parsed_df = parse_user_file(temp_path)

    # remove the temporary file after parsing
    temp_path.unlink()

    return parsed_df

# call on the parser
user_data = run_parser(uploaded_file)


# quick validation of the parsed user data to check for required columns and display a message if any are missing
required_cols = {"rsid", "chromosome", "position", "allele1", "allele2", "zygosity"}
missing_cols = required_cols - set(user_data.columns)
if missing_cols:
    st.error(f"Error: Parsed user data is missing required columns: {missing_cols}\nPlease check your input file and try again.")
    st.stop()

st.success("Your data has been successfully parsed!")

# show preview of the parsed user data in an expander
with st.expander("See parsed data"):
    st.dataframe(user_data.head(50))


#%%
###############################################################################
# Performing ClinVar and AADR Matching.
###############################################################################

st.header("Step 3: Identifying Disease Variants and Any Ancient DNA Matches")

@st.cache_resource(show_spinner="Matching your variants against ClinVar...")

def clinvar_match(user_data, clinvar_data):
    user_data = pd.read_json(user_data, orient="split")
    return match_user_to_clinvar(user_data, clinvar_data)

@st.cache_resource(show_spinner="Matching your variants against ancient individuals...")

def aadr_match(clinvar_matches_json, aadr_data):
    aadr_matches = pd.read_json(clinvar_matches_json, orient="split")
    return compare_aadr_to_user(aadr_df=aadr_data, user_df=clinvar_matches)
    

user_data_json = user_data.to_json(orient="split")
clinvar_matches = clinvar_match(user_data_json, clinvar_data)
clinvar_matches_json = clinvar_matches.to_json(orient="split")
aadr_matches = aadr_match(clinvar_matches_json, aadr_data)

#%%
###############################################################################
# Showing results.
###############################################################################

st.header("Step 4: Results")

tab1, tab2 = st.tabs(["Your Disease Susceptibility", "Ancient Individual Matches"])

###############################################################################

with tab1:
    if clinvar_matches.empty:
        st.warning("No ClinVar disease variants were found in your DNA file.")
    else:
        total     = len(clinvar_matches)
        affected  = (clinvar_matches["disease_state"].str.startswith("affected")).sum()
        carrier   = (clinvar_matches["disease_state"].str.startswith("carrier")).sum()
        unknown   = (clinvar_matches["disease_state"].str.contains("unknown")).sum()
        n_genes   = clinvar_matches["GeneSymbol"].nunique()

        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Total matches",  total)
        m2.metric("Affected",       affected)
        m3.metric("Carrier",        carrier)
        m4.metric("Unknown",        unknown)
        m5.metric("Genes involved", n_genes)

        st.divider()
    
    
    st.download_button(
        label = "Download your matched ClinVar disease variants table here.",
        data = clinvar_matches.to_csv(sep="\t", index=False),
        file_name = "user_clinvar_matches.tsv"
    )

###############################################################################

with tab2:
    if aadr_matches.empty:
        st.warning("No shared SNPs found between your matches and the AADR dataset.")
    else:
        n_individuals = aadr_matches["individual_id"].nunique()
        n_snps        = aadr_matches["rsid"].nunique()
        n_genes       = aadr_matches["GeneSymbol"].nunique()
        both_affected = (
            aadr_matches["aadr_disease_state"].str.startswith("affected") &
            aadr_matches["user_disease_state"].str.startswith("affected")
        ).sum()

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Shared SNPs",           n_snps)
        m2.metric("Ancient individuals",   n_individuals)
        m3.metric("Genes involved",        n_genes)
        m4.metric("Both you & AADR affected", both_affected)

        st.divider()

    st.download_button(
        label = "Download your ancient individual disease variant matches here.",
        data = aadr_matches.to_csv(sep="\t", index=False),
        file_name = "AADR_user_comparison.tsv"
    )

