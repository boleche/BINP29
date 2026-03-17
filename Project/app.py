#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
app.py

Description: Streamlit app script. Run this to open the app locally.

User-defined functions:  
    1. from user_parser import parse_user_file
    2. from clinvar_user_match import load_clinvar, match_user_to_clinvar
    3. from compare import compare_aadr_to_user

Non-standard modules: 
    1. sys
    2. streamlit
    3. pandas
    4. numpy
    5. pathlib
    
    
Procedure:
    1. 


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
import io


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
SAMPLES_DIR = Path(__file__).parent / "samples"
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

# Sample files available for demo use
SAMPLE_FILES = {"23andMe v2 (test sample)": "23andme_v2_test.txt", "Ancestry v2 (test sample)": "Ancestry.com_V2.txt"}

# Define _SampleFile once here so it's available everywhere below
class _SampleFile(io.BytesIO):
    """Wraps a bytes buffer so it behaves like a Streamlit UploadedFile."""
    def __init__(self, data, name):
        super().__init__(data)
        self.name = name
    def getbuffer(self):
        return memoryview(self.getvalue())

sample_file_obj = None


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


    st.divider()

    # adding in the sample file section
    st.header("🧪 Try a Sample File")
    st.markdown(
        "Don't have your own genetic data file? "
        "Select a sample below to load it directly into the app."
    )
    
    # adding a dropdown menu
    sample_choice = st.selectbox(
        "Choose a sample file",
        options=["— None —"] + list(SAMPLE_FILES.keys()),
        index=0,
    )
    
    
    if sample_choice != "— None —":
        sample_filename = SAMPLE_FILES[sample_choice]
        sample_path     = SAMPLES_DIR / sample_filename

        if sample_path.exists():
            # Preview expander (unchanged)
            with st.expander("Preview sample file"):
                try:
                    preview_df = pd.read_csv(
                        sample_path,
                        sep     = "\t",
                        comment = "#",
                        header  = 0,
                        nrows   = 20,
                        dtype   = str,
                    )
                    st.dataframe(preview_df, use_container_width=True)
                except Exception:
                    raw_lines = sample_path.read_text(errors="replace").splitlines()
                    visible   = [l for l in raw_lines if not l.startswith("#")][:20]
                    st.code("\n".join(visible), language="text")
                    
            # Button to confirm loading
            if st.button("▶ Load this sample", type="primary"):
                st.session_state["sample_confirmed"] = sample_filename
        else:
            st.error(...)

    # After the sidebar block, resolve sample_file_obj from session state
    sample_file_obj = None
    if "sample_confirmed" in st.session_state:
        confirmed_path = SAMPLES_DIR / st.session_state["sample_confirmed"]
        if confirmed_path.exists():
            sample_bytes    = confirmed_path.read_bytes()
            sample_file_obj = _SampleFile(sample_bytes, st.session_state["sample_confirmed"])
            st.sidebar.success(f"Sample active: **{st.session_state['sample_confirmed']}**")

# file input section for user to upload their genetic data file
st.header("Step 1: Upload Your Genetic Data")

# adding in a file uploader for the user to upload their genetic data file
uploaded_file = st.file_uploader("Choose a file", type=["txt", "csv", "tsv"], help="Accepted formats: 23andMe raw data, AncestryDNA raw data, .txt, .csv files.")

if uploaded_file and "sample_confirmed" in st.session_state:
    del st.session_state["sample_confirmed"]

active_file = uploaded_file or sample_file_obj

if not active_file:
    st.info("Waiting for file upload, or select a sample file from the sidebar.")
    st.stop()

#%%
###############################################################################
# Parsing the user data file.
###############################################################################

st.header("Step 2: Standardizing Your Genetic Data")

@st.cache_data(show_spinner="Parsing your genetic data...")

# function to run the existing parser functions on the uploaded file and return a standardized dataframe
def run_parser(file_name, file_bytes):
    temp_path = SCRIPTS_DIR / "temp_user_file.txt"
    with open(temp_path, "wb") as f:
        f.write(file_bytes)

    parsed_df = parse_user_file(temp_path)
    temp_path.unlink()
    return parsed_df

# call on the parser
user_data = run_parser(active_file.name, active_file.getvalue())


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

# show preview of the parsed clinvar data used for match
with st.expander("See ClinVar data used for matching"):
    st.dataframe(clinvar_data.head(50))

# show preview of the AADR matched SNP data
with st.expander("See AADR ClinVar matches used for comparison"):
    st.dataframe(aadr_data.head(50))

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
    

    # show preview of user matched data
    with st.expander("See preview of user matches."):
        st.dataframe(clinvar_matches.head(50))
    
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

    # show preview of AADR matched data
    with st.expander("See preview of user matches."):
        st.dataframe(aadr_matches.head(50))

    st.download_button(
        label = "Download your ancient individual disease variant matches here.",
        data = aadr_matches.to_csv(sep="\t", index=False),
        file_name = "AADR_user_comparison.tsv"
    )

