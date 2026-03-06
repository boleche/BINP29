#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extract_busco_sequences.py

Description: Extract BUSCO protein sequences from multiple species and create
one fasta file per BUSCO id. BUSCOs with status "Complete" or "Duplicated"
are used. If duplicated, only the first gene is selected so that each output
file contains one sequence per organism.

Procedure:
    1. Define functions
    2. Establish working directory
    3. Load BUSCO tables and fasta files
    4. Find shared BUSCO ids across species
    5. Write one fasta file per BUSCO

Input: BUSCO full_table.tsv files and species protein fasta files
Output: one fasta file per BUSCO id

Version: 1.00
Date: 2026-03-05
Name: Emma Bolech
"""

#%%
###############################################################################
# Importing modules.
###############################################################################

from pathlib import Path
from collections import defaultdict


#%%
###############################################################################
# Defining functions.
###############################################################################

# read fasta file and return dictionary of gene -> sequence
def read_fasta(file):
    seqs = {}
    header = None
    seq = []
    with open(file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    gene = header.split()[0] # take first word as gene id
                    seqs[gene] = "".join(seq)   # join sequence lines
                header = line[1:]   # remove ">" from header
                seq = []     # reset sequence list for new gene
            else:
                seq.append(line)   # add sequence line to list
        # add last gene after loop
        if header:
            gene = header.split()[0]
            seqs[gene] = "".join(seq)

    return seqs #   return dictionary of gene -> sequence


# read BUSCO full_table.tsv
def read_busco_table(file):
    # return dictionary of busco_id -> list of gene_ids with status "Complete" or "Duplicated"
    buscos = defaultdict(list)
    with open(file) as f:
        for line in f:
            if line.startswith("#"): # skip header lines
                continue
            parts = line.strip().split("\t") # split line into parts by tab
            if len(parts) < 3: # skip lines that don't have at least 3 columns (busco_id, status, gene_id)
                continue
            # extract busco_id, status, and gene_id from the line
            busco_id = parts[0]
            status = parts[1]
            gene_id = parts[2]
            # if status is "Complete" or "Duplicated", add gene_id to the list for this busco_id
            if status in ("Complete", "Duplicated"):
                buscos[busco_id].append(gene_id)
    return buscos


#%%
###############################################################################
# Establish working directory.
###############################################################################

workingdirectory = Path.cwd()
print(f"working directory: {workingdirectory}")


#%%
###############################################################################
# Hard coded paths and species.
###############################################################################

busco_base = Path("/home/inf-25-2025/Desktop/binp29_git/BINP29/Malaria/orthos/busco")
fasta_base = Path("/home/inf-25-2025/Desktop/binp29_git/BINP29/Malaria/orthos/fasta")

output_dir = Path("busco_fasta_files_all")
output_dir.mkdir(exist_ok=True)

species = ["Ht","Pb","Pc","Pf","Pk","Pv","Py","Tg"]


#%%
###############################################################################
# Loading busco tables and fasta sequences.
###############################################################################

# opening busco tables and fasta files for each species and storing in dictionaries
busco_data = {}
sequence_data = {}

for sp in species:
    # construct file paths for busco table and fasta file for this species
    busco_file = busco_base / sp / "run_apicomplexa_odb12" / "full_table.tsv"
    fasta_file = fasta_base / f"{sp}.faa"
    # read busco table and fasta file and store in dictionaries
    busco_data[sp] = read_busco_table(busco_file)
    sequence_data[sp] = read_fasta(fasta_file)


#%%
###############################################################################
# Finding shared busco ids across all species.
###############################################################################

# start with the set of busco ids from the first species and then take the intersection with the sets from the other species to find shared busco ids
shared_buscos = set(busco_data[species[0]])

# loop through the remaining species and take the intersection of shared_buscos with the set of busco ids from each species
for sp in species[1:]:
    shared_buscos &= set(busco_data[sp])



#%%
###############################################################################
# Writing one fasta file per busco.
###############################################################################
# starting with a counter for the number of files created
created = 0
# loop through shared busco ids and write one fasta file per busco with one sequence per species
for busco in sorted(shared_buscos):
    # construct output file path for this busco
    outfile = output_dir / f"{busco}.faa"
    # open output file for writing
    with open(outfile, "w") as out:
        # loop through species and write the sequence for this busco if it exists
        for sp in species:
            # get the gene id for this busco and species from the busco_data dictionary
            gene = busco_data[sp][busco][0]
            seq = sequence_data[sp].get(gene)
            # if the sequence exists, write it to the output file in fasta format
            if seq:
                out.write(f">{sp}\n{seq}\n")
    # increment the counter for created files
    created += 1


#%%
###############################################################################
# Final message.
###############################################################################

print(f"created {created} files in {output_dir}")
print("success")

