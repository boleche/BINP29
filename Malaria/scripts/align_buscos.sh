#!/bin/bash

# Script to align all BUSCO fasta files using clustalo
# Usage: ./align_buscos.sh

# create output directory for aligned files
mkdir -p busco_aligned_all

# loop through each BUSCO fasta file
for input_file in ../busco_fasta_files_all/*.faa; do
    # extract filename without path and extension
    filename=$(basename "$input_file" .faa)
    # define output file
    output_file="busco_aligned_all/${filename}_aligned.faa"
    # run clustalo
    clustalo -i "$input_file" -o "$output_file" -v --threads=10
    
done

