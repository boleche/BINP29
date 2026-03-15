
# checking parsed clinvar dataset

# Usage: clinvar_check.sh ClinVar_parsed.tsv

PARSED=$1

echo "ClinVar Parsed .tsv QC"
echo "File: $PARSED"

# checking all columns necessary are included
echo "-- Checking column headers."
head -1 $PARSED | tr '\t' '\n'


# checking all variants are SNPs
echo "-- Checking for only SNPs."
cut -f1 $PARSED | tail -n +2 | sort | uniq -c

# checking clinical significant 
echo "-- Checking Clinical Significant (Pathogenic or Likely Pathogenic only."
cut -f3 $PARSED | tail -n +2 | sort | uniq -c

# checking for missing rsID values
echo "-- Checking for missing rsID values."
cut -f4 $PARSED | tail -n +2 | grep -cE "^(-|na|NA|)$"

# checking for only Chr37 assembly
echo "-- Checking Assembly."
cut -f6 $PARSED | tail -n +2 | sort | uniq -c

# checking refVCF 
echo "-- Checking Reference alleles."
cut -f14 $PARSED | tail -n +2 | sort | uniq -c

# checking altVCF 
echo "-- Checking Alternate alles."
cut -f15 $PARSED | tail -n +2 | sort | uniq -c 
