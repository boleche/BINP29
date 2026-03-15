# Disease Susceptibility

> SQlite database (SQlite3 python package)


## Summary of Contents





## 1. Downloading datasets.

```

# downloading clinvar dataset

# this is updated weekly: my version was downloaded on March 10th, 2026

# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

# variant_summary.txt

# readme
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README

# downloading given ancient DNA datasets
1. AADR Annotations 2025.xlsx
2. Ancient_samples.txt
3. Modern_samples.txt
4. v54.1_1240K_public.bed
5. v54.1_1240K_public.bim
6. v54.1_1240K_public.fam

# downloading given test user datasets
1. Test1_DNA + txt
2. Test2_DNA + txt
3. Test3_DNA + txt
4. Test4_DNA + txt
5. Test5_DNA + txt


```

## 2. Parsing ClinVar dataset. 

```bash

# check column #
zcat variant_summary.txt.gz | head -1 | tr '\t' '\n' | wc -l
# 43

# call on python script for parsing
python scripts/clinvar_parser.py -f resources/variant_summary.txt -d results/

```
## 3. checking parsed clinvar dataset.


## 4. Parsing AADR dataset.

https://www.ncbi.nlm.nih.gov/snp/rs33451#variant_details

checking SNV # rs33451 to verify reference assembly ---> 

# Variants with ambiguous alternate alleles (IUPAC codes N, R, Y) 
# and missing alternate alleles (na) were excluded from analysis (n=34).


I HAD to switch bim a0 to alt and a1 to ref

The EIGENSTRAT .snp file format — columns 5 and 6 are reference and variant alleles respectively, so EIGENSTRAT lists ref first, alt second. Countingchromosomes
The PLINK .bim file format — PLINK .bim does the opposite where ALT is listed in column 5 and REF in column 6. The convertf tool is not aware of this PLINK .bim format, so after converting from EIGENSTRAT to PLINK, ALT is still listed in column 6 and REF in column 5 of the .bim file. ResearchGate This is exactly your situation.

```bash

conda create -n "binp29_project"

conda install pandas-plink
conda install pandas

python3 scripts/AADR_parsing.py -p resources/AADR/v54.1_1240K_public -c results/ClinVar_parsed.tsv -m resources/AADR/Modern_samples.txt -o AADR_clinvar_matches.tsv -d results/


grep "rs33451" v54.1_1240K_public.bim 
# 3       rs33451 0.643451        42401360        2       4

# matches ncbi position in GRCh37.p13 chr 3	NC_000003.11:g.42401360T>C


# peek at the 61 matched rsIDs to see what alleles look like on both sides
test4 = pd.read_sql_query("""
    SELECT b.rsid, b.ref_allele_nuc, b.alt_allele_nuc, 
           c.ReferenceAlleleVCF, c.AlternateAlleleVCF
    FROM bim b
    INNER JOIN clinvar c ON b.rs_num = c."RS# (dbSNP)"
""", connection)

pd.set_option('display.max_rows', None)
print(test4)
pd.reset_option('display.max_rows')

rsID + ref + alt matches:    count
0     31




```

## 5. Parsing user input files.


```bash

HAD TO remove snps where alt and ref are the same 

# parsing all to same format and including zygosity 
python3 scripts/user_parser.py -i resources/TestUsers/Test1_User.txt -o Test1_parsed.tsv -d results/users_parsed/

# repeat for all test users

# parsing all with clinvar data to output a table with disease state information
python3 scripts/clinvar_user_match.py -i results/users_parsed/Test1_parsed.tsv -c results/ClinVar_parsed.tsv -o Test1_cv_matches.tsv -d results/user_clinvar_matches/


python3 scripts/clinvar_user_match.py -i results/users_parsed/Test5_parsed.tsv -c results/ClinVar_parsed.tsv -o Test5_cv_matches.tsv -d results/user_clinvar_matches/


```

## 6. Matching User Input SNPs to AADR data SNPs.

```bash

python3 scripts/compare.py -i results/user_clinvar_matches/Test4_cv_matches.tsv -a results/AADR_clinvar_matches.tsv -o Test4_comparison.tsv -d results/comparison/

```


## 7. Streamlit app.

```bash

conda activate binp29_project
pip install streamlit


python3 scripts/compare.py -i results/user_clinvar_matches/Test4_cv_matches.tsv -a results/AADR_clinvar_matches.tsv -o Test4_comparison.tsv -d results/comparison/

```
