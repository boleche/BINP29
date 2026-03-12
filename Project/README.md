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


python scripts/clinvar_parser.py -f resources/variant_summary.txt -d results/


```

## 3. Parsing AADR dataset.

https://www.ncbi.nlm.nih.gov/snp/rs33451#variant_details

checking SNV # rs33451 to verify reference assembly ---> 


```bash

grep "rs33451" v54.1_1240K_public.bim 
# 3       rs33451 0.643451        42401360        2       4

# matches ncbi position in GRCh37.p13 chr 3	NC_000003.11:g.42401360T>C

```