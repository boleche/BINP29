# GenMore Application

> GenMore is a locally hosted application that allows users to upload their DNA SNP information to recieve ClinVar pathogenic match results. Users can then compare their disease risk associated markers with ancient individuals. 

## Resource Information

> Original ClinVar data can be found at: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/\
> variant_summary.txt.gz : working version was built with the 2026-03-10 downloaded version.
> 
> AADR data taken from the Reich Lab at: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW\
> v54.1_1240K_public\
> The data was originally downloaded and converted to PLINK format by Eran Elhaik using the following:\
> #Get the data\
> wget https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1/SHARE/public.dir/v54.1_1240K_public.tar\
> #Convertf\
> /home/eelhaik/Tools/EIG-5.0.2/bin/convertf -p /home/eelhaik/Tools/EIG-6.0.1/CONVERTF/par.EIGENSTRAT.PED.aDNA
>


## Project Structure Information




## Scripts Information

#### Pre-Filtering Scripts

```
1. clinvar_parser.py
    * Used to filter the variant_summary.txt file from ClinVar.
    * SNPs only
    * Pathogenic and Likely Pathogenic
    * Remove missing rsIDs
    * Keep only GRCh37 assembly (AADR dataset used GRCh37 assembly for SNPs)
    * Derived dominance information from phenotypic description

2. clinvar_check.sh
    * Used to check that the applied filters worked as expected.
    * Checks that all necessary columns are included.
    * Checks all variants are SNPs.
    * Checks only pathogenic and likely pathogenic SNPs are included.
    * Checks there are no missing rsID entries.
    * Checks for only Chr37 Assembly entries. 
    
3. AADR_parsing.py
    * Used to parse through the PLINK AADR files and match pathogenic ancient DNA SNPs with the parsed ClinVar data.
    * Ref allele cannot = alt allele.
    * Matches on rsID first then on alt allele.
    * Missing genotypes are dropped.
    * Carrier vs. affected info is ascertained from the phenotype list content (dominance is determined and disease status reflects zygosity at the SNP).
    * Modern individuals are excluded based on the included Modern_samples.txt individual names.
    * Non-exact alleles dropped (ex. R).



#### App Functionality Scripts
1. user_parser.py
    * Used to parse user input files into standardized input file format. 
    * Handles 7 different input file types.

2. clinvar_user_match.py
    * Used to match user SNPs to ClinVar filtered pathogenic SNPs.
    * Ref allele cannot = alt allele.
    * Matches on rsID first then on alt allele.
    * Carrier vs. affected info is ascertained from the phenotype list content (dominance is determined and disease status reflects zygosity at the SNP).

3. compare.py
    * Used to compare user output SNPs and ancient DNA AADR SNPs.
    * Load AADR_parsing output and clinvar_user_match output.
    * Merges the two datasets on rsid and individual_id.
    * Compares disease state determinations between the two datasets.
    * Outputs a TSV file with the comparison results.

4. app.py
    * Full script that calls on necessary functions for local app output.


```




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

#
G_sample = G[:1000, :].compute()

total = G_sample.size
n0 = np.sum(G_sample == 0)
n1 = np.sum(G_sample == 1)
n2 = np.sum(G_sample == 2)
nan = np.sum(np.isnan(G_sample))

print(f"0 (ref):  {n0} ({n0/total:.2%})")
print(f"1 (het):  {n1} ({n1/total:.2%})")
print(f"2 (alt):  {n2} ({n2/total:.2%})")
print(f"NaN:      {nan} ({nan/total:.2%})")


0 (ref):  1898020 (11.53%)
1 (het):  940995 (5.71%)
2 (alt):  8449309 (51.31%)
NaN:      5177676 (31.44%)


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




```
