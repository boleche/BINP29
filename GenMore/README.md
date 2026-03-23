# GenMore Application

> GenMore is a locally hosted application that allows users to upload their DNA SNP information to recieve ClinVar pathogenic match results. Users can then compare their disease risk associated markers with ancient individuals. 

## Resource Information

> Original ClinVar data can be found at: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/\
> variant_summary.txt.gz : working version was built with the 2026-03-10 downloaded version.\
> README: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README
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

```


#### App Functionality Scripts

```
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

## Dependencies

> Dependencies are included in the envs/ folder in "plink_parsing_env.yml" and "requirements.txt".\
> Notable packages and versions:

1. pandas-plink v.2.3.2
2. python v.3.13.11
3. streamlit v.1.55.0


## Pre-parsed ClinVar dataset. 

> The ClinVar dataset was parsed prior to app running with the included script "clinvar_parser.py." The parsed output is included in the results/ folder.\
> The successful filtering was confirmed with clinvar_check.sh.


```bash


# call on python script for parsing
python scripts/clinvar_parser.py -f resources/variant_summary.txt -d results/

```
## Pre-parsed AADR dataset.

> The AADR dataset was parsed using pandas-plink and the python script "AADR_parsing.py".\
> The following files downloaded from the AADR database are needed to parse the ancient DNA.\
> The AADR original files were converted to PLINK friendly files so the alternate and reference allele order were swapped.\
> The parser handles this swap to ensure accurate alternate allele matching with ClinVar SNPs.

1. Modern_samples.txt - used to filter out modern individuals.
2. v54.1_1240K_public.bed - used for all SNP X individual information.
3. v54.1_1240K_public.bim - used to SNP and genotype information.
4. v54.1_1240K_public.fam - used for all individual ID index information.


```bash

# call on python script for parsing
python3 scripts/AADR_parsing.py -p resources/AADR/v54.1_1240K_public -c results/ClinVar_parsed.tsv -m resources/AADR/Modern_samples.txt -o AADR_clinvar_matches.tsv -d results/

```


## Running GenMore.

> After downloading this repository to your local computer, from the project root run the following command to open GenMore locally.\
> Ensure proper project structure as shown in the beginning of the README.\
> The pre-parsed ClinVar and AADR datasets as well as two sample inputs are included in the project directory so no\
> pre-filtering is needed as the user.

```bash

streamlit run app.py

```


