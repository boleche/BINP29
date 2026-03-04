# Malaria Project

> 

## Summary of Contents:
1.
2.
3.
4.
5.
6.
7.

## 1. Downloading project materials from the server.

```bash
#!/usr/bin/bash
cp -r /home2/resources/binp29/Data/malaria .

```
> 


## 2. Gene Prediction for Plasmodium Vivax genome using GeneMark.
```bash
#!/usr/bin/bash

gmes_petap.pl --ES --cores 10 --sequence resources/plasmodium_genomes/Plasmodium_vivax.genome --work_dir gene_prediction/gene_prediction_vivax/

# move to the /tmp/PutGenomesHere

#####################################################

# inspecting other gene predictions / genomes and filling out info table ex.

# number of genes
cut -f3 knowlesi.gtf | grep "gene" | wc -l
# 4953

# size of genome
grep -v '^>' Plasmodium_cynomolgi.genome | tr -d '\n' | wc -c
# 26181343

# counting genomic GC ex.
grep -v '^>' Plasmodium_cynomolgi.genome | tr -d '\n' | awk '{s=$0; gc=gsub(/[GgCc]/,"",s); printf "%.6f\n", gc/length($0)*100}'
# 39.078305

```


## 3. Filtering the H. tartakovsky genome.
> Filtered with a GC content threshold of 30% per scaffold (based on the distribution of GC content between the bird and parasite genomes). Filtered at a minimum of 3000 nucleotides per scaffold based on the exercise.

```bash
#!/usr/bin/bash

python3 removeScaffold.py Haemoproteus_tartakovskyi.raw.genome 30 Ht_genome_filtered.genome 3000

# original contig #
grep "^>" Haemoproteus_tartakovskyi.raw.genome | wc -l
# 15048

# filtered contig #
grep "^>" Ht_genome_filtered.genome | wc -l
# 2222

# 12826 excluded

```

## 4. Gene Prediction for filtered H. tartakovsky using GeneMark.
> Using a minimum of 10,000 contigs for training as most genes are not smaller than this. 

```bash
#!/usr/bin/bash

# total contigs
grep "^>" Ht_genome_filtered.genome | wc -l
# 2222

# how many below 10,000
grep ">" Ht_genome_filtered.genome | grep -o "Length=[0-9]*" | cut -d'=' -f2 | awk '$1 < 10000' | wc -l
# 1771

grep ">" Ht_genome_filtered.genome | grep -o "Length=[0-9]*" | cut -d'=' -f2 | awk 'BEGIN{min=999999999; max=0} {if($1>max) max=$1;
 if($1<min) min=$1} END{print "Min:", min, "\nMax:", max}'
# Min: 3003 
# Max: 64494

# do the gene prediction 
gmes_petap.pl --ES --min_contig 10000 --cores 10 --sequence resources/plasmodium_genomes/Ht_genome_filtered.genome --work_dir gene_prediction/gene_prediction_tart/

```

## 5. Making a FASTA file using filtered H. tartakovsky gene prediction.

```bash
#!/usr/bin/bash

# need to clean up first column of our gtf first (remove GC and length info after contig #)
cat ../Ht.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht2.gtf

# use the parser to make a fasta file (-c and -p options used to output an amino acid fasta as well)
gffParse.pl -i ../../resources/plasmodium_genomes/Ht_genome_filered -g ../Ht.gtf -b Ht_fasta -c -p

# INFO: The gff or gtf file ../Ht2.gtf has successfully been parsed.
#       There were 2120 scaffolds containing genes.
#       The scaffolds contained 3683 genes.
#       The genes contained the feature CDS 8289 times.

# INFO: The scaffold/genome file ../../resources/Ht_genome_filered was successfully parsed.
#       There were 2222 scaffolds. This number may be higher than the one above.

```

## 6. Removing genes of avian origin using BLASTx and BLASTp. 

```bash
#!/usr/bin/bash

# confirming blast pathway
echo $BLASTDB
# /home2/resources/blastdb

# running blastx on the generated fna file
nohup blastx -db SwissProt -query ../gene_prediction/gene_prediction_tart/fasta/Ht_fasta.fna -out Ht.blastx -num_threads 10

# running blastp on the generated faa file 
nohup blastp -db SwissProt -query ../gene_prediction/gene_prediction_tart/fasta/Ht_fasta.faa -out Ht.blastp -num_threads 10

# confirming all queries were run
grep "Query=" Ht.blastp | wc -l
# 3683
grep "Query=" Ht.blastx | wc -l
# 3683

######################################################

# soft linking uniprot tax file 
ln -s /resources/binp29/Data/malaria/taxonomy.dat taxonomy.dat
ln -s /resources/binp29/Data/malaria/uniprot_sprot.dat uniprot_sprot.dat

# running parser script which prints out host scaffolds (from the blastx output and within the fna fasta file)
python3 scripts/datParser.py host_scaffold_cleaning/blast_ht/Ht.blastx gene_prediction/gene_prediction_tart/fasta/
Ht_fasta.fna host_scaffold_cleaning/blast_ht/taxonomy.dat host_scaffold_cleaning/blast_ht/uniprot_sprot.dat > host_scaffold_cleaning/blastx_scaffolds.txt

# running parser script which prints out host scaffolds (from the blastp output and within the faa fasta file)
python3 scripts/datParser.py host_scaffold_cleaning/blast_ht/Ht.blastp gene_prediction/gene_prediction_tart/fasta/
Ht_fasta.faa host_scaffold_cleaning/blast_ht/taxonomy.dat host_scaffold_cleaning/blast_ht/uniprot_sprot.dat > host_scaffold_cleaning/blastp_scaffolds.txt

######################################################

# comparing the host scaffolds from blastx output and blastp output
wc -l blastp_scaffolds.txt 
# 23 
wc -l blastx_scaffolds.txt 
# 32
comm -12 blastx_scaffolds.txt blastp_scaffolds.txt | wc -l
# 16 shared scaffolds between the two 

# scaffolds only in blastp output
comm -23 blastp_scaffolds.txt blastx_scaffolds.txt 
# contig00021
# contig00142
# contig00338
# contig00352
# contig00356
# contig00466
# contig01393

# scaffolds only in blastx output
comm -13 blastp_scaffolds.txt blastx_scaffolds.txt 
# contig00042
# contig00044
# contig00080
# contig00232
# contig00235
# contig00238
# contig00239
# contig00621
# contig00677
# contig00970
# contig01030
# contig01132
# contig01529
# contig01562
# contig01643
# contig01899

# because blastx will blast against many different reading frames i will take a strict filtering approach and filter based on all of the blastx host identified scaffolds and all of the remaining blastp host identified scaffolds

######################################################

# calling on the remove_host_scaffolds.py script to remove the id'ed conitgs
python3 scripts/remove_host_scaffolds.py -f resources/plasmodium_genomes/Ht_genome_filtered.genome -i host_scaffold_cleaning/combined_scaffolds.txt -o resources/plasmodium_genomes/Ht_malaria_only.genome

# confirm final contig #
grep "^>" Ht_malaria_only.genome | wc -l
# 2183

```

## 7. Re-doing gene prediction on host-filtered Ht genome. 

```bash
#!/usr/bin/bash

nohup gmes_petap.pl --ES --min_contig 10000 --cores 10 --sequence ../../resources/plasmodium_genomes/Ht_malaria_only.ge
nome --work_dir .

```

## 8. Creating .faa amino acid fasta files for each genome gtf file.

```bash
#!/usr/bin/bash

# ht
cat ht_nohost.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht_nohost2.gtf
gffParse.pl -i ../../resources/plasmodium_genomes/Ht_malaria_only.genome -g Ht_nohost2.gtf -b fasta/Ht -c -p

# repeated for all other species

# copying the relevant .faa files to the orthos directory

# should have
ls
# Ht.faa  Pb.faa  Pc.faa  Pf.faa  Pk.faa  Pv.faa  Py.faa  Tg.faa

```

## 9. Finding orthologs with proteinortho and BUSCO.

```bash
#!/usr/bin/bash

# setting up a conda env for necessary packages
conda create -n busco_env -c conda-forge -c bioconda python=3.10 busco -y
conda install -c bioconda proteinortho

# running protein ortho 
nohup proteinortho6.pl -cpus=10 {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa > proteinortho.log 2>&1 &

# for proteinortho you must remove stop *
for f in {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa; do
  sed -i -E '/^>/! s/[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]//g; /^$/d' "$f"
done

```