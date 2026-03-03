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

gmes_petap.pl --ES --cores 10 --sequence resources/plasmodium_genomes/Plasmodium_vivax.genome --work_dir gene_prediction_vivax/

# move to the /tmp/PutGenomesHere

```


## 3. Filtering the H. tartakovsky genome.
> Filtered with a GC content threshold of 30% per scaffold (based on the distribution of GC content between the bird and parasite genomes). Filtered at a minimum of 3000 nucleotides per scaffold based on the exercise.

```bash
#!/usr/bin/bash

python3 removeScaffold.py Haemoproteus_tartakovskyi.raw.genome 30 Ht_genome_filered 3000

# original contig #
grep "^>" Haemoproteus_tartakovskyi.raw.genome | wc -l
# 15048

# filtered contig #
grep "^>" Ht_genome_filered | wc -l
# 2222

# 12826 excluded

```

## 4. Gene Prediction for filtered H. tartakovsky using GeneMark.
> Using a minimum of 10,000 contigs for training as most genes are not smaller than this. 


```bash
#!/usr/bin/bash

# total contigs
grep "^>" Ht_genome_filered | wc -l
# 2222

# how many below 10,000
grep ">" Ht_genome_filered | grep -o "Length=[0-9]*" | cut -d'=' -f2 | awk '$1 < 10000' | wc -l
# 1771

grep ">" Ht_genome_filered | grep -o "Length=[0-9]*" | cut -d'=' -f2 | awk 'BEGIN{min=999999999; max=0} {if($1>max) max=$1;
 if($1<min) min=$1} END{print "Min:", min, "\nMax:", max}'
# Min: 3003 
# Max: 64494

# do the gene prediction 
gmes_petap.pl --ES --min_contig 10000 --cores 10 --sequence resources/Ht_genome_filered --work_dir gene_prediction_tart/

```

## 5. Making a FASTA file using filtered H. tartakovsky gene prediction.

```bash
#!/usr/bin/bash

# need to clean up first column of our gtf first (remove GC and length info after contig #)
cat ../Ht.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht2.gtf

# use the parser to make a fasta file (-c and -p options used to output an amino acid fasta as well)
gffParse.pl -i ../../resources/Ht_genome_filered -g ../Ht.gtf -b Ht_fasta -c -p

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
nohup blastx -db SwissProt -query ../gene_prediction_tart/fasta/Ht_fasta.fna -out Ht.blastx -num_threads 10

# running blastp on the generated faa file 
nohup blastp -db SwissProt -query ../gene_prediction_tart/fasta/Ht_fasta.faa -out Ht.blastp -num_threads 10

# confirming all queries were run
grep "Query=" Ht.blastp | wc -l
# 3683
grep "Query=" Ht.blastx | wc -l
# 3683




```