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
```bash
#!/usr/bin/bash

python3 removeScaffold.py Haemoproteus_tartakovskyi.raw.genome 30 Ht_genome_filered 3000

```

## 4. 