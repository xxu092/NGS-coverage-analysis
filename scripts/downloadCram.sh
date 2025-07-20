#!/usr/bin/bash -l
#SBATCH -p epyc -n 8 -N 1 --mem 64gb --out logs/downloadCram%A.log 
cd ../dataset
curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/NA12878.hg38.cram