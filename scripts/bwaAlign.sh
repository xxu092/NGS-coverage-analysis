#!/usr/bin/bash -l 
#SBATCH -p epyc -n 1 -N 1 -c 16 --mem 96gb --out logs/bwa%A.log

module load bwa 

bwa index dataset/GCF_000001405.26_GRCh38_genomic.fna.gz

bwa mem -t 16 dataset/GCF_000001405.26_GRCh38_genomic.fna.gz dataset/SRR622461.fastq.gz > SRR622461.sam 
