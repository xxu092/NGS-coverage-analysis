#!/usr/bin/bash -l
#SBATCH -p epyc -n 1 -N 1 -c 16 --mem 96gb --out logs/mosdepth%A.log

module load mosdepth
mosdepth -t 16 -n -b dataset/h38_genomic.genes1.chr1.bed -f dataset/GCF_000001405.26_GRCh38_genomic.fna NA12878chr1 dataset/NA12878.hg38.cram 