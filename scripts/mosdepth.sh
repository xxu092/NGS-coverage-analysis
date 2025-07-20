#!/usr/bin/bash -l
#SBATCH -p epyc -n 1 -N 1 -c 16 --mem 96gb --out logs/mosdepth%A.log

module load mosdepth
mosdepth -t 16 -f dataset/GCF_000001405.26_GRCh38_genomic.fna NA12878 dataset/NA12878.hg38.cram