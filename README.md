# NGS-coverage-analysis
To calculate coverage of sequencing reads I used mosdepth.I took 2 approaches.\
1.Start from a BAM/CRAM file\
2.Start from fastq raw reads.**\
scripts I used are in `scripts/` folder. 

## Coverage analysis from BAM/CRAM file
### 1. First I downloaded cram file and its index file (mapped to hg38 reference genome) 
```
curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/NA12878.hg38.cram
curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/NA12878.hg38.cram.crai
```

### 2. Download reference genome and unzip
```
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
gunzip -k GCF_000001405.26_GRCh38_genomic.fna.gz
```

### 3. Index reference genome (needed for cram file)
```
module load samtools
samtools faidx GCF_000001405.26_GRCh38_genomic.fna
```

### 4.1 Run Mosdepth for the whole genome
```
module load mosdepth
mosdepth -t 16 -f dataset/GCF_000001405.26_GRCh38_genomic.fna NA12878 dataset/NA12878.hg38.cram
```
`NA12878.mosdepth.summary.txt` file includes mean, min, max coverage of each chromosome. For the whole genome we can see mean coverage is at 14.36
`NA12878.mosdepth.per-base.bed.gz` contains coverage at specific locations. Below is the snippet of the result
```
chr1    0       10000   0
chr1    10000   10001   11
chr1    10001   10002   15
chr1    10002   10004   16
chr1    10004   10006   17
chr1    10006   10010   18
chr1    10010   10016   19
chr1    10016   10018   20
chr1    10018   10019   21
chr1    10019   10021   22
chr1    10021   10022   20
chr1    10022   10024   23
chr1    10024   10026   22
chr1    10026   10027   23
chr1    10027   10028   21
chr1    10028   10029   24
chr1    10029   10031   22
chr1    10031   10036   24
```

### 4.2 Run Mosdepth for targetted regions.

I want to calculate coverage in gene regions of chromosome 1. so I downloaded genome annotation file(gff) for reference genome. 
```
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gff.gz

#convert gff.gz file to bed file

module load bedops
zcat GCF_000001405.26_GRCh38_genomic.gff.gz | gff2bed > GCF_000001405.26_GRCh38_genomic.bed

#only extract rows annotated as gene

awk '$8 == "gene"{print}' GCF_000001405.26_GRCh38_genomic.bed > GCF_000001405.26_GRCh38_genomic.genes.bed

#only extract chromosome 1 data 

grep "NC_000001.11" GCF_000001405.26_GRCh38_genomic.genes.bed > h38_genomic.genes.chr1.bed

#change chromosome name to match cram file

awk 'BEGIN {OFS="\t"} $1 == "NC_000001.11" {$1 = "chr1"} {print}' h38_genomic.genes.chr1.bed > h38_genomic.genes1.chr1.bed
```
Run mosdepth with BED file
```
module load mosdepth
mosdepth -t 16 -n -b dataset/h38_genomic.genes1.chr1.bed -f dataset/GCF_000001405.26_GRCh38_genomic.fna NA12878chr1 dataset/NA12878.hg38.cram 
```
In `NA12878chr1.regions.bed.gz` file we can see the calculated coverage for gene regions.Below is snippet of the result.
```
chr1    10953   11523   .       32.76
chr1    11873   14409   .       36.79
chr1    14361   29370   .       40.90
chr1    17368   17436   .       50.69
chr1    30365   30503   .       26.23
chr1    34610   36081   .       25.70
chr1    69090   70008   .       3.99
chr1    120711  133748  .       15.71
chr1    134772  140566  .       8.50
chr1    142436  174392  .       5.35
chr1    184915  199860  .       21.14
chr1    285208  297520  .       6.13
chr1    359407  373177  .       11.02
chr1    450739  451678  .       1.03
chr1    490755  495445  .       12.27
chr1    501066  527676  .       25.15
chr1    627379  629009  .       23.09
chr1    632324  632413  .       15.81
chr1    685715  686654  .       14.00
chr1    725758  730351  .       28.70
chr1    735575  736753  .       34.66
chr1    738986  740800  .       31.08
chr1    764864  778688  .       19.46
chr1    785719  787127  .       16.96
chr1    817370  819834  .       12.37
```

### 5. Plot the whole genome coverage using `plot-dist.py` script from mosdepth package
```
python scripts/plot-dist.py NA12878.mosdepth.global.dist.txt
```

It generates a `dist.html` file where we can visualize coverage 


## Coverage analysis from fastq file 
When we need to generate BAM file from fastq file we can follow the following steps
1. Download fastq raw reads and reference genome
```
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461.fastq.gz
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
```
2. Index the reference genome and perform read mapping using bwa,
script I used is in scripts folder under the name `bwaAlign.sh`
```
module load bwa
bwa index GCF_000001405.26_GRCh38_genomic.fna.gz

bwa mem -t 16 dataset/GCF_000001405.26_GRCh38_genomic.fna.gz dataset/SRR622461.fastq.gz > SRR622461.sam 
```
3. convert SAM to BAM file and sort and index BAM file
```
module load samtools
samtools view -Sb SRR622461.sam > SRR622461.bam
samtools sort -@ 8 -o SRR622461.sorted.bam SRR622461.bam
samtools index SRR622461.sorted.bam
```
Run Mosdepth on targetted region Chr1 genes with script mosdepthBED2.sh

