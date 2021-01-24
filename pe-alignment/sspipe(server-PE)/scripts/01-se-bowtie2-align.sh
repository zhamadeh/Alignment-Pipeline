#!/bin/bash

for sample in `ls /home/zhamadeh/for_ERIBA/fastq/*.fastq.gz`
do
dir="/home/zhamadeh/for_ERIBA/fastq/"
base=$(basename $sample ".fastq.gz")
out="/home/zhamadeh/for_ERIBA/sam/"
bowtie2 -x './refseq/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index' -U ${dir}/${base}.fastq.gz -S ${out}/${base}.sam
done
