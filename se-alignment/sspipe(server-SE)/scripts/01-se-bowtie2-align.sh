#!/bin/bash

for sample in `ls /home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/trimmed/*.trimmed.fastq.gz`
do
dir="/home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/trimmed/"
base=$(basename $sample ".trimmed.fastq.gz")
out="/home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/sam/"
bowtie2 -x './refseq/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index' -U ${dir}/${base}.trimmed.fastq.gz -S ${out}/${base}.sam
done
