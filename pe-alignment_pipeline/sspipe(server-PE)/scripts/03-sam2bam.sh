#Processing SAM files and sorting

#####################################################################################################################################

#!/bin/bash

#go to the directory where this script will write output
cd /home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/sorted

#assign input/output directories and process filenames
dir="/home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/sam/"
base=$(basename $1 ".sam")
out="/home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/sorted/"




#Convert from SAM to BAM and sort.
samtools view -bS ${dir}/${base}.sam | samtools sort - -o ${dir}/${base}.sorted.bam

#Index the bam file (this is done so that non-standard chromosomes can be removed, below)
samtools index ${dir}/${base}.sorted.bam

#Removing non-standard chromosomes from both the reads and the header. This process leaves some reads mis-paired, so we also need to run samtools fixmate, which requires special sorting. 
#We also require that reads be paired by setting the flag -f1, and we filter mapq < 10 with -q10
#CHANGE THE FLAG IF USING SINGLE END READS!!!

samtools view -q10 -h ${dir}/${base}.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY |
        grep -v -E '@SQ.*chrUn|@SQ.*chrM|@SQ.*random|@SQ.*chrEBV' | samtools sort -n - | samtools fixmate -O bam - - | samtools sort - | samtools view -bh -f1 -o ${out}/${base}.sorted






cd ~/sspipe

