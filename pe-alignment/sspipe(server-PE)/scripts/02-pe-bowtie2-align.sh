#Adam's bowtie2 alignment script for hg38. with additional parallelizing (for paired end reads---must be edited for SE reads)

######################################################################################################################

#!/bin/bash

#Loop over all files in the "trimmed" directory, which contains the cutadapt output
for sample in `ls /home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/trimmed/*.1.fastq.gz`
do
	#assign input/output directories and process the file name
	dir="/home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/trimmed/"
	base=$(basename $sample ".1.fastq.gz")
	out="/home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/sam/"
	

	bowtie2 -x 'refseq/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index' -p $1 -1 ${dir}/${base}.1.fastq.gz -2 ${dir}/${base}.2.fastq.gz -S ${out}/${base}.sam

done
