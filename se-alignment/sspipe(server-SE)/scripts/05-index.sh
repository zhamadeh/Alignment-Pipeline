#Parallelized indexing of BAM files

######################################################################################################################################

#!/bin/bash

#Set output directory
dir="/home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/mdup/"

for sample in `ls /home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/mdup/*.bam`
do
	#process file name
	base=$(basename $sample ".bam")

	#parallelized indexing
	samtools index -@$1 ${dir}/${base}.bam 

done




