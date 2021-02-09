#Run GATK's MarkDuplicates and store the metrics somewhere useful

########################################################################################################################################################

#!/bin/bash

#Set input/output/metrics directories, and process the file name
dir="/home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/sorted/"
dirout="/home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/mdup/"
diroutmetrics="./metrics/mdup/"
base=$(basename $1 ".bam")



java -jar '/home/zhamadeh/Packages/Picard/picard-2.20.5/picard.jar' MarkDuplicates \
	I= ${dir}/${base}.sorted \
	O= ${dirout}/${base}.mdup.bam \
	M= ${diroutmetrics}/${base}.txt \
	VALIDATION_STRINGENCY=LENIENT
