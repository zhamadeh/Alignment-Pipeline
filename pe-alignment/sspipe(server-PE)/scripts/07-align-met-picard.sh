#Extracting basic alignment metrics

########################################################################################################################################################################

#!/bin/bash

#assign input/output directores and process file name. Note that this looks back to the relatively unprocessed ".sorted.bam" files produced immediately after alignment.
dir="/home/zhamadeh/for_ERIBA/sam/"
out="./metrics/CollectAlignmentMetrics/"
base=$(basename $1 ".sorted.bam")

java -jar /home/zhamadeh/Packages/Picard/picard-2.20.5/picard.jar CollectAlignmentSummaryMetrics \
R=./refseq/GRCh38.fasta \
I=${dir}/${base}.sorted.bam \
O=${out}/${base}.colalmet.txt

