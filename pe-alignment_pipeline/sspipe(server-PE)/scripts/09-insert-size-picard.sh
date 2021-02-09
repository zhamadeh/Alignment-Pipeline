#Extract and plot insert size distributions
#For PE reads only, of course
#########################################################################################################################

#!/bin/bash


#setting input/ouput and filenames
dir="/home/zhamadeh/for_ERIBA/mdup/"
out="./metrics/CollectInsertSizeMetrics/"
out2="./figs/CollectInsertSizeMetrics/"
base=$(basename $1 ".mdup.bam")

java -jar /home/zhamadeh/Packages/Picard/picard-2.20.5/picard.jar CollectInsertSizeMetrics \
I=${dir}/${base}.mdup.bam \
O=${out}/${base}.colinsert.txt \
H=${out2}/${base}.insert_size_histogram.pdf \
M=0.5

