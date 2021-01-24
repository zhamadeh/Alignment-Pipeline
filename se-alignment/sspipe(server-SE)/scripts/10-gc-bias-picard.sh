#Extract and plot metrics about GC content

########################################################################################################################################################

#!/bin/bash

#settings directories and filenames
dir="/home/zhamadeh/for_ERIBA/mdup/"
out="./metrics/CollectGCBiasMetrics/"
out2="./figs/CollectGCBiasMetrics/"
base=$(basename $1 ".mdup.bam")

java -jar '/home/zhamadeh/Packages/Picard/picard-2.20.5/picard.jar' CollectGcBiasMetrics \
	I=${dir}/${base}.mdup.bam \
	O=${out}/${base}.gc_bias_metrics.txt \
	CHART=${out2}/${base}.gc_bias_metrics.pdf \
	S=${out}/${base}.summary_metrics.txt \
	R=./refseq/GRCh38.fasta
