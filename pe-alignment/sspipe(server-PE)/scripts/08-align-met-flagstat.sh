#More basic alignment metrics. No need to parallelize here because it's so fast.
#Most of the other alignment metrics are for immediate post-alingment BAM files, but this flagstat script is used to createpost-alignment, post-filtering stats (really just read number). 
######################################################################################################################################

#!/bin/bash

printf "Library\tPostfiltering_reads_aligned\n" >/home/zhamadeh/for_ERIBA/sspipe/metrics/postfiltering_reads.txt

for sample in `ls /home/zhamadeh/for_ERIBA/mdup/*.bam`
do
	#assing input/output directories and process file name
	dir="/home/zhamadeh/for_ERIBA/mdup"
	out="./metrics/flagstat/"
	base=$(basename $sample ".mdup.bam")

	samtools flagstat ${dir}/${base}.mdup.bam > ${out}/${base}.flagstat.txt
		
	#extract post-filtering number reads
	num=$(head -5 ${out}/${base}.flagstat.txt | tail -n 1 | cut -f1 -d+)
	name=$(echo $base | sed 's/\.trimmed//g')
	printf "$name\t$num\n" >> /home/zhamadeh/for_ERIBA/sspipe/metrics/postfiltering_reads.txt	

done




