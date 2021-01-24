#Trims the strand-seq index adapters that we use from both reads (set up for PE data)
#Trims low-quality bases from the 3' ends and removes very short reads
#Structure borrowed from Adam's alignment pipeline
#Vincent Hanlon 2019

##################################################################################################################################################


#process the file name
base=$(basename $1 "1.fastq.gz")

#assign input and output directories
dir="/home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/fastq/" 
dirout="/home/zhamadeh/NEW_PAPER_OLD_DATA/pe_alignment_pipeline/trimmed"	


cutadapt \
	-b file:nieks_adapters.fasta \
	-B file:nieks_adapters.fasta \
	-o $dirout/$base.trimmed.1.fastq.gz \
	-p $dirout/$base.trimmed.2.fastq.gz \
	$dir/$base"1.fastq.gz" \
	$dir/$base"2.fastq.gz" \
	-m 30 \
	--nextseq-trim 15
	#-q 15 
	#For the Nextseq this should be --nextseq-trim 15 
       	#For the Miseq this should be -q 15

