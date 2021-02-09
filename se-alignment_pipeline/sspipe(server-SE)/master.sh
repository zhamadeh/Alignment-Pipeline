#This master script and the scripts it calls are a mix of Adam's original pipeline and my custom scripts and alterations.
#Vincent Hanlon 2020

########################################################################################################################################################

#!/bin/bash


#Set number threads to use here:
p=20



#demultiplexing (ignore if starting with fastq)
#bcl2fastq --runfolder-dir /projects/lansdorp/for_ERIBA/200225_NB551308_0076_AH3HLTAFX2 -p $p --output-dir /projects/lansdorp/for_ERIBA/fastq --no-lane-splitting --barcode-mismatches 1

#moving some libraries around after demultiplexing (ignore if starting with fastq)
#mv /projects/lansdorp/for_ERIBA/fastq/pl*/* /projects/lansdorp/for_ERIBA/fastq


#for NextSeq concatenate 4 lanes into one fastq file
#cd /home/zhamadeh/for_ERIBA/uncat_fastq/
#for i in $(ls *gz | cut -f1-2 -d_ | sort -u); do cat "$i"*_R1_*gz > ../fastq/"$i"_R1_001.fastq.gz; cat "$i"*_R2_*gz > ../fastq/"$i"_R2_001.fastq.gz; done​

cd /home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/sspipe/


#breakpointR doesn't like empty libraries, so I'm creating a file to record which ones they are before I delete them.
#>../mdup/libraries_with_zero_reads.txt


parallel --citation

#aligning paired-end reads, processing, removing duplicates, and indexing
#ls /home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/fastq/*fastq.gz | parallel --jobs $p bash "./scripts/01-trim-adapters.sh" {}
./scripts/01-se-bowtie2-align.sh $p
ls /home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/sam/*.sam | parallel --jobs $p bash "./scripts/03-sam2bam.sh" {}




#I've arbitrarily reduced the parallelization for MarkDuplicates to 4-fold, because it seems to use more CPUs than I give it and slow down the server
ls /home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/sorted/*.sorted | rev | cut -f1 -d/ | rev | sed 's/\.sorted//g' | parallel --jobs 2 bash "./scripts/04-mark-dup-picard.sh" {}
./scripts/05-index.sh $p


#creating ideograms with breakpoints
#Rscript ./scripts/06-breakpointR.R $p



#extracting basic metrics
#ls /home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/sam/*.sorted.bam| rev | cut -f1 -d/ | rev | parallel --jobs $p bash "./scripts/07-align-met-picard.sh" {}
#./scripts/08-align-met-flagstat.sh
#ls /home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/mdup/*.bam | rev | cut -f1 -d/ | rev | parallel --jobs $p bash "./scripts/09-insert-size-picard.sh" {}
#ls /home/zhamadeh/NEW_PAPER_OLD_DATA/se_alignment_pipeline/mdup/*.bam | rev | cut -f1 -d/ | rev | parallel --jobs $p bash "./scripts/10-gc-bias-picard.sh" {}



#estimating and plotting complexity
#bash ./scripts/11-preseq.sh "/projects/lansdorp/for_ERIBA/mdup"


#extracting and compiling some more basic metrics
#Rscript ./scripts/12-extract-BPR-stats.R
#bash ./scripts/13-extract-gc-insertsize-stats.sh /projects/lansdorp/for_ERIBA/sspipe/metrics/CollectInsertSizeMetrics/ /projects/lansdorp/for_ERIBA/sspipe/metrics/CollectGCBiasMetrics/ /projects/lansdorp/for_ERIBA/sspipe/metrics

#compiling basic metrics into a big table
#Rscript ./scripts/14-dataman.R

