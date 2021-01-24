#This script implements preseq gc_extrap, which prints complexity curves for scWGS libraries. It takes a bam file as input and plots the predicted coverage for a given sequencing effort for the library. Note that "sequencing effort," measured in Mb, refers 
#to the total number of MAPPED bases after sequencing, including duplicates. Similarly, "genome fraction" is the proportion of the haploid genome that is covered at least one read, including sites with high depth (i.e. the algorithm is naive to alignment errors). 
#If we compare libraries of different read lengths (or SE vs PE), this calculation will need to be revisited by dividing the number of sequenced bases by the total read length (e.g. 76 + 76) to compare the number of unique reads.


#More notes on preseq:

#The complexity curve is a function that predicts the number of unique reads that would be obtained if given number of reads were sequences, i.e, if a library is sequenced with 100,000 reads, the complexity curve might predict 80,000 unique reads.

#One of the related publications makes the point that complexity curves for different libraries have different shapes and can cross: presumably because of different types and amount of amplification bias (e.g. by GC content) as compared with locus dropout 
#(e.g.due to pre-amplification bead cleans).

#This means that the best way to compare two libraries is to plot their complexity curves side-by-side; but it also means that we can simplify by comparing their complexity curves at a fixed point (i.e. a fixed sequencing effort). 
#This point of comparison should be chosen carefully.

#############################################################################################################################################################################################

#process file name
name=$(echo $1| sed 's/\.mdup\.bam//')

#convert BAM files to preseq's preferred input format: Mapped Read format
/home/vhanlon/programs/preseq_v2.0/bam2mr $name.mdup.bam > $name.mr

#run gc_extrap to print a complexity curve
/home/vhanlon/programs/preseq_v2.0/preseq gc_extrap -o $name.preseq.txt $name.mr


#if the preseq analysis was successful, proceed; otherwise, clean up
if test -f $name.preseq.txt; then

	#set column names, which will become axis labels.
	printf "total_Gb\tgenome_fraction\tlower\tupper\n" > $name.preseq.plot

	#process the gc_extrap output in preparation for plotting
	awk '{print $1/1000000000,$2/30990000,$3/30990000,$4/30990000}' OFS="\t" $name.preseq.txt | tail -n +2 | head -21 >> $name.preseq.plot

	#calculate genome fraction of input sequence data
	gfraction=$(samtools depth $name.mdup.bam | awk '$3>0 {sum+=1}END{print sum/30990000}')
	
	#call an R script to plot the complexity curve and calculate how much sequencing would need to be done to obtain a genome fraction of 5%
	Rscript /projects/lansdorp/for_ERIBA/sspipe/scripts/11-plot-one-preseq.R $name $gfraction

else
	#If preseq fails because of defects in the curve, remove the mr file and add the $name to a list of libraries that could not be processed.
	rm $name.mr
	printf $name"\n" >> preseq_analysis_failed.list
fi

rm $name.mr
