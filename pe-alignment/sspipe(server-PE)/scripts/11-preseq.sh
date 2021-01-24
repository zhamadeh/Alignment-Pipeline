#This is a wrapper script that analyzes a batch of libraries using preseq
#It takes a path to a directory as an argument, so run it like this: bash preseq.sh /path/to/directory_with_libraries
#It will generate a complexity curve for each library and plot them individually and all together. It will also estimate the amount of sequencing required to achieve a genome fraction of 5% for each library (see "all_libraries_5perc_genome_fraction.txt")
#If we compare libraries of different read lengths (or SE vs PE), this calculation will need to be revisited by dividing the number of sequenced bases by the total read length (e.g. 76 + 76) to compare the number of unique reads.
#The file all_libraries_5perc_genome_fraction.txt gives the sequencing effort to achieve a genome fraction of 5 percent

#######################################################################################################################################

cd $1

ls *.bam > bam.list
printf "Library\tSequencing_for_5_percent_coverage\n" > all_libraries_5perc_genome_fraction.txt

#intialize a list of libraries that preseq couldn't analyze
printf "#preseq couldn't analyze the libraries below, possibly because they don't contain enough reads for a high-confidence curve\n" > preseq_analysis_failed.list

#Run preseq analysis and generate individual plots
while read name; do bash /projects/lansdorp/for_ERIBA/sspipe/scripts/11-generate-preseq.sh $name | sed 's/\.trimmed//g'| sed 's/\t$//g' >> all_libraries_5perc_genome_fraction.txt; done<bam.list

#Generate a joint plot for all libraries
bash /projects/lansdorp/for_ERIBA/sspipe/scripts/11-plot-all-preseq.sh ./


