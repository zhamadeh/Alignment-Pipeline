#This plots the complexity curves for all libraries in a directory (assuming they have already been generated), along with the std. dev and mean of the complexity curves

#######################################################################################################################################################################################

cd $1

>preseq_all_files.plot


#loop over all libraries for which the preseq analysis worked to create one large file containing preseq data for all libraries
for file in `ls *preseq.plot` 
do
	#process file name
	name=$(echo $file| sed 's/\.preseq\.plot//') 
	
	#process the library file and truncate it so that it isn't bigger than the plot size
	printf "total_Gb\tlib_"$name"\n" > tmp
	tail -n +2 $file | head -21 >> tmp
	mv tmp $file.tmp

	#Assemble preseq data into a single file for plotting
	paste preseq_all_files.plot <(awk '{print $2}' $file.tmp) > tmp
	mv tmp preseq_all_files.plot
	mv $file.tmp last.tmp
done


#Add a first column for the preseq file
paste  <(awk '{print $1}' last.tmp) preseq_all_files.plot > tmp
mv tmp preseq_all_files.plot

#Call a plotting script that uses the all-library preseq file as input
Rscript /projects/lansdorp/for_ERIBA/sspipe/scripts/11-plot-all-preseq.R

rm last.tmp

