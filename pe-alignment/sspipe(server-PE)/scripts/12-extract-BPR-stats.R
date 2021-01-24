#This script manipulates the breakpointR output and extracts some of the calculated metrics.
#Some processing is sample name specific, namely the "Library" line which processes filenames by removing 17 characters. May need to be edited.

#######################################################################################################################################

#Move to the data/ subdirectory of the breakpointR (BPR) output.
setwd("/projects/lansdorp/for_ERIBA/mdup/BPR_output/data")

library(breakpointR)


#Searching for libraries
list <- list.files(pattern='RData$')

#Loading BreakpointR RData files
init <- loadFromFiles(list)



#extracting some BPR stats with reasonable formatting from the weird BreakpointR output objects
names_init <- unname(unlist(lapply(init, function(x)x$ID[1])))
#to get reasonable library names, I truncate them----but this varies by naming scheme
Library <- substr(names_init, 1, nchar(names_init)-17)
Background <- unname(unlist(lapply(init, function(x)x$lib.metrics[1])))
Reads_per_Mb <- unname(unlist(lapply(init, function(x)x$lib.metrics[2])))
Num_regions <- unname(unlist(lapply(init,function(x)length(x$counts))))
Percent_WC <- unname(unlist(lapply(init,function(x)length(grep("wc",mcols(x$counts)$states)))) / Num_regions)
BPR_coverage <- unname(unlist(lapply(init, function(x)x$lib.metrics[3])))/100

#mashing it all into a dataframe and printing it to the screen 
t <- data.frame(Library,BPR_coverage,Background,Reads_per_Mb,Num_regions,Percent_WC)
print(t)

#changing to a new directory where the BPR stats can be stored. 
setwd("/projects/lansdorp/for_ERIBA/sspipe/metrics")


#Writing the BPR stats to a text file
write.table(t,file="BPR_stats.txt",sep="\t",quote=FALSE,row.names=FALSE)

