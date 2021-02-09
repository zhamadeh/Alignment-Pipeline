#This script compiles all metrics calculated so far into two big tables, one of which is a subset of the other. Script largely courtesy of Adam.

#################################################################################################################################################################

# load packages
library(dplyr)
library(tidyverse)
library(data.table)

# Defining variables
genomesize = as.numeric("3099706404")

#set path to directory
setwd("./metrics")

## Mark duplicates
# List all txt files including sub-folders
mdup.list = list.files(path = "./mdup/", 
                       recursive = TRUE,
                       pattern = "\\.txt$", 
                       full.names = TRUE)

# Read all the files and create a ID column to store filenames
mdup = rbindlist(sapply(mdup.list, 
                        fread, 
                        simplify = FALSE),
                 use.names = TRUE, 
                 idcol = "ID" )

# Keep only filenames (not full path in ID varaible)
mdup[, ID := factor(ID, labels = basename(mdup.list))]

# Remove extensions (filename)
mdup$ID = gsub("\\..*", "", mdup$ID)

# Remove columns
mdup = select(mdup, -c(LIBRARY))


##CollectAlignmentMetrics
# List all txt files including sub-folders
colal.list = list.files(path = "./CollectAlignmentMetrics/", 
                        recursive = TRUE,
                        pattern = "\\.txt$", 
                        full.names = TRUE)

# Read all the files and create a ID column to store filenames
colal = rbindlist(sapply(colal.list, 
                         fread, 
                         simplify = FALSE),
                  use.names = TRUE, 
                  idcol = "ID" )

# Keep only filenames (not full path in ID varaible)
colal[, ID := factor(ID, labels = basename(colal.list))]

# Remove extensions (filename)
colal$ID = gsub("\\..*", "", colal$ID)

# Only select paired infomation
colal = filter(colal, CATEGORY == "PAIR")

# Data joining
alignmentMetrics = inner_join(mdup, colal, by = c("ID","ID"))

# Select columns
alignmentMetrics = select(alignmentMetrics, c(ID, 
                                              MEAN_READ_LENGTH, 
                                              TOTAL_READS, 
                                              UNMAPPED_READS, 
                                              PF_READS, 
                                              PCT_PF_READS, 
                                              PF_READS_ALIGNED, 
                                              PCT_PF_READS_ALIGNED, 
                                              READ_PAIR_DUPLICATES, 
                                              UNPAIRED_READ_DUPLICATES, 
                                              PERCENT_DUPLICATION))


# Convert data.frame to numeric matrix
alignmentMetrics[2:11] <- lapply(alignmentMetrics[2:11], as.numeric)

alignmentMetrics$DUPLICATIONS <- (alignmentMetrics$READ_PAIR_DUPLICATES*2 + alignmentMetrics$UNPAIRED_READ_DUPLICATES)

alignmentMetrics %>% rename(Library=ID,Mean_read_length=MEAN_READ_LENGTH,Total_reads=TOTAL_READS,Unmapped_reads=UNMAPPED_READS,PF_reads=PF_READS,Percent_PF_reads=PCT_PF_READS, Initial_reads_aligned=PF_READS_ALIGNED,Percent_reads_aligned=PCT_PF_READS_ALIGNED,Paired_duplicates=READ_PAIR_DUPLICATES,Unpaired_duplicates=UNPAIRED_READ_DUPLICATES,Duplication_rate=PERCENT_DUPLICATION,Total_duplicates=DUPLICATIONS) -> alignmentMetricsRenamed


#add more stats
bpr <- read.table("BPR_stats.txt",sep="\t",header=TRUE)
gc_is <- read.table("gc_insertsize_stats.txt",sep="\t",header=TRUE)
complexity <- read.table("../../mdup/all_libraries_5perc_genome_fraction.txt",sep="\t",header=TRUE)
postfiltering_reads <- read.table("postfiltering_reads.txt",sep="\t",header=TRUE)

#merging dataframes
tmp1 <- merge(alignmentMetricsRenamed, bpr, all=TRUE)
tmp2 <- merge(tmp1, gc_is, all=TRUE)
tmp3 <- merge(tmp2,postfiltering_reads, all=TRUE)
alignmentMetrics2 <- merge(tmp3, complexity, all=TRUE)

#Making the GC content a percentage
alignmentMetrics2$Mode_GC <- alignmentMetrics2$Mode_GC/100

#Calculating a naive coverage
alignmentMetrics2 <- within(alignmentMetrics2, Naive_coverage <- ((alignmentMetrics2$Postfiltering_reads_aligned - alignmentMetrics2$Total_duplicates) * alignmentMetrics2$Mean_read_length)/genomesize)

#Reordering columns
alignmentMetricsFinal <- alignmentMetrics2[,c(1,2,3,4,5,6,7,8,21,9,10,12,11,22,13,14,15,16,17,18,19,20)]


# Export alignment metrics as data table for further characterization and visualization
write.table(x = alignmentMetricsFinal, file = "./metrics_full.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = select(alignmentMetricsFinal,c(Library,Initial_reads_aligned,Postfiltering_reads_aligned,Duplication_rate,Naive_coverage,BPR_coverage,Background,Reads_per_Mb,Percent_WC,Mode_GC,Mode_insert_size,Sequencing_for_5_percent_coverage)), file = "./metrics_summary.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
