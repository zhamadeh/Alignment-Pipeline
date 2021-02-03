BiocManager::install("breakpointR")
library(breakpointR)

for (window_size in c(1000,5000,10000,50000,200000,1000000)){
	breakpointR::breakpointr(inputfolder = "pe_alignment_pipeline/mdup",outputfolder = paste0("pe_alignment_pipeline/bpr_MultipleWind_",window_size,"size/"),
							 pairedEndReads=TRUE, numCPU=8,windowsize=window_size,binMethod="size",peakTh=0.3875,
							 min.mapq=7.75,trim=6.5,background=0.15)
}

