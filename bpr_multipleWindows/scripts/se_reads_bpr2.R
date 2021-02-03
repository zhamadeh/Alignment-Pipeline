BiocManager::install("breakpointR")
library(breakpointR)

for (window_size in c(800,1600)){
	breakpointR::breakpointr(inputfolder = "se_alignment_pipeline/mdup",outputfolder = paste0("se_alignment_pipeline/bpr_MultipleWind_",window_size,"reads/"),
							 pairedEndReads=FALSE, numCPU=8,windowsize=window_size,binMethod="reads",peakTh=0.3875,
							 min.mapq=7.75,trim=6.5,background=0.15)
}
