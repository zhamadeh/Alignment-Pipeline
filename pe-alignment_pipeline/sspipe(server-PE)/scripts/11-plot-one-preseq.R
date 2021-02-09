#This is the plotting script called by preseq.sh to plot an individual library's complexity curve and to append the sequencing effort required to achieve a genome fraction of 5% to stdout.

##################################################################################################################################################################################################################################################

library(ggplot2)
library(stringr)


#The arguments here are (1) the truncated name of the input library and (2) the obverved genome fraction calculated by preseq.sh
args <- commandArgs(T)



#This reads the preseq output for 200 IGSR libraries (the benchmark dataset) and then calculates the mean and std. dev. of the complexity curves for comparison.
all_table <- read.table("/projects/lansdorp/for_ERIBA/sspipe/refseq/benchmark_preseq_all_files.plot", header=T)
mean_table <- data.frame(total_Gb=all_table[,1], Mean=rowMeans(all_table[,-1]), SD1=rowMeans(all_table[,-1]) - apply(all_table[,-1],1,sd), SD2=rowMeans(all_table[,-1]) + apply(all_table[,-1],1,sd))







#Reading the preseq output into a dataframe
complexity_curve <- read.table(paste(args[1],".preseq.plot",sep=""), header=T,sep="\t",col.names=c("total_Gb","genome_fraction","lower","upper"))





#Interpolating between the points given in the preseq output to estimate the sequencing effort required to achieve a genome fraction of 5%. This is printed to stdout
approx_line <- approxfun(complexity_curve$total_Gb, complexity_curve$genome_fraction)
genome_fraction_threshold <- optimize(function(t0) abs(approx_line(t0) - 5), interval = range(complexity_curve$total_Gb))
genome_fraction_threshold$minimum <- round(genome_fraction_threshold$minimum, digits=2)


#If more than 2 Gb of sequencing is needed to get 5% coverage, it's a lost cause: just print 2 (altered to >2 later) and also replace the  genome_fraction_threshold$minimum value
if(genome_fraction_threshold$minimum>=2){ 
	cat(c(args[1],">=2","\n"),sep="\t")	
}else{
	cat(c(args[1],genome_fraction_threshold$minimum),"\n",sep="\t")
}

#process the caption so that it wraps nicely and updates library-specific info
info <- str_wrap(paste("Complexity curve and 95% CI for ", args[1], " (blue line and ribbon). This library requires ", genome_fraction_threshold$minimum," Gb sequencing for a 5% genome fraction. The grey line and ribbon represent the mean and standard deviation 
	of the IGSR benchmark dataset.", sep=""), width=105)



#open a PDF file 
pdf(paste(args[1],".preseq.plot.pdf",sep="")) 

#Plot the complexity curve. Also plotting a horizontal line representing the observed coverage. The scale on the x axis is set to match one full haploid human genome. 
#The last line plots the mean and std. dev. of the benchmark datset. 

ggplot(complexity_curve,aes(x=total_Gb, y=genome_fraction)) + theme_classic() + 
geom_ribbon(data=mean_table, aes(x = total_Gb, ymin = SD1, ymax = SD2), linetype = 2, alpha = 0.3, inherit.aes=FALSE, fill="grey85") + geom_line(data = mean_table, aes(x = total_Gb, y = Mean), colour = "grey75") + 
labs(caption=info, y="predicted haploid genome fraction (%)", x="total sequencing effort (Gb)") + theme(legend.position="none") +
geom_ribbon(aes(ymin=complexity_curve$lower, ymax=complexity_curve$upper), linetype=2, alpha=0.5,  fill="slategray2") + geom_line(colour="steelblue3") +
geom_hline(yintercept=as.numeric(args[2]), linetype=2, size=0.1) + scale_x_continuous(breaks=seq(0,2,0.2),limits=c(0,2)) +
geom_text(aes(2,as.numeric(args[2]),label = "observed genome fraction", vjust = -1, hjust=1)) +
theme(plot.caption = element_text(hjust = 0))


#close PDF
graphics.off()


