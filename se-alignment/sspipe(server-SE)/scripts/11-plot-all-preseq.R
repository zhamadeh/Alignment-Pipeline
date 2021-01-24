#If a file exists in the current directory with preseq data for many libraries, then this script will plot the complexity curves with mean and std. dev.

##########################################################################################################################################################################################################################

library(ggplot2)
library(reshape2)

#Read in table of all preseq data
all_table <- read.table("preseq_all_files.plot", header=T)

#extract the mean and std. dev. of the complexity curves 
mean_table <- data.frame(total_Gb=all_table[,1], Mean=rowMeans(all_table[,-1]), SD1=rowMeans(all_table[,-1]) - apply(all_table[,-1],1,sd), SD2=rowMeans(all_table[,-1]) + apply(all_table[,-1],1,sd))

write.table(mean_table, file="all_mean_with_stdev.preseq.plot",quote=FALSE,sep="\t",row.names=FALSE)

#Use melt to convert a table with many columns into one with 2 columns plus a variable field, which is what ggplot handles best
all_table_long <- melt(all_table, id="total_Gb")
 

#plot both the raw complexity curves and their mean and standard deviation
g <- ggplot(all_table_long,aes(x=total_Gb, y=value, colour=variable)) + theme_classic() + geom_line() + theme(legend.position="none") +
labs(caption="Complexity curves for all libraries (light grey lines) and their mean and standard deviation (blue line and ribbon)", y="predicted haploid genome fraction (%)", x="total sequencing effort (Gb)") + 
scale_x_continuous(breaks=seq(0,2,0.2),limits=c(0,2)) +
scale_colour_manual(values=rep("grey95",ncol(all_table)-1)) + geom_line(data = mean_table, aes(x = total_Gb, y = Mean), color = "steelblue3", size=1) + theme(legend.position="none") +
geom_ribbon(data=mean_table, aes(x = total_Gb, ymin = SD1, ymax = SD2), linetype = 2, alpha = 0.5, inherit.aes=FALSE, fill="slategray2") +
theme(plot.caption = element_text(hjust = 0))

ggsave(filename="all_plots.pdf", device="pdf", plot=g)
