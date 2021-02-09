#server script to wget NV's data from the web

#load packages and data
metadata <- read.table("E-MTAB-5976.sdrf.txt",header=T,fill=T,sep="\t")
library(tidyverse)

#filter out metadata for humans
human<- filter(metadata,metadata$Characteristics.organism.=="Homo sapiens")

#Survey data
human%>%group_by(Characteristics.disease.,Characteristics.cell.type.)%>%summarize(n())
# 1232 BS libraries, 1595 WT
# 567 EBV BS, 665 Fibrorblast BS, 901 EBV WT, 694 Fibroblast WT
human%>%group_by(Characteristics.disease.,human$Comment.LIBRARY_LAYOUT.)%>%summarize(n())
# BS: 1232 single-end, WT: 434 PE + 1161 SE
human%>%group_by(Characteristics.disease.,Factor.Value.protocol.)%>%summarize(n())
#need to remove 8 RNA-seq libraries
humanStrandSeq <- filter(human, human$Factor.Value.protocol. =="Strand-Seq")


#get fastq URLs
humanStrandSeqFastqURLs <- select(humanStrandSeq,Comment.FASTQ_URI.)
humanStrandSeqFastqURLs$Comment.FASTQ_URI.<- paste("wget",humanStrandSeqFastqURLs$Comment.FASTQ_URI.)
write.table(humanStrandSeqFastqURLs,"humanStrandSeqFastqURLs.txt",quote=F,row.names = F,col.names = F)
#then add shebang line to above and save as shell script, scp onto server and run

#ran two alignment pipelines on the data, one for PE fastq and one for SE fastq files
#PE needed the following adaptors trimmed:


