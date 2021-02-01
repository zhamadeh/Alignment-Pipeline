####### Niek's bed files and metadata ########
bpr_se = read.table("breakPointSummary_se.txt",header = T)  %>% select(-c(CI.start,CI.end,genoT))
bpr_pe = read.table("breakPointSummary_pe.txt",header = T) %>% select(-c(CI.start,CI.end,genoT))
bedfile <- read.table("nieks_bed_files/SCE_GM02085_enrichment.BED",sep="\t",header=T)
metadata <- read.table("old_data_analysis/E-MTAB-5976.sdrf.txt",sep="\t",header=T)
metadata <-  select(metadata,c(Comment.ENA_RUN., Extract.Name,Factor.Value.protocol., Characteristics.genotype., Characteristics.organism.))
#packages
library(tidyverse)
library(bedr)
library(GenomicRanges)

#  Put all bed files in one file
bed <- data.frame()
for (file in list.files("nieks_bed_files/",full.names = T)){
	tmp <- read.table(file,sep="\t")
	if (ncol(tmp)==4){
		tmp$V5=tmp$V3 - tmp$V2
		bed <- rbind(tmp,bed)
	} else {bed <- rbind(tmp,bed)}
	#message("columns: ",ncol(tmp)," row namese  ", nrow(tmp) )
}
#merge with metadata
merge <- merge(bed,metadata,by.x="V4",by.y="Extract.Name")

#sort into BLM and WT SCEs
wt <- filter(merge,merge$Characteristics.genotype.=="wild type genotype") %>% select(-c(Characteristics.organism., Characteristics.genotype., Factor.Value.protocol.,Comment.ENA_RUN.)) %>% select(c(V1,V2,V3,V5,V4))
blm <- filter(merge,merge$Characteristics.genotype.=="homozygous Blm knockout")%>% select(-c(Characteristics.organism., Characteristics.genotype., Factor.Value.protocol.,Comment.ENA_RUN.))%>% select(c(V1,V2,V3,V5,V4))
colnames(blm) = c("chr"   ,  "start"     ,   "end"  , "width"    ,    "Extract.Name" )
blm <- blm %>% select(c(chr,start,end,width,Extract.Name))
colnames(wt) = c("chr"   ,  "start"     ,   "end"  , "width"    ,    "Extract.Name" )
wt <- wt %>% select(c(chr,start,end,width,Extract.Name))



bind<- rbind(bpr_pe,bpr_se)
for (i in 1:5){
	bind$filenames<-tools::file_path_sans_ext(bind$filenames)
}
bind$filenames <- gsub( "_", "", as.character(bind$filenames) )
bind$filenames <- as.factor(bind$filenames)
str(metadata$Comment.ENA_RUN.)

bpr_merge <- merge(bind,metadata,by.x="filenames",by.y="Comment.ENA_RUN.") %>% select(-c(filenames, Factor.Value.protocol., , Characteristics.organism.))
bpr_wt <- filter(bpr_merge,bpr_merge$Characteristics.genotype.=="wild type genotype") %>% select(-c(Characteristics.genotype. ))
bpr_blm <- filter(bpr_merge,bpr_merge$Characteristics.genotype.=="homozygous Blm knockout")  %>% select(-c(Characteristics.genotype. ))
bpr_blm$width=bpr_blm$end-bpr_blm$start
bpr_wt$width=bpr_wt$end-bpr_wt$start
colnames(bpr_blm) = c("chr"   ,  "start"     ,   "end"    ,    "Extract.Name", "width"   )
colnames(bpr_wt) = c("chr"   ,  "start"     ,   "end"  ,    "Extract.Name", "width"     )
bpr_blm <- select(bpr_blm,c(chr,start,end,width,Extract.Name))
bpr_wt <- select(bpr_wt,c(chr,start,end,width,Extract.Name))
sces_per_lib_niek <- blm %>% group_by(Extract.Name) %>% dplyr::summarize(n())
sces_per_lib_niek$Extract.Name = "BAIT_BS"
sces_per_lib_bpr <- bpr_blm %>% group_by(Extract.Name) %>% dplyr::summarize(n())
sces_per_lib_bpr$Extract.Name = "BPR_BS"
sces_per_lib_niek_wt <- wt %>% group_by(Extract.Name) %>% dplyr::summarize(n())
sces_per_lib_niek_wt$Extract.Name = "BAIT_WT"
sces_per_lib_bpr_wt <- bpr_wt %>% group_by(Extract.Name) %>% dplyr::summarize(n())
sces_per_lib_bpr_wt$Extract.Name = "BPR_WT"



BAIT_BS <- GRanges(blm)
BAIT_WT <- GRanges(wt)
BPR_BS <- GRanges(bpr_blm)
BPR_WT <- GRanges(bpr_wt)



TP <- length(GenomicRanges::intersect(BAIT_BS,BPR_BS))
FP <- length(GenomicRanges::setdiff(BAIT_BS,BPR_BS))
FN <- length(GenomicRanges::setdiff(BPR_BS,BAIT_BS))

test[-queryHits(findOverlaps(test, centroGRange, type="any")),]




sces_per_lib <- rbind(sces_per_lib_bpr,sces_per_lib_niek,sces_per_lib_bpr_wt,sces_per_lib_niek_wt)

ggplot(sces_per_lib)+geom_jitter(aes(Extract.Name,`n()`,alpha=0.1))+
	geom_boxplot(aes(Extract.Name,`n()`,fill=Extract.Name))+
	theme_classic()


blm$package <- "BAIT_BS"
wt$package <- "BAIT_WT"
bpr_blm$package <- "BPR_BLM"
bpr_wt$package <- "BPR_WT"

join <- rbind(blm,wt,bpr_blm,bpr_wt)

ggplot(join, aes(width)) + stat_ecdf(geom = "step",aes(group=package,color=package)) +
	scale_x_log10() +
	theme_classic()

