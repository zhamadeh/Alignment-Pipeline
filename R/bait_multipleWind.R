################################################
# Packages #
################################################
library(tidyverse)
library(bedr)
library(GenomicRanges)
library(data.table)
library(tools)

################################################
# Data #
################################################

bait_se_200 = read.table("BAIT/BAIT_TOTAL_LIBRARY_SCE_200Kb-window.txt",header = T)  %>% select(-c(transition,callType,aneuploidy))
bait_se_1e6 = read.table("BAIT/BAIT_TOTAL_LIBRARY_SCE_1Mb-window.txt",header = T) %>% select(-c(transition,callType,aneuploidy))


metadata <- read.table("Anxilliary/Metadata/E-MTAB-5976.sdrf.txt",sep="\t",header=T)
metadata <-  select(metadata,c(Comment.ENA_RUN., Extract.Name,Factor.Value.protocol., Characteristics.genotype., Characteristics.organism.))

#  Put all bed files in one file
bed <- data.frame()
for (file in list.files("Benchmark_sce_files/",full.names = T)){
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


bait_se_200$window="200Kb"
bait_se_1e6$window="1Mb"
bind<- rbind(bait_se_200,bait_se_1e6)

bind$filenames<-tools::file_path_sans_ext(bind$index)
bind$filenames <- as.factor(bind$filenames)
str(metadata$Comment.ENA_RUN.)

bait_merge <- merge(bind,metadata,by.x="filenames",by.y="Comment.ENA_RUN.") %>% select(-c(filenames, Factor.Value.protocol., Characteristics.organism.))

bait_wt <- filter(bait_merge,bait_merge$Characteristics.genotype.=="wild type genotype") %>% select(-c(Characteristics.genotype. ))
bait_blm <- filter(bait_merge,bait_merge$Characteristics.genotype.=="homozygous Blm knockout")  %>% select(-c(Characteristics.genotype. ))

bait_blm$width=bait_blm$end - bait_blm$start
bait_wt$width=bait_wt$end - bait_wt$start

colnames(bait_blm) = c("Chromosome"   ,  "start"     ,   "end"  ,"index"  , "window",   "Extract.Name", "width"   )
colnames(bait_wt) = c("Chromosome"   ,  "start"     ,   "end"  ,"index"  ,  "window",   "Extract.Name", "width"   )
bait_blm <- select(bait_blm,c(Chromosome,start,end,width,Extract.Name,window,index))
bait_wt <- select(bait_wt,c(Chromosome,start,end,width,Extract.Name,window,index))

blm$package <- "niek_BS"
wt$package <- "niek_WT"
bait_blm$package <- "bait_BLM"
bait_wt$package <- "bait_WT"

bait_blm$window <- as.factor(bait_blm$window)
bait_wt$window <- as.factor(bait_wt$window)

bait <- rbind(bait_blm,bait_wt)
bait <- select(bait,-c(index))


benchmark <- rbind(wt,blm)
benchmark$window=NA


colnames(bait) = c("chr"   ,  "start"     ,   "end" , "width" , "Extract.Name"  , "window" ,"package")
bait <- select(bait,c(chr,start,end,width,Extract.Name,package,window))

benchmark=benchmark[benchmark$Extract.Name %in% levels(bait$Extract.Name),]

fullDataset <- rbind(bait,benchmark)

##################### dataset complete #####################

##################### dataset complete #####################


BAIT <- GRanges(bait)
BENCHMARK = GRanges(benchmark)


prec_recall=data.frame()
level=levels(bait$window)[1]
i=0

for (level in levels(bait$window)){
	bait_level <- filter(bait,window==level)

	BAIT <- GRanges(bait_level)


	for (i in c(0,10000,100000,500000)){

		(TP <- length(countOverlaps(BENCHMARK,BAIT)[countOverlaps(BENCHMARK,BAIT,maxgap = i) !=0]))
		(FN <- length(BENCHMARK[-queryHits(findOverlaps(BENCHMARK,BAIT, type="any",maxgap = i)),]))
		(FP <-  length(BAIT[-queryHits(findOverlaps(BENCHMARK,BAIT, type="any",maxgap = i)),]))
		FNpTP=FN+TP
		TPpFP=TP+FP
		precision=(TP/(TP+FP))
		recall=(TP/(TP+FN))
		tmp = data.frame("TP"=TP,"FN"=FN,"FP"=FP,"FNpTP"=FNpTP,"TPpFP"=TPpFP,"level"=level,"gap"=i,"precision"=precision,"recall"=recall)
		prec_recall <- rbind(prec_recall,tmp)

	}
}


prec_recall$gap<- as.factor(as.character(prec_recall$gap))
prec_recall$level<- as.factor(as.character(prec_recall$level))
prec_recall$sum=prec_recall$TP+prec_recall$FN+prec_recall$FP
str(prec_recall)

ggplot(prec_recall)+geom_point(aes(recall,precision,shape=gap,color=approach),size=4)+
	geom_line(aes(recall,precision,group=level,color=approach),size=1)+
	#scale_colour_manual(values=c("100_reads.0"="#b8200f","1000_size.0"="#10917c","175_reads.0"="#64167d" ,"50_reads.0"="#116b26","100_reads.10"="#b8200f" ,"1000_size.10"="#10917c" ,"175_reads.10"="#64167d","50_reads.10"="#116b26"))+
	xlim(0,1)+ylim(0,1)+
	#stat_smooth(data=filter(prec_recall,gap=="0"),aes(recall,precision),method="lm",fullrange=TRUE,se=F,color="black")+
	theme_bw()+
	ggsave("prec-recall.png")



ggplot(fullDataset, aes(width)) + stat_ecdf(geom = "step",aes(group=package,color=package)) +
	#geom_density(aes(group=package,color=package))+
	scale_x_log10() + labs(x="Resolution",y="Cumulative density") +
	theme_classic()



prec_recall$approach="BreapointR"
prec_recall1<-select(prec_recall1,-c(trim))
prec_recall2<-prec_recall
prec_recall2$approach="BAIT"

prec_recall=rbind(prec_recall1,prec_recall2)
