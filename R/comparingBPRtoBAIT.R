#packages
library(tidyverse)
library(bedr)
library(GenomicRanges)
library(data.table)


####### Niek's bed files and metadata ########
list.files("bpr_multipleWindows/breakpoints/",full.names = T)

bpr_se = read.table("Breakpoint_analysis/bpr_multipleWindows/bpr-standard/breakPointSummary_se.txt",header = T)  %>% select(-c(CI.start,CI.end,genoT))
bpr_pe = read.table("Breakpoint_analysis/bpr_multipleWindows/bpr-standard//breakPointSummary_pe.txt",header = T) %>% select(-c(CI.start,CI.end,genoT))

bpr_se_1000 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_1000_size_se.txt",header = T)  %>% select(-c(CI.start,CI.end,genoT))
bpr_pe_1000 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_1000_size_pe.txt",header = T) %>% select(-c(CI.start,CI.end,genoT))

bpr_se_100 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_100_reads_pe.txt",header = T)  %>% select(-c(CI.start,CI.end,genoT))
bpr_pe_100 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_100_reads_se.txt",header = T) %>% select(-c(CI.start,CI.end,genoT))

bpr_se_50 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_50_reads_se.txt",header = T)  %>% select(-c(CI.start,CI.end,genoT))
bpr_pe_50 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_50_reads_pe.txt",header = T) %>% select(-c(CI.start,CI.end,genoT))

bpr_se_200 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_200_reads_se.txt",header = T)  %>% select(-c(CI.start,CI.end,genoT))
bpr_pe_200 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_200_reads_pe.txt",header = T) %>% select(-c(CI.start,CI.end,genoT))

bpr_se_2e05 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_2e+05_size_se.txt",header = T)  %>% select(-c(CI.start,CI.end,genoT))
bpr_pe_2e05 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_2e+05_size_pe.txt",header = T) %>% select(-c(CI.start,CI.end,genoT))

bpr_se_5000 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_5000_size_se.txt",header = T)  %>% select(-c(CI.start,CI.end,genoT))
bpr_pe_5000 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_5000_size_pe.txt",header = T) %>% select(-c(CI.start,CI.end,genoT))

bpr_se_800 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_800_reads_se.txt",header = T)  %>% select(-c(CI.start,CI.end,genoT))
bpr_pe_800 = read.table("Breakpoint_analysis/bpr_multipleWindows/breakpoints/breakpoints_800_reads_pe.txt",header = T) %>% select(-c(CI.start,CI.end,genoT))


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



bind1<- rbind(bpr_pe,bpr_se)
bind1$window="175_reads"
bind2<- rbind(bpr_pe_1000,bpr_se_1000)
bind2$window="1000_size"
bind3<- rbind(bpr_pe_100,bpr_se_100)
bind3$window="100_reads"
bind4<-rbind(bpr_se_50,bpr_pe_50)
bind4$window="50_reads"

bind5<-rbind(bpr_se_200,bpr_pe_200)
bind5$window="200_reads"
bind6<-rbind(bpr_se_2e05,bpr_pe_2e05)
bind6$window="2e05_size"
bind7<-rbind(bpr_se_5000,bpr_pe_5000)
bind7$window="5000_size"
bind8<-rbind(bpr_se_800,bpr_pe_800)
bind8$window="800_reads"

bind<- rbind(bind1,bind2,bind3,bind4,bind5,bind6,bind7,bind8)


for (i in 1:5){
	bind$filenames<-tools::file_path_sans_ext(bind$filenames)
}
bind$filenames <- gsub( "_", "", as.character(bind$filenames) )
bind$filenames <- as.factor(bind$filenames)
str(metadata$Comment.ENA_RUN.)

bpr_merge <- merge(bind,metadata,by.x="filenames",by.y="Comment.ENA_RUN.") %>% select(-c(filenames, Factor.Value.protocol., Characteristics.organism.))

bpr_wt <- filter(bpr_merge,bpr_merge$Characteristics.genotype.=="wild type genotype") %>% select(-c(Characteristics.genotype. ))
bpr_blm <- filter(bpr_merge,bpr_merge$Characteristics.genotype.=="homozygous Blm knockout")  %>% select(-c(Characteristics.genotype. ))

bpr_blm$width=bpr_blm$end-bpr_blm$start
bpr_wt$width=bpr_wt$end-bpr_wt$start
colnames(bpr_blm) = c("chr"   ,  "start"     ,   "end"    , "window",   "Extract.Name", "width"   )
colnames(bpr_wt) = c("chr"   ,  "start"     ,   "end"  ,  "window",   "Extract.Name", "width"     )

bpr_blm <- select(bpr_blm,c(chr,start,end,width,Extract.Name,window))
bpr_wt <- select(bpr_wt,c(chr,start,end,width,Extract.Name,window))

blm$package <- "BAIT_BS"
wt$package <- "BAIT_WT"
bpr_blm$package <- "BPR_BLM"
bpr_wt$package <- "BPR_WT"

bpr_blm$window <- as.factor(bpr_blm$window)
bpr_wt$window <- as.factor(bpr_wt$window)

bait <- rbind(wt,blm)
bait$window=NA
bait$trim=NA

fullDataset <- rbind(bait,bpr_blm,bpr_wt)
##################### dataset complete #####################


#inspecting SCEs/library
sces_per_lib_blm_bait <- blm %>% group_by(Extract.Name) %>% summarize(n())
sces_per_lib_blm_bait$package  = "BAIT_BS"
sces_per_lib_wt_bait <- wt %>% group_by(Extract.Name) %>% dplyr::summarize(n())
sces_per_lib_wt_bait$package  = "BAIT_WT"
sces_per_lib_BAIT <- rbind(sces_per_lib_blm_bait,sces_per_lib_wt_bait)
sces_per_lib_BAIT$window <- 0
sces_per_lib_BAIT$package<-as.factor(sces_per_lib_BAIT$package)
sces_per_lib_BAIT$window <- as.factor(as.character(sces_per_lib_BAIT$window))
sces_per_lib_BAIT <- select(sces_per_lib_BAIT,c(Extract.Name,window,`n()`,package))
sces_per_lib_BAIT <- as.data.frame(sces_per_lib_BAIT)


sces_per_lib_blm_bpr <- bpr_blm %>% group_by(Extract.Name,window) %>% dplyr::summarize(n())
sces_per_lib_blm_bpr$package  = "BPR_BS_untrim"
sces_per_lib_wt_bpr <- bpr_wt %>% group_by(Extract.Name,window) %>% dplyr::summarize(n())
sces_per_lib_wt_bpr$package = "BPR_WT_untrim"
sces_per_lib_BPR <- rbind(sces_per_lib_blm_bpr,sces_per_lib_wt_bpr)
sces_per_lib_BPR$package<-as.factor(sces_per_lib_BPR$package)
sces_per_lib_BPR<-as.data.frame(sces_per_lib_BPR)



#collect list libraries with SCE/lib counts above 90th perctentile
blm_bpr_libsToRemove <- setDT(sces_per_lib_blm_bpr)[,.SD[quantile(`n()`, 0.90) < 	`n()`]]$Extract.Name
wt_bpr_libsToRemove <- setDT(sces_per_lib_wt_bpr)[,.SD[quantile(`n()`, 0.90) < 	`n()`]]$Extract.Name
#retabulate SCE/s per
filtered_sces_bpr_blm <- sces_per_lib_blm_bpr[!sces_per_lib_blm_bpr$Extract.Name %in% blm_bpr_libsToRemove, ]
filtered_sces_bpr_blm$package  = "BPR_BS_filt"
filtered_sces_bpr_wt<- sces_per_lib_wt_bpr[!sces_per_lib_wt_bpr$Extract.Name %in% wt_bpr_libsToRemove, ]
filtered_sces_bpr_wt$package = "BPR_WT_filt"
sces_per_lib_filt <- rbind(filtered_sces_bpr_blm,filtered_sces_bpr_wt)
sces_per_lib_filt$package<-as.factor(sces_per_lib_filt$package)
sces_per_lib_filt <- as.data.frame(sces_per_lib_filt)




mergeScesPerLib <- rbind(sces_per_lib_BAIT,sces_per_lib_BPR,sces_per_lib_filt)
mergeScesPerLib$level <- paste0(as.character(mergeScesPerLib$package),"_",as.character(mergeScesPerLib$window))
mergeScesPerLib$level<-as.factor(mergeScesPerLib$level)

mergeOneWindowSize_filt <- rbind(filter(sces_per_lib_filt,window=="175_reads"),sces_per_lib_BAIT)
mergeOneWindowSize <- rbind(filter(sces_per_lib_BPR,window=="175_reads"),sces_per_lib_BAIT)

ggplot(mergeScesPerLib)+geom_jitter(aes(level,`n()`,alpha=0.1))+
	geom_boxplot(aes(level,`n()`,fill=level))+
	theme_classic() +
	labs(x="Package",y="SCEs/library")+
	#scale_fill_manual(values=c("#1f69a1", "#19943c","#318fd6","#31d660")) +
	theme(legend.position= "none",axis.text.x = element_text(angle = 90, hjust = 1))+
	ggsave("SCE_lib_all_window_sizes_filt.png")

ggplot(mergeOneWindowSize_filt)+geom_jitter(aes(package,`n()`,alpha=0.1))+
	geom_boxplot(aes(package,`n()`,fill=package))+
	theme_classic() +
	labs(x="Package",y="SCEs/library")+
	scale_fill_manual(values=c("#1f69a1", "#19943c","#318fd6","#31d660")) +
	theme(legend.position= "none")  + ggsave("SCE_lib_filt_175_reads.png")

ggplot(mergeOneWindowSize)+geom_jitter(aes(package,`n()`,alpha=0.1))+
	geom_boxplot(aes(package,`n()`,fill=package))+
	theme_classic() +
	labs(x="Package",y="SCEs/library")+
	scale_fill_manual(values=c("#1f69a1", "#19943c","#318fd6","#31d660")) +
	theme(legend.position= "none")  + ggsave("SCE_lib.png")

######################## SCE/lib complete ########################

######################## breakpoint resolution #################

#use filtered dataset
bpr_blm_V2 = bpr_blm[!bpr_blm$Extract.Name %in% blm_bpr_libsToRemove, ]
bpr_wt_V2 = bpr_wt[!bpr_wt$Extract.Name %in% wt_bpr_libsToRemove, ]

#try using only overlapping elements
BPR <- GRanges(bpr)
BAIT<- GRanges(bait)

bpr_validated=as.data.frame(subsetByOverlaps(BPR,BAIT)) %>% select(-c(strand,window))
colnames(bpr_validated) = c("chr" ,    "start"  ,      "end"   ,       "width"    ,    "Extract.Name" ,"package" )

bpr_blm_V2$package <- "BPR_BLM"
bpr_wt_V2$package <- "BPR_WT"

bpr_blm_V3 = bpr_blm_V2[bpr_blm_V2$width < quantile(bpr_blm_V2$width , 0.95), ]
bpr_wt_V3 = bpr_wt_V2[bpr_wt_V2$width < quantile(bpr_wt_V2$width , 0.95), ]

joinV1 <- rbind(bait,bpr_blm,bpr_wt)
joinV1 <- rbind(bait,bpr_validated)
joinV2 <- rbind(bait,bpr_blm_V2,bpr_wt_V2)
joinV3 <- rbind(bait,bpr_blm_V3,bpr_wt_V3)


ggplot(joinV1, aes(width)) + stat_ecdf(geom = "step",aes(group=package,color=package)) +
	#geom_density(aes(group=package,color=package))+
	scale_x_log10() +labs(x="Resolution",y="Cumulative density") +
	theme_classic() + ggsave("resolution_V1.png")
ggplot(joinV2, aes(width)) + stat_ecdf(geom = "step",aes(group=package,color=package)) +
	#geom_density(aes(group=package,color=package))+
	scale_x_log10() +labs(x="Resolution",y="Cumulative density") +
	theme_classic() + ggsave("resolution_V2.png")
ggplot(joinV3, aes(width)) + stat_ecdf(geom = "step",aes(group=package,color=package)) +
	#geom_density(aes(group=package,color=package))+
	scale_x_log10() + labs(x="Resolution",y="Cumulative density") +
	theme_classic() + ggsave("resolution_V3")

fullDataset$package<- as.factor(fullDataset$package)
fullDataset$window<-as.factor(fullDataset$window)
fullDataset$trim <-as.factor(fullDataset$trim )

ggplot(fullDataset, aes(width)) + stat_ecdf(geom = "step",aes(group=interaction(package,window),color=interaction(package,window))) +
	#geom_density(aes(group=package,color=package))+
	scale_x_log10() + labs(x="Resolution",y="Cumulative density") +
	theme_classic()# + ggsave("resolution_V3")







bait <- rbind(blm,wt)
BAIT <- GRanges(bait)

bpr=rbind(bpr_blm,bpr_wt)
bprTim <- rbind(bpr_blm_V3,bpr_wt_V3)

bpr$trim <- "0"
bprTim$trim <- "10"
bpr_untrimmed_backup=bpr

bpr=rbind(bpr_untrimmed_backup,bprTim)
bpr$trim<-as.factor(bpr$trim)

prec_recall=data.frame()
trim="0"
bpr_trim <- filter(bpr,trim=="0")

for (level in levels(bpr$window)){
		bpr_level <- filter(bpr_trim,window==level)

		BPR <- GRanges(bpr_level)


		for (i in c(0,10000,100000,500000)){

			(TP <- length(countOverlaps(BAIT,BPR)[countOverlaps(BAIT,BPR,maxgap = i) !=0]))
			(FN <- length(BAIT[-queryHits(findOverlaps(BAIT,BPR, type="any",maxgap = i)),]))
			(FP <-  length(BPR[-queryHits(findOverlaps(BAIT,BPR, type="any",maxgap = i)),]))
			FNpTP=FN+TP
			TPpFP=TP+FP
			precision=(TP/(TP+FP))
			recall=(TP/(TP+FN))
			tmp = data.frame("TP"=TP,"FN"=FN,"FP"=FP,"FNpTP"=FNpTP,"TPpFP"=TPpFP,"trim"=trim,"level"=level,"gap"=i,"precision"=precision,"recall"=recall)
			prec_recall <- rbind(prec_recall,tmp)

		}
	}
trim="10"
bpr_trim <- filter(bpr,trim=="10")

for (level in levels(bpr$window)){
		bpr_level <- filter(bpr_trim,window==level)

		BPR <- GRanges(bpr_level)


		for (i in c(0,10000,100000,500000)){

			(TP <- length(countOverlaps(BAIT,BPR)[countOverlaps(BAIT,BPR,maxgap = i) !=0]))
			(FN <- length(BAIT[-queryHits(findOverlaps(BAIT,BPR, type="any",maxgap = i)),]))
			(FP <-  length(BPR[-queryHits(findOverlaps(BAIT,BPR, type="any",maxgap = i)),]))
			FNpTP=FN+TP
			TPpFP=TP+FP
			precision=(TP/(TP+FP))
			recall=(TP/(TP+FN))
			tmp = data.frame("TP"=TP,"FN"=FN,"FP"=FP,"FNpTP"=FNpTP,"TPpFP"=TPpFP,"trim"=trim,"level"=level,"gap"=i,"precision"=precision,"recall"=recall)
			prec_recall <- rbind(prec_recall,tmp)

		}
	}


prec_recall$gap<- as.factor(as.character(prec_recall$gap))
prec_recall$level<- as.factor(as.character(prec_recall$level))
prec_recall$sum=prec_recall$TP+prec_recall$FN+prec_recall$FP
str(prec_recall)

prec_recall<-prec_recall %>% filter(grepl("reads$", level)) %>% filter(trim=="10")
ggplot(prec_recall)+geom_point(aes(recall,precision,shape=gap,color=interaction(level,trim)),size=4)+
	geom_line(aes(recall,precision,group=interaction(level,trim),color=interaction(level,trim)),size=1)+
	#scale_colour_manual(values=c("100_reads.0"="#b8200f","1000_size.0"="#10917c","175_reads.0"="#64167d" ,"50_reads.0"="#116b26","100_reads.10"="#b8200f" ,"1000_size.10"="#10917c" ,"175_reads.10"="#64167d","50_reads.10"="#116b26"))+
	xlim(0,1)+ylim(0,1)+
	#stat_smooth(data=filter(prec_recall,gap=="0"),aes(recall,precision),method="lm",fullrange=TRUE,se=F,color="black")+
	theme_bw()+
	ggsave("prec-recall.png")




reads<-filter(prec_recall,grepl("reads$", level))

size<-filter(prec_recall,grepl("size$", level))

ggplot(reads)+geom_point(aes(recall,precision,shape=gap,color=interaction(level,trim)),size=4)+
	geom_line(aes(recall,precision,group=interaction(level,trim),color=interaction(level,trim)),size=1)+
	scale_colour_manual(values=c("100_reads.0"="#b8200f" ,"200_reads.0"="#10917c" ,"175_reads.0"="#64167d" ,"50_reads.0"="#116b26" ,"800_reads.0"="#3248a8",
								 "100_reads.10"="#b8200f","200_reads.10"="#10917c","175_reads.10"="#64167d","50_reads.10"="#116b26","800_reads.10"="#3248a8" ))+
	xlim(0,1)+ylim(0,1)+
	#stat_smooth(data=filter(reads,gap=="0"),aes(recall,precision),fullrange=TRUE,se=F,color="black")+
	theme_bw()+
	ggsave("Plots/reads_precRecall.png")
ggplot(size)+geom_point(aes(recall,precision,shape=gap,color=interaction(level,trim)),size=4)+
	geom_line(aes(recall,precision,group=interaction(level,trim),color=interaction(level,trim)),size=1)+
	scale_colour_manual(values=c("1000_size.0"="#10917c","2e05_size.0"="#64167d" ,"5000_size.0"="#116b26" ,"1000_size.10"="#10917c","2e05_size.10"="#64167d" ,"5000_size.10"="#116b26"))+
	xlim(0,1)+ylim(0,1)+
	#stat_smooth(data=filter(prec_recall,gap=="0"),aes(recall,precision),method="lm",fullrange=TRUE,se=F,color="black")+
	theme_bw()+
	ggsave("Plots/size_precRecall.png")




