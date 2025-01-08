# Used to extracted peaks from macs calling (calibrated and overlapping ones are iteratively removed) from multiple samples. Then the peaks are normalized for peak scores within each sample. Peaks are later undergoing iterative removal process (for overlapping ones)


library(optparse)

option_list = list(
 make_option(c("-s", "--sample_info_path"), 
	type="character",
	default=NULL,
	help = "samples_for_merging,use \",\" to separate multiple input", 
	metavar="character"),
  make_option(c("-M","--merging_name"), 
	type="character",
	default="",
	help = "merging  name",
	metavar="character"),
 make_option(c("-o","--output_path"),
	type="character",
	default="./",
	help="path to output directory",
	metavar="characer")
)



library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(tibble)
library(data.table)
library(reshape)

library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(GenomicRanges)
library(future)

###get input parameters
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)


## specify the additional options:
opt$sample_info_path <- "...../scripts_for_sharing/data_preprocessing/input/sample_info.txt"
opt$merging_name <- "NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed"
opt$output_path <- "results/"

out_path <- opt$output_path
if (!file.exists(out_path)){dir.create(out_path)}
if (!file.exists(paste0(out_path,"/peaks"))){dir.create(paste0(out_path,"/peaks"))}

sample_info <- opt$sample_info_path %>% read.table(head=TRUE,sep="\t")
merging <- opt$merging_name

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 200 * 1024^3) # for 100 Gb RAM

samples <- sample_info$sample_id %>% as.character

######################################################## Experiment with different method:
## run normalization and clustering directly
all_peaks <- NULL
for (i in 1:nrow(sample_info)){
	sample_id <- sample_info[i,"sample_id"]
	cat ("read in recentered peaks for ",sample_id,"...\n")
	filtered_centered_peaks_path <- sample_info[i,"Path"] %>% 
	                        gsub(paste0(sample_id,"_processed_multiomic.rds"),paste0("recentered_final.filtered",sample_id,".tsv"),.)

	peaks <- read.table(filtered_centered_peaks_path,head=TRUE,sep="\t")
	
	peaks$Sample <- sample_id
	sum_score=sum(peaks$neg_log10qvalue_summit)/1000000
	peaks$norm_score=peaks$neg_log10qvalue_summit/sum_score
	all_peaks=rbind(all_peaks,peaks)
}

p <- all_peaks

peaks_gr <- StringToGRanges(p$new_peak, sep = c("-", "-"))
seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,peaks_gr) #extract fasta sequence
names(seq) <- p$new_peak
peaks.match.pattern <- vmatchPattern("N", seq) #match peak sequence with N in them
peaks.withN <- names(peaks.match.pattern)[elementNROWS(peaks.match.pattern)>0] # these are peaks that contain "N"
p <- p %>% dplyr::filter(!(new_peak %in% peaks.withN)) %>% dplyr::filter(seqnames!="chrY") ## remove peaks where 'N' is contained and peaks from Y chromosomes

write.table(p,paste0(out_path,"/peaks/MACS2_peaks.",merging,".BySample.tsv",sep=""),sep='\t',quote=FALSE,row.names=FALSE)



p <- read.table(paste0(out_path,"/peaks/MACS2_peaks.",merging,".BySample.tsv",sep=""),sep='\t',header=TRUE)
p1 <- p
recentered_p <- StringToGRanges(p1$new_peak, sep = c("-", "-"))

olap <- findOverlaps(recentered_p,recentered_p) %>% as.data.frame
olap1 <- olap %>% dplyr::filter(queryHits!=olap$subjectHits)

recentered_non_olap <- p1[-olap1$queryHits,]

###Use normalized score instead (as described in Science paper), column #19
pairs <- cbind(p1[olap1$queryHits,c(1:3,12)],olap1$queryHits,p1[olap1$subjectHits,c(1:3,12)],olap1$subjectHits)
colnames(pairs) <- c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')

###removing this step will allow to some more weaker peaks, as suggested here: https://www.archrproject.com/bookdown/the-iterative-overlap-peak-merging-procedure.html)
#pairs=pairs[pairs$score_1>=pairs$score_2,]
pairs <- pairs %>% arrange(desc(score_1))
pairs_all <- pairs

library(doParallel)
registerDoParallel(cores=30)

p$row=as.numeric(rownames(p))
p=as.data.table(p)
setkey(p,row)
all_st=NULL
all_st<-foreach(chr_n=c(1:19,"X","Y")) %dopar% {
	chr=paste("chr",chr_n,sep='')
	pairs=pairs_all[pairs_all$chr_1==chr,]
	pairs=pairs[,c(4,5,9,10)]
	key_list=unique(pairs$row_1)
	pairs=as.data.table(pairs)
	setkey(pairs,row_1)
	all_st_chr=NULL
	for (i in 1:length(key_list)){
		if (length(key_list)>0){
		p_add=p[key_list[1]]
		all_st_chr=rbind(all_st_chr,p_add)

		p_del=pairs[.(key_list[1])]
		del_row_1=c(p_del$row_1,p_del$row_2)
		key_list=key_list[!(key_list %in% del_row_1)]
	}
	print(paste(chr,length(key_list),sep=' '))
}
return(all_st_chr)
}

all_st_f=NULL
for (i in 1:length(all_st)){
	all_st_1=as.data.frame(all_st[[i]])
	all_st_1=all_st_1[!duplicated(all_st_1),]
	all_st_f=rbind(all_st_f,all_st_1)
}

all_st_f <- all_st_f %>% mutate(row=NULL)

recentered_final=rbind(recentered_non_olap,all_st_f)
write.table(recentered_final,paste0(out_path,"peaks/recentered_final.filtered.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(recentered_non_olap,paste0(out_path,"peaks/recentered_nonOverlapping.filtered.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(all_st_f,paste0(out_path,"peaks/recentered_Overlapping.filtered.tsv"),sep='\t',quote=FALSE,row.names=FALSE)


