# 
# library(optparse)
# 
# option_list = list(
#  make_option(c("-s", "--sample_info_path"), ## /diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snATAC/Merging/NMK3_and_old_male/sample_info.txt
# 	type="character",
# 	default=NULL,
# 	help = "samples_for_merging,use \",\" to separate multiple input", # e.g. HTHK1301462,HTHK1301463
# 	metavar="character"),
#   make_option(c("-M","--merging_name"), ## case_id <- "2418"
# 	type="character",
# 	default="",
# 	help = "merging  name",
# 	metavar="character"),
#  make_option(c("-o","--output_path"),
# 	type="character",
# 	default="./",
# 	help="path to output directory",
# 	metavar="characer"),
#  make_option(c("--pc_num"),
#         type="integer",
#         default=30,
#         help = "number of principal components to use",
#         metavar="integer"),
#  make_option(c("--pc_first"),
#         type="integer",
#         default=1,
#         help = "first principal components to use (should be 1 or 2)",
#         metavar="integer")
# )





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

library(future)

# ###get input parameters
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser)


## specify the additional options:
opt <- list()
opt$sample_info_path <- "...../scripts_for_sharing/data_preprocessing/input/sample_info.txt"
opt$merging_name <- "NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed"
opt$output_path <- "results/"
opt$pc_first <- NULL

opt$meta.data_file <- "...../scripts_for_sharing/SuppTable_barcode_QC.txt"

merging <- opt$merging_name
meta.data <- read.table(opt$meta.data_file,head=TRUE,sep="\t")


sample_info <- opt$sample_info_path %>% read.table(head=TRUE,sep="\t")

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 200 * 1024^3) # for 100 Gb RAM

samples <- sample_info$sample_id %>% as.character

if (!file.exists(paste0(out_path,merging,"_snATAC_Merged.rds.gz"))){
	atac <- vector(mode = "list", length = length(samples))
	for (i in 1:nrow(sample_info)){
		cat ("read in ",sample_info[i,"sample_id"] %>% as.character,"...\n")
		atac[[i]]=readRDS(sample_info[i,"Path"] %>% as.character)
		DefaultAssay(atac[[i]]) <- "ATAC_MACS2"
		for (assay in Assays(atac[[i]])){
			if (assay!="ATAC_MACS2"){
			atac[[i]][[assay]] <- NULL	
		}
	}

	#####To obtain the best results - use ALL peaks!
	combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
	peakwidths <- width(combined.peaks)
	combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
	combined.peaks
	peaks.use <- combined.peaks


	#We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.
	matrix.counts=vector(mode = "list", length = length(samples))
	for (i in 1:length(samples)){
	    matrix.counts[[i]] <- FeatureMatrix(	
		fragments = Fragments(atac[[i]]@assays$ATAC_MACS2),
		features = peaks.use,
		sep = c("-","-"),
		cells = colnames(atac[[i]])
		)
	}

	for (i in 1:length(samples)){
		atac[[i]][['peaksinters']] <- CreateChromatinAssay(counts = matrix.counts[[i]],fragments=Fragments(atac[[i]]@assays$ATAC_MACS2))
		atac[[i]]$dataset=samples[i]
		DefaultAssay(atac[[i]])<-'peaksinters'
		atac[[i]]$ATAC_MACS2 <- NULL
	}

	####Merging:
	combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)
	saveRDS(combined, paste0(out_path,merging,"_snATAC_Merged.rds.gz"),compress=TRUE)
	rm(combined);rm(atac);gc()
}
 


##########################################################
## load recentered peaks
recentered_final=read.table(paste0(out_path,"peaks/recentered_final.filtered.tsv"),sep='\t',header=T)

## load atac object, all libraries merged
combined <- readRDS(paste0(out_path,merging,"_snATAC_Merged.rds.gz"))
atac <- combined %>% subset(cells=meta.data$barcode)
atac@meta.data$Nuclei_Dissociation_Date <- plyr::mapvalues(atac@meta.data$orig.ident,from=sample_info$sample_id,to=sample_info$Nuclei_Dissociation_Date)

recentered_p=StringToGRanges(recentered_final$new_peak, sep = c("-", "-"))
matrix.counts <- FeatureMatrix(
	fragments = Fragments(atac@assays$peaksinters),
	features = recentered_p,
	sep = c("-","-"),
	cells = colnames(atac)
)

###next is here:
atac[['peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
		fragments=Fragments(atac@assays$peaksinters))
DefaultAssay(atac)<-'peaksMACS2'

atac[['peaksinters']]<-NULL

###Overlapping ranges supplied -- check this later:

###Remove some assays
saveRDS(atac, paste0(out_path,"peaks/",merging,"_snATAC_Merged_BasedOnSelectedPeaks.rds.gz"),compress=TRUE)


library(harmony)
###################
#ends here for now#
###################
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunTFIDF(atac)
atac <- RunSVD(atac,assay = "peaksMACS2", reduction.key = 'pca_',reduction.name = 'pca') # this is actually an LSI reduction called "pca"
#atac <- RunSVD(atac)
atac <- RunHarmony(atac, "Nuclei_Dissociation_Date", reduction="pca",plot_convergence = TRUE, assay.use = "peaksMACS2",project.dim=F)
#atac <- RunUMAP(object = atac, reduction = 'harmony', dims = opt$pc_first:50)
#atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = opt$pc_first:50)

atac <- RunUMAP(object = atac, reduction = 'harmony', dims = 1:50)
atac <- FindNeighbors(object = atac, reduction = 'harmony', dims = 1:50)

###this helps:
options(future.globals.maxSize= 891289600,future.seed=TRUE)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3,resolution=seq(0.5,1,0.1))
saveRDS(atac, paste0(out_path,"peaks/",merging,"_snATAC_Merged_BasedOnSelectedPeaks_Normalized.rds.gz"),compress=TRUE)


### Add ATAC Gene activity
options(future.globals.maxSize = 200 * 1024^3) # for 100 Gb RAM
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

Annotation(atac) <- annotations

gene.activities <- GeneActivity(atac)
atac[['ATACGeneActivity']] <- CreateAssayObject(counts = gene.activities)

atac <- NormalizeData(object = atac,assay = "ATACGeneActivity",normalization.method = 'LogNormalize',scale.factor = median(atac$nCount_RNA))

saveRDS(atac, paste0(out_path,"peaks/",merging,"_snATAC_Merged_BasedOnSelectedPeaks_Normalized_With_GeneActivity.rds.gz"),compress=TRUE)

