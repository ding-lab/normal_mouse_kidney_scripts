#conda activate scenic_regulation

library(Seurat)
library(SCENIC)
library(dplyr)
library(tibble)
library(Signac)
#library("Signac", lib.loc = "/diskmnt/Projects/Users/rliu/Software/miniconda/envs/Seurat_RNA_ATAC_ST_20220615/lib/R/library")

merging_name <- "NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed"
subset_name <- "PT_S1_S2_S3_100cells_per_sample"

working_dir <- "..../scripts_for_sharing/Fig4_Senic/"
setwd(working_dir)
setwd(paste0(working_dir,"results/"))

############################## Load data ##############################
# loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom"
# loom <- open_loom(loomPath)
# exprMat <- get_dgem(loom)
# cellInfo <- get_cell_annotation(loom)
# close_loom(loom)

#obj_joint <- readRDS(paste0("/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/new_analysis_20240608/mouse/sample_integration_RNA_ATAC_merged/",merging_name,"/results/snRNA_ATAC_jointly_analyzed_for_",merging_name,".rds"))
meta.data <- read.table("..../scripts_for_sharing/SuppTable_barcode_QC.txt",head=TRUE,sep="\t")

## subset PT populations
selected_BC_all <- NULL
n_sample <- 100
selected_populations <- c("PT(S1)","PT(S2)","PT(S3)")
for (s in unique(meta.data$sample_id %>% as.character)){
	for (p in selected_populations){
		BC <- meta.data %>% rownames_to_column %>% filter(cell_type==p) %>% filter(sample_id==s) %>% .$rowname
		if (length(BC)>n_sample){selected_BC <- sample(BC,100)}else{selected_BC <- BC}
		selected_BC_all <- c(selected_BC_all,selected_BC)
	}
}
selected_BC_meta.data <- meta.data[selected_BC_all,] %>% select(sample_id,Age,Sex,cell_type)
write.table(selected_BC_meta.data,file="meta.data_for_selected_BC_for_SCENIC_analysis.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=NA)

#rm(obj_joint)

## load and save cellInfo
cellInfo <- meta.data[selected_BC_all,c("cell_type","Age","Sex","sample_id","nCount_RNA","nFeature_RNA")]
saveRDS(cellInfo, file="int/cellInfo.Rds") ## save cell info


source("..../scripts_for_sharing/colorset.R")
colVars <- list(Age=col_age[c("E16.5","P0","W3","W12","W52","W92")],Sex=col_sex[c("F","M")],cell_type=col_cell_type[c("PT(S1)","PT(S2)","PT(S3)")])
saveRDS(colVars, file="int/colVars.Rds")

############################## Initialize settings ##############################
## for the files within dbDir, instructions on how to obtain them is in this tutorial: 
## https://rdrr.io/github/aertslab/SCENIC/f/vignettes/SCENIC_Setup.Rmd
## wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather
## wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather
dbDir <- "/diskmnt/Projects/Users/rliu/Tools/scenic/Resources"
scenicOptions <- initializeScenic(org="mgi", dbDir=dbDir, nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds" ### specify path for cell info
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds" #### specify path for colvars
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


############################## Co-expression network ##############################
sample_info <- read.table(paste0("scripts_for_sharing/data_preprocessing/input/sample_info.txt"),head=TRUE,sep="\t")
exprMat_list <- list()
for (i in 1:nrow(sample_info)){
	sample_id <- as.character(sample_info$sample_id)[i]
	cat(sample_id,"...\n")
	obj <- readRDS(as.character(sample_info$Path)[i])
	DefaultAssay(obj) <- "RNA"
	obj <- RenameCells(obj,new.names=paste0(sample_id,"_",obj %>% colnames))
	exprMat_list[[sample_id]] <- GetAssayData(obj,"counts")
	exprMat_list[[sample_id]] <- exprMat_list[[sample_id]][,exprMat_list[[sample_id]] %>% colnames %>% intersect(selected_BC_all)]
	cat(exprMat_list[[sample_id]] %>% ncol," cells found\n")
	cat("Done!\n")
}
exprMat <- do.call("cbind",exprMat_list)[,selected_BC_all]

rm(exprMat_list)


genesKept <- geneFiltering(exprMat %>% as.matrix, scenicOptions) ## by default, minCountsPerGene=3*.01*ncol(exprMat);minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered %>% as.matrix, scenicOptions)
exprMat_filtered_log <- log2(as.matrix(exprMat_filtered)+1) 
runGenie3(exprMat_filtered_log, scenicOptions)


############################## Build and score the GRN ##############################
exprMat_log <- log2(exprMat+1)
## For each species, we used two gene-motif rankings (10kb around the TSS or 500bp upstream the TSS), which determine the search space around the transcription start site.

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
# scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=NULL)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log %>% as.matrix)


############################## Binarize activity ##############################
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


############################## Exploring output ##############################
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org




