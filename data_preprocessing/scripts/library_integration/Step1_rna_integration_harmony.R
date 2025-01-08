library(dplyr)
library(tibble)
library(Seurat)
library(Signac)
library(data.table)
library(cowplot)
library(ggplot2)


merging_name <- "NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed"
working_dir <- "......"
setwd(working_dir)

sample_info_summary <- read.table("...../scripts_for_sharing/data_preprocessing/input/sample_info.txt",head=TRUE,sep="\t")
out_path <- "results/"


#================= Load datasets and create Seurat objects =================
n_sample <- sample_info_summary %>% nrow

## load the cell annotation sheet from Simon - onlly the column where the cell type (cell_type_harmony_v1) has a value (not Unknown) will be used here
cell_annotation_sheet <- read.table("...../scripts_for_sharing/SuppTable_barcode_QC.txt",head=TRUE,sep="\t")
sample_ids <- c()
obj_list <- list()
for (i in 1:n_sample){
	sample_id<-sample_info_summary[i,"sample_id"] %>% as.character
	tmp<-readRDS(sample_info_summary[i,"Path"] %>% as.character)

	BC_to_extract <- cell_annotation_sheet %>%  
			filter(sample_id==!!sample_id) %>% 
			.$barcode %>% gsub(paste0(sample_id,"_"),"",.)
	tmp <- tmp %>% subset(cells=BC_to_extract)

	tmp$sample_id<-sample_id
	tmp$original_barcode<-rownames(tmp@meta.data)
	DefaultAssay(tmp) <- "RNA"
	tmp$SCT <- NULL
	tmp$ATAC_MACS2 <- NULL

	obj_list[[sample_id]]<-tmp
	sample_ids<-c(sample_ids,sample_id)
}



#================= Merge =================
#To merge more than two Seurat objects, simply pass a vector of multiple Seurat objects to the y parameter for merge.
merged <- merge(obj_list[[sample_ids[1]]],y=unlist(obj_list)[-1],
		add.cell.ids = sample_ids,
		project = merging_name)

merged@meta.data$Age <- plyr::mapvalues(merged@meta.data$sample_id,from=sample_info_summary$sample_id %>% as.character,to=sample_info_summary$Age %>% as.character)
merged@meta.data$Sex <- plyr::mapvalues(merged@meta.data$sample_id,from=sample_info_summary$sample_id %>% as.character,to=sample_info_summary$Sex %>% as.character)
merged@meta.data$pair_name <- plyr::mapvalues(merged@meta.data$sample_id,from=sample_info_summary$sample_id %>% as.character,to=sample_info_summary$pair_name %>% as.character)
merged@meta.data$Nuclei_Dissociation_Date <- plyr::mapvalues(merged@meta.data$sample_id,from=sample_info_summary$sample_id %>% as.character,to=sample_info_summary$Nuclei_Dissociation_Date %>% as.character)

merged


#subsetting high quality cells
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(merged)
merged <- ScaleData(merged, features = all.genes)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))

#============= Run Harmony =============
library(harmony)
library(RColorBrewer)

source("...../scripts_for_sharing/colorset.R")

#Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
pdf(file = paste0(out_path,"Harmony.pdf"), height=2.5, width=6)
merged <- merged %>% RunHarmony("Nuclei_Dissociation_Date", plot_convergence = TRUE)
dev.off()

harmony_embeddings <- Embeddings(merged, 'harmony')
harmony_embeddings[1:5, 1:5]

#============= Downstream analysis =============
merged <- merged %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = seq(0.5,1,0.1), save.SNN=TRUE) 
  #%>% identity()

saveRDS(merged,paste0(out_path,"Harmony_integration_for_",merging_name,".rds"))

#merged_unknown_removed <- merged %>% subset(cells=merged %>% colnames %>% setdiff(BC_unknown))
p0_1 <- DimPlot(merged,group.by="RNA_snn_res.0.5",label=TRUE)+ggtitle("RNA_snn_res.0.5") + theme(aspect.ratio=1)
p0_2 <- DimPlot(merged,group.by="Age")+ggtitle("Age") + scale_color_manual(values=brewer.pal(6,"Set3") %>% setNames(c("E16.5","P0","W3","W12","W52","W92")))+theme(aspect.ratio=1)

to_plot <- p0_1$data %>% rownames_to_column %>% 
		merge(merged@meta.data %>% rownames_to_column %>% select(rowname,Age,Sex,pair_name,sample_id)) %>% 
		mutate(Age=Age %>% factor(levels=c("E16.5","P0","W3","W12","W52","W92")))
p1 <- p0_1;p1$data <- to_plot
p1 <- p1 + facet_grid(Sex~Age) + theme(aspect.ratio=1)

merged@meta.data$cell_type <- plyr::mapvalues(merged@meta.data %>% rownames,from=cell_annotation_sheet$barcode,to=cell_annotation_sheet$cell_type)

p2 <- DimPlot(merged,group.by="cell_type",label=TRUE)+ggtitle("cell_type") + theme(aspect.ratio=1) #+ scale_color_manual(values=col_cell_type)
p2_1 <- DimPlot(merged,group.by="cell_type",label=FALSE)+ggtitle("cell_type") + theme(aspect.ratio=1) #+ scale_color_manual(values=col_cell_type)
to_plot <- p2$data %>% rownames_to_column %>% 
		merge(merged@meta.data %>% rownames_to_column %>% select(rowname,Age,Sex,pair_name,sample_id)) %>%
		mutate(Age=Age %>% factor(levels=c("E16.5","P0","W3","W12","W52","W92")))
p3 <- p2_1;p3$data <- to_plot
p3 <- p3 + facet_grid(Sex~Age) + theme(aspect.ratio=1)

p3_1 <- p3;p3_1$data <- p3_1$data %>% filter(Age %in% c("E16.5","P0"));p3_1 <- p3_1 + facet_grid(Sex~pair_name) + ggtitle("E16.5/P0")
p4 <- p3;p4$data <- p3$data %>% filter(Age=="E16.5");p4 <- p4 + facet_grid(Sex~pair_name) + ggtitle("E16.5")
p5 <- p3;p5$data <- p3$data %>% filter(Age=="P0");p5 <- p5 + facet_grid(Sex~pair_name) + ggtitle("P0")
p6 <- p3;p6$data <- p3$data %>% filter(Age=="W3");p6 <- p6 + facet_grid(Sex~pair_name) + ggtitle("W3")
p7 <- p3;p7$data <- p3$data %>% filter(Age=="W12");p7 <- p7 + facet_grid(Sex~pair_name) + ggtitle("W12")
p8 <- p3;p8$data <- p3$data %>% filter(Age=="W52");p8 <- p8 + facet_grid(Sex~pair_name) + ggtitle("W52")
p9 <- p3;p9$data <- p3$data %>% filter(Age=="W92");p9 <- p9 + facet_grid(Sex~pair_name) + ggtitle("W92")


pdf(paste0(out_path,"UMAP_for_cell_distribution_overview.pdf"),useDingbats=FALSE,width=35,height=15)
print(p0_1);print(p0_2);print(p1);print(p2);print(p3);
plot_grid(p2,p0_2,ncol=2) %>% print()
print(p3_1)
print(p4);print(p5);print(p6);print(p7);print(p8);print(p9)
dev.off()


