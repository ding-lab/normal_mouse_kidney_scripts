# this script is used to visualize scenic output in customized ways, including - tSNE plot, heatmap for PT-S1/S3/S3-specific regulons (average binary activity), and also use heatmap to visualize the expression of the downstream targets of sex-biased reuglons


library(SCopeLoomR)
library(Seurat)
library(dplyr)
library(tibble)

library(ggplot2)
library(RColorBrewer)


################################ step 0: set up, loading essential files ################################
merging_name <- "NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed"
subset_name <- "PT_S1_S2_S3_100cells_per_sample"
analysis_name <- "scenic_with_all_genes"
working_dir <- "..../scripts_for_sharing/Fig4_Senic/"
setwd(working_dir)

sample_info <- read.table("..../scripts_for_sharing/data_preprocessing/input/sample_info.txt",head=TRUE,sep="\t") %>% select(sample_id,Age,Sex,pair_name)
source("..../scripts_for_sharing/colorset.R")

################################ step 1: visualize tSNE plot ################################
out_path <- "results_plotting/regulon_visualization_in_tSNE/"
if (!file.exists(out_path)){dir.create(out_path)}

## load embeddings
loom <- open_loom("results/output/scenic.loom")
embeddings <- get_embeddings(loom)

## load meta.data
meta.data_for_selected_BC <- read.table("results/meta.data_for_selected_BC_for_SCENIC_analysis.txt",head=TRUE,sep="\t",row.names=1)

library(cowplot)
tSNE_embedding <- embeddings[["SCENIC t-SNE: AUC 191 regulons (50PCs, 30 perplexity)"]]
colnames(tSNE_embedding) <- c("tSNE_1","tSNE_2")
to_plot <- tSNE_embedding %>% as.data.frame %>% rownames_to_column %>% 
		mutate(sample_id=rowname %>% strsplit("_") %>% lapply("[[",1) %>% unlist) %>% 
		merge(sample_info) %>% 
		merge(meta.data_for_selected_BC %>% rownames_to_column %>% select(rowname,cell_type),all.x=TRUE)
p1 <- ggplot(to_plot)+geom_point(aes(x=tSNE_1,y=tSNE_2,color=Age),shape=16,size=0.5)+theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+theme(aspect.ratio=1)+scale_color_manual(values=col_age[c("E16.5","P0","W3","W12","W52","W92")])+guides(colour = guide_legend(override.aes = list(size=2)))
p2 <- ggplot(to_plot)+geom_point(aes(x=tSNE_1,y=tSNE_2,color=Sex),shape=16,size=0.5)+theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+theme(aspect.ratio=1)+scale_color_manual(values=col_sex)+guides(colour = guide_legend(override.aes = list(size=2)))
p3 <- ggplot(to_plot)+geom_point(aes(x=tSNE_1,y=tSNE_2,color=cell_type),shape=16,size=0.5)+theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+theme(aspect.ratio=1)+scale_color_manual(values=col_cell_type[c("PT(S1)","PT(S2)","PT(S3)")])+guides(colour = guide_legend(override.aes = list(size=2)))
pdf(paste0(out_path,"scenic_tsne_visualization_colored_by_age_sex_and_celltype.pdf"),useDingbats=FALSE,width=20)
plot_grid(p1,p2,p3,ncol=3)
dev.off()

################################ step 2: visualize segment-specific/sex-specific regulons ################################
out_path <- "results_plotting/cell_type_specific_regulons/"
if (!file.exists(out_path)){dir.create(out_path)}

######## statistical test to find the regulons of interest
library(rstatix)
## load binary activity
binaryRegulonActivity <- readRDS("results/int/4.1_binaryRegulonActivity.Rds")
cell_info <- readRDS("results/int/cellInfo.Rds")
binaryRegulonActivity_with_annotation <- binaryRegulonActivity %>% as.data.frame %>%
				rownames_to_column("gene_sets") %>%
				tidyr::pivot_longer(!gene_sets,names_to="barcode",values_to="binary_activity") %>%
				merge(cell_info %>% rownames_to_column("barcode") %>% select(barcode,cell_type,Age,Sex),by="barcode",all.x=TRUE)
## for cell type specific regulons, using binary activities and wilcox test
binary_test_res_PT_all <- NULL
to_exclude <- binaryRegulonActivity_with_annotation %>% group_by(gene_sets,Age) %>% summarise(binary_activity=mean(binary_activity)) %>% filter(binary_activity==0 | binary_activity==1) %>% select(gene_sets,Age)
cell_types <- c("PT(S1)","PT(S2)","PT(S3)")
for (cell_type in cell_types){
    binary_test_res_PT <- binaryRegulonActivity_with_annotation %>%
        anti_join(to_exclude,by=c("gene_sets","Age")) %>%
        mutate(cell_type=plyr::mapvalues(cell_type,from=setdiff(cell_types,cell_type),to=paste0("not_",cell_type) %>% rep(2))) %>%
        group_by(gene_sets,Age) %>%
        wilcox_test(binary_activity~cell_type)
    binary_test_res_PT_all <- rbind(binary_test_res_PT_all,binary_test_res_PT)
}
binaryRegulonActivity_average_scores <- binaryRegulonActivity_with_annotation %>% group_by(gene_sets,cell_type,Age) %>% summarise(binary_activity=mean(binary_activity)) %>%
                                        tidyr::pivot_wider(names_from="cell_type",values_from="binary_activity")
binary_test_res_PT_all_with_score_annotations <- binary_test_res_PT_all %>%
                                        merge(binaryRegulonActivity_average_scores,by=c("Age","gene_sets"),all.x=TRUE)
write.table(binary_test_res_PT_all_with_score_annotations,paste0(out_path,"cell_type_specific_regulons_binary_activity_test_for_each_age_and_cell_type.txt"),sep="\t",quote=FALSE,row.names=FALSE)
rm(binary_test_res_PT_all_with_score_annotations)





## select regulons of interest, then plotting for the regulons of interest, per sample, annotated by cell type, age and sex
## load binary test results
cell_type_specific_binary_regulons_test <- read.table(paste0(out_path,"/cell_type_specific_regulons_binary_activity_test_for_each_age_and_cell_type.txt"),head=TRUE,sep="\t")

## select PT-S1/S2/S3 regulons, combine the list
selected_regulons_for_annotation <- readRDS("results/int/4.3_regulonSelections.Rds")
S1_specific_binary_regulons <- cell_type_specific_binary_regulons_test %>% filter(p<1e-2) %>% filter(Age %in% c("W12","W52","W92")) %>%
			filter(group2=="PT(S1)") %>% filter(`PT.S1.`>`PT.S2.` & `PT.S1.`>`PT.S3.`) %>%
			group_by(gene_sets) %>% count() %>% filter(n==3) %>%
			.$gene_sets %>% intersect(c(selected_regulons_for_annotation$corr,selected_regulons_for_annotation$onePercent))
S2_specific_binary_regulons <- cell_type_specific_binary_regulons_test %>% filter(p<1e-2) %>% filter(Age %in% c("W12","W52","W92")) %>%
			filter(group2=="PT(S2)") %>% filter(`PT.S2.`>`PT.S1.` & `PT.S2.`>`PT.S3.`) %>%
			group_by(gene_sets) %>% count() %>% filter(n==3) %>%
			.$gene_sets %>% intersect(c(selected_regulons_for_annotation$corr,selected_regulons_for_annotation$onePercent))
S3_specific_binary_regulons <- cell_type_specific_binary_regulons_test %>% filter(p<1e-2) %>% filter(Age %in% c("W12","W52","W92")) %>%
			filter(group2=="PT(S3)") %>% filter(`PT.S3.`>`PT.S1.` & `PT.S3.`>`PT.S2.`) %>%
			group_by(gene_sets) %>% count() %>% filter(n==3) %>%
			.$gene_sets %>% intersect(c(selected_regulons_for_annotation$corr,selected_regulons_for_annotation$onePercent))
S1_S2_S3_specific_binary_regulons_annotation <- do.call("rbind",
		list(data.frame(regulon=S1_specific_binary_regulons,cell_type="PT(S1)"),                                                                
			data.frame(regulon=S2_specific_binary_regulons,cell_type="PT(S2)"),                                                                
			data.frame(regulon=S3_specific_binary_regulons,cell_type="PT(S3)"))) %>%
		column_to_rownames("regulon")

## get the binaryRegulonActivity matrix for each sample and cell type
binaryRegulonActivity_long_with_annotations <- binaryRegulonActivity %>% as.data.frame %>% rownames_to_column("regulon") %>%
						tidyr::pivot_longer(!regulon,names_to="barcode",values_to="binary_activity") %>%
						merge(cell_info %>% rownames_to_column("barcode"),by="barcode",all.x=TRUE)
average_binaryRegulonActivity_for_each_sample_and_celltype <- binaryRegulonActivity_long_with_annotations %>%
                                                group_by(sample_id,Age,Sex,cell_type,regulon) %>%
                                                summarise(binary_activity=mean(binary_activity)) %>%
                                                mutate(group=paste0(sample_id,"__",Age,"__",Sex,"__",cell_type)) %>% as.data.frame %>%
                                                select(group,regulon,binary_activity) %>%
                                                tidyr::pivot_wider(names_from="group",values_from="binary_activity") %>%
                                                column_to_rownames("regulon")
column_annotations <- average_binaryRegulonActivity_for_each_sample_and_celltype %>% colnames %>% as.data.frame %>% setNames("group") %>%
                                                tidyr::separate("group",c("sample_id","Age","Sex","cell_type"),"__",remove=FALSE) %>%
                                                mutate(cell_type=cell_type %>% factor(levels=c("PT(S1)","PT(S2)","PT(S3)"))) %>%
                                                mutate(Age=Age %>% factor(levels=c("E16.5","P0","W3","W12","W52","W92"))) %>%
                                                mutate(Sex=Sex %>% factor(levels=c("F","M"))) %>%
                                                arrange(cell_type,Sex,Age) %>% mutate(sample_id=NULL) %>%
                                                column_to_rownames("group")

## a few steps to re-order the regulons
x <- average_binaryRegulonActivity_for_each_sample_and_celltype %>% dist %>% hclust
regulon_ordered <- x$labels[x$order]
regulon_ordered <- S1_S2_S3_specific_binary_regulons_annotation %>% rownames_to_column %>%
                mutate(rowname=rowname %>% factor(levels=regulon_ordered)) %>%
                arrange(cell_type,rowname) %>%
                .$rowname %>% as.character

## for S3-specific regulons, do another round of re-ordering (then we will see the pattern of sex-biased regulons)
S3_specific_binary_regulons_average_binaryActivity_for_each_sample_in_PT_S3 <- average_binaryRegulonActivity_for_each_sample_and_celltype[S3_specific_binary_regulons,column_annotations %>% rownames_to_column %>% filter(cell_type=="PT(S3)") %>% .$rowname]
x <- S3_specific_binary_regulons_average_binaryActivity_for_each_sample_in_PT_S3 %>% dist %>% hclust
regulon_ordered_2_S3 <- x$labels[x$order]
regulon_ordered_2_S3 <- x$labels[x$order]

regulon_ordered_2 <- c(S1_S2_S3_specific_binary_regulons_annotation %>% rownames_to_column %>% ## modify the order of regulons again (due to S3 regulon reordering)
			filter(cell_type!="PT(S3)") %>%
			mutate(rowname=rowname %>% factor(levels=regulon_ordered)) %>%
			arrange(cell_type,rowname) %>%
			.$rowname %>% as.character,
		regulon_ordered_2_S3)

## for a more clear visualization, we only keep the regulons whose max average binary activity (across samples and cell types) are > 0.4
average_binary_regulon_activity_threshold <- 0.4
groups_adult <- colnames(average_binaryRegulonActivity_for_each_sample_and_celltype)[grepl("W12|W52|W92",average_binaryRegulonActivity_for_each_sample_and_celltype %>% colnames)]

regulons_to_select_2 <- regulon_ordered_2 %>% intersect(which(apply(average_binaryRegulonActivity_for_each_sample_and_celltype[,groups_adult],1,max)>=average_binary_regulon_activity_threshold) %>% names)
pdf(paste0(out_path,"PT_segment_specific_regulons_based_on_binary_activity_test_average_binary_activity_heatmap_for_each_sample_minimum_regulon_activity_as_",average_binary_regulon_activity_threshold,"_regulon_ordered_2.pdf"),useDingbats=FALSE)
NMF::aheatmap(average_binaryRegulonActivity_for_each_sample_and_celltype[regulons_to_select_2,column_annotations %>% rownames],
		annCol=column_annotations %>% select(cell_type,Age,Sex),annRow=S1_S2_S3_specific_binary_regulons_annotation[regulons_to_select_2,],
		annColor=colVars,main="average binary activity",
		Colv=NA,Rowv=NA)
dev.off()



################################ step 3: visualize some sex-biased regulons and its downstream targets ################################



################################ step 4: downstream of STAT5 and BCL6 ################################

