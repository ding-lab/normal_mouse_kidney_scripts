############## This script is used for plotting downstream targets of selected sex-biased regulons ##############

#library(SCENIC)
library(AUCell)
library(dplyr)
library(tibble)
library(ggplot2)

## setup
merging_name <- "NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed"
subset_name <- "PT_S1_S2_S3_100cells_per_sample"
analysis_name <- "scenic_with_all_genes"
working_dir <- "..../scripts_for_sharing/Fig4_Senic/"
setwd(working_dir)
out_path <- "results_plotting/sex_specific_regulons_visualization/"
if (!file.exists(out_path)){dir.create(out_path)}


## load selected regulons
selected_regulons_for_annotation <- readRDS("results/int/4.3_regulonSelections.Rds")
regulonAUC <- readRDS("results/int/3.4_regulonAUC.Rds")
regulonAUC_scores <- getAUC(regulonAUC)

binaryRegulonActivity <- readRDS("results/int/4.1_binaryRegulonActivity.Rds")
regulons_asGeneSet <- readRDS("results/int/2.6_regulons_asGeneSet.Rds")

cellInfo <- readRDS("results/int/cellInfo.Rds")
colVars <- readRDS("results/int/colVars.Rds");colVars$cell_type <- colVars$cell_type

AUCellThresholds_Info <- read.table("results/int/3.5_AUCellThresholds_Info.tsv",head=TRUE,sep="\t")


## bianry activities with annotations
average_binary_activity_for_each_sample_and_celltype <- binaryRegulonActivity %>% as.data.frame %>% rownames_to_column("regulon") %>% tidyr::pivot_longer(!regulon,names_to="barcode",values_to="binary_activity") %>% merge(cellInfo %>% rownames_to_column("barcode") %>% select(barcode,Age,Sex,sample_id,cell_type)) %>% group_by(regulon,Age,Sex,cell_type,sample_id) %>% summarize(binary_activity=mean(binary_activity))


## heatmap for sex specific regulons
library(ComplexHeatmap)
source("..../scripts_for_sharing/colorset.R")
color_schemes_list <- list("Sex"=col_sex,"Age"=col_age,"platform"=col_platform,"cell_type"=col_cell_type)


############## plotting for regulon scores
to_plot.1 <- regulonAUC_scores %>% as.data.frame %>% rownames_to_column("regulon") %>% 
		tidyr::pivot_longer(!regulon,names_to="barcode",values_to="binary_activity") %>% 
		merge(cellInfo %>% rownames_to_column("barcode"),by="barcode",all.x=TRUE)
ViolinPlot_for_regulon_of_interest <- function(to_plot,regulon,cell_types,Ages,AUCellThresholds_Info){
	to_plot <- to_plot %>% filter(regulon==!!regulon) %>% 
			filter(cell_type %in% cell_types) %>% mutate(cell_type=cell_type %>% factor(levels=cell_types)) %>% 
			filter(Age %in% Ages) %>% mutate(Age=Age %>% factor(levels=Ages))
	threshold <- AUCellThresholds_Info %>% filter(regulon==!!regulon) %>% .$threshold
	p <- ggplot(to_plot)+
		geom_violin(aes(x=sample_id,y=binary_activity),fill=NA)+
		geom_jitter(aes(x=sample_id,y=binary_activity,color=Sex),height=0,width=0.2,shape=16)+
		geom_hline(yintercept=threshold,linetype="dotted")+
		scale_color_manual(values=c("F"="red","M"="blue"))+
#		scale_fill_manual(values=c("F"="red","M"="blue"))+
		facet_grid(.~cell_type+Sex+Age,space="free",scales="free")+
		theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+
		theme(axis.text.x=element_text(angle=90))+
		ggtitle(regulon)
	return (p)
}

BoxPlot_for_regulon_of_interest <- function(to_plot,regulon,cell_types,Ages){
	to_plot <- to_plot %>% filter(regulon==!!regulon) %>%
		filter(cell_type %in% cell_types) %>% mutate(cell_type=cell_type %>% factor(levels=cell_types)) %>%
		filter(Age %in% Ages) %>% mutate(Age=Age %>% factor(levels=Ages))
	p <- ggplot(to_plot)+
		geom_boxplot(aes(x=sample_id,y=binary_activity,fill=Sex))+
		#scale_color_manual(values=c("F"="red","M"="blue"))+
		scale_fill_manual(values=c("F"="red","M"="blue"))+
		facet_grid(.~cell_type+Sex+Age,space="free",scales="free")+
		theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+
		theme(axis.text.x=element_text(angle=90))+
		ggtitle(regulon)
	return (p)
}

############## load expression for plotting (for downstream gene expression)
sample_info <- cellInfo %>% select(sample_id,Age,Sex) %>% unique
average_expression_for_each_sample_and_celltype <- readRDS(paste0("/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/new_analysis_20240608/mouse/sample_integration_RNA_ATAC_merged/",merging_name,"/analysis/average_expression_calculation/results/Average_Expression_for_each_Sample_and_Cell_20240727.rds"))

average_expression_for_each_sample_and_celltype_mtx <- average_expression_for_each_sample_and_celltype$RNA %>% as.data.frame %>% 
						rownames_to_column("gene_symbol") %>% 
						tidyr::pivot_longer(!gene_symbol,names_to="group",values_to="expression") %>% 
						tidyr::separate("group",c("sample_id","cell_type"),"__",remove=FALSE) %>% 
						merge(cellInfo %>% select(sample_id,Age,Sex) %>% unique,by="sample_id",all.x=TRUE) %>% 
						mutate(group=paste0(sample_id,"__",Age,"__",Sex,"__",cell_type)) %>% 
						mutate(Age=Age %>% factor(levels=c("E16.5","P0","W3","W12","W52","W92"))) %>% 
						mutate(cell_type=cell_type %>% factor(levels=c("PT(S1)","PT(S2)","PT(S3)"))) %>% 
						filter(!is.na(Age)) %>% filter(!is.na(cell_type)) %>% 
						arrange(cell_type,Age,Sex) %>% 
						select(gene_symbol,group,expression) %>% 
						tidyr::pivot_wider(names_from="group",values_from="expression") %>% 
						column_to_rownames("gene_symbol")

column_annotation <- average_expression_for_each_sample_and_celltype_mtx %>% colnames %>% as.data.frame %>% setNames("group") %>% 
						tidyr::separate("group",c("sample_id","Age","Sex","cell_type"),"__",remove=FALSE)

get_regulon_downstream_targets <- function(regulons_asGeneSet,regulon){
	regulon <- regulon %>% strsplit(" ") %>% lapply("[[",1) %>% unlist
	targets <- regulons_asGeneSet[[regulon]]
	return (targets)
}

get_heatmap_mtx_and_annotations <- function(expression_mtx,sample_info,regulon,cell_types,Ages){
	genes <- get_regulon_downstream_targets(regulons_asGeneSet=regulons_asGeneSet,regulon=regulon)
	to_plot_mtx <- expression_mtx %>% as.data.frame %>% rownames_to_column("gene_symbol") %>% 
			filter(gene_symbol %in% genes) %>% 
			tidyr::pivot_longer(!gene_symbol,names_to="group",values_to="expression") %>% 
			tidyr::separate("group",c("sample_id","cell_type"),"__",remove=FALSE) %>%
			merge(sample_info,by="sample_id",all.x=TRUE) %>%
			mutate(group=paste0(sample_id,"__",Age,"__",Sex,"__",cell_type)) %>%
			mutate(Age=Age %>% factor(levels=Ages)) %>%
			mutate(cell_type=cell_type %>% factor(levels=cell_types)) %>%
			mutate(Sex=Sex %>% factor(levels=c("F","M"))) %>% 
			filter(!is.na(Age)) %>% filter(!is.na(cell_type)) %>%
			arrange(cell_type,Sex,Age) %>%
			select(gene_symbol,group,expression) %>%
			tidyr::pivot_wider(names_from="group",values_from="expression") %>%
			column_to_rownames("gene_symbol")
	column_annotation <- to_plot_mtx %>% colnames %>% as.data.frame %>% setNames("group") %>%
			tidyr::separate("group",c("sample_id","Age","Sex","cell_type"),"__",remove=FALSE) %>% 
			select(Age,Sex,cell_type,group) %>% 
			column_to_rownames("group")
	column_annotation <- column_annotation[colnames(to_plot_mtx),]
	to_return <- list("exprs_mtx"=to_plot_mtx,"column_annotations"=column_annotation)
	return (to_return)
}


############## plot only selected genes downstream of a regulon based on correlation ##############
Heatmap_for_top_targets_correlated_with_regulon <- function(regulon_activity=average_binary_activity_for_each_sample_and_celltype,
							expression=average_expression_for_each_sample_and_celltype$RNA,
							sample_info=sample_info,
							regulon,cell_type,TF=NULL,other.highlighted.genes=NULL,top_n_correlated=15){
		regulon_activity <- regulon_activity %>% 
					filter(cell_type==cell_type) %>% 
					filter(grepl(paste0("^",!!regulon," "),regulon)) %>% 
					as.data.frame %>% 
					select(sample_id,binary_activity) %>% 
					column_to_rownames("sample_id")
		all_downstream_genes <- get_regulon_downstream_targets(regulons_asGeneSet=regulons_asGeneSet,regulon=regulon)
		mtx_for_all_genes <- average_expression_for_each_sample_and_celltype$RNA[c(all_downstream_genes,TF),] %>% 
					as.data.frame %>% 
					rownames_to_column("gene_symbol") %>% 
					tidyr::pivot_longer(!gene_symbol,names_to="group",values_to="expression") %>% 
					tidyr::separate("group",c("sample_id","cell"),"__") %>% 
					filter(cell==cell_type) %>% mutate(cell=NULL) %>% 
					tidyr::pivot_wider(names_from="gene_symbol",values_from="expression") %>% 
					column_to_rownames("sample_id")
		mtx_for_all_genes <- mtx_for_all_genes[rownames(regulon_activity),]

		top_correlated_genes <- cor(regulon_activity,mtx_for_all_genes) %>% t %>% as.data.frame %>% 
					rownames_to_column("gene_symbol") %>% 
					setNames(c("gene_symbol","correlation")) %>% 
					arrange(desc(correlation)) %>% 
					filter(!gene_symbol %in% c(TF,other.highlighted.genes)) %>% 
					arrange(desc(correlation)) %>% head(top_n_correlated) %>% .$gene_symbol
		
		genes_to_plot <- c(top_correlated_genes,TF,other.highlighted.genes) %>% unique %>% setdiff(paste0(TF,".1"))
		
		Ages_ordered <- c("E16.5","P0","W3","W12","W52","W92")
		samples_ordered <- sample_info %>% mutate(Age=Age %>% factor(levels=Ages_ordered)) %>% arrange(Sex,Age) %>% .$sample_id %>% as.character
		to_plot_mtx <- mtx_for_all_genes[samples_ordered,genes_to_plot] %>% t
		columns_annotation <- sample_info %>% rownames_to_column %>% mutate(rowname=NULL) %>% mutate(sample_id=sample_id %>% factor(levels=samples_ordered)) %>% arrange(sample_id) %>% column_to_rownames("sample_id")

		to_return <- list()
		to_return$mtx <- to_plot_mtx
		to_return$columns_annotation <- columns_annotation
		return (to_return)
}



######### Female regulons
load("/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/new_analysis_20240608/mouse/sample_integration_RNA_ATAC_merged/NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed/analysis/average_expression_calculation/results/expression_plotting/PT_DEG_heatmap/DEGs_res.rda")

M_biased <- DEGs_res %>% filter(trend=="M" & cell_type=="PT(S3)") %>% .$gene_symbol %>% as.character %>% unique
F_biased <- DEGs_res %>% filter(trend=="F" & cell_type=="PT(S3)") %>% .$gene_symbol %>% as.character %>% unique

regulon <- "Cebpd";cell_type <- "PT(S3)";TF <- "Cebpd";
"Socs2" %in% get_regulon_downstream_targets(regulons_asGeneSet,regulon="Cebpd")  # check if Socs2 is downstream of Cebpd - yes - we will manually add Socs2 in the plot for highlight
other.highlighted.genes <- "Socs2"
to_plot <- Heatmap_for_top_targets_correlated_with_regulon(regulon_activity=average_binary_activity_for_each_sample_and_celltype,expression=average_expression_for_each_sample_and_celltype$RNA,sample_info=sample_info,
		regulon=regulon,cell_type=cell_type,TF=TF,other.highlighted.genes=other.highlighted.genes,top_n_correlated=15)
column_ha <- HeatmapAnnotation(Sex=to_plot$columns_annotation$Sex,Age=to_plot$columns_annotation$Age,col = color_schemes_list)
pdf(paste0(out_path,"selected_downstream_targets_expression_for_regulon_",regulon,".pdf"),useDingbats=FALSE)
p <- Heatmap(to_plot$mtx %>% t %>% scale %>% t,top_annotation=column_ha,cluster_columns=FALSE);p %>% print
dev.off()
to_plot$mtx %>% rownames %>% intersect(F_biased)

#[1] "Egf"     "Cyfip2"  "Zbtb44"  "Trim2"   "Glul"    "Jak2"    "Gramd1b"
#[8] "Cebpd"

p1 <- ViolinPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"),AUCellThresholds_Info=AUCellThresholds_Info)
p2 <- BoxPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"))
pdf(paste0(out_path,"sex_specific_regulon_high_in_F_",regulon,"_overview_in_violinPlot.pdf"),useDingbats=FALSE,width=10)
print(p1);print(p2);dev.off()



regulon <- "Creb3l1";cell_type <- "PT(S3)";TF <- "Creb3l1"
"Socs2" %in% get_regulon_downstream_targets(regulons_asGeneSet,regulon="Creb3l1")  # check if Socs2 is downstream of Creb3l1
"Jak2" %in% get_regulon_downstream_targets(regulons_asGeneSet,regulon="Creb3l1")  # check if Jak2 is downstream of Creb3l1
other.highlighted.genes <- "Jak2"
to_plot <- Heatmap_for_top_targets_correlated_with_regulon(regulon_activity=average_binary_activity_for_each_sample_and_celltype,expression=average_expression_for_each_sample_and_celltype$RNA,sample_info=sample_info,
		regulon=regulon,cell_type=cell_type,TF=TF,other.highlighted.genes=other.highlighted.genes,top_n_correlated=15)
column_ha <- HeatmapAnnotation(Sex=to_plot$columns_annotation$Sex,Age=to_plot$columns_annotation$Age,col = color_schemes_list)
pdf(paste0(out_path,"selected_downstream_targets_expression_for_regulon_",regulon,".pdf"),useDingbats=FALSE)
p <- Heatmap(to_plot$mtx %>% t %>% scale %>% t,top_annotation=column_ha,cluster_columns=FALSE);p %>% print
dev.off()
to_plot$mtx %>% rownames %>% intersect(F_biased)
#[1] "Elovl6" "Elovl2" "Sfxn2"  "Stk32b" "Itgb6"  "Grhl2"  "Zfhx3"  "Dsg2"
#[9] "Cdh6"

p1 <- ViolinPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"),AUCellThresholds_Info=AUCellThresholds_Info)
p2 <- BoxPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"))
pdf(paste0(out_path,"sex_specific_regulon_high_in_F_",regulon,"_overview_in_violinPlot.pdf"),useDingbats=FALSE,width=10)
print(p1);print(p2);dev.off()



regulon <- "Foxq1";cell_type <- "PT(S3)";TF <- "Foxq1"
"Socs2" %in% get_regulon_downstream_targets(regulons_asGeneSet,regulon="Foxq1")  # check if Socs2 is downstream of Foxq1
other.highlighted.genes <- "Socs2"
to_plot <- Heatmap_for_top_targets_correlated_with_regulon(regulon_activity=average_binary_activity_for_each_sample_and_celltype,expression=average_expression_for_each_sample_and_celltype$RNA,sample_info=sample_info,
		regulon=regulon,cell_type=cell_type,TF=TF,other.highlighted.genes=other.highlighted.genes,top_n_correlated=15)
column_ha <- HeatmapAnnotation(Sex=to_plot$columns_annotation$Sex,Age=to_plot$columns_annotation$Age,col = color_schemes_list)
pdf(paste0(out_path,"selected_downstream_targets_expression_for_regulon_",regulon,".pdf"),useDingbats=FALSE)
p <- Heatmap(to_plot$mtx %>% t %>% scale %>% t,top_annotation=column_ha,cluster_columns=FALSE);p %>% print
dev.off()
to_plot$mtx %>% rownames %>% intersect(F_biased)
# [1] "Fam107a" "Il17re"  "Zbtb44"  "Dclk2"   "Ston2"   "Foxq1"

p1 <- ViolinPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"),AUCellThresholds_Info=AUCellThresholds_Info)
p2 <- BoxPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"))
pdf(paste0(out_path,"sex_specific_regulon_high_in_F_",regulon,"_overview_in_violinPlot.pdf"),useDingbats=FALSE,width=10)
print(p1);print(p2);dev.off()



######### Male regulons
regulon <- "Bach2_extended";cell_type <- "PT(S3)";TF <- "Bach2";other.highlighted.genes <- NULL
to_plot <- Heatmap_for_top_targets_correlated_with_regulon(regulon_activity=average_binary_activity_for_each_sample_and_celltype,expression=average_expression_for_each_sample_and_celltype$RNA,sample_info=sample_info,
		regulon=regulon,cell_type=cell_type,TF=TF,other.highlighted.genes=other.highlighted.genes,top_n_correlated=15)
column_ha <- HeatmapAnnotation(Sex=to_plot$columns_annotation$Sex,Age=to_plot$columns_annotation$Age,col = color_schemes_list)
pdf(paste0(out_path,"selected_downstream_targets_expression_for_regulon_",regulon,".pdf"),useDingbats=FALSE)
p <- Heatmap(to_plot$mtx %>% t %>% scale %>% t,top_annotation=column_ha,cluster_columns=FALSE);p %>% print
dev.off()
to_plot$mtx %>% rownames %>% intersect(M_biased)
#[1] "Zbtb20"  "Nr1h4"   "Smarca2" "Chst11"  "Errfi1"  "Col19a1" "Uty"
#[8] "Dgkg"    "Cadm1"

p1 <- ViolinPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"),AUCellThresholds_Info=AUCellThresholds_Info)
p2 <- BoxPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"))
pdf(paste0(out_path,"sex_specific_regulon_high_in_F_",regulon,"_overview_in_violinPlot.pdf"),useDingbats=FALSE,width=10)
print(p1);print(p2);dev.off()


regulon <- "Bcl6";cell_type <- "PT(S3)";TF <- "Bcl6";other.highlighted.genes <- NULL
to_plot <- Heatmap_for_top_targets_correlated_with_regulon(regulon_activity=average_binary_activity_for_each_sample_and_celltype,expression=average_expression_for_each_sample_and_celltype$RNA,sample_info=sample_info,
		regulon=regulon,cell_type=cell_type,TF=TF,other.highlighted.genes=other.highlighted.genes,top_n_correlated=15)
column_ha <- HeatmapAnnotation(Sex=to_plot$columns_annotation$Sex,Age=to_plot$columns_annotation$Age,col = color_schemes_list)
pdf(paste0(out_path,"selected_downstream_targets_expression_for_regulon_",regulon,".pdf"),useDingbats=FALSE)
p <- Heatmap(to_plot$mtx %>% t %>% scale %>% t,top_annotation=column_ha,cluster_columns=FALSE);p %>% print
dev.off()
to_plot$mtx %>% rownames %>% intersect(M_biased)
# [1] "Chrm3"   "Chst7"   "Runx1t1" "Mcc"     "Mgmt"    "Ctnna2"  "Zbtb20"
# [8] "Ankrd6"  "Ptn"     "C1qtnf3" "Tc2n"    "Chst11"  "Abcc2"   "Mpped1"
#[15] "Bcl6"

p1 <- ViolinPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"),AUCellThresholds_Info=AUCellThresholds_Info)
p2 <- BoxPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"))
pdf(paste0(out_path,"sex_specific_regulon_high_in_F_",regulon,"_overview_in_violinPlot.pdf"),useDingbats=FALSE,width=10)
print(p1);print(p2);dev.off()


regulon <- "Zbtb20";cell_type <- "PT(S3)";TF <- "Zbtb20";other.highlighted.genes <- NULL
to_plot <- Heatmap_for_top_targets_correlated_with_regulon(regulon_activity=average_binary_activity_for_each_sample_and_celltype,expression=average_expression_for_each_sample_and_celltype$RNA,sample_info=sample_info,
                regulon=regulon,cell_type=cell_type,TF=TF,other.highlighted.genes=other.highlighted.genes,top_n_correlated=15)
column_ha <- HeatmapAnnotation(Sex=to_plot$columns_annotation$Sex,Age=to_plot$columns_annotation$Age,col = color_schemes_list)
pdf(paste0(out_path,"selected_downstream_targets_expression_for_regulon_",regulon,".pdf"),useDingbats=FALSE)
p <- Heatmap(to_plot$mtx %>% t %>% scale %>% t,top_annotation=column_ha,cluster_columns=FALSE);p %>% print
dev.off()
to_plot$mtx %>% rownames %>% intersect(M_biased)
# [1] "Ctnna2"     "Chst11"     "Pdzd2"      "Cntnap5a"   "Nr1h4"
# [6] "Csgalnact1" "Col19a1"    "Tmem178"    "Slco3a1"    "Arhgap42"
#[11] "Dgkg"       "Cadm1"      "Plxdc2"     "Zbtb20"

p1 <- ViolinPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"),AUCellThresholds_Info=AUCellThresholds_Info)
p2 <- BoxPlot_for_regulon_of_interest(to_plot=to_plot.1,regulon=regulon,cell_types="PT(S3)",Ages=c("E16.5","P0","W3","W12","W52","W92"))
pdf(paste0(out_path,"sex_specific_regulon_high_in_F_",regulon,"_overview_in_violinPlot.pdf"),useDingbats=FALSE,width=10)
print(p1);print(p2);dev.off()

