## first, find DA peaks:

library(Seurat)
library(Signac)
library(dplyr)
library(tibble)

library(ggplot2)
library(RColorBrewer)


merging_name <- "NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed"
working_dir <- "..../scripts_for_sharing/Fig6_Motif/"
setwd(working_dir)
out_path <- "results/"

obj_joint <- readRDS(paste0("...../scripts_for_sharing/data_preprocessing/results/snRNA_ATAC_jointly_analyzed_for_",merging_name,".rds"))
meta.data <- read.table("...../scripts_for_sharing/SuppTable_barcode_QC.txt",head=TRUE,sep="\t")

DefaultAssay(obj_joint) <- "peaksMACS2"
cell_types <- c("PT","PT(S1)","PT(S2)","PT(S3)")
Ages <- c("E16.5","P0","W3","W12","W52","W92")
DE_all_peaksMACS2_list <-
foreach(i=seq(1,length(cell_types))) %dopar%{
        DE_all <- NULL
        cell_type <- cell_types[i]
        for (Age in Ages){
                BC_1 <- obj_joint@meta.data %>% rownames_to_column %>% filter(Age==!!Age) %>% filter(cell_type==!!cell_type) %>% filter(Sex=="F") %>% .$barcode ## select cells from female
                BC_2 <- obj_joint@meta.data %>% rownames_to_column %>% filter(Age==!!Age) %>% filter(cell_type==!!cell_type) %>% filter(Sex=="M") %>% .$barcode ## select cells from male

                if (length(BC_1)>=10 & length(BC_2)>=10){
                cat("comparing M and F for ",cell_type," for Age ",Age,"...\n")
                DE <- FindMarkers(obj_joint,ident.1=BC_1,ident.2=BC_2,test.use="LR",latent.vars="peak_region_fragments_500MACS2",logfc.threshold = 0.1,min.pct = 0) %>% rownames_to_column("peaks") %>% mutate(ident.1="F",ident.2="M",Age=Age,cell_type=cell_type)
                DE_all <- rbind(DE_all,DE)
                }
        }
        return (DE_all)
}
DE_all_peaksMACS2 <- do.call("rbind",DE_all_peaksMACS2_list)
write.table(DE_all_peaksMACS2,paste0(out_path,"PT_DE_genes_between_M_and_F_peaksMACS2_assay.txt"),sep="\t",quote=FALSE,row.names=FALSE)
rm(obj_joints)



obj_joint <- readRDS(paste0("...../scripts_for_sharing/data_preprocessing/results/",merging_name,"_chromVAR.rds"))
obj_joint@meta.data$cell_type <- column_to_rownames(meta.data,"barcode")[obj_joint@meta.data %>% rownames,"cell_type"]
DefaultAssay(obj_joint) <- "peaksMACS2"

## specify ordered ages
Ages_ordered <- c("E16.5","P0","W3","W12","W52","W92")
cell_types <- c("PT","PT(S1)","PT(S2)","PT(S3)")

## read in DA peaks
PT_DA_peaks <- read.table("results/PT_DE_genes_between_M_and_F_peaksMACS2_assay.txt",head=TRUE,sep="\t")
PT_DA_peaks$trend <- ifelse(PT_DA_peaks$avg_log2FC>0,"F","M")


# find enriched motifs for peaks of interest - set p = 0.05 as the cutoff
enriched.motifs_all <- NULL
	for (Age in Ages_ordered){
		for (cell_type in cell_types){
			for (trend in c("F","M")){
				DA_peaks <- PT_DA_peaks %>% filter(Age==!!Age) %>% filter(cell_type==!!cell_type) %>% filter(p_val_adj<0.05) %>% filter(trend==!!trend) %>% .$peaks %>% as.character
				if (length(DA_peaks)>=10){
					open.peaks <- AccessiblePeaks(obj_joint,cells=obj_joint@meta.data %>% rownames_to_column %>% filter(Age==!!Age) %>% filter(cell_type==!!cell_type) %>% filter(Sex==trend) %>% .$rowname)

					# match the overall GC content in the peak set
					meta.feature <- GetAssayData(obj_joint, assay = "peaksMACS2", slot = "meta.features")
					peaks.matched <- MatchRegionStats(
							meta.feature = meta.feature[open.peaks,],
							query.feature = meta.feature[DA_peaks,],
							n = 50000)

					# test enrichment
					enriched.motifs <- FindMotifs(object = obj_joint,
									features=DA_peaks,
									background = peaks.matched) %>%
							mutate(Age=Age,cell_type=cell_type,trend=trend)

					enriched.motifs_all <- rbind(enriched.motifs_all,enriched.motifs)
			}
		}
	}
}
write.table(enriched.motifs_all,paste0(out_path,"Motif_enrichment_analysis_for_DARs_for_each_age_sex_and_cell_peak_p_val_cutoff_0.05.txt"),sep="\t")

## plotting
enriched.motifs_all <- read.table(paste0(out_path,"Motif_enrichment_analysis_for_DARs_for_each_age_sex_and_cell_peak_p_val_cutoff_0.05.txt"),head=TRUE,sep="\t")
## visualize motifs of interest, using p value 0.05 as the cutoff
for (cell_type in cell_types){
	enriched_motifs <- enriched.motifs_all %>% filter(p.adjust<0.01) %>% filter(fold.enrichment>0) %>% filter(cell_type==!!cell_type) %>% group_by(Age,trend) %>% arrange(p.adjust,desc(fold.enrichment)) %>% slice(1:20)
	motifs_to_visualize <- enriched_motifs %>% filter(cell_type==!!cell_type) %>% filter(motif.name %in% (enriched_motifs$`motif.name` %>% as.character %>% unique)) %>% mutate(Age=Age %>% factor(levels=Ages_ordered))
	if (nrow(motifs_to_visualize)>0){
		p <- ggplot(motifs_to_visualize)+
			geom_point(aes(x=Age,y=motif.name,color=fold.enrichment,size=-log10(p.adjust)),shape=16)+
			scale_color_gradientn(colors=brewer.pal(9,"YlOrRd")[3:9])+
			facet_grid(.~trend)+
			theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+
			ggtitle(cell_type)
		pdf(paste0(out_path,"plotting/enriched_motifs_summary_for_sex_DAR_from_",cell_type,"_peaks_using_0.05_as_cutoff.pdf"),useDingbats=FALSE,width=6,height=14);print(p);dev.off()
	}
}



###### manually select motifs for plotting, for both PT(S2) and PT(S3)
cell_type <- "PT(S2)"
#female_specific_motifs <- c("ZBTB14","TEF","NRF1","KLF15","HLF","HINFP","E2F6","DBP","CTCFL")
female_specific_motifs <- c("CTCFL","DBP","E2F6","HINFP","HLF","HNF1A","HNF1B","KLF14","KLF15","NFIL3","NRF1","SP4","SP9","TCFL5","TEF","TFDP1","ZBTB14","ZNF460")
male_specific_motifs <- c("POU6F1","POU4F3","POU4F2","POU4F1","NR3C2","NR3C1","Lhx3","HOXB3","HNF4G","HNF4A(var.2)","HNF4A","HMBOX1","GBX2","ESX1","Arid3b","Ar","Alx4","Alx1")
female_male_shared_motifs <- NULL
motifs_to_visualize_S2 <- data.frame(motif_name=c(female_specific_motifs,male_specific_motifs,female_male_shared_motifs),
				sex_specificity_annotation=c("female_specific" %>% rep(female_specific_motifs %>% length),
							"male_specific" %>% rep(male_specific_motifs %>% length),
							"female_male_shared" %>% rep(female_male_shared_motifs %>% length)))

cell_type <- "PT(S3)"
female_specific_motifs <- c("TEAD4","TEAD3","TEAD2","TEAD1","POU5F1B","POU3F4","POU3F2","POU2F1","POU1F1","NFIX(var.2)","NFIC(var.2)","NFIC","NFIB","HLF","GCM2","GCM1")
male_specific_motifs <- c("THRB","RXRG","RXRB","Rxra","PPARD","NR3C2","NR3C1","Nr2f6(var.2)","Nr2f6","NR2F1(var.2)","NR1H2::RXRA","HNF4G","HNF4A(var.2)","Ar")
female_male_shared_motifs <- c("HNF4A","HNF1B","HNF1A")
motifs_to_visualize_S3 <- data.frame(motif_name=c(female_specific_motifs,male_specific_motifs,female_male_shared_motifs),
				sex_specificity_annotation=c("female_specific" %>% rep(female_specific_motifs %>% length),
							"male_specific" %>% rep(male_specific_motifs %>% length),
							"female_male_shared" %>% rep(female_male_shared_motifs %>% length)))
motifs_to_visualize_S2_S3 <- rbind(motifs_to_visualize_S2 %>% mutate(cell_type="PT(S2)"),motifs_to_visualize_S3 %>% mutate(cell_type="PT(S3)"))

to_plot <- motifs_to_visualize_S2_S3 %>%
	merge(all_motifs_annotation,by.x="motif_name",by.y="Name",all.x=TRUE) %>%
	merge(enriched.motifs_all %>% filter(Age %in% c("W3","W12","W52","W92")),by.x=c("motif","cell_type"),by.y=c("motif","cell_type"),all.x=TRUE) %>%
	mutate(Age=Age %>% factor(levels=Ages_ordered %>% rev)) %>%
	unique
p <- ggplot(to_plot)+
        geom_point(aes(x=motif_name,y=Age,color=fold.enrichment,size=-log10(pvalue)))+
        scale_color_gradientn(colors=brewer.pal(9,"YlOrRd")[3:9])+
        facet_grid(trend~cell_type+sex_specificity_annotation+Family,scales="free",space="free")+
        theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+
        theme(axis.text.x=element_text(angle=90))
pdf(paste0(out_path,"plotting/manually_highlight_motifs_summary_for_sex_DARs_for_PT_S2_S3_peaks_using_0.05_as_cutoff.pdf"),useDingbats=FALSE,width=20,height=5)
print(p);dev.off()



## plotting for NR motifs:
all_motifs_annotation <- read.table("..../scripts_for_sharing/Fig6_Motif/motif_annotation_JASPAR_vertebrate_core.txt",head=TRUE,sep="\t")
nuclear_receptor_motifs <- all_motifs_annotation %>% filter(grepl("NR",Family))

to_plot_NR_motifs <- enriched.motifs_all %>%
                                filter(motif %in% as.character(nuclear_receptor_motifs$motif %>% unique)) %>%
                                filter(cell_type %in% c("PT(S1)","PT(S2)","PT(S3)")) %>%
                                merge(nuclear_receptor_motifs %>% select(motif,Family),all.x=TRUE) %>%
                                mutate(Age=Age %>% factor(levels=Ages_ordered)) %>%
                                mutate(Family=Family %>% strsplit("\\(|\\)") %>% lapply("[[",2) %>% unlist)
p <- ggplot(to_plot_NR_motifs)+
	geom_point(aes(x=Age,y=motif.name,color=fold.enrichment,size=-log10(p.adjust)),shape=16)+
	scale_color_gradientn(colors=brewer.pal(9,"YlOrRd")[3:9])+
	facet_grid(Family~trend+cell_type,scales="free",space="free")+
	theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+
	theme(axis.text.x=element_text(angle=90))
pdf(paste0(out_path,"NR_motifs_summary_for_sex_DARs_DA_peaks_p_val_cutoff_as_0.05.pdf"),useDingbats=FALSE,height=10)
print(p);dev.off()


########################## plotting for NR1/2/3 only
highlighted_NR_motifs <- to_plot_NR_motifs %>% filter(Age %in% c("W12","W52","W92")) %>% filter(p.adjust<1e-10) %>% filter(fold.enrichment>=1.5) %>% .$motif %>% as.character %>% unique

nuclear_receptor_1_2_3_motifs <- all_motifs_annotation %>% filter(grepl("NR",Family)) %>% filter(grepl(c("NR1","NR2","NR3") %>% paste0(collapse="|"),Family))
to_plot_NR_1_2_3_motifs <- enriched.motifs_all %>%
                                filter(motif %in% as.character(nuclear_receptor_1_2_3_motifs$motif %>% unique)) %>%
                                filter(cell_type %in% c("PT(S1)","PT(S2)","PT(S3)")) %>%
                                merge(nuclear_receptor_1_2_3_motifs %>% select(motif,Family),all.x=TRUE) %>%
                                mutate(Age=Age %>% factor(levels=Ages_ordered)) %>%
                                mutate(Family=Family %>% strsplit("\\(|\\)") %>% lapply("[[",2) %>% unlist) %>%
                                mutate(annotation=ifelse(motif %in% highlighted_NR_motifs,"highlighted","not_highlighted"))
p <- ggplot(to_plot_NR_1_2_3_motifs)+
        geom_point(aes(x=Age,y=motif.name,color=fold.enrichment,size=-log10(p.adjust)),shape=16)+
        scale_color_gradientn(colors=brewer.pal(9,"YlOrRd")[3:9])+
        facet_grid(annotation+Family~trend+cell_type,scales="free",space="free")+
        theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+
        theme(axis.text.x=element_text(angle=90))
pdf(paste0(out_path,"NR_1_2_3_motifs_summary_for_sex_DARs_DA_peaks_p_val_cutoff_as_0.05.pdf"),useDingbats=FALSE,height=10)
print(p);dev.off()

############# only show the highlighted motifs within NR1/2/3, and also the not-highlighted motifs within NR3 (following the script above)
p2 <- p
p2$data <- p2$data %>% filter(annotation=="highlighted" | grepl("NR3",Family))
motif_ordered <- p2$data %>%
                mutate(Family=Family %>% factor(levels=c("NR1","NR2","NR3"))) %>%
                mutate(annotation=annotation %>% factor(levels=c("highlighted","not_highlighted"))) %>%
                arrange(Family,annotation) %>%
                .$`motif.name` %>% unique
p2$data$`motif.name` <- p2$data$`motif.name` %>% factor(levels=motif_ordered)
p2 <- p2 + facet_grid(Family~trend+cell_type,scales="free",space="free")
pdf(paste0(out_path,"NR_1_2_3_highlighted_motifs_summary_plus_ESR_motifs_for_sex_DARs_DA_peaks_p_val_cutoff_as_0.05.pdf"),useDingbats=FALSE,height=6)
print(p2);dev.off()


pdf(paste0(out_path,"highlighted_NR_motif_plus_ESR_motifs_sequences.pdf"),useDingbats=FALSE,width=5,height=40)
MotifPlot_modified(object = obj_joint,motifs=p2$data %>% group_by(Family) %>% arrange(Family,desc(motif)) %>% .$motif %>% unique,ncol=1)
dev.off()

