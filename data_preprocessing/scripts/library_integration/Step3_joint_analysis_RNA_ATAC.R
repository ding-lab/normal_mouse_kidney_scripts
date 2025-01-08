library(Seurat)
library(Signac)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(clustree)

merging_name <- "NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed"
out_path <- "results/"

obj_RNA <- readRDS(paste0(out_path,"Harmony_integration_for_",merging_name,".rds"))
obj_ATAC <- readRDS(paste0(output_path,"/peaks/",merging_name,"_snATAC_Merged_BasedOnSelectedPeaks_Normalized_With_GeneActivity.rds.gz"))


obj_joint <- obj_ATAC
obj_joint[["RNA"]] <- CreateAssayObject(counts=GetAssayData(obj_RNA,assay="RNA"))
#obj_joint[["SCT"]] <- CreateAssayObject(data=GetAssayData(obj_RNA,assay="SCT"))
obj_joint[["umap"]] <- NULL

DefaultAssay(obj_joint) <- "RNA"
obj_joint[["rna_pca"]] <- obj_RNA[["pca"]]
obj_joint[["rna_harmony"]] <- obj_RNA[["harmony"]]
obj_joint[["rna.umap"]] <- obj_RNA[["umap"]]

DefaultAssay(obj_joint) <- "peaksMACS2"
obj_joint[["atac.umap"]]  <- obj_ATAC[["umap"]]

for (res in seq(0.5,1,0.1)){
	obj_joint@meta.data[,paste0("RNA_snn_res.",res)] <- obj_RNA@meta.data[obj_joint@meta.data %>% rownames,paste0("RNA_snn_res.",res)]
}

####Now we ready to perform joint Dim-reduction:
# Joint clustering
message("[Joint] cluster cells based on joint RNA and ATAC...")
DefaultAssay(obj_joint) = "peaksMACS2"
# build a joint neighbor graph using both assays
obj_joint <- FindMultiModalNeighbors(object = obj_joint,reduction.list = list("rna_harmony", "harmony"),
		#dims.list = list(1:20, 1:50),verbose = TRUE)
		dims.list = list(1:50, 1:50),verbose = TRUE)

# build a joint UMAP visualization
obj_joint <- RunUMAP(object = obj_joint,nn.name = "weighted.nn", #weighted nearest neighbor
		reduction.name = "wnn.umap",reduction.key = "wnnUMAP_",verbose = TRUE)

####Run this later, once you are happy with the results:
obj_joint <- FindClusters(obj_joint,graph.name = "wsnn",algorithm = 3,  # SLM
		resolution = seq(0.5,1,0.1),
		verbose = TRUE
)



###Now make plots:
p1 = DimPlot(obj_joint, label = TRUE, repel = TRUE, reduction = "wnn.umap") + NoLegend() + ggtitle('Joint WNN') + theme(aspect.ratio=1)
p2 = DimPlot(obj_joint, reduction = 'rna.umap', label = TRUE,repel = TRUE, label.size = 2.5) +NoLegend() + ggtitle('RNA') + theme(aspect.ratio=1)
p3 = DimPlot(obj_joint, reduction = 'atac.umap', label = TRUE,repel = TRUE, label.size = 2.5) +NoLegend() +ggtitle('ATAC') + theme(aspect.ratio=1)
p_wsnn = (p1 + p2 + p3) + plot_annotation(title="wsnn based cluster")

p1 = DimPlot(obj_joint,group.by="peaksMACS2_snn_res.0.8",label = TRUE, repel = TRUE, reduction = "wnn.umap") + NoLegend() + ggtitle('Joint WNN') + theme(aspect.ratio=1)
p2 = DimPlot(obj_joint,group.by="peaksMACS2_snn_res.0.8",reduction = 'rna.umap', label = TRUE,repel = TRUE, label.size = 2.5) +NoLegend() + ggtitle('RNA') + theme(aspect.ratio=1)
p3 = DimPlot(obj_joint,group.by="peaksMACS2_snn_res.0.8",reduction = 'atac.umap', label = TRUE,repel = TRUE, label.size = 2.5) +NoLegend() +ggtitle('ATAC') + theme(aspect.ratio=1)
p_peaksMACS2 = (p1 + p2 + p3) + plot_annotation(title="peaksMACS2 based cluster")

p1 = DimPlot(obj_joint,group.by="RNA_snn_res.0.5",label = TRUE, repel = TRUE, reduction = "wnn.umap") + NoLegend() + ggtitle('Joint WNN') + theme(aspect.ratio=1)
p2 = DimPlot(obj_joint,group.by="RNA_snn_res.0.5",reduction = 'rna.umap', label = TRUE,repel = TRUE, label.size = 2.5) +NoLegend() + ggtitle('RNA') + theme(aspect.ratio=1)
p3 = DimPlot(obj_joint,group.by="RNA_snn_res.0.5",reduction = 'atac.umap', label = TRUE,repel = TRUE, label.size = 2.5) +NoLegend() +ggtitle('ATAC') + theme(aspect.ratio=1)
p_RNA = (p1 + p2 + p3) + plot_annotation(title="RNA based cluster")

pdf(paste0(out_path,"UMAP_for_",merging_name,".pdf",sep=""),height=6,width=20)
print(p_wsnn);print(p_peaksMACS2);print(p_RNA)
dev.off()

saveRDS(obj_joint,paste0(out_path,"snRNA_ATAC_jointly_analyzed_for_",merging_name,".rds"))



