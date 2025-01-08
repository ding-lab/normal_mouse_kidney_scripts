## add motif annotation to the joint object

library(optparse)

library(BSgenome.Mmusculus.UCSC.mm10)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)
library(motifmatchr)
library(patchwork)
library(pheatmap)
library(viridis)
library(dplyr)
library(ggplot2)

library(Matrix)
library(stringr)

option_list <- list( 
    make_option(c("-s","--merging_name"),type="character",default=NULL,help="merging_name",metavar="character"),
    make_option(c("-i", "--input"), type="character",default=NULL,help="path to the rds object",metavar="character"),
    make_option(c("-a","--peak_assay"),type="character",default=NULL,help="peak assay, either peaks or peaksinters",metavar="character"),
    make_option(c("-o","--output_path"),type="character",default=NULL,help="path to the output  directory",metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

opt$merging_name <- "NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed"
opt$input <- "results/snRNA_ATAC_jointly_analyzed_for_NMK_E16.5X2_P0X2_W3X3_W12X3_W52X3_W92X3_20240719_doublet_removed.rds"
opt$peak_assay <- "peaksMACS2"
opt$output_path <- "results/"

merging_name <- opt$merging_name

## read in atac object (could be joint object)
input <- opt$input
if (!file.exists(input)){stop(paste0("file ",input," does not exist!"))}
ATAC <- readRDS(input)

## specify assay
peak_assay <- opt$peak_assay

out_path <- opt$output_path
dir.create(out_path)

###PFMatrixList loading...
cat("getting JASPAR2020 matrix...")
library(JASPAR2020)
library(TFBSTools)
jaspar <- function (collection = "CORE", ...)
{
	opts <- list()
	opts["tax_group"] <- "vertebrates"
	opts["collection"] <- collection
	opts <- c(opts, list(...))
	out <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
	if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
		names(out) <- paste(names(out), TFBSTools::name(out),sep = "_")
	return(out)
}
pfm <- jaspar()

## add motifs
DefaultAssay(ATAC) <- peak_assay
motif.matrix <- CreateMotifMatrix(
#  features = granges(ATAC),
	features = StringToGRanges(rownames(ATAC),sep=c("-","-")),
	pwm = pfm,
	genome = BSgenome.Mmusculus.UCSC.mm10, 
	sep = c(":", "-"),use.counts = FALSE)

motif <- CreateMotifObject(
	data = motif.matrix,
	pwm = pfm)

# Add the Motif object to the assay
# ATAC[[peak_assay]] <- AddMotifObject(object = ATAC[[peak_assay]],motif.object = motif)
ATAC <- AddMotifs(ATAC,genome = BSgenome.Mmusculus.UCSC.mm10,pfm = pfm)

# Add the Motif object to the assay
ATAC <- RegionStats(object = ATAC,genome = BSgenome.Mmusculus.UCSC.mm10,sep = c("-", "-"))
ATAC <- RunChromVAR(object = ATAC,genome = BSgenome.Mmusculus.UCSC.mm10)
DefaultAssay(ATAC) <- 'chromvar'

#saveRDS(ATAC,paste('Motif_analysis/chromVAR/',sample,'_chromVAR.rds',sep=""))
saveRDS(ATAC,paste0(out_path,merging_name,'_chromVAR.rds',sep=""))

