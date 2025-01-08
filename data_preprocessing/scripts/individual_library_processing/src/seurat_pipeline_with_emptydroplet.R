## use Seurat to process RNA and ATAC data, for each individual object

library(optparse)
set.seed(1234)
library(future)
plan("multicore", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)

option_list = list(
 make_option(c("-s", "--sample"),
  	type="character", 
  	default=NULL, 
  	help = "sample_name", 
  	metavar="character"),
 make_option(c("-d","--data"),
  	type="character",
  	default=NULL,
  	help = "path to Cellranger-arc data folder (e.g. cellranger output's raw matrices folder)",
  	metavar="character"),
 make_option(c("-m","--macs2_path"),
  	type="character",
  	default=NULL,
  	help = "path to installed MACS2",
  	metavar="character"),
 make_option(c("-o","--output_folder"),
     type="character",
     default=NULL,
     help = "output folder where a sample subfolder will be created",
     metavar="character"),
 make_option(c("-c","--chrom_size"), #/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/Resources/genome_reference/mm10.chrom.sizes.txt
     type="character",
     default=NULL,
     help = "path to mm10.chrom.sizes.txt file",
     metavar="character"),
#CellRanger ATAC QC metrics
 make_option(c("--prf_min"),
	type="integer", 
	default=3000, 
	help = "peak_region_fragments_minimum value for filtering", 
	metavar="integer"),
 make_option(c("--prf_max"),
	type="integer", 
	default=20000, 
	help = "peak_region_fragments_maximum value for filtering", 
	metavar="integer"),
 make_option(c("--pct_min"),
	type="integer", 
	default=15, 
	help = "pct_reads_in_peaks_minimum value for filtering", 
	metavar="integer"),
  make_option(c("--bl_ratio"),
 	type="double", 
 	default=0.05, 
 	help = "blacklist_ratio_minimum value for filtering", 
 	metavar="double"),
#Changed to default=4, based on the latest Signac-vignette
 make_option(c("--ns_max"),
	type="integer", 
	default=4, 
	help = "nucleosome_signal_maximum value for filtering", 
	metavar="integer"),
 make_option(c("--tss"),
	type="integer", 
	default=4, 
	help = "tss_enrichment_minimum value for filtering", 
	metavar="integer"),
 make_option(c("--pc_num"),
	type="integer", 
	default=30, 
	help = "number of principal components to use", 
	metavar="integer"),
 make_option(c("--pc_first"),
	type="integer", 
	default=1, 
	help = "first principal components to use (should be 1 or 2)", 
	metavar="integer"),
make_option(c("--pre_filter_ATAC"),
            type="integer",
            default=500,
            help="min number of ATAC fragments per cell",
            metavar="integer"),

#### RNA QC metrics
make_option(c("--pre_filter"),
            type="integer",
            default=300,
            help="min number of reads per cell to prefilter",
            metavar="integer"),
make_option(c("--fdr_emptydrops"),
            type="double",
            default=0.01,
            help="fdr cutoff for emptDrops model",
            metavar="double"),
make_option(c("--nfeature_min"),
            type="integer",
            default=200,
            help="nFeature_RNA min value for filtering",
            metavar="integer"),
make_option(c("--nfeature_max"),
            type="integer",
            default=10000,
            help="nFeature_RNA max value for filtering",
            metavar="integer"),
make_option(c("--ncount_min"),
            type="integer",
            default=1000,
            help="nCount_RNA min value for filtering",
            metavar="integer"),
make_option(c("--ncount_max"),
            type="integer",
            default=80000,
            help="nCount_RNA max value for filtering",
            metavar="integer"),
make_option(c("--mito_max"),
            type="double",
            default=10,
            help="maximum allowed mitochondrial fraction from 0 to 100",
            metavar="double"),
make_option(c("--scrublet"),
            type="character",
            default=NULL,
            help="path to scrublet folder output",
            metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# an example below:
# opt$sample="20210218-NMKE16_5M"
# opt$data=paste0("/diskmnt/Projects/Mouse_kidney_development/combo_snRNA_snATAC/cellranger_arc_2.0/",opt$sample)
# opt$macs2_path="/diskmnt/Projects/Users/rliu/Software/miniconda/envs/Seurat_RNA_ATAC_ST/bin/macs2"
# opt$output_folder="/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snRNA_ATAC_combo/Seurat/with_SoupX_ambientRNAremoval/out/"
# opt$chrom_size="/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/Resources/genome_reference/mm10.chrom.sizes.txt"
# opt$prf_min=1000
# opt$pct_min=15
# opt$ns_max=5
# opt$pc_first=2
# opt$pc_num=50
# opt$tss=3
# opt$scrublet=paste0("/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snRNA_ATAC_combo/doublet_detection/doublet_auto_detection/out/input_1/combined/",opt$sample,"/",opt$sample,"_combo_scrublet_output_table.csv")


if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (sample_name,atac_data).n", call.=FALSE)
}

#####################################
####### FUNCTIONS ##################
####################################

iterative_removal <- function(all_peaks.f) {
  
  #just load existing peaks if any
  recentered_p=StringToGRanges(regions = all_peaks.f$new_peak, sep = c("-", "-"))
  
  cat(paste0('finding overlapping peaks\n'))
  overlapping=as.data.table(x = findOverlaps(query = recentered_p, 
                                             subject = recentered_p)) # find which peaks overlap
  overlapping=overlapping[queryHits!=subjectHits,]
  overlapping.peak.number <- unique(x = overlapping$queryHits) #these are numbers of overlapping peaks that denote their position in all_peaks.f table
  recentered_non_overlapping=all_peaks.f[-overlapping.peak.number,] # select peaks that are not overlapping as non-overlapping peaks
  # fwrite(recentered_non_overlapping,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_nonOverlapping.',add_filename,'.tsv'),
  #        sep='\t',row.names=FALSE)
  if (length(overlapping.peak.number)>0) {
    tmp <- data.table(chr = all_peaks.f$seqnames[overlapping.peak.number], 
                      num = overlapping.peak.number)
    overlapping.peak.number.split <- split(tmp, by = 'chr', keep.by = T) #split peaks by chromosome 
    registerDoParallel(cores=25)
    #this is where iterative removal of peaks is done
    best_in_overlapping_num <- foreach(peak.numbers=overlapping.peak.number.split) %dopar% {
      cat('removing overlapping peaks in each chromosome\n')
      iterative_removal_core (peak.numbers = peak.numbers, overlapping.f = overlapping)
    }
    stopImplicitCluster()
    best_in_overlapping_num <- do.call('c', best_in_overlapping_num) #combine best peak numbers from all chromosomes
    best_in_overlapping_cancer <- all_peaks.f[best_in_overlapping_num,] #extract peaks themselves
    # fwrite(best_in_overlapping_cancer,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_Overlapping.',add_filename,'.tsv'),
    #        sep='\t',row.names=FALSE)
    recentered_final.f=rbindlist(list(recentered_non_overlapping,best_in_overlapping_cancer))
  } else {
    recentered_final.f=recentered_non_overlapping
  }
  final.overlaps <-  recentered_final.f$new_peak %>% 
    unique %>% 
    StringToGRanges %>% 
    countOverlaps
  if (sum(final.overlaps>1)>0) {
    stop("Execution stopped. Overlapping peaks remained")
  }
  # fwrite(recentered_final.f,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.',add_filename,'.tsv'),sep='\t',
  #        row.names=FALSE)
  
  return(recentered_final.f)
}

# this works like a charm
iterative_removal_core <- function(peak.numbers, overlapping.f) {
  chr = peak.numbers$chr[1]
  running.vector <- peak.numbers$num
  peaks.to.trash <- NULL
  peaks.to.keep <- NULL
  while (length(running.vector) != 0) {
    n <- running.vector[1] # this is the first and the best peak since peaks are sorted by scores
    neighbor.peaks.num.discard <- overlapping.f[queryHits==n, subjectHits] #find positions of other peaks overlapping with the first one 
    running.vector <- setdiff(running.vector, neighbor.peaks.num.discard) # remove them from the list of peaks
    running.vector <- setdiff(running.vector, n)
    peaks.to.keep <- c(peaks.to.keep, n) # add this peak to the keeping list
    peaks.to.trash <- unique(c(peaks.to.trash, neighbor.peaks.num.discard)) # add neighbors to the list of peaks to discard
  }
  cat('done\n')
  return(peaks.to.keep)
}

filter_N_peaks <- function(peak.dt) {

  gr <- StringToGRanges(peak.dt$new_peak, sep = c("-", "-")) #get GRanges object from peaks
  seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,gr) #extract fasta sequence
  names(seq) <- peak.dt$new_peak
  peaks.match.pattern <- vmatchPattern("N", seq) #match peak sequence with N in them
  peaks.withN <- names(peaks.match.pattern)[elementNROWS(peaks.match.pattern)>0] # these are peaks that contain N in their sequence
  toreturn <- peak.dt[! new_peak %in% peaks.withN,]
  # fwrite(toreturn,paste0('peaks/',length(samples.id),'_',cancer.type, '_recentered_final.reproducible.filtered.',add_filename,'.tsv'),sep='\t',
  #        row.names=FALSE)
  return(toreturn)
}

get_upper_outlier_cutoff <- function(vect) {
  suspected.outliers <- vect[vect >= quantile(vect, 0.97, na.rm=TRUE)]
  test <- rosnerTest(vect,
                     k = length(suspected.outliers)
  )
  cutoff <- test$all.stats %>% filter(Outlier) %>% filter(Value > mean(vect)) %>% pull(Value) %>% tail(1) 
  return(cutoff)
}

get_lower_outlier_cutoff <- function(vect) {
  suspected.outliers <- vect[vect >= quantile(vect, 0.97, na.rm=TRUE)]
  test <- rosnerTest(vect,
                     k = length(suspected.outliers)
  )
  cutoff <- test$all.stats %>% filter(Outlier) %>% filter(Value < mean(vect)) %>% pull(Value) %>% head(1) 
  return(cutoff)
}
####################################################

print("Input parameters")
print(paste("--peak_region_fragments_min:",opt$prf_min,sep=""))
print(paste("--peak_region_fragments_max:",opt$prf_max,sep=""))
print(paste("--pct_reads_in_peaks_minimum:",opt$pct_min,sep=""))
print(paste("--blacklist_ratio_minimum:",opt$bl_ratio,sep=""))
print(paste("--nucleosome_signal_maximum:",opt$ns_max,sep=""))
print(paste("--tss_enrichment_minimum:",opt$tss,sep=""))
print(paste("--pc_first:",opt$pc_first,sep=""))
print(paste("--pc_num:",opt$pc_num,sep=""))

##input data
sample=opt$sample
data_folder=opt$d
path.to.chrom.size <- opt$chrom_size
outputpath <- opt$output_folder

print(paste("Cellranger-arc data:",data_folder,sep=""))
##output data
print(sample)

#####LOAD REQUIRED PACKAGES##########
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(EnvStats)
library(outliers)
suppressMessages(library(DropletUtils))
suppressMessages(library(doParallel))

filter <- dplyr::filter
select <- dplyr::select
###########################
########LOAD IN DATA#######
###########################


outputpath=paste(outputpath, '/', sample,"/",sep="")
dir.create(outputpath, showWarnings = F)
setwd(outputpath)

cat('Reading in input matrices\n')
counts <- Read10X_h5(paste(data_folder,"/outs/raw_feature_bc_matrix.h5",sep=""))
rna_counts <- counts$`Gene Expression`
atac_counts <- counts$Peaks
## in this metadata atac_fragments columns is eqvivalent to passed_filters column from singlecell.csv file of cellranger-atac output
#more info on this file https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/per_barcode_metrics#header
metadata <-read.csv(file=paste(data_folder,"/outs/per_barcode_metrics.csv",sep=""), header = TRUE, row.names = 1)
fragment.path <- paste(data_folder,"/outs/atac_fragments.tsv.gz",sep="")


# filter only fragments in the standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# this is practically does not remove anything if you using filtered matrices 
barcodes.non.zero.atac <- Matrix::colSums(atac_counts) > 0
barcodes.non.zero.rna <- Matrix::colSums(rna_counts) > 0
barcodes.non.zero.both <- barcodes.non.zero.atac & barcodes.non.zero.rna

atac_counts.filtered <- atac_counts[,barcodes.non.zero.both]
rna_counts.filtered <- rna_counts[,barcodes.non.zero.both]

# do empty drops detection
out <- emptyDrops(rna_counts.filtered)
fdr.cut <- opt$fdr_emptydrops
out <- out %>% 
  data.frame(check.rows = F, check.names = F) %>% 
  mutate(is.cell = FDR <= fdr.cut)
is.cell <- out$FDR <= fdr.cut #Might need to mess with this cutoff

# Check if p-values are lower-bounded by 'niters'
# (increase 'niters' if any Limited==TRUE and Sig==FALSE)
table(Sig=is.cell, Limited=out$Limited)
write("Number of cells/nuclei:",stdout())
write(sum(is.cell, na.rm=TRUE), stdout())

trueCells = rownames(subset(out,FDR <= fdr.cut))
emptyCells = rownames(subset(out,FDR > fdr.cut))

write("Number of trueCells:",stdout())
write(length(trueCells), stdout())
write("Number of emptyCells:",stdout())
write(length(emptyCells), stdout())


options(repr.plot.width=6, repr.plot.height=5)
ggplot(data.frame(out), aes(Total, -LogProb, color = is.cell)) +
  geom_point(size = 0.5) +
  cowplot::theme_cowplot() +
  labs(x="Total UMI count", y="-Log Probability")
ggsave(paste0(outputpath,"/Scatterplot_emptyDrop_",sample,".pdf"), width = 6, height = 4.5)
dev.off()

br.out <- barcodeRanks(rna_counts.filtered)

# Making a knee plot.

ggplot(data=data.frame(br.out), aes(x=rank, y=total)) +
  geom_point() +
  scale_y_log10() + 
  scale_x_log10() +
  xlab('Rank') + ylab('Total') +
  geom_hline(aes(yintercept = metadata(br.out)$knee, color="knee")) +
  geom_hline(aes(yintercept = metadata(br.out)$inflection, color="inflection")) +
  scale_color_manual(values = c('knee'="dodgerblue", 'inflection'="forestgreen")) +
  theme_classic()
ggsave(paste0(outputpath,"/Kneeplot_emptyDrop_",sample,".pdf"), width = 6, height = 4.5)
dev.off()

#keep only true cells
atac_counts.filtered2 <- atac_counts.filtered[,intersect(colnames(atac_counts.filtered),trueCells)]
rna_counts.filtered2 <- rna_counts.filtered[,intersect(colnames(rna_counts.filtered),trueCells)]

write("Dims of rna_counts.filtered2:",stdout())
write(dim(rna_counts.filtered2), stdout())


#remove cells with cellranger excluded reasons 1 and 3 
#1 = more than 1 gel bead and a nuclei in a partition; 
#2 = low fraction of fragments overlapping peaks; 
#3 = barcode excluded because it was a multiplet
cells.to.keep <- rownames(metadata %>% 
                            filter(excluded_reason != 1) %>%
                            filter(excluded_reason != 3))
length(cells.to.keep)
atac_counts.filtered2 <- atac_counts.filtered2[,intersect(cells.to.keep, colnames(atac_counts.filtered2))]
rna_counts.filtered2 <- rna_counts.filtered2[,intersect(cells.to.keep, colnames(rna_counts.filtered2))]


#cat('prefilter\n')
#keep cells with at least 50 genes detected
#rna.prefiltered.cells <- Matrix::colSums(rna_counts.filtered>0) >= 50
#atac_counts.filtered <- atac_counts.filtered[,rna.prefiltered.cells]
#rna_counts.filtered <- rna_counts.filtered[,rna.prefiltered.cells]


# create object based on RNA data
panc <- CreateSeuratObject(
  counts = rna_counts.filtered2, 
  project = sample
)
panc <- AddMetaData(panc,metadata[-1:-2])
panc[["percent.mt"]] <- PercentageFeatureSet(panc, pattern = "^mt-")
panc$log10GenesPerUMI <- log10(panc$nFeature_RNA) / log10(panc$nCount_RNA)
panc@meta.data <- panc@meta.data %>% mutate(nUMI = nCount_RNA,
                                            nGene = nFeature_RNA)

#Add gene annotations mm10
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

#Add ATAC assay
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts.filtered2,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = fragment.path,
   min.cells = -1, 
   annotation = annotations
)
dim(chrom_assay)

panc[["ATAC"]] <- chrom_assay

## add scrublet if applicable
if (!is.null(opt$scrublet)){
  cat('add scrublet\n')
  scrub.path <- opt$scrublet
  scrub.table <- fread(scrub.path, header = T) %>% data.frame(row.names = 1)
  panc <- AddMetaData(panc, scrub.table)
  panc$predicted_doublet <- panc$predicted_doublet_rna & panc$predicted_doublet_atac
  panc$predicted_doublet[is.na(panc$predicted_doublet_rna)] <- NA
}


if(file.exists(paste0('recentered_final.filtered',sample,'.tsv'))) {
  recentered_final <- fread(paste0('recentered_final.filtered',sample,'.tsv'), header = TRUE)
} else {
  ####2021-03-20: Change Peaks to MACS2
  ###########################################################
  ############MACS2 peak calling#############################
  ###########################################################
  DefaultAssay(panc) <- 'ATAC'
  peaks <- CallPeaks(
    object = panc,
    macs2.path=opt$macs2_path
  )
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
  
  all_peaks=as.data.table(peaks)
  fwrite(all_peaks,paste0(outputpath,'MACS2_peaks.',sample,'.tsv'),sep='\t')
  
  # recenter peaks
  all_peaks[,peak_center:=start+relative_summit_position]
  all_peaks[,recentered_start:=peak_center-250]
  all_peaks[,recentered_end:=peak_center+250]
  all_peaks[,length:=recentered_end-recentered_start+1]
  all_peaks[,new_peak:=paste0(seqnames,"-", recentered_start, '-',recentered_end)]
  
  ####Now check that new start and end don't go beyond the chromosome boundaries
  chr_size=read.table(path.to.chrom.size,sep='\t',header=FALSE)
  colnames(chr_size)=c('seqnames','chr_length')
  all_peaks=merge(all_peaks,chr_size,all.x=TRUE)
  all_peaks=all_peaks[recentered_end<=chr_length && recentered_start>=0,]
  
  ### do iterative removal of overlapping peaks
  all_peaks <- all_peaks[order(neg_log10qvalue_summit, decreasing = T), ]
  recentered_final <- iterative_removal(all_peaks)
  recentered_final <- filter_N_peaks(recentered_final)
  fwrite(recentered_final,paste0(outputpath,'recentered_final.filtered',sample,'.tsv'),sep='\t')
  
}


recentered_p=StringToGRanges(recentered_final$new_peak, sep = c(":", "-"))
matrix.counts <- FeatureMatrix(
  fragments = Fragments(panc@assays$ATAC),
  features = recentered_p,
  sep = c("-","-"),
  cells = colnames(panc)
)

panc[['ATAC_MACS2']] <- CreateChromatinAssay(counts = matrix.counts,
                                             annotation = annotations,
                                             genome = 'mm10',
                                             fragments = fragment.path, min.features = -1)

DefaultAssay(panc)<-'ATAC_MACS2'

# remove ATAC assay
panc[['ATAC']] <- NULL

#add some more QC stuff
peak.data <- GetAssayData(object = panc, assay = 'ATAC_MACS2', slot = "counts")
total_fragments_cell <- panc$atac_fragments
peak.counts <- colSums(x = peak.data)
frip <- peak.counts *100 / total_fragments_cell
panc <- AddMetaData(object = panc, metadata = frip, col.name = 'pct_read_in_peaks_500MACS2')
panc <- AddMetaData(object = panc, metadata = peak.counts, col.name = 'peak_region_fragments_500MACS2')

print(colnames(panc@meta.data))

# plot pre-filter metadata
#### RNA QC
#pdf(paste(outputpath,"/QC_in_sample_",sample, "_RNA.pdf", sep=""), width=15, height=9)
#VlnPlot(object = panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()

## ATAC
# pdf(paste(outputpath,"/QC_in_sample_",sample, "_ATAC.pdf", sep=""), width=10, height=9)
# VlnPlot(object = panc, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)
# dev.off()

# plot metadata associations
pdf(paste0(outputpath,"/FeatureScatter_in_",sample,".pdf",sep=""),width=17,height=7)
plot1 <- FeatureScatter(object = panc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = panc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

ggplot(panc@meta.data, aes_string("log10GenesPerUMI")) +
    geom_density(fill ='#cdb4db') +
    cowplot::theme_cowplot() +
    scale_x_log10() +
    geom_vline(xintercept = 0.8)
ggsave(paste0(outputpath,"/Density_log10GenesPerUMI_",sample,".pdf"), width = 6, height = 4.5)
dev.off()

nUMI.upper.cutoff <- get_upper_outlier_cutoff(panc$nUMI)
nGene.upper.cutoff <- get_upper_outlier_cutoff(panc$nGene)

ggplot(panc@meta.data, aes_string("nUMI", 'nGene', color = 'percent.mt')) +
  geom_point(size = 0.5) +
  cowplot::theme_cowplot() +
  scale_x_log10() + scale_y_log10() +
  geom_vline(xintercept = nUMI.upper.cutoff) +
  geom_vline(xintercept = opt$ncount_min, color = "firebrick") +
  geom_hline(yintercept = opt$nfeature_min, color = "firebrick")
  geom_hline(yintercept = nGene.upper.cutoff)
ggsave(paste0(outputpath,"/FeatureScatter_outlier_cutoffs_in_",sample,".pdf"), width = 6, height = 4.5)
dev.off()


###########################################################
############Quality Control OF SC-ATAC MACS2 DATA################
###########################################################
#https://satijalab.org/signac/articles/panc_vignette.html

# compute nucleosome signal score per cell
DefaultAssay(object = panc) <- 'ATAC_MACS2'
panc <- NucleosomeSignal(object = panc)

# compute TSS enrichment score per cell
panc <- TSSEnrichment(object = panc, fast = FALSE)

#add blacklist ratio and fraction of reads in peaks - this metric is not reported in cellranger-arc output
blacklist_counts <- CountsInRegion(panc, assay = 'ATAC_MACS2', region = blacklist_mm10) %>% data.frame(); colnames(blacklist_counts) <- 'blacklist_region_fragments'
panc <- AddMetaData(panc, blacklist_counts)
panc$blacklist_ratio <- panc$blacklist_region_fragments / panc$peak_region_fragments_500MACS2

# inspecting TSS-enrichment scores
panc$high.tss <- ifelse(panc$TSS.enrichment > opt$tss, 'High', 'Low')
tss_plot=TSSPlot(panc, group.by = 'high.tss') + NoLegend()

# inspecting fragment length periodicity
panc$nucleosome_group <- ifelse(panc$nucleosome_signal > opt$ns_max, paste('NS >', opt$ns_max),  paste('NS <', opt$ns_max))
fragment_period_plot=FragmentHistogram(object = panc, group.by = 'nucleosome_group')	

QCplot <- VlnPlot(object = panc, 
                  features = c('pct_read_in_peaks_500MACS2', 'peak_region_fragments_500MACS2','TSS.enrichment','nucleosome_signal'), 
                  ncol =4)

pdf(paste(outputpath,"/",sample,"_0_ATAC_QC.pdf",sep=""),height=6,width=12)
print(tss_plot)
try(print(fragment_period_plot))
print(QCplot)
dev.off()



tss.peak <- ggplot(panc@meta.data, aes_string("atac_TSS_fragments", 'peak_region_fragments_500MACS2', color ='TSS.enrichment')) +
  geom_point(size = 0.5) +
  cowplot::theme_cowplot() +
  scale_x_log10() + scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "firebrick")

frag.vs.pct <- ggplot(panc@meta.data, aes_string("atac_fragments", 'pct_read_in_peaks_500MACS2', color ='TSS.enrichment')) +
  geom_point(size = 0.5) +
  cowplot::theme_cowplot() +
  scale_x_log10() +
  geom_hline(yintercept = opt$pct_min, color = "firebrick") +
  geom_hline(yintercept = get_lower_outlier_cutoff(panc$pct_read_in_peaks_500MACS2)) +
  geom_vline(xintercept = 1000, color = "firebrick")+
  geom_vline(xintercept = get_upper_outlier_cutoff(panc$atac_fragments))

peak.frag_vs.pct <- ggplot(panc@meta.data, aes_string("peak_region_fragments_500MACS2", 'pct_read_in_peaks_500MACS2', color ='TSS.enrichment')) +
  geom_point(size = 0.5) +
  cowplot::theme_cowplot() +
  scale_x_log10() +
  geom_hline(yintercept = opt$pct_min, color = "firebrick") +
  geom_hline(yintercept = get_lower_outlier_cutoff(panc$pct_read_in_peaks_500MACS2)) +
  geom_vline(xintercept = opt$prf_min, color = "firebrick")+
  geom_vline(xintercept = get_upper_outlier_cutoff(panc$peak_region_fragments_500MACS2))


pdf(paste(outputpath,"/",sample,"_ATAC_QC_scatter.pdf",sep=""),height=6,width=6)
print(tss.peak)
print(frag.vs.pct)
print(peak.frag_vs.pct)
dev.off()


#####
cat('filtering\n')
panc$remove_based_on_RNA <- case_when(panc$nUMI >= nUMI.upper.cutoff & panc$nGene >= nGene.upper.cutoff ~ 'Remove',
                                      panc$nGene < opt$nfeature_min & panc$nUMI < opt$ncount_min ~ 'Remove',
                                      panc$percent.mt >= opt$mito_max ~ 'Remove',
                                      TRUE ~ 'Keep')
panc$remove_based_on_ATAC <- case_when(panc$peak_region_fragments_500MACS2 >= get_upper_outlier_cutoff(panc$peak_region_fragments_500MACS2) ~ 'Remove',
                                       panc$atac_fragments >= get_upper_outlier_cutoff(panc$atac_fragments) ~ 'Remove',
                                       panc$peak_region_fragments_500MACS2 < opt$prf_min ~ 'Remove',
                                       panc$pct_read_in_peaks_500MACS2 < opt$pct_min ~ 'Remove',
                                       panc$TSS.enrichment < opt$tss~ 'Remove',
                                       panc$blacklist_ratio >= opt$bl_ratio ~ 'Remove',
                                       panc$nucleosome_signal >= opt$ns_max  ~ 'Remove',
                                      TRUE ~ 'Keep')

panc$remove <- panc$remove_based_on_RNA == 'Remove' | panc$remove_based_on_ATAC == 'Remove'
#remove cells that are outliers for these QC metrics
panc <- subset(
      x = panc,
      subset = remove, invert = TRUE
)

#### after QC
pdf(paste(outputpath,"/After_QC_in_sample_",sample, "_ATAC_MACS2_RNA.pdf", sep=""), width=25, height=9)
VlnPlot(object = panc, features = c('pct_read_in_peaks_500MACS2', 'peak_region_fragments_500MACS2',  'TSS.enrichment', 'nucleosome_signal', "nFeature_RNA", "nCount_RNA"), ncol = 6)
dev.off()

##################################################
##Normalization and linear dimensional reduction##
##################################################
# RNA analysis
DefaultAssay(panc) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
panc <- NormalizeData(panc, assay = 'RNA')
panc <- CellCycleScoring(panc, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

panc <- SCTransform(panc, 
                    vars.to.regress = c("nCount_RNA","percent.mt","S.Score", "G2M.Score"),
                    return.only.var.genes = F, verbose = FALSE) %>% 
  RunPCA(npcs = opt$pc_num, verbose = FALSE) %>% 
  RunUMAP(dims = 1:opt$pc_num, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(panc) <- "ATAC_MACS2"
panc <- RunTFIDF(panc)
panc <- FindTopFeatures(panc, min.cutoff = 'q0')
panc <- RunSVD(panc)
panc <- RunUMAP(panc, reduction = 'lsi', dims = opt$pc_first:opt$pc_num, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#Check if first LSI-component correlated with the sequencibg depth. If it is, then re-run using LSI components starting from 2 (for exaample, 2:30 instead of 1:30)
depth_corr_plot=DepthCor(panc)
pdf(paste(outputpath,"/",sample,"_DepthCorrelation_1_QC.pdf",sep=""),height=6,width=12)
print(depth_corr_plot)
dev.off()

##################################################
##Non-linear dimension reduction and clustering###
##################################################

# perform graph-based clustering and non-linear dimension reduction for visualization
panc <- FindMultiModalNeighbors(panc, 
                                reduction.list = list("pca", "lsi"), 
                                dims.list = list(1:opt$pc_num, opt$pc_first:opt$pc_num))
panc <- RunUMAP(panc, nn.name = "weighted.nn", 
                reduction.name = "wnn.umap", 
                reduction.key = "wnnUMAP_")
panc <- FindClusters(panc, graph.name = "wsnn", algorithm = 3, verbose = T)


p1 <- DimPlot(panc, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(panc, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(panc, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")


pdf(paste(outputpath,"/",sample,"_2_Dimplots.pdf",sep=""),height=6,width=18)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

p1 <- FeaturePlot(panc, reduction = "wnn.umap", features = c("nCount_RNA","nFeature_RNA")) 
p2 <- FeaturePlot(panc, reduction = "umap.atac", features = c("TSS.enrichment","pct_read_in_peaks_500MACS2")) 
pdf(paste0(sample,"_3_Featureplots.pdf"),height=6,width=14)
print(p1)
print(p2)
dev.off()

# #Linking peaks to genes
# DefaultAssay(object = panc) <- 'ATAC_MACS2'
# # first compute the GC content for each peak
# panc <- RegionStats(panc, genome = BSgenome.Mmusculus.UCSC.mm10)
# 
# # link peaks to genes
# panc <- LinkPeaks(
#   object = panc,
#   peak.assay = 'ATAC_MACS2',
#   expression.assay = "SCT"
# )

#Save object
saveRDS(panc,file = paste(outputpath,"/",sample, "_processed_multiomic.rds", sep=""))

# plot scrublet
if (!is.null(opt$scrublet)){
  pdf(paste0("DimPlot_predicted_doublet_",sample,".pdf"),useDingbats=FALSE)
  DimPlot(object = panc, group.by = 'predicted_doublet',cols = c('black', 'yellow'), 
          reduction = "wnn.umap",label=TRUE,label.size=6)
  dev.off()
}

