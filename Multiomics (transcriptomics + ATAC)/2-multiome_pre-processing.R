# Pre-process multiome data, where data is extracted from the same cells
# Script adapted by Caroline He from https://stuartlab.org/signac/articles/pbmc_multiomic.html
# 2023

#STEP 1: SET UP ENVIORNMENT
library(Seurat)
#setRepositories(ind=1:3) 
#install.packages("Signac")
library(Signac)
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
#library(Herper)
#install_CondaTools(tools="macs2", env="PeakCalling_analysis", pathToMiniConda="~/miniconda3")

## Load annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

setwd("~/multiome")

#STEP 2: PROCESS MULTIOME PATIENTS INDIVIDUALLY

#_____________________________Change this every time_____________________________
#write (brainID)_C(region)
Patient<-'109_C04'
#_____________________________Change this every time_____________________________

## 1. Load the 10x hdf5 file (contains both data types)
inputdata.10x <- Read10X_h5(paste('/home/c_x_he/multiome/outs/',print(Patient), '/outs/filtered_feature_bc_matrix.h5', sep = ""))

## 2. extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

## 3. Create Seurat object
seurat_one <- CreateSeuratObject(counts = rna_counts, project = Patient)
seurat_one[["percent.mt"]] <- PercentageFeatureSet(seurat_one, pattern = "^MT-")


## 4. Add in the ATAC-seq data
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]


frag.file <- paste('/home/c_x_he/multiome/outs/',print(Patient), '/outs/atac_fragments.tsv.gz', sep = "")
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
seurat_one[["ATAC"]] <- chrom_assay

## 5. Perform ATAC Quality Control
DefaultAssay(seurat_one) <- "ATAC"
seurat_one <- NucleosomeSignal(seurat_one)
seurat_one <- TSSEnrichment(seurat_one)

## 6. Perform ATAC peak calling
## call peaks using MACS2
peaks <- CallPeaks(seurat_one, 
                   macs2.path = "/home/c_x_he/miniconda3/envs/PeakCalling_analysis/bin/macs2") #edit for your directory after installing MACS2

## 7. Remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

## 8. quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(seurat_one),
  features = peaks,
  cells = colnames(seurat_one)
)

## 9. create a new assay using the MACS2 peak set and add it to the Seurat object
seurat_one[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = frag.file,
  annotation = annotations
)
assign(paste('C', print(Patient), sep = ""), seurat_one)

## 10. Visualize
VlnPlot(seurat_one, features = c("nCount_ATAC", "nCount_RNA","percent.mt", "TSS.enrichment", "nucleosome_signal"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

assign(paste('S', print(Patient), sep = ""), seurat_one)


## 11. Define the filtration parameters for QC
#_____________________________Change this every time_____________________________
## Here are example values for each parameter.
nCount_ATAC_filter_1<- 7e3
nCount_ATAC_filter_2<-500
nCount_RNA_filter_1<- 5e3
nCount_RNA_filter_2<- 1e3
nucleosome_signal_filter <-3
TSS_enrichment_filter <-1
percent_mt_filter<- 5
#______________________________________________________________________________

seurat_one <- subset(
  x = seurat_one,
  subset = nCount_ATAC < nCount_ATAC_filter_1 &
    nCount_ATAC > nCount_ATAC_filter_2 &
    nCount_RNA < nCount_RNA_filter_1 &
    nCount_RNA > nCount_RNA_filter_2 &
    nucleosome_signal < nucleosome_signal_filter &
    TSS.enrichment > TSS_enrichment_filter &
    percent.mt < percent_mt_filter
)

## 12. Add metadata for each sample using the following command format
## seurat_obj<- AddMetaData(Seurat_obj, metadata= "metadata value to populate column", col.name = "name of metadata column")
## Example code:
S109_C04<-AddMetaData(S109_C04, metadata= "S109_C04", col.name = 'Patient')
S364_C04<-AddMetaData(S364_C04, metadata= "Control", col.name = 'Diagnosis')

# STEP 2: MERGE OBJECTS AND ANALYZE

seurat_obj = merge(S109_C04, y=c(S118_C04, S147_C04, S160_C04, S201_C04, S299_C04, S321_C04, S364_C04))
#_____________________________objects_are_merged_____________________________

##  1. RNA analysis
## Pipeline is the same as specified in the transcriptomics folder

## 2. ATAC analysis
## Exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q0')
seurat_obj <- RunSVD(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

seurat_obj <- FindMultiModalNeighbors(seurat_obj, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
seurat_obj <- RunUMAP(seurat_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seurat_obj <- FindClusters(seurat_obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
DimPlot(seurat_obj)


## 3. Perform sub-clustering on cluster to find additional structure
seurat_obj <- FindSubCluster(seurat_obj, cluster = 1, graph.name = "wsnn", algorithm = 3)
Idents(seurat_obj) <- "sub.cluster"

## 4. Celltype ID
## use the 3-transcripts_celltype_ID_and_DEG.R script for celltype ID 

p1 <- DimPlot(seurat_obj, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(seurat_obj, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

## 5. link genes to peaks
DefaultAssay(seurat_obj) <- "peaks"

## 6. Compute the GC content for each peak
library(BSgenome.Hsapiens.UCSC.hg38)
seurat_obj <- RegionStats(seurat_obj, genome = BSgenome.Hsapiens.UCSC.hg38)

## 7. Link peaks to genes
## will take a long time
seurat_obj <- LinkPeaks(
  object = seurat_obj,
  peak.assay = "peaks",
  expression.assay = "SCT"
)

# Multiome preprocessing complete, proceed to ATAC_analysis script and additional transcriptomics analysis
