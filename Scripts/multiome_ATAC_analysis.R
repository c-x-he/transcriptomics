
setwd("~/multiome")

#load object
load("~/multiome/batch_1_pre_processing_3_peaks_linked.RData")

#-------------------------------------------------------------------------------
#Differential peak analysis
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)

DefaultAssay(seurat_obj) <- 'peaks'

Idents(seurat_obj)<- "Diagnosis"

da_peaks <- FindMarkers(
  object = seurat_obj,
  ident.1 = "CADASIL",
  ident.2 = "Control",
  test.use = 'LR'
)

head(da_peaks)

#-------------------------------------------------------------------------------
#motif analysis 

library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
seurat_obj <- AddMotifs(
  object = seurat_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

#find overrepresented motifs
Idents(seurat_obj)<- "Diagnosis"

da_peaks <- FindMarkers(
  object = seurat_obj,
  ident.1 = 'CADASIL',
  ident.2 = 'Control',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

#compute motif activities
#run chromVar
seurat_obj <- RunChromVAR(
  object = seurat_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(seurat_obj) <- 'chromvar'

differential.activity <- FindMarkers(
  object = seurat_obj,
  ident.1 = 'CADASIL',
  ident.2 = 'Control',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = seurat_obj,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)

#Batch_1_Analysis_1.RData (1/17/23)

#------------------------------------------------------------------------------
#co-accessability networks
#https://stuartlab.org/signac/articles/cicero.html
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)

# convert to CellDataSet format and make the cicero object
seurat_obj.cds <- as.cell_data_set(seurat_obj, assay= "ATAC")
cicero_obj <- make_cicero_cds(seurat_obj.cds, reduced_coordinates = reducedDims(seurat_obj.cds)$UMAP)

# get the chromosome sizes from the Seurat object
DefaultAssay(seurat_obj) <- "peaks"
genome <- seqlengths(seurat_obj)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(cicero_obj, genomic_coords = genome.df, sample_num = 100)

#Find cis-co-accessable networks (CCANs)
ccans <- generate_ccans(conns)

#add links to seurat obj
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(seurat_obj) <- links

CoveragePlot(bone, region = "chr6-43488757-43490254")


#batch_1_analysis_1.RData

#-------------------------------------------------------------------------------
#Transcription Factor footprinting
#https://stuartlab.org/signac/articles/footprint.html

load("~/multiome/Batch_1_analysis_1.RData")

library(Signac)
library(Seurat)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)


# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
seurat_obj <- AddMotifs(seurat_obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

# gather the footprinting information for sets of motifs
seurat_obj <- Footprint(
  object = seurat_obj,
  motif.name = c("GATA2", "CEBPA", "EBF1"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

#Issue here
# plot the footprint data for each group of cells
p2 <- PlotFootprint(seurat_obj, features = c("GATA2", "CEBPA", "EBF1"))

p2 + patchwork::plot_layout(ncol = 1)


#------------------------------------------------------------------------------


#miscilllanius 
