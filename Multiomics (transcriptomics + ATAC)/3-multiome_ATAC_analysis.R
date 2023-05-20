# ATAC Analysis
# Script adapted by Caroline He from https://stuartlab.org/signac/articles/pbmc_multiomic.html
# 2023

# STEP 1: SET UP ENVIORNMENT
setwd("~/multiome")

## load seurat object from 2-multiome_pre-processing.R
load("~/multiome/batch_1_pre_processing_3_peaks_linked.RData")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)


# STEP 2: PERFORM DIFFERENTIAL PEAK ANALYSIS
Idents(seurat_obj)<- "Diagnosis"

da_peaks <- FindMarkers(
  object = seurat_obj,
  ident.1 = "CADASIL", #define groups to compare here, peaks upregulated in ident.1 will laive positive log2fc values
  ident.2 = "Control",
  test.use = 'LR'
)

head(da_peaks)


# STEP 3: PERFORM MOTIF ANALYSIS
## 1. Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

## 2. add motif information
seurat_obj <- AddMotifs(
  object = seurat_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

## 3. find overrepresented motifs
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

## 4. get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

## 5. compute motif activities
## run chromVar
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
