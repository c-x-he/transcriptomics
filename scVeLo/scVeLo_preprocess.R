# Prepare Seurat object for scVeLo analysis
# Script adapted by Caroline He from https://smorabit.github.io/tutorials/8_velocyto/ (it's basically the same script)
# 2023

# STEP 1: Set up enviornment
library(Seurat)
setwd("~/CdCS/scVeLo/")

# STEP 2: PREPARE SEURAT OBJECT
## 1. Load a pre-processed seurat object with celltypes labeled
CdCS90_filt <- readRDS("~/CdCS/cellranger_outs/CdCS/CdCS2/CdCS90_filt.rds")
seurat_obj<-CdCS90_filt 

## 2. Save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)

## 3. write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('counts.mtx'))

## 4. write dimesnionality reduction matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)

##  5. write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)

## Feed objects into scVeLo_analysis.sh