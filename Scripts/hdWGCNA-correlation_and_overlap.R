#Marker gene overlap analysis
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

seurat_all <- readRDS("~/new_frontal_official/hdWGCNA/Ast-subset/hdWGCNA-seurat_obj-basic-end.rds")
seurat_obj<- subset(x = seurat_all, subset = Diagnosis == c('Sporadic', 'Control'))

# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- seurat_obj$Diagnosis
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold= 0.25
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_all,
  deg_df = markers,
  fc_cutoff = 0.25 # log fold change cutoff for overlap analysis
)

setwd("~/new_frontal_official/hdWGCNA/Ast-subset/DEG_overlap")
write.csv(overlap_df, "Ast-E280A.vs.Control_UP.csv")
