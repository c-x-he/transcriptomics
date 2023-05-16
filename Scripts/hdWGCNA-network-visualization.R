#https://github.com/smorabit/hdWGCNA/blob/dev/vignettes/network_visualizations.Rmd

# single-cell analysis package
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
# network analysis & visualization package:
library(igraph)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)

setwd("~/TREM2/Oli-subset")
seurat_obj <- readRDS("~/TREM2/Oli-subset/hdWGCNA-seurat_obj-basic-end.rds")

## Individual module network plots: using top25 genes in each module
ModuleNetworkPlot(seurat_obj)


## Combined hub gene network plots
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)

g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)

plot(g)

## Applying UMAP to co-expression networks
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)
# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()



#Module UMAP Plot

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)

g <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)



