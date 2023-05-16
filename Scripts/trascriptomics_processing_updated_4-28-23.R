#trascriptomics_processing_updated_4-28-23
#https://satijalab.org/seurat/articles/seurat_obj3k_tutorial.html
#https://swaruplab.bio.uci.edu/tutorial/snRNA-essentials/snRNA-essentials.html
#https://satijalab.org/seurat/articles/integration_introduction.html
#https://www.biostars.org/p/9477922/
#https://github.com/chris-mcginnis-ucsf/DoubletFinder
#https://github.com/satijalab/seurat/issues/1717

library(dplyr)
library(Seurat)
library (Matrix)
library(ggplot2)
library(sctransform)
library (EnhancedVolcano)
library (DoubletFinder)
library (pheatmap)
library(circlize)
library(plyr)
library(tidyverse) 
library(ggpubr)
library(cluster)
library(umap)
library(ggrepel)
library(Hmisc)
library(DoubletFinder)

setwd("~/multiome/CADASIL_transcriptomics_ONLY")

#_____________________________Change this every time_____________________________
#write (brainID)_C(region)
#BATCH1_CADASIL: 109_C04 (0.26), 118_C04 (0.3), 147_C04 (0.27), 160_C04 (0.25), 201_C04 (M, 0.21), 299_C04 (0.15), 321_C04 (0.36)
#BATCH1_Control: 364_C04 (M, 0.11)
#BATCH2_Control: 291_C04 (M, 0.16), 327_C04 (M, 0.1), 342_C04 (0.14), 240_C04 (0.13)

Patient<- '299_C04'
diagnosis<- 'CADASIL'
sex<- "Female"
batch<- "1"
doublet_percent<- 0.15

#_____________________________Change this every time_____________________________

#####Lines 31 - 80 should be run individually for each sample

####Load Dataset---- (Select path to the directory where your matrices are, it is the filtered_feature_bc_matrix_SAMPLE_NAME, and it 
####should contain 3 files - matrix, features and barcodes)
seurat_data <- Read10X(paste('~/multiome/Batch1_CADASIL/transcriptomics_cellranger_outs/run_countpremRNA_C',print(Patient), '/outs/filtered_feature_bc_matrix', sep = ""))
#seurat_data <- Read10X(paste('/home/c_x_he/multiome/Batch1_CADASIL/transcriptomics_cellranger_outs/Batch2_Control_outs/run_countpremRNA_C',print(Patient), '/outs/filtered_feature_bc_matrix', sep = ""))

####Create Seurat Object
seurat_obj <- CreateSeuratObject(counts = seurat_data, project = Patient)

#### Store mt.percent
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")

#add metadata
seurat_obj<-AddMetaData(seurat_obj, metadata= print(Patient), col.name = 'Patient')
seurat_obj<-AddMetaData(seurat_obj, metadata= print(diagnosis), col.name = 'Diagnosis')
seurat_obj<-AddMetaData(seurat_obj, metadata= print(batch), col.name = 'Batch')
seurat_obj<-AddMetaData(seurat_obj, metadata= print(sex), col.name = 'Sex')

assign(paste('C', print(Patient), sep = ""), seurat_obj)

###SUBSET low quality cells
plot1 <- FeatureScatter(seurat_obj , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#_____________________________Change this every time_____________________________
#write desired nFeature max cutoff given plot2
nFeature_max <- 8000
# 7500, 6000, 7500, 10000, 7500, 8000, 6000, 
# 8000 
# 9000, 7000, 6000, 5500
#_____________________________Change this every time_____________________________

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature_max  & percent.mt < 5)



#DoubletFinder (per sample)
####Doublet Finder
## pK Identification (no ground-truth)
sweep.res.list_seurat_obj <- paramSweep_v3(seurat_obj, PCs = 1:10, sct = TRUE)
sweep.stats_seurat_obj <- summarizeSweep(sweep.res.list_seurat_obj, GT = FALSE)
bcmvn_seurat_obj <- find.pK(sweep.stats_seurat_obj)

## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)   #argument to run celltype annotation before doublet removal       
nExp_poi <- round(doublet_percent*nrow(seurat_obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#set pK value as value which results in max BCmetric
pK_BC_max<-droplevels(bcmvn_seurat_obj$pK[bcmvn_seurat_obj$BCmetric==max(bcmvn_seurat_obj$BCmetric)])
pK_value<-as.numeric(as.character(pK_BC_max))

## Run DoubletFinder with varying classification stringencies 
seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

#look at metatdata in seurat_obj and change the string of numbers after pANN
df<-seurat_obj@meta.data %>% dplyr::select(starts_with('pANN'))
pANN<-colnames(df)
new_pANN<-substring(pANN, 5)

####check your metadata for the pANN value and replace the next 3 lines of the script with the value for the sample
seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = paste('pANN', new_pANN, sep = ""), sct = TRUE)
seurat_obj@meta.data[,"DF_hi.lo"] <- seurat_obj@meta.data[[paste('DF.classifications', new_pANN, sep = "")]]
seurat_obj@meta.data$DF_hi.lo[which(seurat_obj@meta.data$DF_hi.lo == "Doublet" & seurat_obj@meta.data[[paste('DF.classifications', new_pANN, sep = "")]] == "Singlet")] <- "Doublet_lo"
seurat_obj@meta.data$DF_hi.lo[which(seurat_obj@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"

#remove doublets
Idents(seurat_obj) <- "DF_hi.lo"
seurat_obj <- subset(seurat_obj, idents = "Singlet")

Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj) 

assign(paste('S', print(Patient), sep = ""), seurat_obj)


save.image("SAVEPOINT1_BATCH1-CADASIL_cellranger.RData")
## SAVEPOINT1

#______________________________End of individual sample processing___________________

#### Merge objects
seurat_combined = merge(S109_C04, y=c(S118_C04, S147_C04, S160_C04, S201_C04, S240_C04, S291_C04, S299_C04, S321_C04, S327_C04, S342_C04, S364_C04))
seurat_old = merge(C109_C04, y=c(C118_C04, C147_C04, C160_C04, C201_C04, C240_C04, C291_C04, C299_C04, C321_C04, C327_C04, C342_C04, C364_C04))

####visualize QC metrics from before pre-processing
# use table function to get the number of cells in each Sample as a dataframe
df <- as.data.frame(rev(table(seurat_old$orig.ident)))
colnames(df) <- c('orig.ident', 'n_cells')

# bar plot of the number of cells in each sample
p <- ggplot(df, aes(y=n_cells, x=reorder(orig.ident, -n_cells), fill=orig.ident)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  NoLegend() + RotatedAxis() +
  ylab(expression(italic(N)[cells])) + xlab('orig.ident') +
  ggtitle(paste('Total cells:', sum(df$n_cells))) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="lightgray", size=0.5),
  )

png('figures/Pre_QC_cells_per_sample.png', width=9, height=4, res=200, units='in')
print(p)
dev.off()

# plot distributions of QC metrics, grouped by orig.ident
png('figures/Pre_QC_qc.png', width=10, height=10, res=200, units='in')
VlnPlot(
  seurat_old,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by='orig.ident',
  ncol = 1, pt.size=0)
dev.off()

save.image("SAVEPOINT0_BATCH1-CADASIL_cellranger.RData")


#Integrate to Perform Batch Correction----
#Integrate Data (Batch effect correction)
seurat_combined.list <- SplitObject(seurat_combined, split.by = "orig.ident")

#save.image("SAVEPOINT1_BATCH1-CADASIL_cellranger.RData")

#---

#Normalize Data
seurat_combined.list <- lapply(X = seurat_combined.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = seurat_combined.list, nfeatures = 3000)
seurat_combined.list <- PrepSCTIntegration(object.list = seurat_combined.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_combined.list, normalization.method = "SCT", 
                                  anchor.features = features)

save.image("SAVEPOINT2_BATCH1-cellranger_integrated.RData")

#load("~/multiome/CADASIL_transcriptomics_ONLY/SAVEPOINT2_BATCH1-cellranger_integrated.RData")
seurat_combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(seurat_combined) <- "integrated"

#SAVEPOINT2
save.image("SAVEPOINT2_BATCH1-cellranger_integrated.RData")

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
#load("~/multiome/CADASIL_transcriptomics_ONLY/SAVEPOINT2_BATCH1-cellranger_integrated.RData")



## normalize and scale
seurat_combined <- NormalizeData(seurat_combined)
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_combined)
seurat_combined <- ScaleData(seurat_combined, verbose = FALSE)

seurat_combined <- RunPCA(seurat_combined, features=VariableFeatures(seurat_combined))

seurat_combined <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:10)
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:10)
DimPlot(seurat_combined, reduction = "umap",label = TRUE) + NoLegend()

Idents(seurat_combined) <- "Patient"
DimPlot(seurat_combined)

###test different resolutions, use high resulution (3-7) to get many clusters
seurat_combined <- FindClusters(seurat_combined, resolution=5, 
                                verbose = FALSE,algorithm=1) 

Idents(seurat_combined) <- "seurat_clusters"
DimPlot (seurat_combined, label=T)

####Change default assay to RNA
DefaultAssay(seurat_combined) <- "RNA"


#SAVEPOINT3, continue to celltype ID
save.image("~/multiome/CADASIL_transcriptomics_ONLY/SAVEPOINT3_BATCH1_cellranger_normalized.RData")

#feed seurat_combined into transcripts_celltype_ID script







