# Script adapted by Caroline He from the following sources:
# https://satijalab.org/seurat/articles/seurat_obj3k_tutorial.html
# https://swaruplab.bio.uci.edu/tutorial/snRNA-essentials/snRNA-essentials.html
# https://satijalab.org/seurat/articles/integration_introduction.html
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
# 2023

#STEP 1: SET UP ENVIRONMENT
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

# STEP 2: DEFINE VARIABLES
## I edited the code in a way where you can define your variables here and it will automatically populate the rest of the script
#_____________________________Change this every time_____________________________

Patient<- '299_C04' #The name of the patient the sample is from
Patient_label<- 'CADASIL1' #The given name of the patient the sample is from
diagnosis<- 'CADASIL' #The diagnosis associated with the sample
sex<- "Female"  #The sex associated with the sample
batch<- "1"  #The batch associated with the sample
doublet_percent<- 0.15 #the percent of doublets expected in this sample. For every 1,000 cells found from cellranger, add 1%. For example, 15,000 cells = 15% doublets, 30,000  cells = 30% doublets

#_____________________________Change this every time_____________________________

# STEP 2: PROCESS SAMPLES INDIVIDUALLY

## 1. Load dataset
## Dataset is the output from cellranger, 
## Select path to the directory where your matrices are, it is the filtered_feature_bc_matrix_SAMPLE_NAME, and it should contain 3 files - matrix, features and barcodes)
## Below is an example of how to input the directory, but change it for your particular folder

seurat_data <- Read10X(paste('~/multiome/Batch1_CADASIL/transcriptomics_cellranger_outs/run_countpremRNA_C',print(Patient), '/outs/filtered_feature_bc_matrix', sep = ""))

## 2. Create Seurat Object
seurat_obj <- CreateSeuratObject(counts = seurat_data, project = Patient)

## 3. Store percent of mitochondrial genes
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")

# 4. Add metadata
seurat_obj<-AddMetaData(seurat_obj, metadata= print(Patient), col.name = 'Patient')
seurat_obj<-AddMetaData(seurat_obj, metadata= print(diagnosis), col.name = 'Diagnosis')
seurat_obj<-AddMetaData(seurat_obj, metadata= print(batch), col.name = 'Batch')
seurat_obj<-AddMetaData(seurat_obj, metadata= print(sex), col.name = 'Sex')
seurat_obj<-AddMetaData(seurat_obj, metadata= print(Patient_label), col.name = 'Patient_label')

assign(paste('C', print(Patient), sep = ""), seurat_obj)

## 5. Subset low quality cells (at the individual cell level)
## Visualize the cell distribution

plot1 <- FeatureScatter(seurat_obj , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## write desired nFeature max cutoff given in plot2
#_____________________________Change this every time_____________________________

nFeature_max <- 8000

#_____________________________Change this every time_____________________________

## Subset low quality cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature_max  & percent.mt < 5)

# STEP 3. RUBN DOUBLETFINDER
## Doubletfinder identifies doublets and removes them

## 1. pK Identification (no ground-truth)
sweep.res.list_seurat_obj <- paramSweep_v3(seurat_obj, PCs = 1:10, sct = TRUE)
sweep.stats_seurat_obj <- summarizeSweep(sweep.res.list_seurat_obj, GT = FALSE)
bcmvn_seurat_obj <- find.pK(sweep.stats_seurat_obj)

## 2. Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)   #the number in here can be optimized for your dataset, we found that 0.25 works just as well for ours      
nExp_poi <- round(doublet_percent*nrow(seurat_obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 3. set pK value as value which results in max BCmetric
pK_BC_max<-droplevels(bcmvn_seurat_obj$pK[bcmvn_seurat_obj$BCmetric==max(bcmvn_seurat_obj$BCmetric)])
pK_value<-as.numeric(as.character(pK_BC_max))

## 4. Run DoubletFinder with varying classification stringencies 
seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

## 5. look at metatdata in seurat_obj and change the string of numbers after pANN
df<-seurat_obj@meta.data %>% dplyr::select(starts_with('pANN'))
pANN<-colnames(df)
new_pANN<-substring(pANN, 5)

seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = paste('pANN', new_pANN, sep = ""), sct = TRUE)
seurat_obj@meta.data[,"DF_hi.lo"] <- seurat_obj@meta.data[[paste('DF.classifications', new_pANN, sep = "")]]
seurat_obj@meta.data$DF_hi.lo[which(seurat_obj@meta.data$DF_hi.lo == "Doublet" & seurat_obj@meta.data[[paste('DF.classifications', new_pANN, sep = "")]] == "Singlet")] <- "Doublet_lo"
seurat_obj@meta.data$DF_hi.lo[which(seurat_obj@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"

## 6. remove doublets
Idents(seurat_obj) <- "DF_hi.lo"
seurat_obj <- subset(seurat_obj, idents = "Singlet")

Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj) 

assign(paste('S', print(Patient), sep = ""), seurat_obj)

## SAVEPOINT1
save.image("SAVEPOINT1_BATCH1-CADASIL_cellranger.RData")

#______________________________End of individual sample processing___________________

# STEP 4: MERGE OBJECTS

seurat_combined = merge(S109_C04, y=c(S118_C04, S147_C04, S160_C04, S201_C04, S240_C04, S291_C04, S299_C04, S321_C04, S327_C04, S342_C04, S364_C04))
seurat_old = merge(C109_C04, y=c(C118_C04, C147_C04, C160_C04, C201_C04, C240_C04, C291_C04, C299_C04, C321_C04, C327_C04, C342_C04, C364_C04))

## 1. visualize QC metrics from before pre-processing
## use table function to get the number of cells in each Sample as a dataframe
df <- as.data.frame(rev(table(seurat_old$orig.ident)))
colnames(df) <- c('orig.ident', 'n_cells')

## bar plot of the number of cells in each sample
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

## plot distributions of QC metrics, grouped by orig.ident
png('figures/Pre_QC_qc.png', width=10, height=10, res=200, units='in')
VlnPlot(
  seurat_old,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by='orig.ident',
  ncol = 1, pt.size=0)
dev.off()

# SAVEPOINT0
save.image("SAVEPOINT0_BATCH1-CADASIL_cellranger.RData")

# STEP 5: INTEGRATE OBJECT FOR BATCH CORRECTION
## this step takes a long time
seurat_combined.list <- SplitObject(seurat_combined, split.by = "orig.ident")

## 1. Normalize and integrate data
seurat_combined.list <- lapply(X = seurat_combined.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = seurat_combined.list, nfeatures = 3000)
seurat_combined.list <- PrepSCTIntegration(object.list = seurat_combined.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_combined.list, normalization.method = "SCT", 
                                  anchor.features = features)

seurat_combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(seurat_combined) <- "integrated"

# SAVEPOINT2
save.image("SAVEPOINT2_BATCH1-cellranger_integrated.RData")


# STEP 6: NORMALIZE, SCALE, AND VISUALIZE DATA
## 1. Normalize, scale, and visualize with UMAP
seurat_combined <- NormalizeData(seurat_combined)
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)
seurat_combined <- ScaleData(seurat_combined, verbose = FALSE)
seurat_combined <- RunPCA(seurat_combined, features=VariableFeatures(seurat_combined))

seurat_combined <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:10)
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:10)
DimPlot(seurat_combined, reduction = "umap",label = TRUE) + NoLegend()

Idents(seurat_combined) <- "Patient"
DimPlot(seurat_combined)

## 2. test different resolutions, use high resulution (3-7) to get many clusters
seurat_combined <- FindClusters(seurat_combined, resolution=5, 
                                verbose = FALSE,algorithm=1) 

Idents(seurat_combined) <- "seurat_clusters"
DimPlot (seurat_combined, label=T)

## 3. Change default assay to RNA, to set up for DGE and celltype ID
DefaultAssay(seurat_combined) <- "RNA"

#SAVEPOINT3, continue to celltype ID
save.image("~/multiome/CADASIL_transcriptomics_ONLY/SAVEPOINT3_BATCH1_cellranger_normalized.RData")








