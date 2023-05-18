# Use AddModule Score for cell type annotation
# I adapted the original script from Dasol, and I don't know what vignette he based it off

## STEP 1: LOAD THE MARKER LISTS FOR EACH CELLTYPE

## Astrocytes
ASTSignature <- read_csv("~/multiome/celltypes/cell_type_module_AST.csv")
ASTSignature <- ASTSignature$AST
ASTSignature <- list (ASTSignature)

## Oligodendrocytes
OLISignature <- read_csv("~/multiome/celltypes/cell_type_module_OLI.csv")
OLISignature <- OLISignature$OLI
OLISignature <- list (OLISignature)

## MICROGLIA
MICSignature <- read_csv("~/multiome/celltypes/cell_type_module_MIC.csv")
MICSignature <- MICSignature$MIC
MICSignature <- list (MICSignature)

## OPCs
OPCSignature <- read_csv("~/multiome/celltypes/cell_type_module_OPC.csv")
OPCSignature <- OPCSignature$OPC
OPCSignature <- list (OPCSignature)

## ENDOTHELIAL
ENDOTSignature <- read_csv("~/multiome/celltypes/cell_type_module_ENDOT.csv")
ENDOTSignature <- ENDOTSignature$ENDOT
ENDOTSignature <- list (ENDOTSignature)

## PERICYTES
PERSignature <- read_csv("~/multiome/celltypes/cell_type_module_PER.csv")
PERSignature <- PERSignature$PER
PERSignature <- list (PERSignature)

## NEURONAL
NEUSignature <- read_csv("~/multiome/celltypes/cell_type_module_NEU.csv")
NEUSignature <- NEUSignature$NEU
NEUSignature <- list (NEUSignature)

## EPEN
EPENSignature <- read.csv("~/multiome/celltypes/cell_type_module_EPEN.csv")
EPENSignature <- EPENSignature$EPEN
EPENSignature <- list (EPENSignature)

## Excitatory
Excsignature <- list(c("CBLN2", "EPHA6", "LAMA2", "CNTN5", "PDZD2", "CUX2", "RASGRF2","FAM19A1", "LINC01378", 
                       "CA10","COL5A2", "FAM19A1", "VAT1L", "COL24A1", "CBLN2", "NRP1", "NRG1", "HOMER1", 
                       "SLC35F3","FSTL4", "CNTN5", "RORB", "FSTL5", "IL1RAPL2", "CHN2","TOX", "CPNE4", 
                       "CADPS2", "POU6F2","TSHZ2", "HTR2C", "ITGA8", "ZNF385D", "ASIC2","CDH6", "CRYM", 
                       "NXPH2", "CPNE4","ADAMTSL1", "KIAA1217", "SORCS1", "HS3ST4", "TRPM3", "TOX", "SEMA3E",
                       "MEIS2", "SEMA5A","PTPRK", "PDZRN4", "CDH9", "THEMIS", "FSTL5", "CDH13", "CDH12", "CBLN2",
                       "MLIP","THEMIS", "RNF152", "NTNG2", "STK32B", "KCNMB2", "GAS2L3", "OLFML2B","POSTN", 
                       "B3GAT2", "NR4A2","HS3ST4", "KCNMB2", "MDFIC", "NTM", "CDH9", "MARCH1","TLE4", "FOXP2",
                       "KIAA1217"))

DefaultAssay(seurat_combined) <- "RNA"

# STEP 2: ADD THE MODULE SCORE ASSOCIATED WITH EACH LIST FOR EACH CELL INTO METADATA
seurat_combined <- AddModuleScore(object = seurat_combined,
                                  features = ASTSignature,
                                  name = 'AST')

seurat_combined <- AddModuleScore(object = seurat_combined,
                                  features = OLISignature,
                                  name = 'OLI')
seurat_combined <- AddModuleScore(object = seurat_combined,
                                  features = MICSignature,
                                  name = 'MIC')
seurat_combined <- AddModuleScore(object = seurat_combined,
                                  features = OPCSignature,
                                  name = 'OPC')
seurat_combined <- AddModuleScore(object = seurat_combined,
                                  features = PERSignature,
                                  name = 'PERI')
seurat_combined <- AddModuleScore(object = seurat_combined,
                                  features = ENDOTSignature,
                                  name = 'ENDOT')
seurat_combined <- AddModuleScore(object = seurat_combined,
                                  features = NEUSignature,
                                  name = 'NEUR')
seurat_combined <- AddModuleScore(object = seurat_combined,
                                  features = EPENSignature,
                                  name = 'EPEN')


Idents (seurat_combined) <- "integrated_snn_res.5"

## 2. Visualize the module scores for each cluster
## Identify the most likely celltype associated with each cluster, and label clusters that display a bimodal distribution or clusters with more than one marker as 'doublet'

VlnPlot(seurat_combined, features = c("AST1", "OLI1", "MIC1", "OPC1", "PERI1", 
                                      "ENDOT1", "NEUR1"), 
        pt.size = 0, stack = T, flip = T, sort=T)+NoLegend()

Idents (seurat_combined) <- "integrated_snn_res.5"

#----------------------------------------------------------------------------
## Input the celltype associated with each cluster in order here, the first slot labels cluster 0, second labels cluster 1, etc
new.cluster.ids<- c("Oli","Oli","Oli","Oli","Oli","Oli","Mic","Ast","Oli","Mic","Oli","Ast","Oli","Neu","OPC","doublet","Neu","Neu","Neu","Neu","Neu","Per","Ast","Oli","Neu","doublet","Neu","Neu","Neu","Neu","doublet","Neu","Oli","doublet","Oli","End","Neu","doublet","doublet","Mic","Neu","Ast","doublet","OPC","Neu","Neu","Neu","Neu","Neu","Neu","Ast","Oli","Neu","doublet","Neu","Neu","OPC","doublet","doublet","Per","Neu","Neu","Neu","doublet","doublet","Neu","doublet","doublet","Mic","Endo","Neu","Neu","Neu","doublet","doublet","Neu","Neu","doublet","doublet","doublet","Neu","Neu","doublet","doublet","End","doublet","doublet","doublet","doublet","Neu","doublet","End","doublet","doublet","doublet","doublet","Neu","doublet","doublet","doublet","doublet","doublet","Oli","doublet","doublet","doublet","doublet","doublet","doublet","doublet","Neu")

#----------------------------------------------------------------------------

## 3. Create a metadata column to store the celltype ID for each cell
names(new.cluster.ids) <- levels(seurat_combined)
seurat_combined <- RenameIdents(seurat_combined, new.cluster.ids)
seurat_combined[["Celltype"]] <- Idents(object = seurat_combined)

DimPlot(seurat_combined, label = F) 

## remove clusters that are likely doublets
Idents (seurat_combined) <- "Celltype"
seurat_combined <- subset(seurat_combined, idents = "doublet", invert = TRUE)

DimPlot(seurat_combined, label = F) 

##  4. remove clusters with low UMI and too high percent.mt
## plot distributions of QC metrics, grouped by cluster

VlnPlot(
  seurat_combined,
  features = c("nCount_RNA", "percent.mt"),
  group.by='integrated_snn_res.5',
  ncol = 1, pt.size=0)


#-------------------------------------------------------------------------------

remove_clusters<- c(5,8,9,22,23,32,39,56,68,69,84, 80,110,91,89)

#-------------------------------------------------------------------------------

## Remove low UMI and high percent.mt clusters
Idents (seurat_combined) <- "integrated_snn_res.5"
seurat_combined <- subset(seurat_combined, idents = print(remove_clusters), invert = TRUE)

Idents (seurat_combined) <- "Celltype"
DimPlot(seurat_combined, label = F) 

## 5. QC metric visualization after pre-processing

## plot distributions of QC metrics, grouped by orig.ident
png('figures/Post_QC_qc.png', width=10, height=10, res=200, units='in')
VlnPlot(
  seurat_combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by='orig.ident',
  ncol = 1, pt.size=0)
dev.off()


## use table function to get the number of cells in each Sample as a dataframe
df <- as.data.frame(rev(table(seurat_combined$orig.ident)))
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

png('figures/Post_QC_cells_per_sample.png', width=9, height=4, res=200, units='in')
print(p)
dev.off()

#SAVEPOINT4
save.image("~/multiome/CADASIL_transcriptomics_ONLY/SAVEPOINT4_BATCH1-cellranger_celltype.RData")

# STEP 6: RUN DIFFERENTIAL GENE EXPRESSION ANALYSIS ON EACH CELLTYPE

## 1. Specify the celltype
Celltype<- "Ast"

DefaultAssay(seurat_combined)<- "RNA"
cells <- subset(seurat_combined, idents = print(Celltype))
Idents (cells) <- "Diagnosis"

## 2. Run one test (here I use DESeq2)
comparison_1 <- FindMarkers(cells, ident.1 = "CADASIL", 
                            ident.2 = "Control", 
                            logfc.threshold = 0, min.pct = 0.1, random.seed = 50, test.use='DESeq2')
write.csv(comparison_1, file = paste(print(Celltype), "_CADASILvsControl_DESeq2.csv", , sep = ""))

## 2. Run another test (here I use wilcoxon)
comparison_2 <- FindMarkers(cells, ident.1 = "CADASIL", 
                            ident.2 = "Control", 
                            logfc.threshold = 0, min.pct = 0.1, random.seed = 50, test.use='wilcox')
write.csv(comparison_2, file = paste(print(Celltype), "_CADASILvsControl_Wilcoxon.csv", , sep = ""))
