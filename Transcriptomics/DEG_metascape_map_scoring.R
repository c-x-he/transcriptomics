
#for Peri up
VEGFA_list<- list (c("CALR","CLIC1","CSRP2","EIF4G2","GAPDH","DNAJA1","HSPA1A","HSPB1","HSP90AA1","IGFBP7","CXCL8","LDHA","RPL10A","NFKBIA","PFN1","PGK1","PTMA","RPL5","RPL7","RPL26","RPL27","RPLP2","RPS6","RPS11","CCL2","SDCBP","SELE","SOD2","SSR4","TMSB4X","TXN","TMSB10","ADAMTS1","RACK1","STIP1","RPL13A","ADAMTS9","APOLD1","SH3BGRL3"))

#for Endo up 
#VEGFA_list<-list (c("ACTG1","CALR","CFL1","CLIC1","FN1","GAPDH","DNAJA1","HSPA1A","HSPB1","HSP90AA1","CXCL8","LDHA","RPL10A","NFKBIA","PFN1","PGK1","PTGS2","PTMA","RAC1","RPL5","RPL7","RPL26","RPL27","RPLP2","RPS6","RPS11","CCL2","SDCBP","SELE","SOD2","SSR4","TMSB4X","TXN","TMSB10","ADAMTS1","RACK1","RPL13A","ADAMTS9","APOLD1","SH3BGRL3"))

Celltype<- "Peri"
cells <- subset(seurat_combined, idents = print(Celltype))
                          
DefaultAssay(cells) <- "RNA"

cells <- AddModuleScore(object = cells,
                                  features = VEGFA_list,
                                  name = 'VEGFA Score')

Idents (cells) <- "Diagnosis"
VlnPlot(cells, features = 'VEGFA.Score1', group.by = "Diagnosis")

VlnPlot(cells, features = 'VEGFA.Score1',group.by = "Patient_label", pt.size=0)+NoLegend()

#add new column with pateints labeled as ctrl or cadasil #
seurat_combined$Patient_label <- NA # or any other initialization value
seurat_combined$Patient_label[seurat_combined$Patient == "201_C04"] <- 'CADASIL 1'
seurat_combined$Patient_label[seurat_combined$Patient == "109_C04"] <- 'CADASIL 2'
seurat_combined$Patient_label[seurat_combined$Patient == "299_C04"] <- 'CADASIL 3'
seurat_combined$Patient_label[seurat_combined$Patient == "321_C04"] <- 'CADASIL 4'
seurat_combined$Patient_label[seurat_combined$Patient == "160_C04"] <- 'CADASIL 5'
seurat_combined$Patient_label[seurat_combined$Patient == "118_C04"] <- 'CADASIL 6'
seurat_combined$Patient_label[seurat_combined$Patient == "147_C04"] <- 'CADASIL 7'
seurat_combined$Patient_label[seurat_combined$Patient == "364_C04"] <- 'Control 1'
seurat_combined$Patient_label[seurat_combined$Patient == "342_C04"] <- 'Control 2'
seurat_combined$Patient_label[seurat_combined$Patient == "240_C04"] <- 'Control 3'
seurat_combined$Patient_label[seurat_combined$Patient == "327_C04"] <- 'Control 4'
seurat_combined$Patient_label[seurat_combined$Patient == "291_C04"] <- 'Control 5'

#SAVEPOINT5, new metadata col
#save.image("~/multiome/CADASIL_transcriptomics_ONLY/SAVEPOINT5_BATCG1-cellranger_new_metadata_col.RData")

#trying camila's code

###Heatmap Proteostasis genes
DefaultAssay(seurat_combined)<- "RNA"
Celltype<- "Ast"
Ast.cells <- subset(seurat_combined, idents = print(Celltype))

Autophagy_list <- c("EEF1A1","GFAP","HSP90AA1","LAMP2","UBB","UBC","VIM","BMPR1A",
                    "ITGB1","ITGB8","SMAD1","SMAD9","WWTR1","TUBA1A","RAB7A","TAB2",
                    "ERBIN","GYG2","ZFP36L1","HSPA1A","HSPB1","APOE","FTL","GJA1","DENND4A",
                    "ANKRD28","RHOQ","MAN1C1","COLEC12","HGF","SDC4","GPC6","PTGES3","ASPH",
                    "STOM","ANO6","ETV6","FGF2","CD44","H3-3B","PLCE1","PTPRK","DNER","CAMK2D",
                    "MRAS","PTTG1IP","ANGPTL4","PEBP1","CCNH","AHCYL1","IGFBP7","IFITM3","DPYSL2",
                    "UTRN","IRS2","PLEC","NPAS2","ST8SIA1","ARHGEF3")


AUTOPHAGYSUBSET <- subset(Ast.cells, features = Autophagy_list)

DimPlot (AUTOPHAGYSUBSET)
Idents (AUTOPHAGYSUBSET) <- "Patient_label"
av.exp <- AverageExpression(AUTOPHAGYSUBSET)$RNA

pheatmap(av.exp)

cal_z_score <- function(x){(x - mean(x)) / sd(x)}

av.exp <- t(apply(av.exp, 1, cal_z_score))

pheatmap (av.exp, cluster_cols = F, fontsize_row = 6)

###Proteostasis gene list (upregulated in E280A)
autophagy <- list(c("EEF1A1","GFAP","HSP90AA1","LAMP2","UBB","UBC","VIM","BMPR1A",
                    "ITGB1","ITGB8","SMAD1","SMAD9","WWTR1","TUBA1A","RAB7A","TAB2",
                    "ERBIN","GYG2","ZFP36L1","HSPA1A","HSPB1","APOE","FTL","GJA1","DENND4A",
                    "ANKRD28","RHOQ","MAN1C1","COLEC12","HGF","SDC4","GPC6","PTGES3","ASPH",
                    "STOM","ANO6","ETV6","FGF2","CD44","PLCE1","PTPRK","DNER","CAMK2D",
                    "MRAS","PTTG1IP","ANGPTL4","PEBP1","CCNH","AHCYL1","IGFBP7","IFITM3","DPYSL2",
                    "UTRN","IRS2","PLEC","NPAS2","ST8SIA1","ARHGEF3"))

Ast.cells <- AddModuleScore(object = Ast.cells,
                            features = autophagy,
                            ctrl = 5,
                            name = 'Autophagy')

VlnPlot (Ast.cells, features = "Autophagy1", group.by = "Patient_label", pt.size=0)+NoLegend()












