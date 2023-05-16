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

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

setwd("~/multiome")
#_____________________________Change this every time_____________________________
#write (brainID)_C(region)
#109_C04, 118_C04, 147_C04, 160_C04, 201_C04, 299_C04, 321_C04, 364_C04
Patient<-'109_C04'
#_____________________________Change this every time_____________________________

# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5(paste('/home/c_x_he/multiome/outs/',print(Patient), '/outs/filtered_feature_bc_matrix.h5', sep = ""))

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
seurat_one <- CreateSeuratObject(counts = rna_counts, project = Patient)
seurat_one[["percent.mt"]] <- PercentageFeatureSet(seurat_one, pattern = "^MT-")


# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- "hg38"

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

#ATAC Quality Control
DefaultAssay(seurat_one) <- "ATAC"
seurat_one <- NucleosomeSignal(seurat_one)
seurat_one <- TSSEnrichment(seurat_one)

#peak calling
# call peaks using MACS2
peaks <- CallPeaks(seurat_one, 
                   macs2.path = "/home/c_x_he/miniconda3/envs/PeakCalling_analysis/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(seurat_one),
  features = peaks,
  cells = colnames(seurat_one)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
seurat_one[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = frag.file,
  annotation = annotations
)
assign(paste('C', print(Patient), sep = ""), seurat_one)

#visualize
VlnPlot(seurat_one, features = c("nCount_ATAC", "nCount_RNA","percent.mt", "TSS.enrichment", "nucleosome_signal"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

assign(paste('S', print(Patient), sep = ""), seurat_one)

#load("~/multiome/batch_1_pre_processing.RData")

#_____________________________Change this every time_____________________________
nCount_ATAC_filter_1<- 7e3
#1e3 5e3 7e3 7e3 8e3 7e3 7e3 1e4
#7e3 most common
nCount_ATAC_filter_2<-500
#300 300 3e3 1e3 900 2e3 1e3 6e3
#make 500 min
nCount_RNA_filter_1<- 5e3
#5e3 5e3 5e3 5e3 5e3 5e3 5e3 5e3
#5e3 most common
nCount_RNA_filter_2<- 1e3
#1e3 1.1e3 1e3 1e3 1e3 1e3 1e3 1e2
#1e3 most common
nucleosome_signal_filter <-3
#3
TSS_enrichment_filter <-1
percent_mt_filter<- 5

#______________________________________________________________________________
load("~/multiome/batch_1_pre_processing_2.RData")

#S109_C04, S118_C04, S147_C04, S160_C04, S201_C04, S299_C04, S321_C04, S364_C04
S109_C04<-AddMetaData(S109_C04, metadata= "S109_C04", col.name = 'Patient')
S109_C04<-AddMetaData(S109_C04, metadata= "CADASIL", col.name = 'Diagnosis')
S118_C04<-AddMetaData(S118_C04, metadata= "S118_C04", col.name = 'Patient')
S118_C04<-AddMetaData(S118_C04, metadata= "CADASIL", col.name = 'Diagnosis')
S147_C04<-AddMetaData(S147_C04, metadata= "S147_C04", col.name = 'Patient')
S147_C04<-AddMetaData(S147_C04, metadata= "CADASIL", col.name = 'Diagnosis')
S160_C04<-AddMetaData(S160_C04, metadata= "S160_C04", col.name = 'Patient')
S160_C04<-AddMetaData(S160_C04, metadata= "CADASIL", col.name = 'Diagnosis')
S201_C04<-AddMetaData(S201_C04, metadata= "S201_C04", col.name = 'Patient')
S201_C04<-AddMetaData(S201_C04, metadata= "CADASIL", col.name = 'Diagnosis')
S299_C04<-AddMetaData(S299_C04, metadata= "S299_C04", col.name = 'Patient')
S299_C04<-AddMetaData(S299_C04, metadata= "CADASIL", col.name = 'Diagnosis')
S321_C04<-AddMetaData(S321_C04, metadata= "S321_C04", col.name = 'Patient')
S321_C04<-AddMetaData(S321_C04, metadata= "CADASIL", col.name = 'Diagnosis')
S364_C04<-AddMetaData(S364_C04, metadata= "S364_C04", col.name = 'Patient')
S364_C04<-AddMetaData(S364_C04, metadata= "Control", col.name = 'Diagnosis')

#S109_C04, S118_C04, S147_C04, S160_C04, S201_C04, S299_C04, S321_C04, S364_C04
seurat_obj = merge(S109_C04, y=c(S118_C04, S147_C04, S160_C04, S201_C04, S299_C04, S321_C04, S364_C04))

save.image("~/multiome/update_batch_1_2.RData")
#_____________________________objects_are_merged_____________________________
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- subset(
  x = seurat_obj,
  subset = nCount_ATAC < nCount_ATAC_filter_1 &
    nCount_ATAC > nCount_ATAC_filter_2 &
    nCount_RNA < nCount_RNA_filter_1 &
    nCount_RNA > nCount_RNA_filter_2 &
    nucleosome_signal < nucleosome_signal_filter &
    TSS.enrichment > TSS_enrichment_filter &
    percent.mt < percent_mt_filter
)

# RNA analysis
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q0')
seurat_obj <- RunSVD(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

seurat_obj <- FindMultiModalNeighbors(seurat_obj, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
seurat_obj <- RunUMAP(seurat_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seurat_obj <- FindClusters(seurat_obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
DimPlot(seurat_obj)

#batch_1_preprocessing_2

# perform sub-clustering on cluster to find additional structure
#seurat_obj <- FindSubCluster(seurat_obj, cluster = 1, graph.name = "wsnn", algorithm = 3)
#Idents(seurat_obj) <- "sub.cluster"

#celltype ID
#make list of celltype markers
Ast<- c("AQP4","GJA1","GJB6","SLC4A4","SLC1A2","F3","BMPR1B","FGFR3",
        "SLC39A12","CLDN10","DIO2","ALDOC","ALDH1L1","SLC1A3","CLU",
        "ATP13A4","SLCO1C1","SLC14A1","CHRDL1","GPR37L1","ACSBG1",
        "ATP1A2","SLC25A18","EDNRB","PPAP2B","GFAP","SOX9","SDC4",
        "PPP1R3C","NCAN","MLC1","GLI3","SLC7A11","ACSL6","RFX4","ID4",
        "AGT","SFXN5","GABRG1","PAX6","RORB","GRM3","PTPRZ1","PSD2",
        "SLC6A11","ATP1B2","NTSR2","S1PR1","SLC15A2","ELOVL2","TRIL",
        "SCARA3","MGST1","KIAA1161","FAM107A","BCAN","SPARCL1","NWD1",
        "NTRK2","SLC7A10","SCG3","ACOT11","KCNN3","MFGE8","RANBP3L",
        "GPC5","EZR","ADHFE1","GABRB1","TMEM47","PAMR1","CPE","FABP7",
        "LIX1","SLC13A5","IL33","SLC7A2","EGFR","PREX2","NDRG2","DTNA",
        "ABCD2","HEPACAM","RGS20","ARHGEF26","GPAM","CHI3L1","ADCYAP1R1",
        "GDPD2","SLC1A4","POU3F2","ETNPPL","MEGF10","MT3","TTYH1",
        "PRODH","PLCD4","DDAH1","LGR4","HTRA1")

Endo<- c("APOLD1","FLT1","RGS5","PTPRB","TM4SF1","ABCB1","ITM2A","SDPR",
         "SLCO1A2","FN1","EMCN","ESAM","NOSTRIN","CD34","SLC38A5","CYYR1",
         "PODXL","CDH5","VWF","MECOM","CD93","ABCG2","TEK","PALMD","ERG",
         "CLDN5","PECAM1","KDR","ITGA1","ICAM2","ATP10A","ANXA3","CA4",
         "MYCT1","GIMAP6","ANXA1","PTRF","KIAA1462","EBF1","HMCN1","ENG",
         "IGFBP7","ARHGAP29","ANXA2","OCLN","HIGD1B","SLC2A1","GNG11",
         "SLC19A3","EPAS1","TBX3","SRGN","SOX7","SLC16A4","CAV1","CLIC5",
         "VIM","HEG1","CCDC141","C10ORF10","EDN1","ROBO4","TMEM204",
         "PROM1","IFITM1","LEF1","COBLL1","WWTR1","HBB","ETS1","SLC39A8",
         "COL4A1","OSMR","ADCY4","TIE1","EDN3","THBD","BSG","AHNAK",
         "MYO1B","IL1R1","CXCL12","CLEC14A","GATA2","SGPP2","SHE","PLTP",
         "SPARC","ACVRL1","MMRN2","NID1","TNFSF10","FOXC1","UACA","CGNL1",
         "MFSD2A","NET1","ABCC9","FLI1","C1ORF54")

Mic<- c("CCL4","CCL3","CSF1R","CX3CR1","P2RY12","C1QB","RGS1","GPR183",
        "GPR34","CTSS","LAPTM5","CD53","IL1A","C3AR1","PLEK","FCGR2A",
        "CD83","ITGAM","P2RY13","CD86","TREM2","TYROBP","FCER1G","NCKAP1L",
        "SELPLG","SLC2A5","CD14","C1QC","C1QA","MPEG1","HAVCR2","PTAFR",
        "LY86","AIF1","ALOX5AP","LPCAT2","SLA","PTPRC","FCGR1A","CCL2",
        "BLNK","IL10RA","BCL2A1","C5AR1","RHOH","CD84","CSF3R","TLR7",
        "TLR2","HPGDS","LCP1","CD300A","FYB","MRC1","FAM105A","IRF8",
        "LCP2","RGS10","CD74","PTPN6","TBXAS1","LYZ","DOCK2","TMEM119",
        "NLRP3","ARHGDIB","CCRL2","IKZF1","ARHGAP25","DOCK8","HEXB",
        "THEMIS2","SAMSN1","HK2","PLD4","APBB1IP","ITGB2","RUNX1",
        "SLCO2B1","TLR1","FGD2","HCLS1","GPR84","OLFML3","MAFB","PIK3CG",
        "SIGLEC7","IL1B","PIK3R5","IL6R","CXCL16","CLEC4A","PTGS1","SUSD3",
        "LYN","VAV1","SLC11A1","RBM47","SYK","C10ORF128")

Neu<- c("RELN","VIP","GAD2","SYNPR","GAD1","CNR1","SYT1","SCG2","TAC3",
        "GABRG2","GABRA1","STMN2","DLX1","KCNC2","TMEM130","RAB3C","SST",
        "VSTM2A","SNAP25","ROBO2","CALB2","KIT","CNTNAP2","GABRB2","FSTL5",
        "NRXN3","SYT4","GRIA1","VSNL1","INA","NPY","GRIN2A","IGF1","PENK",
        "ELAVL2","MYT1L","KCNQ5","MEG3","NRIP3","CHGB","CLSTN2","SCN2A",
        "RAB3B","ZMAT4","NELL1","PNOC","ERBB4","SPHKAP","C11ORF87",
        "ADARB2","SLC4A10","KIAA1324","GRIN2B","BCL11A","CELF4","PNMA2",
        "DISP2","NYAP2","SV2B","SERPINI1","SLC2A13","RGS8","RTN1","NAP1L2",
        "CCK","C8ORF34","DYNC1I1","SRRM4","RBFOX1","SLC12A5","NDRG4",
        "ZNF804A","LPPR4","SLITRK4","GPR158","NDNF","KCNJ3","PCSK2",
        "CADPS","OLFM3","GABBR2","SULT4A1","GLRA2","SYT13","CACNA2D1",
        "GDA","SYNGR3","MAL2","PGM2L1","SLC7A14","GPR83","FRMPD4","NELL2",
        "RGS4","CSRNP3","DCX","ATP1A3","ST8SIA3","UCHL1","GAP43")


Oli<- c("PLP1","MOBP","CLDN11","MBP","UGT8","ERMN","MOG","MAG","OPALIN","CNP",
        "MAL","GPR37","TF","MYRF","GJB1","ASPA","ENPP2","BCAS1","LPAR1","FA2H",
        "ENPP6","APOD","CNTN2","CRYAB","KLK6","ERBB3","ANLN","SEPT4","PLEKHB1",
        "TMEFF2","ST18","PTGDS","PEX5L","SLAIN1","QDPR","PLLP","TMEM125","HHIP",
        "LGI3","TUBB4A","PLEKHH1","S1PR5","MAP6D1","GSN","EVI2A","EDIL3",
        "CMTM5","GJC3","CA14","NFASC","TPPP","TMEM88B","TRIM59","CDH19","APLP1",
        "NIPAL4","ADAMTS4","STMN4","S100B","CA2","PRR18","OLIG1","FOLH1",
        "NINJ2","NDRG1","SLC24A2","SGK2","GALNT6","KCNA1","SH3TC2","TTLL7",
        "SH3GL3","DOCK5","SCD","FEZ1","SLC44A1","RHOU","PPP1R16B","TSPAN2",
        "C10ORF90","TNFAIP6","NKAIN2","MOB3B","PRKCQ","PPP1R14A","PLA2G16",
        "DBNDD2","CDK18","PCDH9","ANO4","AGPAT4","OMG","FGFR2","TMEM63A",
        "GLTP","CCP110","PLEKHG3","RAB33A","PSAT1","ZNF536")

OPCs<-c("PDGFRA","TNR","PCDH15","SHC4","VCAN","LHFPL3","NEU4","GPR17",
        "PTPRZ1","OLIG1","MMP16","DSCAM","C8ORF46","SEMA5A","MATN4",
        "UGT8","GRIA3","CNTN1","BCAS1","SULF2","LUZP2","GJC3","NXPH1",
        "APOD","MEGF11","LRRTM3","BRINP3","GALNT13","GRIA4","MYT1","SUSD5",
        "LRRN1","SOX10","PRKCQ","SOX6","ITGB8","TMEM255A","GFRA1","RLBP1",
        "PNLIP","XYLT1","GPSM2","TMEM255B","SEZ6L","STK32A","C14ORF37",
        "LPPR5","SEMA3D","CSPG4","CSMD3","TMEM132B","SCRG1","KCNH8",
        "CACNG4","UGDH","DPP6","BCAT1","PLLP","ERBB3","RNF43","S100B",
        "SORCS1","OLIG2","CHRNA4","KCNJ16","PPAPDC1A","CSMD1","OPCML",
        "PRKG2","COBL","FIGN","ACAN","TGFA","NLGN1","SLC6A13","EMID1",
        "CHST6","TMEM100","GAL3ST1","EDIL3","KCNJ10","SLITRK3","SNTG1",
        "CSPG5","ERBB4","SLC35F1","B3GAT2","C1QL1","SERINC5","CKAP2",
        "LRRTM4","DPYD","SLITRK1","NCALD","CALCRL","SPP1","ZNF488",
        "ADAM12","SULF1","HAS2")

Peri<-c("KDR","APLNR","MFAP4","MCAM","DDR2","COL5A2","MYO1B","TNS1","LMOD1",
        "COL3A1","IGFBP5","DES","FMOD","PRRX1","DPT","RGS5","HSD11B1","VIM",
        "FAP","FBN1","ANGPTL2","SLC38A11","FRZB","SERPING1","NR1H3","PAMR1",
        "FIBIN","ADAM33","REM1","COX4I2","GNB4","PDE5A","SYNPO2","FABP4",
        "PCDH18","POSTN","P2RY14","ECM1","OLFML3","COL15A1","SVEP1","TEK",
        "COL16A1","PODN","HEYL","ALPL","HSPB7","MXRA8","LDB2","STEAP4","HTRA3",
        "MSX1","SOD3","PDGFRA","SPARCL1","ELN","ABCC9","COL1A2","RARRES2",
        "AQP1","FBLN2","CXCL12","C1S","KCNJ8","COG7","COX7A1","HSPB6","MFGE8",
        "ANPEP","ADM","HTRA1","IFITM1","MRGPRF","COL4A2","COL4A1","EDNRA",
        "INPP4B","CDH11","ANGPT2","MMP2","CDH5","TAGLN","PTH1R","ICAM1","ROBO4",
        "THY1","C1QTNF5","HSPB2","CRYAB","CSPG4","ISLR","RASL12","TBX18","NT5E",
        "ZIC1","PLSCR4","BGN","LAMA2","LAMA4","COL6A2","COL6A1","TIMP3","DCN",
        "LUM","GLI1","NDUFA4L2","SGCD","SLIT3","ADAMTS2","PECAM1","EMID1",
        "SPARC","VTN","COL1A1","RAMP2","AOC3","HIGD1B","CYGB","C1QTNF1",
        "COLEC11","CLEC14A","DLK1","FOXC1","ASPN","ECM2","PDZD2","ANGPT1",
        "COL14A1","GPIHBP1","ACVRL1","ATP13A5","FSTL1","ADAMTS5","FNDC1",
        "THBS2","DAAM2","NOTCH3","VEGFA","PDGFRB","CD248","EFEMP2","ACTA2",
        "COL6A3","ABCA8","MXRA5","C1R","HSPB2-C11ORF52","MEDAG","TNS2","ADGRA2",
        "ADGRL4","ADGRF5","CCN2","PLPP3","LHFPL6")

Epen<- c("PLXNB2","PLTP","C20ORF85","EFHC1","APOBEC4","PCP4L1","PTPRT","ARMC3",
         "ENKUR","MORN5","TTLL9","MEIG1","ANGPTL2","WDR38","HDC","SPEF1",
         "TMEM212","ZBBX","SPAG17","ADH7","STOML3","TM4SF1","CELSR2","FAM166B",
         "TEKT2","CROCC","CHEK2","CCDC60","CRYGN","RARRES2","AQP1","USP18",
         "LRRC23","VWA3A","CXCL17","MYO16","TTC29","NWD1","HYDIN","TMEM231",
         "DYNLRB2","TRIM71","CCDC153","CALML4","ZMYND10","AKAP14","LRRIQ1",
         "MYB","RSPH4A","PPIL6","S100B","TMEM107","EFNB3","TEKT1","GAS2L2",
         "KRT15","GFAP","FOXJ1","CCDC40","AK7","HSPA2","SNTN","CAPSL","SPEF2",
         "PVALB","ODF3B","SPAG6","IQCG","EFCAB1","PACRG","TMEM232","TEKT4",
         "RSPH1","SPATS1","TCTE1","RFX2","SIX3","NME5","AQP4","SCGB1A1",
         "LRRC10B","RIIAD1","CYP2A13","DNAH7","DNAH9","C1ORF87","DNAI1",
         "C2ORF73","RGS22","DNAH12","FAM183A","WFDC6","DNAI2","C1ORF194",
         "DNAH6","BPIFA1","C5ORF49","DNAAF1","DNAH2","FAM81B","MAP3K19",
         "VWA3B","AK8","C9ORF24","C6ORF118","LRRC71","PPP1R42","NME9","BPIFB1",
         "MS4A8","DRC1","PIFO","EPPIN-WFDC6","CATIP","TEX26","FAM216B",
         "C11ORF97","ANKRD66","CFAP44","CFAP70","ERICH3","CFAP43","CFAP126",
         "CFAP61","CFAP54","CFAP221","CFAP52","CFAP99","DRC7","CFAP47","SAXO2",
         "CFAP100","CFAP77","CFAP65","DRC3","CFAP206")

DefaultAssay (seurat_obj) <- "RNA"

#Ast, Endo, Mic, Neu, Oli, OPCs, Peri, Epen

celltype<-'Epen'

BrainSUBSET <- subset(seurat_obj, features=get(celltype))
av.expBrainSUBSET <- AverageExpression(BrainSUBSET)$SCT
cal_z_score <- function(x){(x - mean(x)) / sd(x)}
av.expBrainSUBSET <- t(apply(av.expBrainSUBSET, 1, cal_z_score))

library(pheatmap)
plot<-pheatmap (av.expBrainSUBSET, cluster_rows = F, fontsize_row = 7, main= paste(celltype))
#---------------------------------------------------------------------------------

# add celltype annotations
Idents(seurat_obj)<-"seurat_clusters"
new.cluster.ids<- c("Oli","Oli","Oli","Oli","Oli","Endo/Peri","Mic","Mic","Neu","Oli","Neu","Neu","Oli","Mic","Ast","Neu","Neu","OPCs","Ast","Ast","Ast","Ast","Oli","Neu","Mic","Ast","OPCs","Mic","Mic","OPCs","Neu","Neu","OPCs","OPCs","Ast","OPCs","Ast","Peri","Neu","Neu","Unk")

names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$celltype <- Idents(seurat_obj)

p1 <- DimPlot(seurat_obj, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(seurat_obj, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

#Batch_1_pre_processing_3 up to here 1/10/23

#load("~/multiome/batch_1_pre_processing_3_celltype-IDed.RData")

#link genes to peaks
DefaultAssay(seurat_obj) <- "peaks"

# first compute the GC content for each peak
library(BSgenome.Hsapiens.UCSC.hg38)
seurat_obj <- RegionStats(seurat_obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
#will take a long time
seurat_obj <- LinkPeaks(
  object = seurat_obj,
  peak.assay = "peaks",
  expression.assay = "SCT"
)

#Batch_1_pre_processing_3_peaks_linked up to here 1/10/23
