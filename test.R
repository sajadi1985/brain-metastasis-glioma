library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
######################################Initial Data######################################
setwd("C:/Users/bio/Desktop/data-brain")
###Loading data####
lungtest1 <- read.delim(file= "./GSM6112137_LC01.counts.tsv")
row.names(lungtest1) <- lungtest1[,1]
lungtest1 <- lungtest1[,-1]
dim(lungtest1)
#32233  5537
############################## seurat object #######################################
lungtest <- CreateSeuratObject(counts = lungtest1, project = "lungtest", min.cells = 3, min.features = 200)
#################################################
lungtest[["percent.mt"]] <- PercentageFeatureSet(lungtest, pattern = "^MT-")

######################### quality control ########################################
pdf("./Results/quality1_vlnplot-lungtest.pdf", width=20)
VlnPlot(lungtest, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

lungtest <- subset(lungtest, subset = nFeature_RNA <5500 & percent.mt < 15)
dim(lungtest)
#26020  5438
pdf("./Results/quality2_vlnplot-lung.pdf", width=20)
VlnPlot(lungtest,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
################################  Normalization     ######################################
lungtest<- NormalizeData(lungtest, normalization.method = "LogNormalize", scale.factor = 10000)

#########################################    FindVariableFeatures  ######################################
lungtest <- FindVariableFeatures(lungtest, selection.method = "vst", nfeatures = 2000)

############################################ Dimension Reduction  #####################################################

lungtest <- ScaleData(lungtest, verbose = FALSE)
lungtest <- RunPCA(lungtest, npcs = 50, verbose = FALSE)
dataIDs <- Idents(lungtest)
pdf("./pca-lung-test.pdf", width=10)
Idents(lungtest) <- dataIDs
DimPlot(lungtest, reduction = "pca")
dev.off()
pdf("./elbow-lung-test.pdf", width=10)
ElbowPlot(lungtest, ndims = 40, reduction = "pca")
dev.off()
lungtest <- RunUMAP(lungtest, reduction = "pca", dims = 1:30)
###################### Clustering  ###############################################
lungtest <- FindNeighbors(lungtest, dims = 1:30)
lungtest <- FindClusters(lungtest, resolution = 0.5)
dataClass <- Idents(lungtest)

pdf("./lungtest_UMAP1.pdf", width=10)
Idents(lungtest) <- dataIDs
DimPlot(lungtest, reduction = "umap", group.by = 'ident', repel = TRUE
) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
Idents(lungtest) <- dataClass
DimPlot(lungtest, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()




#########################################  DEG and MARKERS  ###################################

DefaultAssay(lungtest) <- "RNA"

pdf("./lung test/Markers_epithelial.pdf", width=15)
Idents(lungtest) <- dataIDs
FeaturePlot(lungtest, features = c("EPCAM","KRT19","KRT18","FGG","SCGB3A1"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(lungtest, c("EPCAM","KRT19","KRT18","FGG","SCGB3A1"))
DotPlot(lungtest, features =  c("EPCAM","KRT19","KRT18","FGG","SCGB3A1"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./lung test/Markers_CD8 EXHUSTED.pdf", width=15)
Idents(lungtest) <- dataIDs
FeaturePlot(lungtest, features = c("CTLA4","CRTAM","DUSP4","CXCL13"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(lungtest, c("CTLA4","CRTAM","DUSP4","CXCL13"))
DotPlot(lungtest, features =  c("CTLA4","CRTAM","DUSP4","CXCL13"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./lung test/Markers_cytotoxic CD4+ T CELL.pdf", width=15)
Idents(lungtest) <- dataIDs
FeaturePlot(lungtest, features = c("IL7R","LTB","IL32","GZMA","GZMB","NKG7","CCL5","GNLY"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(lungtest, c("IL7R","LTB","IL32","GZMA","GZMB","NKG7","CCL5","GNLY"))
DotPlot(lungtest, features =  c("IL7R","LTB","IL32","GZMA","GZMB","NKG7","CCL5","GNLY"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./lung test/Markers_AT2 epithelial.pdf", width=15)
Idents(lungtest) <- dataIDs
FeaturePlot(lungtest, features = c("GPRC5A",	"SLC34A2"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(lungtest, c("GPRC5A",	"SLC34A2"))
DotPlot(lungtest, features =  c("GPRC5A",	"SLC34A2"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./lung test/Markers_macrophage.pdf", width=15)
Idents(lungtest) <- dataIDs
FeaturePlot(lungtest, features = c("C1QA",	"C1QB",	"FCGR3A","AIF1"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(lungtest, c("C1QA",	"C1QB",	"FCGR3A","AIF1"))
DotPlot(lungtest, features =  c("C1QA",	"C1QB",	"FCGR3A","AIF1"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./lung test/Markers_OPC.pdf", width=15)
Idents(lungtest) <- dataIDs
FeaturePlot(lungtest, features = c("PDGFRA",	"VCAN",	"FGAP","BCAN","CLU"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(lungtest, c("PDGFRA",	"VCAN",	"FGAP","BCAN","CLU"))
DotPlot(lungtest, features =  c("PDGFRA",	"VCAN",	"FGAP","BCAN","CLU"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./lung test/marker-oligo.pdf", width=15)
Idents(lungtest) <- dataIDs
FeaturePlot(lungtest, features = c("CLDN11",	"MOG",	"MOBP","PLP1","MAG","TF"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(lungtest, c("CLDN11",	"MOG",	"MOBP","PLP1","MAG","TF"))
DotPlot(lungtest, features =  c("CLDN11",	"MOG",	"MOBP","PLP1","MAG","TF"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


pdf("./lung test/marker.macroglia.pdf", width=15)
Idents(lungtest) <- dataIDs
FeaturePlot(lungtest, features = c("CD81",	"SPP1",	"TREM2"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(lungtest, c("CD81",	"SPP1",	"TREM2"))
DotPlot(lungtest, features =  c("CD81",	"SPP1",	"TREM2"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

