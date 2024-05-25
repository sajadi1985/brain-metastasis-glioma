library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
setwd(".......")
######################################Initial Data######################################
###Loading data
#ALL-GSE186344
breast1 <- Read10X(data.dir = "./186344/breast1")
dim(breast1)
#33694  6798
breast2 <- Read10X(data.dir ="./186344/breast2")
dim(breast2)
#33694 12003
breast3 <- Read10X(data.dir ="./186344/breast3")
dim(breast3)
#33694 11947
#breast4-GSE143423
breast4<-read.csv(file = "./GSE143423_tnbc_scRNAseq_gene_expression_counts.csv")
row.names(breast4) <- breast4[,1]
breast4 <- breast4[,-1]
dim(breast4)
#19826  7247
#GSE117891-glioma
glioma1<-read.csv(file = "./all_6148.umi.count.matrix.tsv",header = TRUE, sep = "\t")
row.names(glioma1) <- glioma1[,1]
glioma1 <- glioma1[,-1]
glioma1<-data.frame(glioma1)
dim(glioma1)
#24990  6148
glioma1N <- glioma1 %>%  select(starts_with(c("S1P2","S1P3","S1P8","S2P1","S2P3","S3P1","S3P8","S4P1","S4P3","S5P5","S6P5","S7P4","S8P5","S11P2","S11P3","S12P5","S13P1","S13P2","S13P3","S15P1","S15P2")))
dim(glioma1N)
#24990  1262
P1 <- colnames(glioma1N)
glioma1T <- select(glioma1,-P1)
#GSE202371-glima
glioma2 <-read.csv(file = "./GSE202371/GSM6112133_GB01.counts.tsv",header = TRUE, sep = "\t")
row.names(glioma2) <- glioma2[,1]
glioma2 <- glioma2[,-1]
dim(glioma2)
#30431  5470
glioma3<-read.csv(file = "./GSE202371/GSM6112134_GB02.counts.tsv",header = TRUE, sep = "\t")
row.names(glioma3) <- glioma3[,1]
glioma3 <- glioma3[,-1]
dim(glioma3)
#27263  3825
glioma4<-read.csv(file = "./GSE202371/GSM6112135_GB03.counts.tsv",header = TRUE, sep = "\t")
row.names(glioma4) <- glioma4[,1]
glioma4 <- glioma4[,-1]
dim(glioma4)
#32393  5552
glioma5<-read.csv(file = "./GSE202371/GSM6112136_GB04.counts.tsv",header = TRUE, sep = "\t")
row.names(glioma5) <- glioma5[,1]
glioma5 <- glioma5[,-1]
dim(glioma5)
# 30561  5142
#GSE135045 -glioma
glioma6 <- read.delim(file= "./GSE135045/GSM3984317_NO.1.expression_matrix.txt")
dim(glioma6)
#35966  4466
glioma7 <- read.delim(file= "./GSE135045/GSM3984318_NO.2.expression_matrix.txt")
dim(glioma7)
#28928  2643
glioma8 <- read.delim(file= "./GSE135045/GSM3984330_NO.16.expression_matrix.txt")
dim(glioma8)
#31166  2926


############################## seurat object #######################################

breast1 <- CreateSeuratObject(counts = breast1, project = "breast1", min.cells = 3, min.features = 200)
breast2 <- CreateSeuratObject(counts = breast2, project = "breast2", min.cells = 3, min.features = 200)
breast3 <- CreateSeuratObject(counts = breast3, project = "breast3", min.cells = 3, min.features = 200)
breast4 <- CreateSeuratObject(counts = breast4, project = "breast4", min.cells = 3, min.features = 200)
p2 <- c("tnbc1","tnbc6k")
p3 <- c("breast4","breast4")
Idents(breast4) <- plyr::mapvalues(x = Idents(breast4), from = p2, to = p3)
glioma1T <- CreateSeuratObject(counts = glioma1T, project = "glioma1T", min.cells = 3, min.features = 200)
p4<-Idents(glioma1T)
length(p4)
p5<-rep("glioma1T",4886)
Idents(glioma1T) <- plyr::mapvalues(x = Idents(glioma1T), from = p4, to = p5)
glioma2 <- CreateSeuratObject(counts = glioma2, project = "glioma2", min.cells = 3, min.features = 200)
glioma3 <- CreateSeuratObject(counts = glioma3, project = "glioma3", min.cells = 3, min.features = 200)
glioma4 <- CreateSeuratObject(counts = glioma4, project = "glioma4", min.cells = 3, min.features = 200)
glioma5 <- CreateSeuratObject(counts = glioma5, project = "glioma5", min.cells = 3, min.features = 200)
glioma6 <- CreateSeuratObject(counts = glioma6, project = "glioma6", min.cells = 3, min.features = 200)
glioma7 <- CreateSeuratObject(counts = glioma7, project = "glioma7", min.cells = 3, min.features = 200)
glioma8 <- CreateSeuratObject(counts = glioma8, project = "glioma8", min.cells = 3, min.features = 200)
################################################################################################
breast1[["percent.mt"]] <- PercentageFeatureSet(breast1, pattern = "^MT-")
breast2[["percent.mt"]] <- PercentageFeatureSet(breast2, pattern = "^MT-")
breast3[["percent.mt"]] <- PercentageFeatureSet(breast3, pattern = "^MT-")
breast4[["percent.mt"]] <- PercentageFeatureSet(breast4, pattern = "^MT-")
glioma1T[["percent.mt"]] <- PercentageFeatureSet(glioma1T, pattern = "^MT-")
glioma2[["percent.mt"]] <- PercentageFeatureSet(glioma2, pattern = "^MT-")
glioma3[["percent.mt"]] <- PercentageFeatureSet(glioma3, pattern = "^MT-")
glioma4[["percent.mt"]] <- PercentageFeatureSet(glioma4, pattern = "^MT-")
glioma5[["percent.mt"]] <- PercentageFeatureSet(glioma5, pattern = "^MT-")
glioma6[["percent.mt"]] <- PercentageFeatureSet(glioma6, pattern = "^MT-")
glioma7[["percent.mt"]] <- PercentageFeatureSet(glioma7, pattern = "^MT-")
glioma8[["percent.mt"]] <- PercentageFeatureSet(glioma8, pattern = "^MT-")

#################### quality control #####################

pdf("./Results/quality1_vlnplot-breast.pdf", width=20)
lapply(c(breast1,breast2,breast3,breast4,glioma1T,glioma2,glioma3,glioma4,glioma5,glioma6,glioma7,glioma8),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

breast1 <- subset(breast1, subset = nFeature_RNA <6500 & percent.mt < 25)
dim(breast1)
# 21175  6671
breast2 <- subset(breast2, subset = nFeature_RNA <5500 & percent.mt < 25)
dim(breast2)
#  20099 11005
breast3 <- subset(breast3, subset = nFeature_RNA <6500 & percent.mt < 25)
dim(breast3)
#20782 10923
breast4 <- subset(breast4, subset = nFeature_RNA <8000 & percent.mt < 25)
dim(breast4)
#115658  5820
glioma1T <- subset(glioma1T, subset = nFeature_RNA <10000 & percent.mt < 20)
dim(glioma1T)
#22791  4845
glioma2 <- subset(glioma2, subset = nFeature_RNA <6000 & percent.mt < 25)
dim(glioma2)
glioma3 <- subset(glioma3, subset = nFeature_RNA <4000 & percent.mt < 25)
dim(glioma3)
glioma4 <- subset(glioma4, subset = nFeature_RNA <6500 & percent.mt < 25)
dim(glioma4)
glioma5 <- subset(glioma5, subset = nFeature_RNA <6000 & percent.mt < 25)
dim(glioma5)
glioma6 <- subset(glioma6, subset = nFeature_RNA <5000 & percent.mt < 25)
dim(glioma6)
glioma7 <- subset(glioma7, subset = nFeature_RNA <4000 & percent.mt < 20)
dim(glioma7)
glioma8 <- subset(glioma8, subset = nFeature_RNA <4000 & percent.mt < 20)
dim(glioma8)
pdf("./Results/quality2_vlnplot-breast.pdf", width=20)
lapply(c(breast1,breast2,breast3,breast4,glioma1T,glioma2,glioma3,glioma4,glioma5,glioma6,glioma7,glioma8),VlnPlot,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#####################################  Normalization    ####################################################
all.list <- c(breast1,breast2,breast3,breast4,glioma1T,glioma2,glioma3,glioma4,glioma5,glioma6,glioma7,glioma8)
names(all.list) <- c("breast1","breast2","breast3","breast4","glioma1T","glioma2","glioma3","glioma4","glioma5","glioma6","glioma7","glioma8")
for (i in names(all.list)) {
  all.list[[i]] <- NormalizeData(all.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}
#####################################   FindVariableFeatures     ###############################################################
for (i in names(all.list)) {
  all.list[[i]] <- FindVariableFeatures(all.list[[i]], selection.method = "vst", nfeatures = 2000)
}
save.image(file = "FindVariableFeatures-breast.RData")
######################################  integration    ######################################################
n=35
all.anchors <- FindIntegrationAnchors(object.list=all.list,dims = 1:n)
all.integrated <- IntegrateData(anchorset = all.anchors, dims = 1:n)
dim(all.integrated )
# 2000 55749
save.image(file = "integrated-breast-Data.RData")
################################   Dimention Reduction& Clustering   ###################################################

DefaultAssay(all.integrated) <- "integrated"
all.integrated <- ScaleData(all.integrated, verbose = FALSE)
all.integrated <- RunPCA(all.integrated, npcs = 50, verbose = FALSE)
dataIDs <- Idents(all.integrated)
pdf("./pca-breast.pdf", width=10)
Idents(all.integrated) <- dataIDs
DimPlot(all.integrated, reduction = "pca")
dev.off()
pdf("./elbow-breast.pdf", width=10)
ElbowPlot(all.integrated, ndims = 40, reduction = "pca")
dev.off()
all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:35)

###  Clustering   #########################################################
all.integrated <- FindNeighbors(all.integrated, dims = 1:35)
all.integrated <- FindClusters(all.integrated, resolution = 0.2)
table(Idents(all.integrated))

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
classIDs<- Idents(all.integrated)
#########################################################
Idents(all.integrated) <- dataIDs
current.cluster.ids <- c("breast1","breast2","breast3","breast4","glioma1T","glioma2","glioma3","glioma4","glioma5","glioma6","glioma7","glioma8")
new.cluster.ids <- c("MB-breast","MB-breast","MB-breast","MB-breast","GLIOMA","GLIOMA","GLIOMA","GLIOMA","GLIOMA","GLIOMA","GLIOMA","GLIOMA")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
dataIDs2<- Idents(all.integrated)

Idents(all.integrated)<-classIDs
clusdata <-table(Idents(all.integrated), dataIDs2)
write.csv(clusdata, file= "./clusdata_integrated-breast.csv")

pdf("./integrated-breast_UMAP1.pdf", width=10)
Idents(all.integrated) <- dataIDs2
DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE,
        order = c("GLIOMA","MB-LUNG"),
        cols= c( "maroon1","yellow3")
) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
Idents(all.integrated) <- classIDs
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


save.image(file = "integrated-breast_UMAP.RData")


################cell types###############

current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12)
new.cluster.ids <- c("epithelial", " macrophage ", "OPC", "proliferating", "fibroblast", "oligodendrocyte","B cell", "T cell", "astrocyte","endothelial","undeterminde","monocyte")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
cellIDs<-Idents(all.integrated)
my_levels <-  c("epithelial", " macrophage ", "OPC", "proliferating", "fibroblast", "oligodendrocyte","B cell", "T cell", "astrocyte","endothelial","undeterminde","monocyte")
all.integrated@active.ident <- factor(x = all.integrated@active.ident, levels = my_levels) 


pdf("./integrated-breast cell type_UMAP.pdf", width=10)
Idents(all.integrated) <- cellIDs
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))

dev.off()

#########################################  DEG and MARKERS ################################################

DefaultAssay(all.integrated) <- "RNA"

all.markers <- FindAllMarkers(all.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./DEGs_integrated-breast.csv")
##################
pdf("./Markers_epithelial.pdf", width=4)
DotPlot(all.integrated, features = c("KRT19","KRT18","CDH1"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_macrophage.pdf", width=4)
DotPlot(all.integrated, features = c("C1QA","CTSB",	"AIF1"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_OPC.pdf", width=4)
DotPlot(all.integrated, features = c("GPR17","SMOC1","FERMT1"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_prolifirating.pdf", width=4)
DotPlot(all.integrated, features = c("UBE2C","BIRC5","NUF2"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_fibroblast.pdf", width=4)
DotPlot(all.integrated, features = c("COL1A1","THY1","DCN"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_oligodendrocyte.pdf", width=4)
DotPlot(all.integrated, features = c("MOG","OLIG1","OLIG2"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_B cell.pdf", width=4)
DotPlot(all.integrated, features = c("JCHAIN","CD79A","MZB1"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_Tcell.pdf", width=4)
DotPlot(all.integrated, features = c("CD3D","CD3E","TRAC"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_astrocyte.pdf", width=4)
DotPlot(all.integrated, features = c("GFAP","SLC1A2","AQP4"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_endothelial.pdf", width=4)
DotPlot(all.integrated, features = c("PECAM1","CLDN5","RAMP2"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_monocyte.pdf", width=4)
DotPlot(all.integrated, features = c("S100A8","S100A9","CTSS"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

################ cell types assignment #################################################

current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12)
new.cluster.ids <- c("epithelial", " macrophage ", "OPC", "proliferating", "fibroblast", "oligodendrocyte","B cell", "T cell", "astrocyte","endothelial","undeterminde","monocyte/neutrophil")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
cellIDs<-Idents(all.integrated)
my_levels <-  c("epithelial", " macrophage ", "OPC", "proliferating", "fibroblast", "oligodendrocyte","B cell", "T cell", "astrocyte","endothelial","undeterminde","monocyte/neutrophil")
# Re-level object@ident
all.integrated@active.ident <- factor(x = all.integrated@active.ident, levels = my_levels) 


pdf("./integrated-breast cell type_UMAP.pdf", width=10)
Idents(all.integrated) <- cellIDs
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))

dev.off()
################################  macrophage sub-clustering  ################
DefaultAssay(all.integrated) <- "integrated"
Idents(all.integrated) <- classIDs
all.integrated$sub_dataID <- dataIDs2
CLUS2 <- subset(x = all.integrated, idents = c("2"), invert = FALSE)
CLUS2$sub_dataID
CLUS2 <- FindVariableFeatures(CLUS2, selection.method = "vst", nfeatures = 2000)
CLUS2 <- ScaleData(CLUS2, verbose = FALSE)
CLUS2<- RunPCA(CLUS2, npcs = 50, verbose = FALSE)

pdf("./pca-clus2.pdf", width=10)
DimPlot(CLUS2, reduction = "pca")
dev.off()
pdf("./elbow-clus2.pdf", width=10)
ElbowPlot(CLUS2, ndims = 40, reduction = "pca")
dev.off()
CLUS2 <- RunUMAP(CLUS2, reduction = "pca", dims = 1:25)
################
CLUS2 <- FindNeighbors(CLUS2, dims = 1:25)
CLUS2 <- FindClusters(CLUS2, resolution = 0.1)
class2IDs <- Idents(CLUS2)
table(Idents(CLUS2))

current.cluster.ids <- c(0, 1, 2,3)
new.cluster.ids <- c(1, 2, 3,4)
Idents(CLUS2) <- plyr::mapvalues(x = Idents(CLUS2), from = current.cluster.ids, to = new.cluster.ids)
class1IDs<- Idents(CLUS2)

clus1data <-table(Idents(CLUS2), CLUS2$sub_dataID)
write.csv(clus1data, file= "./clus2data_CLUS2.csv")

##############################################
my_cols1 <- c('1'='#B95FBB','2'='#E6C122','3'='#F68282','4'='#4B4BF7','5'='#1FA195','6'='#aeadb3')
my_cols2 <- c('myCAF1'='#B95FBB','myCAF2'='#E6C122','iCAF'='#F68282','dCAF'='#4B4BF7','apCAF'='#1FA195','tpCAF'='#aeadb3')

pdf("./integrated-clus2_UMAP1.pdf", width=10)
Idents(CLUS2) <- CLUS2$sub_dataID

DimPlot(CLUS2, reduction = "umap", group.by = 'ident', repel = TRUE,
        order = c("GLIOMA","MB-breast"),
        cols= c( "darkorange","deepskyblue")
) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
Idents(CLUS2) <- class1IDs
DimPlot(CLUS2, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6+cols=my_cols1)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
save.image(file = "integrated-breast-macrophage_UMAP.RData")


#######   DEG and  markers of macrophage sub-clusters  ################

DefaultAssay(CLUS2) <- "RNA"

all.markers <- FindAllMarkers(CLUS2, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./DEGs_CLUS2.csv")

pdf("./Markers_all macrophage.pdf", width=10)
FeaturePlot(CLUS2, features = c( "MKI67","CX3CR1", "TMEM119", "P2RY12","EMP3","LYZ"),cols=c("lightgrey","darkred"))
VlnPlot(CLUS2, c( "MKI67","CX3CR1", "TMEM119", "P2RY12","EMP3","LYZ"))
DotPlot(CLUS2, features = c( "MKI67","CX3CR1", "TMEM119", "P2RY12","EMP3","LYZ"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_Monocyte derived Macrophage(MDM).pdf", width=15)
FeaturePlot(CLUS2, features = c("EMP3", "TMEM119", "TYROBP", "LGALS1", "CD68","LYZ", "CTSL", "EMP3", "ITGAM", "CD14","IITGA4"),cols=c("lightgrey","darkred"))
VlnPlot(CLUS2, c("EMP3", "TMEM119", "TYROBP", "LGALS1", "CD68","LYZ", "CTSL", "EMP3", "ITGAM", "CD14","IITGA4"))
DotPlot(CLUS2, features = c("EMP3", "TMEM119", "TYROBP", "LGALS1", "CD68","LYZ", "CTSL", "EMP3", "ITGAM", "CD14","IITGA4")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_microgelia(MG).pdf", width=15)
FeaturePlot(CLUS2, features = c("CX3CR1","TMEM119","TREM2","CSF1R","P2RY12"),cols=c("lightgrey","darkred"))
VlnPlot(CLUS2, c("CX3CR1","TMEM119","TREM2","CSF1R","P2RY12"))
DotPlot(CLUS2, features = c("CX3CR1","TMEM119","TREM2","CSF1R","P2RY12"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


pdf("./Markers_prolifirating mac.pdf", width=15)
FeaturePlot(CLUS2, features = c("STMN1","MK167","TOP2A","CDK1"),cols=c("lightgrey","darkred"))
VlnPlot(CLUS2, c("STMN1","MK167","TOP2A","CDK1"))
DotPlot(CLUS2, features = c("STMN1","MK167","TOP2A","CDK1"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

################# cell type clus 2 #######
Idents(CLUS2)<-class1IDs
current.cluster.ids <- c(1, 2, 3, 4)
new.cluster.ids <- c("MG","MDM1","MDM2","Priliferating")
Idents(CLUS2) <- plyr::mapvalues(x = Idents(CLUS2), from = current.cluster.ids, to = new.cluster.ids)
cellIDs<-Idents(CLUS2)

pdf("./integrated-breast-MAC-cell type_UMAP.pdf", width=10)
Idents(CLUS2) <- cellIDs
DimPlot(CLUS2, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6+cols=my_cols2)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))

dev.off()

################################# fibroblast sub-clustering  #############################
load("integrated-breast_UMAP.RData")
DefaultAssay(all.integrated) <- "integrated"
Idents(all.integrated) <- classIDs
all.integrated$sub_dataID <- dataIDs2
CLUS5 <- subset(x = all.integrated, idents = c("5"), invert = FALSE)
CLUS5$sub_dataID
CLUS5 <- FindVariableFeatures(CLUS5, selection.method = "vst", nfeatures = 2000)
CLUS5 <- ScaleData(CLUS5, verbose = FALSE)
CLUS5<- RunPCA(CLUS5, npcs = 50, verbose = FALSE)

pdf("./pca-clus5.pdf", width=10)
DimPlot(CLUS5, reduction = "pca")
dev.off()
pdf("./elbow-clus5.pdf", width=10)
ElbowPlot(CLUS5, ndims = 40, reduction = "pca")
dev.off()
CLUS5 <- RunUMAP(CLUS5, reduction = "pca", dims = 1:20)
#########################################
CLUS5 <- FindNeighbors(CLUS5, dims = 1:20)
CLUS5 <- FindClusters(CLUS5, resolution = 0.2)
class5IDs <- Idents(CLUS5)
table(Idents(CLUS5))

current.cluster.ids <- c(0, 1, 2,3,4,5)
new.cluster.ids <- c(1, 2, 3,4,5,6)
Idents(CLUS5) <- plyr::mapvalues(x = Idents(CLUS5), from = current.cluster.ids, to = new.cluster.ids)
class1IDs<- Idents(CLUS5)

clus1data <-table(Idents(CLUS5), CLUS5$sub_dataID)
write.csv(clus1data, file= "./clus5data_CLUS5.csv")
##########################################################


my_cols1 <- c('1'='#B95FBB','2'='#E6C122','3'='#F68282','4'='#4B4BF7','5'='#1FA195','6'='red')
my_cols2 <- c('myCAF1'='#B95FBB','myCAF2'='#E6C122','iCAF'='#F68282','dCAF'='#4B4BF7','apCAF'='#1FA195','tpCAF'='red')

pdf("integrated-clus5_UMAP1.pdf", width=10)
Idents(CLUS5) <- CLUS5$sub_dataID

DimPlot(CLUS5, reduction = "umap", group.by = 'ident', repel = TRUE,
        order = c("GLIOMA","MB-LUNG"),
        cols= c( "darkorange","deepskyblue")
) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
Idents(CLUS5) <- class1IDs
DimPlot(CLUS5, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols = my_cols1)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
current.cluster.ids <- c(1, 2, 3, 4,5,6)
new.cluster.ids <- c("myCAF1","myCAF2","iCAF","dCAF","apCAF","tpCAF")
Idents(CLUS5) <- plyr::mapvalues(x = Idents(CLUS5), from = current.cluster.ids, to = new.cluster.ids)
cellIDs<-Idents(CLUS5)

Idents(CLUS5) <- cellIDs
DimPlot(CLUS5, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols =my_cols2 )+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
save.image(file = "integrated_FIBROBLAST-breast-UMAP1.RData")

##################################### DEG and markers for fibroblast sub-clustering  #############

DefaultAssay(CLUS5) <- "RNA"
all.markers <- FindAllMarkers(CLUS5, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./DEGs_CLUS5.csv")

pdf("./Markers-myCAF1.pdf", width=10)
FeaturePlot(CLUS5, features = c("MMP11","POSTN","COL1A2","COMP","CDH11"),cols=c("lightgrey","#B95FBB"))
VlnPlot(CLUS5, c("MMP11","POSTN","COL1A2","COMP","CDH11"))
DotPlot(CLUS5, features = c("MMP11","POSTN","COL1A2","COMP","CDH11"),cols=c("lightgrey","#B95FBB")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers-myCAF2.pdf", width=10)
FeaturePlot(CLUS5, features = c("MMP11","POSTN","COL1A2","COMP","CDH11"),cols=c("lightgrey","#E6C122"))
VlnPlot(CLUS5, c("MMP11","POSTN","COL1A2","COMP","CDH11"))
DotPlot(CLUS5, features = c("MMP11","POSTN","COL1A2","COMP","CDH11"),cols=c("lightgrey","#E6C122")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers-i CAF.pdf", width=10)
FeaturePlot(CLUS5, features = c("PLA2G2A","CD34","CXCL12","ALDH1A1","KLF4","LEPR"),cols=c("lightgrey","#F68282"))
VlnPlot(CLUS5, c("PLA2G2A","CD34","CXCL12","ALDH1A1","KLF4","LEPR"))
DotPlot(CLUS5, features = c("PLA2G2A","CD34","CXCL12","ALDH1A1","KLF4","LEPR"),cols=c("lightgrey","#F68282")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


pdf("./Markers-tp CAF.pdf", width=10)
FeaturePlot(CLUS5, features = c("NT5E","HSPH1","CD10","TMEM158","MME","PDPN"),cols=c("lightgrey","red"))
VlnPlot(CLUS5, c("NT5E","HSPH1","CD10","TMEM158","MME","PDPN"))
DotPlot(CLUS5, features = c("NT5E","HSPH1","CD10","TMEM158","MME","PDPN"),cols=c("lightgrey","red")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers-ap CAF.pdf", width=10)
FeaturePlot(CLUS5, features = c("HLA-DRA","HLA-DRB1","CD74","CLU"),cols=c("lightgrey","#1FA195"))
VlnPlot(CLUS5, c("HLA-DRA","HLA-DRB1","CD74","CLU"))
DotPlot(CLUS5, features = c("HLA-DRA","HLA-DRB1","CD74","CLU"),cols=c("lightgrey","#1FA195")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers-d CAF.pdf", width=10)
FeaturePlot(CLUS5, features = c("TUBA1B","MKI67","CD74","CLU"),cols=c("lightgrey","#4B4BF7"))
VlnPlot(CLUS5, c("TUBA1B","MKI67","CD74","CLU"))
DotPlot(CLUS5, features = c("TUBA1B","MKI67","CD74","CLU"),cols=c("lightgrey","#4B4BF7")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
