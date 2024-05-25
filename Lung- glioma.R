library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
setwd(".....")
######################################Initial Data######################################
###Loading data####
lung1 <- Read10X(data.dir = "./186344/lung1")
dim(lung1)
# 33694  9985
lung2 <- Read10X(data.dir = "./186344/lung2")
dim(lung2)
# 33694  3891
lung3 <- Read10X(data.dir = "./186344/lung3")
dim(lung3)
#33694  7337
#lung4-GSE131907
lung4 <- read.csv(file= "./data_lung4.csv")
row.names(lung4) <- lung4[,1]
lung4<- lung4[,-1]
dim(lung4)
# 29634  3222
#lung5-GSE143423
lung5 <- read.csv(file = "./GSE143423_lbm_scRNAseq_gene_expression_counts.csv", header=T)
lung5 <- na.omit(lung5)
row.names(lung5) <- lung5[,1]
lung5<- lung5[,-1]
lung5 <- lung5 %>% select(starts_with("lbm1"))
dim(lung5)
#13375  4610
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
lung1 <- CreateSeuratObject(counts = lung1, project = "lung1", min.cells = 3, min.features = 200)
lung2 <- CreateSeuratObject(counts = lung2, project = "lung2", min.cells = 3, min.features = 200)
lung3 <- CreateSeuratObject(counts = lung3, project = "lung3", min.cells = 3, min.features = 200)
lung4 <- CreateSeuratObject(counts = lung4, project = "lung4", min.cells = 3, min.features = 200)
lung5 <- CreateSeuratObject(counts = lung5, project = "lung5", min.cells = 3, min.features = 200)
p2 <- c("lbm1")
p3 <- c("lung5")
Idents(lung5) <- plyr::mapvalues(x = Idents(lung5), from = p2, to = p3)
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
####################################################################################
lung1[["percent.mt"]] <- PercentageFeatureSet(lung1, pattern = "^MT-")
lung2[["percent.mt"]] <- PercentageFeatureSet(lung2, pattern = "^MT-")
lung3[["percent.mt"]] <- PercentageFeatureSet(lung3, pattern = "^MT-")
lung4[["percent.mt"]] <- PercentageFeatureSet(lung4, pattern = "^MT-")
lung5[["percent.mt"]] <- PercentageFeatureSet(lung5, pattern = "^MT-")
glioma1T[["percent.mt"]] <- PercentageFeatureSet(glioma1T, pattern = "^MT-")
glioma2[["percent.mt"]] <- PercentageFeatureSet(glioma2, pattern = "^MT-")
glioma3[["percent.mt"]] <- PercentageFeatureSet(glioma3, pattern = "^MT-")
glioma4[["percent.mt"]] <- PercentageFeatureSet(glioma4, pattern = "^MT-")
glioma5[["percent.mt"]] <- PercentageFeatureSet(glioma5, pattern = "^MT-")
glioma6[["percent.mt"]] <- PercentageFeatureSet(glioma6, pattern = "^MT-")
glioma7[["percent.mt"]] <- PercentageFeatureSet(glioma7, pattern = "^MT-")
glioma8[["percent.mt"]] <- PercentageFeatureSet(glioma8, pattern = "^MT-")

###################################### quality control ######################################

pdf("./Results/quality1_vlnplot-lung.pdf", width=20)
lapply(c(lung1,lung2,lung3,lung4,lung5,glioma1T,glioma2,glioma3,glioma4,glioma5,glioma6,glioma7,glioma8),VlnPlot,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

lung1 <- subset(lung1, subset = nFeature_RNA <7500 & percent.mt < 15)
dim(lung1)
lung2 <- subset(lung2, subset = nFeature_RNA <8000 & percent.mt < 25)
dim(lung2)
#20804  3476
lung3 <- subset(lung3, subset = nFeature_RNA <8000 & percent.mt < 25)
dim(lung3)
# 20813  7202
lung4 <- subset(lung4, subset = nFeature_RNA <8000 & percent.mt < 20)
dim(lung4)
#20504  3202
lung5 <- subset(lung5, subset = nFeature_RNA <6000 & percent.mt < 5)
dim(lung5)
#11683  4602
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
pdf("./Results/quality2_vlnplot-lung.pdf", width=20)
lapply(c(lung1,lung2,lung3,lung4,lung5,glioma1T,glioma2,glioma3,glioma4,glioma5,glioma6,glioma7,glioma8),VlnPlot,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
########################################### Normalization  ##############################################
all.list <- c(lung1,lung2,lung3,lung4,lung5,glioma1T,glioma2,glioma3,glioma4,glioma5,glioma6,glioma7,glioma8)
names(all.list) <- c("lung1","lung2","lung3","lung4","lung5","glioma1T","glioma2","glioma3","glioma4","glioma5","glioma6","glioma7","glioma8")
for (i in names(all.list)) {
  all.list[[i]] <- NormalizeData(all.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}
############################################ FindVariableFeatures ########################################################
for (i in names(all.list)) {
  all.list[[i]] <- FindVariableFeatures(all.list[[i]], selection.method = "vst", nfeatures = 2000)
}
save.image(file = "FindVariableFeatures-lung.RData")
################################## integration #######################################
n=35
all.anchors <- FindIntegrationAnchors(object.list=all.list,dims = 1:n)
all.integrated <- IntegrateData(anchorset = all.anchors, dims = 1:n)
dim(all.integrated )
# 2000 55749
save.image(file = "integrated-lung-Data.RData")
################################ Dimention Reduction ######################################
DefaultAssay(all.integrated) <- "integrated"
all.integrated <- ScaleData(all.integrated, verbose = FALSE)
all.integrated <- RunPCA(all.integrated, npcs = 50, verbose = FALSE)
dataIDs <- Idents(all.integrated)
pdf("./pca-lung.pdf", width=10)
Idents(all.integrated) <- dataIDs
DimPlot(all.integrated, reduction = "pca")
dev.off()
pdf("./elbow-lung.pdf", width=10)
ElbowPlot(all.integrated, ndims = 40, reduction = "pca")
dev.off()
all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:35)

####################### Clustering  #################################
all.integrated <- FindNeighbors(all.integrated, dims = 1:35)
all.integrated <- FindClusters(all.integrated, resolution = 0.3)
table(Idents(all.integrated))

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8,9,10,11,12,13,14)
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15)
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
classIDs<- Idents(all.integrated)
#########################################################
Idents(all.integrated) <- dataIDs
current.cluster.ids <- c("lung1","lung2","lung3","lung4","lung5","glioma1T","glioma2","glioma3","glioma4","glioma5","glioma6","glioma7","glioma8")
new.cluster.ids <- c("MB-LUNG","MB-LUNG","MB-LUNG","MB-LUNG","MB-LUNG","GLIOMA","GLIOMA","GLIOMA","GLIOMA","GLIOMA","GLIOMA","GLIOMA","GLIOMA")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
dataIDs2<- Idents(all.integrated)

Idents(all.integrated)<-classIDs
clusdata <-table(Idents(all.integrated), dataIDs2)
write.csv(clusdata, file= "./clusdata_integrated-lung.csv")
###################################
my_cols <- c('3'='#F68282','15'='#31C53F','5'='#1FA195','1'='#B95FBB','13'='#D4D915',
             '14'='#28CECA','9'='#ff9a36','8'='#2FF18B','11'='#aeadb3','6'='#faf4cf',
             '2'='#CCB1F1','12'='#25aff5','7'='#A4DFF2','4'='#4B4BF7','16'='#AC8F14',
             '10'='#E6C122')
pdf("./integrated-lung_UMAP2.pdf", width=10)
Idents(all.integrated) <- dataIDs2
DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE,
        order = c("GLIOMA","MB-LUNG"),
        cols= c( "maroon1","deepskyblue")
) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
Idents(all.integrated) <- classIDs
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols=my_cols )+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


save.image(file = "integrated-lung_UMAP.RData")


############################################# DEG and MARKERS ############################################
DefaultAssay(all.integrated) <- "RNA"

all.markers <- FindAllMarkers(all.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./DEGs_integrated-lung.csv")

pdf("./Markers_all-lung.pdf", width=15)
DotPlot(all.integrated, features =  c("EPCAM","KRT19","KRT18","AIF1","CTSB","C1QB","MKI67","IGFBPL1","GPRC5A", "NAPSA", "SLC34A2","MOG","CLDN11","MOBP","IL7R","LTB","IQCG","PIFO","NME5","CD8A","GZMA","CCL5","JCHA1N","MZB1","IGHG3","DCN","THY1","COL1A1","GFAP","SLC1A3","AQP4","S100A9","CXCL8","PECAM1","CLDN5","FLT1","HLA-DQB1","HLA-DRB1","HLA-DPB1"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_epithelial.pdf", width=15)
FeaturePlot(all.integrated, features = c("EPCAM","KRT19","KRT18","CDH1"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(all.integrated, c("EPCAM","KRT19","KRT18","CDH1"))
DotPlot(all.integrated, features =  c("EPCAM","KRT19","KRT18","CDH1"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_macrophage.pdf", width=15)
FeaturePlot(all.integrated, features = c("AIF1","LYZ",	"LGMN",	"CTSB",	"CD14",	"HLA-DRA","HLA-DQA","HLA-DPB","C1QA","C1QB","FCGR3A"),cols=c("lightgrey","darkred"))
VlnPlot(all.integrated, c("AIF1","LYZ",	"LGMN",	"CTSB",	"CD14",	"HLA-DRA","HLA-DQA","HLA-DPB","C1QA","C1QB","FCGR3A"))
DotPlot(all.integrated, features = c("AIF1","LYZ",	"LGMN",	"CTSB",	"CD14",	"HLA-DRA","HLA-DQA","HLA-DPB","C1QA","C1QB","FCGR3A"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_Proliferation CELL.pdf", width=12)
FeaturePlot(all.integrated, features = c("MKI67"),cols=c("lightgrey","darkred"))
VlnPlot(all.integrated, c("MKI67"))
DotPlot(all.integrated, features = c("MKI67"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_epithelial(AT2).pdf", width=12)
FeaturePlot(all.integrated, features = c("SFTPB","SFTPC","ABCA3","GPRC5A","SFTPA1", "NAPSA", "SLC34A2", "SLP1"),cols=c("lightgrey","darkred"))
VlnPlot(all.integrated, c("SFTPB","SFTPC","ABCA3","GPRC5A","SFTPA1", "NAPSA", "SLC34A2", "SLP1"	))
DotPlot(all.integrated, features = c("SFTPB","SFTPC","ABCA3","GPRC5A","SFTPA1", "NAPSA", "SLC34A2", "SLP1"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_oligodendrocytes.pdf", width=15)
FeaturePlot(all.integrated, features = c("OLIG1","MOBP","OLIG2","MOG","CLDN11"),cols=c("lightgrey","pink"))
VlnPlot(all.integrated, c("OLIG1","MOBP","OLIG2","MOG","CLDN11"))
DotPlot(all.integrated, features = c("OLIG1","MOBP","OLIG2","MOG","CLDN11"),cols=c("lightgrey","pink")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_CD4 T cell.pdf", width=15)
FeaturePlot(all.integrated, features = c("LTB", "IL7R"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(all.integrated, c("LTB", "IL7R"))
DotPlot(all.integrated, features =  c("LTB", "IL7R"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_astrosyte.pdf", width=15)
FeaturePlot(all.integrated, features = c("IQCG","PIFO","NME5","RSPH1","CAPS","S100B","CRYAB","GFAP"),cols=c("lightgrey","darkred"))
VlnPlot(all.integrated, c("IQCG","PIFO","NME5","RSPH1","CAPS","S100B","CRYAB","GFAP"))
DotPlot(all.integrated, features = c("IQCG","PIFO","NME5","RSPH1","CAPS","S100B","CRYAB","GFAP"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_CD8 T cell.pdf", width=12)
FeaturePlot(all.integrated, features = c("CD8D","GZMA","CD3G","CCL5" ),cols=c("lightgrey","blue"))
VlnPlot(all.integrated, c("CD8D","GZMA","CD3G","CCL5" ))
DotPlot(all.integrated, features = c("CD8D","GZMA","CD3G","CCL5" ),cols=c("lightgrey","blue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_B cells.pdf", width=15)
FeaturePlot(all.integrated, features = c("JCHAIN","MZB1","CD79A","IGHM","IGHG3","IGHA2"),cols=c("lightgrey","darkred"))
VlnPlot(all.integrated, c("JCHAIN","MZB1","CD79A","IGHM","IGHG3","IGHA2"))
DotPlot(all.integrated, features = c("JCHAIN","MZB1","CD79A","IGHM","IGHG3","IGHA2"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


pdf("./Markers_fibroblast.pdf", width=15)
FeaturePlot(all.integrated, features = c("DCN",	"THY1",	"COL1A1","COL1A2"),cols=c("lightgrey","darkgreen"))
VlnPlot(all.integrated, c("DCN","THY1",	"COL1A1","COL1A2"))
DotPlot(all.integrated, features = c("DCN","THY1","COL1A1","COL1A2"),cols=c("lightgrey","darkgreen")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_OPC.pdf", width=15)
FeaturePlot(all.integrated, features = c("GFAP","SLC1A3","AQP4"),cols=c("lightgrey","darkgreen"))
VlnPlot(all.integrated, c("GFAP","SLC1A3","AQP4"))
DotPlot(all.integrated, features = c("GFAP","SLC1A3","AQP4"),cols=c("lightgrey","darkgreen")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


pdf("./Markers_monocyte.pdf", width=15)
FeaturePlot(all.integrated, features = c("S100A9","CXCL8","LYZ","VCAN","LGMN","CTSB","CD14"),cols=c("lightgrey","darkred"))
VlnPlot(all.integrated, c("S100A9","CXCL8","LYZ","VCAN","LGMN","CTSB","CD14"))
DotPlot(all.integrated, features = c("S100A9","CXCL8","LYZ","VCAN","LGMN","CTSB","CD14"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_endothelial.pdf", width=15)
FeaturePlot(all.integrated, features = c("PECAM1","CLDN5","FLT1","RAMP2"),cols=c("lightgrey","red"))
VlnPlot(all.integrated, c("PECAM1","CLDN5","FLT1","RAMP2"))
DotPlot(all.integrated, features = c("PECAM1","CLDN5","FLT1","RAMP2"),cols=c("lightgrey","red")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_Dentritic cell.pdf", width=15)
FeaturePlot(all.integrated, features = c("HLA-DQB1","HLA-DRB1","HLA-DPB1"),cols=c("lightgrey","darkred"))
VlnPlot(all.integrated, c("HLA-DQB1","HLA-DRB1","HLA-DPB1"))
DotPlot(all.integrated, features = c("HLA-DQB1","HLA-DRB1","HLA-DPB1"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()





#######################  cell types   assignment  ######################################################

current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15)
new.cluster.ids <- c("epithelial","macrophage", "proliferating", "epithelial(AT2)", "oligodendrocyte", "CD4 T cell","astrocyte", "CD8 T cell", "B cell","fibroblast","OPC","monocyte","undetermined","endothelial","dentritic cell")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
cellIDs<-Idents(all.integrated)
my_cols <- c('epithelial'='#B95FBB','macrophage'='#CCB1F1','proliferating'='#F68282','epithelial(AT2)'='#4B4BF7','oligodendrocyte'='#1FA195','CD4 T cell'='#faf4cf','astrocyte'='#A4DFF2','CD8 T cell'='#2FF18B','B cell'='#ff9a36','fibroblast'='#E6C122','OPC'='#aeadb3','monocyte'='#25aff5','undetermined'='#D4D915','endothelial'='#28CECA','dentritic cell'='#31C53F')
pdf("./integrated-lung cell type_UMAP2.pdf", width=10)
Idents(all.integrated) <- cellIDs
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols=my_cols)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))

dev.off()



################################  macrophage sub-clusteing  ################

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
CLUS2 <- RunUMAP(CLUS2, reduction = "pca", dims = 1:35)
###############################################################
CLUS2 <- FindNeighbors(CLUS2, dims = 1:35)
CLUS2 <- FindClusters(CLUS2, resolution = 0.2)
class2IDs <- Idents(CLUS2)
table(Idents(CLUS2))

current.cluster.ids <- c(0, 1, 2, 3, 4)
new.cluster.ids <- c(1, 2, 3, 4, 5)
Idents(CLUS2) <- plyr::mapvalues(x = Idents(CLUS2), from = current.cluster.ids, to = new.cluster.ids)
class1IDs<- Idents(CLUS2)

clus1data <-table(Idents(CLUS2), CLUS2$sub_dataID)
write.csv(clus1data, file= "./clus2data_CLUS2.csv")
##########################################################


my_cols1 <- c('1'='#B95FBB','2'='#E6C122','3'='#F68282','4'='#4B4BF7','5'='#1FA195')
my_cols2 <- c('MDM1'='#B95FBB','MG1'='#E6C122','MG2'='#F68282','MG3'='#4B4BF7','MDM2'='#1FA195')

pdf("integrated-clus2_UMAP1.pdf", width=10)
Idents(CLUS2) <- CLUS2$sub_dataID

DimPlot(CLUS2, reduction = "umap", group.by = 'ident', repel = TRUE,
        order = c("GLIOMA","MB-LUNG"),
        cols= c( "maroon1","deepskyblue")
) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
Idents(CLUS2) <- class1IDs
DimPlot(CLUS2, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols = my_cols1)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
current.cluster.ids <- c(1, 2, 3, 4,5)
new.cluster.ids <- c("MDM1","MG1","MG2","MG3","MDM2")
Idents(CLUS2) <- plyr::mapvalues(x = Idents(CLUS2), from = current.cluster.ids, to = new.cluster.ids)
cellIDs<-Idents(CLUS2)

Idents(CLUS2) <- cellIDs
DimPlot(CLUS2, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols =my_cols2 )+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
save.image(file = "integrated_Macrophage-lung-UMAP1.RData")

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



#############################CLUSTERING  CD4 T CELL###############################################
DefaultAssay(all.integrated) <- "integrated"
Idents(all.integrated) <- classIDs
all.integrated$sub_dataID <- dataIDs2
CLUS6 <- subset(x = all.integrated, idents = c("6"), invert = FALSE)
CLUS6$sub_dataID
CLUS6 <- FindVariableFeatures(CLUS6, selection.method = "vst", nfeatures = 2000)
CLUS6 <- ScaleData(CLUS6, verbose = FALSE)
CLUS6<- RunPCA(CLUS6, npcs = 50, verbose = FALSE)

pdf("./pca-clus6.pdf", width=10)
DimPlot(CLUS6, reduction = "pca")
dev.off()
pdf("./elbow-clus6.pdf", width=10)
ElbowPlot(CLUS6, ndims = 40, reduction = "pca")
dev.off()
CLUS6 <- RunUMAP(CLUS6, reduction = "pca", dims = 1:25)
################
CLUS6 <- FindNeighbors(CLUS6, dims = 1:25)
CLUS6 <- FindClusters(CLUS6, resolution = 0.2)
class6IDs <- Idents(CLUS6)
table(Idents(CLUS6))

current.cluster.ids <- c(0, 1, 2,3,4)
new.cluster.ids <- c(1, 2, 3,4,5)
Idents(CLUS6) <- plyr::mapvalues(x = Idents(CLUS6), from = current.cluster.ids, to = new.cluster.ids)
class6IDs<- Idents(CLUS6)

clus6data <-table(Idents(CLUS6), CLUS6$sub_dataID)
write.csv(clus6data, file= "./clus6data_CLUS6.csv")


############################ DEG and markers for CD4+ T cell sub-clustering################################

DefaultAssay(CLUS6) <- "RNA"

all.markers <- FindAllMarkers(CLUS6, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./DEGs_CLUS6.csv")
save.image(file = "integrated_lung-cd4-tcell-UMAP1.RData")

pdf("./Markers_CD4+ Cytotoxic 1.pdf", width=15)
FeaturePlot(CLUS6, features = c("IFNG", "CD4", "FGFBP2", "ITGB1", "GZMA", "CST7", "GNLY"),cols=c("lightgrey","#E6C122"))
VlnPlot(CLUS6, c("IFNG", "CD4", "FGFBP2", "ITGB1", "GZMA", "CST7", "GNLY"))
DotPlot(CLUS6, features = c("IFNG", "CD4", "FGFBP2", "ITGB1", "GZMA", "CST7", "GNLY"),cols=c("lightgrey","#E6C122")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_CD4+ Cytotoxic 2.pdf", width=15)
FeaturePlot(CLUS6, features = c("IFNG", "CD4", "FGFBP2", "ITGB1", "GZMA", "CST7", "GNLY"),cols=c("lightgrey","#F68282"))
VlnPlot(CLUS6, c("IFNG", "CD4", "FGFBP2", "ITGB1", "GZMA", "CST7", "GNLY"))
DotPlot(CLUS6, features = c("IFNG", "CD4", "FGFBP2", "ITGB1", "GZMA", "CST7", "GNLY"),cols=c("lightgrey","#F68282")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_CD4+ Proliferating T.pdf", width=15)
FeaturePlot(CLUS6, features = c("MKI67", "TOP2A", "CDK1", "CCNB1", "KIF4A" , "E2F1", "RRM2"),cols=c("lightgrey","#1FA195"))
VlnPlot(CLUS6, c("MKI67", "TOP2A", "CDK1", "CCNB1", "KIF4A" , "E2F1", "RRM2"))
DotPlot(CLUS6, features = c("MKI67", "TOP2A", "CDK1", "CCNB1", "KIF4A" , "E2F1", "RRM2"),cols=c("lightgrey","#1FA195")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


pdf("./Markers_CD4+ Naive  T.pdf", width=15)
FeaturePlot(CLUS6, features = c("TCF7", "CCR7", "FHIT", "LEF1", "MAL", "NOSIP","CCNB1"),cols=c("lightgrey","#B95FBB"))
VlnPlot(CLUS6, c("TCF7", "CCR7", "FHIT", "LEF1", "MAL", "NOSIP","CCNB1"))
DotPlot(CLUS6, features = c("TCF7", "CCR7", "FHIT", "LEF1", "MAL", "NOSIP","CCNB1"),cols=c("lightgrey","#B95FBB")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()



pdf("./Markers_CD4 TCM.pdf", width=15)
FeaturePlot(CLUS6, features = c("IL7R", "TOB1","ITGB1","AQP3", "PDE4D"),cols=c("lightgrey","darkred"))
VlnPlot(CLUS6, c("IL7R", "TOB1","ITGB1","AQP3", "PDE4D"))
DotPlot(CLUS6, features = c("IL7R", "TOB1","ITGB1","AQP3", "PDE4D"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_CD4 TEM.pdf", width=15)
FeaturePlot(CLUS6, features = c("CCL5","LTB", "GZMK", "GZMA", "KLRB1", "AQP3","DUSP2"),cols=c("lightgrey","darkred"))
VlnPlot(CLUS6, c("CCL5","LTB", "GZMK", "GZMA", "KLRB1", "AQP3","DUSP2"))
DotPlot(CLUS6, features = c("CCL5","LTB", "GZMK", "GZMA", "KLRB1", "AQP3","DUSP2"),cols=c("lightgrey","darkred")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


pdf("./Markers_Regulatory T.pdf", width=15)
FeaturePlot(CLUS6, features = c( "RTKN2", "FOXP3", "AC133644.2", "CD4", "IL2RA", "TIGIT", "CTLA4"),cols=c("lightgrey","#4B4BF7"))
VlnPlot(CLUS6, c( "RTKN2", "FOXP3", "AC133644.2", "CD4", "IL2RA", "TIGIT", "CTLA4"))
DotPlot(CLUS6, features = c( "RTKN2", "FOXP3", "AC133644.2", "CD4", "IL2RA", "TIGIT", "CTLA4"),cols=c("lightgrey","#4B4BF7")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
save.image(file = "integrated_lung-CD4-UMAP1.RData")


############################  cell type assignment CCD4+ T cell sub-clusters ######################
Idents(CLUS6)<-class6IDs
current.cluster.ids <- c(1, 2, 3, 4, 5)
new.cluster.ids <- c("naive","cytotoxic1","cytotoxic2","regulatory(Treg)","Proliferating")
Idents(CLUS6) <- plyr::mapvalues(x = Idents(CLUS6), from = current.cluster.ids, to = new.cluster.ids)
cellIDs<-Idents(CLUS6)



my_cols1 <- c('naive'='#B95FBB','cytotoxic1'='#E6C122','cytotoxic2'='#F68282','regulatory(Treg)'='#4B4BF7','Proliferating'='#1FA195')
my_cols2 <- c('1'='#B95FBB','2'='#E6C122','3'='#F68282','4'='#4B4BF7','5'='#1FA195')

pdf("./integrated-clus6_UMAP1.pdf", width=10)
Idents(CLUS6) <- CLUS6$sub_dataID

DimPlot(CLUS6, reduction = "umap", group.by = 'ident', repel = TRUE,
        order = c("GLIOMA","MB-LUNG"),
        cols= c( "maroon1","deepskyblue")
) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))

Idents(CLUS6) <- class6IDs
DimPlot(CLUS6, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols=my_cols2)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))

Idents(CLUS6) <- cellIDs
DimPlot(CLUS6, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols=my_cols1)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))

dev.off()

############################# CD8 T cell sub-clustering ###############################################
DefaultAssay(all.integrated) <- "integrated"
Idents(all.integrated) <- classIDs
all.integrated$sub_dataID <- dataIDs2
CLUS8 <- subset(x = all.integrated, idents = c("8"), invert = FALSE)
CLUS8$sub_dataID
CLUS8 <- FindVariableFeatures(CLUS8, selection.method = "vst", nfeatures = 2000)
CLUS8 <- ScaleData(CLUS8, verbose = FALSE)
CLUS8<- RunPCA(CLUS8, npcs = 50, verbose = FALSE)

pdf("./pca-clus8.pdf", width=10)
DimPlot(CLUS8, reduction = "pca")
dev.off()
pdf("./elbow-clus8.pdf", width=10)
ElbowPlot(CLUS8, ndims = 40, reduction = "pca")
dev.off()
CLUS8 <- RunUMAP(CLUS8, reduction = "pca", dims = 1:35)
################
CLUS8 <- FindNeighbors(CLUS8, dims = 1:35)
CLUS8 <- FindClusters(CLUS8, resolution = 0.2)
class8IDs <- Idents(CLUS8)
table(Idents(CLUS8))

current.cluster.ids <- c(0, 1, 2,3,4)
new.cluster.ids <- c(1, 2, 3,4,5)
Idents(CLUS8) <- plyr::mapvalues(x = Idents(CLUS8), from = current.cluster.ids, to = new.cluster.ids)
class8IDs<- Idents(CLUS8)

clus8data <-table(Idents(CLUS8), CLUS8$sub_dataID)
write.csv(clus8data, file= "./clus8data_CLUS8.csv")


############################### DEG and markers for CD8+ T cell sub-clustering #############################################
DefaultAssay(CLUS8) <- "RNA"

all.markers <- FindAllMarkers(CLUS8, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./DEGs_CLUS8.csv")

save.image(file = "integrated_lung-cd8-tcell-UMAP1.RData")

pdf("./Markers_CD8+ cytotoxic.pdf", width=15)
FeaturePlot(CLUS8, features = c("GNLY", "ANXA1", "CD8A", "KRT1", "LINC02446", "YBX3","LDHB"),cols=c("lightgrey","#F68282"))
VlnPlot(CLUS8, c("GNLY", "ANXA1", "CD8A", "KRT1", "LINC02446", "YBX3","LDHB"))
DotPlot(CLUS8, features = c("GNLY", "ANXA1", "CD8A", "KRT1", "LINC02446", "YBX3","LDHB"),cols=c("lightgrey","#F68282")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_CD8+ Naive1  T.pdf", width=15)
FeaturePlot(CLUS8, features = c("TCF7", "CD4", "CCR7","IL7R", "TRAC", "NELL2", "LDHB"),cols=c("lightgrey","#B95FBB"))
VlnPlot(CLUS8, c("TCF7", "CD4", "CCR7" ,"IL7R", "TRAC", "NELL2", "LDHB"))
DotPlot(CLUS8, features = c("TCF7", "CD4", "CCR7","IL7R", "TRAC", "NELL2", "LDHB"),cols=c("lightgrey","#B95FBB")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Markers_CD8+ Naive2  T.pdf", width=15)
FeaturePlot(CLUS8, features = c("TCF7", "CD4", "CCR7","IL7R", "TRAC", "NELL2", "LDHB"),cols=c("lightgrey","#E6C122"))
VlnPlot(CLUS8, c("TCF7", "CD4", "CCR7" ,"IL7R", "TRAC", "NELL2", "LDHB"))
DotPlot(CLUS8, features = c("TCF7", "CD4", "CCR7","IL7R", "TRAC", "NELL2", "LDHB"),cols=c("lightgrey","#E6C122")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_CD8+ Proliferating T.pdf", width=15)
FeaturePlot(CLUS8, features = c("MKI67", "TOP2A", "CDK1", "CCNB1", "KIF4A" , "E2F1", "RRM2"),cols=c("lightgrey","#1FA195"))
VlnPlot(CLUS8, c("MKI67", "TOP2A", "CDK1", "CCNB1", "KIF4A" , "E2F1", "RRM2"))
DotPlot(CLUS8, features = c("MKI67", "TOP2A", "CDK1", "CCNB1", "KIF4A" , "E2F1", "RRM2"),cols=c("lightgrey","#1FA195")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_CD8+ exhausted T.pdf", width=15)
FeaturePlot(CLUS8, features = c("CTLA4","CDK1","NKG7", "PD1", "CST7", "CD38", "TCF1"),cols=c("lightgrey","#4B4BF7"))
VlnPlot(CLUS8, c("CTLA4","CDK1","NKG7", "PD1", "CST7", "CD38", "TCF1"))
DotPlot(CLUS8, features = c("CTLA4","CDK1","NKG7", "GZMK", "CST7", "CD38", "TCF1"),cols=c("lightgrey","#4B4BF7")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
######################### cell type assignment CCD8+ T cell sub-clusters ########################

Idents(CLUS8)<-class8IDs
current. cluster.ids <- c(1, 2, 3, 4, 5)
new. cluster.ids <- c("naive1", "naive2", "cytotoxic", "exhausted", "proliferating")
Idents(CLUS8) <- plyr::mapvalues(x = Idents(CLUS8), from = current.cluster.ids, to = new.cluster.ids)
cellIDs<-Idents(CLUS8)
my_cols1 <- c('naive1'='#B95FBB','naive2'='#E6C122','cytotoxic'='#F68282','exhausted'='#4B4BF7','proliferating'='#1FA195')
my_cols2 <- c('1'='#B95FBB','2'='#E6C122','3'='#F68282','4'='#4B4BF7','5'='#1FA195')
pdf("./integrated-clus8_UMAP1.pdf", width=10)
Idents(CLUS8) <- CLUS8$sub_dataID

DimPlot(CLUS8, reduction = "umap", group.by = 'ident', repel = TRUE,
        order = c("GLIOMA","MB-LUNG"),
        cols= c( "maroon1","deepskyblue")
) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))

Idents(CLUS8) <- class8IDs
DimPlot(CLUS8, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols=my_cols2)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))

Idents(CLUS8) <- cellIDs
DimPlot(CLUS8, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6,cols=my_cols1)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))

dev.off()
