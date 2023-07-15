library(Seurat)
library(dplyr)
library(clusterProfiler)
library(monocle3) 
library(reticulate)
library(ggplot2)
library(org.Hs.eg.db)

library(scales) # 用于Col函数查看并更改颜色

rm(list = ls())
gc()

setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
wd <- getwd()
time <- as.character(Sys.Date())
dir.create(time) 

# 载入organoid的原始数据------
setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
KPC.data <- Read10X(data.dir = "D:/Zuo Lab/Data/20210325 sequence Data/Result-X101SC19090611-Z01-F001-37-B1/2.Summary/KPC-org/filtered_feature_bc_matrix")
organoid <- CreateSeuratObject(counts = KPC.data, project = "KPC_organoid", min.cells = 3, min.features = 200)
organoid[["percent.mt"]] <- PercentageFeatureSet(organoid, pattern = "^MT-")
organoid <- subset(organoid, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
VlnPlot(organoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

organoid <- NormalizeData(organoid, normalization.method = "LogNormalize", scale.factor = 10000)
organoid <- FindVariableFeatures(organoid, selection.method = "vst", nfeatures = 2000)

rm(KPC.data)

# 载入Fetal kidney 3的测序数据集-----
  #读取Kidney数据
  # 读取txt格式的矩阵
  data1 <- read.table("2021-06-03/Fetal-Kidney3/Fetal-Kidney3_dge.txt", header = T, sep = ",")
  
  # txt矩阵转换成数据框
  datan = data.frame(data1)
  row.names(datan) <- datan$X
  datan <- datan[,-1]
  # 数据框装转换成稀疏矩阵matrix
  dataan <- as(as.matrix(datan),"dgCMatrix")
  Fetal <- CreateSeuratObject(counts = dataan, project = "KPC", min.cells = 3, min.features = 200)
  
  
  #标准化处理
  Fetal[["percent.mt"]] <- PercentageFeatureSet(Fetal, pattern = "^MT-") 
  Fetal <- NormalizeData(Fetal, normalization.method = "LogNormalize", scale.factor = 10000)
  Fetal <- FindVariableFeatures(Fetal, selection.method = "vst", nfeatures = 2000)
  
  rm(data1, datan, dataan)
# 载入Fetal kidney 4的测序数据集-----
#读取Kidney数据
# 读取txt格式的矩阵
  data1 <- read.table("2021-06-03/Fetal-Kidney4/Fetal-Kidney4_dge.txt", header = T, sep = ",")
  
  # txt矩阵转换成数据框
  datan = data.frame(data1)
  row.names(datan) <- datan$X
  datan <- datan[,-1]
  # 数据框装转换成稀疏矩阵matrix
  dataan <- as(as.matrix(datan),"dgCMatrix")
  Fetal <- CreateSeuratObject(counts = dataan, project = "fetal", min.cells = 3, min.features = 200)

  
  #标准化处理
  Fetal[["percent.mt"]] <- PercentageFeatureSet(Fetal, pattern = "^MT-") 
  Fetal <- NormalizeData(Fetal, normalization.method = "LogNormalize", scale.factor = 10000)
  Fetal <- FindVariableFeatures(Fetal, selection.method = "vst", nfeatures = 2000)
  rm(data1, datan, dataan)
# 载入Fetal kidney 5的测序数据集-----
  setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
  wd <- getwd()
  time <- as.character(Sys.Date())
  dir.create(time) 
  
  data1 <- read.table("2021-06-03/Fetal-Kidney5/Fetal-Kidney5_dge.txt", header = T, sep = ",")
  
  # txt矩阵转换成数据框
  datan = data.frame(data1)
  row.names(datan) <- datan$X
  datan <- datan[,-1]
  # 数据框装转换成稀疏矩阵matrix
  dataan <- as(as.matrix(datan),"dgCMatrix")
  Fetal <- CreateSeuratObject(counts = dataan, project = "fetal", min.cells = 3, min.features = 200)
  
  
  #标准化处理
  Fetal[["percent.mt"]] <- PercentageFeatureSet(Fetal, pattern = "^MT-") 
  Fetal <- NormalizeData(Fetal, normalization.method = "LogNormalize", scale.factor = 10000)
  Fetal <- FindVariableFeatures(Fetal, selection.method = "vst", nfeatures = 2000)
  
  rm(data1, datan, dataan)

# 载入Fetal kidney 6的测序数据集-----
  
data1 <- read.table("2021-06-03/Fetal-Kidney6/Fetal-Kidney6_dge.txt", header = T, sep = ",")
datan = data.frame(data1)
row.names(datan) <- datan$X
datan <- datan[,-1]
dataan <- as(as.matrix(datan),"dgCMatrix")
Fetal <- CreateSeuratObject(counts = dataan, project = "Fetal", min.cells = 3, min.features = 200)

Fetal[["percent.mt"]] <- PercentageFeatureSet(Fetal, pattern = "^MT-") 
Fetal <- NormalizeData(Fetal, normalization.method = "LogNormalize", scale.factor = 10000)
Fetal <- FindVariableFeatures(Fetal, selection.method = "vst", nfeatures = 2000)

rm(data1, datan, dataan)



# 载入Adult kidney 2的测序数据集-----
  #读取Kidney数据
  # 读取txt格式的矩阵
  data1 <- read.table("2021-06-03/Adult-Kidney2/Adult-Kidney2_dge.txt", header = T, sep = ",")
  
  # txt矩阵转换成数据框
  datan = data.frame(data1)
  row.names(datan) <- datan$X
  datan <- datan[,-1]
  # 数据框装转换成稀疏矩阵matrix
  dataan <- as(as.matrix(datan),"dgCMatrix")
  Adult <- CreateSeuratObject(counts = dataan, project = "KPC", min.cells = 3, min.features = 200)
  
  
  #标准化处理
  Adult[["percent.mt"]] <- PercentageFeatureSet(Adult, pattern = "^MT-") 
  Adult <- NormalizeData(Adult, normalization.method = "LogNormalize", scale.factor = 10000)
  Adult <- FindVariableFeatures(Adult, selection.method = "vst", nfeatures = 2000)
  
  rm(data1, datan, dataan)
# 载入D16 kidney organoid的测序数据集-----
  setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
  Ora.data <- Read10X(data.dir = "D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC/2022-10-07/HuOrg_D16_2")
  Ora <- CreateSeuratObject(counts = Ora.data, project = "Ora_referance", min.cells = 3, min.features = 200)
  Ora[["percent.mt"]] <- PercentageFeatureSet(Ora, pattern = "^MT-")
  Ora <- subset(Ora, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  VlnPlot(Ora, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  Ora <- NormalizeData(Ora, normalization.method = "LogNormalize", scale.factor = 10000)
  Ora <- FindVariableFeatures(Ora, selection.method = "vst", nfeatures = 2000)
  
  rm(Ora.data)
  

# 整合分析-----
  scell.anchors <- FindIntegrationAnchors(object.list = c(organoid,Fetal), dims = 1:20)
  scell.anchors <- FindIntegrationAnchors(object.list = c(organoid,Adult), dims = 1:20)
  scell.anchors <- FindIntegrationAnchors(object.list = c(organoid,Ora), dims = 1:20)
  
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(scell)
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  merge <- scell
  rm(scell,scell.anchors)

# 载入Monocle定义下的UMAP-----

  setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
  wd <- getwd()
  time <- as.character(Sys.Date())
  dir.create(time)
  load(file = "2022-07-24/organoid_monocle3_cds.Robj")
  
  cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.0003)
  cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.00045)
  plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
  
  load(file = "2022-07-17/KPC.Robj")
  DimPlot(KPC)
  KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
  colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
  KPC@active.ident <- KPC$Monocle.Cluster
  DimPlot(KPC, pt.size = 1)
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"LH-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Cycling cell"
  
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"LH-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Epithelial cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "6"))<-"Cycling cell"
  
  levels(KPC@active.ident)
  #仅区分分化与未分化的分群
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"Tubule-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"Tubule-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"Tubule-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Cycling cell"
  DimPlot(KPC, pt.size = 1)+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  KPC.markers <- FindAllMarkers(KPC, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(KPC.markers, file = paste(wd,"/",time,"/","differ_undiffer_DEG.csv", sep = ""), sep=",")
  
  organoid_rename <- KPC
  DimPlot(organoid_rename, pt.size = 1)
  rm(cds,KPC)
# 读取对应的参考UMAP数据提取定义-----
  load(file = "2021-06-05/Adult kidney2(redefine).Robj")
  Adult_rename <- scell
  DimPlot(Adult_rename,pt.size = 1)
  rm(scell)
  
  load(file = "2021-06-05/fetal kidney5(redefine).Robj")
  Fetal_rename <- scell
  DimPlot(Fetal_rename,pt.size = 1)
  rm(scell)
  
  load(file = "2021-06-05/fetal kidney4(redefine).Robj")
  Fetal_rename <- scell                                                                                   
  DimPlot(Fetal_rename,pt.size = 1)
  rm(scell)
  
  # 提取细胞
  
  load(file = "2021-06-05/fetal kidney3(redefine).Robj")
  Fetal_rename <- scell
  DimPlot(Fetal_rename,pt.size = 1)
  rm(scell)

# 将两种定义合并（通过as.factor转换为因子类型后加入到metadata中）-----
a <- as.data.frame(Fetal_rename@active.ident)
a[,2] <- rownames(a)
colnames(a) <- c("idents", "cell_id")

a <- as.data.frame(Adult_rename@active.ident)
a[,2] <- rownames(a)
colnames(a) <- c("idents", "cell_id")

b <- as.data.frame(organoid_rename@active.ident)
b[,2] <- rownames(b)
colnames(b) <- c("idents","cell_id")
c <- rbind(b,a)

rm(a,b)

metadata <- as.data.frame(merge@active.ident)
metadata[,2] <- rownames(metadata)
colnames(metadata) <- c("orig.ident", "cell_id")
metadata <- full_join(c,metadata, by=c("cell_id"))
rownames(metadata) <- metadata$cell_id

metadata$idents <- factor(metadata$idents)
merge@meta.data[["idents"]] <- metadata$idents
Idents(merge) <- merge@meta.data$idents

DimPlot(merge,pt.size = 1)
DimPlot(merge,pt.size = 1,group.by = "idents",split.by = "orig.ident")
DimPlot(merge,pt.size = 1,group.by = "idents",split.by = "orig.ident",cols = c("grey","grey","grey","grey","grey",
                                                                               "grey","grey","grey","grey","grey",
                                                                               "grey","grey","grey","grey","grey",
                                                                               "grey","grey"))
DimPlot(merge,pt.size = 1,group.by = "idents",split.by = "orig.ident",cells.highlight = rownames(merge@meta.data[which(merge$idents %in% c("LH-like","Loop of Henle progenitor cell")),]))


DimPlot(merge,pt.size = 0.7)
DimPlot(merge,pt.size = 0.7, split.by = "orig.ident")



# Fetal3
DefaultAssay(merge) <- "RNA"
levels(merge@active.ident)
Idents(merge,cells=WhichCells(merge,idents = c("Endothelial cell_PLVAP high","Endothelial cell_EMCN high","Endothelial cell_CCL21 high","Endothelial cell_GJA5 high"))) <- "Endothelial cell"
Idents(merge,cells=WhichCells(merge,idents = c("Mesangial cell","Nephrogenic mesenchyme cell_DAPL1 high"))) <- "Mesangial cell"
Idents(merge,cells=WhichCells(merge,idents = c("Interstitial progenitor cell","Interstitial cell_POSTN high"))) <- "Mesangial cell"
Idents(merge,cells=WhichCells(merge,idents = c("Lymphatic endothelial cell","Neutrophil","M2 Macrophage"))) <- "Immunocytes"

Idents(merge) <- merge@meta.data$idents

# Fetal4
DefaultAssay(merge) <- "RNA"
levels(merge@active.ident)
Idents(merge,cells=WhichCells(merge,idents = c("Endothelial cell_PLVAP high","Endothelial cell_GJA5 high"))) <- "Endothelial cell"
Idents(merge,cells=WhichCells(merge,idents = c("Interstitial cell_PTN high"))) <- "Mesangial cell"
Idents(merge,cells=WhichCells(merge,idents = c("Interstitial cell_POSTN high","Nephrogenic mesenchyme cell_DAPL1 high"))) <- "Mesangial cell"
Idents(merge,cells=WhichCells(merge,idents = c("S-shaped body cell_LINC01158 high","S-shaped body medial cell"))) <- "S-shaped body cell"
Idents(merge,cells=WhichCells(merge,idents = c("S-shaped body cell_CFAP126 high"))) <- "S-shaped body cell"
Idents(merge,cells=WhichCells(merge,idents = c("Collecting duct cell_CRABP1 high","Collecting duct cell_CALB1 high"))) <- "Collecting duct cell"

Idents(merge) <- merge@meta.data$idents
show_col(hue_pal()(17))
  # LH
    DimPlot(merge,pt.size = 0.7,cols = c("grey","grey","grey","grey","grey",
                                         "grey","#00BC51","grey","grey","#00BCD6",
                                         "grey","grey","grey","grey","grey",
                                         "grey","grey"))
  # PT
    DimPlot(merge,pt.size = 0.7, cols = c("grey","grey","grey","grey","grey",
                                          "grey","grey","#00C087","grey","grey",
                                          "grey","grey","grey","grey","#F166E8",
                                          "grey","grey"))
  # CD
  DimPlot(merge,pt.size = 0.7, cols = c("#F8766D","grey","grey","grey","grey",
                                        "grey","grey","grey","#00C0B2","grey",
                                        "grey","grey","grey","grey","grey",
                                        "grey","grey"))
  # Progenitor 
  DimPlot(merge,pt.size = 0.7, cols = c("grey","grey","grey","grey","grey",
                                        "#45B500","grey","grey","grey","grey",
                                        "grey","grey","grey","grey","grey",
                                        "grey","grey"))
  # Fetal5
  DefaultAssay(merge) <- "RNA"
  levels(merge@active.ident)
  Idents(merge,cells=WhichCells(merge,idents = c("Endothelial cell_PLVAP high","Endothelial cell_CCL21 high"))) <- "Endothelial cell"
  Idents(merge,cells=WhichCells(merge,idents = c("Interstitial cell_POSTN high","Nephrogenic mesenchyme cell_DAPL1 high"))) <- "Mesangial cell"
  Idents(merge,cells=WhichCells(merge,idents = c("S-shaped body cell","S-shaped body medial cell"))) <- "S-shaped body cell"
  Idents(merge,cells=WhichCells(merge,idents = c("S-shaped body cell_CFAP126 high"))) <- "S-shaped body cell"
  Idents(merge,cells=WhichCells(merge,idents = c("Collecting duct cell_CRABP1 high","Collecting duct cell_CALB1 high"))) <- "Collecting duct cell"
  
  Idents(merge) <- merge@meta.data$idents
  DimPlot(merge,pt.size = 0.7)
  DimPlot(merge,pt.size = 0.7,split.by = "orig.ident")
  DimPlot(merge,pt.size = 0.7, group.by = "orig.ident")
  show_col(hue_pal()(19))
  # LH
  DimPlot(merge,pt.size = 0.7,cols = c("grey","grey","grey","grey","grey",
                                       "#6FB000","grey","grey","#00C08E","grey",
                                       "grey","grey","grey","grey","grey",
                                       "grey","grey","grey","grey","grey"))
  # PT
  DimPlot(merge,pt.size = 0.7, cols = c("grey","grey","grey","grey","grey",
                                        "grey","#00B813","grey","grey","grey",
                                        "grey","grey","#00A7FF","grey","grey",
                                        "grey","grey","grey","grey","grey"))
  # CD
  DimPlot(merge,pt.size = 0.7, cols = c("grey","grey","grey","grey","grey",
                                        "grey","grey","#00BD61","grey","#00C0B4",
                                        "grey","grey","grey","grey","grey",
                                        "grey","grey","grey","grey","grey"))
  # Progenitor 
  DimPlot(merge,pt.size = 0.7, cols = c("grey","grey","grey","grey","#9CA700",
                                        "grey","grey","grey","grey","grey",
                                        "grey","grey","grey","grey","grey",
                                        "grey","grey","grey","grey","grey"))
  
  # Adult2
  DefaultAssay(merge) <- "RNA"
  levels(merge@active.ident)
  Idents(merge,cells=WhichCells(merge,idents = c("Loop of Henle(Thick ascending limb)","Loop of Henle_SPP1 high","Loop of Henle_SFN high"))) <- "Loop of Henle"
  Idents(merge,cells=WhichCells(merge,idents = c("Intercalated cell_SPINK1 high","Intercalated cell_SLC26A4 high"))) <- "Intercalated cell"
  Idents(merge,cells=WhichCells(merge,idents = c("T cell","B cell","Macrophage"))) <- "Immunocytes"
  Idents(merge,cells=WhichCells(merge,idents = c("Fenestrated endothelial cell_EMCN high","Fenestrated endothelial cell_SELE high"))) <- "Fenestrated endothelial cell"
  Idents(merge,cells=WhichCells(merge,idents = c("Proximal tubule cell_MT1G high","Proximal tubule cell_ALDOB high"))) <- "Proximal tubule"
  Idents(merge,cells=WhichCells(merge,idents = c("Intercalated cell","Principle cell"))) <- "Collecting duct"
  
  Idents(merge) <- merge@meta.data$idents
  DimPlot(merge,pt.size = 0.7)
  DimPlot(merge,pt.size = 0.7,split.by = "orig.ident")
  DimPlot(merge,pt.size = 0.7, group.by = "orig.ident")
  show_col(hue_pal()(20))
  # LH
  DimPlot(merge,pt.size = 0.7,cols = c("grey","grey","grey","grey","#A3A500",
                                       "grey","grey","#39B600","grey","grey",
                                       "grey","grey","grey","grey","grey",
                                       "grey","grey","grey","grey","grey"))
  # PT
  DimPlot(merge,pt.size = 0.7, cols = c("grey","#EA8331","grey","grey","grey",
                                        "grey","grey","grey","#00BF7D","grey",
                                        "grey","grey","grey","grey","grey",
                                        "grey","grey","grey","grey","grey"))
  # CD
  DimPlot(merge,pt.size = 0.7, cols = c("#F8766D","grey","grey","grey","grey",
                                        "grey","grey","grey","grey","#00C1A3",
                                        "grey","grey","grey","grey","grey",
                                        "grey","grey","grey","grey","grey"))
  # Progenitor 
  DimPlot(merge,pt.size = 0.7, cols = c("grey","grey","grey","grey","grey",
                                        "grey","#39B600","grey","grey","grey",
                                        "grey","grey","grey","grey","grey",
                                        "grey","grey","grey","grey","grey"))


# 以下为与空间转组拟合的数据分析-----
  setwd("D:/Zuo Lab/Data/spatial/KIDNEY/seurat/")
  wd <- getwd()
  time <- as.character(Sys.Date())
  dir.create(time) 
  
# 空间转录组数据-----
  spacial <-Seurat::Load10X_Spatial("data/spacial kidney/GSE171406_Human_Nephrectomy/")
  spacial <- SCTransform(spacial, assay = "Spatial", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  SpatialDimPlot(spacial)
  
# 载入Monocle定义下的UMAP-----
  
  setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
  wd <- getwd()
  time <- as.character(Sys.Date())
  dir.create(time)
  load(file = "2022-07-24/organoid_monocle3_cds.Robj")
  
  cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.0003)
  plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
  
  load(file = "2022-07-17/KPC.Robj")
  DimPlot(KPC)
  KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
  colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
  KPC@active.ident <- KPC$Monocle.Cluster
  DimPlot(KPC, pt.size = 1)
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"LH-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Cell cycle"
  
  KPC.markers <- FindAllMarkers(KPC, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(KPC.markers, file = paste(wd,"/",time,"/","DEG.csv", sep = ""), sep=",")
  colo <- RColorBrewer::brewer.pal(10,"RdYlBu")[1:5]
  colo2 <- RColorBrewer::brewer.pal(10,"RdYlBu")[6:10]
  n <- "IGFBP5"
  FeaturePlot(KPC,features = n,sort.cell = T,pt.size =1.5,min.cutoff = "q40")+
    scale_colour_gradient(low = rev(colo2),high = rev(colo))
  
  
  organoid_rename <- KPC
  DimPlot(organoid_rename, pt.size = 1)
  rm(cds,KPC)
  
  spacial_KPC_organoid <- spacial
  
  Idents(organoid_rename) 
  organoid_rename <- SCTransform(organoid_rename, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
  
  anchors <- FindTransferAnchors(reference = organoid_rename, query = spacial_KPC_organoid, normalization.method = "SCT", dims = 1:10)
  predictions.assay <- TransferData(anchorset = anchors, refdata = organoid_rename@active.ident, prediction.assay = TRUE, 
                                    weight.reduction = spacial_KPC_organoid[["pca"]],dims = 1:15)
  spacial_KPC_organoid[["predictions"]] <- predictions.assay 
  
  DefaultAssay(spacial_KPC_organoid) <- "predictions"

  
  data <- spacial_KPC_organoid@assays$predictions@data
  data  
  for(i in 1:ncol(data)){
    data[6,i] <- rownames(data)[which.max(data[,i])]
  }
  data1 <- t(data)
  all(colnames(data)==colnames(spacial_KPC_organoid))
  all(rownames(data1) == colnames(spacial_KPC_organoid))
  spacial_KPC_organoid[["label"]] <- data1[,6]
  Idents(spacial_KPC_organoid) <- spacial_KPC_organoid$label  
  SpatialDimPlot(spacial_KPC_organoid)
  
  SpatialFeaturePlot(spacial_KPC_organoid, features = c("Progenitor cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q80")
  SpatialFeaturePlot(spacial_KPC_organoid, features = c("LH-like"),alpha = c(0.1,1), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q80", max.cutoff = "q99")+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "top")

  SpatialFeaturePlot(spacial_KPC_organoid, features = c("PT-like"),alpha = c(0.1,1), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q80")+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "right")
  
  SpatialFeaturePlot(spacial_KPC_organoid, features = c("CD-like"),alpha = c(0.1,1), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q50")+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "right")
  
  SpatialDimPlot(spacial_KPC_organoid,group.by = "Idents")
  
  SpatialFeaturePlot(i, features = cell.type, 
                     pt.size.factor = 2,  crop = TRUE)+
    scale_fill_gradient( limits=c(0.01,0.1),breaks=c(0.01,0.05,0.1),
                         low = rev(colo2),high = rev(colo))
  
  
  
  
# Fetal kidney 5-----
  load(file = "2021-06-05/fetal kidney5(redefine).Robj")
  Fetal_rename <- scell
  DimPlot(Fetal_rename,pt.size = 1)
  rm(scell)
  
  Idents(Fetal_rename) 
  Fetal_rename <- SCTransform(Fetal_rename, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
  
  spacial_Fetal5 <- spacial
  anchors <- FindTransferAnchors(reference = Fetal_rename, query = spacial_Fetal5, normalization.method = "SCT", dims = 1:10)
  predictions.assay <- TransferData(anchorset = anchors, refdata = Fetal_rename@active.ident, prediction.assay = TRUE, 
                                    weight.reduction = spacial_Fetal5[["pca"]],dims = 1:15)
  spacial_Fetal5[["predictions"]] <- predictions.assay 
  
  DefaultAssay(spacial_Fetal5) <- "predictions"
  
  
  data <- spacial_Fetal5@assays$predictions@data
  data  
  for(i in 1:ncol(data)){
    data[19,i] <- rownames(data)[which.max(data[,i])]
  }
  data1 <- t(data)
  all(colnames(data)==colnames(spacial_Fetal5))
  all(rownames(data1) == colnames(spacial_Fetal5))
  spacial_Fetal5[["label"]] <- data1[,19]
  Idents(spacial_Fetal5) <- spacial_Fetal5$label  
  SpatialDimPlot(spacial_Fetal5)
  
  rownames(data)
  
  SpatialFeaturePlot(spacial_Fetal5, features = c("Interstitial progenitor cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
  SpatialFeaturePlot(spacial_Fetal5, features = c("Loop of Henle progenitor cell"),alpha = c(0.1,1), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q80")+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "right")
  
  SpatialFeaturePlot(spacial_Fetal5, features = c("Proximal tubule progenitor cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q80")
  SpatialFeaturePlot(spacial_Fetal5, features = c("Podocyte"),alpha = c(0.1,1), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q50", max.cutoff = "q90")+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "right")
  
  SpatialFeaturePlot(spacial_Fetal5, features = c("Proliferating cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
  
  
  
  
  
  n <- "NPHS1"
  FeaturePlot(organoid_rename,features = n,sort.cell = T,pt.size =1.5,min.cutoff = "q40")

# Fetal kidney 4-----
  load(file = "2021-06-05/fetal kidney4(redefine).Robj")
  Fetal_rename <- scell
  DimPlot(Fetal_rename,pt.size = 1)
  rm(scell)
  
  Idents(Fetal_rename) 
  Fetal_rename <- SCTransform(Fetal_rename, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
  
  spacial_Fetal4 <- spacial
  anchors <- FindTransferAnchors(reference = Fetal_rename, query = spacial_Fetal4, normalization.method = "SCT", dims = 1:10)
  predictions.assay <- TransferData(anchorset = anchors, refdata = Fetal_rename@active.ident, prediction.assay = TRUE, 
                                    weight.reduction = spacial_Fetal4[["pca"]],dims = 1:15)
  spacial_Fetal4[["predictions"]] <- predictions.assay 
  
  DefaultAssay(spacial_Fetal4) <- "predictions"
  
  
  data <- spacial_Fetal4@assays$predictions@data
  data  
  for(i in 1:ncol(data)){
    data[20,i] <- rownames(data)[which.max(data[,i])]
  }
  data1 <- t(data)
  all(colnames(data)==colnames(spacial_Fetal4))
  all(rownames(data1) == colnames(spacial_Fetal4))
  spacial_Fetal4[["label"]] <- data1[,20]
  Idents(spacial_Fetal4) <- spacial_Fetal4$label  
  SpatialDimPlot(spacial_Fetal4)
  
  rownames(data)
  
  SpatialFeaturePlot(spacial_Fetal4, features = c("S-shaped body medial cell"),alpha = c(0.1,1),min.cutoff = "q60", pt.size.factor = 1.6, ncol = 2, crop = TRUE)+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "right")
  
  SpatialFeaturePlot(spacial_Fetal4, features = c("Collecting duct cell-CALB1 high"),alpha = c(0.1,1), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q80")+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "right")
  
  SpatialFeaturePlot(spacial_Fetal4, features = c("Proximal tubule progenitor cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q80")
  SpatialFeaturePlot(spacial_Fetal4, features = c("Podocyte"),alpha = c(0.1,1), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q50", max.cutoff = "q90")+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "right")
  
  SpatialFeaturePlot(spacial_Fetal4, features = c("Proliferating cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
  
  
  
  
  
  n <- "NPHS1"
  FeaturePlot(organoid_rename,features = n,sort.cell = T,pt.size =1.5,min.cutoff = "q40")
  
# Fetal kidney 6-----
  load(file = "2021-06-05/fetal kidney6(redefine).Robj")
  Fetal_rename <- scell
  DimPlot(Fetal_rename,pt.size = 1)
  rm(scell)
  
  Idents(Fetal_rename) 
  Fetal_rename <- SCTransform(Fetal_rename, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
  
  spacial_Fetal6 <- spacial
  anchors <- FindTransferAnchors(reference = Fetal_rename, query = spacial_Fetal6, normalization.method = "SCT", dims = 1:10)
  predictions.assay <- TransferData(anchorset = anchors, refdata = Fetal_rename@active.ident, prediction.assay = TRUE, 
                                    weight.reduction = spacial_Fetal6[["pca"]],dims = 1:15)
  spacial_Fetal6[["predictions"]] <- predictions.assay 
  
  DefaultAssay(spacial_Fetal6) <- "predictions"
  
  
  data <- spacial_Fetal6@assays$predictions@data
  data  
  for(i in 1:ncol(data)){
    data[16,i] <- rownames(data)[which.max(data[,i])]
  }
  data1 <- t(data)
  all(colnames(data)==colnames(spacial_Fetal6))
  all(rownames(data1) == colnames(spacial_Fetal6))
  spacial_Fetal6[["label"]] <- data1[,16]
  Idents(spacial_Fetal6) <- spacial_Fetal6$label  
  SpatialDimPlot(spacial_Fetal6)
  
  rownames(data)
  
  SpatialFeaturePlot(spacial_Fetal6, features = c("S-shaped body cell"),min.cutoff = "q60", pt.size.factor = 1.6, ncol = 2, crop = TRUE)
  SpatialFeaturePlot(spacial_Fetal6, features = c("Collecting duct cell-CALB1 high"),alpha = c(0.1,1), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q80")+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "right")
  
  SpatialFeaturePlot(spacial_Fetal6, features = c("Interstitial progenitor cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q80")
  SpatialFeaturePlot(spacial_Fetal6, features = c("Podocyte"),alpha = c(0.1,1), pt.size.factor = 1.6, ncol = 2, crop = TRUE,min.cutoff = "q50", max.cutoff = "q90")+
    theme(legend.text=element_text(color="black",size=rel(1.5)),
          legend.title=element_text(color="black",size=rel(1.5)),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.7, "cm"),
          legend.position = "right")
  
  SpatialFeaturePlot(spacial_Fetal6, features = c("Proliferating cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
  
  
  
  
  
  n <- "NPHS1"
  FeaturePlot(organoid_rename,features = n,sort.cell = T,pt.size =1.5,min.cutoff = "q40")
  
# Adult kidney -----
  load(file = "2021-06-05/Adult")
  Adult_rename <- scell
  DimPlot(Adult_rename,pt.size = 1)
  rm(scell)
  
  
  
# GSEA分析-----
  library(tidyverse)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggplot2)
  library(enrichplot)
  
  setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
  msigdb_GMTs <- "msigdb_v7.0_GMTs"
  msigdb_H <- "2022-10-07/h.all.v2022.1.Hs.symbols.gmt"    #根据个人实际需求下载或编辑
  msigdb_GO <- "2022-10-07/c5.all.v2022.1.Hs.symbols.gmt"    #根据个人实际需求下载或编辑
  kegmt <- read.gmt(file.path(msigdb_GO))

  #准备差异基因列表
  setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
  DEG <- read.csv(file = "2022-10-07/DEG.csv", sep = ",")
  LH_like <- read.csv(file = "2022-10-07/LH-like.CSV",sep = ",", row.names = 1)
  Progenitor <- read.csv(file = "2022-10-07/progenitor_cell.CSV",sep = ",", row.names = 1)
  Tubule_like <-  read.csv(file = "2022-10-09/tubule_like_deg.CSV",sep = ",", row.names = 1)
  
  geneList = Tubule_like$avg_log2FC
  names(geneList) = Tubule_like$gene
  head(geneList)
  geneList = sort(geneList, decreasing = TRUE)
  
  #GSEA分析
  set.seed(1)
  egmt<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
  #转换成数据框
  egmt_result_df <- as.data.frame(egmt)
  write.csv(egmt_result_df, file = "Tubule_like_GSAE.CSV")
  write.table(egmt_result_df,file="GSEA_MSigDb_h.all_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
  save(egmt,egmt_result_df,file = "GSEA_deg_h.all.rda")  
  
  gseaplot2(egmt,7,color="red",pvalue_table = T)
  
  #结果汇总
  gseaplot2(egmt, geneSetID = c(1,3), subplots = 1:3,pvalue_table = T)
  
  #气泡图 展示geneset被激活还是抑制
  egmt2<- setReadable(egmt,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
  dotplot(egmt2,split=".sign",showCategory = 10,font.size = 10, title = "",
          label_format = 30,)+facet_grid(~.sign)
  #edit legends
  # +guides(
  #reverse color order (higher value on top)
  #color = guide_colorbar(reverse = TRUE))
  #reverse size order (higher diameter on top) 
  #size = guide_legend(reverse = TRUE))
  # Title 可以添加标题
  
  #### 经典的GSEA图 
  egmt$ID[]
  a <- which(egmt_result_df$ID %in% c("GOBP_ORGANIC_ACID_METABOLIC_PROCESS"))
  egmt_result_df$ID["GOBP_ORGANIC_ACID_METABOLIC_PROCESS"]
  which(x = egmt_result_df$ID,"GOBP_ORGANIC_ACID_METABOLIC_PROCESS")
  egmt@result[["ID"]][1]
  exp1=dplyr::filter(egmt_result_df, grepl("GOBP_ORGANIC_ACID_METABOLIC_PROCESS", ID))
  which(dplyr::filter(egmt_result_df, grepl("GOBP_ORGANIC_ACID_METABOLIC_PROCESS", ID)))
  egmt$Description
  i <-  which(egmt_result_df$ID %in% c("GOBP_ORGANIC_ACID_METABOLIC_PROCESS"))
  gseap1 <- gseaplot2(egmt,
                      egmt$ID[i],#富集的ID编号
                      title = egmt$Description[i],#标题
                      color = "red", #GSEA线条颜色
                      base_size = 23,#基础字体大小
                      rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                      subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                      ES_geom = "line", #enrichment score用线还是用点"dot"
                      pvalue_table = T) #显示pvalue等信息
  ggsave(gseap1, filename = 'GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION.pdf', width =14, height =10)
  
  #### 合并 GSEA通路 
  gseap2 <- gseaplot2(kk_gse,
                      up_gsea$ID,#富集的ID编号
                      title = "UP_GSEA_all",#标题
                      color = "red",#GSEA线条颜色
                      base_size = 20,#基础字体大小
                      rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                      subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                      ES_geom = "line",#enrichment score用线还是用点"dot"
                      pvalue_table = T) #显示pvalue等信息
  ggsave(gseap2, filename = "GSEA_up_all.pdf",width =12,height =12)

  
# 20221029 将单细胞测序的数据中小管细胞分离出来-----
  
  # Fetal4
  load(file = "2021-06-05/fetal kidney4(redefine).Robj")
  Fetal_rename <- scell                                                                                   
  DimPlot(Fetal_rename,pt.size = 1)
  rm(scell)
  
  DimPlot(Fetal_rename)
  Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Endothelial cell_PLVAP high","Endothelial cell_GJA5 high"))) <- "Endothelial cell"
  Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Interstitial cell_PTN high"))) <- "Mesangial cell"
  Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Interstitial cell_POSTN high","Nephrogenic mesenchyme cell_DAPL1 high"))) <- "Mesangial cell"
  Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("S-shaped body cell_LINC01158 high","S-shaped body medial cell"))) <- "S-shaped body cell"
  Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("S-shaped body cell_CFAP126 high"))) <- "S-shaped body cell"
  Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Collecting duct cell_CRABP1 high","Collecting duct cell_CALB1 high"))) <- "Collecting duct cell"
  levels(Fetal_rename)
  Fetal_rename <- subset(Fetal_rename, idents  = c("Collecting duct cell","Loop of Henle progenitor cell","Proximal tubule progenitor cell"))
  
  # KPC organoid
  setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
  wd <- getwd()
  time <- as.character(Sys.Date())
  dir.create(time)
  load(file = "2022-07-24/organoid_monocle3_cds.Robj")
  
  cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.0003)
  plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
  
  load(file = "2022-07-17/KPC.Robj")
  DimPlot(KPC)
  KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
  colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
  KPC@active.ident <- KPC$Monocle.Cluster
  DimPlot(KPC, pt.size = 1)
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"LH-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Cell cycle"
  levels(KPC@active.ident)
  KPC <- subset(KPC,idents = c("CD-like","PT-like","LH-like"))
  
  # 整合
  scell.anchors <- FindIntegrationAnchors(object.list = c(KPC,Fetal_rename), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(scell)
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  merge <- scell
  rm(scell,scell.anchors)
  
  
  a <- as.data.frame(Fetal_rename@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents", "cell_id")
  
  b <- as.data.frame(KPC@active.ident)
  b[,2] <- rownames(b)
  colnames(b) <- c("idents","cell_id")
  c <- rbind(b,a)
  
  rm(a,b)
  
  metadata <- as.data.frame(merge@active.ident)
  metadata[,2] <- rownames(metadata)
  colnames(metadata) <- c("orig.ident", "cell_id")
  metadata <- full_join(c,metadata, by=c("cell_id"))
  rownames(metadata) <- metadata$cell_id
  
  metadata$idents <- factor(metadata$idents)
  merge@meta.data[["idents"]] <- metadata$idents
  Idents(merge) <- merge@meta.data$idents
  DimPlot(merge,pt.size = 1)
  
  # 绘图
  DimPlot(merge,pt.size = 0.7)
  DimPlot(merge,pt.size = 0.7,split.by = "orig.ident")
  DimPlot(merge,pt.size = 0.7, group.by = "orig.ident")
  show_col(hue_pal()(6))
  # #F8766D #B79F00 #00BA3B #00BFC4 #619CFF #F564E3 
  # LH
  DimPlot(merge,pt.size = 0.7,cols = c("#F8766D","grey","grey","grey","#619CFF","grey"))
  # PT
  DimPlot(merge,pt.size = 0.7, cols = c("grey","#B79F00","grey","grey","grey","#F564E3"))
  # CD
  DimPlot(merge,pt.size = 0.7, cols = c("grey","grey","#00BA3B","#00BFC4","grey","grey"))
  
  
  DefaultAssay(merge) <- "RNA"
  merge@meta.data[["batch"]] <- merge$orig.ident
  merge@meta.data[["celltype"]] <- merge@active.ident
  
  # 转换成为loom文件
  library(loomR)
  library(SeuratDisk)
  main.loom <- as.loom(x = merge, filename = "merge.loom", verbose = FALSE)
  main.loom$close_all()
  write.csv(merge@meta.data,'2022-10-30/merge_metadata.csv')

  Convert("2022-10-29/pbmc.h5ad", dest = "h5seurat", overwrite = F)
  pbmc <-  LoadH5Seurat("2022-10-29/pbmc.h5seurat",meta.data = T) 
  
  # Fetal5

  load(file = "2021-06-05/fetal kidney5(redefine).Robj")
  Fetal_rename <- scell                                                                                   
  DimPlot(Fetal_rename,pt.size = 1)
  rm(scell)
  
  DimPlot(Fetal_rename)
  DefaultAssay(merge) <- "RNA"
  Idents(merge,cells=WhichCells(merge,idents = c("Endothelial cell_PLVAP high","Endothelial cell_CCL21 high"))) <- "Endothelial cell"
  Idents(merge,cells=WhichCells(merge,idents = c("Interstitial cell_POSTN high","Nephrogenic mesenchyme cell_DAPL1 high"))) <- "Mesangial cell"
  Idents(merge,cells=WhichCells(merge,idents = c("S-shaped body cell","S-shaped body medial cell"))) <- "S-shaped body cell"
  Idents(merge,cells=WhichCells(merge,idents = c("S-shaped body cell_CFAP126 high"))) <- "S-shaped body cell"
  Idents(merge,cells=WhichCells(merge,idents = c("Collecting duct cell_CRABP1 high","Collecting duct cell_CALB1 high"))) <- "Collecting duct cell"
  levels(Fetal_rename)
  Fetal_rename <- subset(Fetal_rename, idents  = c("Collecting duct cell","Loop of Henle progenitor cell","Proximal tubule progenitor cell"))
  
  # KPC organoid
  setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
  wd <- getwd()
  time <- as.character(Sys.Date())
  dir.create(time)
  load(file = "2022-07-24/organoid_monocle3_cds.Robj")
  
  cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.00045)
  plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
  
  load(file = "2022-07-17/KPC.Robj")
  DimPlot(KPC)
  KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
  colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
  KPC@active.ident <- KPC$Monocle.Cluster
  DimPlot(KPC, pt.size = 1)
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"LH-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Cell cycle"
  
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"LH-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Epithelial cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "6"))<-"Cycling cell"
  
  levels(KPC@active.ident)
  KPC <- subset(KPC,idents = c("CD-like","PT-like","LH-like"))
  
  # 整合
  DimPlot(KPC)
  scell.anchors <- FindIntegrationAnchors(object.list = c(KPC,Fetal_rename), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(scell)
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  merge <- scell
  rm(scell,scell.anchors)
  
  
  a <- as.data.frame(Fetal_rename@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents", "cell_id")
  
  b <- as.data.frame(KPC@active.ident)
  b[,2] <- rownames(b)
  colnames(b) <- c("idents","cell_id")
  c <- rbind(b,a)
  
  rm(a,b)
  
  metadata <- as.data.frame(merge@active.ident)
  metadata[,2] <- rownames(metadata)
  colnames(metadata) <- c("orig.ident", "cell_id")
  metadata <- full_join(c,metadata, by=c("cell_id"))
  rownames(metadata) <- metadata$cell_id
  
  metadata$idents <- factor(metadata$idents)
  merge@meta.data[["idents"]] <- metadata$idents
  Idents(merge) <- merge@meta.data$idents
  DimPlot(merge,pt.size = 1)
  
  # 绘图
  DimPlot(merge,pt.size = 1.5,cols = c("#BC8589","#2A7EC6","#00ABA9","#D19D88","#88D1AC","#6FB2FF"))+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  DimPlot(merge,pt.size = 0.7,split.by = "orig.ident")
  DimPlot(merge,pt.size = 0.7, group.by = "orig.ident")
  show_col(hue_pal()(6))
  # 配色gray
  # (#F8766D, #d19d88) 粉棕
  # (#BC8589, #D19D88)棕色系
  # (#2A7EC6, #6FB2FF) 蓝色系
  # (#00ABA9， #88D1AC)绿色系
  # #F8766D #B79F00 #00BA3B #00BFC4 #619CFF #F564E3 
  
  # 拆分三组细胞组后在PS中组合
  LH <- subset(merge,idents = c("Loop of Henle progenitor cell","LH-like")) # "#BC8589","#D19D88"
  PT <- subset(merge,idents = c("Proximal tubule progenitor cell","PT-like")) # "#2A7EC6","#6FB2FF"
  CD <- subset(merge,idents = c("Collecting duct cell","CD-like"))# "#00ABA9", "#88D1AC"
  DimPlot(CD,pt.size = 1.5,cols = c("#00ABA9", "#88D1AC"))+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm")) #该方法不行，除非将分群的名称定为等长，否则会因为长度不一样导致图片比例不同
  
  
  Idents(merge,cells=WhichCells(KPC,idents = "Collecting duct cell"))<-""
  Idents(merge,cells=WhichCells(KPC,idents = "Loop of Henle progenitor cell"))<-""
  Idents(merge,cells=WhichCells(KPC,idents = "Proximal tubule progenitor cell"))<-""
  
  # LH
  DimPlot(merge,pt.size = 1.5,cols = c("#BC8589","grey","grey","#D19D88","grey","grey"))+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  
  # PT
  DimPlot(merge,pt.size = 1.5, cols = c("grey","#2A7EC6","grey","grey","grey","#6FB2FF"))+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  # CD
  DimPlot(merge,pt.size = 1.5, cols = c("grey","grey","#00ABA9","grey","#88D1AC","grey"))+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  
  #读取Kidney数据
  load(file = "2021-06-05/Adult kidney2(redefine).Robj")
  Adult_rename <- scell
  DimPlot(Adult_rename,pt.size = 1)
  rm(scell)
  Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Loop of Henle(Thick ascending limb)","Loop of Henle_SPP1 high","Loop of Henle_SFN high"))) <- "Loop of Henle"
  Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Intercalated cell_SPINK1 high","Intercalated cell_SLC26A4 high"))) <- "Intercalated cell"
  Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("T cell","B cell","Macrophage"))) <- "Immunocytes"
  Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Fenestrated endothelial cell_EMCN high","Fenestrated endothelial cell_SELE high"))) <- "Fenestrated endothelial cell"
  Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Proximal tubule cell_MT1G high","Proximal tubule cell_ALDOB high"))) <- "Proximal tubule"
  Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Intercalated cell"))) <- "Collecting duct" #因为principle cell不是主细胞的意思，且与Intercalated cell相差较远，故移除
  levels(Adult_rename)
  Adult_rename <- subset(Adult_rename, idents  = c("Collecting duct","Loop of Henle","Proximal tubule"))
  
  # KPC organoid
  setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
  wd <- getwd()
  time <- as.character(Sys.Date())
  dir.create(time)
  load(file = "2022-07-24/organoid_monocle3_cds.Robj")
  
  cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.0003)
  plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
  
  load(file = "2022-07-17/KPC.Robj")
  DimPlot(KPC)
  KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
  colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
  KPC@active.ident <- KPC$Monocle.Cluster
  DimPlot(KPC, pt.size = 1)
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"LH-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Cycling cell"
  levels(KPC@active.ident)
  KPC <- subset(KPC,idents = c("CD-like","PT-like","LH-like"))
  
  # 整合
  scell.anchors <- FindIntegrationAnchors(object.list = c(KPC,Adult_rename), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(scell)
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  merge <- scell
  rm(scell,scell.anchors)
  
  
  a <- as.data.frame(Adult_rename@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents", "cell_id")
  
  b <- as.data.frame(KPC@active.ident)
  b[,2] <- rownames(b)
  colnames(b) <- c("idents","cell_id")
  c <- rbind(b,a)
  
  rm(a,b)
  
  metadata <- as.data.frame(merge@active.ident)
  metadata[,2] <- rownames(metadata)
  colnames(metadata) <- c("orig.ident", "cell_id")
  metadata <- full_join(c,metadata, by=c("cell_id"))
  rownames(metadata) <- metadata$cell_id
  
  metadata$idents <- factor(metadata$idents)
  merge@meta.data[["idents"]] <- metadata$idents
  Idents(merge) <- merge@meta.data$idents
  DimPlot(merge,pt.size = 1)
  
  # 绘图
  DimPlot(merge,pt.size = 0.7)
  DimPlot(merge,pt.size = 0.7,split.by = "orig.ident")
  DimPlot(merge,pt.size = 0.7, group.by = "orig.ident")
  show_col(hue_pal()(6))
  # #F8766D #B79F00 #00BA3B #00BFC4 #619CFF #F564E3 
  
  DimPlot(merge,pt.size = 1.2,cols = c("#bc8589","#2a7ec6","#00aba9","#88d1ac","#6fb2ff","#d19d88"))+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  
  # LH
  DimPlot(merge,pt.size = 1.2,cols = c("#bc8589","grey","grey","grey","grey","#d19d88"))+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  
  # PT
  DimPlot(merge,pt.size = 1.2, cols = c("grey","#2a7ec6","grey","grey","#6fb2ff","grey"))+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  # CD
  DimPlot(merge,pt.size = 1.2, cols = c("grey","grey","#00aba9","#88d1ac","grey","grey"))+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  
# 与IPSC数据进行比对（数据来源GSE114802）-----
  setwd("D:/Zuo Lab/Data/iPSC single cell organoid")
  ipsc.data <- Read10X(data.dir = "GSE114802_org")
  ipsc_organoid <- CreateSeuratObject(counts = ipsc.data, project = "ipsc_organoid", min.cells = 3, min.features = 200)
  ipsc_organoid[["percent.mt"]] <- PercentageFeatureSet(ipsc_organoid, pattern = "^MT-")
  ipsc_organoid <- subset(ipsc_organoid, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  VlnPlot(ipsc_organoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  ipsc_organoid <- NormalizeData(ipsc_organoid, normalization.method = "LogNormalize", scale.factor = 10000)
  ipsc_organoid <- FindVariableFeatures(ipsc_organoid, selection.method = "vst", nfeatures = 2000)
  
  ipsc.data <- Read10X(data.dir = "GSE114802_org4")
  ipsc_organoid2 <- CreateSeuratObject(counts = ipsc.data, project = "ipsc_organoid2", min.cells = 3, min.features = 200)
  ipsc_organoid2[["percent.mt"]] <- PercentageFeatureSet(ipsc_organoid2, pattern = "^MT-")
  ipsc_organoid2 <- subset(ipsc_organoid2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  VlnPlot(ipsc_organoid2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  ipsc_organoid2 <- NormalizeData(ipsc_organoid2, normalization.method = "LogNormalize", scale.factor = 10000)
  ipsc_organoid2 <- FindVariableFeatures(ipsc_organoid2, selection.method = "vst", nfeatures = 2000)
  
  scell.anchors <- FindIntegrationAnchors(object.list = c(ipsc_organoid,ipsc_organoid2), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  rm(ipsc.data)
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(scell)
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  FeaturePlot(scell, features = c("SOX9"), pt.size = 1)  
  FeaturePlot(scell, features = c("PODXL"), pt.size = 1)  # Cluster3
  FeaturePlot(scell, features = c("NPHS1"), pt.size = 1)  # Cluster3 podocyte
  
  FeaturePlot(scell, features = c("PECAM1"), pt.size = 1)  # Cluster5 endothelium
  
  FeaturePlot(scell, features = c("GATA3"), pt.size = 1)  # Cluster8 gata3 UE/proximal tubule?
  FeaturePlot(scell, features = c("MEIS1"), pt.size = 1)  # Cluster8

  
  nphs1 <- rownames(scell)[match("NPHS1",scell$Symbol)]
  FeaturePlot(comb, nphs1,cols.use = c("grey","blue"),no.legend = FALSE)
  
  
  # GSE109718数据
  a <- read.table(file = "GSE109718/GSM2949337_75_NLW44_B_NoVegf.txt", sep = "\t",header = T,row.names = 1)
  b <- read.table(file = "GSE109718/GSM2949339_77_NLW44_D_NoVegf.txt", sep = "\t",header = T,row.names = 1)

  ipsc_organoid <- CreateSeuratObject(counts = a, project = "ipsc_organoid", min.cells = 3, min.features = 200)
  ipsc_organoid[["percent.mt"]] <- PercentageFeatureSet(ipsc_organoid, pattern = "^MT-")
  ipsc_organoid <- subset(ipsc_organoid, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  VlnPlot(ipsc_organoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  ipsc_organoid <- NormalizeData(ipsc_organoid, normalization.method = "LogNormalize", scale.factor = 10000)
  ipsc_organoid <- FindVariableFeatures(ipsc_organoid, selection.method = "vst", nfeatures = 2000)
  

  ipsc_organoid2 <- CreateSeuratObject(counts = b, project = "ipsc_organoid2", min.cells = 3, min.features = 200)
  ipsc_organoid2[["percent.mt"]] <- PercentageFeatureSet(ipsc_organoid2, pattern = "^MT-")
  ipsc_organoid2 <- subset(ipsc_organoid2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  VlnPlot(ipsc_organoid2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  ipsc_organoid2 <- NormalizeData(ipsc_organoid2, normalization.method = "LogNormalize", scale.factor = 10000)
  ipsc_organoid2 <- FindVariableFeatures(ipsc_organoid2, selection.method = "vst", nfeatures = 2000)
  
  scell.anchors <- FindIntegrationAnchors(object.list = c(ipsc_organoid,ipsc_organoid2), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  rm(ipsc.data)
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(scell)
  scell <- RunPCA(scell, npcs = 20)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:20)
  scell <- FindClusters(object = scell, resolution = 0.15) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:20)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  FeaturePlot(scell, features = c("WFDC2"), pt.size = 1)  # Cluster0 Early Tubular
  FeaturePlot(scell, features = c("NPHS1"), pt.size = 1)  # Cluster1 Early podocyte
  FeaturePlot(scell, features = c("COL3A1"), pt.size = 1)  # Cluster2,5 stromal
  
  FeaturePlot(scell, features = c("TAGLN"), pt.size = 1)  # Cluster2 stromal
  
  FeaturePlot(scell, features = c("SLC3A1"), pt.size = 1)  # Cluster2 Proximal tubular
  FeaturePlot(scell, features = c("NPHS2"), pt.size = 1)  # Cluster2 stromal
  FeaturePlot(scell, features = c("PODXL"), pt.size = 1)  # Cluster2 stromal
  
  rm(list = ls())
  gc()
  
  # GSE119561
  ipsc.data <- Read10X(data.dir = "GSE119561/E6")
  ipsc_organoid <- CreateSeuratObject(counts = ipsc.data, project = "ipsc_organoid", min.cells = 3, min.features = 200)
  ipsc_organoid[["percent.mt"]] <- PercentageFeatureSet(ipsc_organoid, pattern = "^MT-")
  ipsc_organoid <- subset(ipsc_organoid, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  VlnPlot(ipsc_organoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  ipsc_organoid <- NormalizeData(ipsc_organoid, normalization.method = "LogNormalize", scale.factor = 10000)
  ipsc_organoid <- FindVariableFeatures(ipsc_organoid, selection.method = "vst", nfeatures = 2000)
  
  ipsc.data <- Read10X(data.dir = "GSE119561/SIX2_CRE")
  ipsc_organoid2 <- CreateSeuratObject(counts = ipsc.data, project = "ipsc_organoid2", min.cells = 3, min.features = 200)
  ipsc_organoid2[["percent.mt"]] <- PercentageFeatureSet(ipsc_organoid2, pattern = "^MT-")
  ipsc_organoid2 <- subset(ipsc_organoid2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  VlnPlot(ipsc_organoid2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  ipsc_organoid2 <- NormalizeData(ipsc_organoid2, normalization.method = "LogNormalize", scale.factor = 10000)
  ipsc_organoid2 <- FindVariableFeatures(ipsc_organoid2, selection.method = "vst", nfeatures = 2000)
  
  scell.anchors <- FindIntegrationAnchors(object.list = c(ipsc_organoid,ipsc_organoid2), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  rm(ipsc.data)
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(scell)
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  FeaturePlot(scell, features = c("NPHS1"), pt.size = 1)  # Cluster9 Podocytes
  
  FeaturePlot(scell, features = c("CUBN"), pt.size = 1)  # Cluster9 Proximal Tubule
  FeaturePlot(scell, features = c("MT1G"), pt.size = 1)  # Cluster9 Proximal Tubule
  
  FeaturePlot(scell, features = c("SLC12A1"), pt.size = 1)  # Cluster9 Proximal Tubule
  FeaturePlot(scell, features = c("POU3F3"), pt.size = 1)  # Cluster9 Proximal Tubule
  
  FeaturePlot(scell, features = c("TMEM213"), pt.size = 1)  # Cluster9 Proximal Tubule
  FeaturePlot(scell, features = c("GATA3"), pt.size = 1)  # Cluster9 Proximal Tubule
  
  FeaturePlot(scell, features = c("CCND1"), pt.size = 1)  # Cluster9 Proximal Tubule
  
  FeaturePlot(scell, features = c("SIX1"), pt.size = 1)  # Cluster9 Proximal Tubule
  FeaturePlot(scell, features = c("CRABP1"), pt.size = 1)  # Cluster9 Proximal Tubule
  
  
  sub <- subset(scell, idents  = c(9))
  DimPlot(sub,pt.size = 1)

  DefaultAssay(sub) <- "RNA"
  sub[['percent.mito']] <- PercentageFeatureSet(sub, pattern = "^MT-")
  sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000)
  sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000)
  sub <- ScaleData(sub,  vars.to.regress = c("nCount_RNA", "percent.mito"))
  #Integrated
  DefaultAssay(sub) <- "integrated"
  sub <- ScaleData(sub, vars.to.regress = c("nCount_RNA", "percent.mito"))
  sub <- RunPCA(sub, npcs = 50)
  sub <- ProjectDim(object = sub)
  
  sub <- FindNeighbors(object = sub, dims = 1:30)
  sub <- FindClusters(object = sub, resolution = 0.8) 
  sub <- RunUMAP(sub, reduction = "pca", dims = 1:30)
  DimPlot(sub, pt.size = 1)
  
  
  FeaturePlot(sub, features = c("CUBN"), pt.size = 1)  # Cluster1,2,4 Proximal Tubule
  FeaturePlot(sub, features = c("SLC12A1"), pt.size = 1)  # Cluster0,5 Distal Tubule & Loop of Henle
  FeaturePlot(sub, features = c("GATA3"), pt.size = 1)  # Cluster6 Collecting Duct
  
  FeaturePlot(sub, features = c("CCND1"), pt.size = 1)  # Cluster3,7 Epithelium(cycling)
  
  sub_rename <- RenameIdents(sub, `0` = "Distal Tubule & Loop of Henle", `1` = "Proximal Tubule", `2` = "Proximal Tubule", 
                        `3` = "Epithelium(cycling)", `4` = "Proximal Tubule", `5` = "Distal Tubule & Loop of Henle", `6` = "Collecting Duct", `7` = "Epithelium(cycling)")
  DimPlot(sub_rename, pt.size = 1)
  
  
  
# 单独分析-----
  setwd(dir = "D:/Zuo Lab/Data/iPSC single cell organoid")
  ipsc.data <- Read10X(data.dir = "GSE119561/E6")
  ipsc.data <- Read10X(data.dir = "GSE119561/SIX2_CRE")
  
  
  ipsc_organoid <- CreateSeuratObject(counts = ipsc.data, project = "ipsc_organoid", min.cells = 3, min.features = 200)
  ipsc_organoid[["percent.mt"]] <- PercentageFeatureSet(ipsc_organoid, pattern = "^MT-")
  ipsc_organoid <- subset(ipsc_organoid, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  VlnPlot(ipsc_organoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  ipsc_organoid <- NormalizeData(ipsc_organoid, normalization.method = "LogNormalize", scale.factor = 10000)
  ipsc_organoid <- FindVariableFeatures(ipsc_organoid, selection.method = "vst", nfeatures = 2000)

  ipsc_organoid <- ScaleData(ipsc_organoid, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ipsc_organoid <- RunPCA(ipsc_organoid, npcs = 50)
  ipsc_organoid <- ProjectDim(object = ipsc_organoid)
  
  ipsc_organoid <- FindNeighbors(object = ipsc_organoid, dims = 1:30)
  ipsc_organoid <- FindClusters(object = ipsc_organoid, resolution = 0.5) 
  
  ipsc_organoid <- RunUMAP(ipsc_organoid, reduction = "pca", dims = 1:30)
  DimPlot(ipsc_organoid,pt.size = 1)
  DimPlot(ipsc_organoid,group.by = "orig.ident",pt.size = 1)
  DimPlot(ipsc_organoid,split.by = "orig.ident",pt.size = 1,ncol = 3)
  # E6
  FeaturePlot(ipsc_organoid, features = c("NPHS1"), pt.size = 1)  # Cluster9 Podocytes
  FeaturePlot(ipsc_organoid, features = c("MAFB"), pt.size = 1)
  
  FeaturePlot(ipsc_organoid, features = c("CUBN"), pt.size = 1)  # Cluster4 Pt-Tubule cells
  FeaturePlot(ipsc_organoid, features = c("MT1G"), pt.size = 1)  
  FeaturePlot(ipsc_organoid, features = c("SLC12A1"), pt.size = 1) # DTLH
  FeaturePlot(ipsc_organoid, features = c("POU3F3"), pt.size = 1)  
  FeaturePlot(ipsc_organoid, features = c("TMEM213"), pt.size = 1)  # Collecting duct
  FeaturePlot(ipsc_organoid, features = c("GATA3"), pt.size = 1)
  
  FeaturePlot(ipsc_organoid, features = c("MYOD1"), pt.size = 1)  # Cluster8,12 Muscle
  FeaturePlot(ipsc_organoid, features = c("MYOG"), pt.size = 1)
  
  FeaturePlot(ipsc_organoid, features = c("SOX2"), pt.size = 1)  # Cluster10 Progenitor cells
  
  FeaturePlot(ipsc_organoid, features = c("CRABP1"), pt.size = 1)  # Cluster0,2 Stroma
  FeaturePlot(ipsc_organoid, features = c("CRABP2"), pt.size = 1) 
  FeaturePlot(ipsc_organoid, features = c("SIX1"), pt.size = 1) 
  
  FeaturePlot(ipsc_organoid, features = c("PTN"), pt.size = 1)  # Cluster9 Neural
  FeaturePlot(ipsc_organoid, features = c("FABP7"), pt.size = 1)  
  
  FeaturePlot(ipsc_organoid, features = c("CCND1"), pt.size = 1)  # cycling
  FeaturePlot(ipsc_organoid, features = c("EPCAM"), pt.size = 1)  
  
  FeaturePlot(ipsc_organoid, features = c("CCNB1"), pt.size = 1)  # cycling
  FeaturePlot(ipsc_organoid, features = c("CCNB1"), pt.size = 1) 
  
  FeaturePlot(ipsc_organoid, features = c("LYPD1"), pt.size = 1)  # # Cluster13 Early Nephron
  FeaturePlot(ipsc_organoid, features = c("SOX2"), pt.size = 1) 
  
  sub <- subset(ipsc_organoid, idents  = c(2,3,4,5,6,7,8,9,10,11,12,13,14))
  DimPlot(sub,pt.size = 1)
  
  
  sub[['percent.mito']] <- PercentageFeatureSet(sub, pattern = "^MT-")
  sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000)
  sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000)
  sub <- ScaleData(sub,  vars.to.regress = c("nCount_RNA", "percent.mito"))
  sub <- RunPCA(sub, npcs = 50)
  sub <- ProjectDim(object = sub)
  
  sub <- FindNeighbors(object = sub, dims = 1:30)
  sub <- FindClusters(object = sub, resolution = 0.6) 
  sub <- RunUMAP(sub, reduction = "pca", dims = 1:30)
  DimPlot(sub, pt.size = 1)
  
  FeaturePlot(sub, features = c("NPHS1"), pt.size = 1)  # Cluster3, 7 Podocytes
  FeaturePlot(sub, features = c("MAFB"), pt.size = 1)
  
  FeaturePlot(sub, features = c("CUBN"), pt.size = 1)  # Cluster9 Pt-Tubule cells
  FeaturePlot(sub, features = c("MT1G"), pt.size = 1)  
  FeaturePlot(sub, features = c("SLC12A1"), pt.size = 1) # DTLH
  FeaturePlot(sub, features = c("POU3F3"), pt.size = 1)  
  FeaturePlot(sub, features = c("TMEM213"), pt.size = 1)  # Collecting duct
  FeaturePlot(sub, features = c("GATA3"), pt.size = 1)
  
  FeaturePlot(sub, features = c("MYOD1"), pt.size = 1)  # Cluster8,11 Muscle
  FeaturePlot(sub, features = c("MYOG"), pt.size = 1)
  
  FeaturePlot(sub, features = c("SOX2"), pt.size = 1)  # Cluster10 Progenitor cells
  
  FeaturePlot(sub, features = c("CRABP1"), pt.size = 1)  # Cluster0,6,13 Stroma
  FeaturePlot(sub, features = c("CRABP2"), pt.size = 1) 
  FeaturePlot(sub, features = c("SIX1"), pt.size = 1) 
  
  FeaturePlot(sub, features = c("PTN"), pt.size = 1)  # Cluster2,5,14 Neural
  FeaturePlot(sub, features = c("FABP7"), pt.size = 1)  
  
  FeaturePlot(sub, features = c("CCND1"), pt.size = 1)  # cycling
  FeaturePlot(sub, features = c("EPCAM"), pt.size = 1)  
  
  FeaturePlot(sub, features = c("LYPD1"), pt.size = 1)  # # Cluster12 Early Nephron
  
  
  
  sub_rename <- RenameIdents(sub, `0` = "iPSC_Stroma", `1` = "iPSC_Undifferentiated cell", `2` = "iPSC_Neural", 
                             `3` = "iPSC_Podocytes", `4` = "iPSC_Epithelium(cycling)", `5` = "iPSC_Neural", `6` = "iPSC_Stroma", `7` = "iPSC_Podocytes",
                             `8` = "iPSC_Muscle", `9` = "iPSC_Proximal Tubule",`10` = "iPSC_Nephron Progenitors",
                             `11` = "iPSC_Muscle",`12` = "iPSC_Early Nephron",`13` = "iPSC_Stroma", `14` = "iPSC_Neural")
  
  Idents(sub_rename, cells=WhichCells(sub_tub_rename,idents = "Epithelium(cycling)"))<-"iPSC_Epithelium(cycling)"
  Idents(sub_rename, cells=WhichCells(sub_tub_rename,idents = "Distal Tubule & Loop of Henle"))<-"iPSC_Distal Tubule & Loop of Henle"
  Idents(sub_rename, cells=WhichCells(sub_tub_rename,idents = "Proximal Tubule"))<-"iPSC_Proximal Tubule"
  Idents(sub_rename, cells=WhichCells(sub_tub_rename,idents = "Collecting Duct"))<-"iPSC_Collecting Duct"
  DimPlot(sub_rename, pt.size = 1)
  DimPlot(sub_rename, pt.size = 1 )+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  
  scell.markers <- FindAllMarkers(sub_rename, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  write.table(scell.markers, file = paste(wd,"/",time,"/","DEG.csv", sep = ""), sep=",")
  
  FeaturePlot(sub, features = c("MAFB"), pt.size = 1)  # cycling
  FeaturePlot(sub, features = c("KIAA0101"), pt.size = 1) 
  

  FeaturePlot(sub, features = c("SOX2"), pt.size = 1) 

  
  # 分离E6中的肾前体细胞
  FeaturePlot(ipsc_organoid, features = c("DAPL1"), pt.size = 1)
  FeaturePlot(ipsc_organoid, features = c("LYPD1"), pt.size = 1)
  FeaturePlot(ipsc_organoid, features = c("CITED1"), pt.size = 1)
  FeaturePlot(ipsc_organoid, features = c("LHX1"), pt.size = 1)
  FeaturePlot(ipsc_organoid, features = c("JAG1"), pt.size = 1) # Cluster9,13 Nephron progenitors/Early Nephron
  sub_tub <- subset(ipsc_organoid, idents  = c(4,9,13))
  
  sub_tub <- subset(ipsc_organoid, idents  = c(4))
  DimPlot(sub_tub,pt.size = 1)
  

  sub_tub[['percent.mito']] <- PercentageFeatureSet(sub_tub, pattern = "^MT-")
  sub_tub <- NormalizeData(sub_tub, normalization.method = "LogNormalize", scale.factor = 10000)
  sub_tub <- FindVariableFeatures(sub_tub, selection.method = "vst", nfeatures = 2000)
  sub_tub <- ScaleData(sub_tub,  vars.to.regress = c("nCount_RNA", "percent.mito"))
  sub_tub <- RunPCA(sub_tub, npcs = 50)
  sub_tub <- ProjectDim(object = sub_tub)
  
  sub_tub <- FindNeighbors(object = sub_tub, dims = 1:30)
  sub_tub <- FindClusters(object = sub_tub, resolution = 0.6) 
  sub_tub <- RunUMAP(sub_tub, reduction = "pca", dims = 1:30)
  DimPlot(sub_tub, pt.size = 1)
  
  
  FeaturePlot(sub, features = c("CUBN"), pt.size = 1)  # Cluster2,4 Proximal Tubule
  FeaturePlot(sub, features = c("SLC12A1"), pt.size = 1)  # Cluster1 Distal Tubule & Loop of Henle
  FeaturePlot(sub, features = c("GATA3"), pt.size = 1)  # Cluster8 Collecting Duct
  
  FeaturePlot(sub, features = c("CCND1"), pt.size = 1)  # Cluster0,3,5,6,7 Epithelium(cycling)
  
  sub_tub_rename <- RenameIdents(sub_tub, `0` = "Epithelium(cycling)", `1` = "Distal Tubule & Loop of Henle", `2` = "Proximal Tubule", 
                             `3` = "Epithelium(cycling)", `4` = "Proximal Tubule", `5` = "Epithelium(cycling)", `6` = "Epithelium(cycling)", `7` = "Epithelium(cycling)",
                             `8` = "Collecting Duct")
  DimPlot(sub_tub_rename, pt.size = 1)

    
# 与iPSC数据整合-----
  scell.anchors <- FindIntegrationAnchors(object.list = c(sub_rename,Fetal_rename), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1)
  
  merge <- scell
  rm(scell,scell.anchors)

  
  a <- as.data.frame(sub_rename@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents", "cell_id")
  
  b <- as.data.frame(Fetal_rename@active.ident)
  b[,2] <- rownames(b)
  colnames(b) <- c("idents","cell_id")
  c <- rbind(a,b)
  
  rm(a,b)
  
  metadata <- as.data.frame(merge@active.ident)
  metadata[,2] <- rownames(metadata)
  colnames(metadata) <- c("orig.ident", "cell_id")
  metadata <- full_join(c,metadata, by=c("cell_id"))
  rownames(metadata) <- metadata$cell_id
  
  metadata$idents <- factor(metadata$idents)
  merge@meta.data[["idents"]] <- metadata$idents # 这里的顺序不一定是对的
  Idents(merge) <- merge@meta.data$idents  
  
  DefaultAssay(merge) <- "RNA"
  DimPlot(merge,pt.size = 0.7,group.by = "seurat_clusters")
  DimPlot(merge,pt.size = 1)
  DimPlot(merge,pt.size = 1,split.by = "orig.ident")
  DimPlot(merge,pt.size = 1, group.by = "orig.ident")
  show_col(hue_pal()(20))
  # LH
  DimPlot(merge,pt.size = 0.7,cols = c("grey","#C49200","grey","grey","#00B6EB",
                                       "grey","grey"))
  # PT
  DimPlot(merge,pt.size = 0.7,cols = c("grey","grey","#53B400","grey","grey",
                                       "grey","#FB61D7"))
  # CD
  DimPlot(merge,pt.size = 0.7,cols = c("grey","grey","grey","#00C094","grey",
                                       "#A58AFF","grey"))

  
  setwd("D:/Zuo Lab/Data/iPSC single cell organoid/ROBJ/")
  save(Fetal_rename, file = "Fetal5_tubule-only.Robj")
  save(sub_rename, file = "ipsc_tubule-only.Robj")
  save(KPC, file = "KPC_tubule-only.Robj")
  save(merge, file = "merge_ipsc_fetal.Robj")
  
  
# iPSC和KPC数据整合 ------
  scell.anchors <- FindIntegrationAnchors(object.list = c(KPC, sub_rename), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1)
  
  merge <- scell
  rm(scell,scell.anchors)
  
  
  a <- as.data.frame(sub_rename@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents", "cell_id")
  
  b <- as.data.frame(KPC@active.ident)
  b[,2] <- rownames(b)
  colnames(b) <- c("idents","cell_id")
  c <- rbind(b,a)
  
  rm(a,b)
  
  metadata <- as.data.frame(merge@active.ident)
  metadata[,2] <- rownames(metadata)
  colnames(metadata) <- c("orig.ident", "cell_id")
  metadata <- full_join(c,metadata, by=c("cell_id"))
  rownames(metadata) <- metadata$cell_id
  
  metadata$idents <- factor(metadata$idents)
  merge@meta.data[["idents"]] <- metadata$idents # 这里的顺序不一定是对的
  Idents(merge) <- merge@meta.data$idents  
  
  DefaultAssay(merge) <- "RNA"
  DimPlot(merge,pt.size = 0.7,group.by = "seurat_clusters")
  DimPlot(merge,pt.size = 1)
  DimPlot(merge,pt.size = 1,split.by = "orig.ident")
  DimPlot(merge,pt.size = 1, group.by = "orig.ident")
  show_col(hue_pal()(7))
  # LH
  DimPlot(merge,pt.size = 0.7,cols = c("#F8766D","grey","grey","grey","#00B6EB",
                                       "grey","grey"))
  # PT
  DimPlot(merge,pt.size = 0.7,cols = c("grey","#C49A00","grey","grey","grey",
                                       "#A58AFF","grey"))
  # CD
  DimPlot(merge,pt.size = 0.7,cols = c("grey","grey","#53B400","grey","grey",
                                       "grey","#FB61D7"))
  
# iPSC的拟时序-----
  # 创建CDS文件
  data <- as(as.matrix(sub_rename@assays$RNA@counts), 'sparseMatrix')
  pd <-  sub_rename@meta.data
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  colnames(pd)
  
  cds <- new_cell_data_set(data,
                           cell_metadata  = pd,
                           gene_metadata  = fData)
  #Pre-process the data
  cds = preprocess_cds(cds, method = c("PCA"),num_dim = 100)
  plot_pc_variance_explained(cds)
  
  #Reduce dimensionality and visualize the cells
  cds = reduce_dimension(cds) #Monocle uses UMAP by default 
  plot_cells(cds) 

  cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.0003)  
  plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
  
  cds@int_colData$reducedDims$UMAP <- sub_rename@reductions$umap@cell.embeddings
  cds@clusters$UMAP$clusters <- sub_rename@active.ident

  DimPlot(sub, pt.size = 1 )
  
  
  # 定义细胞群（根据之前分析的小提琴图）
  #Annotate your cells according to type

  plot_cells(cds, group_cells_by="cluster", color_cells_by="assigned_cell_type",group_label_size=4,cell_size=1)
  plot_cells(cds,group_cells_by="cluster", color_cells_by="assigned_cell_type",
             graph_label_size=2,cell_size=2,label_cell_groups=FALSE,show_trajectory_graph = F, )+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  plot_cells(cds, color_cells_by = "cluster", cell_size=1, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE,show_trajectory_graph = F)+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.5))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.5))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(2)))+
    theme(axis.title.y = element_text(color="black", size=rel(2)))+
    theme(axis.text.x = element_text(color="black", size=rel(2)))+
    theme(axis.text.y = element_text(color="black", size=rel(2)))+
    theme(legend.text=element_text(color="black",size=rel(2)))+
    theme(legend.title=element_text(color="black",size=rel(2)))+
    theme(legend.key.size = unit(0.5, "cm"))
  colnames(colData(cds))[8] <- "Cell type"
  
  save(cds,file = paste(wd,"/",time,"/","organoid_monocle3_cds.Robj", sep = ""))
  
  
  # 轨迹学习Learn the trajectory graph（使用learn_graph()函数）
  cds <- learn_graph(cds,use_partition = F)
  plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE)
  
  # 2. 细胞按拟时排序
  cds <- order_cells(cds)
  plot_cells(cds, color_cells_by = "pseudotime", cell_size=1, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
  
  plot_cells(cds, color_cells_by = "pseudotime", cell_size=1, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.5))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.5))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(2)))+
    theme(axis.title.y = element_text(color="black", size=rel(2)))+
    theme(axis.text.x = element_text(color="black", size=rel(2)))+
    theme(axis.text.y = element_text(color="black", size=rel(2)))+
    theme(legend.text=element_text(color="black",size=rel(2)))+
    theme(legend.title=element_text(color="black",size=rel(2)))+
    theme(legend.key.size = unit(0.5, "cm"))
  
# iPSC和KPC数据整合(20221211) ------
  load(file = "2022-07-17/KPC.Robj")
  DimPlot(KPC)
  KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
  colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
  KPC@active.ident <- KPC$Monocle.Cluster
  DimPlot(KPC, pt.size = 1)
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"KPC_CD-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"KPC_PT-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"KPC_LH-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"KPC_Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"KPC_Epithelial cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "6"))<-"KPC_Cell cycle"
  
  
  scell.anchors <- FindIntegrationAnchors(object.list = c(sub_rename, KPC), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1)
  
  merge <- scell
  rm(scell,scell.anchors)
  
  
  a <- as.data.frame(sub_rename@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents", "cell_id")
  
  b <- as.data.frame(KPC@active.ident)
  b[,2] <- rownames(b)
  colnames(b) <- c("idents","cell_id")
  c <- rbind(a,b)
  
  rm(a,b)
  
  metadata <- as.data.frame(merge@active.ident)
  metadata[,2] <- substr(rownames(metadata),1,18)
  colnames(metadata) <- c("orig.ident", "cell_id")
  

  metadata[13607,2] <- "TGTTCCGCAAGCGAGT-1_1"
  c[13607,2] <- "TGTTCCGCAAGCGAGT-1_1"
  
  metadata <- full_join(c,metadata, by=c("cell_id"))
  rownames(metadata) <- metadata$cell_id
  
  metadata$idents <- factor(metadata$idents)
  merge@meta.data[["idents"]] <- metadata$idents # 这里的顺序不一定是对的
  Idents(merge) <- merge@meta.data$idents  
  
  DefaultAssay(merge) <- "RNA"
  DimPlot(merge,pt.size = 0.7,group.by = "seurat_clusters")
  DimPlot(merge,pt.size = 1)
  DimPlot(merge,pt.size = 1,split.by = "orig.ident")
  DimPlot(merge,pt.size = 1, group.by = "orig.ident")
  DimPlot(merge, pt.size = 0.7 )+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  
  show_col(hue_pal()(7))
  # LH
  DimPlot(merge,pt.size = 0.7,cols = c("#F8766D","grey","grey","grey","#00B6EB",
                                       "grey","grey"))
  # PT
  DimPlot(merge,pt.size = 0.7,cols = c("grey","#C49A00","grey","grey","grey",
                                       "#A58AFF","grey"))
  # CD
  DimPlot(merge,pt.size = 0.7,cols = c("grey","grey","#53B400","grey","grey",
                                       "grey","#FB61D7"))
  
# 20230308 以下代码是为了毕业论文以及补充文章跑的代码，将上述的代码重新进行规整-----
  # KPC类器官单细胞数据集以及
  {
    # KPC organoid
    setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
    wd <- getwd()
    time <- as.character(Sys.Date())
    dir.create(time)
    load(file = "2022-07-24/organoid_monocle3_cds.Robj")
    
    cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.00045)
    plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
    
    load(file = "2022-07-17/KPC.Robj")
    DimPlot(KPC)
    KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
    colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
    KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
    KPC@active.ident <- KPC$Monocle.Cluster
    DimPlot(KPC, pt.size = 1)

    Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
    Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
    Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"LH-like"
    Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
    Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Epithelial cell"
    Idents(KPC,cells=WhichCells(KPC,idents = "6"))<-"Cycling cell"
    
    levels(KPC@active.ident)
    KPC_Tub <- subset(KPC,idents = c("CD-like","PT-like","LH-like"))
    
    FeaturePlot(KPC,features = "JAG1")
    colo <- RColorBrewer::brewer.pal(10,"RdYlBu")[1:5]
    colo2 <- RColorBrewer::brewer.pal(10,"RdYlBu")[6:10]
    n <- "POU3F3"
    FeaturePlot(KPC,features = n,sort.cell = T,pt.size =1.5,min.cutoff = "q1",max.cutoff = "q95")+ggtitle(n)+
      scale_colour_gradient(low = rev(colo2),high = rev(colo))+ 
      theme(plot.title = element_text(hjust = 0.5,vjust =-0.6,size=rel(2.7)))+
      theme(axis.line = element_line(size=1, colour = "black"))+
      theme(legend.text=element_text(color="black",size=rel(1.3)))+
      theme(axis.title.x = element_text(color="black", size=rel(1.6)))+
      theme(axis.title.y = element_text(color="black", size=rel(1.6)))
    VlnPlot(object = KPC,features = n)
  }
  { # 20230423 提取出除Cycling cell的分群，进行DEG分析，用于下一步与iPSC的肾细胞分群的GO进行对比看差异
    # KPC organoid
    setwd("D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC")
    wd <- getwd()
    time <- as.character(Sys.Date())
    dir.create(time)
    load(file = "2022-07-24/organoid_monocle3_cds.Robj")
    
    cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.00045)
    plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
    
    load(file = "2022-07-17/KPC.Robj")
    DimPlot(KPC)
    KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
    colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
    KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
    KPC@active.ident <- KPC$Monocle.Cluster
    DimPlot(KPC, pt.size = 1)
    
    Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
    Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
    Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"LH-like"
    Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
    Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Epithelial cell"
    Idents(KPC,cells=WhichCells(KPC,idents = "6"))<-"Cycling cell"
    
    levels(KPC@active.ident)
    KPC_nephron <- subset(KPC,idents = c("CD-like","PT-like","LH-like","Progenitor cell"))
    
    FeaturePlot(KPC,features = "JAG1")
    colo <- RColorBrewer::brewer.pal(10,"RdYlBu")[1:5]
    colo2 <- RColorBrewer::brewer.pal(10,"RdYlBu")[6:10]
    n <- "POU3F3"
    FeaturePlot(KPC,features = n,sort.cell = T,pt.size =1.5,min.cutoff = "q1",max.cutoff = "q95")+ggtitle(n)+
      scale_colour_gradient(low = rev(colo2),high = rev(colo))+ 
      theme(plot.title = element_text(hjust = 0.5,vjust =-0.6,size=rel(2.7)))+
      theme(axis.line = element_line(size=1, colour = "black"))+
      theme(legend.text=element_text(color="black",size=rel(1.3)))+
      theme(axis.title.x = element_text(color="black", size=rel(1.6)))+
      theme(axis.title.y = element_text(color="black", size=rel(1.6)))
    VlnPlot(object = KPC,features = n)
  }
  # iPSC数据集只将肾前体细胞和小管细胞分离出来进行GO分析
  {#iPSC原始定义的数据在下面，此处只做提取后差异分析
    iPSC_nephron <- subset(sub_rename,idents = c("iPSC_Podocytes","iPSC_Proximal Tubule","iPSC_Nephron Progenitors","iPSC_Early Nephron",
                                                 "iPSC_Distal Tubule & Loop of Henle",
                                                 "iPSC_Collecting Duct"))
   DimPlot(iPSC_nephron)
   
   #计算差异基因
   DEG_iPSC_nephron <- FindAllMarkers(iPSC_nephron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
   write.table(DEG_iPSC_nephron, file = paste(wd,"/",time,"/","DEG_iPSC_nephron.csv", sep = ""), sep=",")
   #分析iPSC_nephron GO富集分析结果
   x <- "iPSC_Podocytes"
  for (x in c("iPSC_Podocytes","iPSC_Proximal Tubule","iPSC_Nephron Progenitors","iPSC_Early Nephron",
             "iPSC_Distal Tubule & Loop of Henle",
             "iPSC_Collecting Duct")){
   DEG <- DEG_iPSC_nephron[which(DEG_iPSC_nephron$cluster %in% x),]
   GO <- enrichGO(DEG$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05,keyType = "SYMBOL")
   write.table(GO, file = paste(wd,"/",time,"/",x,".csv", sep = ""), sep=",")}
   
   # 将active.idents作为metadata导入
   iPSC_nephron[["Cluster"]] <- iPSC_nephron@active.ident
   # 绘制热图
   colo <- RColorBrewer::brewer.pal(10,"RdYlBu")[1:5]
   colo2 <- RColorBrewer::brewer.pal(10,"RdYlBu")[6:10]
   top30 <- DEG_iPSC_nephron %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
   DoHeatmap(iPSC_nephron, features = top30$gene)
   DoHeatmap(iPSC_nephron, disp.min = -2,disp.max = 2.5, #调整区间使得颜色更明显
             features = as.character(unique(top30$gene)),   
             group.by = "Cluster",  
             assay = "RNA",  
             group.colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#AB82FF","#90EE90","#00CD00","#008B8B"))+  #设置组别颜色
     scale_fill_gradientn(colors = c("black", "#34366b","#ffe889"))+
     theme(legend.text=element_text(color="black",size=rel(1.5)), # 图例字体大小及颜色
           legend.title=element_text(color="black",size=rel(1.5)), # 图例标题大小及颜色
           legend.key.heigh = unit(0.5, "cm"), # 图例图示部分的高度
           legend.key.width = unit(0.7, "cm")) # 图例图示部分的宽度
  
  }
  # 胎儿肾脏的重新定义单细胞数据集
  {
    # Fetal4
    load(file = "2021-06-05/fetal kidney4(redefine).Robj")
    Fetal_rename <- scell                                                                                   
    DimPlot(Fetal_rename,pt.size = 1)
    rm(scell)
    
    DimPlot(Fetal_rename)
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Endothelial cell_PLVAP high","Endothelial cell_GJA5 high"))) <- "Endothelial cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Interstitial cell_PTN high"))) <- "Mesangial cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Interstitial cell_POSTN high","Nephrogenic mesenchyme cell_DAPL1 high"))) <- "Mesangial cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("S-shaped body cell_LINC01158 high","S-shaped body medial cell"))) <- "S-shaped body cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("S-shaped body cell_CFAP126 high"))) <- "S-shaped body cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Collecting duct cell_CRABP1 high","Collecting duct cell_CALB1 high"))) <- "Collecting duct cell"
    levels(Fetal_rename))
  }
  # 胎儿肾脏单独将肾小管细胞分离出来的数据集
  { # Fetal4
    load(file = "2021-06-05/fetal kidney4(redefine).Robj")
    Fetal_rename <- scell                                                                                   
    DimPlot(Fetal_rename,pt.size = 1)
    rm(scell)
    
    DimPlot(Fetal_rename)
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Endothelial cell_PLVAP high","Endothelial cell_GJA5 high"))) <- "Endothelial cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Interstitial cell_PTN high"))) <- "Mesangial cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Interstitial cell_POSTN high","Nephrogenic mesenchyme cell_DAPL1 high"))) <- "Mesangial cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("S-shaped body cell_LINC01158 high","S-shaped body medial cell"))) <- "S-shaped body cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("S-shaped body cell_CFAP126 high"))) <- "S-shaped body cell"
    Idents(Fetal_rename,cells=WhichCells(Fetal_rename,idents = c("Collecting duct cell_CRABP1 high","Collecting duct cell_CALB1 high"))) <- "Collecting duct cell"
    levels(Fetal_rename)
    Fetal_rename <- subset(Fetal_rename, idents  = c("Collecting duct cell","Loop of Henle progenitor cell","Proximal tubule progenitor cell"))}
  # 成人肾脏的重新定义单细胞数据集(同时将小管分离出来)
  {
    load(file = "2021-06-05/Adult kidney2(redefine).Robj")
    Adult_rename <- scell
    DimPlot(Adult_rename,pt.size = 1)
    rm(scell)
    Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Loop of Henle(Thick ascending limb)","Loop of Henle_SPP1 high","Loop of Henle_SFN high"))) <- "Loop of Henle"
    Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Intercalated cell_SPINK1 high","Intercalated cell_SLC26A4 high"))) <- "Intercalated cell"
    Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("T cell","B cell","Macrophage"))) <- "Immunocytes"
    Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Fenestrated endothelial cell_EMCN high","Fenestrated endothelial cell_SELE high"))) <- "Fenestrated endothelial cell"
    Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Proximal tubule cell_MT1G high","Proximal tubule cell_ALDOB high"))) <- "Proximal tubule"
    Idents(Adult_rename,cells=WhichCells(Adult_rename,idents = c("Intercalated cell"))) <- "Collecting duct" #因为principle cell不是主细胞的意思，且与Intercalated cell相差较远，故移除
    levels(Adult_rename)
    Adult_rename <- subset(Adult_rename, idents  = c("Collecting duct","Loop of Henle","Proximal tubule"))
  }
  # iPSC肾类器官重新定义的单细胞数据集（使用的是E6的数据）
  { # 原始数据分析
    setwd(dir = "D:/Zuo Lab/Data/iPSC single cell organoid")
    ipsc.data <- Read10X(data.dir = "GSE119561/E6")

    ipsc_organoid <- CreateSeuratObject(counts = ipsc.data, project = "ipsc_organoid", min.cells = 3, min.features = 200)
    ipsc_organoid[["percent.mt"]] <- PercentageFeatureSet(ipsc_organoid, pattern = "^MT-")
    ipsc_organoid <- subset(ipsc_organoid, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
    VlnPlot(ipsc_organoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    ipsc_organoid <- NormalizeData(ipsc_organoid, normalization.method = "LogNormalize", scale.factor = 10000)
    ipsc_organoid <- FindVariableFeatures(ipsc_organoid, selection.method = "vst", nfeatures = 2000)
    
    ipsc_organoid <- ScaleData(ipsc_organoid, vars.to.regress = c("nCount_RNA", "percent.mito"))
    ipsc_organoid <- RunPCA(ipsc_organoid, npcs = 50)
    ipsc_organoid <- ProjectDim(object = ipsc_organoid)
    
    ipsc_organoid <- FindNeighbors(object = ipsc_organoid, dims = 1:30)
    ipsc_organoid <- FindClusters(object = ipsc_organoid, resolution = 0.5) 
    
    ipsc_organoid <- RunUMAP(ipsc_organoid, reduction = "pca", dims = 1:30)
    DimPlot(ipsc_organoid,pt.size = 1)
    DimPlot(ipsc_organoid,group.by = "orig.ident",pt.size = 1)
    DimPlot(ipsc_organoid,split.by = "orig.ident",pt.size = 1,ncol = 3)
    
    sub <- subset(ipsc_organoid, idents  = c(2,3,4,5,6,7,8,9,10,11,12,13,14))
    DimPlot(sub,pt.size = 1)
    
    
    sub[['percent.mito']] <- PercentageFeatureSet(sub, pattern = "^MT-")
    sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000)
    sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000)
    sub <- ScaleData(sub,  vars.to.regress = c("nCount_RNA", "percent.mito"))
    sub <- RunPCA(sub, npcs = 50)
    sub <- ProjectDim(object = sub)
    
    sub <- FindNeighbors(object = sub, dims = 1:30)
    sub <- FindClusters(object = sub, resolution = 0.6) 
    sub <- RunUMAP(sub, reduction = "pca", dims = 1:30)
    DimPlot(sub, pt.size = 1)
    
    sub_rename <- RenameIdents(sub, `0` = "iPSC_Stroma", `1` = "iPSC_Undifferentiated cell", `2` = "iPSC_Neural", 
                               `3` = "iPSC_Podocytes", `4` = "iPSC_Epithelium(cycling)", `5` = "iPSC_Neural", `6` = "iPSC_Stroma", `7` = "iPSC_Podocytes",
                               `8` = "iPSC_Muscle", `9` = "iPSC_Proximal Tubule",`10` = "iPSC_Nephron Progenitors",
                               `11` = "iPSC_Muscle",`12` = "iPSC_Early Nephron",`13` = "iPSC_Stroma", `14` = "iPSC_Neural")
    
    Idents(sub_rename, cells=WhichCells(sub_tub_rename,idents = "Epithelium(cycling)"))<-"iPSC_Epithelium(cycling)"
    Idents(sub_rename, cells=WhichCells(sub_tub_rename,idents = "Distal Tubule & Loop of Henle"))<-"iPSC_Distal Tubule & Loop of Henle"
    Idents(sub_rename, cells=WhichCells(sub_tub_rename,idents = "Proximal Tubule"))<-"iPSC_Proximal Tubule"
    Idents(sub_rename, cells=WhichCells(sub_tub_rename,idents = "Collecting Duct"))<-"iPSC_Collecting Duct"
    DimPlot(sub_rename, pt.size = 1)
    DimPlot(sub_rename, pt.size = 1 )+
      theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
      theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
      theme(axis.ticks.length = unit(0.2, "cm"))+
      theme(axis.ticks = element_line(size = 1.2))+
      theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
      theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
      theme(legend.text=element_text(color="black",size=rel(1.5)))+
      theme(legend.title=element_text(color="black",size=rel(1.5)))+
      theme(legend.key.size = unit(0.5, "cm"))
    
    
    # 分离E6中的肾前体细胞
    
    sub_tub <- subset(ipsc_organoid, idents  = c(4))
    DimPlot(sub_tub,pt.size = 1)
    
    
    sub_tub[['percent.mito']] <- PercentageFeatureSet(sub_tub, pattern = "^MT-")
    sub_tub <- NormalizeData(sub_tub, normalization.method = "LogNormalize", scale.factor = 10000)
    sub_tub <- FindVariableFeatures(sub_tub, selection.method = "vst", nfeatures = 2000)
    sub_tub <- ScaleData(sub_tub,  vars.to.regress = c("nCount_RNA", "percent.mito"))
    sub_tub <- RunPCA(sub_tub, npcs = 50)
    sub_tub <- ProjectDim(object = sub_tub)
    
    sub_tub <- FindNeighbors(object = sub_tub, dims = 1:30)
    sub_tub <- FindClusters(object = sub_tub, resolution = 0.6) 
    sub_tub <- RunUMAP(sub_tub, reduction = "pca", dims = 1:30)
    DimPlot(sub_tub, pt.size = 1)
    
    sub_tub_rename <- RenameIdents(sub_tub, `0` = "Epithelium(cycling)", `1` = "Distal Tubule & Loop of Henle", `2` = "Proximal Tubule", 
                                   `3` = "Epithelium(cycling)", `4` = "Proximal Tubule", `5` = "Epithelium(cycling)", `6` = "Epithelium(cycling)", `7` = "Epithelium(cycling)",
                                   `8` = "Collecting Duct")
    DimPlot(sub_tub_rename, pt.size = 1)
    
    save(sub_rename, file = "iPSC_organoid(E6).Robj")}
  # iPSC的小管细胞与KPC类器官的小管细胞整合
  {
    load(file = "ROBJ/ipsc_tubule-only.Robj")
    DimPlot(sub_rename, pt.size = 1)
    iPSC_Tub <- subset(sub_rename,idents = c("Distal Tubule & Loop of Henle","Proximal Tubule","Collecting Duct"))
    DimPlot(iPSC_Tub, pt.size = 1)
    
    rm(sub_tub_rename)
    
    scell.anchors <- FindIntegrationAnchors(object.list = c(iPSC_Tub, KPC_Tub), dims = 1:20)
    scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors)) 
    
    DefaultAssay(scell) <- "RNA"
    scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
    scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
    scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
    scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
    #Integrated（分析逻辑）
    DefaultAssay(scell) <- "integrated"
    scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
    scell <- RunPCA(scell, npcs = 50)
    scell <- ProjectDim(object = scell)
    
    scell <- FindNeighbors(object = scell, dims = 1:30)
    scell <- FindClusters(object = scell, resolution = 0.5) 
    
    scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
    DimPlot(scell,pt.size = 1)
    DimPlot(scell,group.by = "orig.ident",pt.size = 1)
    DimPlot(scell,split.by = "orig.ident",pt.size = 1)
    
    merge <- scell
    rm(scell,scell.anchors)
    
    a <- as.data.frame(iPSC_Tub@active.ident) 
    a[,2] <- rownames(a)
    colnames(a) <- c("idents", "cell_id") #提取数据集A的细胞定义
    
    b <- as.data.frame(KPC_Tub@active.ident)
    b[,2] <- rownames(b)
    colnames(b) <- c("idents","cell_id") #提取数据集B的细胞定义
    c <- rbind(a,b) #将两个数据集的定义合并
    
    rm(a,b)
    
    metadata <- as.data.frame(merge@active.ident) #提取整合后的数据集的metadata信息，用于后续与上面提取的细胞分群信息进行合并
    metadata[,2] <- substr(rownames(metadata),1,18) #将细胞信息作为列表的第二列用于索引
    colnames(metadata) <- c("orig.ident", "cell_id") #对列名进行重命名
    
    
    metadata <- full_join(c,metadata, by=c("cell_id")) #按照cell_id将原始的分群信息与整合后的分群进行替换
    rownames(metadata) <- metadata$cell_id
    
    metadata$idents <- factor(metadata$idents) #将新的metadata数据以factor数据类型重新导入回整合后的单细胞测序数据中
    merge@meta.data[["idents"]] <- metadata$idents # 这里的顺序不一定是对的
    Idents(merge) <- merge@meta.data$idents  
    
    DefaultAssay(merge) <- "RNA"
    DimPlot(merge,pt.size = 0.7,group.by = "seurat_clusters")
    DimPlot(merge,pt.size = 1)
    DimPlot(merge,pt.size = 1,split.by = "orig.ident")
    DimPlot(merge,pt.size = 1, group.by = "orig.ident")
    DimPlot(merge, pt.size = 0.7 )+
      theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
      theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
      theme(axis.ticks.length = unit(0.2, "cm"))+
      theme(axis.ticks = element_line(size = 1.2))+
      theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
      theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
      theme(legend.text=element_text(color="black",size=rel(1.5)))+
      theme(legend.title=element_text(color="black",size=rel(1.5)))+
      theme(legend.key.size = unit(0.5, "cm"))
    
    #计算差异基因
    DEG <- FindAllMarkers(merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.table(DEG, file = paste(wd,"/",time,"/","DEG.csv", sep = ""), sep=",")
    
    # 将active.idents改为orig.ident再进行基因差异分析
    merge_test <- merge
    merge_test@active.ident <- as.factor(merge_test$orig.ident) #orig.ident不是因子类型，所以需要as.factor进行转换后再导入到active,idents中
    DEG_ori <- FindAllMarkers(merge_test, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.table(DEG_ori, file = paste(wd,"/",time,"/","DEG_ori.csv", sep = ""), sep=",")
    
    #展示不同类型的基因表达
    gene_name <- c("SOX9","LHX1") #创建基因列表
    genelist = data.frame(gene_name)
    matrix <- AverageExpression(merge)
    a <- matrix$RNA
    b <- a[which(rownames(a) %in% genelist$gene_name),]
    c <- a[genelist$gene_name,]
    color=colorRampPalette(c("black", "#34366b","#ffe889"))(100)
    d <- c[,c(5,4,3,2,1)] #调整列的顺序
    pheatmap(main="",c,scale="row",cluster_rows=F,cluster_cols=F,
                   color=colorRampPalette(c("black", "#34366b","#ffe889"))(100),
                   border_color="white",angle_col=45,
                   cellwidth = 20, cellheight = 15,gaps_col = c(0))
    
    #GO分析每组的细胞功能，筛选基因用于基因热图
    ##分析iPSC衍生的Collecting duct分群的GO富集分析
    DEG_Collecting_Duct <- DEG[which(DEG$cluster %in% "Collecting Duct"),]
    GO_Collecting_Duct <- enrichGO(DEG_Collecting_Duct$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.05,keyType = "SYMBOL")
    write.table(GO_Collecting_Duct, file = paste(wd,"/",time,"/","GO_Collecting_Duct.csv", sep = ""), sep=",")
    
    ##分析iPSC衍生的Proximal Tubule分群的GO富集分析
    DEG_Proximal_Tubule <- DEG[which(DEG$cluster %in% "Proximal Tubule"),]
    GO_Proximal_Tubule <- enrichGO(DEG_Proximal_Tubule$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                                   qvalueCutoff = 0.05,keyType = "SYMBOL")
    write.table(GO_Proximal_Tubule, file = paste(wd,"/",time,"/","GO_Proximal_Tubule.csv", sep = ""), sep=",")
    
    ##分析iPSC衍生的Distal Tubule & Loop of Henle分群的GO富集分析
    DEG_DTLH <- DEG[which(DEG$cluster %in% "Distal Tubule & Loop of Henle"),]
    GO_DTLH <- enrichGO(DEG_DTLH$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                                   qvalueCutoff = 0.05,keyType = "SYMBOL")
    write.table(GO_DTLH, file = paste(wd,"/",time,"/","GO_DTLH.csv", sep = ""), sep=",")
    
    ##分析KPC衍生的PT分群的GO富集分析
    DEG_PT_like <- DEG[which(DEG$cluster %in% "PT-like"),]
    GO_PT_like <- enrichGO(DEG_PT_like$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05,keyType = "SYMBOL")
    write.table(GO_PT_like, file = paste(wd,"/",time,"/","GO_PT_like.csv", sep = ""), sep=",")
    
    ##分析orig.idents下的DEG GO富集分析结果
    DEG_iPSC <- DEG[which(DEG_ori$cluster %in% "ipsc_organoid"),]
    GO_iPSC <- enrichGO(DEG_iPSC$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05,keyType = "SYMBOL")
    write.table(GO_iPSC, file = paste(wd,"/",time,"/","GO_iPSC.csv", sep = ""), sep=",")
    
    DEG_KPC <- DEG[which(DEG_ori$cluster %in% "KPC"),]
    GO_KPC <- enrichGO(DEG_KPC$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05,keyType = "SYMBOL")
    write.table(GO_KPC, file = paste(wd,"/",time,"/","GO_KPC.csv", sep = ""), sep=",")
  }
  # iPSC与KPC类器官数据进行整合
  {
    load(file = "2022-07-17/KPC.Robj")
    DimPlot(KPC)
    KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
    colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
    KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
    KPC@active.ident <- KPC$Monocle.Cluster
    DimPlot(KPC, pt.size = 1)
    Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"KPC_CD-like"
    Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"KPC_PT-like"
    Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"KPC_LH-like"
    Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"KPC_Progenitor cell"
    Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"KPC_Epithelial cell"
    Idents(KPC,cells=WhichCells(KPC,idents = "6"))<-"KPC_Cell cycle"
    
    
    scell.anchors <- FindIntegrationAnchors(object.list = c(sub_rename, KPC), dims = 1:20)
    scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
    
    DefaultAssay(scell) <- "RNA"
    scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
    scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
    scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
    scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
    #Integrated（分析逻辑）
    DefaultAssay(scell) <- "integrated"
    scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
    scell <- RunPCA(scell, npcs = 50)
    scell <- ProjectDim(object = scell)
    
    scell <- FindNeighbors(object = scell, dims = 1:30)
    scell <- FindClusters(object = scell, resolution = 0.5) 
    
    scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
    DimPlot(scell,pt.size = 1)
    DimPlot(scell,group.by = "orig.ident",pt.size = 1)
    DimPlot(scell,split.by = "orig.ident",pt.size = 1)
    
    merge <- scell
    rm(scell,scell.anchors)
    
    
    a <- as.data.frame(sub_rename@active.ident)
    a[,2] <- rownames(a)
    colnames(a) <- c("idents", "cell_id")
    
    b <- as.data.frame(KPC@active.ident)
    b[,2] <- rownames(b)
    colnames(b) <- c("idents","cell_id")
    c <- rbind(a,b)
    
    rm(a,b)
    
    metadata <- as.data.frame(merge@active.ident)
    metadata[,2] <- substr(rownames(metadata),1,18)
    colnames(metadata) <- c("orig.ident", "cell_id")
    
    
    metadata[13607,2] <- "TGTTCCGCAAGCGAGT-1_1"
    c[13607,2] <- "TGTTCCGCAAGCGAGT-1_1"
    
    metadata <- full_join(c,metadata, by=c("cell_id"))
    rownames(metadata) <- metadata$cell_id
    
    metadata$idents <- factor(metadata$idents)
    merge@meta.data[["idents"]] <- metadata$idents # 这里的顺序不一定是对的
    Idents(merge) <- merge@meta.data$idents  
    
    DefaultAssay(merge) <- "RNA"
    DimPlot(merge,pt.size = 0.7,group.by = "seurat_clusters")
    DimPlot(merge,pt.size = 1)
    DimPlot(merge,pt.size = 1,split.by = "orig.ident")
    DimPlot(merge,pt.size = 1, group.by = "orig.ident")
    DimPlot(merge, pt.size = 0.7 )+
      theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
      theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
      theme(axis.ticks.length = unit(0.2, "cm"))+
      theme(axis.ticks = element_line(size = 1.2))+
      theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
      theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
      theme(legend.text=element_text(color="black",size=rel(1.5)))+
      theme(legend.title=element_text(color="black",size=rel(1.5)))+
      theme(legend.key.size = unit(0.5, "cm"))
    
    save(merge, file = "Merge_iPSC&KPC.Robj")
  }
  # iPSC与KPC类器官整合后跑拟时序
  {
    # 创建CDS文件
    data <- as(as.matrix(merge@assays$RNA@counts), 'sparseMatrix')
    pd <-  merge@meta.data
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    colnames(pd)
    
    cds <- new_cell_data_set(data,
                             cell_metadata  = pd,
                             gene_metadata  = fData)
    #Pre-process the data
    cds = preprocess_cds(cds, method = c("PCA"),num_dim = 100)
    plot_pc_variance_explained(cds)
    
    #Reduce dimensionality and visualize the cells
    cds = reduce_dimension(cds) #Monocle uses UMAP by default 
    plot_cells(cds) 
    
    cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.0003)  
    plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
    
    cds@int_colData$reducedDims$UMAP <- merge@reductions$umap@cell.embeddings
    cds@clusters$UMAP$clusters <- merge@active.ident
    
    
    # 定义细胞群（根据之前分析的小提琴图）
    #Annotate your cells according to type
    
    plot_cells(cds, group_cells_by="cluster", color_cells_by="assigned_cell_type",group_label_size=4,cell_size=1)
    plot_cells(cds,group_cells_by="cluster", color_cells_by="assigned_cell_type",
               graph_label_size=2,cell_size=2,label_cell_groups=FALSE,show_trajectory_graph = F, )+
      theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
      theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
      theme(axis.ticks.length = unit(0.2, "cm"))+
      theme(axis.ticks = element_line(size = 1.2))+
      theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
      theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
      theme(legend.text=element_text(color="black",size=rel(1.5)))+
      theme(legend.title=element_text(color="black",size=rel(1.5)))+
      theme(legend.key.size = unit(0.5, "cm"))
    plot_cells(cds, color_cells_by = "cluster", cell_size=1, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE,show_trajectory_graph = F)+
      theme(axis.line.x=element_line(linetype=1,color="black",size=1.5))+
      theme(axis.line.y=element_line(linetype=1,color="black",size=1.5))+
      theme(axis.ticks.length = unit(0.2, "cm"))+
      theme(axis.ticks = element_line(size = 1.2))+
      theme(axis.title.x = element_text(color="black", size=rel(2)))+
      theme(axis.title.y = element_text(color="black", size=rel(2)))+
      theme(axis.text.x = element_text(color="black", size=rel(2)))+
      theme(axis.text.y = element_text(color="black", size=rel(2)))+
      theme(legend.text=element_text(color="black",size=rel(2)))+
      theme(legend.title=element_text(color="black",size=rel(2)))+
      theme(legend.key.size = unit(0.5, "cm"))
    colnames(colData(cds))[8] <- "Cell type"
    
    save(cds,file = paste(wd,"/",time,"/","organoid_monocle3_cds.Robj", sep = ""))
    
    
    # 轨迹学习Learn the trajectory graph（使用learn_graph()函数）
    cds <- learn_graph(cds,use_partition = F)
    plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE)
    
    # 2. 细胞按拟时排序
    cds <- order_cells(cds)
    plot_cells(cds, color_cells_by = "pseudotime", cell_size=1, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE,)
    
    plot_cells(cds, color_cells_by = "pseudotime", cell_size=1, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)+
      theme(axis.line.x=element_line(linetype=1,color="black",size=1.5))+
      theme(axis.line.y=element_line(linetype=1,color="black",size=1.5))+
      theme(axis.ticks.length = unit(0.2, "cm"))+
      theme(axis.ticks = element_line(size = 1.2))+
      theme(axis.title.x = element_text(color="black", size=rel(2)))+
      theme(axis.title.y = element_text(color="black", size=rel(2)))+
      theme(axis.text.x = element_text(color="black", size=rel(2)))+
      theme(axis.text.y = element_text(color="black", size=rel(2)))+
      theme(legend.text=element_text(color="black",size=rel(2)))+
      theme(legend.title=element_text(color="black",size=rel(2)))+
      theme(legend.key.size = unit(0.5, "cm"))
    
  }
# GO富集分析-----
  library(DOSE)
  library(enrichplot)
  library("clusterProfiler")
  library("org.Hs.eg.db")
  data(geneList)
  de <- names(geneList)[abs(geneList) > 2]
  
  edo <- enrichDGN(de)
  
  DimPlot(KPC)
  Idents(KPC,cells = WhichCells(KPC, idents = "KPC_Epithelial cell")) <- "KPC_Differentiated cell"
  Idents(KPC,cells = WhichCells(KPC, idents = "KPC_LH-like")) <- "KPC_Differentiated cell"
  Idents(KPC,cells = WhichCells(KPC, idents = "KPC_PT-like")) <- "KPC_Differentiated cell"
  Idents(KPC,cells = WhichCells(KPC, idents = "KPC_CD-like")) <- "KPC_Differentiated cell"
  Idents(KPC,cells = WhichCells(KPC, idents = "KPC_Cell cycle")) <- "KPC_Cycling cell"
  
  KPC2 <- subset(KPC,idents = c("KPC_Differentiated cell","KPC_Progenitor cell"))
  DimPlot(KPC2)
  
  DEG_LPP <- FindAllMarkers(KPC2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  DEG_Differentiated <- DEG_LPP[which(DEG_LPP$cluster %in% "KPC_Differentiated cell"),]
  DEG_Progenitor <- DEG_LPP[which(DEG_LPP$cluster %in% "KPC_Progenitor cell"),]
  rownames(DEG_Differentiated) <-  DEG_Differentiated$gene
  rownames(DEG_Progenitor) <-  DEG_Progenitor$gene
  
  GO_Differentiated <- enrichGO(DEG_Differentiated$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                                qvalueCutoff = 0.05,keyType = "SYMBOL")
  write.table(GO_Differentiated, file = paste(wd,"/",time,"/","GO_Differentiated.csv", sep = ""), sep=",")
  
  barplot(GO_Differentiated)
  genelist <- as.numeric(DEG_Differentiated[,2]) 
  names(genelist) <- row.names(DEG_Differentiated)  
  
  cnetp1 <- cnetplot(GO_Differentiated,  
                     foldChange = genelist,
                     showCategory = 6,
                     colorEdge = T,
                     node_label = 'all',
                     color_category ='steelblue')
  
  cnetplot(GO_Differentiated, node_label="all") 
  
# 20230323 将KPC单细胞测序与不同时期胎儿肾脏数据进行整合-----
  #读取不同时期胎儿肾脏
  Fetal_w9.data <- Read10X(data.dir = "D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC/Fetal kidney single cell data/GSE114530_RAW/w9/")
  Fetal_w11.data <- Read10X(data.dir = "D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC/Fetal kidney single cell data/GSE114530_RAW/w11")
  Fetal_w13.data <- Read10X(data.dir = "D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC/Fetal kidney single cell data/GSE114530_RAW/w13/")
  Fetal_w16.data <- Read10X(data.dir = "D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC/Fetal kidney single cell data/GSE114530_RAW/w16/")
  Fetal_w18.data <- Read10X(data.dir = "D:/Zuo Lab/Data/20210325 sequence Data/KPC 3D/KPC/Fetal kidney single cell data/GSE114530_RAW/w18/")
  Fetal_w9 <- CreateSeuratObject(counts = Fetal_w9.data, project = "Fetal_w9", min.cells = 3, min.features = 200)
  Fetal_w11 <- CreateSeuratObject(counts = Fetal_w11.data, project = "Fetal_w11", min.cells = 3, min.features = 200)
  Fetal_w13 <- CreateSeuratObject(counts = Fetal_w13.data, project = "Fetal_w13", min.cells = 3, min.features = 200)
  Fetal_w16 <- CreateSeuratObject(counts = Fetal_w16.data, project = "Fetal_w16", min.cells = 3, min.features = 200)
  Fetal_w18 <- CreateSeuratObject(counts = Fetal_w18.data, project = "Fetal_w18", min.cells = 3, min.features = 200)
  rm(Fetal_w9.data,Fetal_w11.data,Fetal_w13.data,Fetal_w16.data,Fetal_w18.data)
  
  w9_idents <- read.csv(file = "Fetal kidney single cell data/GSE114530_RAW/w9/w9_barcodes_celltypes.csv",header = T)
  w11_idents <- read.csv(file = "Fetal kidney single cell data/GSE114530_RAW/w11/w11_barcodes_celltypes.csv",header = T)
  w13_idents <- read.csv(file = "Fetal kidney single cell data/GSE114530_RAW/w13/w13_barcodes_celltypes.csv",header = T)
  w16_idents <- read.csv(file = "Fetal kidney single cell data/GSE114530_RAW/w16/w16_barcodes_celltypes.csv",header = T)
  w18_idents <- read.csv(file = "Fetal kidney single cell data/GSE114530_RAW/w18/w18_barcodes_celltypes.csv",header = T)
  Fetal_w9 <- subset(Fetal_w9, cells = WhichCells(Fetal_w9, cells = w9_idents$cell.barcode))
  Fetal_w11 <- subset(Fetal_w11, cells = WhichCells(Fetal_w11, cells = w11_idents$cell.barcode))
  Fetal_w13 <- subset(Fetal_w13, cells = WhichCells(Fetal_w13, cells = w13_idents$cell.barcode))
  Fetal_w16 <- subset(Fetal_w16, cells = WhichCells(Fetal_w16, cells = w16_idents$cell.barcode))
  Fetal_w18 <- subset(Fetal_w18, cells = WhichCells(Fetal_w18, cells = w18_idents$cell.barcode))
  rm(w9_idents,w11_idents,w13_idents,w16_idents,w18_idents)
  
  
  Fetal_w9 <- NormalizeData(Fetal_w9, normalization.method = "LogNormalize", scale.factor = 10000)
  Fetal_w9 <- FindVariableFeatures(Fetal_w9, selection.method = "vst", nfeatures = 2000)
  Fetal_w11 <- NormalizeData(Fetal_w11, normalization.method = "LogNormalize", scale.factor = 10000)
  Fetal_w11 <- FindVariableFeatures(Fetal_w11, selection.method = "vst", nfeatures = 2000)
  Fetal_w13 <- NormalizeData(Fetal_w13, normalization.method = "LogNormalize", scale.factor = 10000)
  Fetal_w13 <- FindVariableFeatures(Fetal_w13, selection.method = "vst", nfeatures = 2000)
  Fetal_w16 <- NormalizeData(Fetal_w16, normalization.method = "LogNormalize", scale.factor = 10000)
  Fetal_w16 <- FindVariableFeatures(Fetal_w16, selection.method = "vst", nfeatures = 2000)
  Fetal_w18 <- NormalizeData(Fetal_w18, normalization.method = "LogNormalize", scale.factor = 10000)
  Fetal_w18 <- FindVariableFeatures(Fetal_w18, selection.method = "vst", nfeatures = 2000)
  
  scell.anchors <- FindIntegrationAnchors(object.list = c(KPC,Fetal_w16,Fetal_w11,Fetal_w18,Fetal_w13,Fetal_w9), dims = 1:20)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20,features.to.integrate = rownames(scell.anchors))
  scell.anchors <- FindIntegrationAnchors(object.list = c(scell,Fetal_w13), dims = 1:20)
  
  scell.anchors <- FindIntegrationAnchors(object.list = c(scell,Fetal_w16), dims = 1:20)
  
  scell.anchors <- FindIntegrationAnchors(object.list = c(Fetal_w9,Fetal_w11), dims = 1:20)
  scell.anchors <- FindIntegrationAnchors(object.list = c(Fetal_w9,Fetal_w11), dims = 1:20)
  
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:20, features.to.integrate = rownames(scell.anchors))
  
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  #Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(scell)
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)
  DimPlot(scell,group.by = "orig.ident",pt.size = 1)
  DimPlot(scell,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  merge <- scell
  rm(scell,scell.anchors)
  
  
# 20230423 GO富集分析，分别分析KPC和iPSC类器官的肾源性细胞的GO（sub信息在上面）-------------
  # KPC的GO分析
  DimPlot(KPC_nephron)
  KPC_nephron[["Cluster"]] <- KPC_nephron@active.ident
  #修改
  #计算差异基因
  DEG_KPC_nephron <- FindAllMarkers(KPC_nephron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(DEG_KPC_nephron, file = paste(wd,"/",time,"/","DEG_KPC_nephron.csv", sep = ""), sep=",")
  #分析KPC_nephron GO富集分析结果
  DEG_KPC_Pro <- DEG_KPC_nephron[which(DEG_KPC_nephron$cluster %in% "Progenitor cell"),]
  GO_KPC_pro <- enrichGO(DEG_KPC_Pro$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05,keyType = "SYMBOL")
  write.table(GO_KPC_pro, file = paste(wd,"/",time,"/","GO_KPC_pro.csv", sep = ""), sep=",")
  
  DEG_KPC_PT <- DEG_KPC_nephron[which(DEG_KPC_nephron$cluster %in% "PT-like"),]
  GO_KPC_PT <- enrichGO(DEG_KPC_PT$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05,keyType = "SYMBOL")
  write.table(GO_KPC_PT, file = paste(wd,"/",time,"/","GO_KPC_PT.csv", sep = ""), sep=",")
  
  DEG_KPC_LHDT <- DEG_KPC_nephron[which(DEG_KPC_nephron$cluster %in% "LH-like"),]
  GO_KPC_LHDT <- enrichGO(DEG_KPC_LHDT$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05,keyType = "SYMBOL")
  write.table(GO_KPC_LHDT, file = paste(wd,"/",time,"/","GO_KPC_LHDT.csv", sep = ""), sep=",")
  
  
  DEG_KPC_CD <- DEG_KPC_nephron[which(DEG_KPC_nephron$cluster %in% "CD-like"),]
  GO_KPC_CD <- enrichGO(DEG_KPC_CD$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05,keyType = "SYMBOL")
  write.table(GO_KPC_CD, file = paste(wd,"/",time,"/","GO_KPC_CD.csv", sep = ""), sep=",")
  
  # 绘制热图
  colo <- RColorBrewer::brewer.pal(10,"RdYlBu")[1:5]
  colo2 <- RColorBrewer::brewer.pal(10,"RdYlBu")[6:10]
  top30 <- DEG_KPC_nephron %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
  DoHeatmap(KPC_nephron, features = top30$gene,group.colors = )
  DoHeatmap(KPC_nephron, disp.min = -2,disp.max = 2.5, #调整区间使得颜色更明显
            features = as.character(unique(top30$gene)),   
            group.by = "Cluster",  
            assay = "RNA",  
            group.colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#AB82FF","#90EE90","#00CD00","#008B8B"))+  #设置组别颜色
    scale_fill_gradientn(colors = c("black", "#34366b","#ffe889"))+
    theme(legend.text=element_text(color="black",size=rel(1.5)), # 图例字体大小及颜色
          legend.title=element_text(color="black",size=rel(1.5)), # 图例标题大小及颜色
          legend.key.heigh = unit(0.5, "cm"), # 图例图示部分的高度
          legend.key.width = unit(0.7, "cm")) # 图例图示部分的宽度
  
# 20230425 下载AKI患者的肾穿刺样本和正常人的穿刺/活检样本-----
  cat(readLines("AKI kidney and control/GSE210622_metadata.txt"), sep ="\
") #查看txt文本数据的数据类型
  metadata <- read.table(file = "AKI kidney and control/GSE210622_metadata.txt", header = F, fill = TRUE, sep =" ")
  metadata[,8] <- substr(metadata[,1],nchar(metadata[,1])-17,nchar(metadata[,1])) #将细胞信息作为列表的第二列用于索引
  metadata[,9] <- paste(metadata[,8], "1",sep = "_")
  colnames(metadata) <- c("cell_name","orig.ident","nCount_RNA","nFeature_RNA","percent.mt","celltype","group","cell_name2","cell_name_1")
  metadata <- metadata[-1,]
  
  # 读取三个对照肾脏核单细胞测序数据
  sc.data <- Read10X(data.dir = "AKI kidney and control/GSE210622_RAW/CT1/")
  sc.data2 <- Read10X(data.dir = "AKI kidney and control/GSE210622_RAW/CT2/")
  sc.data3 <- Read10X(data.dir = "AKI kidney and control/GSE210622_RAW/CT3/")
  CT1 <- CreateSeuratObject(counts = sc.data, project = "CT1", min.cells = 3, min.features = 200)
  CT2 <- CreateSeuratObject(counts = sc.data, project = "CT2", min.cells = 3, min.features = 200)
  CT3 <- CreateSeuratObject(counts = sc.data, project = "CT3", min.cells = 3, min.features = 200)
  
  CT1[["percent.mt"]] <- PercentageFeatureSet(CT1, pattern = "^MT-")  
  CT2[["percent.mt"]] <- PercentageFeatureSet(CT2, pattern = "^MT-")  
  CT3[["percent.mt"]] <- PercentageFeatureSet(CT3, pattern = "^MT-")  
  VlnPlot(CT1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(CT2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(CT3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # 根据小提琴图选择筛选阈值
  CT1 <- subset(CT1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 2.5)
  CT2 <- subset(CT2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 2.5)
  CT3 <- subset(CT3, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 2.5)
  
  CT1 <- NormalizeData(CT1, normalization.method = "LogNormalize", scale.factor = 10000)
  CT2 <- NormalizeData(CT2, normalization.method = "LogNormalize", scale.factor = 10000)
  CT3 <- NormalizeData(CT3, normalization.method = "LogNormalize", scale.factor = 10000)
  
  CT1 <- FindVariableFeatures(CT1, selection.method = "vst", nfeatures = 2000)
  CT2 <- FindVariableFeatures(CT2, selection.method = "vst", nfeatures = 2000)
  CT3 <- FindVariableFeatures(CT3, selection.method = "vst", nfeatures = 2000)  
  rm(sc.data, sc.data2, sc.data3)
  

  CT1 <- ScaleData(CT1, vars.to.regress = c("nCount_RNA", "percent.mito"))
  CT1 <- RunPCA(CT1, npcs = 50)
  CT1 <- ProjectDim(object = CT1)
  
  CT1 <- FindNeighbors(object = CT1, dims = 1:30)
  CT1 <- FindClusters(object = CT1, resolution = 0.5) 
  
  CT1 <- RunUMAP(CT1, reduction = "pca", dims = 1:30)
  DimPlot(CT1,pt.size = 1)
  DimPlot(CT1,group.by = "orig.ident",pt.size = 1)
  DimPlot(CT1,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  CT1@meta.data[["Cluster"]] <- metadata[which(rownames(CT1@meta.data) %in% metadata[,8]),6]
  metadata[which(rownames(CT@meta.data) %in% metadata[,9]),6]
  metadata_nona <- na.omit(metadata)
  metadata_nona <- metadata[-which( metadata[,1] =="COVID"), ]
  metadata_nona <- metadata_nona[-which( metadata_nona[,1] =="Non-COVID"), ]
  metadata_Control_TN1 <- metadata_nona[which( metadata_nona[,2] =="Control-TN1"), ]
  
  # 根据metadata挑选单细胞数据中的细胞
  CT1_sub <- subset(CT1, cells = metadata_Control_TN1[,8])
  DimPlot(CT1_sub,pt.size = 1)
  
  # 导入已知的细胞定义
  #CT1_sub@meta.data[["Cluster"]] <- metadata_Control_TN1[which(rownames(CT1_sub@meta.data) %in% metadata_Control_TN1[,8]),6]
  #CT1_sub@active.ident <- as.factor(CT1_sub$Cluster)
  #CT1_sub@active.ident <- 
  #DimPlot(CT1_sub)
  
  a <- as.data.frame(CT1_sub@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents","cell_id")
  
  b <- as.data.frame(metadata_Control_TN1)
  rownames(b) <- b[,8]
  b <- b[,c(6,8)]
  colnames(b) <- c("idents","cell_id")
  d <- rbind(a,b)
  c <- left_join(a, b, by=c("cell_id")) 
  
  rownames(c) <- c[,2]

  metadata <- as.data.frame(CT1_sub@active.ident)
  metadata[,2] <- substr(rownames(metadata),1,18)
  colnames(metadata) <- c("orig.ident", "cell_id")
  
  metadata <- full_join(c,metadata, by=c("cell_id"))
  rownames(metadata) <- metadata$cell_id
  
  metadata$idents.y <- factor(metadata$idents.y)
  CT1_sub@meta.data[["Cluster"]] <- metadata$idents.y # 这里的顺序不一定是对的
  Idents(CT1_sub) <- CT1_sub@meta.data$Cluster
  DimPlot(CT1_sub,pt.size = 1)
  save(CT1_sub, file = paste(wd,"/",time,"/","CT1_sub.Robj", sep = ""))
  
  #Integrate 使用“锚”关联各组细胞，消除批次效应
  scell.anchors <- FindIntegrationAnchors(object.list = c(CT1,CT2,CT3), dims = 1:20,k.anchor = 5,k.filter = 10)
  CT <- IntegrateData(anchorset = scell.anchors, dims = 1:50,features.to.integrate = rownames(scell.anchors)) 
  rm(scell.anchors)
  # 该处理完后会自动生成一个“integrated”在Assay中，整合使用的RNA并非所有的RNA，因此作图时需要切换，integrated只用于跑出UMAP，即逻辑分析
  
  # 将默认的Assay定义到“RNA”用于作图（Featureplot，Vlnplot）
  DefaultAssay(CT) <- "RNA"
  CT[['percent.mito']] <- PercentageFeatureSet(CT, pattern = "^MT-")
  CT <- NormalizeData(object = CT, normalization.method = "LogNormalize", scale.factor = 10000)
  CT <- FindVariableFeatures(object = CT, selection.method = "vst", nfeatures = 2000)
  CT <- ScaleData(CT,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  # Integrated（分析逻辑）
  DefaultAssay(CT) <- "integrated"
  CT <- ScaleData(CT, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(CT)
  CT <- RunPCA(CT, npcs = 50)
  CT <- ProjectDim(object = CT)
  
  CT <- FindNeighbors(object = CT, dims = 1:30)
  CT <- FindClusters(object = CT, resolution = 0.5) 
  
  CT <- RunUMAP(CT, reduction = "pca", dims = 1:30)
  DimPlot(CT,pt.size = 1)
  DimPlot(CT,group.by = "orig.ident",pt.size = 1)
  DimPlot(CT,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  CT@meta.data[["Cluster"]] <- metadata[which(rownames(CT@meta.data) %in% metadata[,9]),6]
  metadata[which(rownames(CT@meta.data) %in% metadata[,9]),6]
  CT_sub <- subset(CT, )
  
  a <- as.data.frame(CT@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents", "cell_id")
  
  b <- as.data.frame(metadata)
rownames(b) <- b[,9]
  colnames(b) <- c("idents","cell_id")
  c <- rbind(a,b)
  
  rm(a,b)
  
  metadata <- as.data.frame(merge@active.ident)
  metadata[,2] <- substr(rownames(metadata),1,18)
  colnames(metadata) <- c("orig.ident", "cell_id")
  
  
  metadata[13607,2] <- "TGTTCCGCAAGCGAGT-1_1"
  c[13607,2] <- "TGTTCCGCAAGCGAGT-1_1"
  
  metadata <- full_join(c,metadata, by=c("cell_id"))
  rownames(metadata) <- metadata$cell_id
  
  metadata$idents <- factor(metadata$idents)
  merge@meta.data[["idents"]] <- metadata$idents # 这里的顺序不一定是对的
  Idents(merge) <- merge@meta.data$idents 
# 20230430 AKI患者的肾脏样本分析-----
  # 读取三个对照肾脏核单细胞测序数据
  sc.data <- Read10X(data.dir = "AKI kidney and control/GSE210622_RAW/AKI5/")
  AKI5 <- CreateSeuratObject(counts = sc.data, project = "AKI5", min.cells = 3, min.features = 200)
  
  AKI5[["percent.mt"]] <- PercentageFeatureSet(AKI5, pattern = "^MT-")  
  VlnPlot(AKI5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # 根据小提琴图选择筛选阈值
  AKI5 <- subset(AKI5, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
  AKI5 <- NormalizeData(AKI5, normalization.method = "LogNormalize", scale.factor = 10000)
  AKI5 <- FindVariableFeatures(AKI5, selection.method = "vst", nfeatures = 2000)
  rm(sc.data, sc.data2, sc.data3)
  
  AKI5 <- ScaleData(AKI5, vars.to.regress = c("nCount_RNA", "percent.mito"))
  AKI5 <- RunPCA(AKI5, npcs = 50)
  AKI5 <- ProjectDim(object = AKI5)
  
  AKI5 <- FindNeighbors(object = AKI5, dims = 1:30)
  AKI5 <- FindClusters(object = AKI5, resolution = 0.5) 
  
  AKI5 <- RunUMAP(AKI5, reduction = "pca", dims = 1:30)
  DimPlot(AKI5,pt.size = 1)
  DimPlot(AKI5,group.by = "orig.ident",pt.size = 1)
  DimPlot(AKI5,split.by = "orig.ident",pt.size = 1,ncol = 3)
  
  metadata_nona <- metadata[-which( metadata[,1] =="COVID"), ]
  metadata_nona <- metadata_nona[-which( metadata_nona[,1] =="Non-COVID"), ]
  metadata_AKI <- metadata_nona[which( metadata_nona[,2] =="AKI"), ]
  metadata_AKI5 <- metadata_AKI[which( metadata_AKI[,3] =="5"), ]
  
  AKI5_sub <- subset(AKI5, cells = metadata_AKI5[,8])
  DimPlot(AKI5_sub,pt.size = 1)
  
  
  a <- as.data.frame(AKI5_sub@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents","cell_id")
  
  b <- as.data.frame(metadata_AKI5)
  rownames(b) <- b[,8]
  b <- b[,c(7,8)]
  colnames(b) <- c("idents","cell_id")
  d <- rbind(a,b)
  c <- left_join(a, b, by=c("cell_id")) 
  
  rownames(c) <- c[,2]
  
  metadata <- as.data.frame(AKI5_sub@active.ident)
  metadata[,2] <- substr(rownames(metadata),1,18)
  colnames(metadata) <- c("orig.ident", "cell_id")
  
  metadata <- full_join(c,metadata, by=c("cell_id"))
  rownames(metadata) <- metadata$cell_id
  
  metadata$idents.y <- factor(metadata$idents.y)
  AKI5_sub@meta.data[["Cluster"]] <- metadata$idents.y 
  Idents(AKI5_sub) <- AKI5_sub@meta.data$Cluster
  DimPlot(AKI5_sub,pt.size = 1)


# 20230430 分析与肾发育相关的基因在ipsc organoid和KPC organoid之间有无差异-----
  ipsc_gene <- read.table(file = "2023-04-27/iPSC kidney development.txt", sep = "/")
  ipsc_gene <- t(ipsc_gene)
  duplicated(ipsc_gene)
  ipsc_gene <- as.data.frame(ipsc_gene[!duplicated(ipsc_gene),])
  colnames(ipsc_gene) <- "gene"
  
  KPC_gene <- read.table(file = "2023-04-23/KPC_kidney development.txt", sep = "/")
  KPC_gene <- t(KPC_gene)
  duplicated(KPC_gene)
  KPC_gene <- as.data.frame(KPC_gene[!duplicated(KPC_gene),])
  colnames(KPC_gene) <- "gene"
  
  # 绘制韦恩图展示是否有差异基因
  library(ggVennDiagram)
  inner_join(ipsc_gene, KPC_gene,"gene")
  gene_KPC_only <- anti_join(KPC_gene, ipsc_gene, "gene")  
  gene_KPC_only <- as.data.frame(gene_KPC_only[-2,])
  x <- list(Y = KPC_gene$gene, 
            Z = ipsc_gene$gene)
  ggVennDiagram(x)
  
  DimPlot(AKI5_sub,pt.size = 0.8)+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  DimPlot(CT1_sub,pt.size = 0.8)+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "CD-IC-A"))<-"Collecting duct"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "CD-IC-B"))<-"Collecting duct"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "CD-PC"))<-"Collecting duct"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "CNT"))<-"Collecting duct"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "DCT"))<-"Distal tubule"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "tL"))<-"Proximal tubule"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "TAL"))<-"Loop of Henle"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "IntC"))<-"Interstitial cells"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "Leuk"))<-"leukocytes"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "EC"))<-"Epithelial cell"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "PT"))<-"Proximal tubule"
  Idents(AKI5_sub,cells=WhichCells(AKI5_sub,idents = "Podo"))<-"Podocytes"
  
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "CD-IC-A"))<-"Collecting duct"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "CD-IC-B"))<-"Collecting duct"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "CD-PC"))<-"Collecting duct"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "CNT"))<-"Collecting duct"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "DCT"))<-"Distal tubule"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "tL"))<-"Proximal tubule"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "TAL"))<-"Loop of Henle"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "IntC"))<-"Interstitial cells"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "Leuk"))<-"leukocytes"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "EC"))<-"Epithelial cell"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "PT"))<-"Proximal tubule"
  Idents(CT1_sub,cells=WhichCells(CT1_sub,idents = "Podo"))<-"Podocytes"
  
  
  FeaturePlot(AKI5_sub,features = gene_KPC_only$gene)
  FeaturePlot(CT1_sub,features = gene_KPC_only$gene)
  
  colo <- RColorBrewer::brewer.pal(10,"RdYlBu")[1:5]
  colo2 <- RColorBrewer::brewer.pal(10,"RdYlBu")[6:10]
  for (n in gene_KPC_only$gene) {
    p1 <- FeaturePlot(AKI5_sub,features = n,sort.cell = T,pt.size =1,min.cutoff = "q1",max.cutoff = "q95")+ggtitle(n)+
      scale_colour_gradient(low = rev(colo2),high = rev(colo))+ 
      theme(plot.title = element_text(hjust = 0.5,vjust =-0.6,size=rel(2.5)))+
      theme(axis.line = element_line(size=1, colour = "black"))+
      theme(legend.text=element_text(color="black",size=rel(1.2)))+
      theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.title.y = element_text(color="black", size=rel(1.5)))
    
    CairoPDF(paste(wd,"/",time,"/","AKI5_KPC_only_",n,".pdf", sep = ""), width=8, height=8)
    plot(p1)
    dev.off()
  }

  for (n in gene_KPC_only$gene) {
    p1 <- FeaturePlot(CT1_sub,features = n,sort.cell = T,pt.size =1,min.cutoff = "q1",max.cutoff = "q95")+ggtitle(n)+
      scale_colour_gradient(low = rev(colo2),high = rev(colo))+ 
      theme(plot.title = element_text(hjust = 0.5,vjust =-0.6,size=rel(2.5)))+
      theme(axis.line = element_line(size=1, colour = "black"))+
      theme(legend.text=element_text(color="black",size=rel(1.2)))+
      theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
      theme(axis.title.y = element_text(color="black", size=rel(1.5)))
    
    CairoPDF(paste(wd,"/",time,"/","CT1_KPC_only_",n,".pdf", sep = ""), width=8, height=8)
    plot(p1)
    dev.off()
  }
# 20230502 整合CT1和AKI5来探究基因表达情况
  #Integrate 使用“锚”关联各组细胞，消除批次效应
  scell.anchors <- FindIntegrationAnchors(object.list = c(AKI5_sub,CT1_sub), dims = 1:20,k.anchor = 5,k.filter = 10)
  scell <- IntegrateData(anchorset = scell.anchors, dims = 1:50,features.to.integrate = rownames(scell.anchors)) 
  rm(scell.anchors)
  DefaultAssay(scell) <- "RNA"
  scell[['percent.mito']] <- PercentageFeatureSet(scell, pattern = "^MT-")
  scell <- NormalizeData(object = scell, normalization.method = "LogNormalize", scale.factor = 10000)
  scell <- FindVariableFeatures(object = scell, selection.method = "vst", nfeatures = 2000)
  scell <- ScaleData(scell,  vars.to.regress = c("nCount_RNA", "percent.mito")) 
  # Integrated（分析逻辑）
  DefaultAssay(scell) <- "integrated"
  scell <- ScaleData(scell, vars.to.regress = c("nCount_RNA", "percent.mito"))
  ElbowPlot(scell)
  scell <- RunPCA(scell, npcs = 50)
  scell <- ProjectDim(object = scell)
  
  scell <- FindNeighbors(object = scell, dims = 1:30)
  scell <- FindClusters(object = scell, resolution = 0.5) 
  
  scell <- RunUMAP(scell, reduction = "pca", dims = 1:30)
  DimPlot(scell,pt.size = 1)  
  
  
  ##展示具有代表性的基因
  library(pheatmap)
  matrix <- AverageExpression(scell)
  a <- matrix$RNA
  b <- a[which(rownames(a) %in% gene_KPC_only[,1]),]
  c <- a[gene_KPC_only$gene_name,]
  color=colorRampPalette(c("black", "#34366b","#ffe889"))(100)
  d <- c[,c(5,4,3,2,1)]
  pheatmap(main="",b,scale="row",cluster_rows=F,cluster_cols=F,
                 color=colorRampPalette(c("black", "#34366b","#ffe889"))(100),
                 border_color="white",angle_col=45,
                 cellwidth = 20, cellheight = 15,gaps_col = c(0))  
  
  # 定义细胞类型
  metadata_AKI5_CT1 <- rbind(metadata_AKI5,metadata_Control_TN1)
  
  a <- as.data.frame(scell@active.ident)
  a[,2] <- rownames(a)
  colnames(a) <- c("idents","cell_id")
  

  for (i in c(1:length(metadata_AKI5_CT1$group))) {
    if (metadata_AKI5_CT1[i,7] == "Control") {
      metadata_AKI5_CT1[i,9] <- paste(metadata_AKI5_CT1[i,8], "2",sep = "_")
    }
  }
  
  b <- as.data.frame(metadata_AKI5_CT1)
  
  rownames(b) <- b[,9]
  b <- b[,c(6,7,8,9)]
  colnames(b) <- c("idents_AKI","idents_Control","cell_id","cell_id-n")
  c <- left_join(a, b, by=c("cell_id" = "cell_id-n")) 
  
  rownames(c) <- c[,2]
  
  metadata <- as.data.frame(scell@active.ident)
  metadata[,2] <- rownames(metadata)
  colnames(metadata) <- c("orig.ident", "cell_id")
  
  metadata <- full_join(c,metadata, by=c("cell_id"))

  metadata[1:7241,7] <- metadata[1:7241,4]
  metadata[7242:17540,7] <- metadata[7242:17540,3]
  metadata[1:7241,8] <- "AKI5"
  metadata[7242:17540,8] <- "Control"
  rownames(metadata) <- metadata$cell_id
  colnames(metadata)[7] <- "Cluster"
  colnames(metadata)[8] <- "Sample"
  
  metadata$Cluster <- factor(metadata$Cluster)
  metadata$Sample <- factor(metadata$Sample)
  scell@meta.data[["Cluster"]] <- metadata$Cluster 
  scell@meta.data[["Sample"]] <- metadata$Sample 
  Idents(scell) <- scell@meta.data$Cluster
  Idents(scell) <- scell@meta.data$Sample
  DimPlot(scell,pt.size = 1,group.by = "Sample")+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.3))+
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.3))+
    theme(axis.ticks.length = unit(0.2, "cm"))+
    theme(axis.ticks = element_line(size = 1.2))+
    theme(axis.title.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.x = element_text(color="black", size=rel(1.5)))+
    theme(axis.text.y = element_text(color="black", size=rel(1.5)))+
    theme(legend.text=element_text(color="black",size=rel(1.5)))+
    theme(legend.title=element_text(color="black",size=rel(1.5)))+
    theme(legend.key.size = unit(0.5, "cm"))
  
  DimPlot(scell,pt.size = 1,group.by = "Idents",split.by = "Sample")
  
  FeaturePlot(scell,features = "SOX9",sort.cell = T,pt.size =1,min.cutoff = "q30",max.cutoff = "q95",split.by = "Sample")
  FeaturePlot(scell,features = "ITGA3",sort.cell = T,pt.size =1,min.cutoff = "q50",max.cutoff = "q99",split.by = "Sample")
  
  colo <- RColorBrewer::brewer.pal(10,"RdYlBu")[1:5]
  colo2 <- RColorBrewer::brewer.pal(10,"RdYlBu")[6:10]
  FeaturePlot(scell,features = "SOX9",sort.cell = T,pt.size =1.5,min.cutoff = "q30",max.cutoff = "q95")+
    scale_colour_gradient(low = rev(colo2),high = rev(colo))+ 
    theme(plot.title = element_text(hjust = 0.5,vjust =-0.6,size=rel(2.7)))+
    theme(axis.line = element_line(size=1, colour = "black"))+
    theme(legend.text=element_text(color="black",size=rel(1.3)))+
    theme(axis.title.x = element_text(color="black", size=rel(1.6)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.6)))
  FeaturePlot(scell,features = "IRX3",sort.cell = T,pt.size =1,min.cutoff = "q1",max.cutoff = "q95",split.by = "Sample") #差异不明显
  FeaturePlot(scell,features = "VEGFA",sort.cell = T,pt.size =1,min.cutoff = "q1",max.cutoff = "q95",split.by = "Sample") #差异不明显
  FeaturePlot(scell,features = "NDUFS6",sort.cell = T,pt.size =1,min.cutoff = "q1",max.cutoff = "q95",split.by = "Sample") #差异不明显
  FeaturePlot(scell,features = "ITGB4",sort.cell = T,pt.size =1,min.cutoff = "q1",max.cutoff = "q95",split.by = "Sample") #差异不明显
  FeaturePlot(scell,features = "PLAT",sort.cell = T,pt.size =1,min.cutoff = "q1",max.cutoff = "q95",split.by = "Sample") #差异不明显
  
  # subset整合后的数据分别跑Featureplot-
  Idents(scell) <- scell@meta.data$Sample
  DimPlot(scell)
  Merge_AKI5 <- subset(scell, idents = "AKI5")
  Merge_Control <- subset(scell, idents = "Control")
  
  FeaturePlot(Merge_Control,features = "ITGA3",sort.cell = T,pt.size =1.5,min.cutoff = "q30",max.cutoff = "q98")+
    scale_colour_gradient(low = rev(colo2),high = rev(colo))+ 
    theme(plot.title = element_text(hjust = 0.5,vjust =-0.6,size=rel(2.7)))+
    theme(axis.line = element_line(size=1, colour = "black"))+
    theme(legend.text=element_text(color="black",size=rel(1.3)))+
    theme(axis.title.x = element_text(color="black", size=rel(1.6)))+
    theme(axis.title.y = element_text(color="black", size=rel(1.6)))
  
  Idents(scell,cells=WhichCells(scell,idents = "CD-IC-A"))<-"Collecting duct"
  Idents(scell,cells=WhichCells(scell,idents = "CD-IC-B"))<-"Collecting duct"
  Idents(scell,cells=WhichCells(scell,idents = "CD-PC"))<-"Collecting duct"
  Idents(scell,cells=WhichCells(scell,idents = "CNT"))<-"Collecting duct"
  Idents(scell,cells=WhichCells(scell,idents = "DCT"))<-"Distal tubule"
  Idents(scell,cells=WhichCells(scell,idents = "tL"))<-"Proximal tubule"
  Idents(scell,cells=WhichCells(scell,idents = "TAL"))<-"Loop of Henle"
  Idents(scell,cells=WhichCells(scell,idents = "IntC"))<-"Interstitial cells"
  Idents(scell,cells=WhichCells(scell,idents = "Leuk"))<-"leukocytes"
  Idents(scell,cells=WhichCells(scell,idents = "EC"))<-"Epithelial cell"
  Idents(scell,cells=WhichCells(scell,idents = "PT"))<-"Proximal tubule"
  Idents(scell,cells=WhichCells(scell,idents = "Podo"))<-"Podocytes"
  scell[["Idents"]] <- scell@active.ident
  
  relate_gene <- read.table(file = "2023-05-04/network_predictions_gene.txt",quote = ",")
  colnames(relate_gene) <- "gene"
  GO <- enrichGO(relate_gene$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                                qvalueCutoff = 0.05,keyType = "SYMBOL")
  
  write.table(GO, file = paste(wd,"/",time,"/","GO.csv", sep = ""), sep=",")
  barplot(GO)
  genelist <- as.numeric(DEG_Differentiated[,2]) 
  names(genelist) <- row.names(DEG_Differentiated)  
  
  goselect <- read.csv(file = "2023-05-04/GO_select.csv")
  rownames(goselect) <- goselect$ID
  
  GO@result <- goselect  
  
# 20230519 GSEA分析，只放分化后的-----
  levels(KPC@active.ident)
  KPC_nephron <- subset(KPC,idents = c("CD-like","PT-like","LH-like","Progenitor cell", "Epithelial cell"))
  
  DimPlot(KPC_nephron, pt.size = 1)
  Idents(KPC_nephron,cells=WhichCells(KPC_nephron,idents = "CD-like"))<-"Tubular cell"
  Idents(KPC_nephron,cells=WhichCells(KPC_nephron,idents = "PT-like"))<-"Tubular cell"
  Idents(KPC_nephron,cells=WhichCells(KPC_nephron,idents = "LH-like"))<-"Tubular cell"
  Idents(KPC_nephron,cells=WhichCells(KPC_nephron,idents = "Epithelial cell"))<-"Tubular cell"
  
  Tubular.markers <- FindAllMarkers(KPC_nephron, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  
  geneList = Tubular.markers$avg_log2FC
  names(geneList) = Tubular.markers$gene
  head(geneList)
  geneList = sort(geneList, decreasing = TRUE)
  
  #GSEA分析
  set.seed(1)
  egmt<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
  #转换成数据框
  egmt_result_df <- as.data.frame(egmt)
  write.csv(egmt_result_df, file = "Tubule_like_GSAE.CSV")
  write.table(egmt_result_df,file="GSEA_MSigDb_h.all_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
  save(egmt,egmt_result_df,file = "GSEA_deg_h.all.rda")  
  
  gseaplot2(egmt,7,color="red",pvalue_table = T)
  
  #结果汇总
  gseaplot2(egmt, geneSetID = c(1,3), subplots = 1:3,pvalue_table = T)
  
  #气泡图 展示geneset被激活还是抑制
  egmt2<- setReadable(egmt,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
  dotplot(egmt2,split=".sign",showCategory = 10,font.size = 10, title = "",
          label_format = 30,)+facet_grid(~.sign)
  #edit legends
  # +guides(
  #reverse color order (higher value on top)
  #color = guide_colorbar(reverse = TRUE))
  #reverse size order (higher diameter on top) 
  #size = guide_legend(reverse = TRUE))
  # Title 可以添加标题
  
  #### 经典的GSEA图 
  egmt$ID[]
  a <- which(egmt_result_df$ID %in% c("GOBP_ORGANIC_ACID_METABOLIC_PROCESS"))
  egmt_result_df$ID["GOBP_ORGANIC_ACID_METABOLIC_PROCESS"]
  which(x = egmt_result_df$ID,"GOBP_ORGANIC_ACID_METABOLIC_PROCESS")
  egmt@result[["ID"]][1]
  exp1=dplyr::filter(egmt_result_df, grepl("GOBP_ORGANIC_ACID_METABOLIC_PROCESS", ID))
  which(dplyr::filter(egmt_result_df, grepl("GOBP_ORGANIC_ACID_METABOLIC_PROCESS", ID)))
  egmt$Description
  i <-  which(egmt_result_df$ID %in% c("GOBP_POSITIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION"))
  gseaplot2(egmt,
                      egmt$ID[i],#富集的ID编号
                      title = egmt$Description[i],#标题
                      color = "red", #GSEA线条颜色
                      base_size = 23,#基础字体大小
                      rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                      subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                      ES_geom = "line", #enrichment score用线还是用点"dot"
                      pvalue_table = T) #显示pvalue等信息
  ggsave(gseap1, filename = 'GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION.pdf', width =14, height =10)
  
  #### 合并 GSEA通路 
  gseap2 <- gseaplot2(kk_gse,
                      up_gsea$ID,#富集的ID编号
                      title = "UP_GSEA_all",#标题
                      color = "red",#GSEA线条颜色
                      base_size = 20,#基础字体大小
                      rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                      subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                      ES_geom = "line",#enrichment score用线还是用点"dot"
                      pvalue_table = T) #显示pvalue等信息
  ggsave(gseap2, filename = "GSEA_up_all.pdf",width =12,height =12)