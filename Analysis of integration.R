# iPSC肾类器官重新定义的单细胞数据集（使用的是E6的数据）-----
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
# iPSC与KPC类器官数据进行整合-----
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