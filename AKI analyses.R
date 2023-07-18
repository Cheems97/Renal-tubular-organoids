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
