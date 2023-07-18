library(Seurat)
library(dplyr)
library(clusterProfiler)
library(monocle3) 
library(reticulate)
library(ggplot2)

# Seurat was used for quality control----
  KPC.data <- Read10X(data.dir = "D:/Zuo Lab/Data/20210325 sequence Data/Result-X101SC19090611-Z01-F001-37-B1/2.Summary/KPC-org/filtered_feature_bc_matrix")
  organoid <- CreateSeuratObject(counts = KPC.data, project = "KPC_organoid", min.cells = 3, min.features = 200)
  organoid[["percent.mt"]] <- PercentageFeatureSet(organoid, pattern = "^MT-")
  organoid <- subset(organoid, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  VlnPlot(organoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  organoid <- NormalizeData(organoid, normalization.method = "LogNormalize", scale.factor = 10000)
  organoid <- FindVariableFeatures(organoid, selection.method = "vst", nfeatures = 2000)
  
  rm(KPC.data)
  
  
# Seurat was used to generate S4 files for building cds----
  organoid <- ScaleData(organoid, vars.to.regress = c("nCount_RNA", "percent.mito"))
  organoid <- RunPCA(organoid, npcs = 50)
  organoid <- ProjectDim(object = KPC)
  
  organoid <- RunUMAP(organoid, reduction = "pca", dims = 1:50)
  organoid <- FindNeighbors(organoid, reduction = "pca", dims = 1:50)
  organoid <- FindClusters(organoid, resolution = 0.2)


# Create cds files-----
  data <- as(as.matrix(organoid@assays$RNA@counts), 'sparseMatrix')
  pd <-  organoid@meta.data
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  colnames(pd)
  
  cds <- new_cell_data_set(data,
                           cell_metadata  = pd,
                           gene_metadata  = fData)
  
 # Pre-process the data
  cds = preprocess_cds(cds, method = c("PCA"),num_dim = 100)
  plot_pc_variance_explained(cds)

 # Reduce dimensionality and visualize the cells
  cds = reduce_dimension(cds) #Monocle uses UMAP by default 
  plot_cells(cds) 
  
  cds = cluster_cells(cds, python_home = "E:\\Work\\Python\\python.exe",verbose = T,resolution = 0.00045)
  plot_cells(cds,graph_label_size=0.5,cell_size=0.5,label_cell_groups=FALSE)
  
# The single cell data information after monocle3 analysis was transferred to seurat's single cell file for presentation-----
  KPC@reductions$umap@cell.embeddings <- cds@int_colData$reducedDims$UMAP
  colnames(KPC@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  KPC@meta.data[["Monocle.Cluster"]] <- cds@clusters$UMAP$clusters
  KPC@active.ident <- KPC$Monocle.Cluster
  DimPlot(KPC, pt.size = 1)
  
  Idents(KPC,cells=WhichCells(KPC,idents = "1"))<-"CD-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "2"))<-"PT-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "3"))<-"DTLH-like"
  Idents(KPC,cells=WhichCells(KPC,idents = "4"))<-"Progenitor cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "5"))<-"Epithelial cell"
  Idents(KPC,cells=WhichCells(KPC,idents = "6"))<-"Cycling cell"
  

