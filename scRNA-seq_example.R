library(R.matlab)
library(Seurat)
library(ggplot2)

#Import Data and format for Seurat 
mtx <- readMat("/data/singlecell.mat")
dge <- mtx$DGE
genes <- as.vector(mtx$genes)
genes <- trimws(genes, "r")
genes <- gsub('.{6}$', '', genes)
samples <- mtx$sample.type
colnames(dge) <- genes
rownames(dge) <- make.names(samples, unique = T)
dge_t <- t(dge)

#Create Seurat Object and Add MetaData
cbl_dev <- CreateSeuratObject(counts = dge_t, project = "80k")
cbl_dev$samples <- samples

#Filter outlier cells - Remove 4SD above median for nGenes (nFeature_RNA) & nUMIs (nCount_RNA)
nGenes_outlier <- (4*sd(cbl_dev@meta.data$nFeature_RNA)) + median(cbl_dev@meta.data$nFeature_RNA)
nUMIs_outlier <- (4*sd(cbl_dev@meta.data$nCount_RNA)) + median(cbl_dev@meta.data$nCount_RNA)
cbl_dev <- subset(x = cbl_dev, subset = nFeature_RNA > 200 & nFeature_RNA < nGenes_outlier & nCount_RNA < nUMIs_outlier & percent.mt < 1)
VlnPlot(object = cbl_dev, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", pt.size = 0)

#Clustering Pipeline
#Use either:
cbl_dev <- SCTransform(cbl_dev, vars.to.regress = "CC.Difference") #Regularized negative binomial regression to normalize UMI count data
#or: 
cbl_dev <- NormalizeData(cbl_dev, normalization.method = 'LogNormalize', scale.factor = 10000) #Natural-log tranformed using log1p counts
cbl_dev <- FindVariableFeatures(cbl_dev, selection.method = 'vst') #Determines most variable features using local polynomial regression
cbl_dev <- ScaleData(cbl_dev, vars.to.regress = 'CC.Difference') #Scales and centers dataset, regresses out cell cycling
#Rest of pipeline
cbl_dev <- RunPCA(cbl_dev, npcs = 50) #Principal Component Analysis
ElbowPlot(cbl_dev, ndims = 50)
cbl_dev <- RunUMAP(cbl_dev, dims = 1:50) #Uniform Manifold Approximation and Projection
cbl_dev <- FindNeighbors(cbl_dev, dims = 1:50) #Shared Nearest Neighbor Graph
cbl_dev <- FindClusters(cbl_dev, resolution = 1) #Louvain Clustering
DimPlot(cbl_dev, reduction = 'umap', label = T)

#Differential Gene Test - Wilcoxon Rank Sum test to determine differentially expressed genes for each cluster
markers <- FindAllMarkers(cbl_dev, only.pos = T, test.use = 'wilcox')  
