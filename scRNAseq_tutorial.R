# Tutorial Single Cell 
#Autor: Nathalia Portilla (2024)

#A Set up
#1.Install and load libraries 
#srun --pty --partition=medium --mem=20GB --time=04:00:00 bash in bash terminal in case of using interative session

install.packages("tidyverse")
install.packages("pak")
install.packages("BiocManager")
BiocManager::install("hdf5r")
library(pak)
pak::pkg_install("satijalab/seurat")
#Note: maybe you need to install the hdf5r library so update using sudo apt install libhdf5-dev on terminal
install.packages("/hpcfs/home/ciencias_biologicas/na.portilla10/journal/single_cell/hdf5r_1.3.11.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(Seurat)
library(hdf5r)

#2.Read dataset,using Read function reads in the output of the cellranger pipeline from 10X, 
#returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number 
#of molecules for each feature (i.e. gene; row) that are detected in each cell (column). Note that more recent 
#versions of cellranger now also output using the h5 file format,
h5_data<-Read10X_h5(filename ="/hpcfs/home/ciencias_biologicas/na.portilla10/journal/single_cell/NSCLC.h5")
str(h5_data)
#For Gene expression analysis 
h5_genes<-h5_data$`Gene Expression`
h5_genes
#3. Initialize the Seurat object with the raw (non-normalized data).
h5<-CreateSeuratObject(counts = h5_genes,assay = "RNA",project="H5",min.cells = 3, min.features = 200)

#B QC Metrics 
#1 Map mithocondrial genes on each cell
h5[["percent.mt"]] <- PercentageFeatureSet(h5, pattern = "^MT-")
#2 Viz to QC metrics
# Visualize QC metrics as a violin plot
violin<-VlnPlot(h5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("VlnPlot.png", plot = violin, width = 10, height = 5, dpi = 300)
#Feature scatter 
CountVsFeauture <- FeatureScatter(h5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm') +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )
ggsave("CountVSFeature.png", plot = CountVsFeauture, width = 10, height = 5, dpi = 300)

#C Filtering 
h5 <- subset(h5, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#D Normalize data 
#normalizes the feature expression measurements for each cell by the total expression, 
#multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
#In Seurat v5, Normalized values are stored in 
h5<-pbmc <- NormalizeData(h5, normalization.method = "LogNormalize", scale.factor = 10000)
str(h5)
#E Feature Selection
#1. features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).
h5 <- FindVariableFeatures(h5, selection.method = "vst", nfeatures = 2000)
# 2.Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(h5), 10)
# 3. plot variable features with and without labels
plot1 <- VariableFeaturePlot(h5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )
ggsave("High_Features.png", plot = plot2, width = 10, height = 5, dpi = 300)

#F SCALE DATA 
#Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. 
all.genes <- rownames(h5)
h5 <- ScaleData(h5, features = all.genes)
str(h5)

#G Lineal dimensionality reduction
h5 <-RunPCA(h5, features = VariableFeatures(object = h5))
print(h5[["pca"]], dims = 1:5, nfeatures = 5)
pca<-DimPlot(h5, reduction = "pca") + NoLegend()
ggsave("pca.png", plot = pca)
#heatmap<-DimHeatmap(h5, dims = 1:15, cells = 500, balanced = TRUE)
#ggsave("Heatmap_hetero.png", plot = heatmap, width = 10, height = 5, dpi = 300)

#H Determinate dimensionality 
elbowz<-ElbowPlot(h5)+
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )
ggsave("elbow.png", plot = elbowz)

#I Clustering 
h5<- FindNeighbors(h5, dims = 1:15)
h5<- FindClusters(h5, resolution = c(0.1,0.3, 0.5, 0.7, 1))

clusters<-DimPlot(h5, group.by = "RNA_snn_res.0.5", label = TRUE)
ggsave("clusters.png", plot = clusters)

# setting identity of clusters
Idents(h5)
Idents(h5) <- "RNA_snn_res.0.1"

#J NON LINEAL DIMESIONALITY REDUCTION 
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
h5 <- RunUMAP(h5, dims = 1:8)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
umap<-DimPlot(h5, reduction = "umap")
ggsave("umap.png", plot = umap)
saveRDS(h5, file = "h5_tutorial.rds")

#K Biomarkers
# find all markers of cluster 2
cluster2.markers <- FindMarkers(h5, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(h5, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
h5.markers <- FindAllMarkers(h5, only.pos = TRUE)
h5.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#plot
  h5.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
h5.markers.heatmap<-DoHeatmap(h5, features = top10$gene) 
ggsave("gene_expression_heatmap.png", plot = h5.markers.heatmap)
