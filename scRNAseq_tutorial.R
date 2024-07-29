# Tutorial Single Cell 
#Autor: Nathalia Portilla (2024)

#A Set up
#1.Install and load libraries 
install.packages("tidyverse")
install.packages("pak")
install.packages("BiocManager")
pak::pkg_install("satijalab/seurat")
#Note: maybe you need to install the hdf5r library so update using sudo apt install libhdf5-dev on terminal
install.packages("/home/nathalia/Documentos/Journal/Single_cell/hdf5r_1.3.11.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(Seurat)
library(pak)
library(hdf5r)

#2.Read dataset,using Read function reads in the output of the cellranger pipeline from 10X, 
#returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number 
#of molecules for each feature (i.e. gene; row) that are detected in each cell (column). Note that more recent 
#versions of cellranger now also output using the h5 file format,
h5_data<-Read10X_h5(filename ="/home/nathalia/Documentos/Journal/Single_cell/PBMC_matrix.h5")
summary(h5_data)
#For Gene expression analysis 
h5_genes<-h5_data$`Gene Expression`
#3. Initialize the Seurat object with the raw (non-normalized data).
h5<-CreateSeuratObject(counts = h5_genes,assay = "RNA",project="H5",min.cells = 3, min.features = 200)

#B QC Metrics 
#1 Map mithocondrial genes on each cell
h5[["percent.mt"]] <- PercentageFeatureSet(h5, pattern = "^MT-")
#2 Viz to QC metrics
# Visualize QC metrics as a violin plot
VlnPlot(h5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Feature scatter 
CountVsFeauture <- FeatureScatter(h5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ geom_smooth(method = 'lm')
CountVsFeauture

