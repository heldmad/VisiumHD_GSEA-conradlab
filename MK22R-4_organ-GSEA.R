# try subsetting whole slide into whole organs and comparing those!
# instead of comparing clusters within single organs
# eg lungs vs all other organs, etc
# use the same Loupe cluster labeling that you used previously, 
# just change clustering method to be whole organs instead of Loupe cluster numbers
# this file is for one slide; explain that this is to be replicated for each slide


# Install/load required packages

install.packages("tidyverse")
install.packages("readxl")
install.packages("writexl")
install.packages("utils")
install.packages("Seurat")
install.packages("Matrix")
install.packages("hdf5r")
install.packages("patchwork")
install.packages("future")
install.packages("sctransform")
install.packages("arrow")
install.packages("png")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EBImage")
BiocManager::install("clusterProfiler")
#BiocManager::install("celldex")
#BiocManager::install("SingleR")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
library(tidyverse)
library(readxl)
library(writexl)
library(utils)
library(Seurat)
library(Matrix)
library(hdf5r)
library(patchwork)
library(future)
library(sctransform)
library(arrow)
library(png)
library(EBImage)
library(clusterProfiler)
library(celldex)
library(SingleR)
library("org.Hs.eg.db", character.only = TRUE)
library(enrichplot)


# Loading in VisiumHD data (for Seurat purposes)
localdir <- '~/Desktop/spatial_analysis/MK22R-4_combined_out/008um/binned_outputs/square_008um/'
list.files(localdir)
MK22R4_combined <- Load10X_Spatial(data.dir = localdir, 
                                   filename = 'filtered_feature_bc_matrix.h5', 
                                   assay = "spatial", 
                                   slice = "slice1", 
                                   filter.matrix = TRUE)

# Set default assays to 'spatial'
Assays(MK22R4_combined)
DefaultAssay(MK22R4_combined) <- 'spatial'

# Visualize initial data
vln.plot <- VlnPlot(MK22R4_combined, features = "nCount_spatial", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(MK22R4_combined, features = "nCount_spatial", image.scale = "hires") + theme(legend.position = "right") #image.scale = "hires"
vln.plot | count.plot

# Update the default image to the hires one (for easier visualization in the future)
hires_image <- EBImage::readImage("~/Desktop/spatial_analysis/MK22R-4_combined_out/008um/spatial/tissue_hires_image.png")
# Assign it to the image slot in your Seurat object
DefaultAssay(MK22R4_combined) <- "spatial"
MK22R4_combined@images$slice1@image <- hires_image
MK22R4_combined@images$slice1@scale.factors

# Normalize the data
MK22R4 <- NormalizeData(MK22R4_combined)
MK22R4 <- FindVariableFeatures(MK22R4)
MK22R4 <- ScaleData(MK22R4)
MK22R4 <- subset(MK22R4, subset = nCount_spatial > 0) #remove NA to do SCT transform
