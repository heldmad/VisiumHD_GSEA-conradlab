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