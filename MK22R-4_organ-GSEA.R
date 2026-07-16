# Install/load required packages

install.packages("tidyverse")
install.packages("readxl")
install.packages("writexl")
install.packages("utils")
install.packages("png")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
library(tidyverse)
library(readxl)
library(writexl)
library(utils)
library(png)
library(clusterProfiler)
library("org.Hs.eg.db", character.only = TRUE)
library(enrichplot)


# Create an object with ALL genes expressed in the UNCLUSTERED LUNG (DEGs in lung only regardless of cluster)
# data from LOUPE BROWSER
all_genes_wholelung <- read.csv("~/Desktop/scratchwork/MK22R-4_lungs_whole_DEGs.csv", header = TRUE)
all_genes_wholelung <- all_genes_wholelung %>%
  dplyr::select(FeatureName, lungs_QCd.Log2.Fold.Change, lungs_QCd.P.Value)

# Upload NEW term2gene lists from G2P and ClinGen
clingen <- read_xlsx("~/Downloads/all_clingen_GCEPs_forR.xlsx")
clingen <- clingen %>%
  dplyr::select(GENE, TERM)
G2P <- read_xlsx("~/Downloads/G2P_withlabels_forR.xlsx")

# Combine into one data frame and make a TERM2GENE list with gene as entrez IDs
new_terms <- rbind(clingen, G2P)
new_TERM2GENE_entrez <- select(org.Hs.eg.db,
                               keys = new_terms$GENE,
                               columns = c("ENTREZID", "ENSEMBL"),
                               keytype = "SYMBOL")
new_terms <- new_terms %>%
  dplyr::select(TERM = TERM, SYMBOL = GENE)
new_TERM2GENE <- left_join(new_terms, new_TERM2GENE_entrez, by = "SYMBOL")
new_TERM2GENE <- new_TERM2GENE[!duplicated(new_TERM2GENE[,1:3]), ]
new_TERM2GENE <- new_TERM2GENE %>%
  dplyr::select(TERM = TERM, GENE = ENTREZID)
new_TERM2GENE <- new_TERM2GENE %>%
  filter(TERM != "NA")


## Whole Lung DGEs GSEA
# filtering out high p-value and low L2FC genes

#set seed
set.seed(1)

#create a vector for whole lung DEGs
filtered_wholelung <- all_genes_wholelung %>%
  filter(lungs_QCd.Log2.Fold.Change > 0.5 | lungs_QCd.Log2.Fold.Change < -0.5,
         lungs_QCd.P.Value < 0.25)
filtered_wholelung <- filtered_wholelung %>%
  dplyr::select(SYMBOL = FeatureName, 
                L2FC = lungs_QCd.Log2.Fold.Change, 
                PVal_Adj = lungs_QCd.P.Value)
filtered_wholelung_entrez <- select(org.Hs.eg.db,
                                    keys = filtered_wholelung$SYMBOL,
                                    columns = c("ENTREZID", "ENSEMBL"),
                                    keytype = "SYMBOL")
filtered_wholelunglist <- left_join(filtered_wholelung, filtered_wholelung_entrez, by = "SYMBOL")
filtered_wholelunglist <- filtered_wholelunglist[!duplicated(filtered_wholelunglist[,1]), ]
filtered_wholelunglist <- na.omit(filtered_wholelunglist)
allgenesList_wlf <- filtered_wholelunglist[,2]
names(allgenesList_wlf) <- as.character(filtered_wholelunglist[,4])
allgenesList_wlf <- sort(allgenesList_wlf, decreasing = TRUE)

wholelung_gsea <- GSEA(geneList = allgenesList_wlf,
                       TERM2GENE = new_TERM2GENE,
                       minGSSize = 1,
                       maxGSSize = 7000,
                       pvalueCutoff = 0.25,
                       seed = TRUE,
                       scoreType = "pos")
readable_wholelung_gsea <- setReadable(wholelung_gsea, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
df_wholelung_gsea <- readable_wholelung_gsea@result
gseaplot2(wholelung_gsea, df_wholelung_gsea[, 1])
cnetplot(readable_wholelung_gsea)


# Repeat the process above with all organs in the slide
# Download CSV with DGEs for that organ and read it into this file


# Heart as another example
all_genes_wholeheart <- read.csv("~/Desktop/MK22R-4_heartvsall.csv", header = TRUE)
all_genes_wholeheart <- all_genes_wholeheart %>%
  dplyr::select(FeatureName, heart_QCd.Log2.Fold.Change, heart_QCd.P.Value)

# filtering out high p-value and low L2FC genes
#create a vector for whole lung DEGs
filtered_wholeheart <- all_genes_wholeheart %>%
  filter(heart_QCd.Log2.Fold.Change > 0.5 | heart_QCd.Log2.Fold.Change < -0.5,
         heart_QCd.P.Value < 0.25)
filtered_wholeheart <- filtered_wholeheart %>%
  dplyr::select(SYMBOL = FeatureName, 
                L2FC = heart_QCd.Log2.Fold.Change, 
                PVal_Adj = heart_QCd.P.Value)
filtered_wholeheart_entrez <- select(org.Hs.eg.db,
                                    keys = filtered_wholeheart$SYMBOL,
                                    columns = c("ENTREZID", "ENSEMBL"),
                                    keytype = "SYMBOL")
filtered_wholeheartlist <- left_join(filtered_wholeheart, filtered_wholeheart_entrez, by = "SYMBOL")
filtered_wholeheartlist <- filtered_wholeheartlist[!duplicated(filtered_wholeheartlist[,1]), ]
filtered_wholeheartlist <- na.omit(filtered_wholeheartlist)
allgenesList_whf <- filtered_wholeheartlist[,2]
names(allgenesList_whf) <- as.character(filtered_wholeheartlist[,4])
allgenesList_whf <- sort(allgenesList_whf, decreasing = TRUE)

wholeheart_gsea <- GSEA(geneList = allgenesList_whf,
                       TERM2GENE = new_TERM2GENE,
                       minGSSize = 1,
                       maxGSSize = 7000,
                       pvalueCutoff = 0.25,
                       seed = TRUE,
                       scoreType = "pos")
readable_wholeheart_gsea <- setReadable(wholeheart_gsea, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
df_wholeheart_gsea <- readable_wholeheart_gsea@result
gseaplot2(wholeheart_gsea, df_wholeheart_gsea[, 1])
cnetplot(readable_wholeheart_gsea)
