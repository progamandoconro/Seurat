################## Requiered packages ##########################################

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)

############################# Exp_1 #################################################

##Read, normalize and find features in the data

RO48_1 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/RO48_1/")%>%
  CreateSeuratObject(min.cells = 3, min.features = 200,project = "RO48 Exp.1")%>%
  NormalizeData(verbose = FALSE)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

DMSO_1 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/DMSO_1/")%>%
  CreateSeuratObject(min.cells = 3, min.features = 200,project = "DMSO Exp.1")%>%
  NormalizeData(verbose = FALSE)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

# Create variables for grouping later

RO48_1$RO48_1 <- 'RO48_1'
DMSO_1$DMSO_1 <- 'DMSO_1'

# Integrate data

exp_1 <- FindIntegrationAnchors(object.list = list(DMSO_1,RO48_1), dims = 1:20)%>%
IntegrateData(dims = 1:20)

#Remove objects to save RAM

rm(DMSO_1,RO48_1)
DefaultAssay(exp_1) <- "integrated"

# PCA, t-SNE and Cluster in one pipe

exp_1 <-ScaleData(exp_1,verbose = FALSE)%>%
RunPCA(npcs = 30, verbose = FALSE)%>%
RunUMAP(reduction = "pca", dims = 1:20)%>%
FindNeighbors( reduction = "pca", dims = 1:20)%>%
FindClusters( resolution = 0.5)

# Visualization

g1 <- DimPlot(exp_1, reduction = "umap",group.by = 'DMSO_1')
    
g1 + scale_color_discrete(labels=c('DMSO_1','RO48_1'))
 
DimPlot(exp_1, reduction = "umap", label = TRUE)

rm(exp_1)
################################### Exp_2 ########################################

##Read, normalize and find features in the data

RO48_2 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/RO48_2/")%>%
  CreateSeuratObject(min.cells = 30, min.features = 2000,project = "RO48 Exp.2")%>%
  NormalizeData(verbose = FALSE)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)


DMSO_2 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/DMSO_2/")%>%
  CreateSeuratObject(min.cells = 30, min.features = 2000,project = "DMSO Exp.2")%>%
  NormalizeData(verbose = FALSE)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

# Create variables for grouping later

RO48_2$RO48_2 <- 'RO48_2'
DMSO_2$DMSO_21 <- 'DMSO_2'

# Integrate data

exp_2 <- FindIntegrationAnchors(object.list = list(DMSO_2,RO48_2), dims = 1:20)%>%
  IntegrateData(dims = 1:20)

#Remove objects to save RAM

rm(DMSO_2,RO48_2)
DefaultAssay(exp_2) <- "integrated"

# PCA, t-SNE and Cluster in one pipe

exp_2 <-ScaleData(exp_2,verbose = FALSE)%>%
  RunPCA(npcs = 30, verbose = FALSE)%>%
  RunUMAP(reduction = "pca", dims = 1:20)%>%
  FindNeighbors( reduction = "pca", dims = 1:20)%>%
  FindClusters( resolution = 0.5)

# Visualization

g2 <- DimPlot(exp_2, reduction = "umap",group.by = 'DMSO_2')

g2 + scale_color_discrete(labels=c('DMSO_2','RO48_2'))

DimPlot(exp_2, reduction = "umap", label = TRUE)

rm(exp_2)
################################# Exp1 / Exp2 ###############################

##Read, normalize and find features in the data
RO48_1 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/RO48_1/")%>%
  CreateSeuratObject(min.cells = 3, min.features = 200,project = "RO48 Exp.1")%>%
  NormalizeData(verbose = FALSE)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

DMSO_1 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/DMSO_1/")%>%
  CreateSeuratObject(min.cells = 3, min.features = 200,project = "DMSO Exp.1")%>%
  NormalizeData(verbose = FALSE)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

RO48_2 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/RO48_2/")%>%
  CreateSeuratObject(min.cells = 3, min.features = 200,project = "RO48 Exp.2")%>%
  NormalizeData(verbose = FALSE)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

DMSO_2 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/DMSO_2/")%>%
  CreateSeuratObject(min.cells = 3, min.features = 200,project = "DMSO Exp.2")%>%
  NormalizeData(verbose = FALSE)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

# Create variables for grouping later

RO48_1$RO48_1 <- 'RO48_1'
DMSO_1$DMSO_1 <- 'DMSO_1'
RO48_2$RO48_2 <- 'RO48_2'
DMSO_2$DMSO_2 <- 'DMSO_2'

# Integrate data

exp_full <- FindIntegrationAnchors(object.list = list(DMSO_1,RO48_1,DMSO_2,RO48_2), dims = 1:20)%>%
  IntegrateData(dims = 1:20)

DefaultAssay(exp_full) <- "integrated"

#Remove objects to save RAM

rm(DMSO_1,RO48_1,DMSO_2,RO48_2)

# PCA, t-SNE and Cluster in one pipe

exp_full <-ScaleData(exp_full,verbose = FALSE)%>%
  RunPCA(npcs = 30, verbose = FALSE)%>%
  RunUMAP(reduction = "pca", dims = 1:20)%>%
  FindNeighbors( reduction = "pca", dims = 1:20)%>%
  FindClusters( resolution = 0.5)

# Visualization

g2 <- DimPlot(exp_full, reduction = "umap",group.by = 'DMSO_1')

g2 + scale_color_discrete(labels=c('DMSO_1','RO48_1,DMSO_2','RO48_2'))

DimPlot(exp_full, reduction = "umap", label = TRUE)

####################################################################################################################
