      ################## Requiered packages ##########################################
     
      library(dplyr)
      library(Seurat)
      library(cowplot)
      
          ########################## Exp_1 #################################################
        
          RO48_1 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/RO48_1/")%>%
          CreateSeuratObject(min.cells = 3, min.features = 200,project = "RO48 Exp.1")
        
        FeatureScatter(object = RO48_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
    
      DMSO_1 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/DMSO_1/")%>%
      CreateSeuratObject(min.cells = 3, min.features = 200,project = "DMSO Exp.1")
      
      FeatureScatter(object = DMSO_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
      
      mix_1 <- merge (x= DMSO_1,y=RO48_1)%>%
        NormalizeData(verbose = FALSE)%>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
      
      FeatureScatter(object = mix_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 1,smooth = T)
     
    rm(RO48_1,DMSO_1)#remove objects to reduce RAM load
   
  ########################## Exp_2 ##################################################
   
    RO48_2 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/RO48_2/")%>%
      CreateSeuratObject(min.cells = 3, min.features = 200,project = "RO48 Exp.2")
    
    FeatureScatter(object = RO48_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
    
    DMSO_2 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/DMSO_2/")%>%
      CreateSeuratObject(min.cells = 3, min.features = 200,project = "DMSO Exp.2")
    
    FeatureScatter(object = DMSO_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
    
    mix_2 <- merge (x= DMSO_2,y=RO48_2)%>%
      NormalizeData(verbose = FALSE)%>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
    
    FeatureScatter(object = mix_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 1,smooth = T)
      
    rm(RO48_2,DMSO_2)
    
   ######################### PCA exp_1 #################################################
  
    all.genes <- rownames(x = mix_1)
    
    mix_1 <- ScaleData(object = mix_1, features = all.genes)
    mix_1 <- RunPCA(object = mix_1, features = VariableFeatures(object = mix_1))

    VizDimLoadings(object = mix_1, dims = 1:2, reduction = 'pca')
    DimPlot(object = mix_1, reduction = 'pca')
    
    DimHeatmap(object =mix_1, dims = 1, cells = 500, balanced = TRUE)
    DimHeatmap(object =mix_1, dims = 1:15, cells = 500, balanced = TRUE)
    
    ElbowPlot(object = mix_1)
    
    ################### Clusters exp_1 ################################################
    mix_1 <- FindNeighbors(object = mix_1, dims = 1:10)
    mix_1 <- FindClusters(object = mix_1, resolution = 0.5)
    
    
    ######################### PCA exp_2 #################################################
    
    all.genes <- rownames(x = mix_2)
    
    mix_2 <- ScaleData(object = mix_1, features = all.genes)
    mix_2 <- RunPCA(object = mix_2, features = VariableFeatures(object = mix_2))
    
    VizDimLoadings(object = mix_2, dims = 1:2, reduction = 'pca')
    DimPlot(object = mix_2, reduction = 'pca')
    
    DimHeatmap(object =mix_2, dims = 1, cells = 500, balanced = TRUE)
    DimHeatmap(object =mix_2, dims = 1:15, cells = 500, balanced = TRUE)
    
    ElbowPlot(object = mix_2)
    
    ################### Clusters exp_2 ################################################
    mix_2 <- FindNeighbors(object = mix_2, dims = 1:10)
    mix_2 <- FindClusters(object = mix_2, resolution = 0.5)
    
   
    
