     ########################## Requiered packages ##########################################
     
      library(dplyr)
      library(Seurat)
      library(cowplot)
      
      ############################# Exp_1 #################################################
        
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
        
        ####################### UMAP exp_1 #######################################################
        #pip install umap-learn
      
      mix_1 <- RunUMAP(object = mix_1, dims = 1:10)
      DimPlot(object = mix_1, reduction = 'umap')
      
      ######################### cluster biomarkers exp_1 #######################################
      
      # find all markers of cluster 1
      cluster1.markers <- FindMarkers(object = mix_1, ident.1 = 1, min.pct = 0.25)
      head(x = cluster1.markers, n = 5)
      # find all markers distinguishing cluster 5 from clusters 0 and 3
      cluster5.markers <- FindMarkers(object = mix_1, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
      head(x = cluster5.markers, n = 5)
      # find markers for every cluster compared to all remaining cells, report only the positive ones
      pbmc.markers <- FindAllMarkers(object = mix_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
      
   ######################## Exp_2 ##################################################
   
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
    
    ####################### UMAP exp_2 #######################################################
    
    mix_2 <- RunUMAP(object = mix_2, dims = 1:10)
    DimPlot(object = mix_2, reduction = 'umap')
    
    ######################### cluster biomarkers exp_2 #######################################
    
    # find all markers of cluster 1
    cluster1.markers <- FindMarkers(object = mix_2, ident.1 = 1, min.pct = 0.25)
    head(x = cluster1.markers, n = 5)
    # find all markers distinguishing cluster 5 from clusters 0 and 3
    cluster5.markers <- FindMarkers(object = mix_2, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
    head(x = cluster5.markers, n = 5)
    # find markers for every cluster compared to all remaining cells, report only the positive ones
    pbmc.markers <- FindAllMarkers(object = mix_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
    
    ##################################################################################
    ###################### Full_data #################################################
    ##################################################################################
   
    rm(mix_1,mix_2)
    
    RO48_1 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/RO48_1/")%>%
      CreateSeuratObject(min.cells = 30, min.features = 2000,project = "RO48 Exp.1")
    DMSO_1 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/DMSO_1/")%>%
      CreateSeuratObject(min.cells = 30, min.features = 2000,project = "DMSO Exp.1")
    
    mix_1 <- merge(x=RO48_1,y= DMSO_1)%>%
      NormalizeData(verbose = FALSE)%>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
    
    rm(RO48_1,DMSO_1)
    
    RO48_2 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/RO48_2/")%>%
      CreateSeuratObject(min.cells = 30, min.features = 2000,project = "RO48 Exp.2")
    DMSO_2 <- Read10X(data.dir = "~/Dropbox/DataScience/Fiver/DMSO_2/")%>%
      CreateSeuratObject(min.cells = 30, min.features = 2000,project = "DMSO Exp.2")
    
    mix_2 <- merge(x=RO48_2,y= DMSO_2)%>%
      NormalizeData(verbose = FALSE)%>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
    
    rm(RO48_2,DMSO_2)
    
    mix_full <- merge(x=mix_1,y=mix_2)%>%
      NormalizeData(verbose = FALSE)%>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  
    rm (mix_1,mix_2)
    
    ##################################################################################
    
    FeatureScatter(object = mix_full, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 1,smooth = F)
    
    ######################### PCA exp_full #################################################
    
    all.genes <- rownames(x = mix_full)
    
    mix_full <- ScaleData(object = mix_full, features = all.genes)
    mix_full <- RunPCA(object = mix_full, features = VariableFeatures(object = mix_full))
    
    VizDimLoadings(object = mix_full, dims = 1:2, reduction = 'pca')
      DimPlot(object = mix_full, reduction = 'pca')
    
    DimHeatmap(object =mix_full, dims = 1, cells = 500, balanced = TRUE)
    DimHeatmap(object =mix_full, dims = 1:15, cells = 500, balanced = TRUE)
    
    ElbowPlot(object = mix_full)
    
    ################### Clusters exp_full ################################################
    
    mix_full <- FindNeighbors(object = mix_full, dims = 1:10)
    mix_full <- FindClusters(object = mix_full, resolution = 0.5)
    
    ####################### UMAP exp_full #######################################################
    #Need to install pip on Mac and run in console: sudo apt-get install pip 
    #Install the package with: pip install umap-learn
    
    mix_full <- RunUMAP(object = mix_full, dims = 1:10)
    DimPlot(object = mix_full, reduction = 'umap')
    
    ###########################################################################################
