## Creating Fibroblast species object (Fig7) with our Fibroblast scRNA-Seq data (7/22/82 wo ChP 4V&LV) and human snRNA-seq data (Yang et al.)
## Script attempted integration of the datasets via Seurat (CCA)  and subsetted to common homologs
## CCA integration was finally not used, but this object was used as starting point for the BBKNN workflow

## Load Packages
library(Seurat)
library(ggplot2)
library(gridExtra)

## Set parameters
# setwd("G:/VIB_G_drive/")
setwd("/run/media/clintdn/CN1465-DATA/VIB_G_drive") #Linux

sampleName_Human<-"Human_snRNAseq_CP_COVID"
sampleFolder<-paste0(sampleName_Human,"/")
sampleName <- "CCA_integrated_Human_and_Mouse_CP_W2"

## Read human dataset
seuratObj_Human_subset_converted<-readRDS(file = paste0(sampleFolder,"Robjects/seuratObj_Human_subset_converted.rds"))
Original_human_gene_symbols<-readRDS(file=paste0(sampleFolder,"Robjects/Original_human_gene_symbols.rds"))
M.genes<-readRDS(file=paste0(sampleFolder,"Robjects/Converted_human_to_mouse_gene_symbols.rds"))
M.genes_df<-readRDS(file=paste0(sampleFolder,"Robjects/Converted_human_to_mouse_gene_symbols_df.rds"))

##### Read mouse object Urvb (scRNAseq)
seuratObj_Mouse <- readRDS(file=paste0("/home/clintdn/VIB/DATA/Roos/Daan 1/FB_datasets/Merge_Urvb_extra/Robjects/seuratObj_subset_Urvb_FB_6datasets_harmony_RNA.rds"))

seuratObj_Mouse <- UpdateSeuratObject(seuratObj_Mouse)
DimPlot(seuratObj_Mouse, reduction = "umap", label=T,repel = T, group.by="annotated_clusters_new", pt.size = 1)

## Defaultassay
DefaultAssay(seuratObj_Mouse)
DefaultAssay(seuratObj_Human_subset_converted)


############################################################################################################
############################################################################################################

#######################
###### Workflow 2 ##### 
#######################

## First set the genes to the same for both seuratobjects -> increase overlap?
## Drawback is it will remove certain genes which would be in RNA slot

##### Extract count matrix and metadata
Counts_human <- GetAssayData(seuratObj_Human_subset_converted[["RNA"]], slot = "counts")
Metadata_human <- seuratObj_Human_subset_converted@meta.data

Counts_mouse <- GetAssayData(seuratObj_Mouse[["RNA"]], slot = "counts")
Metadata_mouse <- seuratObj_Mouse@meta.data

Common_genes<-intersect(rownames(Counts_human),rownames(Counts_mouse))

Counts_human<-Counts_human[Common_genes,] #11367
Counts_mouse<-Counts_mouse[Common_genes,] #11367

## Recreate seuratObjects
seuratObj_Human_final <- CreateSeuratObject(counts = Counts_human, project = "Human_CP_subset", assay = "RNA")
seuratObj_Mouse_final <- CreateSeuratObject(counts = Counts_mouse, project = "Mouse_CP_subset", assay = "RNA")

## Add metaData back 
seuratObj_Human_final@meta.data <-cbind(seuratObj_Human_final@meta.data, Metadata_human)
seuratObj_Mouse_final@meta.data <-cbind(seuratObj_Mouse_final@meta.data, Metadata_mouse)

# Remove duplicate columns (slight diff due to removal certain genes. Remove the old ones!)
seuratObj_Human_final@meta.data[c(4,5,6)]<-NULL
seuratObj_Mouse_final@meta.data[c(1,5,6)]<-NULL #Keep orig ident original

## Normalize
seuratObj_Human_final <- NormalizeData(object = seuratObj_Human_final, normalization.method = "LogNormalize", scale.factor = 10000)
seuratObj_Mouse_final <- NormalizeData(object = seuratObj_Mouse_final, normalization.method = "LogNormalize", scale.factor = 10000)

## Find Variable Features first for Human data
seuratObj_Human_final <- FindVariableFeatures(seuratObj_Human_final, assay = "RNA", selection.method = "vst", nfeatures = 2000)
seuratObj_Mouse_final <- FindVariableFeatures(seuratObj_Mouse_final, assay = "RNA", selection.method = "vst", nfeatures = 2000)

# select features that are repeatedly variable across datasets for integration
FB_list_2 <- list(seuratObj_Mouse_final, seuratObj_Human_final)
features_2 <- SelectIntegrationFeatures(object.list = FB_list_2)

## Integration
# We then identify anchors using the FindIntegrationAnchors() function, 
# which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
FB.anchors_2 <- FindIntegrationAnchors(object.list = FB_list_2, anchor.features = features_2, dims = 1:30)

FB.combined_2 <- IntegrateData(anchorset = FB.anchors_2)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(FB.combined_2) <- "integrated"

# Run the standard workflow for visualization and clustering
FB.combined_2 <- ScaleData(FB.combined_2, verbose = FALSE)
FB.combined_2 <- RunPCA(FB.combined_2, npcs = 30, verbose = FALSE)
FB.combined_2 <- RunUMAP(FB.combined_2, reduction = "pca", dims = 1:30)
FB.combined_2 <- FindNeighbors(FB.combined_2, reduction = "pca", dims = 1:30)
FB.combined_2 <- FindClusters(FB.combined_2, resolution = 0.5)

## Add metadata column (Species)
FB.combined_2@meta.data$Species<-"Human"
FB.combined_2@meta.data[which(FB.combined_2@meta.data$Lab == "Vandenbroucke"),"Species"] <- "Mouse"

# Visualization
p3 <- DimPlot(FB.combined_2, reduction = "umap", group.by = "Species") + ggtitle("Species_workflow2")
p4 <- DimPlot(FB.combined_2, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Numbered_annotation_workflow2")
p3+p4

FeaturePlot(FB.combined_2, features = "Igfbp6", min.cutoff = "q10", max.cutoff = "q90")
FB.combined_2$NewClusters_combo<-paste0(FB.combined_2@active.ident,"_",FB.combined_2$Species)
table(FB.combined_2$NewClusters_combo)

####################################################################################################################

########## Via PCelbowplot ##########
ElbowPlot(object = FB.combined_2, ndims = 30)

## Perform clustering
dimsToTry<-c(seq(10,30,by=5))

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  FB.combined_2 <- FindNeighbors(object = FB.combined_2, reduction = "pca", dims = dimsToUse)
  FB.combined_2 <- FindClusters(object = FB.combined_2,  graph.name = "integrated_snn", resolution = resToUse)

  ##### Create UMAP plot
  FB.combined_2 <- RunUMAP(FB.combined_2, dims = dimsToUse, n_neighbors = 30, assay = "integrated", reduction ="pca",
                           reduction.name = "umap", reduction.key = "integratedUMAP_")
  umapPlot<-DimPlot(FB.combined_2, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(FB.combined_2, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(sampleFolder,"CCA_FB_combined_2/UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

### Final
dimsToTry<-c(20)
resToUse<-0.8
diagnostics<-list()
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  FB.combined_2 <- FindNeighbors(object = FB.combined_2, reduction = "pca", dims = dimsToUse)
  FB.combined_2 <- FindClusters(object = FB.combined_2,  graph.name = "integrated_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  FB.combined_2 <- RunUMAP(FB.combined_2, dims = dimsToUse, n_neighbors = 30, assay = "integrated", reduction ="pca",
                           reduction.name = "umap", reduction.key = "integratedUMAP_")
  umapPlot<-DimPlot(FB.combined_2, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(FB.combined_2, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(sampleFolder,"CCA_FB_combined_2/UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

names(FB.combined_2)

### Clustering: trying out clusTree
Perplexity<-20
Resolution<-0.8

FB.combined_2 <- FindNeighbors(object = FB.combined_2, reduction = "pca", dims = 1:Perplexity)
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  FB.combined_2 <- FindClusters(object = FB.combined_2,  graph.name = "integrated_snn", resolution = res)
}

library("clustree")
pdf(file=paste0(sampleFolder,"CCA_FB_combined_2/Clustree.pdf"))
clustree(FB.combined_2, prefix = "integrated_snn_res.")
dev.off()

# 0.8 seems reasonable
# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
FB.combined_2$Integrated_RNA_clusters <- FB.combined_2$integrated_snn_res.0.8
Idents(FB.combined_2)<-FB.combined_2@meta.data$Integrated_RNA_clusters

############################################################################################

##### Save object
saveRDS(FB.combined_2, file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))

##### Read object
FB.combined_2<- readRDS(file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))
