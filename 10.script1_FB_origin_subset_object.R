## Creating Fibroblast origin subset object with only the Fibroblast clusters from the Fibroblast origin complete object
## Workflow: don't take along prolif FBs or PVM or vascular cells and only take cells from clusters 1,2,5,8,14,20
## Only took along cells from those clusters which mapped on the left side of the UMAP (performed manual tracing of the clusters)
## Initial basic Seurat and CCA workflow script to create the Fibroblast origin subset object in Figure 1 manuscript
## Rebuttal: Performed CCA instead of harmony and now includes the new Betsholtz lab data

library(Seurat)
library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(gridExtra)
library(clustree)

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")

source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

################################################################################
################################################################################

sampleName <- "Rebuttal7_Full_Merge_subset_CCA_v7_FB_datasets" 
sampleFolder<-paste0("Rebuttal_mouse_data/Merge_subset_Final","/")
output.dir<-"Rebuttal_mouse_data/Merge_subset_Final/"

## Prepare folders
dir.create("Rebuttal_mouse_data/Merge_subset_Final/")
dir.create("Rebuttal_mouse_data/Merge_subset_Final/Plots")
dir.create("Rebuttal_mouse_data/Merge_subset_Final/Plots/RNA")
dir.create("Rebuttal_mouse_data/Merge_subset_Final/Robjects")

### Read RDS object
seuratObjAll <- readRDS(file="Rebuttal_mouse_data/Full_merge_Final/Robjects/seuratObj_Rebuttal4_Full_merge_Final_FB_datasets_harmony_RNA.rds")

### Subset data ### (non-prolif FBs)
## Cluster data combined with manual subselection (to remove outlier cells of those clusters!)
clusterMatrix<-seuratObjAll@meta.data
umapTable<-as.data.frame(seuratObjAll[['umap']]@cell.embeddings, stringsAsFactors = F)

Idents(seuratObjAll)<-seuratObjAll@meta.data$harmony_clusters
DimPlot(seuratObjAll, label = T)

U1 <- DimPlot(seuratObjAll, reduction = "umap", label = T, label.size = 4)
seuratObjAll <- CellSelector(U1, object=seuratObjAll, ident="Selected")
p1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObjAll, idents = c("Selected"))) #11639
wantedCells<-intersect(rownames(seuratObjAll@meta.data[which(seuratObjAll$harmony_clusters %in% c(1,2,5,8,14,20)),]), WhichCells(seuratObjAll, idents = c("Selected")))
p2<-colorSomeCells(clusterMatrix, umapTable, wantedCells) #11312 

pdf(file=paste0(output.dir,"/Plots/RNA/16a_Cell_selection_subset_CCA.pdf"), width = 16, height = 10)
print(p1)
print(p2)
dev.off()

seuratObj<-subset(seuratObjAll, cells = wantedCells)

dim(seuratObj)
# 9651

diagnostics<-list()

unique(sapply(X = strsplit(colnames(seuratObj), split = "_"), FUN = "[", 1))
# "FB1"   "EP1"   "EP2"   "EP3"   "VASC1" "VASC2" "EP4"   "URVB" "FB2"   "FB3"  "FB4"  

########## Get HVG ##########
seuratObj <- FindVariableFeatures(object = seuratObj, assay = "RNA", selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(seuratObj, assay = "RNA"))
# 2000

### Add to diagnostics
diagnostics[['varGenes']]<-length(VariableFeatures(seuratObj, assay = "RNA"))

########## Scale ##########
seuratObj <- ScaleData(seuratObj, assay = "RNA")

# Run PCA on rna normalized through scran/scater
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj, assay = "RNA"), 
                    npcs = 150, ndims.print = 1:5, nfeatures.print = 10, assay = "RNA")

## Use RNA going forward 
DefaultAssay(object = seuratObj)<-"RNA"

## Split up object for CCA
## Remark!!
## seuratObj2 and seuratObj5 do not contain enough cells! From ependymal datasets. Not able to combine with other dataset. Left out of subset!
## seuratObj3 also causing errors, but important to keep. Also 10X from Betsholtz lab. Combine with seuratObj7
## seuratObj8 also 10X from Betsholtz lab. Combine with seuratObj7

c("FB1","EP1","EP2","EP3","VASC1","VASC2","EP4","URVB","FB2","FB3","FB4")

seuratObj1<-subset(seuratObj, cells = colnames(seuratObj)[grep("FB1",colnames(seuratObj))])
# seuratObj2<-subset(seuratObj, cells = colnames(seuratObj)[c(grep("EP1",colnames(seuratObj)),
#                                                             grep("EP2",colnames(seuratObj)),
#                                                             grep("EP3",colnames(seuratObj)))])
# seuratObj3<-subset(seuratObj, cells = colnames(seuratObj)[grep("VASC1",colnames(seuratObj))]) #Not enough cells BFB1/BFB6 -> put with object 7 and 8. ALso 10X
seuratObj4<-subset(seuratObj, cells = colnames(seuratObj)[grep("VASC2",colnames(seuratObj))])
# seuratObj5<-subset(seuratObj, cells = colnames(seuratObj)[grep("EP4",colnames(seuratObj))])
seuratObj6<-subset(seuratObj, cells = colnames(seuratObj)[grep("URVB",colnames(seuratObj))])
seuratObj7<-subset(seuratObj, cells = colnames(seuratObj)[c(grep("FB2",colnames(seuratObj)),grep("FB3",colnames(seuratObj)),grep("VASC1",colnames(seuratObj)))])
# seuratObj8<-subset(seuratObj, cells = colnames(seuratObj)[grep("FB3",colnames(seuratObj))]) #Put with object 8. Also 10X
seuratObj9<-subset(seuratObj, cells = colnames(seuratObj)[grep("FB4",colnames(seuratObj))])

# select features that are repeatedly variable across datasets for integration
FB_list_2 <- list(seuratObj6,seuratObj1,seuratObj4,seuratObj7, seuratObj9) #seuratObj2,seuratObj5,seuratObj3,seuratObj8
features_2 <- SelectIntegrationFeatures(object.list = FB_list_2)

## Run loop to try out different parameters
dimsToTryCCA<-c(seq(15,30,by=5))

for(maxCCA in dimsToTryCCA){
  
  ## Integration
  # We then identify anchors using the FindIntegrationAnchors() function, 
  # which takes a list of Seurat objects as input, and use these anchors to integrate the four datasets together with IntegrateData().
  FB.anchors_2 <- FindIntegrationAnchors(object.list = FB_list_2, anchor.features = features_2, dims = 1:maxCCA, normalization.method = "LogNormalize", reduction = "cca")
  
  FB.combined_2 <- IntegrateData(anchorset = FB.anchors_2, normalization.method = "LogNormalize", dims = 1:maxCCA, k.weight = 100)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(FB.combined_2) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  FB.combined_2 <- ScaleData(FB.combined_2, verbose = FALSE)
  FB.combined_2 <- RunPCA(FB.combined_2, npcs = maxCCA, verbose = FALSE)

  # ########## Via PCelbowplot ##########
  # ElbowPlot(object = FB.combined_2, ndims = 30)
  
  ## Perform clustering
  dimsToTry<-c(seq(10,maxCCA,by=5))
  
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
    umapPlotSplit<-DimPlot(FB.combined_2, reduction = "umap", label = T, label.size = 2, group.by="New_clusters")
    
    ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
           file=paste0(output.dir,"Plots/RNA/16c_UMAP_",maxCCA,"_CCA_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
    
    clusterMatrix2<-FB.combined_2@meta.data
    umapTable2<-as.data.frame(FB.combined_2[['umap']]@cell.embeddings, stringsAsFactors = F)
    
    wantedCells<-rownames(FB.combined_2@meta.data[which(FB.combined_2$New_clusters == "Fibroblasts Type 2"),])
    p2<-colorSomeCells(clusterMatrix2, umapTable2, wantedCells) # 
    
    ggsave(p2, file=paste0(output.dir,"Plots/RNA/16c_extra_UMAP_",maxCCA,"_CCA_BBCs_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 12, height=10)
    
    resolutions <- seq(0.3,0.8,by=0.1)
    
    for(res in resolutions){
      FB.combined_2 <- FindClusters(object = FB.combined_2,  graph.name = "integrated_snn", resolution = res)
    }
    
    pdf(file=paste0(output.dir,"Plots/RNA/16c_extra_",maxCCA,"_CCA_overview_res_",min(dimsToUse),"-",max(dimsToUse),".pdf"))
    print(DimPlot(FB.combined_2, reduction = "umap", label = T, group.by = "integrated_snn_res.0.3", label.size = 2, repel = F))
    print(DimPlot(FB.combined_2, reduction = "umap", label = T, group.by = "integrated_snn_res.0.4", label.size = 2, repel = F))
    print(DimPlot(FB.combined_2, reduction = "umap", label = T, group.by = "integrated_snn_res.0.5", label.size = 2, repel = F))
    print(DimPlot(FB.combined_2, reduction = "umap", label = T, group.by = "integrated_snn_res.0.6", label.size = 2, repel = F))
    print(DimPlot(FB.combined_2, reduction = "umap", label = T, group.by = "integrated_snn_res.0.7", label.size = 2, repel = F))
    dev.off()
    
  }
}


################################################################################
################################################################################
### MANUAL PART -> choose final PCA dims and resolution
################################################################################
################################################################################

## Integration
# We then identify anchors using the FindIntegrationAnchors() function, 
# which takes a list of Seurat objects as input, and use these anchors to integrate the four datasets together with IntegrateData().
FB.anchors_2 <- FindIntegrationAnchors(object.list = FB_list_2, anchor.features = features_2, dims = 1:15, normalization.method = "LogNormalize", reduction = "cca")

FB.combined_2 <- IntegrateData(anchorset = FB.anchors_2, normalization.method = "LogNormalize", dims = 1:15, k.weight = 100)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(FB.combined_2) <- "integrated"

# Run the standard workflow for visualization and clustering
FB.combined_2 <- ScaleData(FB.combined_2, verbose = FALSE)
FB.combined_2 <- RunPCA(FB.combined_2, npcs = 15, verbose = FALSE)

### Final
dimsToTry<-c(10)
resToUse<-0.4
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
  umapPlotSplit<-DimPlot(FB.combined_2, reduction = "umap", label = T, group.by="New_clusters")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/16c_UMAP_final_check_15_CCA_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

## Check location of all cells original dataset annotation
clusterMatrix2<-FB.combined_2@meta.data
umapTable2<-as.data.frame(FB.combined_2[['umap']]@cell.embeddings, stringsAsFactors = F)

pdf(file=paste0(output.dir,"Plots/RNA/16c_UMAP_final_check_15_CCA_",min(dimsToUse),"-",max(dimsToUse),"_Color_old_annotation.pdf"), width = 12)
for (i in levels(as.factor(FB.combined_2@meta.data$New_clusters))) {
  C1<-colorSomeCells(clusterMatrix2, umapTable2, WhichCells(FB.combined_2, cells = rownames(FB.combined_2@meta.data[which(FB.combined_2@meta.data$New_clusters==i),])))
  C1<-C1+ggtitle(i)
  print(C1)
}
dev.off()

### Clustering: trying out clusTree
Perplexity<-10
Resolution<-0.4

FB.combined_2 <- FindNeighbors(object = FB.combined_2, reduction = "pca", dims = 1:Perplexity)
resolutions <- seq(0.1,1,by=0.1)

for(res in resolutions){
  FB.combined_2 <- FindClusters(object = FB.combined_2,  graph.name = "integrated_snn", resolution = res)
}

library("clustree")
pdf(file=paste0(output.dir,"Plots/RNA/16d_Clustree_final_check_15_CCA_",min(dimsToUse),"-",max(dimsToUse),".pdf"))
clustree(FB.combined_2, prefix = "integrated_snn_res.")
dev.off()

# Final Resolution and final clusters
res <- 0.4
diagnostics[['res']]<-res

DimPlot(FB.combined_2, reduction = "umap", label = T, group.by = "integrated_snn_res.0.4", label.size = 2, repel = F)

FB.combined_2$Integrated_RNA_clusters <- FB.combined_2$integrated_snn_res.0.4
Idents(FB.combined_2)<-FB.combined_2@meta.data$Integrated_RNA_clusters

## Save plot original annotation
pdf(file=paste0(output.dir,"Plots/RNA/16e_annotation_original_datasets_final_check_15_CCA_",min(dimsToUse),"-",max(dimsToUse),".pdf"), width = 20, height = 15)
DimPlot(FB.combined_2, reduction = "umap", label = T, group.by = "New_clusters", label.size = 3, repel = F)
dev.off()

########################################################################

## Save plot new unsupervised clustering
umapPlot<-DimPlot(FB.combined_2, reduction = "umap", label = T, group.by= "Integrated_RNA_clusters", label.size = 6)
FB.combined_2@active.assay

pdf(file=paste0(output.dir,"Plots/RNA/16f_tSNE_UMAP.pdf"), width = 17*0.45, height = 12.4*0.45)
umapPlot
dev.off()

########################################################################

experiment<-"Rebuttal7_Full_Merge_subset_CCA_v7_FB_datasets"

# Save objects
saveRDS(FB.combined_2, file = paste0(output.dir,"Robjects/seuratObj_",experiment,".rds"))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment, ".rds"))

# Read objects
FB.combined_2<-readRDS(file = paste0(output.dir,"Robjects/seuratObj_",experiment,".rds"))
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment, ".rds"))
