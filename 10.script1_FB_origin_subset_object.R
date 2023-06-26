## Creating Fibroblast origin subset object with only the Fibroblast clusters from the Fibroblast origin complete object
## Workflow: don't take along prolif FBs (4) or PVM (19) or vascular cells and only take cells from cluster 1/3/5/8/11 
## Only took along cells from those clusters which mapped on the left side of the UMAP (performed manual tracing of the clusters)
## Initial basic Seurat and Harmony workflow script to create the Fibroblast origin subset object in Figure 1 manuscript
## Rebuttal: Stronger harmony settings

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

## Prepare folders
dir.create("Rebuttal_mouse_data/Merge_subset/")
dir.create("Rebuttal_mouse_data/Merge_subset/Plots")
dir.create("Rebuttal_mouse_data/Merge_subset/Plots/RNA")
dir.create("Rebuttal_mouse_data/Merge_subset/Robjects")

output.dir<-"Rebuttal_mouse_data/Merge_subset/"

### Read RDS object
seuratObjAll <- readRDS(file="Rebuttal_mouse_data/Full_merge/Robjects/seuratObj_Rebuttal2_Full_merge_FB_datasets_harmony_RNA.rds")

### Subset data ### (non-prolif FBs)
clusterMatrix<-seuratObjAll@meta.data
umapTable<-as.data.frame(seuratObjAll[['umap']]@cell.embeddings, stringsAsFactors = F)

Idents(seuratObjAll)<-seuratObjAll@meta.data$harmony_clusters
DimPlot(seuratObjAll, label = T)

U1 <- DimPlot(seuratObjAll, reduction = "umap", label = T, label.size = 4)
seuratObjAll <- CellSelector(U1, object=seuratObjAll, ident="Selected")
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObjAll, idents = c("Selected"))) #9806
wantedCells<-intersect(rownames(seuratObjAll@meta.data[which(seuratObjAll$harmony_clusters %in% c(1,3,5,8,11)),]), WhichCells(seuratObjAll, idents = c("Selected")))
colorSomeCells(clusterMatrix, umapTable, wantedCells) #9651

seuratObj<-subset(seuratObjAll, cells = wantedCells)

dim(seuratObj)
# 9651

diagnostics<-list()

unique(sapply(X = strsplit(colnames(seuratObj), split = "_"), FUN = "[", 1))
# "FB1"   "EP1"   "EP2"   "EP3"   "VASC1" "VASC2" "EP4"   "URVB"

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

################################################################################
########## RUN HARMONY
################################################################################
library('cowplot')
library("harmony")

########## Create vlnPlot before running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "pca", pt.size = 0.2, group.by = "sample_origin")
p2 <- VlnPlot(object = seuratObj, features = "PC_1", pt.size = 0.2, group.by = "sample_origin")
plot_grid(p1,p2)
ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/RNA/1a_vlnPlot_beforeAlignment.png"), width = 12, height = 8)

########## Run Harmony ##########
### Increase theta parameter in case of bad overlap!
options(repr.plot.height = 3, repr.plot.width = 6)
seuratObj<-RunHarmony(seuratObj, group.by.vars = c("sample_origin","Dataset"), theta = c(2,4), plot_convergence = TRUE, nclust = 50,
                      max.iter.cluster = 100, max.iter.harmony = 20, dims.use=1:40) #V1

### Get embeddings
harmony_embeddings <- Embeddings(seuratObj, 'harmony')
harmony_embeddings[1:5, 1:5]


########## Create vlnPlot after running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "harmony", pt.size = 0.2, group.by = "sample_origin")
p2 <- VlnPlot(object = seuratObj, features = "harmony_1", pt.size = 0.2, group.by = "sample_origin")
plot_grid(p1,p2)
ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/RNA/1b_vlnPlot_afterAlignment.png"), width = 12, height = 8)


########################################
########## Choose dims
########################################

########## Via heatmap ##########

pdf(file=paste0(output.dir,"/Plots/RNA/1c_heatmapHarmony.pdf"))
PCHeatmap(seuratObj, reduction = "harmony", dims = 1:12, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, reduction = "harmony", dims = 13:24, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, reduction = "harmony", dims = 25:36, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, reduction = "harmony", dims = 37:40, cells = 500, balanced = TRUE)
dev.off()

########## Via PCelbowplot ##########
ElbowPlot(object = seuratObj, ndims = 40)


dimsToTry<-c(seq(10,30,by=5))

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "harmony", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="harmony")
  umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

### Final
dimsToTry<-c(25)
resToUse<-0.8
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "harmony", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="harmony")
  umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

names(seuratObj)

### Clustering: trying out clusTree

Perplexity<-25
Resolution<-0.8
Perplexity_UMAP<-25

seuratObj <- FindNeighbors(object = seuratObj, reduction = "harmony", dims = 1:Perplexity)
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  seuratObj <- FindClusters(object = seuratObj,  resolution = res)
}

pdf(file=paste0(output.dir,"Plots/RNA/10c_Clustree.pdf"))
clustree(seuratObj, prefix = "RNA_snn_res.")
dev.off()

# 0.8 seems reasonable
# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
seuratObj$harmony_clusters_subset <- seuratObj$RNA_snn_res.0.8


################################################################################
################################################################################
### AUTOMATIC PART
################################################################################
################################################################################

umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, group.by= "harmony_clusters_subset", label.size = 6)
seuratObj@active.assay

pdf(file=paste0(output.dir,"Plots/RNA/11_tSNE_UMAP.pdf"), width = 17*0.45, height = 12.4*0.45)
umapPlot
dev.off()

experiment<-"Rebuttal2_Full_merge_subset_FB_datasets"

# Save objects
saveRDS(seuratObj, file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_harmony_RNA.rds"))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment, "_harmony_RNA.rds"))

# Read objects
seuratObj<-readRDS(file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_harmony_RNA.rds"))
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment, "_harmony_RNA.rds"))

