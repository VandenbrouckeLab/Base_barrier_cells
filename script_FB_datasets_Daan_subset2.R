## FB analysis Daan
## v2: Different type of subsetting -> don't take along prolif FBs and only take cells from cluster1, 9 and 16 in bottom half!!
## Due to issue with strange clusters in subset v1!!

library(Seurat)
library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(gridExtra)
library(clustree)

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")


################################################################################
################################################################################

## Prepare folders
dir.create("Merge_subset2")
dir.create("Merge_subset2/Plots")
dir.create("Merge_subset2/Plots/RNA")
dir.create("Merge_subset2/Robjects")

output.dir<-"Merge_subset2/"

### Read RDS object
seuratObjAll <- readRDS(file="Merge/Robjects/seuratObj_Merge_FB_datasets_harmony_RNA.rds")

### Subset data ### (non-prolif FBs)
clusterMatrix<-seuratObjAll@meta.data
tsneTable<-as.data.frame(seuratObjAll[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObjAll[['umap']]@cell.embeddings, stringsAsFactors = F)

Idents(seuratObjAll)<-seuratObjAll@meta.data$sliced_clusters
DimPlot(seuratObjAll, label = T)

umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_2 < -5)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObjAll, idents = c(1,9,16)))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-subset(seuratObjAll, cells = wantedCells)

dim(seuratObj)
# 11694 (subset 1) -> 9676 (subset 2 without prolif!)

diagnostics<-list()

unique(sapply(X = strsplit(colnames(seuratObj), split = "_"), FUN = "[", 1))
# [1] "FB1"   "EP1"   "EP3"   "VASC1" "VASC2" "EP4"   "URVB"
#No cells from EP2!!

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

## Use RNA going forward (to avoid mistakes!!!!!!!)
DefaultAssay(object = seuratObj)<-"RNA"

################################################################################
########## RUN HARMONY
################################################################################
library('cowplot')
library("harmony")

########## Create vlnPlot before running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "pca", pt.size = 0.2, group.by = "orig.ident")
p2 <- VlnPlot(object = seuratObj, features = "PC_1", pt.size = 0.2, group.by = "orig.ident")
plot_grid(p1,p2)
# ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/RNA/1a_vlnPlot_beforeAlignment.png"))

########## Run Harmony ##########
### Increase theta parameter in case of bad overlap!
options(repr.plot.height = 3, repr.plot.width = 6)
seuratObj<-RunHarmony(seuratObj, group.by.vars = "orig.ident", theta = 2, plot_convergence = TRUE, nclust = 50,
                      max.iter.cluster = 100, max.iter.harmony = 20, dims.use=1:40)


### Get embeddings
harmony_embeddings <- Embeddings(seuratObj, 'harmony')
harmony_embeddings[1:5, 1:5]


########## Create vlnPlot after running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "harmony", pt.size = 0.2, group.by = "orig.ident")
p2 <- VlnPlot(object = seuratObj, features = "harmony_1", pt.size = 0.2, group.by = "orig.ident")
plot_grid(p1,p2)
# ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/RNA/1b_vlnPlot_afterAlignment.png"))


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


dimsToTry<-c(seq(15,35,by=5))

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "harmony", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "RNA", reduction = "harmony")
  tsnePlot<-DimPlot(seuratObj, reduction = "tsne", label=T, label.size = 8)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "tsne", label=F, group.by="orig.ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
  
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
dimsToTry<-c(30)
resToUse<-0.8
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "harmony", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "RNA", reduction = "harmony")
  tsnePlot<-DimPlot(seuratObj, reduction = "tsne", label=T, label.size = 8)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "tsne", label=F, group.by="orig.ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="harmony")
  umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

names(seuratObj)

### Clustering: trying out clusTree

Perplexity<-30
Resolution<-0.8
Perplexity_UMAP<-30

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
seuratObj$harmony_clusters <- seuratObj$RNA_snn_res.0.8


################################################################################
################################################################################
### AUTOMATIC PART
################################################################################
################################################################################

# seuratObj <- RunTSNE(seuratObj, reduction = "SCT_pca", dims = 1:Perplexity, assay = "SCT")
TSNEPlot(seuratObj)

# seuratObj <- RunUMAP(seuratObj, dims = 1:Perplexity_UMAP, reduction = "SCT_pca", assay = "SCT")
umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, group.by= "harmony_clusters", label.size = 6)
tsnePlot<-TSNEPlot(seuratObj, reduction = "tsne", label = T, group.by= "harmony_clusters", label.size = 6)
seuratObj@active.assay

pdf(file=paste0(output.dir,"Plots/RNA/11_tSNE_UMAP.pdf"), width = 17*0.45, height = 12.4*0.45)
umapPlot
tsnePlot
dev.off()

experiment<-"Merge_FB_datasets_subset2"

# Save objects
saveRDS(seuratObj, file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_harmony_RNA.rds"))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment, "_harmony_RNA.rds"))

# Read objects
seuratObj<-readRDS(file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_harmony_RNA.rds"))
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment, "_harmony_RNA.rds"))

