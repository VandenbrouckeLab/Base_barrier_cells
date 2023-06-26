## Script for merging our datasets with public datasets to investigate the origin of our FBs
## Initial basic Seurat and Harmony workflow script to create Fibroblast origin complete object in Figure 1 manuscript
## Rebuttal: Stronger harmony settings

library(Seurat)
library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(gridExtra)
library(clustree)

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")


## Desisto fibroblasts (No filter necessary!!)
Desisto_FB<-read.csv(file = "GSE150219_Desisto/Foxc1_1.2.3_fibroblast.data.csv") # Log counts!!!!!!!!!!!!!!!!!!!!!!!!
Desisto_FB_metadata<-read.csv(file = "GSE150219_Desisto/Foxc1_1.2.3_fibroblast.metadata.csv")

rownames(Desisto_FB)<-Desisto_FB[,1]
Desisto_FB<-Desisto_FB[,-1]

colnames(Desisto_FB)<-gsub("X1","1",colnames(Desisto_FB))
colnames(Desisto_FB)<-gsub("X2","2",colnames(Desisto_FB))
colnames(Desisto_FB)<-gsub("X3","3",colnames(Desisto_FB))

rownames(Desisto_FB_metadata)<-Desisto_FB_metadata[,1]
Desisto_FB_metadata<-Desisto_FB_metadata[,-1]

Desisto_FB[1:5,1:5]

seuratObj1 <- CreateSeuratObject(counts = Desisto_FB, project = "Desisto_FB", min.cells = 3, min.features = 200) #Already log transformed!
seuratObj1@meta.data<-Desisto_FB_metadata

####################################

## Shah seuratobjects
seuratObj2<-readRDS(file="GSE100320_Shah/GSM_2677817/Robjects/seuratObj_final_GSM_2677817.rds")
seuratObj3<-readRDS(file="GSE100320_Shah/GSM_2677818/Robjects/seuratObj_final_GSM_2677818.rds")
seuratObj4<-readRDS(file="GSE100320_Shah/GSM_2677819/Robjects/seuratObj_final_GSM_2677819.rds")


####################################

## Vanlandewijck

seuratObj5<-readRDS(file="GSE98816_Vanlandewijck/Robjects/seuratObj_final_GSE98816_Vanlandewijck.rds")


####################################

## Atlas seuratobjects

seuratObj6<-readRDS(file="Atlas_Zeisel/Robjects/seuratObj_Vascular_cells_atlas.rds")
seuratObj7<-readRDS(file="Atlas_Zeisel/Robjects/seuratObj_Ependymal_cells_atlas.rds")

####################################

## Seuratobject Daan and Nina (Updated!! Now no immune cells or CPE cells!!)

seuratObj8<-readRDS(file="Urvb_datasets_rebuttal//Robjects/seuratObj_final_Urvb_datasets_rebuttal.rds")

################################################################################

## Update Seuratobjects (Now in R 4.0.5)
seuratObj2<-UpdateSeuratObject(seuratObj2)
seuratObj3<-UpdateSeuratObject(seuratObj3)
seuratObj4<-UpdateSeuratObject(seuratObj4)
seuratObj5<-UpdateSeuratObject(seuratObj5)
seuratObj6<-UpdateSeuratObject(seuratObj6)
seuratObj7<-UpdateSeuratObject(seuratObj7)
seuratObj8<-UpdateSeuratObject(seuratObj8)

################################################################################
################################################################################


## Prepare folders
dir.create("Rebuttal_mouse_data/Full_merge/")
dir.create("Rebuttal_mouse_data/Full_merge/Plots")
dir.create("Rebuttal_mouse_data/Full_merge/Plots/RNA")
dir.create("Rebuttal_mouse_data/Full_merge/Robjects")

output.dir<-"Rebuttal_mouse_data/Full_merge/"

## Merge all seuratobjects
seuratObj <- merge(seuratObj1, y = c(seuratObj2,seuratObj3,seuratObj4,seuratObj5,seuratObj6,seuratObj7,seuratObj8),
                   add.cell.ids = c("FB1","EP1","EP2","EP3","VASC1","VASC2","EP4","URVB"), project = "FB_merge", merge.data = T)

dim(seuratObj)
# [1] 29623 30800

diagnostics<-list()

unique(sapply(X = strsplit(colnames(seuratObj), split = "_"), FUN = "[", 1))


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

## Extra split for different thetas
levels(as.factor(seuratObj$orig.ident))
seuratObj@meta.data$sample_origin<-as.factor(seuratObj@meta.data$orig.ident)
levels(seuratObj@meta.data$sample_origin)<-c("FB_DeSisto_et_al_1","FB_DeSisto_et_al_2","FB_DeSisto_et_al_3",
                                             "Ependymal_Zeisel_et_al","Vascular_Vanlandewijck_et_al",
                                             "Ependymal_Shah_et_al_1","Ependymal_Shah_et_al_2","Ependymal_Shah_et_al_3",
                                             "CP_LpsNeg_4V_Urvb","CP_LpsNeg_LV_Urvb","CP_Young_4V_Urvb","CP_Young_LV_Urvb",
                                             "Vascular_Zeisel_et_al")
seuratObj$Dataset <- as.factor(seuratObj@meta.data$orig.ident)
levels(seuratObj@meta.data$Dataset)<-c("FB_DeSisto_et_al","FB_DeSisto_et_al","FB_DeSisto_et_al",
                                       "Ependymal_Zeisel_et_al","Vascular_Vanlandewijck_et_al",
                                       "Ependymal_Shah_et_al","Ependymal_Shah_et_al","Ependymal_Shah_et_al",
                                       "Verhaege_et_al","Verhaege_et_al","Verhaege_et_al","Verhaege_et_al",
                                       "Vascular_Zeisel_et_al")

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


dimsToTry<-c(seq(15,40,by=5))

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
dimsToTry<-c(15)
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

Perplexity<-15
Resolution<-0.8
Perplexity_UMAP<-15

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

umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, group.by= "harmony_clusters", label.size = 6)
seuratObj@active.assay

pdf(file=paste0(output.dir,"Plots/RNA/11_UMAP.pdf"), width = 17*0.45, height = 12.4*0.45)
umapPlot
dev.off()

experiment<-"Rebuttal2_Full_merge_FB_datasets"

# Save objects
saveRDS(seuratObj, file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_harmony_RNA.rds"))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment, "_harmony_RNA.rds"))

# Read objects
seuratObj<-readRDS(file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_harmony_RNA.rds"))
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment, "_harmony_RNA.rds"))

