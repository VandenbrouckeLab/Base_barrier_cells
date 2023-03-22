## Human single nucleus CP -> check FB expression

# options(repos = c(getOption("repos"), BiocManager::repositories()))
# getOption("repos")
# packrat::get_opts()
# packrat::status()
# packrat::snapshot()

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library("clustree")

# BiocManager::install("impute") version = "3.8")
# BiocManager::install("preprocessCore") version = "3.8")
# BiocManager::install("GO.db") version = "3.8")
# BiocManager::install("AnnotationDbi") version = "3.8")
# devtools::install_github('dambi/DisneyTool@seuratv3', host="github.ugent.be/api/v3", auth_token = 'e5ca75c8c2f815aa7f1195cb0b6b6a3190064707')
#cutils::download.file("https://github.ugent.be/api/v3/repos/dambi/DisneyTool/tarball/master", destfile = "test.zip", method = "curl")
library('DisneyTools')

source('/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/RAW_DATA/script_functions_COVID.R') #KEVIN

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

# setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")
setwd("G:/VIB_G_drive/")
setwd("/run/media/clintdn/CN1465-DATA//VIB_G_drive/")

sampleName<-"Human_snRNAseq_CP_COVID"
sampleFolder<-paste0(sampleName,"/")

experiment <- sampleName
output.dir <- sampleFolder

##add some subfolders
dir.create(paste0(sampleFolder,"results"))
dir.create(paste0(sampleFolder,"Robjects"))

########################################
##### Some variables
########################################

### General variables
diagnostics<-list()

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')


################################################################################
########## LOAD DATA
################################################################################

##### Read object
seuratObj <- readRDS(file=paste0(sampleFolder,"Robjects/COVID-19_brain_snRNA-seq_choroid_plexus_final_seurat_v3.2.3.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,".rds"))

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/COVID-19_brain_snRNA-seq_choroid_plexus_final_seurat_v3.2.3.rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,".rds"))

##### Create new clusters
##new clusters
clusterMatrix<-seuratObj@meta.data
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

## DimPlot UMAP
D1<-DimPlot(seuratObj, reduction = "umap", label=T,repel = T, group.by="cellID", pt.size = 1)

pdf(file=paste0(sampleFolder,"results/UMAP_plot_",sampleName,".pdf"), width = 13, height = 10)
D1
dev.off()

Markers_cellID<-FindMarkers(seuratObj, ident.1 = "Mesenchymal", ident.2 = levels(Idents(seuratObj))[-2], assay = "integrated",
                            min.diff.pct = 0.25)

saveRDS(Markers_cellID, file=paste0(sampleFolder,"Robjects/Markers_cellID_",sampleName,".rds"))
Markers_cellID<-readRDS(file=paste0(sampleFolder,"Robjects/Markers_cellID_",sampleName,".rds"))

F1<-FeaturePlot(object = seuratObj, features = c("IGFBP6","CLDN11","ALPL","CDH11"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T)

pdf(file=paste0(sampleFolder,"results/Feature_plot_",sampleName,".pdf"), width = 15, height = 11)
F1
dev.off()

## Featureplots paper (May 2022)
features<-c("IGFBP6","CLDN11","ALPL","CDH11")

pdf(file=paste0(sampleFolder,"results_paper_Human_alone/Feature_plot_paper_4_markers_check_",sampleName,"_viridisC_ordered.pdf"), height = 11, width = 15)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) +
    scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

## Check for cDCs for Daan (19/04/22)
U1<-DimPlot(seuratObj, reduction = "umap", label=T,repel = T, group.by="cellID", pt.size = 1)
S1<-CellSelector(U1, seuratObj,ident="selectedCells")
FeaturePlot(seuratObj, features = "XCR1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T)
Markers<-FindMarkers(S1, ident.1 = "selectedCells", ident.2 = "Macrophage")
Markers_sig <- Markers %>% dplyr::filter(p_val_adj < 0.05)
Markers2<-FindMarkers(S1, ident.1 = "selectedCells", ident.2 = "Endothelial")
Markers2_sig <- Markers2 %>% dplyr::filter(p_val_adj < 0.05)

######################################################################################################
######################################################################################################

## Extra analysis May 2022: separate FBs and analyze again through pipeline
## Use integrated assay and rerun from HVG

dir.create(paste0(sampleFolder,"results_subset"))

seuratObj_subset<-subset(seuratObj, idents = "Mesenchymal")
DefaultAssay(seuratObj_subset)

################################################################################
########## GET HVG
################################################################################

seuratObj_subset <- FindVariableFeatures(object = seuratObj_subset, selection.method = "vst", nfeatures = 2000)

################################################################################
########## SCALE DATA
################################################################################

seuratObj_subset <- ScaleData(object = seuratObj_subset, features = VariableFeatures(seuratObj_subset))

################################################################################
########## PCA
################################################################################
seuratObj_subset <- RunPCA(object = seuratObj_subset, features = VariableFeatures(seuratObj_subset), npcs = 50, ndims.print = 1:5, nfeatures.print = 10)

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

### Create PCElbowplot
ElbowPlot(object = seuratObj_subset, ndims = 40)

################################################################################
########## CLUSTER THE CELLS (35 dims res 1.0)
################################################################################
# library(reticulate)
# use_condaenv(condaenv = "snowflakes", conda = "~/anaconda3/bin/conda", required = TRUE)
# reticulate::py_install(packages ='umap-learn') #Restart R session and reload packages https://github.com/satijalab/seurat/issues/1760

dimsToTry<-c(10,15,20,25)
resToUse<-0.8

### Final
dimsToTry<-c(20)
resToUse<-0.8
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse


for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj_subset <- FindNeighbors(object = seuratObj_subset, dims = dimsToUse)
  seuratObj_subset <- FindClusters(object = seuratObj_subset, resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj_subset <- RunUMAP(seuratObj_subset, dims = dimsToUse, n.neighbors = 30) #, umap.method = 'umap-learn', metric = 'correlation'
  umapPlot<-DimPlot(seuratObj_subset, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj_subset, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(sampleFolder,"results_subset/10b_UMAP_final_",min(dimsToUse),"-",max(dimsToUse),"-", resToUse,".png"), width = 15, height = 6)
  
}

### Clustering: trying out clusTree

Perplexity<-25
Resolution<-0.8
Perplexity_UMAP<-25

seuratObj_subset <- FindNeighbors(object = seuratObj_subset, reduction = "pca", dims = 1:Perplexity)
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  seuratObj_subset <- FindClusters(object = seuratObj_subset,  resolution = res)
}

pdf(file=paste0(sampleFolder,"results_subset/12c_Clustree.pdf"))
clustree(seuratObj_subset, prefix = "integrated_snn_res.")
dev.off()

DimPlot(seuratObj_subset, reduction = "umap", label = T, group.by = "integrated_snn_res.1", repel = T, label.size = 4) + NoLegend()
DimPlot(seuratObj_subset, reduction = "umap", label = T, group.by = "integrated_snn_res.0.8", repel = T, label.size = 4) + NoLegend()

# 0.8 seems reasonable
# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
seuratObj_subset$seurat_clusters <- seuratObj_subset$integrated_snn_res.0.8

## Check markers
F1<-FeaturePlot(object = seuratObj_subset, features = c("IGFBP6","CLDN11","ALPL","CDH11"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T)

pdf(file=paste0(sampleFolder,"results_subset/Feature_plot_integrated_",sampleName,".pdf"), width = 15, height = 11)
F1
dev.off()

## Read in BBKNN object for annotation
FB.combined_regressed <- readRDS(file = paste0("Human_snRNAseq_CP_COVID/Robjects/Ridge_regressed_seuratObj_new_BBKNN_integrated_Human_and_Mouse_CP.rds")) #new

DimPlot(FB.combined_regressed, reduction = "umap", label = T, group.by = "annotated_BBKNN_clusters", repel = T, label.size = 4) + NoLegend()

seuratObj_subset$BBKNN_annotation<-"Unknown"
seuratObj_subset$BBKNN_annotation <- FB.combined_regressed@meta.data[intersect(colnames(seuratObj_subset), colnames(FB.combined_regressed)),"annotated_BBKNN_clusters"]

D1<-DimPlot(seuratObj_subset, reduction = "umap", label = T, group.by = "BBKNN_annotation", repel = T, label.size = 4)

pdf(file=paste0(sampleFolder,"results_subset/Dimplot_BBKNN_annotation__integrated_",sampleName,".pdf"), width = 10, height = 8)
D1
dev.off()

##### Save object
saveRDS(seuratObj_subset, file=paste0(sampleFolder,"Robjects/Human_snRNA-seq_Mesenchymal_subset_integrated_assay.rds"))

##### Read object
seuratObj_subset <- readRDS(file=paste0(sampleFolder,"Robjects/Human_snRNA-seq_Mesenchymal_subset_integrated_assay.rds"))

######################################################################################################
######################################################################################################

## Extra analysis May 2022: separate FBs and analyze again through pipeline
## Retry with RNA assay

dir.create(paste0(sampleFolder,"results_subset"))

seuratObj_subset<-subset(seuratObj, idents = "Mesenchymal")
DefaultAssay(seuratObj_subset)<-"RNA"

################################################################################
########## NORMALIZE
################################################################################
seuratObj_subset <- NormalizeData(object = seuratObj_subset, normalization.method = "LogNormalize", scale.factor = 10000)

################################################################################
########## GET HVG
################################################################################

seuratObj_subset <- FindVariableFeatures(object = seuratObj_subset, selection.method = "vst", nfeatures = 2000)

################################################################################
########## SCALE DATA
################################################################################

seuratObj_subset <- ScaleData(object = seuratObj_subset, features = VariableFeatures(seuratObj_subset))

################################################################################
########## PCA
################################################################################
seuratObj_subset <- RunPCA(object = seuratObj_subset, features = VariableFeatures(seuratObj_subset), npcs = 50, ndims.print = 1:5, nfeatures.print = 10)

################################################################################
########## RUN HARMONY
################################################################################
library('cowplot')
library("harmony")

########## Create vlnPlot before running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj_subset, reduction = "pca", pt.size = 0.2, group.by = "orig.ident") #Adapted
p2 <- VlnPlot(object = seuratObj_subset, features = "PC_1", pt.size = 0.2, group.by = "orig.ident") #Adapted
plot_grid(p1,p2)
ggsave(plot_grid(p1, p2), file=paste0(output.dir,"results_subset/1a_vlnPlot_beforeAlignment.png"))

########## Run Harmony ##########
### Increase theta parameter in case of bad overlap!
options(repr.plot.height = 3, repr.plot.width = 6)
seuratObj_subset<-RunHarmony(seuratObj_subset, group.by.vars = "orig.ident", theta = 4, plot_convergence = TRUE, nclust = 50, #Adapted
                      max.iter.cluster = 100, max.iter.harmony = 20, dims.use=1:40)


### Get embeddings
harmony_embeddings <- Embeddings(seuratObj_subset, 'harmony')
harmony_embeddings[1:5, 1:5]


########## Create vlnPlot after running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj_subset, reduction = "harmony", pt.size = 0.2, group.by = "orig.ident") #Adapted
p2 <- VlnPlot(object = seuratObj_subset, features = "harmony_1", pt.size = 0.2, group.by = "orig.ident") #Adapted
plot_grid(p1,p2)
ggsave(plot_grid(p1, p2), file=paste0(output.dir,"results_subset/1b_vlnPlot_afterAlignment.png"))


########################################
########## Choose dims
########################################

########## Via PCelbowplot ##########
ElbowPlot(object = seuratObj_subset, ndims = 40)


dimsToTry<-20

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj_subset <- FindNeighbors(object = seuratObj_subset, reduction = "harmony", dims = dimsToUse)
  seuratObj_subset <- FindClusters(object = seuratObj_subset, resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj_subset <- RunUMAP(seuratObj_subset, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="harmony")
  umapPlot<-DimPlot(seuratObj_subset, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj_subset, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"results_subset/10b_UMAP_harmony_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

### Clustering: trying out clusTree

Perplexity<-20
Resolution<-0.8
Perplexity_UMAP<-20

seuratObj_subset <- FindNeighbors(object = seuratObj_subset, reduction = "harmony", dims = 1:Perplexity)
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  seuratObj_subset <- FindClusters(object = seuratObj_subset,  resolution = res)
}

pdf(file=paste0(sampleFolder,"results_subset/12c_Clustree_harmony.pdf"))
clustree(seuratObj_subset, prefix = "RNA_snn_res.")
dev.off()

DimPlot(seuratObj_subset, reduction = "umap", label = T, group.by = "RNA_snn_res.1", repel = T, label.size = 4) + NoLegend()
DimPlot(seuratObj_subset, reduction = "umap", label = T, group.by = "RNA_snn_res.0.8", repel = T, label.size = 4) + NoLegend()

# 0.8 seems reasonable
# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
seuratObj_subset$seurat_clusters <- seuratObj_subset$RNA_snn_res.0.8

## Check markers
F1<-FeaturePlot(object = seuratObj_subset, features = c("IGFBP6","CLDN11","ALPL","CDH11"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T)

pdf(file=paste0(sampleFolder,"results_subset/Feature_plot_Harmony_RNA_",sampleName,".pdf"), width = 15, height = 11)
F1
dev.off()

# Markers_12<-FindMarkers(seuratObj_subset, ident.1 = 12, ident.2 = c(1:11))

## Read in BBKNN object for annotation
FB.combined_regressed <- readRDS(file = paste0("Human_snRNAseq_CP_COVID/Robjects/Ridge_regressed_seuratObj_new_BBKNN_integrated_Human_and_Mouse_CP.rds")) #new

DimPlot(FB.combined_regressed, reduction = "umap", label = T, group.by = "annotated_BBKNN_clusters", repel = T, label.size = 4) + NoLegend()

seuratObj_subset$BBKNN_annotation<-"Unknown"
seuratObj_subset$BBKNN_annotation <- FB.combined_regressed@meta.data[intersect(colnames(seuratObj_subset), colnames(FB.combined_regressed)),"annotated_BBKNN_clusters"]

D1<-DimPlot(seuratObj_subset, reduction = "umap", label = T, group.by = "BBKNN_annotation", repel = T, label.size = 4)

pdf(file=paste0(sampleFolder,"results_subset/Dimplot_BBKNN_annotation_Harmony_RNA_",sampleName,".pdf"), width = 10, height = 8)
D1
dev.off()

##### Save object
saveRDS(seuratObj_subset, file=paste0(sampleFolder,"Robjects/Human_snRNA-seq_Mesenchymal_subset_Harmony_RNA_assay.rds"))

##### Read object
seuratObj_subset <- readRDS(file=paste0(sampleFolder,"Robjects/Human_snRNA-seq_Mesenchymal_subset_Harmony_RNA_assay.rds"))

######################################################################################################
######################################################################################################

## Extra analysis May 2022: separate FBs and analyze again through pipeline
## Retry with simple subset without reanalysis

dir.create(paste0(sampleFolder,"results_subset"))

seuratObj_subset<-subset(seuratObj, idents = "Mesenchymal")
DefaultAssay(seuratObj_subset)<-"integrated"

seuratObj_subset$BBKNN_annotation<-"Unknown"
seuratObj_subset$BBKNN_annotation <- FB.combined_regressed@meta.data[intersect(colnames(seuratObj_subset), colnames(FB.combined_regressed)),"annotated_BBKNN_clusters"]

D1<-DimPlot(seuratObj_subset, reduction = "umap", label = T, pt.size = 1, group.by = "BBKNN_annotation", repel = T, label.size = 4)

pdf(file=paste0(sampleFolder,"results_subset/Dimplot_BBKNN_annotation_simple_subset_",sampleName,".pdf"), width = 13, height = 10)
D1
dev.off()

## DimPlot UMAP
D2<-DimPlot(seuratObj_subset, reduction = "umap", label=T, repel = T, group.by="cellID", pt.size = 1, label.size = 4)

pdf(file=paste0(sampleFolder,"results_subset/UMAP_plot_simple_subset_",sampleName,".pdf"), width = 13, height = 10)
D2
dev.off()

## Featureplots paper (May 2022)
features<-c("IGFBP6","CLDN11","ALPL","CDH11")

pdf(file=paste0(sampleFolder,"results_subset/Feature_plot_paper_4_markers_simple_subset_",sampleName,"_viridisC_ordered.pdf"), height = 10, width = 13)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj_subset, features =feature, cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) +
    scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

##### Save object
saveRDS(seuratObj_subset, file=paste0(sampleFolder,"Robjects/Human_snRNA-seq_Mesenchymal_simple_subset.rds"))

##### Read object
seuratObj_subset <- readRDS(file=paste0(sampleFolder,"Robjects/Human_snRNA-seq_Mesenchymal_simple_subset.rds"))

######################################################################################################
######################################################################################################

## Extra analysis May 2022: separate FBs and analyze again through pipeline
## Retry with subset without reanalysis, but also include some cleaning!!!!

dir.create(paste0(sampleFolder,"results_subset"))

seuratObj_subset<-subset(seuratObj, idents = "Mesenchymal")
DefaultAssay(seuratObj_subset)<-"integrated"

U1<-DimPlot(seuratObj_subset, reduction = "umap")

seuratObj_subset2<-CellSelector(U1, object=seuratObj_subset, ident = "Outliers")
seuratObj_subset<-subset(seuratObj_subset2, idents = "Mesenchymal") #Remove outliers

seuratObj_subset$BBKNN_annotation<-"Unknown"
seuratObj_subset$BBKNN_annotation <- FB.combined_regressed@meta.data[intersect(colnames(seuratObj_subset), colnames(FB.combined_regressed)),"annotated_BBKNN_clusters"]

D1<-DimPlot(seuratObj_subset, reduction = "umap", label = T, pt.size = 1, group.by = "BBKNN_annotation", repel = T, label.size = 4)

pdf(file=paste0(sampleFolder,"results_subset/Dimplot_BBKNN_annotation_subset_and_cleaned_",sampleName,".pdf"), width = 13, height = 10)
D1
dev.off()

## DimPlot UMAP
D2<-DimPlot(seuratObj_subset, reduction = "umap", label=T, repel = T, group.by="cellID", pt.size = 1, label.size = 4)

pdf(file=paste0(sampleFolder,"results_subset/UMAP_plot_subset_and_cleaned_",sampleName,".pdf"), width = 13, height = 10)
D2
dev.off()

## Featureplots paper (May 2022)
features<-c("IGFBP6","CLDN11","ALPL","CDH11")

pdf(file=paste0(sampleFolder,"results_subset/Feature_plot_paper_4_markers_subset_and_cleaned_",sampleName,"_viridisC_ordered.pdf"), height = 10, width = 13)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj_subset, features =feature, cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) +
    scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

##### Save object
saveRDS(seuratObj_subset, file=paste0(sampleFolder,"Robjects/Human_snRNA-seq_Mesenchymal_subset_and_cleaned.rds"))

##### Read object
seuratObj_subset <- readRDS(file=paste0(sampleFolder,"Robjects/Human_snRNA-seq_Mesenchymal_subset_and_cleaned.rds"))