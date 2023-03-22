## Exploration of public human snRNA-Seq data (Yang et al.) and checking Fibroblast markers
## Data obtained from https://twc-stanford.shinyapps.io/scrna_brain_covid19/

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library("clustree")

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