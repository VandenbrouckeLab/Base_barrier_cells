## Other option: 3V removed from Lehtinen data!!!
## CCA used to perform integration between datasets

## FB analysis Daan
library(Seurat)
library(ggplot2)
# library(SingleCellExperiment)
# library(scran)
# library(scater)
library(dplyr)
library(gridExtra)
library(clustree)

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")


####################################

## Seuratobject Daan and Nina 6 datasets

seuratObj1<-readRDS(file="Urvb_6datasets/Robjects/seuratObj_final_Urvb_6datasets.rds")

####################################

## Seuratobject Lehtinen
seuratObj2<-readRDS(file="Lehtinen_CP_scRNAseq//Robjects/seuratObj_mesenchymal_Lehtinen_CP.rds")

## Subset without 3V
seuratObj2@meta.data$Ventricle<-"CP"
seuratObj2@meta.data[grepl("3V", rownames(seuratObj2@meta.data)), "Ventricle"] <- "3V"
seuratObj2@meta.data[grepl("TCP", rownames(seuratObj2@meta.data)), "Ventricle"] <- "LV"
seuratObj2@meta.data[grepl("HCP", rownames(seuratObj2@meta.data)), "Ventricle"] <- "4V"

seuratObj2<-subset(seuratObj2, cells = rownames(seuratObj2@meta.data[which(seuratObj2@meta.data$Ventricle != "3V"),]))

################################################################################
################################################################################

## Prepare folders
dir.create("CCA_Fig2")
dir.create("CCA_Fig2/Plots")
dir.create("CCA_Fig2/Plots/RNA")
dir.create("CCA_Fig2/Robjects")
dir.create("CCA_Fig2/results")

sampleFolder<-"CCA_Fig2/"
sampleName<-"CCA_FB_6datasets_Fig2"

## DefaultAssay
DefaultAssay(seuratObj1)
DefaultAssay(seuratObj2)

## Normalize
seuratObj1 <- NormalizeData(object = seuratObj1, normalization.method = "LogNormalize", scale.factor = 10000)
seuratObj2 <- NormalizeData(object = seuratObj2, normalization.method = "LogNormalize", scale.factor = 10000)

## Find Variable Features first for Human data
seuratObj1 <- FindVariableFeatures(seuratObj1, assay = "RNA", selection.method = "vst", nfeatures = 2000)
seuratObj2 <- FindVariableFeatures(seuratObj2, assay = "RNA", selection.method = "vst", nfeatures = 2000)

# select features that are repeatedly variable across datasets for integration
FB_list <- list(seuratObj1, seuratObj2)
features <- SelectIntegrationFeatures(object.list = FB_list)

## Integration
# We then identify anchors using the FindIntegrationAnchors() function, 
# which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
FB.anchors <- FindIntegrationAnchors(object.list = FB_list, anchor.features = features, dims = 1:30)

FB.combined <- IntegrateData(anchorset = FB.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(FB.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
FB.combined <- ScaleData(FB.combined, verbose = FALSE)
FB.combined <- RunPCA(FB.combined, npcs = 30, verbose = FALSE)
FB.combined <- RunUMAP(FB.combined, reduction = "pca", dims = 1:30)
FB.combined <- FindNeighbors(FB.combined, reduction = "pca", dims = 1:30)
FB.combined <- FindClusters(FB.combined, resolution = 0.5)

## Add metadata column (Lab)
FB.combined@meta.data$Lab<-FB.combined@meta.data$orig.ident
FB.combined@meta.data[which(FB.combined@meta.data$Lab == "RVD1_LpsNegFour" | FB.combined@meta.data$Lab == "RVD2_LpsNegLat" |
                            FB.combined@meta.data$Lab == "RVD5_Y4V" | FB.combined@meta.data$Lab == "RVD6_YLV" |
                            FB.combined@meta.data$Lab == "RVD7_O4V" | FB.combined@meta.data$Lab == "RVD8_OLV"),
                    "Lab"]<-"Vandenbroucke"
FB.combined@meta.data[which(FB.combined@meta.data$Lab == "E.HCP.V1.Rep1" |
                            FB.combined@meta.data$Lab == "E.HCP.V2.Rep1" | FB.combined@meta.data$Lab == "E.HCP.V2.Rep2" |
                            FB.combined@meta.data$Lab == "E.TCP.V1.Rep1" | FB.combined@meta.data$Lab == "E.TCP.V2.Rep1" | 
                            FB.combined@meta.data$Lab == "E.TCP.V2.Rep2" ),
                    "Lab"]<-"Lehtinen"


# Visualization
p3 <- DimPlot(FB.combined, reduction = "umap", group.by = "Lab") + ggtitle("CCA workflow")
p4 <- DimPlot(FB.combined, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Numbered_annotation_CCA")
p3+p4

FeaturePlot(FB.combined, features = c("Igfbp6","Cldn11"), min.cutoff = "q10", max.cutoff = "q90")

##################

## Read other seuratObject for annotation
seuratObjOld<-readRDS("Merge_Fig2/Robjects/seuratObj_Merge_FB_6datasets_Fig2_harmony_RNA_check.rds")

## Adjust cell names
Test_meta_data<-seuratObjOld@meta.data
rownames(Test_meta_data)<-gsub("URVB_","",rownames(Test_meta_data))
rownames(Test_meta_data)<-gsub("Lehtinen_","",rownames(Test_meta_data))

## Get old annotation and put in new seuratobject
FB.combined@meta.data[colnames(FB.combined),"annotated_clusters_old"] <- Test_meta_data[colnames(FB.combined),"annotated_clusters"]

DimPlot(FB.combined, reduction = "umap", group.by = "annotated_clusters_old", label = TRUE, repel = TRUE) + ggtitle("Old_annotation_Harmony")

#######################################################################
######################################################################

dim(FB.combined)
# [1] 2000 4561

diagnostics<-list()

## Add new metadata
FB.combined@meta.data$Author<-FB.combined@meta.data$orig.ident
FB.combined@meta.data$Age<-FB.combined@meta.data$orig.ident
FB.combined@meta.data$Ventricle<-FB.combined@meta.data$orig.ident

FB.combined@meta.data[which(FB.combined@meta.data$Author == "RVD1_LpsNegFour" | FB.combined@meta.data$Author == "RVD2_LpsNegLat"),
                    "Author"]<-"CP_Verhaege_et_al"
FB.combined@meta.data[which(FB.combined@meta.data$Author == "RVD5_Y4V" | FB.combined@meta.data$Author == "RVD6_YLV" |
                            FB.combined@meta.data$Author == "RVD7_O4V" | FB.combined@meta.data$Author == "RVD8_OLV"),
                    "Author"]<-"CP_Gorle_et_al"
FB.combined@meta.data[which(FB.combined@meta.data$Author == "E.HCP.V1.Rep1" |
                            FB.combined@meta.data$Author == "E.HCP.V2.Rep1" | FB.combined@meta.data$Author == "E.HCP.V2.Rep2" |
                            FB.combined@meta.data$Author == "E.TCP.V1.Rep1" | FB.combined@meta.data$Author == "E.TCP.V2.Rep1" | 
                            FB.combined@meta.data$Author == "E.TCP.V2.Rep2" ),
                    "Author"]<-"CP_Dani_et_al"

FB.combined@meta.data[which(FB.combined@meta.data$Age == "RVD1_LpsNegFour" | FB.combined@meta.data$Age == "RVD2_LpsNegLat"),
                    "Age"]<-"7w"
FB.combined@meta.data[which(FB.combined@meta.data$Age == "RVD5_Y4V" | FB.combined@meta.data$Age == "RVD6_YLV"),
                    "Age"]<-"22w"
FB.combined@meta.data[which(FB.combined@meta.data$Age == "RVD7_O4V" | FB.combined@meta.data$Age == "RVD8_OLV"),
                    "Age"]<-"82w"
FB.combined@meta.data[which(FB.combined@meta.data$Age == "E.HCP.V1.Rep1" |
                            FB.combined@meta.data$Age == "E.HCP.V2.Rep1" | FB.combined@meta.data$Age == "E.HCP.V2.Rep2" |
                            FB.combined@meta.data$Age == "E.TCP.V1.Rep1" | FB.combined@meta.data$Age == "E.TCP.V2.Rep1" | 
                            FB.combined@meta.data$Age == "E.TCP.V2.Rep2" ),
                    "Age"]<-"Embryonal"


FB.combined@meta.data[which(FB.combined@meta.data$Ventricle == "RVD1_LpsNegFour" | FB.combined@meta.data$Ventricle == "RVD5_Y4V" |
                            FB.combined@meta.data$Ventricle == "RVD7_O4V" | FB.combined@meta.data$Ventricle == "E.HCP.V1.Rep1" |
                            FB.combined@meta.data$Ventricle == "E.HCP.V2.Rep1" | FB.combined@meta.data$Ventricle == "E.HCP.V2.Rep2"),
                    "Ventricle"]<-"4V"
FB.combined@meta.data[which(FB.combined@meta.data$Ventricle == "RVD2_LpsNegLat" | FB.combined@meta.data$Ventricle == "RVD6_YLV" |
                            FB.combined@meta.data$Ventricle == "RVD8_OLV" | FB.combined@meta.data$Ventricle == "E.TCP.V1.Rep1" | 
                            FB.combined@meta.data$Ventricle == "E.TCP.V2.Rep1" | FB.combined@meta.data$Ventricle == "E.TCP.V2.Rep2"),
                    "Ventricle"]<-"LV"

head(FB.combined@meta.data)

####################################################################################################################

########## Via PCelbowplot ##########
ElbowPlot(object = FB.combined, ndims = 30)

## Perform clustering
dimsToTry<-c(seq(15,30,by=5))

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  FB.combined <- FindNeighbors(object = FB.combined, reduction = "pca", dims = dimsToUse)
  FB.combined <- FindClusters(object = FB.combined,  graph.name = "integrated_snn", resolution = resToUse)
  
  # ##### Create tSNE plot
  # FB.combined <- RunTSNE(object = FB.combined, dims = dimsToUse, assay = "RNA", reduction = "RNA_pca", 
  #                      reduction.name = "RNA_tsne", reduction.key = "RNATSNE_")
  # tsnePlot<-DimPlot(FB.combined, reduction = "RNA_tsne", label=T, label.size = 8)
  # tsnePlotSplit<-DimPlot(FB.combined, reduction = "RNA_tsne", label=F, group.by="ident", pt.size = 2)
  # 
  # ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
  #        file=paste0(sampleFolder,"Plots/RNA/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  # 
  
  ##### Create UMAP plot
  FB.combined <- RunUMAP(FB.combined, dims = dimsToUse, n_neighbors = 30, assay = "integrated", reduction ="pca",
                           reduction.name = "umap", reduction.key = "integratedUMAP_")
  umapPlot<-DimPlot(FB.combined, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(FB.combined, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(sampleFolder,"Plots/RNA/UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

### Final
dimsToTry<-c(30)
resToUse<-0.8
diagnostics<-list()
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  FB.combined <- FindNeighbors(object = FB.combined, reduction = "pca", dims = dimsToUse)
  FB.combined <- FindClusters(object = FB.combined,  graph.name = "integrated_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  FB.combined <- RunUMAP(FB.combined, dims = dimsToUse, n_neighbors = 30, assay = "integrated", reduction ="pca",
                           reduction.name = "umap", reduction.key = "integratedUMAP_")
  umapPlot<-DimPlot(FB.combined, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(FB.combined, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(sampleFolder,"Plots/RNA/UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

names(FB.combined)

### Clustering: trying out clusTree

Perplexity<-30
Resolution<-0.8
Perplexity_UMAP<-30

FB.combined <- FindNeighbors(object = FB.combined, reduction = "pca", dims = 1:Perplexity)
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  FB.combined <- FindClusters(object = FB.combined,  graph.name = "integrated_snn", resolution = res)
}

library("clustree")
pdf(file=paste0(sampleFolder,"Plots/RNA/Clustree.pdf"))
clustree(FB.combined, prefix = "integrated_snn_res.")
dev.off()

# 0.8 seems reasonable
# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
FB.combined$Integrated_RNA_clusters <- FB.combined$integrated_snn_res.0.8
Idents(FB.combined)<-FB.combined@meta.data$Integrated_RNA_clusters

umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, group.by= "Integrated_RNA_clusters", label.size = 6)

pdf(file=paste0(sampleFolder,"Plots/RNA/11_UMAP.pdf"), width = 17*0.45, height = 12.4*0.45)
umapPlot
dev.off()

##### Save object
saveRDS(FB.combined, file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))
# saveRDS(diagnostics, file=paste0("CCA_FB_combined_2/diagnostics_",sampleName,".rds"))

##### Read object
FB.combined<- readRDS(file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))
# diagnostics<- readRDS(file=paste0("CCA_FB_combined_2/diagnostics_",sampleName,".rds"))

##########################################################################################################
##########################################################################################################

### Find RNAmarkers for every Integrated cluster compared to all remaining cells, report only the positive ones
library(future)
plan("multiprocess", workers = 8)

RNAMarkers_RNAclus <- FindAllMarkers(FB.combined, assay = "integrated", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"/Robjects/IntegratedmarkersList_integratedclus_",sampleName,".rds"))

# ### Add to diagnostics
# diagnostics[['RNAmarkersPerIntegratedcluster']]<-paste0(table(RNAMarkers_RNAclus$cluster)," RNA markers for integrated cluster ",rownames(table(RNAMarkers_RNAclus$cluster)))
# saveRDS(diagnostics, file=paste0(sampleFolder,"diagnostics_",sampleName,".rds"))

### Create list with markers
totalNrRNAclusters_RNAclus<-max(as.numeric(names(table(RNAMarkers_RNAclus$cluster))))
totalNrRNAclusters_RNAclusPlusOne<-totalNrRNAclusters_RNAclus+1
RNAmarkersList_RNAclus<-list()

for(i in 1:totalNrRNAclusters_RNAclusPlusOne){
  clusterNr<-i-1
  
  tmp<-RNAMarkers_RNAclus[RNAMarkers_RNAclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  RNAmarkersList_RNAclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(RNAmarkersList_RNAclus)<-paste0("Integratedcluster",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results/IntegratedmarkersList_Integratedclus_",sampleName,".xlsx"))


######################################
##### Annotation
########################################

## Annotate clusters reference
D1 <- DimPlot(FB.combined, reduction = "umap", group.by = "Integrated_RNA_clusters", label = TRUE, repel = TRUE)

D2 <- DimPlot(FB.combined, reduction = "umap", group.by = "annotated_clusters_old", label = TRUE, repel = TRUE) + ggtitle("Old_annotation_Harmony")

pdf(file=paste0(sampleFolder,"results/UMAP_old_annotation_",sampleName,".pdf"), width = 20, height = 10)
D1 + D2
dev.off()

Idents(FB.combined)<-FB.combined$Integrated_RNA_clusters

## Split off Slc1a3+ FBs (even though not separate cluster) -> pretty clear in UMAP space!!!
U1 <- DimPlot(FB.combined, reduction = "umap", group.by = "Integrated_RNA_clusters", label = TRUE, repel = TRUE)

cells.located <- CellSelector(plot = U1)
FB.combined <- SetIdent(FB.combined, cells = cells.located, value =  14)

## Clean object further (consolidate contamination clusters)
U1 <- DimPlot(FB.combined, reduction = "umap", label = TRUE, repel = TRUE)
cells.located <- CellSelector(plot = U1)
FB.combined <- SetIdent(FB.combined, cells = cells.located, value =  12)

U1 <- DimPlot(FB.combined, reduction = "umap", label = TRUE, repel = TRUE)
cells.located <- CellSelector(plot = U1)
FB.combined <- SetIdent(FB.combined, cells = cells.located, value =  8)

## Check dimplot
DimPlot(FB.combined, reduction = "umap", label = TRUE, repel = TRUE)

FB.combined$sliced_clusters<-Idents(FB.combined)

FB.combined@meta.data$sliced_clusters<-factor(FB.combined@meta.data$sliced_clusters,levels=sort(as.numeric(levels(FB.combined@meta.data$sliced_clusters)))) #reorder levels

FB.combined@meta.data$Integrated_annotated_clusters <- FB.combined@meta.data$sliced_clusters
levels(FB.combined@meta.data$Integrated_annotated_clusters) <- c(rep("Stromal Fibroblasts",3),"Other Fibroblasts","Stalk Fibroblasts",
                                                                 "Stromal Fibroblasts",rep("Mural/EC-like Fibroblasts",2),"CPE doublets",
                                                                 rep("Proliferating Fibroblasts",2),"MF doublets", "EC doublets",
                                                                 "ABCs", "Slc1a3+ Fibroblasts")

U_annot<-DimPlot(FB.combined, reduction = "umap", label = T, group.by = "Integrated_annotated_clusters", repel = T, label.size = 4) + NoLegend()
ggsave(U_annot, file=paste0(sampleFolder,"results/UMAP_annotated_cluster_",sampleName,".png"), height = 10, width = 12, dpi = "retina")

Idents(FB.combined)<-FB.combined@meta.data$Integrated_annotated_clusters

###############################################################
###############################################################

FB.combined_subset <- subset(FB.combined, idents = levels(Idents(FB.combined))[-c(5,7,8)])
FB.combined_subset$Integrated_annotated_clusters<-as.factor(as.character(FB.combined_subset$Integrated_annotated_clusters))

##################################################################################################################

## Check Igfbp6 and Cldn11 expression

# Featureplots on full object
# Module score Seurat
library(viridis)

FB_stalk_signature <- list(c("Cldn11","Igfbp6"))

FB.combined_subset <- AddModuleScore(object = FB.combined_subset, assay = "RNA", features = FB_stalk_signature, name = "FB_stalk_signature_score")
F2 <- FeaturePlot(object = FB.combined_subset, features = "FB_stalk_signature_score1", order = T)  + scale_color_viridis(option = "C")
pdf(file = paste0(sampleFolder,"results/2_UMAP_Featureplot_modulescore_Cldn11_Igfbp6_",sampleName,".pdf"), width = 12, height = 10)
F2
dev.off()

Idents(FB.combined_subset)<-FB.combined_subset@meta.data$Integrated_annotated_clusters
F2.5 <- FeaturePlot(object = FB.combined_subset, features = "FB_stalk_signature_score1", order = T, label = T, repel = T, label.size = 3, cols = c("Yellow","Red"))
pdf(file = paste0(sampleFolder,"results/2_UMAP_Featureplot_modulescore_Cldn11_Igfbp6_labeled_",sampleName,".pdf"), width = 12, height = 10)
F2.5
dev.off()

# Extra featureplots
DefaultAssay(FB.combined_subset)<-"RNA"
F1 <- FeaturePlot(object = FB.combined_subset, features = c("Cldn11","Igfbp6","Dpep1","Alpl"), order = T, pt.size = 1, cols = c("Yellow","Red"), combine = F)
pdf(file = paste0(sampleFolder,"results/2_UMAP_Featureplots_FB_markers_RNA_",sampleName,".pdf"), width = 12, height = 10)
F1
dev.off()

DefaultAssay(FB.combined_subset)<-"integrated"
F1.5 <- FeaturePlot(object = FB.combined_subset, features = c("Cldn11","Igfbp6","Dpep1","Alpl"), order = T, pt.size = 1, cols = c("Yellow","Red"), combine = F)
pdf(file = paste0(sampleFolder,"results/2_UMAP_Featureplots_FB_markers_integrated_",sampleName,".pdf"), width = 12, height = 10)
F1.5
dev.off()

## Featureplots paper (May 2022)
features<-c("Dcn", "Dpep1", "Igfbp6" , "Cldn11")

pdf(file=paste0(sampleFolder,"results/Paper/Feature_plot_paper_4_markers_check_",sampleName,"_viridisC_ordered.pdf"), height = 10, width = 12)
for (feature in features) {
  F1<-FeaturePlot(object = FB.combined_subset, features =feature, cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) +
    scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

# Create new clusters: split on source
FB.combined_subset@meta.data$newClustersTmp<-FB.combined_subset@meta.data$Integrated_annotated_clusters
FB.combined_subset@meta.data$newClusters<-paste0(FB.combined_subset@meta.data$newClustersTmp,"_",FB.combined_subset@meta.data$Age)
head(FB.combined_subset@meta.data)

# Subset stalk and violin plot
Idents(FB.combined_subset)<-FB.combined_subset@meta.data$newClusters
FB.combined_double_subset<-subset(FB.combined_subset, ident = c("Stalk Fibroblasts_Embryonal","Stalk Fibroblasts_7w", "Stalk Fibroblasts_22w","Stalk Fibroblasts_82w",
                                                          "Stromal Fibroblasts_Embryonal","Stromal Fibroblasts_7w", "Stromal Fibroblasts_22w","Stromal Fibroblasts_82w"))

# ## Initial order
# Idents(FB.combined_double_subset)<-factor(FB.combined_double_subset@active.ident, levels = levels(Idents(FB.combined_double_subset))[c(8,7,2,1,4,3,6,5)])
# DimPlot(FB.combined_double_subset)
# 
# library(RColorBrewer)
# ColorSet<-brewer.pal(n =8, name = "Paired")
# 
# V1<-VlnPlot(FB.combined_double_subset, features = c("Igfbp6","Cldn11","Dpep1","Alpl"), assay = "RNA", ncol = 4, cols = ColorSet) & theme(plot.margin = margin(0.5,0.5,0.5,1, "cm"))
# V1.5<-VlnPlot(FB.combined_double_subset, features = c("Igfbp6","Cldn11","Dpep1","Alpl"), assay = "integrated", ncol = 4, cols = ColorSet) & theme(plot.margin = margin(0.5,0.5,0.5,1, "cm"))
# 
# pdf(file=paste0(sampleFolder,"results/3_Test_violinplot_RNA_",sampleName,".pdf"), width = 25, height= 10)
# V1
# dev.off()
# 
# pdf(file=paste0(sampleFolder,"results/3_Test_violinplot_integrated_",sampleName,".pdf"), width = 25, height= 10)
# V1.5
# dev.off()

## New order
Idents(FB.combined_double_subset)<-factor(FB.combined_double_subset@active.ident, levels = levels(Idents(FB.combined_double_subset))[c(8,2,4,6,7,1,3,5)])
DimPlot(FB.combined_double_subset)

library(RColorBrewer)
ColorSet<-brewer.pal(n =8, name = "Paired")[c(1,3,5,7,2,4,6,8)]

V1<-VlnPlot(FB.combined_double_subset, features = c("Igfbp6","Cldn11","Dpep1","Alpl"), assay = "RNA", ncol = 4, cols = ColorSet) & theme(plot.margin = margin(0.5,0.5,0.5,1, "cm"))

pdf(file=paste0(sampleFolder,"results/3_Test_violinplot_RNA_order_Roos_",sampleName,".pdf"), width = 25, height= 10)
V1
dev.off()

################################################################################

## Sample distribution across labs

# create a dataset
Sample <- FB.combined_subset@meta.data$Lab
cluster <- FB.combined_subset$Integrated_annotated_clusters
# Aggr <- rep(experiment,length(cluster)) 

data <- data.frame(table(Sample, cluster))
# data2 <- data.frame(table(cluster,Aggr))

# barplotAggr(seuratObj, listLabels)

# Stacked
png(file=paste0(sampleFolder,"results/4_SampleDistribution_ggplot2_annot_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

# png(file=paste0(sampleFolder,"CCA_FB_combined_2/3_SampleDistribution_ggplot2_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
# ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
#   geom_bar(position="fill", stat="identity", colour="white")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
# dev.off()

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split1<-data[which(data$Sample=="Lehtinen"),]
data_split2<-data[which(data$Sample=="Vandenbroucke"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})


data_new<-rbind(data_split1,data_split2)

# png(file=paste0(sampleFolder,"results/4_SampleDistribution_ggplot2_annot_1_",sampleName,"_adjusted.png"), width = 2000, height = 1500, res = 300)
S2<-ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red',"Turquoise")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
# dev.off()

pdf(file=paste0(sampleFolder,"results/4_SampleDistribution_ggplot2_annotated_v2_",sampleName,"_adjusted.pdf"), width =8, height= 6)
S2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #Remove grid lines
dev.off()


################################################################################

## Sample distribution across ventricles

# create a dataset
Sample <- FB.combined_subset@meta.data$Ventricle
cluster <- FB.combined_subset$Integrated_annotated_clusters
# Aggr <- rep(experiment,length(cluster)) 

data <- data.frame(table(Sample, cluster))
# data2 <- data.frame(table(cluster,Aggr))

# barplotAggr(seuratObj, listLabels)

# Stacked
png(file=paste0(sampleFolder,"results/5_SampleDistribution_ggplot2_ventricle_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

# png(file=paste0(sampleFolder,"CCA_FB_combined_2/3_SampleDistribution_ggplot2_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
# ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
#   geom_bar(position="fill", stat="identity", colour="white")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
# dev.off()

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split1<-data[which(data$Sample=="4V"),]
data_split2<-data[which(data$Sample=="LV"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})


data_new<-rbind(data_split1,data_split2)

# png(file=paste0(sampleFolder,"results/4_SampleDistribution_ggplot2_annot_1_",sampleName,"_adjusted.png"), width = 2000, height = 1500, res = 300)
S2<-ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red',"Turquoise")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
# dev.off()

pdf(file=paste0(sampleFolder,"results/5_SampleDistribution_ggplot2_ventricle_v2_",sampleName,"_adjusted.pdf"), width =8, height= 6)
S2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #Remove grid lines
dev.off()


################################################################################

## Sample distribution across timepoints

## Create subset first of Stromal and Stalk FBs
FB.combined_double_subset<-subset(FB.combined_subset, ident = c("Stalk Fibroblasts", "Stromal Fibroblasts"))
DimPlot(FB.combined_double_subset)

# create a dataset
FB.combined_double_subset@meta.data$Age<-factor(FB.combined_double_subset@meta.data$Age, levels = c("Embryonal", "7w", "22w", "82w"))
Sample <- FB.combined_double_subset@meta.data$Age
cluster <- FB.combined_double_subset@active.ident

data <- data.frame(table(Sample, cluster))

# Stacked
png(file=paste0(sampleFolder,"results/7_SampleDistribution_ggplot2_age_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split1<-data[which(data$Sample=="Embryonal"),]
data_split2<-data[which(data$Sample=="7w"),]
data_split3<-data[which(data$Sample=="22w"),]
data_split4<-data[which(data$Sample=="82w"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})
data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-round((x/sum(data_split3$Freq))*100,2)})
data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-round((x/sum(data_split4$Freq))*100,2)})

data_new<-rbind(data_split1,data_split2,data_split3,data_split4)

# png(file=paste0(sampleFolder,"results/4_SampleDistribution_ggplot2_annot_1_",sampleName,"_adjusted.png"), width = 2000, height = 1500, res = 300)
S2<-ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red',"Turquoise")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
# dev.off()

pdf(file=paste0(sampleFolder,"results/7_SampleDistribution_ggplot2_age_v2_",sampleName,"_adjusted.pdf"), width =8, height= 6)
S2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #Remove grid lines
dev.off()

############################################

## Also create UMAP with ventricle info (preference Daan)
U_ventricle<-DimPlot(FB.combined_subset, reduction = "umap", label = T, group.by = "Ventricle", repel = T, label.size = 4) + NoLegend()
ggsave(U_ventricle, file=paste0(sampleFolder,"results/6_UMAP_ventricle_",sampleName,".png"), height = 10, width = 12, dpi = "retina")

## Also create UMAP with age info (preference Daan)
U_age<-DimPlot(FB.combined_subset, reduction = "umap", label = T, group.by = "Age", repel = T, label.size = 4) + NoLegend()
ggsave(U_age, file=paste0(sampleFolder,"results/6_UMAP_age_",sampleName,".png"), height = 10, width = 12, dpi = "retina")

FB.combined_subset@meta.data$Age<-factor(FB.combined_subset@meta.data$Age, levels = c("Embryonal", "7w", "22w", "82w"))
U_age_v2<-DimPlot(FB.combined_subset, reduction = "umap", label = T, group.by = "Age", repel = T, label.size = 4) + NoLegend()
ggsave(U_age_v2, file=paste0(sampleFolder,"results/6_UMAP_age_v2_",sampleName,".png"), height = 10, width = 12, dpi = "retina")

## Remake UMAP for subset
U_annot<-DimPlot(FB.combined_subset, reduction = "umap", label = T, group.by = "Integrated_annotated_clusters", repel = T, label.size = 4) + NoLegend()
ggsave(U_annot, file=paste0(sampleFolder,"results/6_UMAP_annotated_",sampleName,".png"), height = 10, width = 12, dpi = "retina")

############################################
## Paper figures

FB.combined_subset@meta.data$Age<-factor(FB.combined_subset@meta.data$Age, levels = c("Embryonal", "7w", "22w", "82w"))
U_age_v2<-DimPlot(FB.combined_subset, reduction = "umap", label = T, group.by = "Age", repel = T, label.size = 4) 
pdf(file=paste0(sampleFolder,"results/6_UMAP_age_v2_",sampleName,".pdf"), width =13, height= 10)
U_age_v2
dev.off()

U_annot<-DimPlot(FB.combined_subset, reduction = "umap", label = T, group.by = "Integrated_annotated_clusters", repel = T, label.size = 4)
pdf(file=paste0(sampleFolder,"results/6_UMAP_annotated_",sampleName,".pdf"), width =14, height= 10)
U_annot
dev.off()

## Increase dot size for paper
FB.combined_subset@meta.data$Age<-factor(FB.combined_subset@meta.data$Age, levels = c("Embryonal", "7w", "22w", "82w"))
U_age_big<-DimPlot(FB.combined_subset, reduction = "umap", label = T, group.by = "Age", repel = T, label.size = 4, pt.size = 1) 
pdf(file=paste0(sampleFolder,"results/6_UMAP_age_bigger_dots_",sampleName,".pdf"), width =13, height= 10)
U_age_big
dev.off()

U_annot_big<-DimPlot(FB.combined_subset, reduction = "umap", label = T, group.by = "Integrated_annotated_clusters", repel = T, 
                 label.size = 4, pt.size = 1)
pdf(file=paste0(sampleFolder,"results/6_UMAP_annotated_bigger_dots_",sampleName,".pdf"), width =14, height= 10)
U_annot_big
dev.off()
################################################################################

######################################
##### Markers annotated clusters
########################################

### Find RNAmarkers for every Integrated cluster compared to all remaining cells, report only the positive ones
library(future)
plan("multiprocess", workers = 8)

Idents(FB.combined_subset)

RNAMarkers_integratedclus <- FindAllMarkers(FB.combined_subset, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_integratedclus$cluster)

IntegratedMarkers_integratedclus <- FindAllMarkers(FB.combined_subset, assay = "integrated", only.pos = TRUE)
table(IntegratedMarkers_integratedclus$cluster)

saveRDS(RNAMarkers_integratedclus, file=paste0(sampleFolder,"/Robjects/RNAmarkersList_integratedclus_annotated_",sampleName,".rds"))
saveRDS(IntegratedMarkers_integratedclus, file=paste0(sampleFolder,"/Robjects/IntegratedmarkersList_integratedclus_annotated_",sampleName,".rds"))

### Create list with RNA markers
totalNrIntegratedclusters_RNA<-names(table(RNAMarkers_integratedclus$cluster))
RNAmarkersList_Integratedclus<-list()

for(i in totalNrIntegratedclusters_RNA){

  tmp<-RNAMarkers_integratedclus[RNAMarkers_integratedclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  RNAmarkersList_Integratedclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Create list with integrated markers
totalNrIntegratedclusters_Integrated<-names(table(IntegratedMarkers_integratedclus$cluster))
IntegratedmarkersList_Integratedclus<-list()

for(i in totalNrIntegratedclusters_Integrated){

  tmp<-IntegratedMarkers_integratedclus[IntegratedMarkers_integratedclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  IntegratedmarkersList_Integratedclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_Integratedclus, file =paste0(sampleFolder, "results/RNAmarkersList_Integratedclus_annotated_",sampleName,".xlsx"))
write.xlsx(IntegratedmarkersList_Integratedclus, file =paste0(sampleFolder, "results/IntegratedmarkersList_Integratedclus_annotated_",sampleName,".xlsx"))

################################################################################

## Change Stalk to Base for paper (30/01/23)
levels(FB.combined_subset@meta.data[["annotated_clusters_old"]])[5]<-"Base Fibroblasts"
levels(FB.combined_subset@meta.data[["Integrated_annotated_clusters"]])[6]<-"Base Fibroblasts"
levels(FB.combined_subset@meta.data[["newClustersTmp"]])[3]<-"Base Fibroblasts"
FB.combined_subset@meta.data[["newClusters"]]<-as.factor(FB.combined_subset@meta.data[["newClusters"]])
levels(FB.combined_subset@meta.data[["newClusters"]])[20]<-"Base Fibroblasts_22w"
levels(FB.combined_subset@meta.data[["newClusters"]])[21]<-"Base Fibroblasts_7w"
levels(FB.combined_subset@meta.data[["newClusters"]])[22]<-"Base Fibroblasts_82w"
levels(FB.combined_subset@meta.data[["newClusters"]])[23]<-"Base Fibroblasts_Embryonal"
FB.combined_subset$FB_Base_signature_score1<-FB.combined_subset$FB_stalk_signature_score1
FB.combined_subset$FB_stalk_signature_score1<-NULL
Idents(FB.combined_subset)<-FB.combined_subset$Integrated_annotated_clusters

################################################################################

###############

## Markers Daan volcanoplot (Mar 2023)

Idents(FB.combined_subset)<-FB.combined_subset$newClusters

########## 2. GET MARKERS (everything!!) ##########
getDEgenes<-function(ident1, ident2){
  markersDiff <- FindMarkers(FB.combined_subset, ident.1 = ident1, ident.2 = ident2,
                             logfc.threshold = 0.0, min.pct = 0.1) #0.30, 0.10, min.diff.pct = 0.15
  #markersDiff<-markersDiff[markersDiff$p_val_adj < 0.01,]
  markersDiff<-markersDiff[order(markersDiff$avg_logFC, decreasing = T),]
  
  markersDiff$geneSymbol<-rownames(markersDiff)
  markersDiff$pct.1<-markersDiff$pct.1+0.001
  markersDiff$pct.2<-markersDiff$pct.2+0.001
  
  markersDiff<-rbind(markersDiff[markersDiff$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/pct.2*avg_logFC),
                     markersDiff[markersDiff$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/pct.1*avg_logFC))
  markersDiff<-markersDiff[order(markersDiff$score, decreasing = T),]
  return(markersDiff)
}

#### Get diff markers Stalk FBs 88w vs 7w #####
library(future)
plan("multiprocess", workers = 6)

levels(Idents(FB.combined_subset))[grep("Base",levels(Idents(FB.combined_subset)))] 
levels(Idents(FB.combined_subset))[grep("Stromal",levels(Idents(FB.combined_subset)))] 

# Clusters to compare!!!!
"Base Fibroblasts_22w"       
"Base Fibroblasts_7w"       
"Base Fibroblasts_82w"       
"Base Fibroblasts_Embryonal"

Base_82w_vs_Base_7w_volcano<-getDEgenes("Base Fibroblasts_82w","Base Fibroblasts_7w")
Base_82w_vs_Base_7w_volcano<-Base_82w_vs_Base_7w_volcano[order(Base_82w_vs_Base_7w_volcano$avg_logFC,decreasing = T),]
head(Base_82w_vs_Base_7w_volcano)
dim(Base_82w_vs_Base_7w_volcano)

Base_22w_vs_Base_7w_volcano<-getDEgenes("Base Fibroblasts_22w","Base Fibroblasts_7w")
Base_22w_vs_Base_7w_volcano<-Base_22w_vs_Base_7w_volcano[order(Base_22w_vs_Base_7w_volcano$avg_logFC,decreasing = T),]
head(Base_22w_vs_Base_7w_volcano)
dim(Base_22w_vs_Base_7w_volcano)

Base_82w_vs_Base_22w_volcano<-getDEgenes("Base Fibroblasts_82w","Base Fibroblasts_22w")
Base_82w_vs_Base_7w_volcano<-Base_82w_vs_Base_7w_volcano[order(Base_82w_vs_Base_7w_volcano$avg_logFC,decreasing = T),]
head(Base_82w_vs_Base_7w_volcano)
dim(Base_82w_vs_Base_7w_volcano)


##add to list
listDiffMarkers_volcano<-tibble::lst(Base_82w_vs_Base_7w_volcano,Base_22w_vs_Base_7w_volcano,
                                     Base_82w_vs_Base_22w_volcano)

lapply(listDiffMarkers_volcano, dim)

##Add geneSymbol in column (for the export)
listDiffMarkers_volcano<-lapply(listDiffMarkers_volcano,function(x){x<-cbind(x,'gene'=rownames(x))})
# ##Filter on adj.P-value
# listDiffMarkers_volcano<-lapply(listDiffMarkers_volcano, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDiffMarkers_volcano<-lapply(listDiffMarkers_volcano, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                                           x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
# listDiffMarkers_volcano<-lapply(listDiffMarkers_volcano, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
listDiffMarkers_volcano<-lapply(listDiffMarkers_volcano,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDiffMarkers_volcano,file=paste0(sampleFolder,"Robjects/Markers_for_our_Base_FBs_volcanoplot_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDiffMarkers_volcano, paste0(sampleFolder,"results/Age_DE_analysis/Markers_for_our_Base_FBs_volcanoplot_",sampleName,".xlsx"))

#########
library(EnhancedVolcano)

list_names_volc<-c("Base FBs 82w vs 7w","Base FBs 22w vs 7w","Base FBs 82w vs 22w")

for (k in 1:length(listDiffMarkers_volcano)) {
  E1<-EnhancedVolcano(listDiffMarkers_volcano[[k]],
                      lab = listDiffMarkers_volcano[[k]][["gene"]],
                      x = "avg_logFC",
                      y = "p_val_adj",
                      xlim = c(-3,3),
                      title = paste0(list_names_volc[[k]]),
                      pCutoff = 0.01,
                      FCcutoff = 0.25,
                      pointSize = 2.0,
                      labSize = 3.0)
  pdf(file=paste0(sampleFolder,"results/Age_DE_analysis/Volcanoplot_",names(listDiffMarkers_volcano)[k],"_",sampleName,".pdf"), width=7, height=15)
  print(E1)
  dev.off()
}

#################################################################

levels(Idents(FB.combined_subset))[grep("Stromal",levels(Idents(FB.combined_subset)))] 

# Clusters to compare!!!!
"Stromal Fibroblasts_22w"          
"Stromal Fibroblasts_7w"          
"Stromal Fibroblasts_82w"          

Stromal_82w_vs_Stromal_7w_volcano<-getDEgenes("Stromal Fibroblasts_82w","Stromal Fibroblasts_7w")
Stromal_82w_vs_Stromal_7w_volcano<-Stromal_82w_vs_Stromal_7w_volcano[order(Stromal_82w_vs_Stromal_7w_volcano$avg_logFC,decreasing = T),]
head(Stromal_82w_vs_Stromal_7w_volcano)
dim(Stromal_82w_vs_Stromal_7w_volcano)

Stromal_22w_vs_Stromal_7w_volcano<-getDEgenes("Stromal Fibroblasts_22w","Stromal Fibroblasts_7w")
Stromal_22w_vs_Stromal_7w_volcano<-Stromal_22w_vs_Stromal_7w_volcano[order(Stromal_22w_vs_Stromal_7w_volcano$avg_logFC,decreasing = T),]
head(Stromal_22w_vs_Stromal_7w_volcano)
dim(Stromal_22w_vs_Stromal_7w_volcano)

Stromal_82w_vs_Stromal_22w_volcano<-getDEgenes("Stromal Fibroblasts_82w","Stromal Fibroblasts_22w")
Stromal_82w_vs_Stromal_7w_volcano<-Stromal_82w_vs_Stromal_7w_volcano[order(Stromal_82w_vs_Stromal_7w_volcano$avg_logFC,decreasing = T),]
head(Stromal_82w_vs_Stromal_7w_volcano)
dim(Stromal_82w_vs_Stromal_7w_volcano)


##add to list
listDiffMarkers_volcano_stromal<-tibble::lst(Stromal_82w_vs_Stromal_7w_volcano,Stromal_22w_vs_Stromal_7w_volcano,
                                             Stromal_82w_vs_Stromal_22w_volcano)

lapply(listDiffMarkers_volcano_stromal, dim)

names(listDiffMarkers_volcano_stromal)<-c("Stromal_82w_vs_7w_volcano","Stromal_22w_vs_7w_volcano",
                                          "Stromal_82w_vs_22w_volcano")

##Add geneSymbol in column (for the export)
listDiffMarkers_volcano_stromal<-lapply(listDiffMarkers_volcano_stromal,function(x){x<-cbind(x,'gene'=rownames(x))})
# ##Filter on adj.P-value
# listDiffMarkers_volcano_stromal<-lapply(listDiffMarkers_volcano_stromal, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDiffMarkers_volcano_stromal<-lapply(listDiffMarkers_volcano_stromal, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                                                           x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
# listDiffMarkers_volcano_stromal<-lapply(listDiffMarkers_volcano_stromal, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
listDiffMarkers_volcano_stromal<-lapply(listDiffMarkers_volcano_stromal,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDiffMarkers_volcano_stromal,file=paste0(sampleFolder,"Robjects/Markers_for_our_Stromal_FBs_volcanoplot_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDiffMarkers_volcano_stromal, paste0(sampleFolder,"results/Age_DE_analysis/Markers_for_our_Stromal_FBs_volcanoplot_",sampleName,".xlsx"))

#########
library(EnhancedVolcano)

list_names_volc_stromal<-c("Stromal FBs 82w vs 7w","Stromal FBs 22w vs 7w","Stromal FBs 82w vs 22w")

for (k in 1:length(listDiffMarkers_volcano_stromal)) {
  E1<-EnhancedVolcano(listDiffMarkers_volcano_stromal[[k]],
                      lab = listDiffMarkers_volcano_stromal[[k]][["gene"]],
                      x = "avg_logFC",
                      y = "p_val_adj",
                      xlim = c(-3.5,3.5),
                      title = paste0(list_names_volc_stromal[[k]]),
                      pCutoff = 0.01,
                      FCcutoff = 0.25,
                      pointSize = 2.0,
                      labSize = 3.0)
  pdf(file=paste0(sampleFolder,"results/Age_DE_analysis/Volcanoplot_",names(listDiffMarkers_volcano_stromal)[k],"_",sampleName,".pdf"), width=7, height=15)
  print(E1)
  dev.off()
}

### GOEA v3: with background scRNA-seq experiment 26/04/21 ###
library(clusterProfiler)
library(org.Mm.eg.db)
# options(connectionObserver = NULL) #If issue loading package org.MM.eg.db, need to run this line!!!

## All genes object as background
Genes_scRnaseq<-rownames(FB.combined_subset@assays$RNA@counts)

## Filter for DE genes
listDiffMarkers_volcano_stromal_filtered<-listDiffMarkers_volcano_stromal
##Filter on adj.P-value and logFC
listDiffMarkers_volcano_stromal_filtered<-lapply(listDiffMarkers_volcano_stromal_filtered, function(x){dplyr::filter(x, p_val_adj<0.01)})
listDiffMarkers_volcano_stromal_filtered<-lapply(listDiffMarkers_volcano_stromal_filtered, function(x){dplyr::filter(x, abs(avg_logFC)>0.25)})

## Stromal FBs
Test1_background<-enrichGO(
  listDiffMarkers_volcano_stromal_filtered$Stromal_82w_vs_22w_volcano$geneSymbol,
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = Genes_scRnaseq, ##Full background
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D1_background<-dotplot(Test1_background, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Most correct and best results with full background!!!!!!!! Use those!!

### Save results ###
pdf(file=paste0(sampleFolder,"results/Age_DE_analysis/Dotplot_FB1s_82w_vs_22w_full_background_",sampleName,".pdf"), width = 10, height = 10)
D1_background
dev.off()

write.xlsx(Test1_background@result,file=paste0(sampleFolder,"results/Age_DE_analysis/EnrichGO_FB1s_82w_vs_22w_full_background_",sampleName,".xlsx"))

saveRDS(Test1_background,file=paste0(sampleFolder,"Robjects/EnrichGO_FB1s_82w_vs_22w_full_background_",sampleName,".xlsx"))

##############################################################################################################

## Save subset
saveRDS(FB.combined_subset, file = paste0(sampleFolder,"Robjects/seuratObj_subset_",sampleName,".rds"))

##### Read object
FB.combined_subset<- readRDS(file = paste0(sampleFolder,"Robjects/seuratObj_subset_",sampleName,".rds"))
