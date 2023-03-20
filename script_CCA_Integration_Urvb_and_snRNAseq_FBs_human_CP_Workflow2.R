## Integrate snRNAseq Human Covid dataset and scRNAseq Urvb
## Convert Human gene symbols to mouse gene symbols
## Try integration via Seurat (CCA or RPCA??)

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

# pdf(file = paste0(sampleFolder,"Dimplot__",sampleName,"_both_workflows.pdf"), width = 20, height = 15)
# (p1 + p2) / (p3 + p4)
# dev.off()

FeaturePlot(FB.combined_2, features = "Igfbp6", min.cutoff = "q10", max.cutoff = "q90")
FB.combined_2$NewClusters_combo<-paste0(FB.combined_2@active.ident,"_",FB.combined_2$Species)
table(FB.combined_2$NewClusters_combo)

# saveRDS(FB.combined_2, file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))

####################################################################################################################

# ## Perform RPCA instead # https://satijalab.org/seurat/articles/integration_large_datasets.html
# 
# ## Create list of objects:
# Spleen.list<-list(seuratObj_SanVic1,seuratObj_SanVic2_WT,seuratObj_SanVic2_DKO,seuratObj_Laurien,seuratObj_Jessica)
# 
# ## Findvariablefeatures on RNA slot of objects:
# Spleen.list <- lapply(X = Spleen.list, FUN = function(x) {
#   x <- FindVariableFeatures(x, verbose = T)
# })
# 
# ## select features for downstream integration, and run PCA on each object in the list, 
# ## which is required for running the alternative reciprocal PCA workflow
# 
# features <- SelectIntegrationFeatures(object.list = Spleen.list)
# Spleen.list <- lapply(X = Spleen.list, FUN = function(x) {
#   x <- ScaleData(x, features = features, verbose = T)
#   x <- RunPCA(x, features = features, verbose = T)
# })
# 
# ##Choose best dataset to be the reference: choose the KI dataset which wasn't sorted and doesn't have a strong stimulus!
# anchors <- FindIntegrationAnchors(object.list = Spleen.list, reference = c(4), reduction = "rpca",
#                                   dims = 1:50)
# Spleen.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
# 
# ## Perform further analysis
# Spleen.integrated <- ScaleData(Spleen.integrated, verbose = T)
# Spleen.integrated <- RunPCA(Spleen.integrated, verbose = T)
# Spleen.integrated <- RunUMAP(Spleen.integrated, dims = 1:50)
# 
# ## Save files
# saveRDS(anchors, "CCA_FB_combined_2/Anchors.rds")
# saveRDS(Spleen.integrated, "CCA_FB_combined_2/Spleen.integrated.rds")
# 
# ## Visualize
# D1<- DimPlot(Spleen.integrated, group.by = "orig.ident")
# 
# D2<-DimPlot(Spleen.integrated, reduction = "umap", group.by = "annotated_clusters", label = TRUE, repel = TRUE) +
#   NoLegend()
# 
# pdf(file = paste0(sampleFolder,"Dimplot_UMAP_",sampleName,"_annotation.pdf"),
#     height = 15, width = 30)
# D1 + D2
# dev.off()
# 
# D3<-DimPlot(Spleen.integrated, ncol=7, reduction = "umap", split.by = "orig.ident", group.by = "annotated_clusters", label = TRUE, repel = TRUE) +
#   NoLegend()
# 
# pdf(file = paste0(sampleFolder,"Dimplot_UMAP_",sampleName,"_split.pdf"),
#     height = 20, width = 40)
# D3
# dev.off()

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
  
  # ##### Create tSNE plot
  # FB.combined_2 <- RunTSNE(object = FB.combined_2, dims = dimsToUse, assay = "RNA", reduction = "RNA_pca", 
  #                      reduction.name = "RNA_tsne", reduction.key = "RNATSNE_")
  # tsnePlot<-DimPlot(FB.combined_2, reduction = "RNA_tsne", label=T, label.size = 8)
  # tsnePlotSplit<-DimPlot(FB.combined_2, reduction = "RNA_tsne", label=F, group.by="ident", pt.size = 2)
  # 
  # ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
  #        file=paste0(sampleFolder,"Plots/RNA/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  # 
  
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
Perplexity_UMAP<-20

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

##### Save object
saveRDS(FB.combined_2, file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))
# saveRDS(diagnostics, file=paste0("CCA_FB_combined_2/diagnostics_",sampleName,".rds"))

##### Read object
FB.combined_2<- readRDS(file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))
# diagnostics<- readRDS(file=paste0("CCA_FB_combined_2/diagnostics_",sampleName,".rds"))

########################################################################

### Find RNAmarkers for every Integrated cluster compared to all remaining cells, report only the positive ones
library(future)
plan("multiprocess", workers = 8)

RNAMarkers_RNAclus <- FindAllMarkers(FB.combined_2, assay = "integrated", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"CCA_FB_combined_2/Robjects/IntegratedmarkersList_integratedclus_",sampleName,".rds"))

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
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  RNAmarkersList_RNAclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(RNAmarkersList_RNAclus)<-paste0("Integratedcluster",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "CCA_FB_combined_2/IntegratedmarkersList_Integratedclus_",sampleName,".xlsx"))


######################################
##### Markers annotated clusters
########################################

## Annotate clusters reference
DimPlot(FB.combined_2, reduction = "umap", group.by = "Integrated_RNA_clusters", label = TRUE, repel = TRUE)

# FB.combined_2@meta.data$Integrated_RNA_clusters<-factor(FB.combined_2@meta.data$Integrated_RNA_clusters,levels=sort(as.numeric(levels(FB.combined_2@meta.data$Integrated_RNA_clusters)))) #reorder levels

FB.combined_2@meta.data$Integrated_annotated_clusters <- FB.combined_2@meta.data$Integrated_RNA_clusters
levels(FB.combined_2@meta.data$Integrated_annotated_clusters) <- c("Stromal FBs Mouse","Stromal FBs Human","Stromal FBs Human","Stromal FBs Mouse",
                                                                   "Stromal FBs Mouse","Stromal FBs Human","Other FBs Human 1","Other FBs Human 2",
                                                                   "Stromal FBs Mouse","Hsp+ Stromal FBs Human","Stalk FBs","Other FBs Human 2",
                                                                   "Other FBs Human 3 ", "ABCs")
U_annot<-DimPlot(FB.combined_2, reduction = "umap", label = T, group.by = "Integrated_annotated_clusters", repel = T, label.size = 4) + NoLegend()
ggsave(U_annot, file=paste0(sampleFolder,"CCA_FB_combined_2/UMAP_annotated_cluster1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

#############################################################################################################
Idents(FB.combined_2)<-FB.combined_2$Integrated_annotated_clusters

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(FB.combined_2) <- "RNA"
Stalk.markers <- FindConservedMarkers(FB.combined_2, ident.1 = "Stalk FBs", grouping.var = "Species", verbose = FALSE)
head(Stalk.markers)

## Save results
saveRDS(Stalk.markers, file=paste0(sampleFolder,"CCA_FB_combined_2/Robjects/Conserved_stalk_markers_",sampleName,".rds"))
write.xlsx(Stalk.markers, file =paste0(sampleFolder, "CCA_FB_combined_2/Conserved_stalk_markers_",sampleName,".xlsx"), row.names = T)

## Get top markers
Stalk.markers.top<-Stalk.markers[which(Stalk.markers$Mouse_avg_log2FC > 1 & Stalk.markers$Human_avg_log2FC >1),]
Stalk.markers.top<-Stalk.markers.top[order(Stalk.markers.top$Mouse_avg_log2FC, decreasing = T),]

## Plot out data
DefaultAssay(FB.combined_2) <- "RNA"
F1 <- FeaturePlot(FB.combined_2, features = head(rownames(Stalk.markers.top), n =9), pt.size = 1,  min.cutoff = 'q2', max.cutoff = 'q98', ncol = 3)
F2 <- FeaturePlot(FB.combined_2, features = head(rownames(Stalk.markers.top), n =9), order = T,  pt.size = 1,  min.cutoff = 'q2', max.cutoff = 'q98', ncol = 3)

## Make subset of Stalk data
FB.Stalk_only<-subset(FB.combined_2, idents = "Stalk FBs")

## Try dotplots
markers.to.plot <- c("Igfbp6", "Cldn11", "S100a6", "Nbl1", "Vim")

D1 <- DotPlot(FB.Stalk_only, features = markers.to.plot, cols = c("blue", "red"), assay = "RNA", split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot Stalk Only")

D2 <- DotPlot(FB.combined_2, features = markers.to.plot, cols = c("blue", "red"), assay = "RNA", split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot All")

FB.Stalk_only$NewClusters_combo_annotated<-paste0(FB.Stalk_only$Integrated_annotated_clusters,"_",FB.Stalk_only$Species)
FB.combined_2$NewClusters_combo_annotated<-paste0(FB.combined_2$Integrated_annotated_clusters,"_",FB.combined_2$Species)

D3 <- DotPlot(FB.Stalk_only, group.by = "NewClusters_combo_annotated", assay = "RNA", features = markers.to.plot) +
  RotatedAxis() + ggtitle("Dotplot Stalk Only with scale")

D4 <- DotPlot(FB.combined_2, group.by = "NewClusters_combo_annotated", assay = "RNA", features = markers.to.plot) +
  RotatedAxis() + ggtitle("Dotplot All with scale")

## Save plots
pdf(file=paste0(sampleFolder,"CCA_FB_combined_2/1_Featureplots_Stalk_Markers_non_ordered_",sampleName,".pdf"), width = 20, height = 15)
F1
dev.off()

pdf(file=paste0(sampleFolder,"CCA_FB_combined_2/1_Featureplots_Stalk_Markers_ordered_",sampleName,".pdf"), width = 20, height = 15)
F2
dev.off()

pdf(file=paste0(sampleFolder,"CCA_FB_combined_2/1_Dotplots_Stalk_Markers_",sampleName,".pdf"), width = 10, height = 8)
D1
D2
D3
D4
dev.off()

## Second version of markers
# Filter on adj P-value for human and mouse and order according to Max P-value -> top 5

Top_markers_v2<- c("Igfbp6","Clu","S100a6", "Itm2a","Islr")

D2_v2 <- DotPlot(FB.combined_2, features = Top_markers_v2, assay = "RNA", cols = c("blue", "red"), split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot All v2")
D4_v2 <- DotPlot(FB.combined_2, group.by = "NewClusters_combo_annotated", assay = "RNA", features = Top_markers_v2) +
  RotatedAxis() + ggtitle("Dotplot All with scale v2")

pdf(file=paste0(sampleFolder,"CCA_FB_combined_2/1_Dotplots_Stalk_Markers_v2_",sampleName,".pdf"), width = 10, height = 8)
D2_v2
D4_v2
dev.off()

VlnPlot(FB.combined_2, split.by = "Species", assay = "integrated", features = Top_markers_v2) +
  RotatedAxis() + ggtitle("Violinplot All with scale v2")

DotPlot(FB.combined_2, features = Top_markers_v2, assay = "integrated", cols = c("blue", "red"), split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot All v2")

################################################################################################################

## Check Igfbp6 and Cldn11 expression in regressed object

# Featureplots on full object
# Module score Seurat
library(viridis)

FB_stalk_signature <- list(c("Cldn11","Igfbp6"))

FB.combined_2 <- AddModuleScore(object = FB.combined_2, assay = "RNA", features = FB_stalk_signature, name = "FB_stalk_signature_score")
F2 <- FeaturePlot(object = FB.combined_2, features = "FB_stalk_signature_score1", order = T)  + scale_color_viridis(option = "C")
pdf(file = paste0(sampleFolder,"CCA_FB_combined_2/2_UMAP_Featureplot_modulescore_Cldn11_Igfbp6_",sampleName,".pdf"), width = 15, height = 10)
F2
dev.off()

DefaultAssay(FB.combined_2) <- "RNA"
Idents(FB.combined_2)<-FB.combined_2@meta.data$Integrated_annotated_clusters
F2.5 <- FeaturePlot(object = FB.combined_2, features = "FB_stalk_signature_score1", order = T, label = T, repel = T, label.size = 3, cols = c("Yellow","Red"))
pdf(file = paste0(sampleFolder,"CCA_FB_combined_2/2_UMAP_Featureplot_modulescore_Cldn11_Igfbp6_labeled_",sampleName,".pdf"), width = 15, height = 10)
F2.5
dev.off()

# Extra featureplots
F1 <- FeaturePlot(object = FB.combined_2, features = c("Cldn11","Igfbp6","Dpep1","Alpl","Cdh11"), order = T, pt.size = 1, cols = c("Yellow","Red"), combine = F)
pdf(file = paste0(sampleFolder,"CCA_FB_combined_2/2_UMAP_Featureplots_FB_markers_RNA_",sampleName,".pdf"), width = 15, height = 10)
F1
dev.off()

# # Create new clusters: split on source
# FB.combined_2@meta.data$newClustersTmp<-FB.combined_2@meta.data$Integrated_annotated_clusters
# FB.combined_2@meta.data$newClusters<-paste0(FB.combined_2@meta.data$newClustersTmp,"_",FB.combined_2@meta.data$Species)
# head(FB.combined_2@meta.data)

# Subset stalk and violin plot
Idents(FB.combined_2)<-FB.combined_2@meta.data$NewClusters_combo_annotated
FB.stalkandstromal<-subset(FB.combined_2, ident = c("Stromal FBs Mouse_Mouse","ABCs_Mouse","Stalk FBs_Mouse","Stromal FBs Human_Human","Stalk FBs_Human"))

Idents(FB.stalkandstromal)<-factor(FB.stalkandstromal@active.ident, levels = levels(Idents(FB.stalkandstromal))[c(2,5,3,4,1)])
levels(Idents(FB.stalkandstromal))<-c("ABCs_Mouse","Stalk FBs_Human","Stalk FBs_Mouse","Stromal FBs_Human","Stromal FBs_Mouse")
DimPlot(FB.stalkandstromal)

library(RColorBrewer)
ColorSet<-brewer.pal(n = 6, name = "Paired")[2:6]

V1<-VlnPlot(FB.stalkandstromal, features = c("Igfbp6","Cldn11","Cdh11","Alpl"), assay = "RNA", ncol = 4, cols = ColorSet) & theme(plot.margin = margin(0.5,0.5,0.5,1, "cm"))

pdf(file=paste0(sampleFolder,"CCA_FB_combined_2/4_Test_violinplot_RNA_",sampleName,".pdf"), width = 25, height= 10)
V1
dev.off()

## UMAP with ventricle info (Daan prefers that!)
# Update Ventricle metadata (all human data is LV!!!)
FB.combined_2$Ventricle<-FB.combined_2$orig.ident

FB.combined_2@meta.data[which(FB.combined_2@meta.data$Ventricle == "RVD1_LpsNegFour" |
                              FB.combined_2@meta.data$Ventricle == "RVD5_Y4V" | FB.combined_2@meta.data$Ventricle == "RVD7_O4V"),
                      "Ventricle"]<-"4V"

FB.combined_2@meta.data[which(FB.combined_2@meta.data$Ventricle == "ct" | FB.combined_2@meta.data$Ventricle == "RVD2_LpsNegLat" |
                                FB.combined_2@meta.data$Ventricle == "cv" | FB.combined_2@meta.data$Ventricle == "RVD6_YLV" |
                                FB.combined_2@meta.data$Ventricle == "flu" | FB.combined_2@meta.data$Ventricle == "RVD8_OLV"),
                        "Ventricle"]<-"LV"

U_ventricle<-DimPlot(FB.combined_2, reduction = "umap", label = T, group.by = "Ventricle", repel = T, label.size = 4) + NoLegend()
ggsave(U_ventricle, file=paste0(sampleFolder,"CCA_FB_combined_2/UMAP_ventricle_",sampleName,".png"), height = 10, width = 15, dpi = "retina")


################################################################################

## Sample distribution across clusters

# create a dataset
Sample <- FB.combined_2@meta.data$Species
cluster <- FB.combined_2$Integrated_annotated_clusters
# Aggr <- rep(experiment,length(cluster)) 

data <- data.frame(table(Sample, cluster))
# data2 <- data.frame(table(cluster,Aggr))

# barplotAggr(seuratObj, listLabels)

# Stacked
png(file=paste0(sampleFolder,"CCA_FB_combined_2/3_SampleDistribution_ggplot2_annot_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"CCA_FB_combined_2/3_SampleDistribution_ggplot2_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
dev.off()

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split1<-data[which(data$Sample=="Human"),]
data_split2<-data[which(data$Sample=="Mouse"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})

data_new<-rbind(data_split1,data_split2)

# png(file=paste0(sampleFolder,"results/4_SampleDistribution_ggplot2_annot_1_",sampleName,"_adjusted.png"), width = 2000, height = 1500, res = 300)
S2<-ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red',"Turquoise")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
# dev.off()

pdf(file=paste0(sampleFolder,"CCA_FB_combined_2/3_SampleDistribution_ggplot2_annotated_",sampleName,"_adjusted.pdf"), width =8, height= 6)
S2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #Remove grid lines
dev.off()
