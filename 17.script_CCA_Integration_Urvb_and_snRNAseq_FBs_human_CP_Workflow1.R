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
sampleName <- "CCA_integrated_Human_and_Mouse_CP_W1"

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
###### Workflow 1 ##### -> No extra subsetting of seuratobjects (different sets of genes)
#######################

## Find Variable Features first for Human data
seuratObj_Human_subset_converted <- FindVariableFeatures(seuratObj_Human_subset_converted, assay = "RNA", 
                                                         selection.method = "vst", nfeatures = 2000)

# select features that are repeatedly variable across datasets for integration
FB_list <- list(seuratObj_Mouse, seuratObj_Human_subset_converted)
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

## Add metadata column (Species)
FB.combined@meta.data$Species<-"Human"
FB.combined@meta.data[which(FB.combined@meta.data$Lab == "Vandenbroucke"),"Species"] <- "Mouse"

# Visualization
p1 <- DimPlot(FB.combined, reduction = "umap", group.by = "Species") + ggtitle("Species_workflow1")
p2 <- DimPlot(FB.combined, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Numbered_annotation_workflow1")
p1 + p2

FeaturePlot(FB.combined, features = "Igfbp6", min.cutoff = "q2", max.cutoff = "q98")
FB.combined$NewClusters_combo<-paste0(FB.combined@active.ident,"_",FB.combined$Species)
table(FB.combined$NewClusters_combo)

# saveRDS(FB.combined, file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))

####################################################################################################################

########## Via PCelbowplot ##########
ElbowPlot(object = FB.combined, ndims = 30)

## Perform clustering
dimsToTry<-c(seq(10,30,by=5))

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
         file=paste0(sampleFolder,"CCA_FB_combined/UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
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
  FB.combined <- FindNeighbors(object = FB.combined, reduction = "pca", dims = dimsToUse)
  FB.combined <- FindClusters(object = FB.combined,  graph.name = "integrated_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  FB.combined <- RunUMAP(FB.combined, dims = dimsToUse, n_neighbors = 30, assay = "integrated", reduction ="pca",
                               reduction.name = "umap", reduction.key = "integratedUMAP_")
  umapPlot<-DimPlot(FB.combined, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(FB.combined, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(sampleFolder,"CCA_FB_combined/UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

names(FB.combined)

### Clustering: trying out clusTree

Perplexity<-20
Resolution<-0.8
Perplexity_UMAP<-20

FB.combined <- FindNeighbors(object = FB.combined, reduction = "pca", dims = 1:Perplexity)
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  FB.combined <- FindClusters(object = FB.combined,  graph.name = "integrated_snn", resolution = res)
}

library("clustree")
pdf(file=paste0(sampleFolder,"CCA_FB_combined/Clustree.pdf"))
clustree(FB.combined, prefix = "integrated_snn_res.")
dev.off()

# 0.8 seems reasonable
# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
FB.combined$Integrated_RNA_clusters <- FB.combined$integrated_snn_res.0.8
Idents(FB.combined)<-FB.combined@meta.data$Integrated_RNA_clusters

##### Save object
saveRDS(FB.combined, file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))
# saveRDS(diagnostics, file=paste0("CCA_FB_combined/diagnostics_",sampleName,".rds"))

##### Read object
FB.combined<- readRDS(file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))
# diagnostics<- readRDS(file=paste0("CCA_FB_combined/diagnostics_",sampleName,".rds"))

########################################################################

### Find RNAmarkers for every Integrated cluster compared to all remaining cells, report only the positive ones
library(future)
plan("multiprocess", workers = 8)

RNAMarkers_RNAclus <- FindAllMarkers(FB.combined, assay = "integrated", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"CCA_FB_combined/Robjects/IntegratedmarkersList_integratedclus_",sampleName,".rds"))

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
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "CCA_FB_combined/IntegratedmarkersList_Integratedclus_",sampleName,".xlsx"))


######################################
##### Markers annotated clusters
########################################

## Annotate clusters reference
DimPlot(FB.combined, reduction = "umap", group.by = "Integrated_RNA_clusters", label = TRUE, repel = TRUE)

# FB.combined@meta.data$Integrated_RNA_clusters<-factor(FB.combined@meta.data$Integrated_RNA_clusters,levels=sort(as.numeric(levels(FB.combined@meta.data$Integrated_RNA_clusters)))) #reorder levels

FB.combined@meta.data$Integrated_annotated_clusters <- FB.combined@meta.data$Integrated_RNA_clusters
levels(FB.combined@meta.data$Integrated_annotated_clusters) <- c("Stromal FBs Human","Stromal FBs Mouse","Stromal FBs Mouse","Stromal FBs Mouse",
                                                                  "Stromal FBs Human","Stromal FBs Human","Stromal FBs Human","Other FBs Human 1",
                                                                  "Stromal FBs Mouse","Stalk FBs","Hsp+ Stromal FBs Human",
                                                                  "Other FBs Human 2 ","Other FBs Human 1", "CPE contamination", "ABCs")
U_annot<-DimPlot(FB.combined, reduction = "umap", label = T, group.by = "Integrated_annotated_clusters", repel = T, label.size = 4) + NoLegend()
ggsave(U_annot, file=paste0(sampleFolder,"CCA_FB_combined/UMAP_annotated_cluster_",sampleName,".png"), height = 10, width = 15, dpi = "retina")


##################################################################################################################

################################################################################

## Sample distribution across clusters

# create a dataset
Sample <- FB.combined@meta.data$Species
cluster <- FB.combined$Integrated_annotated_clusters
# Aggr <- rep(experiment,length(cluster)) 

data <- data.frame(table(Sample, cluster))
# data2 <- data.frame(table(cluster,Aggr))

# barplotAggr(seuratObj, listLabels)

# Stacked
png(file=paste0(sampleFolder,"CCA_FB_combined/3_SampleDistribution_ggplot2_annot_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"CCA_FB_combined/3_SampleDistribution_ggplot2_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
dev.off()

#############################################################################################################
Idents(FB.combined)<-FB.combined$Integrated_annotated_clusters

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(FB.combined) <- "RNA"
Stalk.markers <- FindConservedMarkers(FB.combined, ident.1 = "Stalk FBs", grouping.var = "Species", verbose = FALSE)
head(Stalk.markers)

## Save results
saveRDS(Stalk.markers, file=paste0(sampleFolder,"CCA_FB_combined/Robjects/Conserved_stalk_markers_",sampleName,".rds"))
write.xlsx(Stalk.markers, file =paste0(sampleFolder, "CCA_FB_combined/Conserved_stalk_markers_",sampleName,".xlsx"), row.names = T)

## Get top markers
Stalk.markers.top<-Stalk.markers[which(Stalk.markers$Mouse_avg_log2FC > 1 & Stalk.markers$Human_avg_log2FC >1),]
Stalk.markers.top<-Stalk.markers.top[order(Stalk.markers.top$Mouse_avg_log2FC, decreasing = T),]

## Plot out data
DefaultAssay(FB.combined) <- "RNA"
F1 <- FeaturePlot(FB.combined, features = head(rownames(Stalk.markers.top), n =5), pt.size = 1,  min.cutoff = 'q2', max.cutoff = 'q98', ncol = 3)
F2 <- FeaturePlot(FB.combined, features = head(rownames(Stalk.markers.top), n =5), order = T,  pt.size = 1,  min.cutoff = 'q2', max.cutoff = 'q98', ncol = 3)

## Make subset of Stalk data
FB.Stalk_only<-subset(FB.combined, idents = "Stalk FBs")

## Try dotplots
markers.to.plot <- c("Igfbp6", "Ptgds", "S100a6", "Nbl1", "Vim","Cldn11")

D1 <- DotPlot(FB.Stalk_only, features = markers.to.plot, cols = c("blue", "red"), assay = "RNA", split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot Stalk Only")

D2 <- DotPlot(FB.combined, features = markers.to.plot, cols = c("blue", "red"), assay = "RNA", split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot All")

FB.Stalk_only$NewClusters_combo_annotated<-paste0(FB.Stalk_only$Integrated_annotated_clusters,"_",FB.Stalk_only$Species)
FB.combined$NewClusters_combo_annotated<-paste0(FB.combined$Integrated_annotated_clusters,"_",FB.combined$Species)

D3 <- DotPlot(FB.Stalk_only, group.by = "NewClusters_combo_annotated", assay = "RNA", features = markers.to.plot) +
  RotatedAxis() + ggtitle("Dotplot Stalk Only with scale")

D4 <- DotPlot(FB.combined, group.by = "NewClusters_combo_annotated", assay = "RNA", features = markers.to.plot) +
  RotatedAxis() + ggtitle("Dotplot All with scale")

## Save plots
pdf(file=paste0(sampleFolder,"CCA_FB_combined/Featureplots_Stalk_Markers_non_ordered_",sampleName,".pdf"), width = 20, height = 15)
F1
dev.off()

pdf(file=paste0(sampleFolder,"CCA_FB_combined/Featureplots_Stalk_Markers_ordered_",sampleName,".pdf"), width = 20, height = 15)
F2
dev.off()

pdf(file=paste0(sampleFolder,"CCA_FB_combined/Dotplots_Stalk_Markers_",sampleName,".pdf"), width = 10, height = 8)
D1
D2
D3
D4
dev.off()

## Second version of markers
# Filter on adj P-value for human and mouse and order according to Max P-value -> top 5 (Without Rik of Rpl genes!)

Top_markers_v2<- c("Igfbp6","Clu","Vim", "Islr","Tpt1")

D2_v2 <- DotPlot(FB.combined, features = Top_markers_v2, cols = c("blue", "red"), assay = "RNA", split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot All v2")
D4_v2 <- DotPlot(FB.combined, group.by = "NewClusters_combo_annotated", assay = "RNA", features = Top_markers_v2) +
  RotatedAxis() + ggtitle("Dotplot All with scale v2")

pdf(file=paste0(sampleFolder,"CCA_FB_combined/Dotplots_Stalk_Markers_v2_",sampleName,".pdf"), width = 10, height = 8)
D2_v2
D4_v2
dev.off()
