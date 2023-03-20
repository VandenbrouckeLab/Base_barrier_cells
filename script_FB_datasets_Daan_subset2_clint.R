## Follow-up FB analysis with atlas datasets to check origin of FB2s
## Subsetted to contain only the FBs (no prolif FBs in subset2!!)

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('openxlsx')

# BiocManager::install("impute") version = "3.8")
# BiocManager::install("preprocessCore") version = "3.8")
# BiocManager::install("GO.db") version = "3.8")
# BiocManager::install("AnnotationDbi") version = "3.8")
# devtools::install_github('dambi/DisneyTool@seuratv3', host="github.ugent.be/api/v3", auth_token = 'e5ca75c8c2f815aa7f1195cb0b6b6a3190064707')
#cutils::download.file("https://github.ugent.be/api/v3/repos/dambi/DisneyTool/tarball/master", destfile = "test.zip", method = "curl")
library('DisneyTools')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")

sampleName <- "Merge_FB_datasets_subset2" #Change for this analysis!!!
sampleFolder<-paste0("Merge_subset2","/")

##add some subfolders
dir.create(paste0(sampleFolder,"results_merge_subset2"))
dir.create(paste0(sampleFolder,"results_merge_subset2/QC"))
dir.create(paste0(sampleFolder,"results_merge_subset2/Robjects"))

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

##### Read object
seuratObj <- readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_",sampleName,"_harmony_RNA.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,"_harmony_RNA.rds"))

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_",sampleName,"_harmony_RNA.rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,"_harmony_RNA.rds"))

##### Create new clusters
##new clusters
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

Idents(seuratObj)<-seuratObj@meta.data$harmony_clusters
seuratObj@meta.data$sliced_clusters<-seuratObj@meta.data$harmony_clusters
# Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters

## Check 1.0 clustering:
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "RNA_snn_res.1", label.size = 8)

## Use 1.0 res clustering!!  
## Doublets5 = cluster 10 in res 1.0!!!
Idents(seuratObj)<-seuratObj@meta.data$RNA_snn_res.1
seuratObj@meta.data$sliced_clusters<-seuratObj@meta.data$RNA_snn_res.1

## Save new clustering
seuratObj@meta.data$sliced_clusters <- seuratObj@active.ident #Sliced clustering
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident #Sliced clustering

################################################################################
########## PLOTS
################################################################################

## New function
drawUMI_mitoPlot_new<-function(coordsTable, reductionType, clusterMatrix, columnName, titleInfo){
  
  columnNr<-which(colnames(clusterMatrix)==columnName)
  
  p <- ggplot()+
    geom_point(aes(x=tSNE_1,y=tSNE_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20) 
  
  if(reductionType=="umap"){
    p <- ggplot()+
      geom_point(aes(x=UMAP_1,y=UMAP_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20) 
  }
  
  p<-p +
    scale_colour_gradientn(colours = c("darkblue","cyan","green","yellow","orange","darkred")) +
    ggtitle(paste0(titleInfo," (",reductionType,")")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  
  return(p)
}

########## UMI plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'nCount_RNA',"UMI")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'nCount_RNA',"UMI")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_merge_subset2/QC/1_UMI_",sampleName,".png"), width = 20)

########## mito.genes plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'subsets_Mito_percent',"mito")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'subsets_Mito_percent',"mito")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_merge_subset2/QC/2_percMito_",sampleName,".png"), width = 20)


########## PCA plot ##########
# pdf(file=paste0(sampleFolder,"results_merge_subset2/QC/13a_PCA_",sampleName,".pdf"), width=10)
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,3))
# dev.off()

#RNA clusters
dir.create(paste0(sampleFolder,"results_merge_subset2/Annotation"))

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/1_annotation_Color_RNA_clusters_on_harmony_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$harmony_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$harmony_clusters==i),])))
  C1<-C1+ggtitle(paste0("Harmony_cluster_",i))
  print(C1)
}
dev.off()

################################################################################
########## CHECK DE GENES
################################################################################
dir.create(paste0(sampleFolder,"results_merge_subset2/Feature_plots"))

##### Epithelial marker ->
Features<-c("Otx2", "Ttr")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"CPE"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Endothelial marker ->
Features<-c("Pecam1", "Flt1","Plvap")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"EC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Vascular associated marker ->
Features<-c("Pdgfrb", "Mylk","Myh11","Tagln")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"VAC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Fibroblast marker ->
Features<-c("Dcn", "Col1a1","Dpep1","Igfbp6")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"FB"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Pdgfra", "Pdgfrb","Lum","Acta2","Rgs5")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"FB_bis"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Neuronal + glial cell marker ->
Features<-c("Tubb3","Slc1a3", "Fabp7","Olig1") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Neuronal_and_glial"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Mitotic cell marker ->
Features<-c("Birc5", "Mki67")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Mitotic_cells"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

#########################################################################################################################################

########################################
##### all clusters vs all clusters
########################################

# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

dir.create(paste0(sampleFolder,"results_merge_subset2/Marker_lists"))

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
RNAMarkers_RNAclus <- FindAllMarkers(seuratObj, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_RNAclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerRNAcluster']]<-paste0(table(RNAMarkers_RNAclus$cluster)," RNA markers for RNA cluster ",rownames(table(RNAMarkers_RNAclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset2/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
names(RNAmarkersList_RNAclus)<-paste0("RNAcluster",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results_merge_subset2/Marker_lists/RNAmarkersList_RNAclus_",sampleName,".xlsx"))

########################################
dir.create(paste0(sampleFolder,"results_merge_subset2/Heatmaps"))

# ## Load in markers
# RNAMarkers_RNAclus<- readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_RNAclus_",sampleName,".rds"))
#
# ## Perform on a subset -> better view of smaller clusters!!
# seuratObj.small <- subset(seuratObj, downsample = 500)
#
# ## Heatmap RNA markers on RNA clusters
# top10 <- RNAMarkers_RNAclus %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# D1<-DoHeatmap(seuratObj.small, features = top10$gene, group.by = "RNA_clusters") + NoLegend()
# ggsave(D1, file=paste0(sampleFolder, "results_merge_subset2/Heatmaps/Heatmap_RNAmarkersList_RNAclus_",sampleName,".png"),
#        height = 20, width = 12, dpi = "retina")
#
# pdf(file=paste0(sampleFolder, "results_merge_subset2/Heatmaps/Heatmap_RNAmarkersList_RNAclus_",sampleName,".pdf"),
#     height = 20, width = 25)
# DoHeatmap(seuratObj.small, features = top10$gene, group.by = "RNA_clusters") + NoLegend()
# dev.off()

########################################################################################################################

## Check other annotation

Idents(seuratObj)<- seuratObj@meta.data$orig.ident

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/1_annotation_Color_datasets_",sampleName,".pdf"), width = 15)
for (i in levels(Idents(seuratObj))) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = i))
  C1<-C1+ggtitle(i)
  print(C1)
}
dev.off()

Idents(seuratObj)<-seuratObj@meta.data$New_clusters #Revert

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/1_annotation_Color_old_idents_",sampleName,".pdf"), width = 15)
for (i in levels(Idents(seuratObj))) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = i))
  C1<-C1+ggtitle(i)
  print(C1)
}
dev.off()

Idents(seuratObj)<-seuratObj@meta.data$harmony_clusters #Revert

U_old_annot<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "New_clusters", label.size = 3)
ggsave(U_old_annot, file=paste0(sampleFolder,"results_merge_subset2/Annotation/1_UMAP_old_annot_",sampleName,".png"), height = 8, width = 15, dpi = "retina")

## Quantify contamination (400-500 cells)
nrow(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Doublets_1"),])
nrow(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Doublets_2"),])
nrow(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Doublets_3"),])
nrow(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Doublets_4"),])
nrow(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Doublets_5"),])
nrow(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Epithelial cells"),])
nrow(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Endothelial cells"),])
nrow(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "NK cells"),])
nrow(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Ccr7+ cDCs"),])

colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "meninges 4 (dura)"),]))
colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "meninges prolif. 1"),]))
colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "meninges prolif. 2"),]))
colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "meninges 1"),]))
colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "meninges 2"),]))
colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "meninges 3 (arachnoid)"),]))
colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "VLMC2"),]))
colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "VLMC1"),]))


FeaturePlot(object = seuratObj, features = "Pdgfra", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Pdgfrb", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Dcn", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Lum", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Dpep1", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Igfbp6", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Col1a1", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)


FeaturePlot(object = seuratObj, features = "Col18a1", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "S100a6", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Ngfr", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Aldh1a2", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)

DimPlot(seuratObj, reduction = "umap", label = T, group.by = "sample_origin", label.size = 4)

########################################################################################################################

########################################
##### RNA clusters post-slice
########################################
########################################
seuratObj@meta.data$sliced_clusters<- factor(seuratObj@meta.data$sliced_clusters,sort(as.numeric(levels(seuratObj@meta.data$sliced_clusters)))) #reorder levels
seuratObj@meta.data$annotated_clusters <- factor(seuratObj@meta.data$annotated_clusters,sort(as.numeric(levels(seuratObj@meta.data$annotated_clusters)))) #reorder levels
Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters

U_annot<-DimPlot(seuratObj, reduction = "umap", label = T, group.by = "sliced_clusters", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_merge_subset2/Annotation/2_UMAP_sliced_clusters_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
RNAMarkers_RNAclus <- FindAllMarkers(seuratObj, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_RNAclus_",sampleName,"_sliced.rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerRNAclustersliced']]<-paste0(table(RNAMarkers_RNAclus$cluster)," RNA markers for RNA cluster sliced ",rownames(table(RNAMarkers_RNAclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset2/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

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
names(RNAmarkersList_RNAclus)<-paste0("RNAclustersliced",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results_merge_subset2/Marker_lists/RNAmarkersList_RNAclus_",sampleName,"_sliced.xlsx"))


########################################
##### Markers annotated clusters
########################################
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident
levels(seuratObj@meta.data$annotated_clusters) <- c("Fibroblasts (Arachnoid, Stromal)",rep("Fibroblasts (Pial)",4),
                                                    "VLMC1","Fibroblasts (Dura, Stalk, ABC)","VLMC2",
                                                    "Endothelial cells","Proliferating Fibroblasts",
                                                    "FB/CPE doublets","Other VLMCs","FB/MF doublets")

U_annot<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_merge_subset2/Annotation/2_UMAP_annotated1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

RNAMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_SCTclus$cluster)
saveRDS(RNAMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_SCTclus_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerSCTclusterannotated']]<-paste0(table(RNAMarkers_SCTclus$cluster)," RNA markers for annotated cluster ",rownames(table(RNAMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset2/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrRNAclusters_SCTclus<-names(table(RNAMarkers_SCTclus$cluster))
# totalNrRNAclusters_SCTclusPlusOne<-totalNrRNAclusters_SCTclus
RNAmarkersList_SCTclus<-list()

for(i in totalNrRNAclusters_SCTclus){
  # clusterNr<-i-1
  
  tmp<-RNAMarkers_SCTclus[RNAMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  RNAmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
# names(RNAmarkersList_SCTclus)<-paste0("SCTclustersliced",0:totalNrRNAclusters_SCTclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_subset2/Marker_lists/RNAmarkersList_SCTclus_",sampleName,"_annotated.xlsx"))

######################################################

### Make heatmap for annotated clusters
RNAMarkers_SCTclus<- readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_SCTclus_",sampleName,"_annotated.rds"))

## Perform on a subset -> better view of smaller clusters!!
seuratObj.small <- subset(seuratObj, downsample = 500)

########## Get HVG ##########
seuratObj.small <- FindVariableFeatures(object = seuratObj.small, selection.method = "vst", nfeatures = nrow(seuratObj.small@assays$RNA))
length(VariableFeatures(seuratObj.small))

########## Scale ##########
seuratObj.small <- ScaleData(seuratObj.small)

## Heatmap RNA markers on RNA clusters
top10 <- RNAMarkers_SCTclus %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
D1<-DoHeatmap(seuratObj.small, features = top10$gene, group.by = "annotated_clusters") + NoLegend()
ggsave(D1, file=paste0(sampleFolder, "results_merge_subset2/Heatmaps/Heatmap_RNAmarkersList_Annotatedclus_",sampleName,".png"), 
       height = 20, width = 12, dpi = "retina")

pdf(file=paste0(sampleFolder, "results_merge_subset2/Heatmaps/Heatmap_RNAmarkersList_Annotatedclus_",sampleName,".pdf"), 
    height = 25, width = 25)
DoHeatmap(seuratObj.small, features = top10$gene, group.by = "annotated_clusters") + NoLegend()
dev.off()


######################################################

# Frequency tables (sliced)
Sample <- seuratObj@meta.data$sample_origin
cluster <- seuratObj@meta.data$annotated_clusters
Aggr <- rep(sampleName,length(cluster))

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

# Stacked
library("ggthemes")

png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/3_SampleDistribution_ggplot2_annot_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/3_SampleDistribution_ggplot2_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/3_SampleDistribution_ggplot2_annot_3_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data2, aes(fill=cluster, y=Freq, x=Aggr)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

# ##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
# data_split1<-data[which(data$Sample=="ABR1"),]
# data_split2<-data[which(data$Sample=="ABR2"),]
# data_split3<-data[which(data$Sample=="ABR3"),]
# data_split4<-data[which(data$Sample=="DV1"),]
# data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
# data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})
# data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-round((x/sum(data_split3$Freq))*100,2)})
# data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-round((x/sum(data_split4$Freq))*100,2)})
# 
# data_new<-rbind(data_split1,data_split2,data_split3,data_split4)
# 
# png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/3_SampleDistribution_ggplot2_annot_1_",sampleName,"_adjusted.png"), width = 2000, height = 1500, res = 300)
# ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + scale_fill_manual(values=c('Blue','Red',"Turquoise","Magenta")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
#   geom_bar(position="fill", stat="identity", colour="white")
# dev.off()

######################################################

## Extra vis
# seuratObj@meta.data$sample_origin<-as.factor(seuratObj@meta.data$orig.ident)
# levels(seuratObj@meta.data$sample_origin)<-c("FB_DeSisto_et_al_1","FB_DeSisto_et_al_2","FB_DeSisto_et_al_3",
#                                              "Ependymal_Zeisel_et_al","Vascular_Vanlandewijck_et_al",
#                                              "Ependymal_Shah_et_al_1","Ependymal_Shah_et_al_2","Ependymal_Shah_et_al_3",
#                                              "CP_LpsNeg_4V_Urvb","CP_LpsNeg_LV_Urvb","CP_Young_4V_Urvb","CP_Young_LV_Urvb",
#                                              "Vascular_Zeisel_et_al")

library(RColorBrewer)
Colorset<-c(brewer.pal(12,"Set3"),"magenta")

U_origin<-DimPlot(seuratObj, reduction = "umap", label = F, group.by = "sample_origin", label.size = 3, cols = Colorset)
ggsave(U_origin, file=paste0(sampleFolder,"results_merge_subset2/Annotation/4_UMAP_dataset_origin1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

U_origin2<-DimPlot(seuratObj, reduction = "umap", label = F, split.by ="sample_origin", group.by = "annotated_clusters", label.size = 3, ncol = 4)
ggsave(U_origin2, file=paste0(sampleFolder,"results_merge_subset2/Annotation/4_UMAP_dataset_origin2_",sampleName,".png"), height = 20, width = 25, dpi = 300)

U_origin3<-DimPlot(seuratObj, reduction = "umap", label = F, split.by ="sample_origin", group.by = "sample_origin", label.size = 3, cols = Colorset, ncol = 4)
ggsave(U_origin3, file=paste0(sampleFolder,"results_merge_subset2/Annotation/4_UMAP_dataset_origin3_",sampleName,".png"), height = 20, width = 25, dpi = 300)

pdf(file = paste0(sampleFolder,"results_merge_subset2/Annotation/4_UMAP_dataset_origin2_",sampleName,".pdf"), width = 15, height = 10)
U_origin2
dev.off()

pdf(file = paste0(sampleFolder,"results_merge_subset2/Annotation/4_UMAP_dataset_origin3_",sampleName,".pdf"), width = 15, height = 10)
U_origin3
dev.off()

########################################################################################################################
########################################################################################################################
########################################################################################################################

# Remove bad clusters!!
seuratObj_clean<-subset(seuratObj, idents =c("Fibroblasts (Arachnoid, Stromal)","Fibroblasts (Pial)",
                                             "VLMC1","Fibroblasts (Dura, Stalk, ABC)","VLMC2","Other VLMCs"))

# Clean outliers!!!!
Idents(seuratObj_clean)<-as.character(Idents(seuratObj_clean))

U1 <- DimPlot(seuratObj_clean, reduction = "umap", label = T, label.size = 4)
seuratObj_clean <- CellSelector(U1, object=seuratObj_clean, ident="Outliers1")

U_outlier<-DimPlot(seuratObj_clean, reduction = "umap", label = T, repel = T, label.size = 4)
ggsave(U_outlier, file=paste0(sampleFolder,"results_merge_subset2/Annotation/5.0_UMAP_clean_outliers_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

seuratObj_clean<-subset(seuratObj_clean, idents =c("Fibroblasts (Arachnoid, Stromal)","Fibroblasts (Pial)",
                                                  "VLMC1","Fibroblasts (Dura, Stalk, ABC)","VLMC2","Other VLMCs"))

# Remake figures
Idents(seuratObj_clean)<-as.factor(as.character(Idents(seuratObj_clean)))
seuratObj_clean@meta.data$annotated_clusters<-Idents(seuratObj_clean)

U_annot<-DimPlot(seuratObj_clean, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_merge_subset2/Annotation/5_UMAP_clean_annotated1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

pdf(file=paste0(sampleFolder, "results_merge_subset2/Annotation/5_UMAP_clean_annotated1_",sampleName,".pdf"), height = 10, width = 15)
U_annot
dev.off()

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

RNAMarkers_SCTclus <- FindAllMarkers(seuratObj_clean, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_SCTclus$cluster)
saveRDS(RNAMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_SCTclus_clean_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerSCTclustercleanannotated']]<-paste0(table(RNAMarkers_SCTclus$cluster)," RNA markers for clean annotated cluster ",rownames(table(RNAMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset2/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrRNAclusters_SCTclus<-names(table(RNAMarkers_SCTclus$cluster))
# totalNrRNAclusters_SCTclusPlusOne<-totalNrRNAclusters_SCTclus
RNAmarkersList_SCTclus<-list()

for(i in totalNrRNAclusters_SCTclus){
  # clusterNr<-i-1
  
  tmp<-RNAMarkers_SCTclus[RNAMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  RNAmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
# names(RNAmarkersList_SCTclus)<-paste0("SCTclustersliced",0:totalNrRNAclusters_SCTclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_subset2/Marker_lists/RNAmarkersList_SCTclus_clean_",sampleName,"_annotated.xlsx"))

### Make heatmap for annotated clusters
RNAMarkers_SCTclus<- readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_SCTclus_clean_",sampleName,"_annotated.rds"))

## Perform on a subset -> better view of smaller clusters!!
seuratObj_clean.small <- subset(seuratObj_clean, downsample = 500)

########## Get HVG ##########
seuratObj_clean.small <- FindVariableFeatures(object = seuratObj_clean.small, selection.method = "vst", nfeatures = nrow(seuratObj_clean.small@assays$RNA))
length(VariableFeatures(seuratObj_clean.small))

########## Scale ##########
seuratObj_clean.small <- ScaleData(seuratObj_clean.small)

## Heatmap RNA markers on RNA clusters
top10 <- RNAMarkers_SCTclus %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
D1<-DoHeatmap(seuratObj_clean.small, features = top10$gene, group.by = "annotated_clusters") + NoLegend()
ggsave(D1, file=paste0(sampleFolder, "results_merge_subset2/Heatmaps/Heatmap_RNAmarkersList_Annotatedclus_clean_",sampleName,".png"), 
       height = 20, width = 12, dpi = "retina")

pdf(file=paste0(sampleFolder, "results_merge_subset2/Heatmaps/Heatmap_RNAmarkersList_Annotatedclus_clean_",sampleName,".pdf"), 
    height = 25, width = 25)
DoHeatmap(seuratObj_clean.small, features = top10$gene, group.by = "annotated_clusters") + NoLegend()
dev.off()


######################################################

# Frequency tables (sliced)
Sample <- seuratObj_clean@meta.data$sample_origin
cluster <- seuratObj_clean@meta.data$annotated_clusters
Aggr <- rep(sampleName,length(cluster))

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

# Stacked
library("ggthemes")

png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/5_SampleDistribution_ggplot2_annot_1_",sampleName,"_clean.png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/5_SampleDistribution_ggplot2_annot_2_",sampleName,"_clean.png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/5_SampleDistribution_ggplot2_annot_3_",sampleName,"_clean.png"), width = 2000, height = 1500, res = 300)
ggplot(data2, aes(fill=cluster, y=Freq, x=Aggr)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

######################################################

## Extra vis
library(RColorBrewer)
Colorset<-c(brewer.pal(12,"Set3"),"magenta")

U_origin<-DimPlot(seuratObj_clean, reduction = "umap", label = F, group.by = "sample_origin", label.size = 3, cols = Colorset)
ggsave(U_origin, file=paste0(sampleFolder,"results_merge_subset2/Annotation/5_UMAP_clean_dataset_origin1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

U_origin2<-DimPlot(seuratObj_clean, reduction = "umap", label = F, split.by ="sample_origin", group.by = "annotated_clusters", label.size = 3, ncol = 4)
ggsave(U_origin2, file=paste0(sampleFolder,"results_merge_subset2/Annotation/5_UMAP_clean_dataset_origin2_",sampleName,".png"), height = 20, width = 25, dpi = 300)

U_origin3<-DimPlot(seuratObj_clean, reduction = "umap", label = F, split.by ="sample_origin", group.by = "sample_origin", label.size = 3, cols = Colorset, ncol = 4)
ggsave(U_origin3, file=paste0(sampleFolder,"results_merge_subset2/Annotation/5_UMAP_clean_dataset_origin3_",sampleName,".png"), height = 20, width = 25, dpi = 300)

pdf(file = paste0(sampleFolder,"results_merge_subset2/Annotation/5_UMAP_dataset_clean_origin2_",sampleName,".pdf"), width = 15, height = 10)
U_origin2
dev.off()

pdf(file = paste0(sampleFolder,"results_merge_subset2/Annotation/5_UMAP_dataset_clean_origin3_",sampleName,".pdf"), width = 15, height = 10)
U_origin3
dev.off()

#############################
### Extra detail clusters
############################
detail1_Pial_vs_Arach<-FindMarkers(seuratObj_clean, ident.1 = "Fibroblasts (Pial)", ident.2 = "Fibroblasts (Arachnoid, Stromal)",
                                   min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

detail2_Pial_vs_Dura<-FindMarkers(seuratObj_clean, ident.1 = "Fibroblasts (Pial)", ident.2 = "Fibroblasts (Dura, Stalk, ABC)",
                                   min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

detail3_Arach_vs_Dura<-FindMarkers(seuratObj_clean, ident.1 = "Fibroblasts (Arachnoid, Stromal)", ident.2 = "Fibroblasts (Dura, Stalk, ABC)",
                                  min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

detail4_Arach_vs_VLMC2<-FindMarkers(seuratObj_clean, ident.1 = "Fibroblasts (Arachnoid, Stromal)", ident.2 = "VLMC2",
                                   min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

detail5_VLMC2_vs_Dura<-FindMarkers(seuratObj_clean, ident.1 = "VLMC2", ident.2 = "Fibroblasts (Dura, Stalk, ABC)",
                                   min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

##### Create list
listDEgenesExtra<-tibble::lst(detail1_Pial_vs_Arach, detail2_Pial_vs_Dura, detail3_Arach_vs_Dura,
                              detail4_Arach_vs_VLMC2,detail5_VLMC2_vs_Dura)

##Add geneSymbol in column (for the export)
listDEgenesExtra<-lapply(listDEgenesExtra,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenesExtra<-lapply(listDEgenesExtra, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDEgenesExtra<-lapply(listDEgenesExtra, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                             x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
# listDEgenesExtra<-lapply(listDEgenesExtra, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
listDEgenesExtra<-lapply(listDEgenesExtra,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDEgenesExtra,file=paste0(sampleFolder,"results_merge_subset2/Robjects/detailClusters_clean_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDEgenesExtra, paste0(sampleFolder,"results_merge_subset2/Marker_lists/detailClusters_clean_",sampleName,".xlsx"))


#########################################################################

## Correlated genes
H1<-corGenesHeatmap(seuratObj_clean.small,"Igfbp6")
H2<-corGenesHeatmap(seuratObj_clean.small,"Dpep1")
H3<-corGenesHeatmap(seuratObj_clean.small,"Fxyd5")
H4<-corGenesHeatmap(seuratObj_clean.small,"Gsn")
H5<-corGenesHeatmap(seuratObj_clean.small,"Slc38a2")
H6<-corGenesHeatmap(seuratObj_clean.small,"Col4a1")
H7<-corGenesHeatmap(seuratObj_clean.small,"Dcn")
H8<-corGenesHeatmap(seuratObj_clean.small,"Pdgfra")
H9<-corGenesHeatmap(seuratObj_clean.small,"Col1a1")

#########################################################################

## Change Stalk to Base for paper (30/01/23)
levels(seuratObj_clean$annotated_clusters)[2]<-"Fibroblasts (Dura, Base, ABC)"
levels(seuratObj_clean$New_clusters_clean)[4]<-"CP Base"
seuratObj_clean$FB_Base_signature_score1<-seuratObj_clean$FB_stalk_signature_score1
seuratObj_clean$FB_stalk_signature_score1<-NULL
Idents(seuratObj_clean)<-seuratObj_clean$annotated_clusters

#########################################################################

##### Read object
seuratObj_clean <- readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_clean_",sampleName,"_harmony_RNA.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,"_harmony_RNA.rds"))

##### Save object
saveRDS(seuratObj_clean, file=paste0(sampleFolder,"Robjects/seuratObj_clean_",sampleName,"_harmony_RNA.rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,"_harmony_RNA.rds"))

#########################################################################

########################################
##### different markers
########################################

# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

######### 1. CREATE NEW CLUSTERS ##########
DimPlot(seuratObj_clean, reduction = "umap", label = T, label.size = 3, group.by = "annotated_clusters")

### Create new clusters: split on source
seuratObjNew<-seuratObj_clean
Idents(seuratObjNew)<-seuratObjNew@meta.data$annotated_clusters
seuratObjNew@meta.data$newClustersTmp<-Idents(seuratObjNew)
seuratObjNew@meta.data$Treatment<-as.character(seuratObjNew@meta.data$New_clusters)
seuratObjNew@meta.data$Combo_annot<-paste0(seuratObjNew@meta.data$newClustersTmp,"_",seuratObjNew@meta.data$Treatment)
head(seuratObjNew@meta.data)

### Use the new clusters
Idents(seuratObjNew)<-seuratObjNew@meta.data$Combo_annot

DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, label.size = 3) + theme(legend.position = "none")

Table_freq_Combo<-as.data.frame(sort(table(seuratObjNew@meta.data$Combo_annot), decreasing = T))
write.xlsx(Table_freq_Combo, paste0(sampleFolder,"results_merge_subset2/Annotation/6_Frequency_table_combo_annotation_",sampleName,".xlsx"))


########## 2. GET MARKERS ##########
# getDEgenes<-function(ident1, ident2){
#   markersDiff <- FindMarkers(seuratObjNew, ident.1 = ident1, ident.2 = ident2,
#                              min.pct = 0.10) # No min diff pct!! , min.diff.pct = 0.15 or logFC logfc.threshold = 0.30,
#   markersDiff<-markersDiff[markersDiff$p_val_adj < 0.01,]
#   markersDiff<-markersDiff[order(markersDiff$avg_logFC, decreasing = T),]
# 
#   markersDiff$geneSymbol<-rownames(markersDiff)
#   markersDiff$pct.1<-markersDiff$pct.1+0.001
#   markersDiff$pct.2<-markersDiff$pct.2+0.001
# 
#   markersDiff<-rbind(markersDiff[markersDiff$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/pct.2*avg_logFC),
#                      markersDiff[markersDiff$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/pct.1*avg_logFC))
#   markersDiff<-markersDiff[order(markersDiff$score, decreasing = T),]
#   return(markersDiff)
# }

########## 2. GET MARKERS (everything!! IPA!!) ##########
getDEgenes<-function(ident1, ident2){
  markersDiff <- FindMarkers(seuratObjNew, ident.1 = ident1, ident.2 = ident2,
                             logfc.threshold = 0, min.pct = 0) #0.30, 0.10, min.diff.pct = 0.15
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

## Choose!!! SCT/RNA
DefaultAssay(seuratObjNew)<-"RNA"

# Clusters to compare!!!!
"Fibroblasts (Arachnoid, Stromal)_meninges 3 (arachnoid)"
"Fibroblasts (Arachnoid, Stromal)_Fibroblasts Type 1"
"Fibroblasts (Arachnoid, Stromal)_VLMC2"
"Fibroblasts (Dura, Stalk, ABC)_meninges 4 (dura)"
"Fibroblasts (Dura, Stalk, ABC)_ABC"
"Fibroblasts (Dura, Stalk, ABC)_Fibroblasts Type 2"

#### Get diff markers Arach cluster #####
Arach_M3_vs_Stromal<-getDEgenes("Fibroblasts (Arachnoid, Stromal)_meninges 3 (arachnoid)","Fibroblasts (Arachnoid, Stromal)_Fibroblasts Type 1")
Arach_M3_vs_Stromal<-Arach_M3_vs_Stromal[order(Arach_M3_vs_Stromal$avg_logFC,decreasing = T),]
head(Arach_M3_vs_Stromal)
dim(Arach_M3_vs_Stromal)

Arach_M3_vs_VLMC2<-getDEgenes("Fibroblasts (Arachnoid, Stromal)_meninges 3 (arachnoid)","Fibroblasts (Arachnoid, Stromal)_VLMC2")
Arach_M3_vs_VLMC2<-Arach_M3_vs_VLMC2[order(Arach_M3_vs_VLMC2$avg_logFC,decreasing = T),]
head(Arach_M3_vs_VLMC2)
dim(Arach_M3_vs_VLMC2)

Arach_Stromal_vs_VLMC2<-getDEgenes("Fibroblasts (Arachnoid, Stromal)_Fibroblasts Type 1","Fibroblasts (Arachnoid, Stromal)_VLMC2")
Arach_Stromal_vs_VLMC2<-Arach_Stromal_vs_VLMC2[order(Arach_Stromal_vs_VLMC2$avg_logFC,decreasing = T),]
head(Arach_Stromal_vs_VLMC2)
dim(Arach_Stromal_vs_VLMC2)

#### Get diff markers Dura cluster #####
Dura_M4_vs_Stalk<-getDEgenes("Fibroblasts (Dura, Stalk, ABC)_meninges 4 (dura)","Fibroblasts (Dura, Stalk, ABC)_Fibroblasts Type 2")
Dura_M4_vs_Stalk<-Dura_M4_vs_Stalk[order(Dura_M4_vs_Stalk$avg_logFC,decreasing = T),]
head(Dura_M4_vs_Stalk)
dim(Dura_M4_vs_Stalk)

Dura_M4_vs_ABC<-getDEgenes("Fibroblasts (Dura, Stalk, ABC)_meninges 4 (dura)","Fibroblasts (Dura, Stalk, ABC)_ABC")
Dura_M4_vs_ABC<-Dura_M4_vs_ABC[order(Dura_M4_vs_ABC$avg_logFC,decreasing = T),]
head(Dura_M4_vs_ABC)
dim(Dura_M4_vs_ABC)

Dura_Stalk_vs_ABC<-getDEgenes("Fibroblasts (Dura, Stalk, ABC)_Fibroblasts Type 2","Fibroblasts (Dura, Stalk, ABC)_ABC")
Dura_Stalk_vs_ABC<-Dura_Stalk_vs_ABC[order(Dura_Stalk_vs_ABC$avg_logFC,decreasing = T),]
head(Dura_Stalk_vs_ABC)
dim(Dura_Stalk_vs_ABC)

##add to list
listDiffMarkers_IPA<-tibble::lst(Arach_M3_vs_Stromal,Arach_M3_vs_VLMC2,Arach_Stromal_vs_VLMC2,
                                 Dura_M4_vs_Stalk,Dura_M4_vs_ABC,Dura_Stalk_vs_ABC)

lapply(listDiffMarkers_IPA, dim)
listDiffMarkers_IPA<-lapply(listDiffMarkers_IPA,function(x){x<-x[order(x$score, decreasing=T),]})

#Check settings
saveRDS(listDiffMarkers_IPA, file=paste0(sampleFolder,"results_merge_subset2/Robjects/markersDiffSamples_clean_IPA_",sampleName,".rds"))

### Write to Excel
library('openxlsx')
write.xlsx(listDiffMarkers_IPA, file = paste0(sampleFolder,"results_merge_subset2/Marker_lists/summaryDiffMarkers_clean_IPA_",sampleName,".xlsx"))

######################################################################################################

# Clusters to compare!!!!
"Fibroblasts (Arachnoid, Stromal)_Fibroblasts Type 1"
"Fibroblasts (Dura, Stalk, ABC)_Fibroblasts Type 2"


detail1_Fib1_in_Arach_vs_all<-FindMarkers(seuratObjNew, ident.1 = "Fibroblasts (Arachnoid, Stromal)_Fibroblasts Type 1",
                                   min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = T)

detail2_Fib2_in_Dura_vs_all<-FindMarkers(seuratObjNew, ident.1 = "Fibroblasts (Dura, Stalk, ABC)_Fibroblasts Type 2",
                                  min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = T)

Idents(seuratObjNew)<-seuratObjNew@meta.data$New_clusters #Revert

detail1_Fib1_vs_all<-FindMarkers(seuratObjNew, ident.1 = "Fibroblasts Type 1",
                                 min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = T)

detail2_Fib2_vs_all<-FindMarkers(seuratObjNew, ident.1 = "Fibroblasts Type 2",
                                 min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = T)

## Cross ref with markers CPE in FB Merge
library(openxlsx)
Markers_Full_Merge<-read.xlsx("Merge/results_merge/Marker_lists/RNAmarkersList_SCTclus_Merge_FB_datasets_annotated.xlsx", sheet = "Epithelial_cells")

Filter1<-intersect(Markers_Full_Merge$gene,rownames(detail1_Fib1_in_Arach_vs_all))
Filter2<-intersect(Markers_Full_Merge$gene,rownames(detail2_Fib2_in_Dura_vs_all))
Filter3<-intersect(Markers_Full_Merge$gene,rownames(detail1_Fib1_vs_all))
Filter4<-intersect(Markers_Full_Merge$gene,rownames(detail2_Fib2_vs_all))

## Filter lists for CPE markers
Filtered_Fib1_in_Arach_vs_all<-detail1_Fib1_in_Arach_vs_all[setdiff(rownames(detail1_Fib1_in_Arach_vs_all),Markers_Full_Merge$gene),]
Filtered_Fib2_in_Dura_vs_all<-detail2_Fib2_in_Dura_vs_all[setdiff(rownames(detail2_Fib2_in_Dura_vs_all),Markers_Full_Merge$gene),]
Filtered_Fib1_vs_all<-detail1_Fib1_vs_all[setdiff(rownames(detail1_Fib1_vs_all),Markers_Full_Merge$gene),]
Filtered_Fib2_vs_all<-detail2_Fib2_vs_all[setdiff(rownames(detail2_Fib2_vs_all),Markers_Full_Merge$gene),]

##### Create list
listDEgenesFBs<-tibble::lst(detail1_Fib1_in_Arach_vs_all, Filtered_Fib1_in_Arach_vs_all, 
                            detail1_Fib1_vs_all,Filtered_Fib1_vs_all,
                            detail2_Fib2_in_Dura_vs_all,Filtered_Fib2_in_Dura_vs_all,
                            detail2_Fib2_vs_all,Filtered_Fib2_vs_all)

##Add geneSymbol in column (for the export)
listDEgenesFBs<-lapply(listDEgenesFBs,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenesFBs<-lapply(listDEgenesFBs, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDEgenesFBs<-lapply(listDEgenesFBs, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                             x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
# listDEgenesFBs<-lapply(listDEgenesFBs, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
listDEgenesFBs<-lapply(listDEgenesFBs,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDEgenesFBs,file=paste0(sampleFolder,"results_merge_subset2/Robjects/Markers_for_our_FBs_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDEgenesFBs, paste0(sampleFolder,"results_merge_subset2/Marker_lists/Markers_for_our_FBs_",sampleName,".xlsx"))


######################################################################################################
## Average expression heatmaps
Idents(seuratObjNew)
# cluster.averages <- AverageExpression(seuratObj)
# head(cluster.averages[["SCT"]][, 1:5])

cluster.averages <- AverageExpression(seuratObjNew, return.seurat = TRUE)
cluster.averages

## Markers to show
FB1_filtered_markers<-as.character(head(rownames(Filtered_Fib1_vs_all),10))
FB2_filtered_markers<-as.character(head(rownames(Filtered_Fib2_vs_all),10))

FB_markers_heatmap<-c(FB1_filtered_markers,FB2_filtered_markers)

## Create heatmap
library("RColorBrewer")
Colorset<-c(brewer.pal(12,"Set3"),"#d3d3d3","magenta4","#00ffff","turquoise4",
            "Orangered", "Blue4")

H1 <- DoHeatmap(cluster.averages, assay = "RNA",
                features = FB_markers_heatmap, size = 3,
                draw.lines = FALSE, group.colors = Colorset) + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")

pdf(file=paste0(sampleFolder,"results_merge_subset2/Heatmaps/Test_heatmap_average_population_expression_",sampleName,".pdf"), width = 15, height = 15)
H1
dev.off()

######################################################################################################

Idents(seuratObjNew)<-seuratObjNew@meta.data$New_clusters #Revert

IPA_Fib1_vs_all<-FindMarkers(seuratObjNew, ident.1 = "Fibroblasts Type 1",
                                 min.pct = 0.01, logfc.threshold = 0.01)

IPA_Fib2_vs_all<-FindMarkers(seuratObjNew, ident.1 = "Fibroblasts Type 2",
                                 min.pct = 0.01, logfc.threshold = 0.01)

## Filter lists for CPE markers
Filtered_IPA_Fib1_vs_all<-IPA_Fib1_vs_all[setdiff(rownames(IPA_Fib1_vs_all),Markers_Full_Merge$gene),]
Filtered_IPA_Fib2_vs_all<-IPA_Fib2_vs_all[setdiff(rownames(IPA_Fib2_vs_all),Markers_Full_Merge$gene),]

##### Create list
list_IPA_FBs<-tibble::lst(IPA_Fib1_vs_all,Filtered_IPA_Fib1_vs_all,
                            IPA_Fib2_vs_all,Filtered_IPA_Fib2_vs_all)

##Add geneSymbol in column (for the export)
list_IPA_FBs<-lapply(list_IPA_FBs,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
# list_IPA_FBs<-lapply(list_IPA_FBs, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
list_IPA_FBs<-lapply(list_IPA_FBs, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                         x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
# list_IPA_FBs<-lapply(list_IPA_FBs, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
list_IPA_FBs<-lapply(list_IPA_FBs,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(list_IPA_FBs,file=paste0(sampleFolder,"results_merge_subset2/Robjects/Markers_IPA_for_our_FBs_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(list_IPA_FBs, paste0(sampleFolder,"results_merge_subset2/Marker_lists/Markers_IPA_for_our_FBs_",sampleName,".xlsx"))

######################################################################################################

## New annotation: combine dataset names
seuratObj_clean@meta.data$sample_origin2<-seuratObj_clean@meta.data$sample_origin
levels(seuratObj_clean@meta.data$sample_origin2)<-c("FB_DeSisto_et_al","FB_DeSisto_et_al",         
"FB_DeSisto_et_al","Ependymal_Zeisel_et_al" ,"Vascular_Vanlandewijck_et_al","Ependymal_Shah_et_al",      
"Ependymal_Shah_et_al","Ependymal_Shah_et_al","CP_Urvb","CP_Urvb","CP_Urvb","CP_Urvb","Vascular_Zeisel_et_al")

# Frequency tables (sliced)
Sample <- seuratObj_clean@meta.data$sample_origin2
cluster <- seuratObj_clean@meta.data$annotated_clusters
Aggr <- rep(sampleName,length(cluster))

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

# Stacked
library("ggthemes")

png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/7_SampleDistribution_ggplot2_annot_1_",sampleName,"_clean.png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/7_SampleDistribution_ggplot2_annot_2_",sampleName,"_clean.png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png(file=paste0(sampleFolder,"results_merge_subset2/Annotation/7_SampleDistribution_ggplot2_annot_3_",sampleName,"_clean.png"), width = 2000, height = 1500, res = 300)
ggplot(data2, aes(fill=cluster, y=Freq, x=Aggr)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

######################################################

## Extra vis
library(RColorBrewer)
Colorset<-c(brewer.pal(8,"Set1")[c(1,8,6,2,3,4)])

U_origin<-DimPlot(seuratObj_clean, reduction = "umap", label = F, group.by = "sample_origin2", label.size = 3, cols = Colorset)
ggsave(U_origin, file=paste0(sampleFolder,"results_merge_subset2/Annotation/7_UMAP_clean_dataset_origin1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

U_origin2<-DimPlot(seuratObj_clean, reduction = "umap", label = F, split.by ="sample_origin2", group.by = "annotated_clusters", label.size = 3, ncol = 3)
ggsave(U_origin2, file=paste0(sampleFolder,"results_merge_subset2/Annotation/7_UMAP_clean_dataset_origin2_",sampleName,".png"), height = 20, width = 25, dpi = 300)

U_origin3<-DimPlot(seuratObj_clean, reduction = "umap", label = F, split.by ="sample_origin2", group.by = "sample_origin2", label.size = 3, cols = Colorset, ncol = 3)
ggsave(U_origin3, file=paste0(sampleFolder,"results_merge_subset2/Annotation/7_UMAP_clean_dataset_origin3_",sampleName,".png"), height = 20, width = 25, dpi = 300)

pdf(file = paste0(sampleFolder,"results_merge_subset2/Annotation/7_UMAP_dataset_clean_origin2_",sampleName,".pdf"), width = 15, height = 10)
U_origin2
dev.off()

pdf(file = paste0(sampleFolder,"results_merge_subset2/Annotation/7_UMAP_dataset_clean_origin3_",sampleName,".pdf"), width = 15, height = 10)
U_origin3
dev.off()

## Visualize detailed metadata Desisto paper
seuratObj_clean@meta.data$annotated_clusters_detail<-as.character(seuratObj_clean@meta.data$annotated_clusters)
seuratObj_clean@meta.data[!is.na(seuratObj_clean$SubCluster),"annotated_clusters_detail"]<-seuratObj_clean@meta.data[!is.na(seuratObj_clean$SubCluster),"SubCluster"]
U1<-DimPlot(seuratObj_clean, reduction = "umap", label = T, label.size = 5, pt.size = 1, repel = T, group.by = "annotated_clusters_detail")
U2<-DimPlot(seuratObj_clean, reduction = "umap", label = T, label.size = 5, repel = T, pt.size = 1, group.by = "SubCluster")

## Load in data to check
clusterMatrix<-seuratObj_clean@meta.data
umapTable<-as.data.frame(seuratObj_clean[['umap']]@cell.embeddings, stringsAsFactors = F)

## Check location three clusters
U3<-colorSomeCells(clusterMatrix, umapTable, 
               WhichCells(seuratObj_clean, cells = rownames(seuratObj_clean@meta.data[which(seuratObj_clean@meta.data$SubCluster=="M2-6"),])))

pdf(file = paste0(sampleFolder,"results_merge_subset2/Annotation/15_UMAP_Desisto_detailed_annotation_",sampleName,".pdf"), width = 15, height = 10)
U1
U2
U3
dev.off()

######################################################################################################

## Update old metadata: change to NewClusters annotation (09/03/21)
seuratObj_clean@meta.data$New_clusters_clean<-as.factor(seuratObj_clean@meta.data$New_clusters)
levels(seuratObj_clean@meta.data$New_clusters_clean)
Idents(seuratObj_clean)<- seuratObj_clean@meta.data$New_clusters_clean

## Load in data to check
clusterMatrix<-seuratObj_clean@meta.data
umapTable<-as.data.frame(seuratObj_clean[['umap']]@cell.embeddings, stringsAsFactors = F)

## Check location three clusters
colorSomeCells(clusterMatrix, umapTable, 
               WhichCells(seuratObj_clean, cells = rownames(seuratObj_clean@meta.data[which(seuratObj_clean@meta.data$New_clusters_clean=="Doublets_3" | 
                                                                                            seuratObj_clean@meta.data$New_clusters_clean=="Doublets_5" |
                                                                                            seuratObj_clean@meta.data$New_clusters_clean=="Unknown"),])))

Extra_cells_Urvb<-WhichCells(seuratObj_clean, cells = rownames(seuratObj_clean@meta.data[which(seuratObj_clean@meta.data$New_clusters_clean=="Doublets_3" | 
                                                                               seuratObj_clean@meta.data$New_clusters_clean=="Doublets_5" |
                                                                               seuratObj_clean@meta.data$New_clusters_clean=="Unknown"),]))
## Divide based on UMAP location: bottom right -> classify as FB Type 2
umapSlice1<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 > 3, UMAP_2 < -2) 
wantedCells1<-setdiff(Extra_cells_Urvb,umapSlice1$cell) #Type 1
wantedCells2<-intersect(umapSlice1$cell, Extra_cells_Urvb) #Type 2

colorSomeCells(clusterMatrix, umapTable, wantedCells1)
colorSomeCells(clusterMatrix, umapTable, wantedCells2)

## Reannotate
seuratObj_clean<-SetIdent(object = seuratObj_clean, cells = wantedCells1, value = "Fibroblasts Type 1")
seuratObj_clean<-SetIdent(object = seuratObj_clean, cells = wantedCells2, value = "Fibroblasts Type 2")
DimPlot(seuratObj_clean, reduction = "umap", label = T, label.size = 4)

seuratObj_clean@meta.data$New_clusters_clean<-Idents(seuratObj_clean)

## Sort to order Daan
seuratObj_clean@meta.data$New_clusters_clean<-factor(seuratObj_clean@meta.data$New_clusters_clean, 
                                                     levels=levels(seuratObj_clean@meta.data$New_clusters_clean)[c(2,4,5,1,3,7,6,14,15,8,9,10,11,12,13)])
Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean
DimPlot(seuratObj_clean, reduction = "umap", label = T, label.size = 4)

###################

## Determine DEG for this annotation
library(future)
plan("multiprocess", workers = 6)
plan()

RNAMarkers_SCTclus <- FindAllMarkers(seuratObj_clean, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_SCTclus$cluster)
saveRDS(RNAMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_Old_annotation_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['RNAmarkersOldAnnotation']]<-paste0(table(RNAMarkers_SCTclus$cluster)," RNA markers for old annotated cluster ",rownames(table(RNAMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset2/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrRNAclusters_SCTclus<-names(table(RNAMarkers_SCTclus$cluster))
# totalNrRNAclusters_SCTclusPlusOne<-totalNrRNAclusters_SCTclus
RNAmarkersList_SCTclus<-list()

for(i in totalNrRNAclusters_SCTclus){
  # clusterNr<-i-1
  
  tmp<-RNAMarkers_SCTclus[RNAMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  RNAmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
# names(RNAmarkersList_SCTclus)<-paste0("SCTclustersliced",0:totalNrRNAclusters_SCTclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_subset2/Marker_lists/RNAmarkersList_Old_annotation_",sampleName,".xlsx"))

###################

## Determine DEG for this annotation (stricter!!)
RNAMarkers_SCTclus_strict <- FindAllMarkers(seuratObj_clean, assay = "RNA", min.pct = 0.6, only.pos = TRUE)
table(RNAMarkers_SCTclus_strict$cluster)
saveRDS(RNAMarkers_SCTclus_strict, file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_Old_annotation_",sampleName,"_strict.rds"))

### Add to diagnostics
diagnostics[['RNAmarkersOldAnnotationStrict']]<-paste0(table(RNAMarkers_SCTclus_strict$cluster)," RNA markers for old annotated cluster strict ",rownames(table(RNAMarkers_SCTclus_strict$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset2/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrRNAclusters_SCTclus_strict<-names(table(RNAMarkers_SCTclus_strict$cluster))
# totalNrRNAclusters_SCTclusPlusOne<-totalNrRNAclusters_SCTclus_strict
RNAmarkersList_SCTclus_strict<-list()

for(i in totalNrRNAclusters_SCTclus_strict){
  # clusterNr<-i-1
  
  tmp<-RNAMarkers_SCTclus_strict[RNAMarkers_SCTclus_strict$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  RNAmarkersList_SCTclus_strict[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
# names(RNAmarkersList_SCTclus_strict)<-paste0("SCTclustersliced",0:totalNrRNAclusters_SCTclus_strict)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_SCTclus_strict, file =paste0(sampleFolder, "results_merge_subset2/Marker_lists/RNAmarkersList_Old_annotation_",sampleName,"_strict.xlsx"))

###################

## Cross ref with markers CPE in FB Merge
library(openxlsx)
Markers_Full_Merge<-read.xlsx("Merge/results_merge/Marker_lists/RNAmarkersList_SCTclus_Merge_FB_datasets_annotated.xlsx", sheet = "Epithelial_cells")

Filter1<-intersect(Markers_Full_Merge$gene,RNAMarkers_SCTclus_strict[which(RNAMarkers_SCTclus_strict$cluster == "Fibroblasts Type 1"),"gene"])
Filter2<-intersect(Markers_Full_Merge$gene,RNAMarkers_SCTclus_strict[which(RNAMarkers_SCTclus_strict$cluster == "Fibroblasts Type 2"),"gene"])

## Filter lists for CPE markers
RNAmarkersList_SCTclus_strict_filtered<-RNAmarkersList_SCTclus_strict
RNAmarkersList_SCTclus_strict_filtered$`Fibroblasts Type 1`<-RNAmarkersList_SCTclus_strict_filtered$`Fibroblasts Type 1`[setdiff(RNAmarkersList_SCTclus_strict_filtered$`Fibroblasts Type 1`$gene,Markers_Full_Merge$gene),]
RNAmarkersList_SCTclus_strict_filtered$`Fibroblasts Type 2`<-RNAmarkersList_SCTclus_strict_filtered$`Fibroblasts Type 2`[setdiff(RNAmarkersList_SCTclus_strict_filtered$`Fibroblasts Type 2`$gene,Markers_Full_Merge$gene),]

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_SCTclus_strict_filtered, file =paste0(sampleFolder, "results_merge_subset2/Marker_lists/RNAmarkersList_Old_annotation_",sampleName,"_strict_filtered.xlsx"))

###################
# top3_logFC <- RNAMarkers_SCTclus_strict %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
top3_score<- as.character()
for (pop in 1:length(RNAmarkersList_SCTclus_strict_filtered)) {
  extra_genes<-head(RNAmarkersList_SCTclus_strict_filtered[[pop]]$gene,3)
  top3_score<-c(top3_score,extra_genes)
}

## Create dotplot of top markers old annotation
Colors_dotplot<-c("#071AE5","#F50635") #030720

wantedGenes<-top3_score
wantedGenes<-rev(unique(wantedGenes))
D1<-DotPlot(seuratObj_clean, features = wantedGenes, cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/08_Dotplot_top3_markers_score_strict_filtered.pdf"), width = 15, height = 10)
DotPlot(seuratObj_clean, features = wantedGenes, cols = Colors_dotplot) + RotatedAxis()
dev.off()

# library("ggthemes")
# library("scales")
# library("ggpubr")
# D1<-DotPlot(seuratObj_clean, features = wantedGenes, cols = Colors_dotplot)
# ggpar(D1, orientation = "horizontal") #ggpubr package to reorient ggplot objects

###################

## Create average heatmap again with new order and 90 degree twist
## Average expression heatmaps
Idents(seuratObj_clean)

cluster.averages <- AverageExpression(seuratObj_clean, return.seurat = TRUE)
cluster.averages

## Markers to show (strict-filtered list!)
Filtered_strict_list_FB1s<-read.xlsx(xlsxFile = paste0(sampleFolder, "results_merge_subset2/Marker_lists/RNAmarkersList_Old_annotation_",sampleName,"_strict_filtered.xlsx"), sheet = 1)
Filtered_strict_list_FB2s<-read.xlsx(xlsxFile = paste0(sampleFolder, "results_merge_subset2/Marker_lists/RNAmarkersList_Old_annotation_",sampleName,"_strict_filtered.xlsx"), sheet = 4)

FB1_filtered_markers<-as.character(head(Filtered_strict_list_FB1s$gene,10))
FB2_filtered_markers<-as.character(head(Filtered_strict_list_FB2s$gene,10))

FB_markers_heatmap<-c(FB1_filtered_markers,FB2_filtered_markers)

## Create heatmap
library("RColorBrewer")
Colorset<-c(brewer.pal(12,"Set3"),"magenta4","#00ffff", "Blue4")

H1 <- DoHeatmap(cluster.averages, assay = "RNA",
                features = FB_markers_heatmap, size = 3,
                draw.lines = FALSE, group.colors = Colorset) + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")

# library("ggthemes")
# library("scales")
library("ggpubr")
H2 <- DoHeatmap(cluster.averages, assay = "RNA",
                features = rev(FB_markers_heatmap), size = 3,
                draw.lines = FALSE, group.colors = Colorset) + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")
H2<-ggpar(H2, orientation = "horizontal") #ggpubr package to reorient ggplot objects

pdf(file=paste0(sampleFolder,"results_merge_subset2/Heatmaps/Heatmap_average_population_expression_",sampleName,"_vertical.pdf"), width = 15, height = 15)
H1
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Heatmaps/Heatmap_average_population_expression_",sampleName,"_horizontal.pdf"), width = 15, height = 10)
H2
dev.off()

## Revert back Idents of object
Idents(seuratObj_clean)<-seuratObj_clean@meta.data$annotated_clusters
DimPlot(seuratObj_clean, reduction = "umap", label = T, label.size = 4)

######################################################################################################

################################################################################
########## EXPORT FOR TSNE TOOL
################################################################################
dir.create(paste0(sampleFolder,"results_merge_subset2/neededFilesOnlineTool"))

listLabels<-list(levels(seuratObj_clean@meta.data$sample_origin))

exportShiny<-DisneyTools::expShiny(seuratObj_clean, conditions = unlist(listLabels), assay = "RNA", clustering.labels = "sliced_clusters") #Use annotated version!
saveRDS(exportShiny, file=paste0(sampleFolder,"results_merge_subset2/neededFilesOnlineTool/SO_",sampleName,".rds"))
annotation <- exportShiny@meta.data$annotated_clusters %>% as.factor()
saveRDS(annotation, file=paste0(sampleFolder,"results_merge_subset2/neededFilesOnlineTool/SO_",sampleName,"_annotation.rds"))

######################################################################################################

## Update meeting 22/03/21 (previously made excel files not yet changed!!)
levels(seuratObj_clean@meta.data$New_clusters_clean)[c(1,4)]<-c("CP Stromal","CP Stalk")

U_old_annot_clean<-DimPlot(seuratObj_clean, reduction = "umap", label = T, repel = T, label.size = 4, group.by = "New_clusters_clean")
ggsave(U_old_annot_clean, file=paste0(sampleFolder,"results_merge_subset2/Annotation/9_UMAP_old_annot_clean_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/9_UMAP_old_annot_clean_",sampleName,".pdf"), width = 15, height = 10)
U_old_annot_clean
dev.off()

##Increase dot size for paper
U_old_annot_clean_v2<-DimPlot(seuratObj_clean, reduction = "umap", label = T, repel = T, label.size = 4,
                              pt.size = 1, group.by = "New_clusters_clean")

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/9_UMAP_old_annot_clean_bigger_dots_",sampleName,".pdf"), width = 15, height = 10)
U_old_annot_clean_v2
dev.off()

###################

### Read markers and perform filtering (workflow above)
RNAMarkers_SCTclus_strict<-readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/RNAmarkersList_Old_annotation_",sampleName,"_strict.rds"))
## Re-perform filtering (see above!!!!!!!!!!!!!!)

# top3_logFC <- RNAMarkers_SCTclus_strict %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
top3_score<- as.character()
for (pop in 1:length(RNAmarkersList_SCTclus_strict_filtered)) {
  extra_genes<-head(RNAmarkersList_SCTclus_strict_filtered[[pop]]$gene,3)
  top3_score<-c(top3_score,extra_genes)
}

## Create dotplot of top markers old annotation
Colors_dotplot<-c("#071AE5","#F50635") #030720

wantedGenes<-top3_score
wantedGenes<-rev(unique(wantedGenes))
DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = wantedGenes, cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/9_Dotplot_top3_markers_score_strict_filtered_new.pdf"), width = 15, height = 10)
DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = wantedGenes, cols = Colors_dotplot) + RotatedAxis()
dev.off()

###################

## Final plots paper May 2022

library(viridis)

## Feature plots panel D
features<-c("Dcn", "Dpep1", "Igfbp6" , "Cldn11")

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Paper/Feature_plot_paper_4_markers_check_",sampleName,"_viridisC_ordered.pdf"), height = 10, width = 15)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj_clean, features =feature, cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) +
  scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

##Modulescore
FB_stalk_signature <- list(c("Cldn11","Igfbp6"))

seuratObj_clean <- AddModuleScore(object = seuratObj_clean, assay = "RNA", features = FB_stalk_signature, name = "FB_stalk_signature_score")
F2 <- FeaturePlot(object = seuratObj_clean, features = "FB_stalk_signature_score1", 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)  + scale_color_viridis(option = "D")
pdf(file = paste0(sampleFolder,"results_merge_subset2/Feature_plots/Paper/Featureplot_modulescore_Cldn11_Igfbp6_",sampleName,".pdf"), width = 15, height = 10)
F2
dev.off()

## Extra plot paper October 2022
library(viridis)

## Feature plots panel D
features<-c("Tcf21")

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Paper/Feature_plot_paper_Tcf21_",sampleName,"_viridisC_ordered.pdf"), height = 10, width = 15)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj_clean, features =feature, cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) +
    scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

###################
## Check Lrat expr
F1<-FeaturePlot(object = seuratObj_clean, features ="Lrat", cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=F)

#No expr!

###################
## Check Nestin (Nes) expr
F1<-FeaturePlot(object = seuratObj_clean, features ="Nes", cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)


###################
## Check Pdgfra en Pdgfrb expr
F1<-FeaturePlot(object = seuratObj_clean, features =c("Pdgfra","Pdgfrb"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)


###################
## Check Adam12 expr
F1<-FeaturePlot(object = seuratObj_clean, features ="Adam12", cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)

###################
## Check ACE2  expr
F1<-FeaturePlot(object = seuratObj_clean, features ="Ace2", cols = c("yellow", "red"), label = T,
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)

pdf(file = paste0(sampleFolder,"results_merge_subset2/Feature_plots/Featureplot_Ace2_",sampleName,".pdf"), width = 15, height = 10)
F1
dev.off()

###################
## DE analysis with new clean metadata! (Including Unknown and doublets!!)
# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)

# Clusters to compare!!!!
"CP Stromal"
"CP Stalk"

Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean

Fib1_All_vs_Fib2_All<-FindMarkers(seuratObj_clean, ident.1 = "CP Stromal", ident.2 = "CP Stalk",
                                  min.pct = 0.10, logfc.threshold = 0.30, only.pos = F)

### Create new clusters: split on source
seuratObjNew<-seuratObj_clean
Idents(seuratObjNew)<-seuratObjNew@meta.data$annotated_clusters
seuratObjNew@meta.data$newClustersTmp<-Idents(seuratObjNew)
seuratObjNew@meta.data$Treatment<-as.character(seuratObjNew@meta.data$New_clusters_clean)
seuratObjNew@meta.data$Combo_annot<-paste0(seuratObjNew@meta.data$newClustersTmp,"_",seuratObjNew@meta.data$Treatment)
head(seuratObjNew@meta.data)

### Use the new clusters
Idents(seuratObjNew)<-seuratObjNew@meta.data$Combo_annot

# Clusters to compare!!!!
"Fibroblasts (Arachnoid, Stromal)_CP Stromal"
"Fibroblasts (Dura, Stalk, ABC)_CP Stalk"

DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, label.size = 3) + theme(legend.position = "none")

Fib1_Arach_vs_Fib2_Dura<-FindMarkers(seuratObjNew, ident.1 = "Fibroblasts (Arachnoid, Stromal)_CP Stromal", ident.2 = "Fibroblasts (Dura, Stalk, ABC)_CP Stalk",
                                     min.pct = 0.10, logfc.threshold = 0.30, only.pos = F)


##### Create list
listDEgenesFBs<-tibble::lst(Fib1_All_vs_Fib2_All,Fib1_Arach_vs_Fib2_Dura)

##Add geneSymbol in column (for the export)
listDEgenesFBs<-lapply(listDEgenesFBs,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenesFBs<-lapply(listDEgenesFBs, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDEgenesFBs<-lapply(listDEgenesFBs, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                         x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
# listDEgenesFBs<-lapply(listDEgenesFBs, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
listDEgenesFBs<-lapply(listDEgenesFBs,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDEgenesFBs,file=paste0(sampleFolder,"results_merge_subset2/Robjects/Markers_for_our_FBs_paper_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDEgenesFBs, paste0(sampleFolder,"results_merge_subset2/Marker_lists/Markers_for_our_FBs_paper_",sampleName,".xlsx"))

## Revert annotation
Idents(seuratObj_clean)<-seuratObj_clean@meta.data$annotated_clusters #Revert

###################

### GOEA ###
library(clusterProfiler)

## Stromal FBs
Test<-enrichGO(
  as.character(listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC > 0),"gene"]),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D1<-dotplot(Test, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Stalk FBs
Test2<-enrichGO(
  as.character(listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC < 0),"gene"]),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D2<-dotplot(Test2, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Combine
D_comp<-merge_result(list(Stromal_CP=Test, Stalk_CP=Test2)) %>%
  dotplot(., split="ONTOLOGY", showCategory=15) + facet_grid(ONTOLOGY~., scale="free")

# ## Combine
# Test_combo<-Test@result
# Test_combo2<-Test2@result
# 
# Test_combo$Type<-"Stromal FBs"
# Test_combo2$Type<-"Stalk FBs"
# 
# Test_combo3<-rbind(Test_combo,Test_combo2)
# Test3<-Test
# Test3@result<-Test_combo3
# Test3@result$Type<-as.character(Test3@result$Type)
# Test3@result$Description<-as.character(Test3@result$Description)
# 
# Test3@gene<-c(Test@gene,Test2@gene)
# Test3@geneSets<-c(Test@geneSets,Test2@geneSets)
# 
# dotplot(Test3, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
# 
# dotplot(Test3, split="Type", showCategory=30) + facet_grid(Type~., scale="free")
# 
# ## from Tommy's code
# p <- ggplot(Test3) + 
#   geom_point(aes(x = Type, y = reorder(Description,GeneRatio), size = GeneRatio, color = p.adjust)) +
#   theme_bw(base_size = 14) +
#   scale_colour_gradient(limits=c(0, 0.10), low="red") +
#   ylab(NULL) +
#   ggtitle("GO pathway enrichment")
# 
# p + facet_grid(Type~., scale="free")

## Stromal FBs Arach
Test_Arach<-enrichGO(
  as.character(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC > 0),"gene"]),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D3<-dotplot(Test_Arach, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Stalk FBs Dura
Test2_Dura<-enrichGO(
  as.character(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC < 0),"gene"]),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D4<-dotplot(Test2_Dura, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Combine
D_comp2<-merge_result(list(Stromal_CP=Test_Arach, Stalk_CP=Test2_Dura)) %>%
  dotplot(., split="ONTOLOGY", showCategory=15) + facet_grid(ONTOLOGY~., scale="free")


## Check genes associated with GO terms
GeneList_Test<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC > 0),"avg_logFC"]
names(GeneList_Test)<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC > 0),"gene"]
H1<-heatplot(Test, foldChange=GeneList_Test)

GeneList_Test2<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC < 0),"avg_logFC"]*(-1) #Convert logFC
names(GeneList_Test2)<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC < 0),"gene"]
H2<-heatplot(Test2, foldChange=GeneList_Test2)

GeneList_Test_Arach<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC > 0),"avg_logFC"]
names(GeneList_Test_Arach)<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC > 0),"gene"]
H3<-heatplot(Test_Arach, foldChange=GeneList_Test_Arach)

GeneList_Test2_Dura<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC < 0),"avg_logFC"]*(-1) #Convert logFC
names(GeneList_Test2_Dura)<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC < 0),"gene"]
H4<-heatplot(Test2_Dura, foldChange=GeneList_Test2_Dura)


#### Adapted comparison plot for FB1 All vs FB2 All (26/03/21)
## Filter
Test_filtered<-Test
Test2_filtered<-Test2
Test_filtered@result<-Test_filtered@result[-which(Test_filtered@result$ONTOLOGY == "MF"),]
Test2_filtered@result<-Test2_filtered@result[-which(Test2_filtered@result$ONTOLOGY == "MF"),]
## Combine
D_comp_filtered<-merge_result(list(Stromal_CP=Test_filtered, Stalk_CP=Test2_filtered)) %>%
  dotplot(., split="ONTOLOGY", showCategory=25) + facet_grid(ONTOLOGY~., scale="free")
## Version Roos (28/04/22)
D_comp_filtered_v2<-merge_result(list(Stalk_CP=Test2_filtered, Stromal_CP=Test_filtered)) %>%
  dotplot(., split="ONTOLOGY", showCategory=25) + facet_grid(ONTOLOGY~., scale="free")

## Version SING (06/05/22)
D_comp_filtered_SING<-dotplot(Test2_filtered, split="ONTOLOGY", showCategory=25) + facet_grid(ONTOLOGY~., scale="free")
D_comp_filtered_SING_v2<-merge_result(list(Stalk_CP=Test2_filtered)) %>%
  dotplot(., split="ONTOLOGY", showCategory=25) + facet_grid(ONTOLOGY~., scale="free")

# ###############
# ## March 2022 update
# Test_filtered_v2<-Test
# Test2_filtered_v2<-Test2
# Test_filtered_v2@result<-Test_filtered_v2@result[-which(Test_filtered_v2@result$ONTOLOGY == "MF" | Test_filtered_v2@result$ONTOLOGY == "CC"),]
# Test2_filtered_v2@result<-Test2_filtered_v2@result[-which(Test2_filtered_v2@result$ONTOLOGY == "MF" | Test2_filtered_v2@result$ONTOLOGY == "CC" ),]
# Test_filtered_v3<-Test
# Test2_filtered_v3<-Test2
# Test_filtered_v3@result<-Test_filtered_v3@result[-which(Test_filtered_v3@result$ONTOLOGY == "MF" | Test_filtered_v3@result$ONTOLOGY == "BP"),]
# Test2_filtered_v3@result<-Test2_filtered_v3@result[-which(Test2_filtered_v3@result$ONTOLOGY == "MF" | Test2_filtered_v3@result$ONTOLOGY == "BP" ),]
# ## Combine
# D_comp_filtered_v2<-merge_result(list(Stromal_CP=Test_filtered_v2, Stalk_CP=Test2_filtered_v2)) %>%
#   dotplot(., split="ONTOLOGY", showCategory=25) + facet_grid(ONTOLOGY~., scale="free")
# D_comp_filtered_v3<-merge_result(list(Stromal_CP=Test_filtered_v3, Stalk_CP=Test2_filtered_v3)) %>%
#   dotplot(., split="ONTOLOGY", showCategory=25) + facet_grid(ONTOLOGY~., scale="free")
# ###############

### Save results ###
pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB1s_All_vs_FB2s_All_",sampleName,".pdf"), width = 10, height = 10)
D1
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB2s_All_vs_FB1s_All_",sampleName,".pdf"), width = 10, height = 10)
D2
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB1s_Arach_vs_FB2s_Dura_",sampleName,".pdf"), width = 10, height = 10)
D3
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB2s_Dura_vs_FB1s_Arach_",sampleName,".pdf"), width = 10, height = 10)
D4
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_",sampleName,".pdf"), width = 10, height = 15)
D_comp
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_Arach_and_FB2s_Dura_",sampleName,".pdf"), width = 10, height = 15)
D_comp2
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_",sampleName,"_without_MF.pdf"), width = 10, height = 15)
D_comp_filtered
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_order_Roos_",sampleName,"_without_MF.pdf"), width = 10, height = 15)
D_comp_filtered_v2
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB2s_All_vs_FB1s_All_SING_",sampleName,"_without_MF.pdf"), width = 15, height = 12)
D_comp_filtered_SING
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB2s_All_vs_FB1s_All_SINGv2_",sampleName,"_without_MF.pdf"), width = 8, height = 11)
D_comp_filtered_SING_v2
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB1s_All_vs_FB2s_All_",sampleName,".pdf"), width = 18, height = 10)
H1
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB2s_All_vs_FB1s_All_",sampleName,".pdf"), width = 18, height = 10)
H2
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB1s_Arach_vs_FB2s_Dura_",sampleName,".pdf"), width = 18, height = 10)
H3
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB2s_Dura_vs_FB1s_Arach_",sampleName,".pdf"), width = 18, height = 10)
H4
dev.off()

####



### Save objects ###
saveRDS(Test,file=paste0(sampleFolder,"results_merge_subset2/Robjects/EnrichGO_FB1s_All_vs_FB2s_All_",sampleName,".rds"))
saveRDS(Test2,file=paste0(sampleFolder,"results_merge_subset2/Robjects/EnrichGO_FB2s_All_vs_FB1s_All_",sampleName,".rds"))
saveRDS(Test_Arach,file=paste0(sampleFolder,"results_merge_subset2/Robjects/EnrichGO_FB1s_Arach_vs_FB2s_Dura_",sampleName,".rds"))
saveRDS(Test2_Dura,file=paste0(sampleFolder,"results_merge_subset2/Robjects/EnrichGO_FB2s_Dura_vs_FB1s_Arach_",sampleName,".rds"))

write.xlsx(Test@result,file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_EnrichGO_FB1s_All_vs_FB2s_All_",sampleName,".xlsx"))
write.xlsx(Test2@result,file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_EnrichGO_FB2s_All_vs_FB1s_All_",sampleName,".xlsx"))
write.xlsx(Test_Arach@result,file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_EnrichGO_FB1s_Arach_vs_FB2s_Dura_",sampleName,".xlsx"))
write.xlsx(Test2_Dura@result,file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_EnrichGO_FB2s_Dura_vs_FB1s_Arach_",sampleName,".xlsx"))

### Read objects ###
Test<-readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/EnrichGO_FB1s_All_vs_FB2s_All_",sampleName,".rds"))
Test2<-readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/EnrichGO_FB2s_All_vs_FB1s_All_",sampleName,".rds"))
Test_Arach<-readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/EnrichGO_FB1s_Arach_vs_FB2s_Dura_",sampleName,".rds"))
Test2_Dura<-readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/EnrichGO_FB2s_Dura_vs_FB1s_Arach_",sampleName,".rds"))

listDEgenesFBs<-readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/Markers_for_our_FBs_paper_",sampleName,".rds"))

###################

### DAVID results ###

## Table 1
DAVID_FB1sAll_vs_FB2sALL<-read.xlsx(paste0(sampleFolder,"results_merge_subset2/DAVID/DAVID_FB1s_All_vs_FB2s_All.xlsx"))
DAVID_FB1sAll_vs_FB2sALL_GO<-DAVID_FB1sAll_vs_FB2sALL[grep("GOTERM",DAVID_FB1sAll_vs_FB2sALL$Category),]
DAVID_FB1sAll_vs_FB2sALL_GO$Type<-"BP"
DAVID_FB1sAll_vs_FB2sALL_GO[grep("CC",DAVID_FB1sAll_vs_FB2sALL_GO$Category),"Type"]<-"CC"
DAVID_FB1sAll_vs_FB2sALL_GO[grep("MF",DAVID_FB1sAll_vs_FB2sALL_GO$Category),"Type"]<-"MF"
colnames(DAVID_FB1sAll_vs_FB2sALL_GO)[4]<-"Percentage"

p1 <- ggplot(head(DAVID_FB1sAll_vs_FB2sALL_GO,50)) +
  geom_point(aes(x = Type, y = reorder(Term,Percentage), size = Percentage, color = Benjamini)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("DAVID GO analysis")

p1 <- p1 + facet_grid(Type~., scale="free")

pdf(file=paste0(sampleFolder,"results_merge_subset2/DAVID/10_Dotplot_DAVID_FB1s_All_vs_FB2s_All_",sampleName,".pdf"), width = 10, height = 15)
p1
dev.off()

## Table 2
DAVID_FB2sAll_vs_FB1sALL<-read.xlsx(paste0(sampleFolder,"results_merge_subset2/DAVID/DAVID_FB2s_All_vs_FB1s_All.xlsx"))
DAVID_FB2sAll_vs_FB1sALL_GO<-DAVID_FB2sAll_vs_FB1sALL[grep("GOTERM",DAVID_FB2sAll_vs_FB1sALL$Category),]
DAVID_FB2sAll_vs_FB1sALL_GO$Type<-"BP"
DAVID_FB2sAll_vs_FB1sALL_GO[grep("CC",DAVID_FB2sAll_vs_FB1sALL_GO$Category),"Type"]<-"CC"
DAVID_FB2sAll_vs_FB1sALL_GO[grep("MF",DAVID_FB2sAll_vs_FB1sALL_GO$Category),"Type"]<-"MF"
colnames(DAVID_FB2sAll_vs_FB1sALL_GO)[4]<-"Percentage"

p2 <- ggplot(head(DAVID_FB2sAll_vs_FB1sALL_GO,50)) +
  geom_point(aes(x = Type, y = reorder(Term,Percentage), size = Percentage, color = Benjamini)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("DAVID GO analysis")

p2 <- p2 + facet_grid(Type~., scale="free")

pdf(file=paste0(sampleFolder,"results_merge_subset2/DAVID/10_Dotplot_DAVID_FB2s_All_vs_FB1s_All_",sampleName,".pdf"), width = 10, height = 15)
p2
dev.off()

## Table 3
DAVID_FB1sArach_vs_FB2sDura<-read.xlsx(paste0(sampleFolder,"results_merge_subset2/DAVID/DAVID_FB1s_Arach_vs_FB2s_Dura.xlsx"))
DAVID_FB1sArach_vs_FB2sDura_GO<-DAVID_FB1sArach_vs_FB2sDura[grep("GOTERM",DAVID_FB1sArach_vs_FB2sDura$Category),]
DAVID_FB1sArach_vs_FB2sDura_GO$Type<-"BP"
DAVID_FB1sArach_vs_FB2sDura_GO[grep("CC",DAVID_FB1sArach_vs_FB2sDura_GO$Category),"Type"]<-"CC"
DAVID_FB1sArach_vs_FB2sDura_GO[grep("MF",DAVID_FB1sArach_vs_FB2sDura_GO$Category),"Type"]<-"MF"
colnames(DAVID_FB1sArach_vs_FB2sDura_GO)[4]<-"Percentage"

p3 <- ggplot(head(DAVID_FB1sArach_vs_FB2sDura_GO,50)) +
  geom_point(aes(x = Type, y = reorder(Term,Percentage), size = Percentage, color = Benjamini)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("DAVID GO analysis")

p3 <- p3 + facet_grid(Type~., scale="free")

pdf(file=paste0(sampleFolder,"results_merge_subset2/DAVID/10_Dotplot_DAVID_FB1s_Arach_vs_FB2s_Dura_",sampleName,".pdf"), width = 10, height = 15)
p3
dev.off()

## Table 4
DAVID_FB2sDura_vs_FB1sArach<-read.xlsx(paste0(sampleFolder,"results_merge_subset2/DAVID/DAVID_FB2s_Dura_vs_FB1s_Arach.xlsx"))
DAVID_FB2sDura_vs_FB1sArach_GO<-DAVID_FB2sDura_vs_FB1sArach[grep("GOTERM",DAVID_FB2sDura_vs_FB1sArach$Category),]
DAVID_FB2sDura_vs_FB1sArach_GO$Type<-"BP"
DAVID_FB2sDura_vs_FB1sArach_GO[grep("CC",DAVID_FB2sDura_vs_FB1sArach_GO$Category),"Type"]<-"CC"
DAVID_FB2sDura_vs_FB1sArach_GO[grep("MF",DAVID_FB2sDura_vs_FB1sArach_GO$Category),"Type"]<-"MF"
colnames(DAVID_FB2sDura_vs_FB1sArach_GO)[4]<-"Percentage"

p4 <- ggplot(head(DAVID_FB2sDura_vs_FB1sArach_GO,50)) +
  geom_point(aes(x = Type, y = reorder(Term,Percentage), size = Percentage, color = Benjamini)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("DAVID GO analysis")

p4 <- p4 + facet_grid(Type~., scale="free")

pdf(file=paste0(sampleFolder,"results_merge_subset2/DAVID/10_Dotplot_DAVID_FB2s_Dura_vs_FB1s_Arach_",sampleName,".pdf"), width = 10, height = 15)
p4
dev.off()

#########################################################################################################
#########################################################################################################

### GOEA v2: with background bulk experiment 26/04/21 ###
library(clusterProfiler)
library(org.Mm.eg.db)
# options(connectionObserver = NULL) #If issue loading package org.MM.eg.db, need to run this line!!!

listDEgenesFBs<- readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/Markers_for_our_FBs_paper_",sampleName,".rds"))
Background_bulk<- read.xlsx("../../DAVID_Background_Genelist_exp1761_BulkRNASEQ_PBSandLPS.xlsx", colNames = F)

## Stromal FBs
Test_background<-enrichGO(
  as.character(listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC > 0),"gene"]),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = Background_bulk$X1,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D1_background<-dotplot(Test_background, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Stalk FBs
Test2_background<-enrichGO(
  as.character(listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC < 0),"gene"]),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = Background_bulk$X1,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D2_background<-dotplot(Test2_background, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Combine
D_comp_background<-merge_result(list(Stromal_CP=Test_background, Stalk_CP=Test2_background)) %>%
  dotplot(., split="ONTOLOGY", showCategory=15) + facet_grid(ONTOLOGY~., scale="free")

# # ## Combine
# # Test_combo<-Test@result
# # Test_combo2<-Test2@result
# # 
# # Test_combo$Type<-"Stromal FBs"
# # Test_combo2$Type<-"Stalk FBs"
# # 
# # Test_combo3<-rbind(Test_combo,Test_combo2)
# # Test3<-Test
# # Test3@result<-Test_combo3
# # Test3@result$Type<-as.character(Test3@result$Type)
# # Test3@result$Description<-as.character(Test3@result$Description)
# # 
# # Test3@gene<-c(Test@gene,Test2@gene)
# # Test3@geneSets<-c(Test@geneSets,Test2@geneSets)
# # 
# # dotplot(Test3, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
# # 
# # dotplot(Test3, split="Type", showCategory=30) + facet_grid(Type~., scale="free")
# # 
# # ## from Tommy's code
# # p <- ggplot(Test3) + 
# #   geom_point(aes(x = Type, y = reorder(Description,GeneRatio), size = GeneRatio, color = p.adjust)) +
# #   theme_bw(base_size = 14) +
# #   scale_colour_gradient(limits=c(0, 0.10), low="red") +
# #   ylab(NULL) +
# #   ggtitle("GO pathway enrichment")
# # 
# # p + facet_grid(Type~., scale="free")
# 
# ## Stromal FBs Arach
# Test_Arach<-enrichGO(
#   as.character(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC > 0),"gene"]),
#   'org.Mm.eg.db',
#   keyType = "SYMBOL",
#   ont = "ALL",
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   universe,
#   qvalueCutoff = 0.2,
#   minGSSize = 10,
#   maxGSSize = 500,
#   readable = FALSE,
#   pool = T
# )
# 
# D3<-dotplot(Test_Arach, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
# 
# ## Stalk FBs Dura
# Test2_Dura<-enrichGO(
#   as.character(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC < 0),"gene"]),
#   'org.Mm.eg.db',
#   keyType = "SYMBOL",
#   ont = "ALL",
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   universe,
#   qvalueCutoff = 0.2,
#   minGSSize = 10,
#   maxGSSize = 500,
#   readable = FALSE,
#   pool = T
# )
# 
# D4<-dotplot(Test2_Dura, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
# 
# ## Combine
# D_comp2<-merge_result(list(Stromal_CP=Test_Arach, Stalk_CP=Test2_Dura)) %>%
#   dotplot(., split="ONTOLOGY", showCategory=15) + facet_grid(ONTOLOGY~., scale="free")


## Check genes associated with GO terms
GeneList_Test_background<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC > 0),"avg_logFC"]
names(GeneList_Test_background)<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC > 0),"gene"]
H1_background<-heatplot(Test_background, foldChange=GeneList_Test_background)

GeneList_Test2_background<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC < 0),"avg_logFC"]*(-1) #Convert logFC
names(GeneList_Test2_background)<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC < 0),"gene"]
H2_background<-heatplot(Test2_background, foldChange=GeneList_Test2_background)

# GeneList_Test_Arach<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC > 0),"avg_logFC"]
# names(GeneList_Test_Arach)<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC > 0),"gene"]
# H3<-heatplot(Test_Arach, foldChange=GeneList_Test_Arach)
# 
# GeneList_Test2_Dura<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC < 0),"avg_logFC"]*(-1) #Convert logFC
# names(GeneList_Test2_Dura)<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC < 0),"gene"]
# H4<-heatplot(Test2_Dura, foldChange=GeneList_Test2_Dura)


#### Adapted comparison plot for FB1 All vs FB2 All (26/03/21)
## Filter
Test_filtered_background<-Test_background
Test2_filtered_background<-Test2_background
Test_filtered_background@result<-Test_filtered_background@result[-which(Test_filtered_background@result$ONTOLOGY == "MF"),]
Test2_filtered_background@result<-Test2_filtered_background@result[-which(Test2_filtered_background@result$ONTOLOGY == "MF"),]
## Combine
D_comp_filtered_background<-merge_result(list(Stromal_CP=Test_filtered_background, Stalk_CP=Test2_filtered_background)) %>%
  dotplot(., split="ONTOLOGY", showCategory=25) + facet_grid(ONTOLOGY~., scale="free")


### Save results ###
pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB1s_All_vs_FB2s_All_background_",sampleName,".pdf"), width = 10, height = 10)
D1_background
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB2s_All_vs_FB1s_All_background_",sampleName,".pdf"), width = 10, height = 10)
D2_background
dev.off()

# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB1s_Arach_vs_FB2s_Dura_",sampleName,".pdf"), width = 10, height = 10)
# D3
# dev.off()
# 
# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB2s_Dura_vs_FB1s_Arach_",sampleName,".pdf"), width = 10, height = 10)
# D4
# dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_background_",sampleName,".pdf"), width = 10, height = 15)
D_comp_background
dev.off()

# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_Arach_and_FB2s_Dura_",sampleName,".pdf"), width = 10, height = 15)
# D_comp2
# dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_background_",sampleName,"_without_MF.pdf"), width = 10, height = 15)
D_comp_filtered_background
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB1s_All_vs_FB2s_All_background_",sampleName,".pdf"), width = 18, height = 10)
H1_background
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB2s_All_vs_FB1s_All_background_",sampleName,".pdf"), width = 18, height = 10)
H2_background
dev.off()

# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB1s_Arach_vs_FB2s_Dura_",sampleName,".pdf"), width = 18, height = 10)
# H3
# dev.off()
# 
# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB2s_Dura_vs_FB1s_Arach_",sampleName,".pdf"), width = 18, height = 10)
# H4
# dev.off()

####


### Did not save objects!!!!!!!! ###


#########################################################################################################
#########################################################################################################

### GOEA v3: with background scRNA-seq experiment 26/04/21 ###
library(clusterProfiler)
library(org.Mm.eg.db)
# options(connectionObserver = NULL) #If issue loading package org.MM.eg.db, need to run this line!!!

listDEgenesFBs<- readRDS(file=paste0(sampleFolder,"results_merge_subset2/Robjects/Markers_for_our_FBs_paper_",sampleName,".rds"))
seuratObj_Urvb<-readRDS(file = "Urvb_datasets/Robjects/seuratObj_final_Urvb_datasets.rds")

Genes_scRnaseq<-rowSums(seuratObj_Urvb@assays$RNA@counts)
Background_scRNAseq<-names(Genes_scRnaseq)

## Stromal FBs
Test_background<-enrichGO(
  as.character(listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC > 0),"gene"]),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = Background_scRNAseq,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D1_background<-dotplot(Test_background, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Stalk FBs
Test2_background<-enrichGO(
  as.character(listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC < 0),"gene"]),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = Background_scRNAseq,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D2_background<-dotplot(Test2_background, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Combine
D_comp_background<-merge_result(list(Stromal_CP=Test_background, Stalk_CP=Test2_background)) %>%
  dotplot(., split="ONTOLOGY", showCategory=15) + facet_grid(ONTOLOGY~., scale="free")

# # ## Combine
# # Test_combo<-Test@result
# # Test_combo2<-Test2@result
# # 
# # Test_combo$Type<-"Stromal FBs"
# # Test_combo2$Type<-"Stalk FBs"
# # 
# # Test_combo3<-rbind(Test_combo,Test_combo2)
# # Test3<-Test
# # Test3@result<-Test_combo3
# # Test3@result$Type<-as.character(Test3@result$Type)
# # Test3@result$Description<-as.character(Test3@result$Description)
# # 
# # Test3@gene<-c(Test@gene,Test2@gene)
# # Test3@geneSets<-c(Test@geneSets,Test2@geneSets)
# # 
# # dotplot(Test3, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
# # 
# # dotplot(Test3, split="Type", showCategory=30) + facet_grid(Type~., scale="free")
# # 
# # ## from Tommy's code
# # p <- ggplot(Test3) + 
# #   geom_point(aes(x = Type, y = reorder(Description,GeneRatio), size = GeneRatio, color = p.adjust)) +
# #   theme_bw(base_size = 14) +
# #   scale_colour_gradient(limits=c(0, 0.10), low="red") +
# #   ylab(NULL) +
# #   ggtitle("GO pathway enrichment")
# # 
# # p + facet_grid(Type~., scale="free")
# 
# ## Stromal FBs Arach
# Test_Arach<-enrichGO(
#   as.character(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC > 0),"gene"]),
#   'org.Mm.eg.db',
#   keyType = "SYMBOL",
#   ont = "ALL",
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   universe,
#   qvalueCutoff = 0.2,
#   minGSSize = 10,
#   maxGSSize = 500,
#   readable = FALSE,
#   pool = T
# )
# 
# D3<-dotplot(Test_Arach, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
# 
# ## Stalk FBs Dura
# Test2_Dura<-enrichGO(
#   as.character(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC < 0),"gene"]),
#   'org.Mm.eg.db',
#   keyType = "SYMBOL",
#   ont = "ALL",
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   universe,
#   qvalueCutoff = 0.2,
#   minGSSize = 10,
#   maxGSSize = 500,
#   readable = FALSE,
#   pool = T
# )
# 
# D4<-dotplot(Test2_Dura, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
# 
# ## Combine
# D_comp2<-merge_result(list(Stromal_CP=Test_Arach, Stalk_CP=Test2_Dura)) %>%
#   dotplot(., split="ONTOLOGY", showCategory=15) + facet_grid(ONTOLOGY~., scale="free")


## Check genes associated with GO terms
GeneList_Test_background<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC > 0),"avg_logFC"]
names(GeneList_Test_background)<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC > 0),"gene"]
H1_background<-heatplot(Test_background, foldChange=GeneList_Test_background)

GeneList_Test2_background<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC < 0),"avg_logFC"]*(-1) #Convert logFC
names(GeneList_Test2_background)<-listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_logFC < 0),"gene"]
H2_background<-heatplot(Test2_background, foldChange=GeneList_Test2_background)

# GeneList_Test_Arach<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC > 0),"avg_logFC"]
# names(GeneList_Test_Arach)<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC > 0),"gene"]
# H3<-heatplot(Test_Arach, foldChange=GeneList_Test_Arach)
# 
# GeneList_Test2_Dura<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC < 0),"avg_logFC"]*(-1) #Convert logFC
# names(GeneList_Test2_Dura)<-listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura[which(listDEgenesFBs$Fib1_Arach_vs_Fib2_Dura$avg_logFC < 0),"gene"]
# H4<-heatplot(Test2_Dura, foldChange=GeneList_Test2_Dura)


#### Adapted comparison plot for FB1 All vs FB2 All (26/03/21)
## Filter
Test_filtered_background<-Test_background
Test2_filtered_background<-Test2_background
Test_filtered_background@result<-Test_filtered_background@result[-which(Test_filtered_background@result$ONTOLOGY == "MF"),]
Test2_filtered_background@result<-Test2_filtered_background@result[-which(Test2_filtered_background@result$ONTOLOGY == "MF"),]
## Combine
D_comp_filtered_background<-merge_result(list(Stromal_CP=Test_filtered_background, Stalk_CP=Test2_filtered_background)) %>%
  dotplot(., split="ONTOLOGY", showCategory=25) + facet_grid(ONTOLOGY~., scale="free")


### Save results ###
pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB1s_All_vs_FB2s_All_background_scRNAseq_",sampleName,".pdf"), width = 10, height = 10)
D1_background
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB2s_All_vs_FB1s_All_background_scRNAseq_",sampleName,".pdf"), width = 10, height = 10)
D2_background
dev.off()

# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB1s_Arach_vs_FB2s_Dura_",sampleName,".pdf"), width = 10, height = 10)
# D3
# dev.off()
# 
# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_FB2s_Dura_vs_FB1s_Arach_",sampleName,".pdf"), width = 10, height = 10)
# D4
# dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_background_scRNAseq_",sampleName,".pdf"), width = 10, height = 15)
D_comp_background
dev.off()

# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_Arach_and_FB2s_Dura_",sampleName,".pdf"), width = 10, height = 15)
# D_comp2
# dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_background_scRNAseq_",sampleName,"_without_MF.pdf"), width = 10, height = 15)
D_comp_filtered_background
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB1s_All_vs_FB2s_All_background_scRNAseq_",sampleName,".pdf"), width = 18, height = 10)
H1_background
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB2s_All_vs_FB1s_All_background_scRNAseq_",sampleName,".pdf"), width = 18, height = 10)
H2_background
dev.off()

# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB1s_Arach_vs_FB2s_Dura_",sampleName,".pdf"), width = 18, height = 10)
# H3
# dev.off()
# 
# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_Heatplot_FB2s_Dura_vs_FB1s_Arach_",sampleName,".pdf"), width = 18, height = 10)
# H4
# dev.off()

####

write.xlsx(Test_background@result,file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_EnrichGO_FB1s_All_vs_FB2s_All_background_scRNAseq_",sampleName,".xlsx"))
write.xlsx(Test2_background@result,file=paste0(sampleFolder,"results_merge_subset2/Annotation/10_EnrichGO_FB2s_All_vs_FB1s_All_background_scRNAseq_",sampleName,".xlsx"))

### Did not save objects!!!!!!!! ###

#########################################################################################################
#########################################################################################################

## Final change to dataset names (2 options -> 2 metadata columns) 26/03/21
## New annotation: combine dataset names
seuratObj_clean@meta.data$sample_origin2_optionA<-seuratObj_clean@meta.data$sample_origin2
seuratObj_clean@meta.data$sample_origin2_optionB<-seuratObj_clean@meta.data$sample_origin2
levels(seuratObj_clean@meta.data$sample_origin2_optionA)<-c("FB_DeSisto_et_al","Ependymal_Zeisel_et_al","Vascular_Vanlandewijck_et_al",
                                                            "Ependymal_Shah_et_al","CP_Verhaege_et_al","Vascular_Zeisel_et_al")
levels(seuratObj_clean@meta.data$sample_origin2_optionB)<-c("FB_DeSisto_et_al","Ependymal_Zeisel_et_al","Vascular_Vanlandewijck_et_al",
                                                            "Ependymal_Shah_et_al","Choroid_Plexus_cells","Vascular_Zeisel_et_al")

library(RColorBrewer)
Colorset<-c(brewer.pal(8,"Set1")[c(1,8,6,2,3,4)])

U_originA<-DimPlot(seuratObj_clean, reduction = "umap", label = F, group.by = "sample_origin2_optionA", label.size = 3, cols = Colorset)
ggsave(U_originA, file=paste0(sampleFolder,"results_merge_subset2/Annotation/11_UMAP_clean_dataset_origin1_optionA_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

U_originB<-DimPlot(seuratObj_clean, reduction = "umap", label = F, group.by = "sample_origin2_optionB", label.size = 3, cols = Colorset)
ggsave(U_originB, file=paste0(sampleFolder,"results_merge_subset2/Annotation/11_UMAP_clean_dataset_origin1_optionB_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

#########################################################################################################

## Check matrixome collagen genes (16/12/21)
Matrixome_df<-read.xlsx("Documentation/MGIBatchReport_matrixome.xlsx", sheet = "Matrixome")
Matrixome_genes<-Matrixome_df$Symbol
intersect(Matrixome_genes,rownames(seuratObj_clean)) #All!

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = Matrixome_genes, cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/12_Dotplot_matrixome_",sampleName,".pdf"), width = 15, height = 10)
D1
dev.off()

# Check overlap with marker genes our FBs (filtered for CPE genes)
Markers_our_FBs<-readRDS("Merge_subset2/results_merge_subset2/Robjects/Markers_for_our_FBs_Merge_FB_datasets_subset2.rds")
intersect(Markers_our_FBs$Filtered_Fib1_vs_all$gene,Matrixome_genes)
intersect(Markers_our_FBs$Filtered_Fib2_vs_all$gene,Matrixome_genes)

#########################################################################################################

## Check genes human BBB two FB types (16/02/22)

Human_FB_genes<-c("LAMA2",'FBLN1','KCNMA1',
                  'SLC4A4','SLC24A3','SLC7A2','SLC39A11','SLC47A1','SLC38A2','SLC22A23','SLC26A2','SLC9B2','SLC2A3','SLC13A3',
                  'SLC35G1','SLC41A2','SLC26A7','SLC1A3','SLC7A11','ABCA9','ABCA10','ABCA8','ABCA6','ABCA11')
## ABCA10 and ABCA11 no mouse genes and ABCA8 is Abca8a according to MGIbatchquery

###First letter upper case
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

Mouse_FB_genes<-firstup(Human_FB_genes)
Mouse_FB_genes[22]<-"Abca8a"
intersect(Mouse_FB_genes,rownames(seuratObj_clean)) #2 lost: ABCA10 and ABCA11!

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = intersect(Mouse_FB_genes,rownames(seuratObj_clean)), 
            cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/14_Dotplot_Human_BBB_paper_",sampleName,".pdf"), width = 15, height = 10)
D1
dev.off()

#########################################################################################################

## FeaturePlot Crym (17/02/22)

Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = "Crym", 
            cols = Colors_dotplot) + RotatedAxis()

F1<-FeaturePlot(object = seuratObj_clean, features = "Crym", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) 

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_Crym_",sampleName,".pdf"), height = 20, width = 15)
F1/D1
dev.off()

## FeaturePlot integrins (11/04/22)
Features <- rownames(seuratObj_clean)[grep("Itg",rownames(seuratObj_clean))]

# PER 36
plots1 <- FeaturePlot(seuratObj_clean, features = Features[1:32],order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots1 <- lapply(X = plots1, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))

pdf(file = paste0(sampleFolder,"results_merge_subset2/Feature_plots/Integrin_Featureplots_q2-q98_per36.pdf"), width = 30, height = 30)
CombinePlots(plots = plots1, ncol=6)
dev.off()

## FeaturePlot Cxcl12 (12/04/22)

Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = "Cxcl12", 
            cols = Colors_dotplot) + RotatedAxis()

F1<-FeaturePlot(object = seuratObj_clean, features = "Cxcl12", cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) 

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plot_Cxcl12_",sampleName,".pdf"), height = 20, width = 15)
F1/D1
dev.off()

#########################################################################################################

## FeaturePlot conference Roos (11/06/22)

Gene_list<-c("Klf5","Tbx18","Lef1","Tjp1","Cdh1","Ctnnb1")

Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = Gene_list, 
            cols = Colors_dotplot) + RotatedAxis()

Idents(seuratObj_clean)<-seuratObj_clean$New_clusters_clean

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plots_and_dotplot_conference_Roos_",sampleName,".pdf"), height = 10, width = 15)
for (feature in Gene_list) {
  F1<-FeaturePlot(object = seuratObj_clean, features =feature, cols = c("yellow", "red"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', label = T, pt.size = 1.5, order=T)
  print(F1)
}
D1
dev.off()

#########################################################################################################

## FeaturePlot request Daan (01/07/22)
# Transporters: Abcb1a (gene for PGP-1, major efflux transporter), Slc7a14, Slc9a2, Slc4a4, Slc4a10, Slc12a2, Slc38a2, Slc26a7, Stra6 (transporter for retinol)
# Trancytosis: Cav1, Cav2, Sdpr (gene for Cavin-2, in human the gene might be CAVIN2)
# Adherens/Tight junctions: Cxadr, Cdh1, Cdh11, Cdh17, Cdh5, Tjp1, Tjp2

Gene_list<-c("Abcb1a", "Slc7a14", "Slc9a2", "Slc4a4", "Slc4a10", "Slc12a2", "Slc38a2", "Slc26a7", "Stra6",
              "Cav1", "Cav2", "Sdpr",
              "Cxadr", "Cdh1", "Cdh11", "Cdh17", "Cdh5", "Tjp1", "Tjp2")

Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = Gene_list, 
            cols = Colors_dotplot) + RotatedAxis()

Idents(seuratObj_clean)<-seuratObj_clean$New_clusters_clean

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plots_and_dotplot_Request_Daan_July2022_",sampleName,".pdf"), height = 10, width = 15)
for (feature in Gene_list) {
  F1<-FeaturePlot(object = seuratObj_clean, features =feature, cols = c("yellow", "red"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', label = T, pt.size = 1.5, order=T)
  print(F1)
}
D1
dev.off()

## Vimentin
Gene_list<-c("Vim")

Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = Gene_list, 
            cols = Colors_dotplot) + RotatedAxis()

Idents(seuratObj_clean)<-seuratObj_clean$New_clusters_clean

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plots_and_dotplot_Request2_Daan_July2022_",sampleName,".pdf"), height = 10, width = 15)
for (feature in Gene_list) {
  F1<-FeaturePlot(object = seuratObj_clean, features =feature, cols = c("yellow", "red"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', label = T, pt.size = 1.5, order=T)
  print(F1)
}
D1
dev.off()

#########################################################################################################

## Check FB marker genes (01/09/22) (Update 01/03/23)
Figure_genes<-c("Steap4", "Heph", "Cp", "Cxcl12", "Dpep1", "Itga8", "Itgbl1", "Tjp1", "Gja1", "Gjb2", "Gjb6", "Cdh1", "Cdh5", "Cldn11")

Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean
levels(Idents(seuratObj_clean))

# seuratObj_subset<-subset(seuratObj_clean, idents = c("CP Stromal", "CP Stalk"))
seuratObj_subset<-subset(seuratObj_clean, idents = c("CP Stromal", "CP Base"))

# seuratObj_subset$New_clusters_clean<-factor(as.character(seuratObj_subset$New_clusters_clean), levels = c("CP Stalk","CP Stromal"))
seuratObj_subset$New_clusters_clean<-factor(as.character(seuratObj_subset$New_clusters_clean), levels = c("CP Base","CP Stromal"))

levels(seuratObj_subset$New_clusters_clean)<-c("Base barrier cells","Stromal")

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_subset, group.by = "New_clusters_clean", features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/16_Dotplot_figure3B_",sampleName,".pdf"), width = 10, height = 6)
D1
dev.off()

# Bigger dotplot (01/03/23)
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_subset, group.by = "New_clusters_clean", features = rev(Figure_genes), cols = Colors_dotplot,
            dot.scale = 14) + RotatedAxis() #dot.scale = 8

# pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/16_Dotplot_figure4B_bigger_dots_v2_",sampleName,".pdf"), width = 7.5, height = 4.5)
pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/16_Dotplot_figure4B_bigger_dots_v1_",sampleName,".pdf"), width = 10, height = 6)
D1
dev.off()

#########################################################################################################

## Check FB marker genes (29/11/22)
Figure_genes<-c("Tbx18", "Wnt1", "Mesp1")

Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean
levels(Idents(seuratObj_clean))

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Annotation/17_Dotplot_extra_29112022_",sampleName,".pdf"), width = 10, height = 6)
D1
dev.off()

#########################################################################################################

## Check Ccl20 expr (15/12/22)
Figure_genes<-c("Ccl20")

Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean
levels(Idents(seuratObj_clean))

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plots_and_dotplot_Ccl20_",sampleName,".pdf"), height = 10, width = 15)
for (feature in Figure_genes) {
  F1<-FeaturePlot(object = seuratObj_clean, features =feature, cols = c("yellow", "red"), 
                  reduction = "umap", label = T, pt.size = 1.5, order=T) # Turned off cutoffs to show the 1 positive cell!!!!!!!!!!
  print(F1)
}
D1
dev.off()

#########################################################################################################

## Check gene list St Louis Daan expr (24/01/23)
Figure_genes<-c("Sema3a", "Sema3d", "Slit2", "Slit3", "Ntn1", "Efna5", "Efnb2", "Dpp4")

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean
levels(Idents(seuratObj_clean))

D2<-DotPlot(seuratObj_clean, features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plots_and_dotplots_St_Louis_",sampleName,".pdf"), height = 10, width = 15)
D1
D2

for (feature in Figure_genes) {
  F1<-FeaturePlot(object = seuratObj_clean, features =feature, cols = c("yellow", "red"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', label = T, pt.size = 1.5, order=T)
  print(F1)
}

dev.off()

#########################################################################################################

## Check gene Daan paper (27/02/23)
Figure_genes<-c("Foxc1")

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean
levels(Idents(seuratObj_clean))

D2<-DotPlot(seuratObj_clean, features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plots_and_dotplots_Foxc1_",sampleName,".pdf"), height = 10, width = 15)
D1
D2

for (feature in Figure_genes) {
  F1<-FeaturePlot(object = seuratObj_clean, features =feature, cols = c("yellow", "red"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', label = T, pt.size = 1.5, order=T)
  print(F1)
}

dev.off()

#########################################################################################################

## Daan paper (01/03/23)
Figure_genes<-c("Aqp1","Aqp4")

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean
levels(Idents(seuratObj_clean))

D2<-DotPlot(seuratObj_clean, features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset2/Feature_plots/Feature_plots_and_dotplots_Aqp_",sampleName,".pdf"), height = 10, width = 15)
D1
D2

for (feature in Figure_genes) {
  F1<-FeaturePlot(object = seuratObj_clean, features =feature, cols = c("yellow", "red"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', label = T, pt.size = 1.5, order=T)
  print(F1)
}

dev.off()
