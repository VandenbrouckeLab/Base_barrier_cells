## Analysis with 6 Urvb dataset (non-LPS)
## To perform check between ages without any other datasets confounding analysis

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')

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

sampleName <- "Urvb_FB_6datasets" #Change for this analysis!!!
sampleFolder<-paste0("Merge_Urvb_extra","/")

##add some subfolders
dir.create(paste0(sampleFolder,"results"))
dir.create(paste0(sampleFolder,"results/QC"))
dir.create(paste0(sampleFolder,"results/Robjects"))

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

## Add extra metaData
metaData_FBs <- readRDS(file = "Urvb_6datasets/Robjects/MetaData_FBs.rds")
colnames(metaData_FBs) <- c("Dataset", "annotated_clusters_original")
seuratObj@meta.data<-cbind(seuratObj@meta.data, metaData_FBs)

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

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results/QC/1_UMI_",sampleName,".png"), width = 20)

########## mito.genes plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'percent.mito',"mito")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'percent.mito',"mito")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results/QC/2_percMito_",sampleName,".png"), width = 20)


########## PCA plot ##########
# pdf(file=paste0(sampleFolder,"results/QC/13a_PCA_",sampleName,".pdf"), width=10)
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,3))
# dev.off()

#RNA clusters
dir.create(paste0(sampleFolder,"results/Annotation"))

pdf(file=paste0(sampleFolder,"results/Annotation/1_annotation_Color_RNA_clusters_on_harmony_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$harmony_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$harmony_clusters==i),])))
  C1<-C1+ggtitle(paste0("Harmony_cluster_",i))
  print(C1)
}
dev.off()

################################################################################
########## CHECK DE GENES
################################################################################
dir.create(paste0(sampleFolder,"results/Feature_plots"))

##### Epithelial marker ->
Features<-c("Otx2", "Ttr")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"CPE"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Fras1 marker ->
Features<-c("Fras1", "Xist","Neat1","Steap2","Folr1","Car2")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"CPE_extra"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order =T)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_ordered_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Endothelial marker ->
Features<-c("Pecam1", "Flt1","Plvap")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"EC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Vascular associated marker ->
Features<-c("Pdgfrb", "Mylk","Myh11","Tagln")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"VAC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Fibroblast marker ->
Features<-c("Dcn", "Col1a1")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"FB"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Macrophage marker ->
Features<-c("Adgre1", "Csf1r","Fcgr1")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"MF"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Microglia marker ->
Features<-c("P2ry12")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Microglia-like"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### NK cell marker ->
Features<-c("Klrb1c", "Gzmb")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"NK"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### cDC marker ->
Features<-c("Cd209a", "Ccr7","Xcr1")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"DC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Neutrophil marker ->
Features<-c("S100a8", "Ngp","Retnlg")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"NF"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Basophil marker ->
FeaturePlot(object = seuratObj, features = "Fcer1a", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Fcer1g", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = c("Fcer1a", "Fcer1g"), cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

##### Monocyte marker ->
Features<-c("Plac8", "Ly6c2","Ccr2")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"MC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Mast cell marker ->
FeaturePlot(object = seuratObj, features = "Ptgds", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

##### Neuronal + glial cell marker ->
Features<-c("Tubb3","Slc1a3", "Fabp7","Olig1") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Neuronal_and_glial"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Mitotic cell marker ->
Features<-c("Birc5", "Mki67")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Mitotic_cells"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Dissociation effect ->
Features<-c("Fos", "Junb","Atf3","Dusp1","Ccl4")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Dissociation"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


#########################################################################################################################################

########################################
##### all clusters vs all clusters
########################################

# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

dir.create(paste0(sampleFolder,"results/Marker_lists"))

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
RNAMarkers_RNAclus <- FindAllMarkers(seuratObj, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"results/Robjects/RNAmarkersList_RNAclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerRNAcluster']]<-paste0(table(RNAMarkers_RNAclus$cluster)," RNA markers for RNA cluster ",rownames(table(RNAMarkers_RNAclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results/Marker_lists/RNAmarkersList_RNAclus_",sampleName,".xlsx"))

########################################
dir.create(paste0(sampleFolder,"results/Heatmaps"))

# ## Load in markers
# RNAMarkers_RNAclus<- readRDS(file=paste0(sampleFolder,"results/Robjects/RNAmarkersList_RNAclus_",sampleName,".rds"))
#
# ## Perform on a subset -> better view of smaller clusters!!
# seuratObj.small <- subset(seuratObj, downsample = 500)
#
# ## Heatmap RNA markers on RNA clusters
# top10 <- RNAMarkers_RNAclus %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# D1<-DoHeatmap(seuratObj.small, features = top10$gene, group.by = "RNA_clusters") + NoLegend()
# ggsave(D1, file=paste0(sampleFolder, "results/Heatmaps/Heatmap_RNAmarkersList_RNAclus_",sampleName,".png"),
#        height = 20, width = 12, dpi = "retina")
#
# pdf(file=paste0(sampleFolder, "results/Heatmaps/Heatmap_RNAmarkersList_RNAclus_",sampleName,".pdf"),
#     height = 20, width = 25)
# DoHeatmap(seuratObj.small, features = top10$gene, group.by = "RNA_clusters") + NoLegend()
# dev.off()


# #############################
# ### Extra detail clusters
# ############################
# Idents(seuratObj)<-seuratObj@meta.data$harmony_clusters
# 
# detail1_16vs9<-FindMarkers(seuratObj, ident.1 = 16, ident.2 = 9, min.pct = 0.10,
#                            min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)
# 
# detail2_1vs9<-FindMarkers(seuratObj, ident.1 = 1, ident.2 = 9, min.pct = 0.10,
#                           min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)
# 
# detail3_8vs1.16.9<-FindMarkers(seuratObj, ident.1 = 8, ident.2 = c(1,9,16), min.pct = 0.10,
#                                min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)
# 
# 
# ##### Create list
# listDEgenesExtra<-tibble::lst(detail1_16vs9, detail2_1vs9, detail3_8vs1.16.9)
# 
# ##Add geneSymbol in column (for the export)
# listDEgenesExtra<-lapply(listDEgenesExtra,function(x){x<-cbind(x,'gene'=rownames(x))})
# ##Filter on adj.P-value
# listDEgenesExtra<-lapply(listDEgenesExtra, function(x){dplyr::filter(x, p_val_adj<0.01)})
# ##Add score
# listDEgenesExtra<-lapply(listDEgenesExtra, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
#                                                              x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
# # listDEgenesExtra<-lapply(listDEgenesExtra, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
# ##Sort on logFC
# listDEgenesExtra<-lapply(listDEgenesExtra,function(x){x<-x[order(x$score, decreasing=T),]})
# 
# saveRDS(listDEgenesExtra,file=paste0(sampleFolder,"results/Robjects/detailClusters_sliced_",sampleName,".rds"))
# 
# ##write to Excel
# library('openxlsx')
# write.xlsx(listDEgenesExtra, paste0(sampleFolder,"results/Marker_lists/detailClusters_sliced_",sampleName,".xlsx"))

########################################################################################################################

## Check other annotation

DimPlot(seuratObj, reduction = "umap", label = T, group.by = "Dataset", label.size = 4)

DimPlot(seuratObj, reduction = "umap", label = T, group.by = "annotated_clusters_original", label.size = 4)

pdf(file=paste0(sampleFolder,"results/Annotation/1_UMAP_split_datasets_",sampleName,".pdf"), width = 40, height = 30)
DimPlot(seuratObj, reduction = "umap", label = T, split.by = "Dataset", group.by = "annotated_clusters_original", label.size = 4, ncol = 3)
dev.off()

U_old_annot<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters_original", label.size = 3)
ggsave(U_old_annot, file=paste0(sampleFolder,"results/Annotation/1_UMAP_old_annot_",sampleName,".png"), height = 8, width = 15, dpi = "retina")


########################################################################################################################

########################################
##### RNA clusters post-slice (after split 10!!)
########################################
########################################

# seuratObj@meta.data$sliced_clusters<- factor(seuratObj@meta.data$sliced_clusters,sort(as.numeric(levels(seuratObj@meta.data$sliced_clusters)))) #reorder levels
# seuratObj@meta.data$annotated_clusters3 <- factor(seuratObj@meta.data$annotated_clusters3,sort(as.numeric(levels(seuratObj@meta.data$annotated_clusters3)))) #reorder levels
# Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters
# 
# U_annot<-DimPlot(seuratObj, reduction = "umap", label = T, group.by = "sliced_clusters", label.size = 4)
# ggsave(U_annot, file=paste0(sampleFolder,"results/Annotation/2_UMAP_sliced_clusters_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

# ### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
# RNAMarkers_RNAclus <- FindAllMarkers(seuratObj, assay = "RNA", only.pos = TRUE)
# table(RNAMarkers_RNAclus$cluster)
# saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"results/Robjects/RNAmarkersList_RNAclus_",sampleName,"_sliced.rds"))
# 
# ### Add to diagnostics
# diagnostics[['RNAmarkersPerRNAclustersliced']]<-paste0(table(RNAMarkers_RNAclus$cluster)," RNA markers for RNA cluster sliced ",rownames(table(RNAMarkers_RNAclus$cluster)))
# saveRDS(diagnostics, file=paste0(sampleFolder,"results/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))
# 
# ### Create list with markers
# totalNrRNAclusters_RNAclus<-max(as.numeric(names(table(RNAMarkers_RNAclus$cluster))))
# totalNrRNAclusters_RNAclusPlusOne<-totalNrRNAclusters_RNAclus+1
# RNAmarkersList_RNAclus<-list()
# 
# for(i in 1:totalNrRNAclusters_RNAclusPlusOne){
#   clusterNr<-i-1
#   
#   tmp<-RNAMarkers_RNAclus[RNAMarkers_RNAclus$cluster==clusterNr,]
#   tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
#   
#   RNAmarkersList_RNAclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
# }
# names(RNAmarkersList_RNAclus)<-paste0("RNAclustersliced",0:totalNrRNAclusters_RNAclus)
# 
# ### Write to Excel
# library('openxlsx')
# write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results/Marker_lists/RNAmarkersList_RNAclus_",sampleName,"_sliced.xlsx"))


########################################
##### Markers annotated clusters
########################################
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident
seuratObj@meta.data$annotated_clusters_detailed <- seuratObj@active.ident
levels(seuratObj@meta.data$annotated_clusters) <- c("Stromal Fibroblasts","Stromal Fibroblasts","Stromal Fibroblasts","Stalk Fibroblasts",
                                                    "Stromal Fibroblasts","CPE contaminated FBs","MF contaminated FBs",
                                                    "Stromal Fibroblasts","EC contaminated FBs","Stromal Fibroblasts", "Stalk Fibroblasts")
levels(seuratObj@meta.data$annotated_clusters_detailed) <- c("Stromal Fibroblasts 1","Klf4+ Stromal Fibroblasts","High Mito Stromal Fibroblasts","Stalk Fibroblasts",
                                                             "Stromal Fibroblasts 2","CPE contaminated FBs","MF contaminated FBs",
                                                             "High Mito Klf4+ Stromal Fibroblasts","EC contaminated FBs","ISG+ Stromal Fibroblasts", "ABCs")

U_annot<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
U_annot2<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters_detailed", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results/Annotation/2_UMAP_annotated1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")
ggsave(U_annot2, file=paste0(sampleFolder,"results/Annotation/2_UMAP_annotated2_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters

# Clean outliers and subset object!!!!
U1 <- DimPlot(seuratObj, reduction = "umap", label = T, label.size = 4)
seuratObjNew <- CellSelector(U1, object=seuratObj, ident="CPE contaminated FBs")

U1 <- DimPlot(seuratObjNew, reduction = "umap", label = T, label.size = 4)
seuratObjNew <- CellSelector(U1, object=seuratObjNew, ident="MF contaminated FBs")

U1 <- DimPlot(seuratObjNew, reduction = "umap", label = T, label.size = 4)
seuratObjNew <- CellSelector(U1, object=seuratObjNew, ident="EC contaminated FBs")

seuratObjNew<-subset(seuratObjNew, idents=c("Stromal Fibroblasts","Stalk Fibroblasts"))

DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters_detailed", label.size = 4)

U_annot_subset<-DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
U_annot2_subset<-DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters_detailed", label.size = 4)
ggsave(U_annot_subset, file=paste0(sampleFolder,"results/Annotation/4_UMAP_annotated1_subset_",sampleName,".png"), height = 10, width = 10, dpi = "retina")
ggsave(U_annot2_subset, file=paste0(sampleFolder,"results/Annotation/4_UMAP_annotated2_subset_",sampleName,".png"), height = 10, width = 12, dpi = "retina")

## Update factors
seuratObjNew@meta.data$annotated_clusters <-as.factor(as.character(seuratObjNew@meta.data$annotated_clusters))
seuratObjNew@meta.data$annotated_clusters_detailed <- as.factor(as.character(seuratObjNew@meta.data$annotated_clusters_detailed))

## New annotation
seuratObjNew@meta.data$annotated_clusters_new <- seuratObjNew@meta.data$annotated_clusters_detailed
levels(seuratObjNew@meta.data$annotated_clusters_new) <- c("ABCs","Stromal Fibroblasts","Stromal Fibroblasts",
                                                           "Stromal Fibroblasts", "Stromal Fibroblasts","Stalk Fibroblasts",                  
                                                           "Stromal Fibroblasts", "Stromal Fibroblasts")
seuratObjNew@meta.data$annotated_clusters_new<-factor(seuratObjNew@meta.data$annotated_clusters_new, levels = levels(seuratObjNew@meta.data$annotated_clusters_new)[c(1,3,2)])

library(RColorBrewer)
Colorset<-c(brewer.pal(3,"Set2"))

U_annot_subset3<-DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters_new", cols = Colorset, label.size = 4)
ggsave(U_annot_subset3, file=paste0(sampleFolder,"results/Annotation/4_UMAP_annotated3_subset_",sampleName,".png"), height = 10, width = 10, dpi = "retina")

###############################################################

## Save subset object
saveRDS(seuratObjNew, file=paste0(sampleFolder,"Robjects/seuratObj_subset_",sampleName,"_harmony_RNA.rds"))

## Read subset object
seuratObjNew <- readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_subset_",sampleName,"_harmony_RNA.rds"))

###############################################################

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

Idents(seuratObjNew)<-seuratObjNew@meta.data$annotated_clusters_detailed

RNAMarkers_SCTclus <- FindAllMarkers(seuratObjNew, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_SCTclus$cluster)
saveRDS(RNAMarkers_SCTclus, file=paste0(sampleFolder,"results/Robjects/RNAmarkersList_SCTclus_subset_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerSCTclusterannotatedSubset']]<-paste0(table(RNAMarkers_SCTclus$cluster)," RNA markers for annotated cluster ",rownames(table(RNAMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(RNAmarkersList_SCTclus, file =paste0(sampleFolder, "results/Marker_lists/RNAmarkersList_SCTclus_subset_",sampleName,"_annotated.xlsx"))

######################################################

### Make heatmap for annotated clusters
RNAMarkers_SCTclus<- readRDS(file=paste0(sampleFolder,"results/Robjects/RNAmarkersList_SCTclus_subset_",sampleName,"_annotated.rds"))

## Perform on a subset -> better view of smaller clusters!!
seuratObjNew.small <- subset(seuratObjNew, downsample = 500)

########## Get HVG ##########
seuratObjNew.small <- FindVariableFeatures(object = seuratObjNew.small, selection.method = "vst", nfeatures = nrow(seuratObjNew.small@assays$RNA))
length(VariableFeatures(seuratObjNew.small))

########## Scale ##########
seuratObjNew.small <- ScaleData(seuratObjNew.small)

## Heatmap RNA markers on RNA clusters
top10 <- RNAMarkers_SCTclus %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
D1<-DoHeatmap(seuratObjNew.small, features = top10$gene, group.by = "annotated_clusters_detailed") + NoLegend()
ggsave(D1, file=paste0(sampleFolder, "results/Heatmaps/Heatmap_RNAmarkersList_Annotatedclus_",sampleName,"_2.png"), 
       height = 20, width = 12, dpi = "retina")

pdf(file=paste0(sampleFolder, "results/Heatmaps/Heatmap_RNAmarkersList_Annotatedclus_",sampleName,"_2.pdf"), 
    height = 25, width = 25)
DoHeatmap(seuratObjNew.small, features = top10$gene, group.by = "annotated_clusters_detailed") + NoLegend()
dev.off()


######################################################

# Frequency tables (sliced)
Sample <- seuratObjNew@meta.data$orig.ident
cluster <- seuratObjNew@meta.data$annotated_clusters_detailed
Aggr <- rep(sampleName,length(cluster))

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

# Stacked
library("ggthemes")

png(file=paste0(sampleFolder,"results/Annotation/3_SampleDistribution_ggplot2_annot_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"results/Annotation/3_SampleDistribution_ggplot2_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png(file=paste0(sampleFolder,"results/Annotation/3_SampleDistribution_ggplot2_annot_3_",sampleName,".png"), width = 2000, height = 1500, res = 300)
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
# png(file=paste0(sampleFolder,"results/Annotation/3_SampleDistribution_ggplot2_annot_1_",sampleName,"_adjusted.png"), width = 2000, height = 1500, res = 300)
# ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + scale_fill_manual(values=c('Blue','Red',"Turquoise","Magenta")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
#   geom_bar(position="fill", stat="identity", colour="white")
# dev.off()

######################################################


## Markers Daan volcanoplot (Jan 2022)

### Create new clusters: split on source
seuratObjNew@meta.data$newClustersTmp<-seuratObjNew@meta.data$annotated_clusters_new
seuratObjNew@meta.data$newClusters<-paste0(seuratObjNew@meta.data$newClustersTmp,"_",seuratObjNew@meta.data$Age)
head(seuratObjNew@meta.data)

### Use the new clusters
seuratObjNew@meta.data$newClusters<- as.factor(seuratObjNew@meta.data$newClusters) #reorder levels
seuratObjNew@meta.data$newClusters<-factor(seuratObjNew@meta.data$newClusters, levels = levels(seuratObjNew@meta.data$newClusters)[c(2,1,3,5,4,6,8,7,9)]) #reorder levels
Idents(seuratObjNew)<-seuratObjNew@meta.data$newClusters

## Create dimplot
U_new_subset<-DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "newClusters", label.size = 4)
ggsave(U_new_subset, file=paste0(sampleFolder,"results/Annotation/4_UMAP_NewClusters_subset_",sampleName,".png"), height = 10, width = 10, dpi = "retina")


########## 2. GET MARKERS (everything!!) ##########
getDEgenes<-function(ident1, ident2){
  markersDiff <- FindMarkers(seuratObjNew, ident.1 = ident1, ident.2 = ident2,
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

levels(Idents(seuratObjNew))[grep("Stalk",levels(Idents(seuratObjNew)))] 
levels(Idents(seuratObjNew))[grep("Stromal",levels(Idents(seuratObjNew)))] 

# Clusters to compare!!!!
"Stalk Fibroblasts_7w"       
"Stalk Fibroblasts_22w"       
"Stalk Fibroblasts_82w"       

Stalk_82w_vs_Stalk_7w_volcano<-getDEgenes("Stalk Fibroblasts_82w","Stalk Fibroblasts_7w")
Stalk_82w_vs_Stalk_7w_volcano<-Stalk_82w_vs_Stalk_7w_volcano[order(Stalk_82w_vs_Stalk_7w_volcano$avg_logFC,decreasing = T),]
head(Stalk_82w_vs_Stalk_7w_volcano)
dim(Stalk_82w_vs_Stalk_7w_volcano)

Stalk_22w_vs_Stalk_7w_volcano<-getDEgenes("Stalk Fibroblasts_22w","Stalk Fibroblasts_7w")
Stalk_22w_vs_Stalk_7w_volcano<-Stalk_22w_vs_Stalk_7w_volcano[order(Stalk_22w_vs_Stalk_7w_volcano$avg_logFC,decreasing = T),]
head(Stalk_22w_vs_Stalk_7w_volcano)
dim(Stalk_22w_vs_Stalk_7w_volcano)

Stalk_82w_vs_Stalk_22w_volcano<-getDEgenes("Stalk Fibroblasts_82w","Stalk Fibroblasts_22w")
Stalk_82w_vs_Stalk_22w_volcano<-Stalk_82w_vs_Stalk_22w_volcano[order(Stalk_82w_vs_Stalk_22w_volcano$avg_logFC,decreasing = T),]
head(Stalk_82w_vs_Stalk_22w_volcano)
dim(Stalk_82w_vs_Stalk_22w_volcano)


##add to list
listDiffMarkers_volcano<-tibble::lst(Stalk_82w_vs_Stalk_7w_volcano,Stalk_22w_vs_Stalk_7w_volcano,
                                     Stalk_82w_vs_Stalk_22w_volcano)

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

saveRDS(listDiffMarkers_volcano,file=paste0(sampleFolder,"results/Robjects/Markers_Stalk_FBs_volcanoplot_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDiffMarkers_volcano, paste0(sampleFolder,"results/Marker_lists/Markers_Stalk_FBs_volcanoplot_",sampleName,".xlsx"))

#########
library(EnhancedVolcano)

list_names_volc<-c("Stalk FBs 82w vs 7w","Stalk FBs 22w vs 7w","Stalk FBs 82w vs 22w")

for (k in 1:length(listDiffMarkers_volcano)) {
  E1<-EnhancedVolcano(listDiffMarkers_volcano[[k]],
                      lab = listDiffMarkers_volcano[[k]][["gene"]],
                      x = "avg_logFC",
                      y = "p_val_adj",
                      xlim = c(-8, 8),
                      title = paste0(list_names_volc[[k]]),
                      pCutoff = 0.01,
                      FCcutoff = 0.25,
                      pointSize = 2.0,
                      labSize = 3.0)
  ggsave(E1,file=paste0(sampleFolder,"results/Annotation/4_Volcanoplot_",names(listDiffMarkers_volcano)[k],"_",sampleName,".png"),dpi = 300, width=10, height=10)
}


## Check Igfbp6 and Cldn11 expression

# Featureplots on full object
# Module score Seurat
library(viridis)

FB_stalk_signature <- list(c("Cldn11","Igfbp6"))

seuratObjNew <- AddModuleScore(object = seuratObjNew, features = FB_stalk_signature, name = "FB_stalk_signature_score")
F2 <- FeaturePlot(object = seuratObjNew, features = "FB_stalk_signature_score1", order = T)  + scale_color_viridis(option = "C")
pdf(file = paste0(sampleFolder,"results/Annotation/6_UMAP_Featureplot_modulescore_Cldn11_Igfbp6_",sampleName,".pdf"), width = 8, height = 10)
F2
dev.off()

Idents(seuratObjNew)<-seuratObjNew@meta.data$annotated_clusters
F2.5 <- FeaturePlot(object = seuratObjNew, features = "FB_stalk_signature_score1", order = T, label = T, repel = T, label.size = 3, cols = c("Yellow","Red"))
pdf(file = paste0(sampleFolder,"results/Annotation/6_UMAP_Featureplot_modulescore_Cldn11_Igfbp6_labeled_",sampleName,".pdf"), width = 8, height = 10)
F2.5
dev.off()

# Nebulosa code (in R v4)
seuratObjNew<-UpdateSeuratObject(seuratObjNew)

F3<-plot_density(seuratObjNew, features =  c("Cldn11","Igfbp6"), slot = "data", reduction = "umap",joint = TRUE)
pdf(file = paste0(sampleFolder,"results/Annotation/6_UMAP_Nebulosa_Cldn11_Igfbp6_",sampleName,".pdf"), width = 20, height = 8)
F3
dev.off()

# Subset stalk and violin plot
Idents(seuratObjNew)<-seuratObjNew@meta.data$newClusters
seuratObj_Stalk_FBs<-subset(seuratObjNew, ident = c("Stalk Fibroblasts_7w", "Stalk Fibroblasts_22w","Stalk Fibroblasts_82w"))

V1<-VlnPlot(seuratObj_Stalk_FBs, features = c("Igfbp6","Cldn11"))

ggsave(V1, file=paste0(sampleFolder,"results/Annotation/6_Test_violinplot_",sampleName,".png"), width =10, height= 8, dpi = "retina")

V2<-VlnPlot(seuratObj_Stalk_FBs, features = c("Dpep1","Alpl"))

ggsave(V2, file=paste0(sampleFolder,"results/Annotation/7_Violinplot_extra_Stalk_",sampleName,".png"), width =10, height= 8, dpi = "retina")

# Subset stromal
seuratObj_Stromal_FBs<-subset(seuratObjNew, ident = c("Stromal Fibroblasts_7w", "Stromal Fibroblasts_22w","Stromal Fibroblasts_82w"))

V3<-VlnPlot(seuratObj_Stromal_FBs, features = c("Dpep1","Alpl"))

ggsave(V3, file=paste0(sampleFolder,"results/Annotation/7_Violinplot_extra_Stromal_",sampleName,".png"), width =10, height= 8, dpi = "retina")

# Extra featureplots
F1 <- FeaturePlot(object = seuratObjNew, features = c("Cldn11","Igfbp6","Dpep1","Alpl"), order = T, pt.size = 1, cols = c("Yellow","Red"), combine = F)
pdf(file = paste0(sampleFolder,"results/Annotation/7_UMAP_Featureplots_FB_markers_",sampleName,".pdf"), width = 8, height = 10)
F1
dev.off()


# Subset and violin plot extra (18/01/22) -> Stromal + Stalk
Idents(seuratObjNew)<-seuratObjNew@meta.data$newClusters
seuratObj_Stalk_and_Stromal<-subset(seuratObjNew, ident = c("Stalk Fibroblasts_7w", "Stalk Fibroblasts_22w","Stalk Fibroblasts_82w",
                                                            "Stromal Fibroblasts_7w", "Stromal Fibroblasts_22w","Stromal Fibroblasts_82w"))

seuratObj_Stalk_and_Stromal@active.ident <- factor(seuratObj_Stalk_and_Stromal@active.ident, levels = levels(Idents(seuratObj_Stalk_and_Stromal))[c(1,4,2,5,3,6)])

library(RColorBrewer)
ColorSet<-brewer.pal(n = 8, name = "Paired")[3:8]
V_final2<-VlnPlot(seuratObj_Stalk_and_Stromal, features = c("Igfbp6","Cldn11","Dpep1","Alpl"), cols = ColorSet, ncol = 4) & theme(plot.margin = margin(0.5,0.5,0.5,1, "cm"))

pdf(file=paste0(sampleFolder,"results/Annotation/8_Violinplot_final_v2_",sampleName,".pdf"), width =37.5, height= 10)
V_final2
dev.off()


#################################################################

levels(Idents(seuratObjNew))[grep("Stromal",levels(Idents(seuratObjNew)))] 

# Clusters to compare!!!!
"Stromal Fibroblasts_7w"          
"Stromal Fibroblasts_22w"          
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
Stromal_82w_vs_Stromal_22w_volcano<-Stromal_82w_vs_Stromal_22w_volcano[order(Stromal_82w_vs_Stromal_22w_volcano$avg_logFC,decreasing = T),]
head(Stromal_82w_vs_Stromal_22w_volcano)
dim(Stromal_82w_vs_Stromal_22w_volcano)


##add to list
listDiffMarkers_volcano_stromal<-tibble::lst(Stromal_82w_vs_Stromal_7w_volcano,Stromal_22w_vs_Stromal_7w_volcano,
                                             Stromal_82w_vs_Stromal_22w_volcano)

lapply(listDiffMarkers_volcano_stromal, dim)

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

saveRDS(listDiffMarkers_volcano_stromal,file=paste0(sampleFolder,"results/Robjects/Markers_Stromal_FBs_volcanoplot_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDiffMarkers_volcano_stromal, paste0(sampleFolder,"results/Marker_lists/Markers_Stromal_FBs_volcanoplot_",sampleName,".xlsx"))

#########
library(EnhancedVolcano)

list_names_volc_stromal<-c("Stromal FBs 82w vs 7w","Stromal FBs 22w vs 7w","Stromal FBs 82w vs 22w")

for (k in 1:length(listDiffMarkers_volcano_stromal)) {
  E1<-EnhancedVolcano(listDiffMarkers_volcano_stromal[[k]],
                      lab = listDiffMarkers_volcano_stromal[[k]][["gene"]],
                      x = "avg_logFC",
                      y = "p_val_adj",
                      xlim = c(-8, 8),
                      title = paste0(list_names_volc_stromal[[k]]),
                      pCutoff = 0.01,
                      FCcutoff = 0.25,
                      pointSize = 2.0,
                      labSize = 3.0)
  ggsave(E1,file=paste0(sampleFolder,"results/Annotation/5_Volcanoplot_",names(listDiffMarkers_volcano_stromal)[k],"_",sampleName,".png"),dpi = 300, width=10, height=10)
}


#################################################################

## Check other dataset for amount of DEGs for similar amount of cells Stalk FBs
table(seuratObjNew@meta.data$newClusters)

# 200 Stalk FBs vs 1k+ Stromal FBs!!
# ABCs_7w                ABCs_22w                ABCs_82w 
# 10                       9                      11 
# Stalk Fibroblasts_7w   Stalk Fibroblasts_22w   Stalk Fibroblasts_82w 
# 70                      68                     104 
# Stromal Fibroblasts_7w Stromal Fibroblasts_22w Stromal Fibroblasts_82w 
# 399                     606                     921 


## More than 25 ABCs -> check DEGs vs stalk FBs
Idents(seuratObjNew)<-seuratObjNew@meta.data$annotated_clusters_new

ABCs_vs_Stalk<-getDEgenes("ABCs","Stalk Fibroblasts")
ABCs_vs_Stalk<-ABCs_vs_Stalk[order(ABCs_vs_Stalk$avg_logFC,decreasing = T),]
head(ABCs_vs_Stalk)
dim(ABCs_vs_Stalk)

##add to list
listDiffMarkers_ABCs_vs_Stalk<-tibble::lst(ABCs_vs_Stalk)

lapply(listDiffMarkers_ABCs_vs_Stalk, dim)

##Add geneSymbol in column (for the export)
listDiffMarkers_ABCs_vs_Stalk<-lapply(listDiffMarkers_ABCs_vs_Stalk,function(x){x<-cbind(x,'gene'=rownames(x))})
# ##Filter on adj.P-value
# listDiffMarkers_ABCs_vs_Stalk<-lapply(listDiffMarkers_ABCs_vs_Stalk, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDiffMarkers_ABCs_vs_Stalk<-lapply(listDiffMarkers_ABCs_vs_Stalk, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                                                           x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
# listDiffMarkers_ABCs_vs_Stalk<-lapply(listDiffMarkers_ABCs_vs_Stalk, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
listDiffMarkers_ABCs_vs_Stalk<-lapply(listDiffMarkers_ABCs_vs_Stalk,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDiffMarkers_ABCs_vs_Stalk,file=paste0(sampleFolder,"results/Robjects/Markers_ABCs_vs_Stalk_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDiffMarkers_ABCs_vs_Stalk, paste0(sampleFolder,"results/Marker_lists/Markers_ABCs_vs_Stalk_",sampleName,".xlsx"))

##################################################################################

## Perform GO enrichment on 82w vs 7w lists (Stalk and Stromal)
listDiffMarkers_volcano_stalk <- readRDS(file=paste0(sampleFolder,"results/Robjects/Markers_Stalk_FBs_volcanoplot_",sampleName,".rds"))
listDiffMarkers_volcano_stromal <- readRDS(file=paste0(sampleFolder,"results/Robjects/Markers_Stromal_FBs_volcanoplot_",sampleName,".rds"))

## Filter lists
# Stalk
listDiffMarkers_volcano_stalk_DE<-lapply(listDiffMarkers_volcano_stalk, function(x){dplyr::filter(x, p_val_adj<0.01)})
listDiffMarkers_volcano_stalk_DE<-lapply(listDiffMarkers_volcano_stalk_DE, function(x){rbind(x[x$avg_logFC > 0.25,],
                                                                                       x[x$avg_logFC < -0.25,])})
# Stromal
listDiffMarkers_volcano_stromal_DE<-lapply(listDiffMarkers_volcano_stromal, function(x){dplyr::filter(x, p_val_adj<0.01)})
listDiffMarkers_volcano_stromal_DE<-lapply(listDiffMarkers_volcano_stromal_DE, function(x){rbind(x[x$avg_logFC > 0.25,],
                                                                                             x[x$avg_logFC < -0.25,])})
## GSEA in R:
#Load packages
library(clusterProfiler)
library(org.Mm.eg.db)
# options(connectionObserver = NULL) #If issue loading package org.MM.eg.db, need to run this line!!!

#Create dir
dir.create(paste0(sampleFolder,"results/GSEA"))

## 3 comparisons Stalk FBs
EnrichGO_Stalk<-tibble::lst()
for (Comparison in 1:length(listDiffMarkers_volcano_stalk_DE)) {
  #Background:
  Background<-as.character(rownames(seuratObjNew)) #Each list use same background!
  #DE genes:
  
  #EnrichGO:
  EnrichGO_Stalk[Comparison]<-enrichGO( #No results for last comparison!!
    as.character(rownames(listDiffMarkers_volcano_stalk_DE[[Comparison]])),
    'org.Mm.eg.db',
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = Background,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = T
  )
}

## Save results
saveRDS(EnrichGO_Stalk, paste0(sampleFolder,"results/GSEA/EnrichGO_Stalk_FBs_results_",sampleName,".rds"))

for (Comparison in 1:length(EnrichGO_Stalk)) { #1:length(listDiffMarkers) -> only 3 because no results for  comp 4!!
  ## Create dotplot top 10 each category
  D_background<-dotplot(EnrichGO_Stalk[[Comparison]], split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
  
  ## Check genes associated with GO terms
  GeneList_Test_background<-listDiffMarkers_volcano_stalk_DE[[Comparison]]$avg_logFC
  names(GeneList_Test_background)<-listDiffMarkers_volcano_stalk_DE[[Comparison]]$geneSymbol
  H_background<-heatplot(EnrichGO_Stalk[[Comparison]], foldChange=GeneList_Test_background)
  
  ## Save files
  pdf(file=paste0(sampleFolder,"results/GSEA/Dotplot_",names(listDiffMarkers_volcano_stalk_DE)[Comparison],"_",sampleName,".pdf"), width = 15, height = 10)
  print(D_background)
  dev.off()
  
  
  pdf(file=paste0(sampleFolder,"results/GSEA/Heatplot_",names(listDiffMarkers_volcano_stalk_DE)[Comparison],"_",sampleName,".pdf"), width = 40, height = 10)
  print(H_background)
  dev.off()
  
  write.xlsx(EnrichGO_Stalk[[Comparison]]@result,file=paste0(sampleFolder,"results/GSEA/EnrichGO_",names(listDiffMarkers_volcano_stalk_DE)[Comparison],"_",sampleName,".xlsx"))
}

## 3 comparisons Stromal FBs
EnrichGO_Stromal<-tibble::lst()
for (Comparison in 1:length(listDiffMarkers_volcano_stromal_DE)) {
  #Background:
  Background<-as.character(rownames(seuratObjNew)) #Each list use same background!
  #EnrichGO:
  EnrichGO_Stromal[Comparison]<-enrichGO( #No results for last comparison!!
    as.character(rownames(listDiffMarkers_volcano_stromal_DE[[Comparison]])),
    'org.Mm.eg.db',
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = Background,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = T
  )
}

## Save results
saveRDS(EnrichGO_Stromal, paste0(sampleFolder,"results/GSEA/EnrichGO_Stromal_FBs_results_",sampleName,".rds"))

for (Comparison in 1:length(EnrichGO_Stromal)) { #1:length(listDiffMarkers) -> only 3 because no results for  comp 4!!
  ## Create dotplot top 10 each category
  D_background<-dotplot(EnrichGO_Stromal[[Comparison]], split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
  
  ## Check genes associated with GO terms
  GeneList_Test_background<-listDiffMarkers_volcano_stromal_DE[[Comparison]]$avg_logFC
  names(GeneList_Test_background)<-listDiffMarkers_volcano_stromal_DE[[Comparison]]$geneSymbol
  H_background<-heatplot(EnrichGO_Stromal[[Comparison]], foldChange=GeneList_Test_background)
  
  ## Save files
  pdf(file=paste0(sampleFolder,"results/GSEA/Dotplot_",names(listDiffMarkers_volcano_stromal_DE)[Comparison],"_",sampleName,".pdf"), width = 15, height = 10)
  print(D_background)
  dev.off()
  
  
  pdf(file=paste0(sampleFolder,"results/GSEA/Heatplot_",names(listDiffMarkers_volcano_stromal_DE)[Comparison],"_",sampleName,".pdf"), width = 40, height = 10)
  print(H_background)
  dev.off()
  
  write.xlsx(EnrichGO_Stromal[[Comparison]]@result,file=paste0(sampleFolder,"results/GSEA/EnrichGO_",names(listDiffMarkers_volcano_stromal_DE)[Comparison],"_",sampleName,".xlsx"))
}

#################################################################################

## FeaturePlot claudins and occludins (11/04/22)
Features <- rownames(seuratObjNew)[c(grep("Cldn",rownames(seuratObjNew)),grep("Ocln",rownames(seuratObjNew)))]

# PER 36
plots1 <- FeaturePlot(seuratObjNew, features = Features[1:13],order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots1 <- lapply(X = plots1, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))

pdf(file = paste0(sampleFolder,"results/Feature_plots/Claudins_Featureplots_q2-q98_per36_",sampleName,".pdf"), width = 15, height = 10)
CombinePlots(plots = plots1, ncol=5)
dev.off()


#################################################################################

## FeaturePlot Ccl2 (02/05/22)
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObjNew, group.by = "newClusters", features = "Ccl2", 
            cols = Colors_dotplot) + RotatedAxis()

F1<-FeaturePlot(object = seuratObjNew, features = "Ccl2", cols = c("yellow", "red"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', label = T, pt.size = 1.5, order=T) 

pdf(file=paste0(sampleFolder,"results/Feature_plots/Feature_plot_Ccl2_",sampleName,".pdf"), height = 20, width = 15)
F1/D1
dev.off()

# ################################################################################
# ########## EXPORT FOR TSNE TOOL
# ################################################################################
# dir.create(paste0(sampleFolder,"results/neededFilesOnlineTool"))
# 
# listLabels<-list('ABR1','ABR2','ABR3','DV1')
# 
# exportShiny<-DisneyTools::expShiny(seuratObj, conditions = unlist(listLabels), assay = "RNA", clustering.labels = "sliced_clusters") #Use annotated version!
# saveRDS(exportShiny, file=paste0(sampleFolder,"results/neededFilesOnlineTool/SO_",sampleName,".rds"))
# annotation <- exportShiny@meta.data$annotated_clusters %>% as.factor()
# saveRDS(annotation, file=paste0(sampleFolder,"results/neededFilesOnlineTool/SO_",sampleName,"_annotation.rds"))
# 
# # exportShinyNew<-DisneyTools::expShiny(seuratObjNew, conditions = unlist(listLabels), assay = "RNA")
# # saveRDS(exportShinyNew, file=paste0(sampleFolder,"results/neededFilesOnlineTool/SO_",sampleName,"_Pos_vs_Neg.rds"))
# # annotationNew <- exportShinyNew@meta.data$newClusters %>% as.factor() 
# # saveRDS(annotationNew, file=paste0(sampleFolder,"results/neededFilesOnlineTool/SO_",sampleName,"_Pos_vs_Neg_annotation.rds"))
# 
