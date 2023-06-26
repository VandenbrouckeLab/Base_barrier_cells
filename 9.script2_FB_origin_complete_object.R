## Script for merging our datasets with public datasets to investigate the origin of our FBs
## Follow-up script to process and explore the Fibroblast origin complete object in Figure 1 manuscript
## Initially the object still contains everything and subsequently it is subsetted to remove ChP Epithelial cells and Immune cells
## Rebuttal: stronger harmony settings

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")

sampleName <- "Rebuttal2_Full_merge_FB_datasets" 
sampleFolder<-paste0("Rebuttal_mouse_data/Full_merge","/")

##add some subfolders
dir.create(paste0(sampleFolder,"results_merge"))
dir.create(paste0(sampleFolder,"results_merge/QC"))
dir.create(paste0(sampleFolder,"results_merge/Robjects"))

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
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

## Check 1.0 clustering: stick to 0.8
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "RNA_snn_res.1", label.size = 8)

# seuratObj@meta.data$harmony_clusters<-seuratObj$RNA_snn_res.1
Idents(seuratObj)<-seuratObj@meta.data$harmony_clusters
seuratObj@meta.data$sliced_clusters<-seuratObj@meta.data$harmony_clusters

DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

dir.create(paste0(sampleFolder,"results_merge/Annotation"))

pdf(file=paste0(sampleFolder,"results_merge/Annotation/1_annotation_numbered_harmony_",sampleName,".pdf"), width = 14, height = 10)
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "harmony_clusters", label.size = 4, repel = T)
dev.off()

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
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'nCount_RNA',"UMI")

ggsave(p2, file=paste0(sampleFolder,"results_merge/QC/1_UMI_",sampleName,".png"), width = 20)

########## mito.genes plot ##########
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'subsets_Mito_percent',"mito")

ggsave( p2, file=paste0(sampleFolder,"results_merge/QC/2_percMito_",sampleName,".png"), width = 20)


########## PCA plot ##########
# pdf(file=paste0(sampleFolder,"results_merge/QC/13a_PCA_",sampleName,".pdf"), width=10)
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,3))
# dev.off()

#RNA clusters
pdf(file=paste0(sampleFolder,"results_merge/Annotation/1_annotation_Color_RNA_clusters_on_harmony_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$harmony_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$harmony_clusters==i),])))
  C1<-C1+ggtitle(paste0("Harmony_cluster_",i))
  print(C1)
}
dev.off()

################################################################################
########## CHECK DE GENES
################################################################################
dir.create(paste0(sampleFolder,"results_merge/Feature_plots"))

##### Epithelial marker ->
Features<-c("Otx2", "Ttr")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"CPE"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Endothelial marker ->
Features<-c("Pecam1", "Flt1","Plvap")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"EC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Vascular associated marker ->
Features<-c("Pdgfrb", "Mylk","Myh11","Tagln")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"VAC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Fibroblast marker ->
Features<-c("Dcn", "Col1a1")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"FB"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Neuronal + glial cell marker ->
Features<-c("Tubb3","Slc1a3", "Fabp7","Olig1") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Neuronal_and_glial"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Mitotic cell marker ->
Features<-c("Birc5", "Mki67")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Mitotic_cells"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

## Shah dataset: Ependymal cells
Features<-c("Meis2","Meg3", "Stmn2","Btg1","Celf4") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Shah_Neuroblasts"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Hmgb2","Hist1h2ap", "Top2a","Cenpf") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Shah_aNSCs_late"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Ccdc153", "Mia", "Rarres2","Tmem212","Tm4sf1") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Shah_Ependymal_cells"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Htra1", "Atp1a2", "Ntm","Aldoc","Mfge8") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Shah_qNSCs"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Cst3", "Fxyd1", "Id3","Fbxo2","Gfap") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Shah_aNSCs"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Plp1", "Olig1", "Mbp","Ptgds","Cldn11") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Shah_OL_OPCs"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Bsg", "Cldn5", "Itm2a","Ly6c1","Cxcl12") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Shah_Undefined1"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

#########################################################################################################################################

########################################
##### all clusters vs all clusters
########################################

# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

dir.create(paste0(sampleFolder,"results_merge/Marker_lists"))

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
RNAMarkers_RNAclus <- FindAllMarkers(seuratObj, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"results_merge/Robjects/RNAmarkersList_RNAclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerRNAcluster']]<-paste0(table(RNAMarkers_RNAclus$cluster)," RNA markers for RNA cluster ",rownames(table(RNAMarkers_RNAclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
names(RNAmarkersList_RNAclus)<-paste0("RNAcluster",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results_merge/Marker_lists/RNAmarkersList_RNAclus_",sampleName,".xlsx"))

########################################################################################################################

## Check other annotation
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "ClusterName", label.size = 4)

DimPlot(seuratObj, reduction = "umap", label = T, group.by = "Clusters", label.size = 4)

pdf(file=paste0(sampleFolder,"results_merge/Annotation/1_UMAP_split_datasets_",sampleName,".pdf"), width = 40, height = 30)
DimPlot(seuratObj, reduction = "umap", label = T, split.by = "orig.ident", ncol = 4)
dev.off()

Idents(seuratObj)<- seuratObj@meta.data$orig.ident
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "1_Foxcl_ctl"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "2_Foxcl_ctl"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "3_Foxcl_ctl"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "GSM_2677817"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "GSM_2677818"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "GSM_2677819"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "GSE98816_Vanlandewijck"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "Vascular_cells_atlas"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "Ependymal_cells_atlas"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "RVD1_LpsNegFour"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "RVD2_LpsNegLat"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "RVD5_Y4V"))
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "RVD6_YLV"))

pdf(file=paste0(sampleFolder,"results_merge/Annotation/1_annotation_Color_datasets_",sampleName,".pdf"), width = 15)
for (i in levels(Idents(seuratObj))) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = i))
  C1<-C1+ggtitle(i)
  print(C1)
}
dev.off()

seuratObj@meta.data$New_clusters <- "Unknown"

## Update clusters URVB
metaData_FBs <- readRDS(file = "Urvb_datasets_rebuttal/Robjects/MetaData_FBs.rds")
rownames(metaData_FBs)<-gsub("_","-",rownames(metaData_FBs))
URVB_cells<-intersect(paste0("URVB_",rownames(metaData_FBs)),rownames(seuratObj@meta.data))
URVB_cells_old<-gsub("URVB_","",URVB_cells)
seuratObj@meta.data[URVB_cells,"New_clusters"]<-as.character(metaData_FBs[URVB_cells_old,"annotated_clusters"])

## Update clusters Vasc atlas
VASC_atlas_cells<-colnames(seuratObj[,grep("VASC2",colnames(seuratObj))])

seuratObj@meta.data[VASC_atlas_cells,"New_clusters"]<-seuratObj@meta.data[VASC_atlas_cells,"clusters3..duplicated.cells3.."]

## Update clusters Ep atlas
EP_atlas_cells<-colnames(seuratObj[,grep("EP4",colnames(seuratObj))])

seuratObj@meta.data[EP_atlas_cells,"New_clusters"]<-seuratObj@meta.data[EP_atlas_cells,"ClusterName"]

## Update clusters Fb Desisto
FB_Desisto_cells<-colnames(seuratObj[,grep("FB1",colnames(seuratObj))])

seuratObj@meta.data[FB_Desisto_cells,"New_clusters"]<-seuratObj@meta.data[FB_Desisto_cells,"Clusters"]

## Update clusters Vasc Vanlandewijck
Vasc_Vldw_cells<-colnames(seuratObj[,grep("VASC1",colnames(seuratObj))])
Vasc_annot<-gsub("VASC1_","",Vasc_Vldw_cells)
Mat_Vasc_annot<-matrix(unlist(strsplit(Vasc_annot, "_")), ncol=2, byrow=TRUE)
Vasc_annot_v2<-Mat_Vasc_annot[,1]

seuratObj@meta.data[Vasc_Vldw_cells,"New_clusters"]<-Vasc_annot_v2

## Save plot
pdf(file=paste0(sampleFolder,"results_merge/Annotation/1_annotation_original_dataset_harmony_",sampleName,".pdf"), width = 17, height = 10)
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "New_clusters", label.size = 2, repel = T)
dev.off()

########################################################################################################################

########################################
##### RNA clusters post-slice 
########################################
########################################

seuratObj@meta.data$sliced_clusters<- factor(seuratObj@meta.data$sliced_clusters,sort(as.numeric(levels(seuratObj@meta.data$sliced_clusters)))) #reorder levels
seuratObj@meta.data$annotated_clusters <- factor(seuratObj@meta.data$annotated_clusters,sort(as.numeric(levels(seuratObj@meta.data$annotated_clusters)))) #reorder levels

########################################
##### Manual annotation
########################################
levels(seuratObj@meta.data$annotated_clusters) <- c("Endothelial and vascular cells","Fibroblasts","Endothelial and vascular cells",
                                                    "Fibroblasts","Proliferating cells","Fibroblasts","Smooth muscle cells",
                                                    "Endothelial and vascular cells","Fibroblasts","Endothelial and vascular cells",
                                                    "Neuronal and glial cells","Fibroblasts","Endothelial and vascular cells",
                                                    "Neuronal and glial cells","Neuronal and glial cells","Endothelial and vascular cells",
                                                    "Neuronal and glial cells","Endothelial and vascular cells",
                                                    "Neuronal and glial cells","Perivascular macrophages")
U_annot<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_merge/Annotation/2_UMAP_annotated_",sampleName,".png"), height = 10, width = 16, dpi = "retina")

Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

RNAMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_SCTclus$cluster)
saveRDS(RNAMarkers_SCTclus, file=paste0(sampleFolder,"results_merge/Robjects/RNAmarkersList_SCTclus_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerSCTclusterannotated']]<-paste0(table(RNAMarkers_SCTclus$cluster)," RNA markers for annotated cluster ",rownames(table(RNAMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrRNAclusters_SCTclus<-names(table(RNAMarkers_SCTclus$cluster))
# totalNrRNAclusters_SCTclusPlusOne<-totalNrRNAclusters_SCTclus
RNAmarkersList_SCTclus<-list()

for(i in totalNrRNAclusters_SCTclus){
  # clusterNr<-i-1
  
  tmp<-RNAMarkers_SCTclus[RNAMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  RNAmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
# names(RNAmarkersList_SCTclus)<-paste0("SCTclustersliced",0:totalNrRNAclusters_SCTclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge/Marker_lists/RNAmarkersList_RNAclus_",sampleName,"_annotated.xlsx"))


######################################################

# Frequency tables (sliced)
Sample <- seuratObj@meta.data$orig.ident
cluster <- seuratObj@meta.data$annotated_clusters
Aggr <- rep(sampleName,length(cluster))

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

# Stacked
library("ggthemes")

png(file=paste0(sampleFolder,"results_merge/Annotation/3_SampleDistribution_ggplot2_annot_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"results_merge/Annotation/3_SampleDistribution_ggplot2_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png(file=paste0(sampleFolder,"results_merge/Annotation/3_SampleDistribution_ggplot2_annot_3_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data2, aes(fill=cluster, y=Freq, x=Aggr)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

######################################################

## Extra metadata

seuratObj@meta.data$sample_origin2<-as.factor(seuratObj@meta.data$orig.ident)
levels(seuratObj@meta.data$sample_origin2)<-c(rep("Meningeal fibroblasts (DeSisto et al)",3),
                                              "Ependymal cells (Zeisel et al)","Brain vascular cells (Vanlandewijck et al)",
                                              rep("Ependymal cells (Shah et al)",3),
                                              rep("ChP vascular cells (Verhaege et al)",4),
                                              "CNS vascular cells (Zeisel et al)")

######################################################

# Frequency tables (sliced)
Sample <- seuratObj@meta.data$sample_origin2
cluster <- seuratObj@meta.data$annotated_clusters
Aggr <- rep(sampleName,length(cluster))

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(Sample,Aggr))

# Stacked
library("ggthemes")
library(RColorBrewer)
Colorset<-c(brewer.pal(8,"Set1")[c(1,7,2,3,4,5)])

p1<-ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white") + scale_fill_manual(values=Colorset)

p2<-ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p3<-ggplot(data2, aes(fill=Sample, y=Freq, x=Aggr)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white") + scale_fill_manual(values=Colorset)

pdf(file=paste0(sampleFolder,"results_merge/Annotation/3_SampleDistribution_ggplot2_annot_4_",sampleName,".pdf"), width = 10, height = 7.5)
print(p1)
dev.off()

pdf(file=paste0(sampleFolder,"results_merge/Annotation/3_SampleDistribution_ggplot2_annot_5_",sampleName,".pdf"), width = 10, height = 7.5)
print(p2)
dev.off()

pdf(file=paste0(sampleFolder,"results_merge/Annotation/3_SampleDistribution_ggplot2_annot_6_",sampleName,".pdf"), width = 10, height = 7.5)
print(p3)
dev.off()

###########################

U_annot_sub<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)

pdf(file = paste0(sampleFolder,"results_merge/Annotation/4_UMAP_annotated_",sampleName,".pdf"), width = 16, height = 10)
U_annot_sub
dev.off()

U_origin<-DimPlot(seuratObj, reduction = "umap", label = F, group.by = "sample_origin2", label.size = 3, cols = Colorset)
pdf(file = paste0(sampleFolder,"results_merge/Annotation/4_UMAP_origin_",sampleName,".pdf"), width = 14, height = 10)
U_origin
dev.off()

##############################################

## ScType Automatic annotation 
## https://github.com/IanevskiAleksandr/sc-type

# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

## Perform with own custom list
gs_list_custom<-readRDS("Rebuttal_mouse_data/Automatic_coarse_annotation_file_rebuttal_scType_filtered.rds") #Filtered list Panglao

# get cell-type by cell matrix
## Adapted to include all genes instead of just the 2000 scaled genes -> better scores!!
es.max = sctype_score(scRNAseqData = seuratObj[["RNA"]]@data, scaled = F, 
                      gs = gs_list_custom, gs2 = NULL, gene_names_to_uppercase = F)

# Please note that sctype_score function (used above) accepts both positive and negative markers through gs and gs2 arguments. 
# In case, there are no negative markers (i.e. markers providing evidence against a cell being of specific cell type) 
# just set gs2 argument to NULL (i.e. gs2 = NULL).

# merge by cluster
cL_results = do.call("rbind", lapply(unique(seuratObj@meta.data$harmony_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seuratObj@meta.data[seuratObj@meta.data$harmony_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seuratObj@meta.data$harmony_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

seuratObj@meta.data$customclassif2 = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seuratObj@meta.data$customclassif2[seuratObj@meta.data$harmony_clusters == j] = as.character(cl_type$type[1])
}

## Save results
Colorset<-c(brewer.pal(12,"Set3")[c(1,3:6,7,8,10)],"maroon2",brewer.pal(12,"Set3")[c(11,12)])

U_automatic<-DimPlot(seuratObj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif2', label.size = 3, cols = Colorset)

pdf(file = paste0(sampleFolder,"results_merge/Annotation/4_UMAP_automatic_filtered_",sampleName,".pdf"), width = 14, height = 10)
U_automatic
dev.off()

write.xlsx(cL_results, file =paste0(sampleFolder, "results_merge/Marker_lists/Automatic_annotation_scores_filtered_",sampleName,".xlsx"))

##############################################

## Extra: plots for rebuttal for quality data and M&M

## Convert to sce
sce<-as.SingleCellExperiment(seuratObj)

##### Get mitochondrial genes #####
is.mito <- grepl("^mt-", rownames(sce), ignore.case = T)
sum(is.mito)
##13
rownames(sce)[is.mito]

##### Calculate QC metrics #####
### => pData(sce) is created
sce <- perCellQCMetrics(sce, subsets=list(Mt=rownames(sce)[is.mito]))

##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("orig.ident"=seuratObj$sample_origin2,
                     "nGene"=sce$detected,"nUMI"=sce$sum,"percent.mito"=sce$subsets_Mt_percent, 
                     stringsAsFactors = F)

pdf(file=paste0(sampleFolder,"results_merge/Annotation/5_QC_stats_rebuttal_",sampleName,".pdf"))
for (i in levels(seuratObj$sample_origin2)[c(1,2,4,5,6)]){
  metaData_subset<-metaData[which(metaData$orig.ident==i),]
  p1<-ggplot(metaData_subset, aes(x = nUMI)) +
    geom_histogram(binwidth=100) +
    xlab(paste0("nUMI ",i)) +
    ylab("Frequency") +
    theme_bw() #Count depth
  
  p2<-ggplot(metaData_subset, aes(x = nGene)) +
    geom_histogram(binwidth=20) +
    xlab(paste0("nGene ",i)) +
    ylab("Frequency") +
    theme_bw() # histogram of nr of genes
  
  p3<-ggplot(metaData_subset, aes(x = percent.mito)) +
    geom_histogram(binwidth=1) +
    xlab(paste0("percent.mito ",i)) +
    ylab("Frequency") +
    theme_bw()
  
  p4<-ggplot(metaData_subset, aes(x = nUMI,y=nGene,colour=percent.mito)) +
    geom_point(size=0.5) +
    scale_color_gradient2(midpoint=10, low="black", mid="white",
                          high="red", space ="Lab" )+
    xlab(paste0("nUMI ",i)) +
    ylab(paste0("Nr of Genes ",i)) +
    geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
    theme_bw()
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
}
dev.off()

pdf(file=paste0(sampleFolder,"results_merge/Annotation/5_QC_stats_rebuttal_extra_",sampleName,".pdf"))
for (i in "Brain vascular cells (Vanlandewijck et al)"){
  metaData_subset<-metaData[which(metaData$orig.ident==i),]
  p1<-ggplot(metaData_subset, aes(x = nUMI)) +
    geom_histogram(binwidth=10000) +
    xlab(paste0("nUMI ",i)) +
    ylab("Frequency") +
    theme_bw() #Count depth
  
  p2<-ggplot(metaData_subset, aes(x = nGene)) +
    geom_histogram(binwidth=20) +
    xlab(paste0("nGene ",i)) +
    ylab("Frequency") +
    theme_bw() # histogram of nr of genes
  
  p4<-ggplot(metaData_subset, aes(x = nUMI,y=nGene)) +
    geom_point(size=0.5) +
    # scale_color_gradient2(midpoint=10, low="black", mid="white",
    #                       high="red", space ="Lab" )+
    xlab(paste0("nUMI ",i)) +
    ylab(paste0("Nr of Genes ",i)) +
    geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
    theme_bw()
  print(p1)
  print(p2)
  print(p4)
  
}
dev.off()

##############################################
## Create diet object for online tool
seuratObj_diet<-DietSeurat(seuratObj, counts = T, data = T, scale.data = F,
                           assays = c("RNA"), dimreducs = "umap", graphs = NULL)

## New metadata names
seuratObj_diet$Numbered_clusters<-seuratObj_diet$harmony_clusters
seuratObj_diet$Automatic_annotation<-seuratObj_diet$customclassif2
seuratObj_diet$Origin<-seuratObj_diet$sample_origin2
seuratObj_diet$Manual_annotation<-seuratObj_diet$annotated_clusters

## Metadata columns
# "orig.ident" "Phase" "SCT_clusters" "ADT_clusters" "annotated_clusters" "final_annotation2021"  
DimPlot(seuratObj_diet, reduction = "umap", label = T, repel = F, group.by = "Automatic_annotation", label.size = 3,
        cols =Colorset2)
# cols =c(brewer.pal(n = 6, name = "Set3")[1],"Yellow",brewer.pal(n = 6, name = "Set3")[c(3:6)]))
library(RColorBrewer)
Colorset1<-c(brewer.pal(8,"Set1")[c(1,7,2,3,4,5)]) #origin
Colorset2<-c(brewer.pal(12,"Set3")[c(1,3:6,7,8,10)],"#EE30A7",brewer.pal(12,"Set3")[c(11,12)]) #auto

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

All<-c(levels(as.factor(seuratObj_diet$Manual_annotation)),levels(as.factor(seuratObj_diet$Origin)),
       levels(as.factor(seuratObj_diet$Numbered_clusters)),levels(as.factor(seuratObj_diet$Automatic_annotation)))
Color_info<-c(gg_color_hue(6),Colorset1,gg_color_hue(20),Colorset2)
Metadata_column<-c(rep("Manual_annotation",6),rep("Origin",6),
                   rep("Numbered_clusters",20),rep("Automatic_annotation",11))
Info_Kevin<-as.data.frame(cbind(All,Color_info,Metadata_column))

write.xlsx(Info_Kevin, file =paste0(sampleFolder, "results_merge/Annotation/Info_Kevin_rebuttal_",sampleName,".xlsx"))

##### Save object
saveRDS(seuratObj_diet, file=paste0(sampleFolder,"Robjects/seuratObj_paper_diet_rebuttal_",sampleName,"_2023.rds"))

##### Read object
seuratObj_diet<-readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_paper_diet_rebuttal_",sampleName,"_2023.rds"))
