## Creating Fibroblast origin subset object with only the Fibroblast clusters from the Fibroblast origin complete object
## Follow-up script to process and explore the Fibroblast origin subset object in Figure 1 manuscript
## Rebuttal: Stronger harmony settings

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('openxlsx')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")

sampleName <- "Rebuttal2_Full_merge_subset_FB_datasets" 
sampleFolder<-paste0("Rebuttal_mouse_data/Merge_subset","/")

##add some subfolders
dir.create(paste0(sampleFolder,"results_merge_subset"))
dir.create(paste0(sampleFolder,"results_merge_subset/QC"))
dir.create(paste0(sampleFolder,"results_merge_subset/Robjects"))

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

Idents(seuratObj)<-seuratObj@meta.data$harmony_clusters_subset
seuratObj@meta.data$sliced_clusters<-seuratObj@meta.data$harmony_clusters_subset

## Check 1.0 clustering:
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "RNA_snn_res.1", label.size = 8)
## Better with 0.8  

## Save new clustering
seuratObj@meta.data$annotated_clusters_subset <- seuratObj@active.ident #Sliced clustering

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

ggsave(p2, file=paste0(sampleFolder,"results_merge_subset/QC/1_UMI_",sampleName,".png"), width = 10)

########## mito.genes plot ##########
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'subsets_Mito_percent',"mito")

ggsave( p2, file=paste0(sampleFolder,"results_merge_subset/QC/2_percMito_",sampleName,".png"), width = 10)


########## PCA plot ##########
# pdf(file=paste0(sampleFolder,"results_merge_subset/QC/13a_PCA_",sampleName,".pdf"), width=10)
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,3))
# dev.off()

#RNA clusters
dir.create(paste0(sampleFolder,"results_merge_subset/Annotation"))

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/1_annotation_Color_RNA_clusters_on_harmony_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$harmony_clusters_subset))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$harmony_clusters_subset==i),])))
  C1<-C1+ggtitle(paste0("Harmony_cluster_",i))
  print(C1)
}
dev.off()

################################################################################
########## CHECK DE GENES
################################################################################
dir.create(paste0(sampleFolder,"results_merge_subset/Feature_plots"))

##### Epithelial marker ->
Features<-c("Otx2", "Ttr")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"CPE"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Endothelial marker ->
Features<-c("Pecam1", "Flt1","Plvap")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"EC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Vascular associated marker ->
Features<-c("Pdgfrb", "Mylk","Myh11","Tagln")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"VAC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Fibroblast marker ->
Features<-c("Dcn", "Col1a1","Dpep1","Igfbp6")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"FB"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Pdgfra", "Pdgfrb","Lum","Acta2","Rgs5")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"FB_bis"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Neuronal + glial cell marker ->
Features<-c("Tubb3","Slc1a3", "Fabp7","Olig1") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Neuronal_and_glial"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Ependymal cell marker ->
Features<-c("Foxj1","Tmem212","Ccdc153","Dnah11") #Tubb3 neuronal
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Ependymal"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Mitotic cell marker ->
Features<-c("Birc5", "Mki67")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Mitotic_cells"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_subset/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

#########################################################################################################################################

########################################
##### all clusters vs all clusters
########################################

# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

dir.create(paste0(sampleFolder,"results_merge_subset/Marker_lists"))

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
RNAMarkers_RNAclus <- FindAllMarkers(seuratObj, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"results_merge_subset/Robjects/RNAmarkersList_RNAclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerRNAcluster']]<-paste0(table(RNAMarkers_RNAclus$cluster)," RNA markers for RNA cluster ",rownames(table(RNAMarkers_RNAclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results_merge_subset/Marker_lists/RNAmarkersList_RNAclus_",sampleName,".xlsx"))

########################################################################################################################

## Check new clustering
## Save plot
pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/1_annotation_new_numbered_",sampleName,".pdf"), width = 15, height = 10)
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "harmony_clusters_subset", label.size = 4, repel = T)
dev.off()

## Check original annotation
## Save plot
pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/1_annotation_original_datasets_",sampleName,".pdf"), width = 25, height = 15)
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "New_clusters", label.size = 2, repel = T)
dev.off()

########################################################################################################################

## Check other annotation

Idents(seuratObj)<- seuratObj@meta.data$orig.ident

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/1_annotation_Color_datasets_",sampleName,".pdf"), width = 15)
for (i in levels(Idents(seuratObj))) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = i))
  C1<-C1+ggtitle(i)
  print(C1)
}
dev.off()

Idents(seuratObj)<-seuratObj@meta.data$New_clusters #Revert

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/1_annotation_Color_old_idents_",sampleName,".pdf"), width = 15)
for (i in levels(Idents(seuratObj))) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = i))
  C1<-C1+ggtitle(i)
  print(C1)
}
dev.off()

Idents(seuratObj)<-seuratObj@meta.data$harmony_clusters_subset #Revert

DimPlot(seuratObj, reduction = "umap", label = F, repel = T, group.by = "New_clusters", label.size = 3)

# Check ependymal datasets cells
colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Unknown"),]))

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
FeaturePlot(object = seuratObj, features = "Ntm", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Tagln", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Wnt4", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Stmn1", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Ptgds", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Ptgds", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Mbp", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Spp1", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Col15a1", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "eGFP", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Pla1a", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)
FeaturePlot(object = seuratObj, features = "Tm4sf1", cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)

########################################################################################################################

########################################
##### RNA clusters post-slice
########################################
########################################
seuratObj@meta.data$sliced_clusters<- factor(seuratObj@meta.data$sliced_clusters,sort(as.numeric(levels(seuratObj@meta.data$sliced_clusters)))) #reorder levels
seuratObj@meta.data$annotated_clusters_subset <- factor(seuratObj@meta.data$annotated_clusters_subset,sort(as.numeric(levels(seuratObj@meta.data$annotated_clusters_subset)))) #reorder levels

########################################
##### Markers annotated clusters
########################################
levels(seuratObj@meta.data$annotated_clusters_subset) <- c("Ntm hi FBs 1","Ntm hi FBs 2", #Also Col18a1!
                                                           "Tagln hi Arachnoid FBs","Wnt4 hi FBs",
                                                           "FBs 1","VLMCs","Col15a1 hi FBs","Barrier cells",
                                                           "EC/Mural like FBs","FBs 2",
                                                           "Contaminating endothelial cells", "Contaminating ependymal cells")

U_annot<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters_subset", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_merge_subset/Annotation/2_UMAP_annotated1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters_subset

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

names(RNAmarkersList_RNAclus)<-levels(as.factor(seuratObj$annotated_clusters_subset))

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results_merge_subset/Marker_lists/RNAmarkersList_RNAclus_",sampleName,"_annotated.xlsx"))

######################################################

# Frequency tables (sliced)
Sample <- seuratObj@meta.data$sample_origin2
cluster <- seuratObj@meta.data$annotated_clusters_subset
Aggr <- rep(sampleName,length(cluster))

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

# Stacked
library("ggthemes")

png(file=paste0(sampleFolder,"results_merge_subset/Annotation/3_SampleDistribution_ggplot2_annot_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"results_merge_subset/Annotation/3_SampleDistribution_ggplot2_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png(file=paste0(sampleFolder,"results_merge_subset/Annotation/3_SampleDistribution_ggplot2_annot_3_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data2, aes(fill=cluster, y=Freq, x=Aggr)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")
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

# Next, let's prepare gene sets from the input cell marker file. 
# By default, we use our in-built cell marker DB, however, feel free to use your own data. 
# Just prepare an input XLSX file in the same format as our DB file. 
# DB file should contain four columns (tissueType - tissue type, cellName - cell type, 
# geneSymbolmore1 - positive marker genes, geneSymbolmore2 - marker genes not expected to be expressed by a cell type)
# In addition, provide a tissue type your data belongs to:

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# Use own custom list
gs_list_custom<-readRDS("Rebuttal_mouse_data/Automatic_coarse_annotation_file_rebuttal_scType_filtered.rds") #Filtered list Panglao -> Auto_annotation_coarse2

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = seuratObj[["RNA"]]@data, scaled = F,
                      gs = gs_list_custom, gs2 = NULL, gene_names_to_uppercase = F)

# Please note that sctype_score function (used above) accepts both positive and negative markers through gs and gs2 arguments. 
# In case, there are no negative markers (i.e. markers providing evidence against a cell being of specific cell type) 
# just set gs2 argument to NULL (i.e. gs2 = NULL).

# merge by cluster
cL_results = do.call("rbind", lapply(unique(seuratObj@meta.data$annotated_clusters_subset), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seuratObj@meta.data[seuratObj@meta.data$annotated_clusters_subset==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seuratObj@meta.data$annotated_clusters_subset==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

seuratObj@meta.data$Auto_annotation_coarse2 = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seuratObj@meta.data$Auto_annotation_coarse2[seuratObj@meta.data$annotated_clusters_subset == j] = as.character(cl_type$type[1])
}

## Save results
Colorset<-brewer.pal(12,"Set3")[c(3:6,10,11,9,12)]

U_automatic<-DimPlot(seuratObj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Auto_annotation_coarse2', label.size = 3, cols = Colorset)

pdf(file = paste0(sampleFolder,"results_merge_subset/Annotation/4_UMAP_automatic_filtered_",sampleName,".pdf"), width = 14, height = 10)
U_automatic
dev.off()

write.xlsx(cL_results, file =paste0(sampleFolder, "results_merge_subset/Marker_lists/Automatic_coarse_annotation_scores_filtered_",sampleName,".xlsx"))

## Check complete origin object results
Colorset<-c(brewer.pal(12,"Set3")[c(5,12)])

U_automatic_old<-DimPlot(seuratObj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif2', label.size = 3, cols = Colorset)

pdf(file = paste0(sampleFolder,"results_merge_subset/Annotation/4_UMAP_automatic_filtered_full_merge_object_",sampleName,".pdf"), width = 14, height = 10)
U_automatic_old
dev.off()

########################################################################################################################
########################################################################################################################
########################################################################################################################

# Remove bad clusters!!
seuratObj_clean<-subset(seuratObj, idents =c("Ntm hi FBs 1","Ntm hi FBs 2", #Also Col18a1!
                                             "Tagln hi Arachnoid FBs","Wnt4 hi FBs",
                                             "FBs 1","VLMCs","Col15a1 hi FBs","Barrier cells",
                                             "EC/Mural like FBs","FBs 2"))

# Clean outliers!!!!
U1 <- DimPlot(seuratObj_clean, reduction = "umap", label = T, label.size = 4)
seuratObj_clean <- CellSelector(U1, object=seuratObj_clean, ident="Outliers1")

U_outlier<-DimPlot(seuratObj_clean, reduction = "umap", label = T, repel = T, label.size = 4)
ggsave(U_outlier, file=paste0(sampleFolder,"results_merge_subset/Annotation/5_UMAP_clean_outliers_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

seuratObj_clean<-subset(seuratObj_clean, idents =c("Ntm hi FBs 1","Ntm hi FBs 2", #Also Col18a1!
                                                   "Tagln hi Arachnoid FBs","Wnt4 hi FBs",
                                                   "FBs 1","VLMCs","Col15a1 hi FBs","Barrier cells",
                                                   "EC/Mural like FBs","FBs 2"))

# Remake figures
Idents(seuratObj_clean)<-as.factor(as.character(Idents(seuratObj_clean)))
seuratObj_clean@meta.data$annotated_clusters_subset<-Idents(seuratObj_clean)

U_annot<-DimPlot(seuratObj_clean, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters_subset", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_merge_subset/Annotation/5_UMAP_clean_annotated1_",sampleName,".png"), height = 10, width = 8, dpi = "retina")

pdf(file=paste0(sampleFolder, "results_merge_subset/Annotation/5_UMAP_clean_annotated1_",sampleName,".pdf"), height = 10, width = 8)
U_annot
dev.off()

### Find RNAmarkers for every RNA cluster compared to all remaining cells, report only the positive ones
# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

RNAMarkers_RNAclus <- FindAllMarkers(seuratObj_clean, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"results_merge_subset/Robjects/RNAmarkersList_RNAclus_clean_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['RNAmarkersPerRNAclustercleanannotated']]<-paste0(table(RNAMarkers_RNAclus$cluster)," RNA markers for clean annotated cluster ",rownames(table(RNAMarkers_RNAclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrRNAclusters_RNAclus<-names(table(RNAMarkers_RNAclus$cluster))
# totalNrRNAclusters_RNAclusPlusOne<-totalNrRNAclusters_RNAclus
RNAmarkersList_RNAclus<-list()

for(i in totalNrRNAclusters_RNAclus){
  # clusterNr<-i-1
  
  tmp<-RNAMarkers_RNAclus[RNAMarkers_RNAclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  RNAmarkersList_RNAclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
# names(RNAmarkersList_RNAclus)<-paste0("RNAclustersliced",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results_merge_subset/Marker_lists/RNAmarkersList_RNAclus_clean_",sampleName,"_annotated.xlsx"))

######################################################

# Frequency tables (sliced)
Sample <- seuratObj_clean@meta.data$sample_origin2
cluster <- seuratObj_clean@meta.data$annotated_clusters_subset
Aggr <- rep(sampleName,length(cluster))

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(Sample,Aggr))

# Stacked
library("ggthemes")
library(RColorBrewer)
Colorset<-c(brewer.pal(8,"Set1")[c(1,7,2,3,4,5)])

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/5_SampleDistribution_ggplot2_annot_1_",sampleName,"_clean.pdf"), width = 10, height = 7.5)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white") + scale_fill_manual(values=Colorset)
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/5_SampleDistribution_ggplot2_annot_2_",sampleName,"_clean.pdf"), width = 10, height = 7.5)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/5_SampleDistribution_ggplot2_annot_3_",sampleName,"_clean.pdf"), width = 10, height = 7.5)
ggplot(data2, aes(fill=Sample, y=Freq, x=Aggr)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white") + scale_fill_manual(values=Colorset)
dev.off()

######################################################

## Extra vis
U_origin<-DimPlot(seuratObj_clean, reduction = "umap", label = F, group.by = "sample_origin2", label.size = 3, cols = Colorset)

U_origin2<-DimPlot(seuratObj_clean, reduction = "umap", label = F, split.by ="sample_origin2", group.by = "annotated_clusters_subset", label.size = 3, ncol = 3)

U_origin3<-DimPlot(seuratObj_clean, reduction = "umap", label = F, split.by ="sample_origin2", group.by = "sample_origin2", label.size = 3, cols = Colorset, ncol = 3)

pdf(file = paste0(sampleFolder,"results_merge_subset/Annotation/5_UMAP_clean_dataset_origin1_",sampleName,".pdf"), width = 15, height = 15)
U_origin
dev.off()

pdf(file = paste0(sampleFolder,"results_merge_subset/Annotation/5_UMAP_dataset_clean_origin2_",sampleName,".pdf"), width = 15, height = 15)
U_origin2
dev.off()

pdf(file = paste0(sampleFolder,"results_merge_subset/Annotation/5_UMAP_dataset_clean_origin3_",sampleName,".pdf"), width = 15, height = 12)
U_origin3
dev.off()

## Check original annotation
## Save plot
pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/5_UMAP_clean_with_annotation_original_datasets_",sampleName,".pdf"), width = 12, height = 12)
DimPlot(seuratObj_clean, reduction = "umap", label = T, group.by = "New_clusters", label.size = 2, repel = T)
dev.off()

## Show numbered clusters on clean dataset for rebuttal
seuratObj_clean$harmony_clusters_subset<-as.factor(as.numeric(seuratObj_clean_diet$harmony_clusters_subset))
levels(seuratObj_clean$harmony_clusters_subset)<-c(0,1,2,3,4,5,6,7,8,9)

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/5_UMAP_clean_with_numbered_clustering_",sampleName,".pdf"), width = 9, height = 12)
DimPlot(seuratObj_clean, reduction = "umap", label = T, group.by = "harmony_clusters_subset", label.size = 5, repel = T)
dev.off()

## Show automatic annotation on clean dataset for rebuttal
## Save results
Colorset<-brewer.pal(12,"Set3")[c(5:6,10,11,9,12)]

U_automatic<-DimPlot(seuratObj_clean, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Auto_annotation_coarse2', label.size = 3, cols = Colorset)

pdf(file = paste0(sampleFolder,"results_merge_subset/Annotation/5_UMAP_clean_automatic_anno_filtered_",sampleName,".pdf"), width = 9, height = 12)
U_automatic
dev.off()

##############################################

## Update old metadata: 
## Clean up Fibroblast annotation of LV (no split between FBs. Split now)
clusterMatrix<-seuratObj_clean@meta.data
umapTable<-as.data.frame(seuratObj_clean[['umap']]@cell.embeddings, stringsAsFactors = F)

Extra_BBCs<-intersect(rownames(seuratObj_clean@meta.data[which(seuratObj_clean$New_clusters == "Fibroblasts"),]),rownames(umapTable[which(umapTable$UMAP_2 > 5),]))
Extra_Stromal<-intersect(rownames(seuratObj_clean@meta.data[which(seuratObj_clean$New_clusters == "Fibroblasts"),]),rownames(umapTable[which(umapTable$UMAP_2 < 5),]))
seuratObj_clean@meta.data$New_clusters_clean<-as.character(seuratObj_clean@meta.data$New_clusters)
seuratObj_clean@meta.data[Extra_BBCs,"New_clusters_clean"]<-"Fibroblasts Type 2"
seuratObj_clean@meta.data[Extra_Stromal,"New_clusters_clean"]<-"Fibroblasts Type 1"
seuratObj_clean@meta.data$New_clusters_clean<-as.factor(seuratObj_clean@meta.data$New_clusters_clean)
levels(seuratObj_clean@meta.data$New_clusters_clean)

levels(seuratObj_clean@meta.data$New_clusters_clean)<-c("ABC","Brain.Mural","Brain.Mural","Brain.Mural","Brain.Mural",
                                                        "Brain.Pdgfra","Brain.Pdgfra","Ependymal","Ependymal","Ependymal",
                                                        "ChP stromal fibroblasts","ChP base fibroblasts",
                                                        "Pial and arachnoid fibroblasts","Pial and arachnoid fibroblasts","Arachnoid fibroblasts",
                                                        "Dura fibroblasts","Prolif. meningeal fibroblasts","Prolif. meningeal fibroblasts","Ependymal",
                                                        "VECA","VLMC1","VLMC2")

seuratObj_clean@meta.data$New_clusters_clean<-as.factor(as.character(seuratObj_clean@meta.data$New_clusters_clean))
levels(seuratObj_clean@meta.data$New_clusters_clean)

## Sort to order Daan
seuratObj_clean@meta.data$New_clusters_clean<-factor(seuratObj_clean@meta.data$New_clusters_clean, 
                                                     levels=levels(seuratObj_clean@meta.data$New_clusters_clean)[c(6,3,4,5,1,7,2,9,10,12,13,11,8)])

Colorset<-c(brewer.pal(12,"Paired"),"magenta")
Colorset[11]<-"yellow"
U_old<-DimPlot(seuratObj_clean, reduction = "umap", group.by = "New_clusters_clean", repel = T, label = T, label.size = 4, cols = Colorset)

pdf(file = paste0(sampleFolder,"results_merge_subset/Annotation/6_UMAP_dataset_old_annot_clean_",sampleName,".pdf"), width = 12, height = 12)
U_old
dev.off()

# Frequency tables (sliced)
Sample <- seuratObj_clean@meta.data$New_clusters_clean
cluster <- seuratObj_clean@meta.data$annotated_clusters_subset
Aggr <- rep(sampleName,length(cluster))

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(Sample,Aggr))

# Stacked
library("ggthemes")

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/6_SampleDistribution_ggplot2_annot_1_",sampleName,"_clean.pdf"), width = 10, height = 7.5)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white") +
  scale_fill_manual(values=Colorset)
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/6_SampleDistribution_ggplot2_annot_3_",sampleName,"_clean.pdf"), width = 10, height = 7.5)
ggplot(data2, aes(fill=Sample, y=Freq, x=Aggr)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white") +
  scale_fill_manual(values=Colorset)
dev.off()

###################
Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean

## Determine DEG for this annotation
library(future)
plan("multiprocess", workers = 6)
plan()

RNAMarkers_RNAclus <- FindAllMarkers(seuratObj_clean, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_RNAclus$cluster)
saveRDS(RNAMarkers_RNAclus, file=paste0(sampleFolder,"results_merge_subset/Robjects/RNAmarkersList_Old_annotation_cleaned_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['RNAmarkersOldAnnotation']]<-paste0(table(RNAMarkers_RNAclus$cluster)," RNA markers for old annotated cluster ",rownames(table(RNAMarkers_RNAclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrRNAclusters_RNAclus<-names(table(RNAMarkers_RNAclus$cluster))
# totalNrRNAclusters_RNAclusPlusOne<-totalNrRNAclusters_RNAclus
RNAmarkersList_RNAclus<-list()

for(i in totalNrRNAclusters_RNAclus){
  # clusterNr<-i-1
  
  tmp<-RNAMarkers_RNAclus[RNAMarkers_RNAclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  RNAmarkersList_RNAclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
# names(RNAmarkersList_RNAclus)<-paste0("RNAclustersliced",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results_merge_subset/Marker_lists/RNAmarkersList_Old_annotation_cleaned_",sampleName,".xlsx"))

###################

## Determine DEG for this annotation (stricter!!)
RNAMarkers_RNAclus_strict <- FindAllMarkers(seuratObj_clean, assay = "RNA", min.pct = 0.6, only.pos = TRUE)
table(RNAMarkers_RNAclus_strict$cluster)
saveRDS(RNAMarkers_RNAclus_strict, file=paste0(sampleFolder,"results_merge_subset/Robjects/RNAmarkersList_Old_annotation_cleaned_",sampleName,"_strict.rds"))

### Add to diagnostics
diagnostics[['RNAmarkersOldAnnotationStrict']]<-paste0(table(RNAMarkers_RNAclus_strict$cluster)," RNA markers for old annotated cluster strict ",rownames(table(RNAMarkers_RNAclus_strict$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_subset/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrRNAclusters_RNAclus_strict<-names(table(RNAMarkers_RNAclus_strict$cluster))
# totalNrRNAclusters_RNAclusPlusOne<-totalNrRNAclusters_RNAclus_strict
RNAmarkersList_RNAclus_strict<-list()

for(i in totalNrRNAclusters_RNAclus_strict){
  # clusterNr<-i-1
  
  tmp<-RNAMarkers_RNAclus_strict[RNAMarkers_RNAclus_strict$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  RNAmarkersList_RNAclus_strict[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
# names(RNAmarkersList_RNAclus_strict)<-paste0("RNAclustersliced",0:totalNrRNAclusters_RNAclus_strict)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus_strict, file =paste0(sampleFolder, "results_merge_subset/Marker_lists/RNAmarkersList_Old_annotation_cleaned_",sampleName,"_strict.xlsx"))

###################

## Cross ref with markers CPE in LPSnegaggr and add Ttr (worst culprit of contamination!)
library(openxlsx)
Markers_CPE<-read.xlsx("../LpsNegAggr/results/QC/groupClusters_LpsNegAggr.xlsx", sheet = "group5_cl0.1.2.3.5.7")

Filter1<-intersect(Markers_CPE$gene,RNAMarkers_RNAclus_strict[which(RNAMarkers_RNAclus_strict$cluster == "ChP stromal fibroblasts"),"gene"])
Filter2<-intersect(Markers_CPE$gene,RNAMarkers_RNAclus_strict[which(RNAMarkers_RNAclus_strict$cluster == " ChP base fibroblasts"),"gene"])

## Filter lists for CPE markers
RNAmarkersList_RNAclus_strict_filtered<-RNAmarkersList_RNAclus_strict
RNAmarkersList_RNAclus_strict_filtered$`ChP stromal fibroblasts`<-RNAmarkersList_RNAclus_strict_filtered$`ChP stromal fibroblasts`[setdiff(RNAmarkersList_RNAclus_strict_filtered$`ChP stromal fibroblasts`$gene,c(Markers_CPE$gene,"Ttr")),]
RNAmarkersList_RNAclus_strict_filtered$`ChP base fibroblasts`<-RNAmarkersList_RNAclus_strict_filtered$`ChP base fibroblasts`[setdiff(RNAmarkersList_RNAclus_strict_filtered$`ChP base fibroblasts`$gene,c(Markers_CPE$gene,"Ttr")),]

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus_strict_filtered, file =paste0(sampleFolder, "results_merge_subset/Marker_lists/RNAmarkersList_Old_annotation_cleaned_",sampleName,"_strict_filtered.xlsx"))

###################
top3_score<- as.character()
for (pop in 1:length(RNAmarkersList_RNAclus_strict_filtered)) {
  extra_genes<-head(RNAmarkersList_RNAclus_strict_filtered[[pop]]$gene,3)
  top3_score<-c(top3_score,extra_genes)
}
write.xlsx(as.data.frame(top3_score, drop =F),paste0(sampleFolder,"results_merge_subset/Annotation/08_Dotplot_top_genes_Fig1F.xlsx"))

## Create dotplot of top markers old annotation
Colors_dotplot<-c("#071AE5","#F50635") #030720

wantedGenes<-top3_score
wantedGenes<-unique(wantedGenes)
D1<-DotPlot(seuratObj_clean, features = wantedGenes, cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/08_Dotplot_top3_markers_score_strict_filtered.pdf"), width = 15, height = 10)
D1
dev.off()

###################

## Revert back Idents of object
Idents(seuratObj_clean)<-seuratObj_clean@meta.data$annotated_clusters_subset
DimPlot(seuratObj_clean, reduction = "umap", label = T, label.size = 4)

######################################################################################################

##Increase dot size for paper
library(RColorBrewer)
Colorset<-c(brewer.pal(12,"Paired"),"magenta")
Colorset[11]<-"yellow"
U_old_annot_clean_v2<-DimPlot(seuratObj_clean, reduction = "umap", label = T, repel = T, label.size = 4, cols = Colorset,
                              pt.size = 1, group.by = "New_clusters_clean")

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/9_UMAP_old_annot_clean_bigger_dots_",sampleName,".pdf"), width = 10, height = 10)
U_old_annot_clean_v2
dev.off()

###################

## Final plots paper
library(viridis)

## Feature plots panel D
features<-c("Dcn", "Dpep1", "Igfbp6" , "Cldn11")

pdf(file=paste0(sampleFolder,"results_merge_subset/Feature_plots/Feature_plot_paper_4_markers_check_",sampleName,"_viridisC_ordered.pdf"), height = 10, width = 8)
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
pdf(file = paste0(sampleFolder,"results_merge_subset/Feature_plots/Featureplot_modulescore_Cldn11_Igfbp6_",sampleName,".pdf"), width = 8, height = 10)
F2
dev.off()

###################
## DE analysis with new clean metadata! 
# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)

# Clusters to compare!!!!
"CP Stromal"
"CP Stalk"

Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean

Fib1_All_vs_Fib2_All<-FindMarkers(seuratObj_clean, ident.1 = "ChP stromal fibroblasts", ident.2 = "ChP base fibroblasts",
                                  min.pct = 0.10, logfc.threshold = 0.30, only.pos = F)

##### Create list
listDEgenesFBs<-tibble::lst(Fib1_All_vs_Fib2_All)

##Add geneSymbol in column (for the export)
listDEgenesFBs<-lapply(listDEgenesFBs,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenesFBs<-lapply(listDEgenesFBs, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDEgenesFBs<-lapply(listDEgenesFBs, function(x){rbind(x[x$avg_log2FC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_log2FC),
                                                         x[x$avg_log2FC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_log2FC))})
# listDEgenesFBs<-lapply(listDEgenesFBs, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_log2FC)})
##Sort on logFC
listDEgenesFBs<-lapply(listDEgenesFBs,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDEgenesFBs,file=paste0(sampleFolder,"results_merge_subset/Robjects/Markers_for_our_FBs_paper_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDEgenesFBs, paste0(sampleFolder,"results_merge_subset/Marker_lists/Markers_for_our_FBs_paper_",sampleName,".xlsx"))

## Revert annotation
Idents(seuratObj_clean)<-seuratObj_clean@meta.data$annotated_clusters_subset #Revert

###################

### GOEA ###
library(clusterProfiler)

## Background
seuratObj8<-readRDS(file="Urvb_datasets_rebuttal/Robjects/seuratObj_final_Urvb_datasets_rebuttal.rds") #Urvb dataset for full merge
Background_scRNAseq<-rownames(seuratObj8@assays$RNA@counts)

## Stromal FBs
Test<-enrichGO(
  as.character(listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_log2FC > 0),"gene"]),
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

D1<-dotplot(Test, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Stalk FBs
Test2<-enrichGO(
  as.character(listDEgenesFBs$Fib1_All_vs_Fib2_All[which(listDEgenesFBs$Fib1_All_vs_Fib2_All$avg_log2FC < 0),"gene"]),
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

D2<-dotplot(Test2, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Combine
D_comp<-merge_result(list(Stromal_CP=Test, Stalk_CP=Test2)) %>%
  dotplot(., split="ONTOLOGY", showCategory=15) + facet_grid(ONTOLOGY~., scale="free")


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


# ###############

### Save results ###
pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/10_Dotplot_FB1s_All_vs_FB2s_All_",sampleName,".pdf"), width = 10, height = 10)
D1
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/10_Dotplot_FB2s_All_vs_FB1s_All_",sampleName,".pdf"), width = 10, height = 10)
D2
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_",sampleName,".pdf"), width = 10, height = 15)
D_comp
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_",sampleName,"_without_MF.pdf"), width = 10, height = 15)
D_comp_filtered
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/10_Dotplot_Comparison_FB1s_and_FB2s_All_order_Roos_",sampleName,"_without_MF.pdf"), width = 10, height = 15)
D_comp_filtered_v2
dev.off()

####

### Save objects ###
saveRDS(Test,file=paste0(sampleFolder,"results_merge_subset/Robjects/EnrichGO_FB1s_All_vs_FB2s_All_",sampleName,".rds"))
saveRDS(Test2,file=paste0(sampleFolder,"results_merge_subset/Robjects/EnrichGO_FB2s_All_vs_FB1s_All_",sampleName,".rds"))

write.xlsx(Test@result,file=paste0(sampleFolder,"results_merge_subset/Annotation/10_EnrichGO_FB1s_All_vs_FB2s_All_",sampleName,".xlsx"))
write.xlsx(Test2@result,file=paste0(sampleFolder,"results_merge_subset/Annotation/10_EnrichGO_FB2s_All_vs_FB1s_All_",sampleName,".xlsx"))


#########################################################################################################
#########################################################################################################

## Check matrixome collagen genes (16/12/21)
Matrixome_df<-read.xlsx("Documentation/MGIBatchReport_matrixome.xlsx", sheet = "Matrixome")
Matrixome_genes<-Matrixome_df$Symbol
intersect(Matrixome_genes,rownames(seuratObj_clean)) #All!

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = Matrixome_genes, cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/12_Dotplot_matrixome_",sampleName,".pdf"), width = 15, height = 10)
D1
dev.off()

# Create new version with a few extra genes
Matrixome_genes_v2<-unique(c(Matrixome_genes,"Lamb2","Fbln1", "Col6a2", "Col12a1"))
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_clean, group.by = "New_clusters_clean", features = Matrixome_genes_v2, cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/12_Dotplot_matrixome_v2_rebuttal_",sampleName,".pdf"), width = 15, height = 10)
D1
dev.off()

#########################################################################################################

## Check FB marker genes 
Figure_genes<-c("Steap4", "Heph", "Cp", "Cxcl12", "Dpep1", "Itga8", "Itgbl1", "Tjp1", "Gja1", "Gjb2", "Gjb6", "Cdh1", "Cdh5", "Cldn11") #Orig
Figure_genes2<-c("Trf","Cp","Steap4", "Heph","Cxcl12", "Dpep1", "Itga8", "Itgbl1", "Tjp1", "Gja1", "Gjb2", "Gjb6", "Cdh1", "Cdh5", "Cldn11",
                 "Pmp22","Tjp2","Pkp4","Perp","Cav1","Cxadr","Lsr","Wnt5a","S100a10","Plec","Kdr","Ntrk2","Sipa1l1","Fam107a","Flrt2") #rebuttal
# S100a10/Wnt5a/Plec/Tjp1/Kdr/Pmp22/Gjb2/Pkp4/Ntrk2/Cldn11/Cdh5/Sipa1l1/Gja1/Fam107a/Gjb6/Flrt2/Cav1/Lsr/Cdh1/Fn1

Idents(seuratObj_clean)<-seuratObj_clean@meta.data$New_clusters_clean
levels(Idents(seuratObj_clean))

seuratObj_subset<-subset(seuratObj_clean, idents = c("ChP base fibroblasts", "ChP stromal fibroblasts"))

seuratObj_subset$New_clusters_clean<-factor(as.character(seuratObj_subset$New_clusters_clean), levels = c("ChP base fibroblasts", "ChP stromal fibroblasts"))

levels(seuratObj_subset$New_clusters_clean)<-c("Base barrier cells","Stromal")

# Create dotplot
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_subset, group.by = "New_clusters_clean", features = Figure_genes, cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/16_Dotplot_figure3B_",sampleName,".pdf"), width = 10, height = 6)
D1
dev.off()

# Bigger dotplot 
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObj_subset, group.by = "New_clusters_clean", features = Figure_genes, cols = Colors_dotplot,
            dot.scale = 14) + RotatedAxis() #dot.scale = 8

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/16_Dotplot_figure4B_bigger_dots_v1_",sampleName,".pdf"), width = 10, height = 6)
D1
dev.off()

# Bigger dotplot new genes rebuttal
D2<-DotPlot(seuratObj_subset, group.by = "New_clusters_clean", features = Figure_genes2, cols = Colors_dotplot,
            dot.scale = 13) + RotatedAxis() #dot.scale = 8

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/16_Dotplot_figure4B_bigger_dots_v2_",sampleName,".pdf"), width = 16, height = 6)
D2
dev.off()

##############

## Extra: plots for rebuttal for quality data and M&M

## Convert to sce
sce<-as.SingleCellExperiment(seuratObj_clean)

##### Get mitochondrial genes #####
is.mito <- grepl("^mt-", rownames(sce))
sum(is.mito)
##13
rownames(sce)[is.mito]

##### Calculate QC metrics #####
### => pData(sce) is created
sce <- perCellQCMetrics(sce, subsets=list(Mt=rownames(sce)[is.mito]))

##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("orig.ident"=seuratObj_clean$sample_origin2,
                     "nGene"=sce$detected,"nUMI"=sce$sum,"percent.mito"=sce$subsets_Mt_percent, 
                     stringsAsFactors = F)

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/17_QC_stats_rebuttal_",sampleName,".pdf"))
for (i in levels(seuratObj_clean$sample_origin2)[c(1,2,4,5,6)]){
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

pdf(file=paste0(sampleFolder,"results_merge_subset/Annotation/17_QC_stats_rebuttal_extra_",sampleName,".pdf"))
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

#########################################################################

##### Read object
seuratObj_clean <- readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_clean_",sampleName,"_harmony_RNA.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,"_harmony_RNA.rds"))

##### Save object
saveRDS(seuratObj_clean, file=paste0(sampleFolder,"Robjects/seuratObj_clean_",sampleName,"_harmony_RNA.rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,"_harmony_RNA.rds"))

#########################################################################
## Create diet object for online tool
seuratObj_clean_diet<-DietSeurat(seuratObj_clean, counts = T, data = T, scale.data = F,
                                 assays = c("RNA"), dimreducs = "umap", graphs = NULL)

## New metadata names
seuratObj_clean_diet$Original_annotation<-seuratObj_clean_diet$New_clusters_clean
seuratObj_clean_diet$Origin<-seuratObj_clean_diet$sample_origin2
seuratObj_clean_diet$New_annotation<-seuratObj_clean_diet$annotated_clusters_subset
seuratObj_clean_diet$Numbered_annotation<-seuratObj_clean_diet$harmony_clusters_subset

## Metadata columns
# "orig.ident" "Phase" "SCT_clusters" "ADT_clusters" "annotated_clusters_subset" "final_annotation2021"  
DimPlot(seuratObj_clean_diet, reduction = "umap", label = T, repel = T, group.by = "Origin", label.size = 5,
        # cols =gg_color_hue(10)) 
        cols =Colorset1)
library(RColorBrewer)
Colorset<-c(brewer.pal(12,"Paired"),"#FF00FF")
Colorset[11]<-"#FFFF00"
Colorset1<-c(brewer.pal(8,"Set1")[c(1,7,2,3,4,5)]) #origin

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

All<-c(levels(as.factor(seuratObj_clean_diet$Original_annotation)),levels(as.factor(seuratObj_clean_diet$Origin)),
       levels(as.factor(seuratObj_clean_diet$New_annotation)),levels(seuratObj_clean_diet$Numbered_annotation))
Color_info<-c(Colorset,Colorset1,
              gg_color_hue(10),gg_color_hue(10))
Metadata_column<-c(rep("Original annotation",13),rep("Origin",6),
                   rep("New annotation",10),rep("Numbered_clusters",10))
Info_Kevin<-as.data.frame(cbind(All,Color_info,Metadata_column))

write.xlsx(Info_Kevin, file =paste0(sampleFolder, "results_merge_subset/Annotation/Info_Kevin_rebuttal_",sampleName,".xlsx"))

##### Save object
saveRDS(seuratObj_clean_diet, file=paste0(sampleFolder,"Robjects/seuratObj_paper_diet_rebuttal_",sampleName,"_2023.rds"))

##### Read object
seuratObj_clean_diet<-readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_paper_diet_rebuttal_",sampleName,"_2023.rds"))

######################################################################################################