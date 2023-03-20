library(Seurat)
library(ggplot2)
library(gridExtra)
library(reticulate)
library(openxlsx)

# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

setwd("/run/media/clintdn/CN1465-DATA/VIB_G_drive") #Linux
sampleName_Human<-"Human_snRNAseq_CP_COVID"
sampleFolder<-paste0(sampleName_Human,"/")
sampleName_old <- "CCA_integrated_Human_and_Mouse_CP_W2"
sampleName <- "BBKNN_integrated_Human_and_Mouse_CP"

#############################################################################################################

## Performed in conda env: env_R_and_Python!!!!!!!

FB.combined_2<- readRDS(file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName_old,".rds"))

## Object contains PCA from integrated data, not from RNA
## Remake object and rerun normalize, scale and PCA
Counts_RNA <- GetAssayData(FB.combined_2[["RNA"]], slot = "counts")
Metadata <- FB.combined_2@meta.data

FB.combined_new <- CreateSeuratObject(counts = Counts_RNA, project = "BBKNN_FB_combined", assay = "RNA")
FB.combined_new@meta.data<-Metadata

FB.combined_new <- NormalizeData(object = FB.combined_new, normalization.method = "LogNormalize", scale.factor = 10000)
FB.combined_new <- FindVariableFeatures(FB.combined_new, assay = "RNA", selection.method = "vst", nfeatures = 2000)
FB.combined_new <- ScaleData(FB.combined_new, verbose = FALSE)
FB.combined_new <- RunPCA(FB.combined_new, npcs = 30, verbose = FALSE)

# Converting the Seurat object to an AnnData file is a two-step process. 
# First, we save the Seurat object as an h5Seurat file. 
# For more details about saving Seurat objects to h5Seurat files, please see this vignette; 
# after the file is saved, we can convert it to an AnnData file for use in Scanpy. 
# Full details about the conversion processes are listed in the manual page for the Convert function

saveRDS(FB.combined_new, file = paste0(sampleFolder,"Robjects/seuratObj_adata_",sampleName,".rds"))

SaveH5Seurat(FB.combined_new, filename = paste0(sampleFolder,"Robjects/seuratObj_adata_",sampleName,".h5Seurat"))
Convert(paste0(sampleFolder,"Robjects/seuratObj_adata_",sampleName,".h5Seurat"), dest = "h5ad")

# At this point, there is no plan to create a BBKNN R package. 
# However, it can be ran quite easily via reticulate. 
# Using the base functions is the same as in python. 
# If you're in possession of a PCA matrix and a batch assignment vector and want to get UMAP coordinates out of it, 
# you can use the following code snippet to do so. 
# The weird PCA computation part and replacing it with your original values is unfortunately necessary 
# due to how AnnData innards operate from a reticulate level. Provide your python path in use_python()

use_python("/home/clintdn/anaconda3/envs/env_R_and_Python/bin/python3")
anndata = import("anndata",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
sc = import("scanpy",convert=FALSE)

# We can view the AnnData file in Scanpy by using the read_h5ad function
adata = sc$read_h5ad(paste0(sampleFolder,"Robjects/seuratObj_adata_",sampleName,".h5ad"))
adata

## Issues following tutorial on bbknn github!!!
# pca<-adata$obsm["X_pca"]
# batch<-adata$obs$Species
# 
# adata_bbknn = anndata$AnnData(X=pca, obs=batch)
# # Error in py_call_impl(callable, dots$args, dots$keywords) : 
# #   ValueError: Cannot convert <class 'pandas.core.series.Series'> to DataFrame
# # 
# # Detailed traceback:
# #   File "/home/clintdn/anaconda3/envs/env_R_and_Python/lib/python3.9/site-packages/anndata/_core/anndata.py", line 308, in __init__
# # self._init_as_actual(
# #   File "/home/clintdn/anaconda3/envs/env_R_and_Python/lib/python3.9/site-packages/anndata/_core/anndata.py", line 502, in _init_as_actual
# #   self._obs = _gen_dataframe(obs, self._n_obs, ["obs_names", "row_names"])
# #   File "/home/clintdn/anaconda3/envs/env_R_and_Python/lib/python3.9/functools.py", line 877, in wrapper
# #   return dispatch(args[0].__class__)(*args, **kw)
# #   File "/home/clintdn/anaconda3/envs/env_R_and_Python/lib/python3.9/site-packages/anndata/_core/anndata.py", line 128, in _
# #   raise ValueError(f"Cannot convert {type(anno)} to DataFrame")
# 
# pca<-FB.combined_new@reductions$pca
# batch<-FB.combined_new@meta.data$Species
# data_bbknn = anndata$AnnData(X=pca, obs=batch)
# # A dimensional reduction object with key PC_ 
# # Number of dimensions: 30 
# # Projected dimensional reduction calculated:  FALSE 
# # Jackstraw run: FALSE 
# # Computed using assay: integrated 
# # Error in py_call_impl(callable, dots$args, dots$keywords) : 
# #   Evaluation error: Unable to convert R object to Python type.

# ## Just work with object and run BBKNN (STANDARD)
# bbknn$bbknn(adata,batch_key="Species") #Adds obsp distances + connectivities and uns neighbors
# sc$tl$umap(adata) #Adds obsm X_umap and uns umap
# umap = py_to_r(adata$obsm[["X_umap"]])

## Just work with object and run BBKNN (EXTRA part for better clustering)
## One may not be armed in a biological grouping at the time of analysis, in which case a coarse clustering 
## (aiming to feature as many batches as possible per cluster) can be used in its place. This results in the following analysis flow:
# bbknn
# clustering
# ridge regression
# pca
# bbknn

bbknn$bbknn(adata,batch_key="Species") #Adds obsp distances + connectivities and uns neighbors
sc$tl$umap(adata)
sc$pl$umap(adata, color='Species')
sc$pl$umap(adata, color='Integrated_annotated_clusters')
sc$tl$leiden(adata, resolution = 0.4) #Leiden clustering
sc$pl$umap(adata, color = 'leiden')

# Armed with this crude biological grouping, let's use it to call bbknn.ridge_regression(), 
# and use the output of that to recompute the PCA. We can then re-run BBKNN and create a new UMAP, which mixes the batches better.

bbknn$ridge_regression(adata, batch_key='Species', confounder_key='leiden')
sc$pp$pca(adata)
bbknn$bbknn(adata,batch_key="Species")
sc$tl$umap(adata)
sc$pl$umap(adata, color='Species')
sc$pl$umap(adata, color='Integrated_annotated_clusters')

umap = py_to_r(adata$obsm[["X_umap"]])
pca = py_to_r(adata$obsm[["X_pca"]])

# saveRDS(umap, file = paste0(sampleFolder,"Robjects/Ridge_regressed_UMAP_",sampleName,".rds")) #old

## Add in rownames otherwise error later
rownames(umap)<-rownames(FB.combined_new[["pca"]])
rownames(pca)<-rownames(FB.combined_new[["pca"]])

## Redo and save pca too! Better clusters?? (January 27th)
saveRDS(umap, file = paste0(sampleFolder,"Robjects/Ridge_regressed_UMAP_new_",sampleName,".rds"))
saveRDS(pca, file = paste0(sampleFolder,"Robjects/Ridge_regressed_PCA_new_",sampleName,".rds"))

## Add into object and finish pipeline
umap_dr<-CreateDimReducObject(umap, key="UMAP_")
pca_dr<-CreateDimReducObject(pca, key="pythonPCA_")

FB.combined_new@reductions$umap<-umap_dr
FB.combined_new@reductions$pca<-pca_dr #overwrite old pca!!!!!!!!!!! #New

# FB.combined_new <- FindNeighbors(FB.combined_new, reduction = "pca", dims = 1:30) #old
FB.combined_new <- FindNeighbors(FB.combined_new, reduction = "pca", dims = 1:30) #new
FB.combined_new <- FindClusters(FB.combined_new, resolution = 0.5)

saveRDS(FB.combined_new, file = paste0(sampleFolder,"Robjects/Ridge_regressed_seuratObj_new_",sampleName,".rds"))

################################################################################################################

## Load objects
FB.combined_basic <- readRDS(file = paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))
# FB.combined_regressed <- readRDS(file = paste0(sampleFolder,"Robjects/Ridge_regressed_seuratObj_",sampleName,".rds")) #old
FB.combined_regressed <- readRDS(file = paste0(sampleFolder,"Robjects/Ridge_regressed_seuratObj_new_",sampleName,".rds")) #new

## Save object
saveRDS(FB.combined_regressed, file = paste0(sampleFolder,"Robjects/Ridge_regressed_seuratObj_new_",sampleName,".rds")) #new

## Create plots
umapPlot<-DimPlot(FB.combined_basic, reduction = "umap", label = T, label.size = 8)
umapPlotSplit<-DimPlot(FB.combined_basic, reduction = "umap", label = F, group.by="orig.ident")

ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
       file=paste0(sampleFolder,"BBKNN/UMAP_1_30.png"), width = 20, height=10)

umapPlot<-DimPlot(FB.combined_regressed, reduction = "umap", label = T, label.size = 8)
umapPlotSplit<-DimPlot(FB.combined_regressed, reduction = "umap", label = F, group.by="orig.ident")

ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
       file=paste0(sampleFolder,"BBKNN/UMAP_regressed_new_1_30.png"), width = 20, height=10)

## Check features
FeaturePlot(FB.combined_basic, features = "Igfbp6", min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(FB.combined_regressed, features = "Igfbp6", min.cutoff = "q10", max.cutoff = "q90")

################################################################################################################

### Find RNAmarkers for every Integrated cluster compared to all remaining cells, report only the positive ones
library(future)
plan("multiprocess", workers = 8)

RNAMarkers_RNAclus_basic <- FindAllMarkers(FB.combined_basic, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_RNAclus_basic$cluster)
saveRDS(RNAMarkers_RNAclus_basic, file=paste0(sampleFolder,"BBKNN/Robjects/RNAmarkersList_RNAclus_basic_",sampleName,".rds"))

### Create list with markers
totalNrRNAclusters_RNAclus<-max(as.numeric(names(table(RNAMarkers_RNAclus_basic$cluster))))
totalNrRNAclusters_RNAclusPlusOne<-totalNrRNAclusters_RNAclus+1
RNAmarkersList_RNAclus_basic<-list()

for(i in 1:totalNrRNAclusters_RNAclusPlusOne){
  clusterNr<-i-1
  
  tmp<-RNAMarkers_RNAclus_basic[RNAMarkers_RNAclus_basic$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  RNAmarkersList_RNAclus_basic[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(RNAmarkersList_RNAclus_basic)<-paste0("Integratedcluster",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus_basic, file =paste0(sampleFolder, "BBKNN/RNAmarkersList_RNAclus_basic_",sampleName,".xlsx"))

################################################################################################################

### Find RNAmarkers for every Integrated cluster compared to all remaining cells, report only the positive ones
RNAMarkers_RNAclus_regressed <- FindAllMarkers(FB.combined_regressed, assay = "RNA", only.pos = TRUE)
table(RNAMarkers_RNAclus_regressed$cluster)
saveRDS(RNAMarkers_RNAclus_regressed, file=paste0(sampleFolder,"BBKNN/Robjects/RNAmarkersList_RNAclus_regressed_new_",sampleName,".rds"))

### Create list with markers
totalNrRNAclusters_RNAclus<-max(as.numeric(names(table(RNAMarkers_RNAclus_regressed$cluster))))
totalNrRNAclusters_RNAclusPlusOne<-totalNrRNAclusters_RNAclus+1
RNAmarkersList_RNAclus_regressed<-list()

for(i in 1:totalNrRNAclusters_RNAclusPlusOne){
  clusterNr<-i-1
  
  tmp<-RNAMarkers_RNAclus_regressed[RNAMarkers_RNAclus_regressed$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  RNAmarkersList_RNAclus_regressed[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(RNAmarkersList_RNAclus_regressed)<-paste0("Integratedcluster",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus_regressed, file =paste0(sampleFolder, "BBKNN/RNAmarkersList_RNAclus_regressed_new_",sampleName,".xlsx"))

################################################################################################################

######################################
##### Markers annotated clusters
########################################

## Annotate clusters reference
DimPlot(FB.combined_basic, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(FB.combined_regressed, reduction = "umap", label = TRUE, repel = TRUE)

DimPlot(FB.combined_basic, reduction = "umap", group.by = "Integrated_annotated_clusters", label = TRUE, repel = TRUE)
D1<-DimPlot(FB.combined_regressed, reduction = "umap", group.by = "Integrated_annotated_clusters", label = TRUE, repel = TRUE)
ggsave(D1, file=paste0(sampleFolder,"BBKNN/5_UMAP_regressed_new_CCA_annotation_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

VlnPlot(FB.combined_basic, features = c("Igfbp6","Cldn11"))
VlnPlot(FB.combined_regressed, features = c("Igfbp6","Cldn11"))
# FB.combined_regressed@meta.data$Integrated_RNA_clusters<-factor(FB.combined_regressed@meta.data$Integrated_RNA_clusters,levels=sort(as.numeric(levels(FB.combined_regressed@meta.data$Integrated_RNA_clusters)))) #reorder levels

FB.combined_regressed@meta.data$BBKNN_clusters <- FB.combined_regressed@active.ident
FB.combined_regressed@meta.data$annotated_BBKNN_clusters <- FB.combined_regressed@active.ident

levels(FB.combined_regressed@meta.data$annotated_BBKNN_clusters) <- c("Stromal FBs Human","Stromal FBs Mouse","Stromal FBs Mouse",
                                                                      "Stromal FBs Human","Stromal FBs Human","Stromal FBs Mouse",
                                                                      "Other FBs Human 1","Stalk FBs","Other FBs Human 2",
                                                                      "Other FBs Human 3","Hsp+ Stromal FBs Human", "ABCs")

U_annot<-DimPlot(FB.combined_regressed, reduction = "umap", label = T, group.by = "annotated_BBKNN_clusters", repel = T, label.size = 4) + NoLegend()
ggsave(U_annot, file=paste0(sampleFolder,"BBKNN/5_UMAP_regressed_new_annotated_clusters_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

################################################################################################################

## Check Igfbp6 and Cldn11 expression in regressed object

# Featureplots on full object
# Module score Seurat
library(viridis)

FB_stalk_signature <- list(c("Cldn11","Igfbp6"))

FB.combined_regressed <- AddModuleScore(object = FB.combined_regressed, assay = "RNA", features = FB_stalk_signature, name = "FB_stalk_signature_score")
F2 <- FeaturePlot(object = FB.combined_regressed, features = "FB_stalk_signature_score1", 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order = T)  + scale_color_viridis(option = "D")
pdf(file = paste0(sampleFolder,"BBKNN/2_UMAP_new_Featureplot_modulescore_Cldn11_Igfbp6_",sampleName,".pdf"), width = 15, height = 10)
F2
dev.off()

# Idents(FB.combined_regressed)<-FB.combined_regressed@meta.data$annotated_BBKNN_clusters
# F2.5 <- FeaturePlot(object = FB.combined_regressed, features = "FB_stalk_signature_score1", order = T, label = T, repel = T, label.size = 3, cols = c("Yellow","Red"))
# pdf(file = paste0(sampleFolder,"BBKNN/2_UMAP_new_Featureplot_modulescore_Cldn11_Igfbp6_labeled_",sampleName,".pdf"), width = 15, height = 10)
# F2.5
# dev.off()

# # Extra featureplots
# F1 <- FeaturePlot(object = FB.combined_regressed, features = c("Cldn11","Igfbp6","Dpep1","Alpl","Cdh11"), order = T, pt.size = 1, cols = c("Yellow","Red"), combine = F)
# pdf(file = paste0(sampleFolder,"BBKNN/2_UMAP_new_Featureplots_FB_markers_RNA_",sampleName,".pdf"), width = 15, height = 10)
# F1
# dev.off()

## Featureplots paper (May 2022)
features<-c("Cldn11","Igfbp6","Dpep1","Alpl","Cdh11","Dcn")

pdf(file=paste0(sampleFolder,"BBKNN/2_Feature_plot_paper_4_markers_",sampleName,"_viridisC_ordered.pdf"), height = 10, width = 15)
for (feature in features) {
  F1<-FeaturePlot(object = FB.combined_regressed, features =feature, cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) +
    scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

# Create new clusters: split on source
FB.combined_regressed@meta.data$newClustersTmp<-FB.combined_regressed@meta.data$annotated_BBKNN_clusters
FB.combined_regressed@meta.data$newClusters<-paste0(FB.combined_regressed@meta.data$newClustersTmp,"_",FB.combined_regressed@meta.data$Species)
head(FB.combined_regressed@meta.data)

# Subset stalk and violin plot
Idents(FB.combined_regressed)<-FB.combined_regressed@meta.data$newClusters
# FB.stalkandstromal<-subset(FB.combined_regressed, ident = c("Stalk FBs_Mouse", "Stalk FBs_Human","ABCs_Mouse","Stromal FBs Human_Human","Stromal FBs Mouse_Mouse"))
FB.stalkandstromal<-subset(FB.combined_regressed, ident = c("Stalk FBs_Mouse", "Stalk FBs_Human","Stromal FBs Human_Human","Stromal FBs Mouse_Mouse")) #11/03/22

# Idents(FB.stalkandstromal)<-factor(FB.stalkandstromal@active.ident, levels = levels(Idents(FB.stalkandstromal))[c(3,5,1,4,2)])
Idents(FB.stalkandstromal)<-factor(FB.stalkandstromal@active.ident, levels = levels(Idents(FB.stalkandstromal))[c(4,1,3,2)]) #110322
# levels(Idents(FB.stalkandstromal))<-c("ABCs_Mouse","Stalk FBs_Human","Stalk FBs_Mouse","Stromal FBs_Human","Stromal FBs_Mouse")
levels(Idents(FB.stalkandstromal))<-c("Stalk FBs_Human","Stalk FBs_Mouse","Stromal FBs_Human","Stromal FBs_Mouse") #110322

DimPlot(FB.stalkandstromal)

library(RColorBrewer)
# ColorSet<-brewer.pal(n = 6, name = "Paired")[2:6]
ColorSet<-brewer.pal(n = 6, name = "Paired")[3:6] #110322

V1<-VlnPlot(FB.stalkandstromal, features = c("Igfbp6","Cldn11","Cdh11","Alpl"), assay = "RNA", ncol = 4, cols = ColorSet) & theme(plot.margin = margin(0.5,0.5,0.5,1, "cm"))

pdf(file=paste0(sampleFolder,"BBKNN/4_Violinplot_new_v2_",sampleName,".pdf"), width = 25, height= 10) #110322
V1
dev.off()

## UMAP with ventricle info (Daan prefers that!)
# Update Ventricle metadata (all human data is LV!!!)
FB.combined_regressed$Ventricle<-FB.combined_regressed$orig.ident

FB.combined_regressed@meta.data[which(FB.combined_regressed@meta.data$Ventricle == "RVD1_LpsNegFour" |
                                FB.combined_regressed@meta.data$Ventricle == "RVD5_Y4V" | FB.combined_regressed@meta.data$Ventricle == "RVD7_O4V"),
                                "Ventricle"]<-"4V"

FB.combined_regressed@meta.data[which(FB.combined_regressed@meta.data$Ventricle == "ct" | FB.combined_regressed@meta.data$Ventricle == "RVD2_LpsNegLat" |
                                FB.combined_regressed@meta.data$Ventricle == "cv" | FB.combined_regressed@meta.data$Ventricle == "RVD6_YLV" |
                                FB.combined_regressed@meta.data$Ventricle == "flu" | FB.combined_regressed@meta.data$Ventricle == "RVD8_OLV"),
                                "Ventricle"]<-"LV"

U_ventricle<-DimPlot(FB.combined_regressed, reduction = "umap", label = T, group.by = "Ventricle", repel = T, label.size = 4) + NoLegend()
ggsave(U_ventricle, file=paste0(sampleFolder,"BBKNN/5_UMAP_new_ventricle_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

################################################################################

## Paper figures
U_annot<-DimPlot(FB.combined_regressed, reduction = "umap", label = T, group.by = "annotated_BBKNN_clusters", repel = T, label.size = 4)

umapPlotSplit<-DimPlot(FB.combined_regressed, reduction = "umap", label = F, pt.size = 0.5, group.by="orig.ident")

umapPlotSplit2<-DimPlot(FB.combined_regressed, reduction = "umap", label = F, pt.size = 1, group.by="Species")

pdf(file=paste0(sampleFolder,"BBKNN/5_UMAP_regressed_new_annotated_clusters_",sampleName,".pdf"), width =16, height= 10)
U_annot
dev.off()

pdf(file=paste0(sampleFolder,"BBKNN/5_UMAP_regressed_origin_",sampleName,".pdf"), width =16, height= 10)
umapPlotSplit
dev.off()

pdf(file=paste0(sampleFolder,"BBKNN/5_UMAP_regressed_Species_",sampleName,".pdf"), width =16, height= 10)
umapPlotSplit2
dev.off()

## Adapt annotation for paper
FB.combined_regressed$annotated_BBKNN_clusters_v2<-FB.combined_regressed$annotated_BBKNN_clusters
levels(FB.combined_regressed$annotated_BBKNN_clusters_v2)[7]<-"Other FBs Human 4"

U_annot_v2<-DimPlot(FB.combined_regressed, reduction = "umap", label = T, group.by = "annotated_BBKNN_clusters_v2", 
                 repel = T, label.size = 4, pt.size = 1)

pdf(file=paste0(sampleFolder,"BBKNN/5_UMAP_regressed_new_annotated_clusters_v2_",sampleName,".pdf"), width =16, height= 10)
U_annot_v2
dev.off()

################################################################################

## Sample distribution across clusters

# create a dataset
Sample <- FB.combined_regressed@meta.data$Species
cluster <- FB.combined_regressed$annotated_BBKNN_clusters
# Aggr <- rep(experiment,length(cluster)) 

data <- data.frame(table(Sample, cluster))
# data2 <- data.frame(table(cluster,Aggr))

# barplotAggr(seuratObj, listLabels)

# Stacked
png(file=paste0(sampleFolder,"BBKNN/3_SampleDistribution_ggplot2_new_annot_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"BBKNN/3_SampleDistribution_ggplot2_new_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
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

pdf(file=paste0(sampleFolder,"BBKNN/3_SampleDistribution_ggplot2_new_annotated_",sampleName,"_adjusted.pdf"), width =8, height= 6)
S2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #Remove grid lines
dev.off()

#############################################################################################################
Idents(FB.combined_regressed)<-FB.combined_regressed$annotated_BBKNN_clusters

## Chekc for conserved markers between species
Stalk.markers <- FindConservedMarkers(FB.combined_regressed, ident.1 = "Stalk FBs", grouping.var = "Species", verbose = FALSE)
head(Stalk.markers)

## Save results
saveRDS(Stalk.markers, file=paste0(sampleFolder,"BBKNN/Robjects/Conserved_stalk_markers_new_",sampleName,".rds"))
write.xlsx(Stalk.markers, file =paste0(sampleFolder, "BBKNN/Conserved_stalk_markers_new_",sampleName,".xlsx"), row.names = T)

## Read results
Stalk.markers<-readRDS(file=paste0(sampleFolder,"BBKNN/Robjects/Conserved_stalk_markers_new_",sampleName,".rds"))

## Get top markers
Stalk.markers.top<-Stalk.markers[which(Stalk.markers$Mouse_avg_log2FC > 1 & Stalk.markers$Human_avg_log2FC >1),]
Stalk.markers.top<-Stalk.markers.top[order(Stalk.markers.top$Mouse_avg_log2FC, decreasing = T),]

## Plot out data
F1 <- FeaturePlot(FB.combined_regressed, features = head(rownames(Stalk.markers.top), n =6), pt.size = 1,  min.cutoff = 'q2', max.cutoff = 'q98', ncol = 3)
F2 <- FeaturePlot(FB.combined_regressed, features = head(rownames(Stalk.markers.top), n =6), order = T,  pt.size = 1,  min.cutoff = 'q2', max.cutoff = 'q98', ncol = 3)

## Make subset of Stalk data
FB.Stalk_only<-subset(FB.combined_regressed, idents = "Stalk FBs")

## Try dotplots
markers.to.plot <- rownames(Stalk.markers.top)[1:6]

D1 <- DotPlot(FB.Stalk_only, features = markers.to.plot, cols = c("blue", "red"), assay = "RNA", split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot Stalk Only")

D2 <- DotPlot(FB.combined_regressed, features = markers.to.plot, cols = c("blue", "red"), assay = "RNA", split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot All")

D3 <- DotPlot(FB.Stalk_only, group.by = "newClusters", assay = "RNA", features = markers.to.plot) +
  RotatedAxis() + ggtitle("Dotplot Stalk Only with scale")

D4 <- DotPlot(FB.combined_regressed, group.by = "newClusters", assay = "RNA", features = markers.to.plot) +
  RotatedAxis() + ggtitle("Dotplot All with scale")

## Save plots
pdf(file=paste0(sampleFolder,"BBKNN/1_Featureplots_new_Stalk_Markers_non_ordered_",sampleName,".pdf"), width = 20, height = 15)
F1
dev.off()

pdf(file=paste0(sampleFolder,"BBKNN/1_Featureplots_new_Stalk_Markers_ordered_",sampleName,".pdf"), width = 20, height = 15)
F2
dev.off()

pdf(file=paste0(sampleFolder,"BBKNN/1_Dotplots_new_Stalk_Markers_",sampleName,".pdf"), width = 10, height = 8)
D1
D2
D3
D4
dev.off()

## Second version of markers
# Filter on adj P-value for human and mouse and order according to Max P-value -> top 5

Top_markers_v2<- c("S100b","Cpe","Vim", "Cystm1","Igfbp6") #"Fth1", 

D2_v2 <- DotPlot(FB.combined_regressed, features = Top_markers_v2, assay = "RNA", cols = c("blue", "red"), split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot paper")
D4_v2 <- DotPlot(FB.combined_regressed, group.by = "newClusters", assay = "RNA", features = Top_markers_v2) +
  RotatedAxis() + ggtitle("Dotplot All with scale v2")

pdf(file=paste0(sampleFolder,"BBKNN/1_Dotplots_new_Stalk_Markers_v2_",sampleName,".pdf"), width = 10, height = 8)
D2_v2
D4_v2
dev.off()

pdf(file=paste0(sampleFolder,"BBKNN/1_Dotplots_new_Stalk_Markers_adjPV_paper_",sampleName,".pdf"), width = 10, height = 8)
D2_v2
dev.off()

####################

## Scatterplot Kia paper (07/07/22)

# A basic scatterplot with color depending on Species
Stalk.markers$name<-rownames(Stalk.markers)
ggplot(Stalk.markers, aes(x=-log(Mouse_p_val_adj), y=-log(Human_p_val_adj), label = name)) + 
  geom_point(size=1)  +geom_text(hjust=0, vjust=-0.1)

## New DE analysis Stalk vs Stromal FBs
Idents(FB.combined_regressed)<-FB.combined_regressed@meta.data[["NewClusters_combo_annotated"]]
Mouse_DE<-FindMarkers(FB.combined_regressed, ident.1 = "Stalk FBs_Mouse", ident.2 = c("ABCs_Mouse","Stromal FBs Mouse_Mouse"))
Human_DE<-FindMarkers(FB.combined_regressed, ident.1 = "Stalk FBs_Human", ident.2 = c("Hsp+ Stromal FBs Human_Human","Other FBs Human 1_Human","Other FBs Human 2_Human",
                                                                                      "Other FBs Human 3 _Human","Stromal FBs Human_Human","Stromal FBs Mouse_Human"))
Mouse_DE_v2<-FindMarkers(FB.combined_regressed, ident.1 = "Stalk FBs_Mouse", ident.2 = c("Stromal FBs Mouse_Mouse"),
                         logfc.threshold = 0, min.pct = 0)
Human_DE_v2<-FindMarkers(FB.combined_regressed, ident.1 = "Stalk FBs_Human", ident.2 = c("Stromal FBs Human_Human"),
                         logfc.threshold = 0, min.pct = 0)
## Save results
saveRDS(Mouse_DE, file=paste0(sampleFolder,"BBKNN/Robjects/Mouse_DE_v1_",sampleName,".rds"))
write.xlsx(Mouse_DE, file =paste0(sampleFolder, "BBKNN/Mouse_DE_v1_",sampleName,".xlsx"), row.names = T)
saveRDS(Mouse_DE_v2, file=paste0(sampleFolder,"BBKNN/Robjects/Mouse_DE_v2_",sampleName,".rds"))
write.xlsx(Mouse_DE_v2, file =paste0(sampleFolder, "BBKNN/Mouse_DE_v2_",sampleName,".xlsx"), row.names = T)
saveRDS(Human_DE, file=paste0(sampleFolder,"BBKNN/Robjects/Human_DE_v1_",sampleName,".rds"))
write.xlsx(Human_DE, file =paste0(sampleFolder, "BBKNN/Human_DE_v1_",sampleName,".xlsx"), row.names = T)
saveRDS(Human_DE_v2, file=paste0(sampleFolder,"BBKNN/Robjects/Human_DE_v2_",sampleName,".rds"))
write.xlsx(Human_DE_v2, file =paste0(sampleFolder, "BBKNN/Human_DE_v2_",sampleName,".xlsx"), row.names = T)

## Read results
Mouse_DE_v2 <- readRDS(file=paste0(sampleFolder,"BBKNN/Robjects/Mouse_DE_v2_",sampleName,".rds"))
Human_DE_v2 <- readRDS(file=paste0(sampleFolder,"BBKNN/Robjects/Human_DE_v2_",sampleName,".rds"))

Mouse_DE_v2$Name<-rownames(Mouse_DE_v2)
Human_DE_v2$Name<-rownames(Human_DE_v2)

# Filter Rps/Rpl genes
Mouse_DE_v2_fil<-Mouse_DE_v2[!grepl("Rps|Rpl", rownames(Mouse_DE_v2)),]
Human_DE_v2_fil<-Human_DE_v2[!grepl("Rps|Rpl", rownames(Human_DE_v2)),]

######
## V1
#Look at top 150 DE genes each
intersect(head(Mouse_DE_v2_fil$Name,150),head(Human_DE_v2_fil$Name,150))
Mouse_top_DE_genes<-head(Mouse_DE_v2_fil$Name,150)
Human_top_DE_genes<-head(Human_DE_v2_fil$Name,150)

#Filter on top 150 and append names to combine tables
Mouse_DE_v2_fil_DE<-Mouse_DE_v2_fil[unique(c(Mouse_top_DE_genes,Human_top_DE_genes)),]
Human_DE_v2_fil_DE<-Human_DE_v2_fil[unique(c(Mouse_top_DE_genes,Human_top_DE_genes)),]
colnames(Mouse_DE_v2_fil_DE)<-paste0(colnames(Mouse_DE_v2_fil_DE),"_Mouse")
colnames(Human_DE_v2_fil_DE)<-paste0(colnames(Human_DE_v2_fil_DE),"_Human")
Scatterplot_table<-cbind(Mouse_DE_v2_fil_DE,Human_DE_v2_fil_DE)

#Get colors for plot depending on adj p-value and logFC
Scatterplot_table$Color<-"Black"

for (row in 1:nrow(Scatterplot_table)){
  if(Scatterplot_table[row,"p_val_adj_Mouse"] < 0.05 && Scatterplot_table[row,"p_val_adj_Human"] < 0.05){
    if(Scatterplot_table[row,"avg_log2FC_Mouse"] < 0 && Scatterplot_table[row,"avg_log2FC_Human"] < 0){
      Scatterplot_table[row,"Color"]<-"Blue"
    } else if(Scatterplot_table[row,"avg_log2FC_Mouse"] > 0 && Scatterplot_table[row,"avg_log2FC_Human"] > 0){
      Scatterplot_table[row,"Color"]<-"Red"
    } else {
      Scatterplot_table[row,"Color"]<-"Pink"
    }
  }
}

#Adjust p-value depending on UP or DOWN regulated
for (row in 1:nrow(Scatterplot_table)){
  Scatterplot_table[row,"p_val_adj_Mouse"]<-(-log(Scatterplot_table[row,"p_val_adj_Mouse"]))
  Scatterplot_table[row,"p_val_adj_Human"]<-(-log(Scatterplot_table[row,"p_val_adj_Human"]))
  if(Scatterplot_table[row,"avg_log2FC_Mouse"] < 0){
    Scatterplot_table[row,"p_val_adj_Mouse"]<-(Scatterplot_table[row,"p_val_adj_Mouse"] * -1)
  }
  if(Scatterplot_table[row,"avg_log2FC_Human"] < 0){
    Scatterplot_table[row,"p_val_adj_Human"]<-(Scatterplot_table[row,"p_val_adj_Human"] * -1)
  }
}

library(ggrepel)
p1<-ggplot(Scatterplot_table, aes(x=p_val_adj_Mouse, y=p_val_adj_Human, label = Name_Mouse, color = Color, size = 0.5)) + 
  geom_point(size=0.5)  + geom_text_repel(max.overlaps = 50, size = 2.5, vjust = 1) +
  # geom_label_repel(size = 0.5,box.padding = unit(0.5, "lines")) +
  scale_color_manual(values = c("Black" = "black",
                                "Pink"="hotpink",
                                "Red"="red")) +
  theme(legend.position="none")

pdf(file=paste0(sampleFolder,"BBKNN/6_Scatterplot_Kia_for_Stalk_FBs_",sampleName,".pdf"), width = 25, height = 20)
p1
dev.off()

saveRDS(Scatterplot_table, file=paste0(sampleFolder,"BBKNN/Robjects/Scatterplot_table_",sampleName,".rds"))
write.xlsx(Scatterplot_table, file =paste0(sampleFolder, "BBKNN/Scatterplot_table_",sampleName,".xlsx"), row.names = T)

table(Scatterplot_table$Color)
# Black  Pink   Red 
# 219    15    48 

#######
## V2
#Look at all significant DE genes
intersect(Mouse_DE_v2_fil[which(Mouse_DE_v2_fil$p_val_adj < 0.05),"Name"],
          Human_DE_v2_fil[which(Human_DE_v2_fil$p_val_adj < 0.05),"Name"])
Mouse_DE_genes<-Mouse_DE_v2_fil[which(Mouse_DE_v2_fil$p_val_adj < 0.05),"Name"]
Human_DE_genes<-Human_DE_v2_fil[which(Human_DE_v2_fil$p_val_adj < 0.05),"Name"]

#Filter on top 150 and append names to combine tables
Mouse_DE_v2_fil_DE<-Mouse_DE_v2_fil[unique(c(Mouse_DE_genes,Human_DE_genes)),]
Human_DE_v2_fil_DE<-Human_DE_v2_fil[unique(c(Mouse_DE_genes,Human_DE_genes)),]
colnames(Mouse_DE_v2_fil_DE)<-paste0(colnames(Mouse_DE_v2_fil_DE),"_Mouse")
colnames(Human_DE_v2_fil_DE)<-paste0(colnames(Human_DE_v2_fil_DE),"_Human")
Scatterplot_table<-cbind(Mouse_DE_v2_fil_DE,Human_DE_v2_fil_DE)

#Get colors for plot depending on adj p-value and logFC
Scatterplot_table$Color<-"Black"

for (row in 1:nrow(Scatterplot_table)){
  if(Scatterplot_table[row,"p_val_adj_Mouse"] < 0.05 && Scatterplot_table[row,"p_val_adj_Human"] < 0.05){
    if(Scatterplot_table[row,"avg_log2FC_Mouse"] < 0 && Scatterplot_table[row,"avg_log2FC_Human"] < 0){
      Scatterplot_table[row,"Color"]<-"Blue"
    } else if(Scatterplot_table[row,"avg_log2FC_Mouse"] > 0 && Scatterplot_table[row,"avg_log2FC_Human"] > 0){
      Scatterplot_table[row,"Color"]<-"Red"
    } else {
      Scatterplot_table[row,"Color"]<-"Pink"
    }
  }
}

#Adjust p-value depending on UP or DOWN regulated
for (row in 1:nrow(Scatterplot_table)){
  Scatterplot_table[row,"p_val_adj_Mouse"]<-(-log(Scatterplot_table[row,"p_val_adj_Mouse"]))
  Scatterplot_table[row,"p_val_adj_Human"]<-(-log(Scatterplot_table[row,"p_val_adj_Human"]))
  if(Scatterplot_table[row,"avg_log2FC_Mouse"] < 0){
    Scatterplot_table[row,"p_val_adj_Mouse"]<-(Scatterplot_table[row,"p_val_adj_Mouse"] * -1)
  }
  if(Scatterplot_table[row,"avg_log2FC_Human"] < 0){
    Scatterplot_table[row,"p_val_adj_Human"]<-(Scatterplot_table[row,"p_val_adj_Human"] * -1)
  }
}

library(ggrepel)
p1<-ggplot(Scatterplot_table, aes(x=p_val_adj_Mouse, y=p_val_adj_Human, label = Name_Mouse, color = Color, size = 0.5)) + 
  geom_point(size=0.5)  + geom_text_repel(max.overlaps = 50, size = 2.5, vjust = 1) +
  # geom_label_repel(size = 0.5,box.padding = unit(0.5, "lines")) +
  scale_color_manual(values = c("Black" = "black",
                                "Pink"="hotpink",
                                "Red"="red",
                                "Blue" = "blue")) +
  theme(legend.position="none")

pdf(file=paste0(sampleFolder,"BBKNN/6_Scatterplot_Kia_for_Stalk_FBs_All_",sampleName,".pdf"), width = 25, height = 20)
p1
dev.off()

saveRDS(Scatterplot_table, file=paste0(sampleFolder,"BBKNN/Robjects/Scatterplot_table_all_",sampleName,".rds"))
write.xlsx(Scatterplot_table, file =paste0(sampleFolder, "BBKNN/Scatterplot_table_all_",sampleName,".xlsx"), row.names = T)

table(Scatterplot_table$Color)
# Black  Blue  Pink   Red 
# 928     2    32    65 

####################

## Extra 10/03/22
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=2)

markers.to.plot_paper <- rownames(Stalk.markers.top)[c(1,2,4,5,6,3)]
FB.combined_regressed@active.ident <- factor(FB.combined_regressed@active.ident,
                                             levels = c("Other FBs Human 3","Other FBs Human 2","Other FBs Human 1","Hsp+ Stromal FBs Human",
                                                        "Stromal FBs Human","Stromal FBs Mouse", "Stalk FBs","ABCs"))

D5 <- DotPlot(FB.combined_regressed, features = markers.to.plot_paper, cols = color_list, assay = "RNA", split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot Conserved Markers logFC paper")

pdf(file=paste0(sampleFolder,"BBKNN/1_Dotplot_Stalk_Markers_logFC_paper_",sampleName,".pdf"), width = 10, height = 8)
D5
dev.off()

D2_v3 <- DotPlot(FB.combined_regressed, features = Top_markers_v2, assay = "RNA", cols = color_list, split.by = "Species") +
  RotatedAxis() + ggtitle("Dotplot Conserved Markers Adj PV paper")

pdf(file=paste0(sampleFolder,"BBKNN/1_Dotplot_Stalk_Markers_adjPV_paper_",sampleName,".pdf"), width = 10, height = 8)
D2_v3
dev.off()

####################
# Extra 05/05/22: Check for human markers for stalk FBs
# Get cell IDs from this object and then perform marker analysis in human object
table(FB.combined_regressed@meta.data$newClusters)
Idents(FB.combined_regressed)<-FB.combined_regressed@meta.data$newClusters

Stalk_human_cell_IDs<-rownames(FB.combined_regressed@meta.data[which(FB.combined_regressed@meta.data$newClusters == "Stalk FBs_Human"),])
# Stalk_human_cell_IDs<-CellSelector(FB.combined_regressed, ident = "Stalk FBs_Human")

##### Read human object
seuratObj <- readRDS(file=paste0("/run/media/clintdn/CN1465-DATA//VIB_G_drive/Human_snRNAseq_CP_COVID/Robjects/COVID-19_brain_snRNA-seq_choroid_plexus_final_seurat_v3.2.3.rds"))

DimPlot(seuratObj, reduction = "umap", label=T,repel = T,  pt.size = 1)

seuratObj@meta.data$newID<-as.character(seuratObj@meta.data$cellID)
seuratObj@meta.data[Stalk_human_cell_IDs,"newID"]<-"Stalk FBs"
Idents(seuratObj)<-seuratObj@meta.data$newID

DimPlot(seuratObj, reduction = "umap", label=T,repel = T,  pt.size = 1)

Markers_cellID<-FindMarkers(seuratObj, ident.1 = "Stalk FBs", ident.2 = levels(Idents(seuratObj))[-8], assay = "integrated",
                            min.diff.pct = 0.25)

Markers_cellID$cluster<-"Stalk_FBs"
Markers_cellID$score<-Markers_cellID$pct.1/(Markers_cellID$pct.2+0.01)*Markers_cellID$avg_log2FC
  
Markers_cellID<-Markers_cellID[order(Markers_cellID$score, decreasing=TRUE),]
Markers_cellID$gene<-rownames(Markers_cellID)

### Write to Excel
library('openxlsx')
write.xlsx(Markers_cellID, file =paste0(sampleFolder,"BBKNN/Stalk_FBs_human_markers_vs_all_other_cells_in_human_dataset.xlsx"))

saveRDS(Markers_cellID, file=paste0(sampleFolder,"Robjects/Markers_cellID_",sampleName,".rds"))
Markers_cellID<-readRDS(file=paste0(sampleFolder,"Robjects/Markers_cellID_",sampleName,".rds"))

##################################################################

## Update annotation Stalk -> Base for paper (30/01/23)
FB.combined_regressed@meta.data[["annotated_clusters"]]<-as.factor(FB.combined_regressed@meta.data[["annotated_clusters"]])
levels(FB.combined_regressed@meta.data[["annotated_clusters"]])[1]<-"Base Fibroblasts"
FB.combined_regressed@meta.data[["annotated_clusters_detailed"]]<-as.factor(FB.combined_regressed@meta.data[["annotated_clusters_detailed"]])
levels(FB.combined_regressed@meta.data[["annotated_clusters_detailed"]])[6]<-"Base Fibroblasts"
levels(FB.combined_regressed@meta.data[["newClustersTmp"]])[4]<-"Base FBs"
FB.combined_regressed@meta.data[["newClusters"]]<-as.factor(FB.combined_regressed@meta.data[["newClusters"]])
levels(FB.combined_regressed@meta.data[["newClusters"]])[6]<-"Base FBs_Human"
levels(FB.combined_regressed@meta.data[["newClusters"]])[7]<-"Base FBs_Mouse"
FB.combined_regressed$FB_Base_signature_score1<-FB.combined_regressed$FB_stalk_signature_score1
FB.combined_regressed$FB_stalk_signature_score1<-NULL
FB.combined_regressed@meta.data[["annotated_clusters_new"]]<-as.factor(FB.combined_regressed@meta.data[["annotated_clusters_new"]])
levels(FB.combined_regressed@meta.data[["annotated_clusters_new"]])[2]<-"Base Fibroblasts"
levels(FB.combined_regressed@meta.data[["Integrated_annotated_clusters"]])[6]<-"Base FBs"
FB.combined_regressed@meta.data[["NewClusters_combo_annotated"]]<-as.factor(FB.combined_regressed@meta.data[["NewClusters_combo_annotated"]])
levels(FB.combined_regressed@meta.data[["NewClusters_combo_annotated"]])[6]<-"Base FBs_Human"
levels(FB.combined_regressed@meta.data[["NewClusters_combo_annotated"]])[7]<-"Base FBs_Mouse"
levels(FB.combined_regressed@meta.data[["annotated_BBKNN_clusters"]])[4]<-"Base FBs Mouse and Human"
levels(FB.combined_regressed@meta.data[["annotated_BBKNN_clusters"]])[8]<-"ABCs Mouse"
levels(FB.combined_regressed@meta.data[["annotated_BBKNN_clusters_v2"]])[4]<-"Base FBs Mouse and Human"
levels(FB.combined_regressed@meta.data[["annotated_BBKNN_clusters_v2"]])[8]<-"ABCs Mouse"
Idents(FB.combined_regressed)<-FB.combined_regressed$annotated_BBKNN_clusters_v2
