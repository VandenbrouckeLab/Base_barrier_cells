## Script for merging our datasets with public datasets to investigate the origin of our FBs
## Follow-up script to process and explore the Fibroblast origin complete object in Figure 1 manuscript
## Initially the object still contains everything and subsequently it is subsetted to remove ChP Epithelial cells and Immune cells

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

sampleName <- "Merge_FB_datasets" #Change for this analysis!!!
sampleFolder<-paste0("Merge","/")

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
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

Idents(seuratObj)<-seuratObj@meta.data$harmony_clusters
seuratObj@meta.data$sliced_clusters<-seuratObj@meta.data$harmony_clusters
# Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters

## Check 1.0 clustering:
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "RNA_snn_res.1", label.size = 8)

## Cleaning
# Split up cluster 10 -> 26
umapSlice1<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 < 0.5, UMAP_2 > 4)
umapSlice2<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 < -1, UMAP_2 > 2)
wantedCells1<-intersect(umapSlice1$cell, WhichCells(seuratObj, idents = 10))
wantedCells2<-intersect(umapSlice2$cell, WhichCells(seuratObj, idents = 10))
wantedCells<-c(wantedCells1,wantedCells2)
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 26)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

## Save new clustering
seuratObj@meta.data$sliced_clusters <- seuratObj@active.ident #Sliced clustering
seuratObj@meta.data$annotated_clusters3 <- seuratObj@active.ident #Sliced clustering

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

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_merge/QC/1_UMI_",sampleName,".png"), width = 20)

########## mito.genes plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'subsets_Mito_percent',"mito")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'subsets_Mito_percent',"mito")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_merge/QC/2_percMito_",sampleName,".png"), width = 20)


########## PCA plot ##########
# pdf(file=paste0(sampleFolder,"results_merge/QC/13a_PCA_",sampleName,".pdf"), width=10)
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,3))
# dev.off()

#RNA clusters
dir.create(paste0(sampleFolder,"results_merge/Annotation"))

pdf(file=paste0(sampleFolder,"results_merge/Annotation/1_annotation_Color_RNA_clusters_on_harmony_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$harmony_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$harmony_clusters==i),])))
  C1<-C1+ggtitle(paste0("Harmony_cluster_",i))
  print(C1)
}
dev.off()

################################################################################
########## CHECK marker GENES
################################################################################
dir.create(paste0(sampleFolder,"results_merge/Feature_plots"))

##### Epithelial marker ->
Features<-c("Otx2", "Ttr")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"CPE"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Fras1 marker ->
Features<-c("Fras1", "Xist","Neat1","Steap2","Folr1","Car2")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"CPE_extra"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order =T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_ordered_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Endothelial marker ->
Features<-c("Pecam1", "Flt1","Plvap")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"EC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Vascular associated marker ->
Features<-c("Pdgfrb", "Mylk","Myh11","Tagln")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"VAC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Fibroblast marker ->
Features<-c("Dcn", "Col1a1")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"FB"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Macrophage marker ->
Features<-c("Adgre1", "Csf1r","Fcgr1")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"MF"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Microglia marker ->
Features<-c("P2ry12")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Microglia-like"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### NK cell marker ->
Features<-c("Klrb1c", "Gzmb")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"NK"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### cDC marker ->
Features<-c("Cd209a", "Ccr7","Xcr1")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"DC"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Neutrophil marker ->
Features<-c("S100a8", "Ngp","Retnlg")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"NF"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

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
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

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
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Mitotic cell marker ->
Features<-c("Birc5", "Mki67")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Mitotic_cells"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Dissociation effect ->
Features<-c("Fos", "Junb","Atf3","Dusp1","Ccl4")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Dissociation"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"),
                reduction = "tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


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
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  RNAmarkersList_RNAclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(RNAmarkersList_RNAclus)<-paste0("RNAcluster",0:totalNrRNAclusters_RNAclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_RNAclus, file =paste0(sampleFolder, "results_merge/Marker_lists/RNAmarkersList_RNAclus_",sampleName,".xlsx"))

########################################

#############################
### Extra detail clusters
############################
Idents(seuratObj)<-seuratObj@meta.data$harmony_clusters

detail1_16vs9<-FindMarkers(seuratObj, ident.1 = 16, ident.2 = 9, min.pct = 0.10,
                           min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

detail2_1vs9<-FindMarkers(seuratObj, ident.1 = 1, ident.2 = 9, min.pct = 0.10,
                          min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

detail3_8vs1.16.9<-FindMarkers(seuratObj, ident.1 = 8, ident.2 = c(1,9,16), min.pct = 0.10,
                               min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)


##### Create list
listDEgenesExtra<-tibble::lst(detail1_16vs9, detail2_1vs9, detail3_8vs1.16.9)

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

saveRDS(listDEgenesExtra,file=paste0(sampleFolder,"results_merge/Robjects/detailClusters_sliced_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDEgenesExtra, paste0(sampleFolder,"results_merge/Marker_lists/detailClusters_sliced_",sampleName,".xlsx"))

########################################################################################################################

## Check other annotation
DimPlot(seuratObj, reduction = "umap", label = T, group.by = "ClusterName", label.size = 4)

DimPlot(seuratObj, reduction = "umap", label = T, group.by = "Clusters", label.size = 4)

pdf(file=paste0(sampleFolder,"results_merge/Annotation/1_UMAP_split_datasets_",sampleName,".pdf"), width = 40, height = 30)
DimPlot(seuratObj, reduction = "umap", label = T, split.by = "orig.ident", group.by = "orig.ident", label.size = 4, ncol = 4)
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

### Transfer old annotation ###

## Update clusters URVB
seuratObj_agg15 <- readRDS(file="../RVD_aggr15/Robjects/seuratObj_sliced_RVD_aggr15_Harmony.rds")
colnames(seuratObj_agg15)
length(colnames(seuratObj[,grep("URVB",colnames(seuratObj))]))

URVB_cells<-intersect(paste0("URVB_",colnames(seuratObj_agg15)),colnames(seuratObj[,grep("URVB",colnames(seuratObj))])) # Different filtering!
URVB_cells_old<-gsub("URVB_","",URVB_cells)
seuratObj@meta.data[URVB_cells,"New_clusters"]<-as.character(seuratObj_agg15@meta.data[URVB_cells_old,"annotated_clusters"])

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

DimPlot(seuratObj, reduction = "umap", label = T, group.by = "New_clusters", label.size = 2)

## Update EP clusters Shah (no original annotation)
EP1_cells<-colnames(seuratObj[,grep("EP1",colnames(seuratObj))])
EP2_cells<-colnames(seuratObj[,grep("EP2",colnames(seuratObj))])
EP3_cells<-colnames(seuratObj[,grep("EP3",colnames(seuratObj))])

colorSomeCells(clusterMatrix, umapTable, EP1_cells)
colorSomeCells(clusterMatrix, umapTable, EP2_cells)
colorSomeCells(clusterMatrix, umapTable, EP3_cells)

umapSlice1<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_2 > 0)
umapSlice2<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_2 < 0)

## Annotate based on expression
seuratObj@meta.data[intersect(c(EP1_cells,EP2_cells,EP3_cells), WhichCells(seuratObj, idents = 19)),"New_clusters"]<-"Neuroblasts"
seuratObj@meta.data[intersect(c(EP1_cells,EP2_cells,EP3_cells), WhichCells(seuratObj, idents = 8)),"New_clusters"]<-"aNSCs-late"
seuratObj@meta.data[intersect(c(EP1_cells,EP2_cells,EP3_cells), WhichCells(seuratObj, idents = 21)),"New_clusters"]<-"Ependymal_cells"
seuratObj@meta.data[intersect(umapSlice1$cell, intersect(c(EP1_cells,EP2_cells,EP3_cells), WhichCells(seuratObj, idents = 14))),"New_clusters"]<-"Ependymal_cells"
seuratObj@meta.data[intersect(umapSlice2$cell, intersect(c(EP1_cells,EP2_cells,EP3_cells), WhichCells(seuratObj, idents = 14))),"New_clusters"]<-"NSCs"
seuratObj@meta.data[intersect(c(EP1_cells,EP2_cells,EP3_cells), WhichCells(seuratObj, idents = 25)),"New_clusters"]<-"OL/OPCs"
seuratObj@meta.data[intersect(c(EP1_cells,EP2_cells,EP3_cells), WhichCells(seuratObj, idents = c(2,6))),"New_clusters"]<-"Undefined1"
seuratObj@meta.data[intersect(c(EP1_cells,EP2_cells,EP3_cells), WhichCells(seuratObj, idents = 7)),"New_clusters"]<-"Undefined2"
seuratObj@meta.data[intersect(c(EP1_cells,EP2_cells,EP3_cells), WhichCells(seuratObj, idents = c(5,24))),"New_clusters"]<-"Microglia"

rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Unknown"),])

# Leftover unknown:
# 300 cells from Ependymal paper which cluster separately
# 1200 cells which aren't annotated in aggr15 URVB (filtered there, but not here!!)
colorSomeCells(clusterMatrix,umapTable,rownames(seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "Unknown"),]))

# Cells from Vanlandewijck paper not fully annotated (no metadata)
# EC cells above, mural cells mid, part of Pdgfra cells at bottom, but also some oligodendrocytes middle!!
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "GSE98816_Vanlandewijck"))

U_old_annot<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "New_clusters", label.size = 3)
ggsave(U_old_annot, file=paste0(sampleFolder,"results_merge/Annotation/1_UMAP_old_annot_",sampleName,".png"), height = 8, width = 15, dpi = "retina")


########################################################################################################################

########################################
##### RNA clusters post-slice (after split 10!!)
########################################
########################################

seuratObj@meta.data$sliced_clusters<- factor(seuratObj@meta.data$sliced_clusters,sort(as.numeric(levels(seuratObj@meta.data$sliced_clusters)))) #reorder levels
seuratObj@meta.data$annotated_clusters3 <- factor(seuratObj@meta.data$annotated_clusters3,sort(as.numeric(levels(seuratObj@meta.data$annotated_clusters3)))) #reorder levels
Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters

U_annot<-DimPlot(seuratObj, reduction = "umap", label = T, group.by = "sliced_clusters", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_merge/Annotation/5_UMAP_sliced_clusters_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

########################################
##### Markers annotated clusters
########################################
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident
levels(seuratObj@meta.data$annotated_clusters) <- c("Epithelial_cells","Fibroblasts","Endothelial cells","Epithelial_cells","Epithelial_cells","Macrophages",
                                                    "Endothelial cells","VAC","Proliferating cells","Fibroblasts","VECC","VSMCA","Epithelial_cells","Endothelial cells",
                                                    rep("Neuronal and glial cells",2),"Fibroblasts","Epithelial_cells","Other immune cells","Neuronal and glial cells",
                                                    "Epithelial_cells","Neuronal and glial cells","Epithelial_cells","Epithelial_cells",
                                                    "Microglia","Neuronal and glial cells")
U_annot<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_merge/Annotation/2_UMAP_annotated1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

seuratObj@meta.data$annotated_clusters2 <- as.character(seuratObj@meta.data$annotated_clusters)
seuratObj@meta.data[which(seuratObj@meta.data$annotated_clusters2 == "Proliferating cells"),"annotated_clusters2"]<-"Proliferating Fibroblasts"
seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "aNSCs-late"),"annotated_clusters2"]<-"Proliferating NSCs"
seuratObj@meta.data$annotated_clusters2 <-as.factor(seuratObj@meta.data$annotated_clusters2)

U_annot2<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters2", label.size = 4)
ggsave(U_annot2, file=paste0(sampleFolder,"results_merge/Annotation/2_UMAP_annotated2_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

U_origin2<-DimPlot(seuratObj, reduction = "umap", label = F, split.by ="orig.ident", group.by = "annotated_clusters2", label.size = 3, cols = Colorset, ncol = 4)
ggsave(U_origin2, file=paste0(sampleFolder,"results_merge/Annotation/1_UMAP_dataset_origin2_",sampleName,".png"), height = 20, width = 25, dpi = 300)

Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters2

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
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  RNAmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
# names(RNAmarkersList_SCTclus)<-paste0("SCTclustersliced",0:totalNrRNAclusters_SCTclus)

### Write to Excel
library('openxlsx')
write.xlsx(RNAmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge/Marker_lists/RNAmarkersList_SCTclus_",sampleName,"_annotated.xlsx"))

######################################################

### Make heatmap for annotated clusters
RNAMarkers_SCTclus<- readRDS(file=paste0(sampleFolder,"results_merge/Robjects/RNAmarkersList_SCTclus_",sampleName,"_annotated.rds"))

## Perform on a subset -> better view of smaller clusters!!
seuratObj.small <- subset(seuratObj, downsample = 500)

########## Get HVG ##########
seuratObj.small <- FindVariableFeatures(object = seuratObj.small, selection.method = "vst", nfeatures = nrow(seuratObj.small@assays$RNA))
length(VariableFeatures(seuratObj.small))

########## Scale ##########
seuratObj.small <- ScaleData(seuratObj.small)

## Heatmap RNA markers on RNA clusters
top10 <- RNAMarkers_SCTclus %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
D1<-DoHeatmap(seuratObj.small, features = top10$gene, group.by = "annotated_clusters2") + NoLegend()
ggsave(D1, file=paste0(sampleFolder, "results_merge/Heatmaps/Heatmap_RNAmarkersList_Annotatedclus_",sampleName,"_2.png"), 
       height = 20, width = 12, dpi = "retina")

pdf(file=paste0(sampleFolder, "results_merge/Heatmaps/Heatmap_RNAmarkersList_Annotatedclus_",sampleName,"_2.pdf"), 
    height = 25, width = 25)
DoHeatmap(seuratObj.small, features = top10$gene, group.by = "annotated_clusters2") + NoLegend()
dev.off()


######################################################

# Frequency tables (sliced)
Sample <- seuratObj@meta.data$orig.ident
cluster <- seuratObj@meta.data$annotated_clusters2
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

## Extra vis

seuratObj@meta.data$sample_origin<-as.factor(seuratObj@meta.data$orig.ident)
levels(seuratObj@meta.data$sample_origin)<-c("FB_DeSisto_et_al_1","FB_DeSisto_et_al_2","FB_DeSisto_et_al_3",
                                             "Ependymal_Zeisel_et_al","Vascular_Vanlandewijck_et_al",
                                             "Ependymal_Shah_et_al_1","Ependymal_Shah_et_al_2","Ependymal_Shah_et_al_3",
                                             "CP_LpsNeg_4V_Urvb","CP_LpsNeg_LV_Urvb","CP_Young_4V_Urvb","CP_Young_LV_Urvb",
                                             "Vascular_Zeisel_et_al")

library(RColorBrewer)
Colorset<-c(brewer.pal(12,"Set3"),"magenta")

U_origin<-DimPlot(seuratObj, reduction = "umap", label = F, group.by = "sample_origin", label.size = 3, cols = Colorset)
ggsave(U_origin, file=paste0(sampleFolder,"results_merge/Annotation/4_UMAP_dataset_origin1_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

U_origin2<-DimPlot(seuratObj, reduction = "umap", label = F, split.by ="sample_origin", group.by = "annotated_clusters2", label.size = 3, ncol = 4)
ggsave(U_origin2, file=paste0(sampleFolder,"results_merge/Annotation/4_UMAP_dataset_origin2_",sampleName,".png"), height = 20, width = 25, dpi = 300)

U_origin3<-DimPlot(seuratObj, reduction = "umap", label = F, split.by ="sample_origin", group.by = "sample_origin", label.size = 3, cols = Colorset, ncol = 4)
ggsave(U_origin3, file=paste0(sampleFolder,"results_merge/Annotation/4_UMAP_dataset_origin3_",sampleName,".png"), height = 20, width = 25, dpi = 300)

pdf(file = paste0(sampleFolder,"results_merge/Annotation/4_UMAP_dataset_origin2_",sampleName,".pdf"), width = 30, height = 20)
U_origin2
dev.off()

pdf(file = paste0(sampleFolder,"results_merge/Annotation/4_UMAP_dataset_origin3_",sampleName,".pdf"), width = 30, height = 20)
U_origin3
dev.off()

#############################################################################################
#############################################################################################

## Split off CPE and Imm (not relevant for paper)
## Also clean up EC <-> CPE (Cluster 10!!!)
seuratObj@meta.data$annotated_clusters3 <- seuratObj@active.ident
levels(seuratObj@meta.data$annotated_clusters3) <- c("Epithelial_cells","Fibroblasts","Endothelial cells","Epithelial_cells","Epithelial_cells","Macrophages",
                                                     "Endothelial cells","VAC","Proliferating cells","Fibroblasts","VECC","VSMCA","Epithelial_cells","Endothelial cells",
                                                     rep("Neuronal and glial cells",2),"Fibroblasts","Epithelial_cells","Other immune cells","Neuronal and glial cells",
                                                     "Epithelial_cells","Neuronal and glial cells","Epithelial_cells","Epithelial_cells",
                                                     "Microglia","Neuronal and glial cells","Epithelial_cells")
seuratObj@meta.data$annotated_clusters3<-as.character(seuratObj@meta.data$annotated_clusters3)
seuratObj@meta.data[which(seuratObj@meta.data$annotated_clusters3 == "Proliferating cells"),"annotated_clusters3"]<-"Proliferating Fibroblasts"
seuratObj@meta.data[which(seuratObj@meta.data$New_clusters == "aNSCs-late"),"annotated_clusters3"]<-"Proliferating NSCs"
seuratObj@meta.data$annotated_clusters3 <-as.factor(seuratObj@meta.data$annotated_clusters3)

U_annot3<-DimPlot(seuratObj, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters3", label.size = 4)
ggsave(U_annot3, file=paste0(sampleFolder,"results_merge/Annotation/5_UMAP_annotated3_",sampleName,".png"), height = 10, width = 15, dpi = "retina")

# Subset object
Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters3
seuratObjNew<-subset(seuratObj, idents=c("Fibroblasts","Endothelial cells","Neuronal and glial cells","Proliferating Fibroblasts",
                                         "Proliferating NSCs","VAC","VECC","VSMCA"))

DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters3", label.size = 4)

# Clean outliers!!!!
Idents(seuratObjNew)<-as.character(seuratObjNew@meta.data$annotated_clusters3)

U1 <- DimPlot(seuratObjNew, reduction = "umap", label = T, label.size = 4)
seuratObjNew <- CellSelector(U1, object=seuratObjNew, ident="Outliers1")

U1 <- DimPlot(seuratObjNew, reduction = "umap", label = T, label.size = 4)
seuratObjNew <- CellSelector(U1, object=seuratObjNew, ident="Outliers2")

seuratObjNew<-subset(seuratObjNew, idents=c("Fibroblasts","Endothelial cells","Neuronal and glial cells","Proliferating Fibroblasts",
                                            "Proliferating NSCs","VAC","VECC","VSMCA"))

######
library(RColorBrewer)
Colorset<-c(brewer.pal(12,"Set3"),"magenta")

U_old_annot_sub<-DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "New_clusters", label.size = 3)

pdf(file=paste0(sampleFolder,"results_merge/Annotation/6_UMAP_old_annot_",sampleName,".pdf"), height = 10, width = 13)
U_old_annot_sub
dev.off()

U_annot_sub<-DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters3", label.size = 4)
ggsave(U_annot_sub, file=paste0(sampleFolder,"results_merge/Annotation/6_UMAP_annotated_",sampleName,"_no_CPE_or_Imm.png"), height = 10, width = 9, dpi = "retina")

pdf(file = paste0(sampleFolder,"results_merge/Annotation/6_UMAP_annotated_",sampleName,"_no_CPE_or_Imm.pdf"), width = 9, height = 10)
U_annot_sub
dev.off()

U_origin_sub<-DimPlot(seuratObjNew, reduction = "umap", label = F, group.by = "sample_origin", label.size = 3, cols = Colorset)
ggsave(U_origin_sub, file=paste0(sampleFolder,"results_merge/Annotation/6_UMAP_dataset_origin1_",sampleName,"_no_CPE_or_Imm.png"), height = 10, width = 10, dpi = "retina")

U_origin2_sub<-DimPlot(seuratObjNew, reduction = "umap", label = F, split.by ="sample_origin", group.by = "annotated_clusters3", label.size = 3, ncol = 4)
ggsave(U_origin2_sub, file=paste0(sampleFolder,"results_merge/Annotation/6_UMAP_dataset_origin2_",sampleName,"_no_CPE_or_Imm.png"), height = 20, width = 16, dpi = 300)

U_origin3_sub<-DimPlot(seuratObjNew, reduction = "umap", label = F, split.by ="sample_origin", group.by = "sample_origin", label.size = 3, cols = Colorset, ncol = 4)
ggsave(U_origin3_sub, file=paste0(sampleFolder,"results_merge/Annotation/6_UMAP_dataset_origin3_",sampleName,"_no_CPE_or_Imm.png"), height = 20, width = 16, dpi = 300)

pdf(file = paste0(sampleFolder,"results_merge/Annotation/6_UMAP_dataset_origin2_",sampleName,"_no_CPE_or_Imm.pdf"), width = 16, height = 20)
U_origin2_sub
dev.off()

pdf(file = paste0(sampleFolder,"results_merge/Annotation/6_UMAP_dataset_origin3_",sampleName,"_no_CPE_or_Imm.pdf"), width = 16, height = 20)
U_origin3_sub
dev.off()

###########################

## Extra figure 23/03/21: Group datasets by paper as with subset!!(13 -> 6)

## New annotation: combine dataset names
seuratObjNew@meta.data$sample_origin2<-seuratObjNew@meta.data$sample_origin
levels(seuratObjNew@meta.data$sample_origin2)<-c("FB_DeSisto_et_al","FB_DeSisto_et_al",         
                                                 "FB_DeSisto_et_al","Ependymal_Zeisel_et_al" ,"Vascular_Vanlandewijck_et_al","Ependymal_Shah_et_al",      
                                                 "Ependymal_Shah_et_al","Ependymal_Shah_et_al","CP_Urvb","CP_Urvb","CP_Urvb","CP_Urvb","Vascular_Zeisel_et_al")

library(RColorBrewer)
Colorset<-c(brewer.pal(8,"Set1")[c(1,8,6,2,3,4)])

U_origin2<-DimPlot(seuratObjNew, reduction = "umap", label = F, group.by = "sample_origin2", label.size = 3, cols = Colorset)
ggsave(U_origin2, file=paste0(sampleFolder,"results_merge/Annotation/7_UMAP_dataset_origin1_",sampleName,"_no_CPE_or_Imm.png"), height = 10, width = 10, dpi = "retina")

pdf(file = paste0(sampleFolder,"results_merge/Annotation/7_UMAP_dataset_origin1_",sampleName,"_no_CPE_or_Imm.pdf"), width = 10, height = 10)
U_origin2
dev.off()

###########################

## Final change to dataset names (2 options -> 2 metadata columns) 26/03/21
## New annotation: combine dataset names
seuratObjNew@meta.data$sample_origin2_optionA<-seuratObjNew@meta.data$sample_origin2
seuratObjNew@meta.data$sample_origin2_optionB<-seuratObjNew@meta.data$sample_origin2
levels(seuratObjNew@meta.data$sample_origin2_optionA)<-c("FB_DeSisto_et_al","Ependymal_Zeisel_et_al","Vascular_Vanlandewijck_et_al",
                                                            "Ependymal_Shah_et_al","CP_Verhaege_et_al","Vascular_Zeisel_et_al")
levels(seuratObjNew@meta.data$sample_origin2_optionB)<-c("FB_DeSisto_et_al","Ependymal_Zeisel_et_al","Vascular_Vanlandewijck_et_al",
                                                            "Ependymal_Shah_et_al","Choroid_Plexus_cells","Vascular_Zeisel_et_al")

library(RColorBrewer)
Colorset<-c(brewer.pal(8,"Set1")[c(1,8,6,2,3,4)])

U_origin2A<-DimPlot(seuratObjNew, reduction = "umap", label = F, group.by = "sample_origin2_optionA", label.size = 3, cols = Colorset)
ggsave(U_origin2A, file=paste0(sampleFolder,"results_merge/Annotation/8_UMAP_dataset_origin1_optionA_",sampleName,"_no_CPE_or_Imm.png"), height = 10, width = 10, dpi = "retina")

U_origin2B<-DimPlot(seuratObjNew, reduction = "umap", label = F, group.by = "sample_origin2_optionB", label.size = 3, cols = Colorset)
ggsave(U_origin2B, file=paste0(sampleFolder,"results_merge/Annotation/8_UMAP_dataset_origin1_optionB_",sampleName,"_no_CPE_or_Imm.png"), height = 10, width = 10, dpi = "retina")

###########################

## Extra plot for Daan (17/12/21)
# Blend featureplot
F1 <- FeaturePlot(seuratObjNew, features = c("Cldn11","Igfbp6"), blend = TRUE, blend.threshold = 0.5) #Default threshold 0.5
pdf(file = paste0(sampleFolder,"results_merge/Annotation/9_UMAP_Featureplot_blend_Cldn11_Igfbp6_",sampleName,"_no_CPE_or_Imm.pdf"), width = 20, height = 8)
F1
dev.off()

# Module score Seurat
library(viridis)

FB_stalk_signature <- list(c("Cldn11","Igfbp6"))

seuratObjNew <- AddModuleScore(object = seuratObjNew, features = FB_stalk_signature, name = "FB_stalk_signature_score")
F2 <- FeaturePlot(object = seuratObjNew, features = "FB_stalk_signature_score1", order = T)  + scale_color_viridis(option = "C")
pdf(file = paste0(sampleFolder,"results_merge/Annotation/9_UMAP_Featureplot_modulescore_Cldn11_Igfbp6_",sampleName,"_no_CPE_or_Imm.pdf"), width = 8, height = 10)
F2
dev.off()

Idents(seuratObjNew)<-seuratObjNew@meta.data$New_clusters
F2.5 <- FeaturePlot(object = seuratObjNew, features = "FB_stalk_signature_score1", order = T, label = T, repel = T, label.size = 3, cols = c("Yellow","Red"))
pdf(file = paste0(sampleFolder,"results_merge/Annotation/9_UMAP_Featureplot_modulescore_Cldn11_Igfbp6_labeled_",sampleName,"_no_CPE_or_Imm.pdf"), width = 8, height = 10)
F2.5
dev.off()

# Nebulosa code (in R v4)
seuratObjNew<-UpdateSeuratObject(seuratObjNew)

F3<-plot_density(seuratObjNew, features =  c("Cldn11","Igfbp6"), slot = "data", reduction = "umap",joint = TRUE)
pdf(file = paste0(sampleFolder,"results_merge/Annotation/9_UMAP_Nebulosa_Cldn11_Igfbp6_",sampleName,"_no_CPE_or_Imm.pdf"), width = 20, height = 8)
F3
dev.off()

###########################

##### Read object
seuratObjNew <- readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_",sampleName,"_harmony_RNA_no_CPE_or_Imm.rds"))

##### Save object
saveRDS(seuratObjNew, file=paste0(sampleFolder,"Robjects/seuratObj_",sampleName,"_harmony_RNA_no_CPE_or_Imm.rds"))

#############################################################################################
#############################################################################################

## Extra Daan 12/07/21
## Check markers specific for Stalk FBs compared to other brain cells

DimPlot(seuratObjNew, reduction = "umap", label = T, repel = T, group.by = "New_clusters", label.size = 3)

Idents(seuratObjNew)<-seuratObjNew@meta.data$New_clusters
levels(Idents(seuratObjNew))
detail_Stalk_FBs_strict<-FindMarkers(seuratObjNew, ident.1 = levels(Idents(seuratObjNew))[41], ident.2 = levels(Idents(seuratObjNew))[-41], 
                           min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

detail_Stalk_FBs_super_strict<-FindMarkers(seuratObjNew, ident.1 = levels(Idents(seuratObjNew))[41], ident.2 = levels(Idents(seuratObjNew))[-41], 
                                     min.pct = 0.50, min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

detail_Stalk_FBs_lenient<-FindMarkers(seuratObjNew, ident.1 = levels(Idents(seuratObjNew))[41], ident.2 = levels(Idents(seuratObjNew))[-41], 
                                           min.pct = 0.10, logfc.threshold = 0.30, only.pos = FALSE)


##### Create list
listDEgenesExtra<-tibble::lst(detail_Stalk_FBs_strict, detail_Stalk_FBs_super_strict, detail_Stalk_FBs_lenient)

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

saveRDS(listDEgenesExtra,file=paste0(sampleFolder,"results_merge/Robjects/detail_Stalk_FBs_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDEgenesExtra, paste0(sampleFolder,"results_merge/Marker_lists/detail_Stalk_FBs_",sampleName,".xlsx"))

## Featureplot super strict markers
as.character(listDEgenesExtra$detail_Stalk_FBs_super_strict$gene[1:9])

##### Stalk FBs
Features<-c(as.character(listDEgenesExtra$detail_Stalk_FBs_super_strict$gene[1:9]))
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Stalk_FBs_super_Strict"
F1<-FeaturePlot(object = seuratObjNew, features = Features, cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order =T)
ggsave(F1, file=paste0(sampleFolder,"results_merge/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

########

## Extra october 2021
## ABC markers
detail_ABCs_lenient<-FindMarkers(seuratObjNew, ident.1 = levels(Idents(seuratObjNew))[25], ident.2 = levels(Idents(seuratObjNew))[-25], 
                                 min.pct = 0.10, logfc.threshold = 0.30, only.pos = FALSE)

detail_ABCs_strict<-FindMarkers(seuratObjNew, ident.1 = levels(Idents(seuratObjNew))[25], ident.2 = levels(Idents(seuratObjNew))[-25], 
                                min.pct = 0.10, min.diff.pct = 0.25, logfc.threshold = 0.30, only.pos = FALSE)

##### Create list
listDEgenesABCs<-tibble::lst(detail_ABCs_lenient, detail_ABCs_strict)

##Add geneSymbol in column (for the export)
listDEgenesABCs<-lapply(listDEgenesABCs,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenesABCs<-lapply(listDEgenesABCs, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDEgenesABCs<-lapply(listDEgenesABCs, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                           x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
# listDEgenesABCs<-lapply(listDEgenesABCs, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
listDEgenesABCs<-lapply(listDEgenesABCs,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDEgenesABCs,file=paste0(sampleFolder,"results_merge/Robjects/detail_ABCs_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDEgenesABCs, paste0(sampleFolder,"results_merge/Marker_lists/detail_ABCs_",sampleName,".xlsx"))
