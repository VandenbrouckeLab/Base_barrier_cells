## Integrate snRNAseq Human Covid dataset and snRNAseq Lehtinen lab
## Convert Human gene symbols to mouse gene symbols
## Then perform Harmony as first check. Perhaps try integration afterwards

################################################################################
########## GENERAL
################################################################################

###### Load packages ######
library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')

source('/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/RAW_DATA/script_functions_COVID.R') #KEVIN

########################################
##### Getwd
########################################

setwd("/run/media/clintdn/CN1465-DATA/VIB_G_drive/")

sampleName_Human<-"Human_snRNAseq_CP_COVID"
sampleFolder<-paste0(sampleName_Human,"/")

##### Read human object
seuratObj_Human <- readRDS(file=paste0(sampleFolder,"Robjects/COVID-19_brain_snRNA-seq_choroid_plexus_final_seurat_v3.2.3.rds"))

DimPlot(seuratObj_Human, reduction = "umap", label=T,repel = T, group.by="cellID", pt.size = 1)

##### Read mouse object Lehtinen mesenchymal cells all ages
seuratObj_Mouse <- readRDS(file=paste0("/home/clintdn/VIB/DATA/Roos/Daan 1/FB_datasets/Merge_Lehtinen_extra/Robjects/seuratObj_Lehtinen_snRNAseq_datasets_harmony_RNA.rds"))

seuratObj_Mouse <- UpdateSeuratObject(seuratObj_Mouse)
DimPlot(seuratObj_Mouse, reduction = "umap", label=T,repel = T, group.by="annotated_clusters", pt.size = 1)

##### Subset object to mesenchymal cells
seuratObj_Human_subset<-subset(seuratObj_Human, idents = "Mesenchymal")

##### Extract count matrix and metadata
Counts_human_mesenchymal_data <- GetAssayData(seuratObj_Human_subset[["SCT"]], slot = "counts") #RNA counts slot doesn't seem correct!!!!!!
Metadata_human_mesenchymal_data <- seuratObj_Human_subset@meta.data

Original_human_gene_symbols <- rownames(Counts_human_mesenchymal_data)

## Convert rownames 
H.genes <- Original_human_gene_symbols

# Convert with BioMart
# Basic function to convert human to mouse gene names 
## Update: changed useMart to useEnsembl!!!!!!!!! (14/10/20)
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  # humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  # print(head(humanx))
  return(genesV2)
}

M.genes <- convertHumanGeneList(H.genes) #humanx
M.genes_df <- convertHumanGeneList(H.genes) #genesV2

## Recreate seuratobject (remove genes which didn't get converted??)
Counts_human_mesenchymal_data_subset_mouse<-Counts_human_mesenchymal_data[M.genes_df$HGNC.symbol,]
rownames(Counts_human_mesenchymal_data_subset_mouse)<-M.genes_df$MGI.symbol

seuratObj_Human_subset_converted <- CreateSeuratObject(counts = Counts_human_mesenchymal_data_subset_mouse, project = "Human_CP_converted", assay = "RNA", min.cells = 3, min.features = 200) #Stick data in RNA slot

## Add metaData back 
seuratObj_Human_subset_converted@meta.data <-cbind(seuratObj_Human_subset_converted@meta.data, Metadata_human_mesenchymal_data)

# Remove duplicate columns (slight diff due to removal certain genes. Remove the old ones!)
seuratObj_Human_subset_converted@meta.data[c(4,5,6)]<-NULL

## Normalize
seuratObj_Human_subset_converted <- NormalizeData(object = seuratObj_Human_subset_converted, normalization.method = "LogNormalize", scale.factor = 10000)

## Evolution number of genes
length(H.genes)
length(M.genes_df$MGI.symbol)
length(M.genes)

## Subset rows of both seuratobjects???
## Convert certain metadata columns to same name??

## Save objects
saveRDS(seuratObj_Human_subset_converted, file = paste0(output.dir,"Robjects/seuratObj_Human_subset_converted.rds"))
saveRDS(Original_human_gene_symbols, file=paste0(output.dir,"Robjects/Original_human_gene_symbols.rds"))
saveRDS(M.genes, file=paste0(output.dir,"Robjects/Converted_human_to_mouse_gene_symbols.rds"))
saveRDS(M.genes_df, file=paste0(output.dir,"Robjects/Converted_human_to_mouse_gene_symbols_df.rds"))

#############################################3
#############################################

## Merge seuratobjects
seuratObj <- merge(seuratObj_Mouse, y = seuratObj_Human_subset_converted, 
                   project = "FB_species_merge", merge.data = T)

dim(seuratObj)
# [1] 23557 12616

diagnostics<-list()

unique(sapply(X = strsplit(colnames(seuratObj), split = "_"), FUN = "[", 1))

## Set directories and names again
sampleName <- "Integrated_Human_and_Mouse_CP"
output.dir <- sampleFolder

dir.create(paste0(sampleFolder,"Plots/"))
dir.create(paste0(sampleFolder,"Plots/RNA/"))

########## Get HVG ##########
seuratObj <- FindVariableFeatures(object = seuratObj, assay = "RNA", selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(seuratObj, assay = "RNA"))
# 2000

### Add to diagnostics
diagnostics[['varGenes']]<-length(VariableFeatures(seuratObj, assay = "RNA"))

########## Scale ##########
seuratObj <- ScaleData(seuratObj, assay = "RNA")

# Run PCA on rna normalized through scran/scater
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj, assay = "RNA"), 
                    npcs = 150, ndims.print = 1:5, nfeatures.print = 10, assay = "RNA")

## Use RNA going forward (to avoid mistakes!!!!!!!)
DefaultAssay(object = seuratObj)<-"RNA"

## Add metadata column (Species)
seuratObj@meta.data$Species<-"Human"
seuratObj@meta.data[which(seuratObj@meta.data$Lab == "Lehtinen"),"Species"] <- "Mouse"

################################################################################
########## RUN HARMONY
################################################################################
library('cowplot')
library("harmony")

########## Create vlnPlot before running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "pca", pt.size = 0.2, group.by = "Species") #Adapted
p2 <- VlnPlot(object = seuratObj, features = "PC_1", pt.size = 0.2, group.by = "Species") #Adapted
plot_grid(p1,p2)
# ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/RNA/1a_vlnPlot_beforeAlignment.png"))

########## Run Harmony ##########
### Increase theta parameter in case of bad overlap!
options(repr.plot.height = 3, repr.plot.width = 6)
seuratObj<-RunHarmony(seuratObj, group.by.vars = "Species", theta = 4, plot_convergence = TRUE, nclust = 50, #Adapted
                      max.iter.cluster = 100, max.iter.harmony = 20, dims.use=1:40)


### Get embeddings
harmony_embeddings <- Embeddings(seuratObj, 'harmony')
harmony_embeddings[1:5, 1:5]


########## Create vlnPlot after running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "harmony", pt.size = 0.2, group.by = "Species") #Adapted
p2 <- VlnPlot(object = seuratObj, features = "harmony_1", pt.size = 0.2, group.by = "Species") #Adapted
plot_grid(p1,p2)
ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/RNA/1b_vlnPlot_afterAlignment.png"))


########################################
########## Choose dims
########################################

########## Via PCelbowplot ##########
ElbowPlot(object = seuratObj, ndims = 40)


dimsToTry<-c(seq(15,30,by=5))

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "harmony", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  # ##### Create tSNE plot
  # seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "RNA", reduction = "harmony")
  # tsnePlot<-DimPlot(seuratObj, reduction = "tsne", label=T, label.size = 8)
  # tsnePlotSplit<-DimPlot(seuratObj, reduction = "tsne", label=F, group.by="orig.ident", pt.size = 2)
  # 
  # ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
  #        file=paste0(output.dir,"Plots/RNA/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="harmony")
  umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

### Final
dimsToTry<-c(30)
resToUse<-0.8
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "harmony", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "RNA", reduction = "harmony")
  tsnePlot<-DimPlot(seuratObj, reduction = "tsne", label=T, label.size = 8)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "tsne", label=F, group.by="orig.ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="harmony")
  umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

# names(seuratObj)
# 
# ### Clustering: trying out clusTree
# 
# Perplexity<-25
# Resolution<-0.8
# Perplexity_UMAP<-25
# 
# seuratObj <- FindNeighbors(object = seuratObj, reduction = "harmony", dims = 1:Perplexity)
# resolutions <- seq(0,1,by=0.1)
# 
# for(res in resolutions){
#   seuratObj <- FindClusters(object = seuratObj,  resolution = res)
# }
# 
# pdf(file=paste0(output.dir,"Plots/RNA/10c_Clustree.pdf"))
# clustree(seuratObj, prefix = "RNA_snn_res.")
# dev.off()
# 
# # 0.8 seems reasonable
# # Final Resolution and final clusters
# res <- 0.8
# diagnostics[['res']]<-res
# seuratObj$harmony_clusters <- seuratObj$RNA_snn_res.0.8
# 
# 
# ################################################################################
# ################################################################################
# ### AUTOMATIC PART
# ################################################################################
# ################################################################################
# 
# # seuratObj <- RunTSNE(seuratObj, reduction = "SCT_pca", dims = 1:Perplexity, assay = "SCT")
# TSNEPlot(seuratObj)
# 
# # seuratObj <- RunUMAP(seuratObj, dims = 1:Perplexity_UMAP, reduction = "SCT_pca", assay = "SCT")
# umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, group.by= "harmony_clusters", label.size = 6)
# tsnePlot<-TSNEPlot(seuratObj, reduction = "tsne", label = T, group.by= "harmony_clusters", label.size = 6)
# seuratObj@active.assay
# 
# pdf(file=paste0(output.dir,"Plots/RNA/11_tSNE_UMAP.pdf"), width = 17*0.45, height = 12.4*0.45)
# umapPlot
# tsnePlot
# dev.off()

experiment<-sampleName

# Save objects
saveRDS(seuratObj, file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_harmony_RNA.rds"))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment, "_harmony_RNA.rds"))

# Read objects
seuratObj<-readRDS(file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_harmony_RNA.rds"))
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment, "_harmony_RNA.rds"))

