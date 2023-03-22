## Creating Fibroblast species object (Fig7) with our Fibroblast scRNA-Seq data (7/22/82 wo ChP 4V&LV) and human snRNA-seq data (Yang et al.)
## Script performs conversion Human gene symbols to mouse gene symbols
## initially part of script which attempted harmony integration

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