## Script for processing 7-, 22-wo and 82-wo ChP 4V and LV samples from our lab
## Subset to only Fibroblasts based on metadata processed objects
## Run until log normalization
## Save seuratobject

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')

source('/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/RAW_DATA/script_functions_COVID.R') #KEVIN

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")

sampleName<-"Urvb_6datasets"
sampleFolder<-paste0(sampleName,"/")

experiment <- sampleName
output.dir <- sampleFolder

##add some subfolders
dir.create(paste0(sampleFolder,"results"))
dir.create(paste0(sampleFolder,"results/QC"))
dir.create(paste0(sampleFolder,"Robjects"))


########################################
##### Some variables
########################################

### General variables
diagnostics<-list()

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')


################################################################################
########## LOAD DATA
################################################################################
rawData_1 <- as.matrix(Read10X("~/VIB/DATA/Roos/Daan 1/LpsNegFourVentr/filtered_gene_bc_matrices/mm10/"))
rawData_2 <- as.matrix(Read10X("~/VIB/DATA/Roos/Daan 1/LpsNegLatVentr/filtered_gene_bc_matrices/mm10/"))
rawData_3 <- as.matrix(Read10X("~/VIB/DATA/Roos/Nina/RVD5_Y4V/filtered_gene_bc_matrices/mm10/"))
rawData_4 <- as.matrix(Read10X("~/VIB/DATA/Roos/Nina/RVD6_YLV/filtered_gene_bc_matrices/mm10/"))
rawData_5 <- as.matrix(Read10X("~/VIB/DATA/Roos/Nina/RVD7_O4V/filtered_gene_bc_matrices/mm10/"))
rawData_6 <- as.matrix(Read10X("~/VIB/DATA/Roos/Nina/RVD8_OLV/filtered_gene_bc_matrices/mm10/"))

seuratObj_1 <- readRDS(file = "~/VIB/DATA/Roos/Daan 1/LpsNegFourVentr/Robjects/seuratObj_final_LpsNegFourVentr.rds")
seuratObj_2 <- readRDS(file = "~/VIB/DATA/Roos/Daan 1/LpsNegLatVentr/Robjects/seuratObj_sliced3_LpsNegLatVentr.rds")
seuratObj_3 <- readRDS(file = "~/VIB/DATA/Roos/Nina/RVD5_Y4V/Robjects/seuratObj_RVD5_Y4V.rds")
seuratObj_4 <- readRDS(file = "~/VIB/DATA/Roos/Nina/RVD6_YLV/Robjects/seuratObj_RVD6_YLV.rds")
seuratObj_5 <- readRDS(file = "~/VIB/DATA/Roos/Nina/RVD7_O4V/Robjects/seuratObj_slicedRVD7_O4V.rds")
seuratObj_6 <- readRDS(file = "~/VIB/DATA/Roos/Nina/RVD8_OLV/Robjects/seuratObj_sliced_RVD8_OLV.rds")

##Annotation
## LpsNeg4V
## Update 04/03/21
seuratObj_1@meta.data$annotated_clusters <- seuratObj_1@active.ident
seuratObj_1@meta.data$annotated_clusters<- factor(seuratObj_1@meta.data$annotated_clusters,levels(seuratObj_1@meta.data$annotated_clusters)[c(2:17,1)]) #reorder levels
levels(seuratObj_1@meta.data$annotated_clusters) <- c(rep("Epithelial cells",2), "Macrophages", "Endothelial cells", rep("Epithelial cells",2),"Fibroblasts Type 1", 
                                                    "Macrophages", "Epithelial cells", "Vascular associated cells", "NK cells",  "Fibroblasts Type 2",
                                                    "Doublets", rep("Other immune cells",3), "Doublets 2")  # based on markers (change 7.11.19 cl12 to doublets)

## LpsNegLV
### Create annotated UMAP ###
seuratObj_2@meta.data$annotated_clusters <- seuratObj_2@active.ident
seuratObj_2@meta.data$annotated_clusters <- factor(seuratObj_2@meta.data$annotated_clusters,levels(seuratObj_2@meta.data$annotated_clusters)[c(5:18,4,3,2,1)]) #reorder levels
levels(seuratObj_2@meta.data$annotated_clusters) <- c(rep("Epithelial cells",4), "Macrophages",rep("Epithelial cells",3), "Endothelial cells", "Fibroblasts", 
                                                    "Microglia-like Macrophages", "Other immune cells", "Vascular associated cells", "Other immune cells","Neuronal cells", 
                                                    "Doublets", "Xist+ Epithelial cells","NK cells")  # based on markers

## Y4V
### Create annotated UMAP
seuratObj_3@meta.data$annotated_clusters <- seuratObj_3@active.ident
# seuratObj_3@meta.data$annotated_clusters <- factor(seuratObj_3@meta.data$annotated_clusters,levels(seuratObj_3@meta.data$annotated_clusters)[c(3:15,2,1)]) #reorder levels
levels(seuratObj_3@meta.data$annotated_clusters) <- c("Epithelial cells 1","Epithelial cells 2","Epithelial cells 1","Endothelial cells","Macrophages",'Fibroblasts Type 1',
                                                    'Epithelial cells 2','Epithelial cells 1',"Vascular associated cells",'Epithelial cells 1',"Doublets_1",
                                                    "Doublets_2",'NK cells','Other Macrophages','Other Immune Cells','Doublets_3',"Fibroblasts Type 2","Doublets_4") #split

## YLV
### Create annotated UMAP
seuratObj_4@meta.data$annotated_clusters <- seuratObj_4@active.ident
# seuratObj_4@meta.data$annotated_clusters <- factor(seuratObj_4@meta.data$annotated_clusters,levels(seuratObj_4@meta.data$annotated_clusters)[c(3:15,2,1)]) #reorder levels
levels(seuratObj_4@meta.data$annotated_clusters) <- c("Epithelial cells","Epithelial cells","Epithelial cells","Epithelial cells","Epithelial cells","Epithelial cells",
                                                    "Epithelial cells","Epithelial cells","Epithelial cells","Endothelial cells","Macrophages","Vascular associated cells",
                                                    'Fibroblasts',"Doublets_1","Doublets_2",'NK cells','Other Immune Cells','Microglia-like Macrophages')

## O4V
### Create annotated UMAP
seuratObj_5@meta.data$annotated_clusters <- seuratObj_5@active.ident
seuratObj_5@meta.data$annotated_clusters <- factor(seuratObj_5@meta.data$annotated_clusters,levels(seuratObj_5@meta.data$annotated_clusters)[c(3:20,2,1)]) #reorder levels
levels(seuratObj_5@meta.data$annotated_clusters) <- c("Epithelial cells","Epithelial cells","Epithelial cells","Fibroblasts Type 1","Macrophages","Epithelial cells",
                                                    "Endothelial cells","Vascular associated cells","Epithelial cells","Doublets_1","Doublets_2","Doublets_3",
                                                    'NK cells','Other Immune Cells','Fibroblasts Type 2',"Doublets_4","Doublets_5",'Other Vascular associated cells',
                                                    "Doublets_6",'T-cells')

## OLV
### Create annotated UMAP
seuratObj_6@meta.data$annotated_clusters <- seuratObj_6@active.ident
seuratObj_6@meta.data$annotated_clusters <- factor(seuratObj_6@meta.data$annotated_clusters,levels(seuratObj_6@meta.data$annotated_clusters)[c(3:18,2,1)]) #reorder levels
levels(seuratObj_6@meta.data$annotated_clusters) <- c(rep("Epithelial cells",7),"Endothelial cells",'Fibroblasts',"Macrophages",
                                                    "Doublets_1",rep("Epithelial cells",2),"Doublets_2",'Other Immune Cells',
                                                    'NK cells',"Doublets_3","Vascular associated cells")

## Update idents
Idents(seuratObj_1)<-seuratObj_1@meta.data$annotated_clusters
Idents(seuratObj_2)<-seuratObj_2@meta.data$annotated_clusters
Idents(seuratObj_3)<-seuratObj_3@meta.data$annotated_clusters
Idents(seuratObj_4)<-seuratObj_4@meta.data$annotated_clusters
Idents(seuratObj_5)<-seuratObj_5@meta.data$annotated_clusters
Idents(seuratObj_6)<-seuratObj_6@meta.data$annotated_clusters

## Check annotation
DimPlot(seuratObj_1, reduction = "umap", label = T, repel = T, label.size = 4)
DimPlot(seuratObj_2, reduction = "umap", label = T, repel = T, label.size = 4)
DimPlot(seuratObj_3, reduction = "umap", label = T, repel = T, label.size = 4)
DimPlot(seuratObj_4, reduction = "umap", label = T, repel = T, label.size = 4)
DimPlot(seuratObj_5, reduction = "umap", label = T, repel = T, label.size = 4)
DimPlot(seuratObj_6, reduction = "umap", label = T, repel = T, label.size = 4)

## Choose FBs
Cells1<-WhichCells(seuratObj_1, idents = c("Fibroblasts Type 1", "Fibroblasts Type 2"))
Cells2<-WhichCells(seuratObj_2, idents = c("Fibroblasts"))
Cells3<-WhichCells(seuratObj_3, idents = c("Fibroblasts Type 1", "Fibroblasts Type 2"))
Cells4<-WhichCells(seuratObj_4, idents = c("Fibroblasts"))
Cells5<-WhichCells(seuratObj_5, idents = c("Fibroblasts Type 1", "Fibroblasts Type 2"))
Cells6<-WhichCells(seuratObj_6, idents = c("Fibroblasts"))

## Subset rawData
rawData_1_subset<-rawData_1[,Cells1]
rawData_2_subset<-rawData_2[,Cells2]
rawData_3_subset<-rawData_3[,Cells3]
rawData_4_subset<-rawData_4[,Cells4]
rawData_5_subset<-rawData_5[,Cells5]
rawData_6_subset<-rawData_6[,Cells6]

### Change colnames
colnames(rawData_1_subset)<-paste0(colnames(rawData_1_subset),"-1")
colnames(rawData_2_subset)<-paste0(colnames(rawData_2_subset),"-2")
colnames(rawData_3_subset)<-paste0(colnames(rawData_3_subset),"-3")
colnames(rawData_4_subset)<-paste0(colnames(rawData_4_subset),"-4")
colnames(rawData_5_subset)<-paste0(colnames(rawData_5_subset),"-5")
colnames(rawData_6_subset)<-paste0(colnames(rawData_6_subset),"-6")

### Check before merge
cbind(head(rownames(rawData_1_subset)),head(rownames(rawData_2_subset)),head(rownames(rawData_3_subset)),
      head(rownames(rawData_4_subset)),head(rownames(rawData_5_subset)),head(rownames(rawData_6_subset)))
cbind(tail(rownames(rawData_1_subset)),tail(rownames(rawData_2_subset)),tail(rownames(rawData_3_subset)),
      tail(rownames(rawData_4_subset)),tail(rownames(rawData_5_subset)),tail(rownames(rawData_6_subset)))

dim(rawData_1_subset)
dim(rawData_2_subset)
dim(rawData_3_subset)
dim(rawData_4_subset)
dim(rawData_5_subset)
dim(rawData_6_subset)

sum(ncol(rawData_1_subset),ncol(rawData_2_subset),ncol(rawData_3_subset),ncol(rawData_4_subset),ncol(rawData_5_subset),ncol(rawData_6_subset))
#2580 cells

### Do merge
rawDataRNA<-cbind(rawData_1_subset, rawData_2_subset,rawData_3_subset,rawData_4_subset,rawData_5_subset,rawData_6_subset)
dim(rawDataRNA)
diagnostics[['dimBeforeSeuratObj']]<-paste0(nrow(rawDataRNA)," genes - ",ncol(rawDataRNA)," cells")

### Remove some variables
rm(rawData_1)
rm(rawData_2)
rm(rawData_3)
rm(rawData_4)
rm(rawData_5)
rm(rawData_6)
gc()

#############################
########## PREP DATA
#############################
message("########################## Preparing Data ##########################")

rownames(rawDataRNA) <- stringr::str_remove(rownames(rawDataRNA),"GRCh38.99_____________")
rownames(rawDataRNA) <- stringr::str_remove(rownames(rawDataRNA),"SARS-CoV-2_cellranger_")

#diagnostics
diagnostics[['dimRawDataRNA']]<-paste0(nrow(rawDataRNA)," genes - ",ncol(rawDataRNA)," cells")
diagnostics[['nrGenes']]<-nrow(rawDataRNA)
diagnostics[['nrCells']]<-ncol(rawDataRNA)

#### adding the batch variable from the hashing
batch <- rep(experiment, colnames(rawDataRNA) %>% length())
cells.use  <- colnames(rawDataRNA)

nrZeros<-sum(rawDataRNA==0)/(nrow(rawDataRNA)*ncol(rawDataRNA))*100

##### In each cell: how many genes are expressed (count > 0) #####
cellCounts<-apply(rawDataRNA,2,function (x){sum(x>0)})
##### For each gene: in how many cells is it expressed? (count > 0) #####
geneCounts<-apply(rawDataRNA,1,function (x){sum(x>0)})

##### Add to diagnostics #####
diagnostics[['dimRawData']]<-paste0(nrow(rawDataRNA)," genes - ",ncol(rawDataRNA)," cells")
diagnostics[['nrGenes']]<-nrow(rawDataRNA)
diagnostics[['nrCells']]<-ncol(rawDataRNA)
diagnostics[['zeroInflation']]<-nrZeros
diagnostics[['minGenesPerCell']]<-min(cellCounts)
diagnostics[['maxGenesPerCell']]<-max(cellCounts)
diagnostics[['meanGenesPerCell']]<-mean(cellCounts)
diagnostics[['medianGenesPerCell']]<-median(cellCounts)
diagnostics[['cellsLess200genes']]<-length(cellCounts[cellCounts<200])
diagnostics[['genesNotExpressed']]<-length(geneCounts[geneCounts==0])
diagnostics[['genesLess3cells']]<-length(geneCounts[geneCounts<3])

### Remove some variables
rm(geneCounts)
rm(cellCounts)

###########################
########## QC: CELLS
###########################

##### Create object #####
sce<-SingleCellExperiment(list(counts=rawDataRNA))
diagnostics[['dimSce']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")

##### Get mitochondrial genes #####
is.mito <- grepl("^MT-", rownames(sce), ignore.case = TRUE)
sum(is.mito)
##13
rownames(sce)[is.mito]


##### Calculate QC metrics #####
### => pData(sce) is created
sce<- addPerCellQC(sce, subsets=list(Mito=is.mito))
dim(colData(sce))
# colnames(colData(sce))

### List samples
listLabels<-c("RVD1_LpsNegFour","RVD2_LpsNegLat","RVD5_Y4V","RVD6_YLV","RVD7_O4V","RVD8_OLV")

##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("staticNr"=colnames(rawDataRNA),"orig.ident"=listLabels[[1]], "nGene"=sce$detected,"nUMI"=sce$sum,
                     "percent.mito"=sce$subsets_Mito_percent,
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1


for(i in 2:length(listLabels)){
  toSearch<-paste0("-",i)
  metaData[grep(toSearch,rownames(metaData)), which(colnames(metaData)=="orig.ident")]<-listLabels[[i]]
}

table(metaData$orig.ident)

##### Add to diagnostics #####
diagnostics[['splitSamples']]<-paste0(table(metaData$orig.ident)," cells of sample ",rownames(table(metaData$orig.ident)))

################################################################################
########## FINALIZE QC
################################################################################

dim(sce)
# saveRDS(sce, file=paste0(sampleFolder,"Robjects/sce.rds"))
# sce <- readRDS(file=paste0(sampleFolder,"Robjects/sce.rds"))

rawDataFiltered<-rawDataRNA[rownames(sce),colnames(sce)]
dim(rawDataFiltered)
# 27998  2580
diagnostics[['dimBeforeSeuratObj']]<-paste0(nrow(rawDataFiltered)," genes - ",ncol(rawDataFiltered)," cells")

### Remove some variables
rm(sce)
rm(rawData)


################################################################################
########## CREATE SEURAT OBJECT
################################################################################

##### Create object #####
seuratObj <- CreateSeuratObject(counts = rawDataFiltered, project = "seuratObj", min.cells = 3, min.features = 200)

dim(seuratObj)
#14582  2580

GetAssayData(seuratObj, assay = "RNA", slot="counts")[1:5,1:5]
seuratObj[['RNA']]@counts[1:5,1:5]

##### Add to diagnostics #####
diagnostics[['dimAfterSeuratObj']]<-paste0(nrow(seuratObj)," genes - ",ncol(seuratObj)," cells")


################################################################################
########## FILTER DATA
################################################################################

seuratObj[["percent.mito"]] <- PercentageFeatureSet(object = seuratObj, pattern = "^mt-")
head(seuratObj@meta.data)

png(file=paste0(sampleFolder,"results/QC/4_vlnPlotSeurat.png"), width = 850, height = 642)
VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"))
dev.off()

##### Add orig ident
metaDataTable<-seuratObj@meta.data

metaDataTable$orig.ident<-as.character(metaDataTable$orig.ident)
for(i in 1:length(listLabels)){
  toSearch<-paste0('-',i)
  metaDataTable[grep(toSearch,rownames(metaDataTable)), which(colnames(metaDataTable)=="orig.ident")]<-listLabels[[i]]
}
seuratObj@meta.data<-metaDataTable

head(metaDataTable)
table(metaDataTable$orig.ident)

##### Add to diagnostics #####
diagnostics[['splitSamplesAfterFiltering']]<-paste0(table(metaDataTable$orig.ident)," cells of sample ",rownames(table(metaDataTable$orig.ident)))


################################################################################
########## NORMALIZE
################################################################################
seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)

##### Check per group #####
metaDataTable<-seuratObj@meta.data
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@data),2,function(x){sum(x>0)})

drawVlnPlotSeurat_split(metaDataTable, paste0(sampleFolder,"results/QC/5_afterNorm_splitted.png"))


################################################################################
########## Save object
################################################################################

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_final_",sampleName,".rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_final_",sampleName,".rds"))


## Save metadata
metaData_FBs1<-seuratObj_1@meta.data[Cells1,c("orig.ident","annotated_clusters")]
rownames(metaData_FBs1)<-paste0(rownames(metaData_FBs1),"_",1)
metaData_FBs2<-seuratObj_2@meta.data[Cells2,c("orig.ident","annotated_clusters")]
rownames(metaData_FBs2)<-paste0(rownames(metaData_FBs2),"_",2)
metaData_FBs3<-seuratObj_3@meta.data[Cells3,c("orig.ident","annotated_clusters")]
rownames(metaData_FBs3)<-paste0(rownames(metaData_FBs3),"_",3)
metaData_FBs4<-seuratObj_4@meta.data[Cells4,c("orig.ident","annotated_clusters")]
rownames(metaData_FBs4)<-paste0(rownames(metaData_FBs4),"_",4)
metaData_FBs5<-seuratObj_5@meta.data[Cells5,c("orig.ident","annotated_clusters")]
rownames(metaData_FBs5)<-paste0(rownames(metaData_FBs5),"_",5)
metaData_FBs6<-seuratObj_6@meta.data[Cells6,c("orig.ident","annotated_clusters")]
rownames(metaData_FBs6)<-paste0(rownames(metaData_FBs6),"_",6)
metaData_FBs<-rbind(metaData_FBs1,metaData_FBs2,metaData_FBs3,metaData_FBs4,metaData_FBs5,metaData_FBs6)
saveRDS(metaData_FBs, file = "Urvb_6datasets/Robjects/MetaData_FBs.rds")

