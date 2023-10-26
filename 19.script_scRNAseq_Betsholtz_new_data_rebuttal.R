## Rebuttal: Processing new Betsholtz lab data (10/2023)
## https://betsholtzlab.org/Publications/BrainFB/Data/BFBdata.html
## Three datasets to read in, normalize and annotate according to metadata Betsholtz lab
## No extra filtering! Just like their count table!

## Raw read counts!!
## Run until log normalization
## Save seuratobject

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('openxlsx')

###################################################################################################################
###################################################################################################################

######################################################################################
#################### Dataset 1: GSE227713_Betsholtz_leptomeninges #################### 
######################################################################################

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")

sampleName<-"GSE227713_Betsholtz_leptomeninges"
sampleFolder<-"Betsholtz_new/Dissected_leptomeninges/"

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

rawDataSparse <- read.table("Betsholtz_new/Dissected_leptomeninges/Raw_data/dataset1_dissected_leptomeninges_raw_counts.txt.gz", header = T)

rawDataRNA<-as.matrix(rawDataSparse)
dim(rawDataRNA)  #18107  1341

### Remove some variables
rm(rawDataSparse)

#############################
########## PREP DATA
#############################
message("########################## Preparing Data ##########################")

rownames(rawDataRNA) <- stringr::str_remove(rownames(rawDataRNA),"GRCh38.99_____________")

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

################################################################################
########## CREATE SEURAT OBJECT
################################################################################

##### Create object #####
seuratObj <- CreateSeuratObject(counts = rawDataRNA, project = "seuratObj") 

dim(seuratObj)

GetAssayData(seuratObj, assay = "RNA", slot="counts")[1:5,1:5]
seuratObj[['RNA']]@counts[1:5,1:5]

##### Add to diagnostics #####
diagnostics[['dimAfterSeuratObj']]<-paste0(nrow(seuratObj)," genes - ",ncol(seuratObj)," cells")

seuratObj[["percent.mito"]] <- PercentageFeatureSet(object = seuratObj, pattern = "^mt-")
head(seuratObj@meta.data)

png(file=paste0(sampleFolder,"results/QC/4_vlnPlotSeurat.png"), width = 850, height = 642)
VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"))
dev.off()

##### Add orig ident
metaDataTable<-seuratObj@meta.data
metaDataTable$orig.ident<-as.character(metaDataTable$orig.ident)
metaDataTable[, which(colnames(metaDataTable)=="orig.ident")]<-sampleName

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
########## ANNOTATE
################################################################################
## Read in annotation Betsholtz paper
Dataset1_annotation<-read.xlsx("Betsholtz_new/Dissected_leptomeninges/Raw_data/dataset1_dissected_leptomeninges_annotation.xlsx")

rownames(Dataset1_annotation)<-Dataset1_annotation$Cells
all(rownames(Dataset1_annotation)==rownames(seuratObj@meta.data))

seuratObj@meta.data$Annotation_Betsholtz_paper<-Dataset1_annotation$annot

################################################################################
########## Save object
################################################################################

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_final_",sampleName,".rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_final_",sampleName,".rds"))

###################################################################################################################
###################################################################################################################

######################################################################################
#################### Dataset 2: GSE233270_Betsholtz_dural #################### 
######################################################################################

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")

sampleName<-"GSE233270_Betsholtz_dural"
sampleFolder<-"Betsholtz_new/Dissected_dural/"

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

rawDataSparse <- read.table("Betsholtz_new/Dissected_dural/Raw_data/dataset2_dissected_dural_raw_counts.txt.gz", header = T)

rawDataRNA<-as.matrix(rawDataSparse)
dim(rawDataRNA) #20685   345

### Remove some variables
rm(rawDataSparse) 

#############################
########## PREP DATA
#############################
message("########################## Preparing Data ##########################")

rownames(rawDataRNA) <- stringr::str_remove(rownames(rawDataRNA),"GRCh38.99_____________")

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

################################################################################
########## CREATE SEURAT OBJECT
################################################################################

##### Create object #####
seuratObj <- CreateSeuratObject(counts = rawDataRNA, project = "seuratObj") 

dim(seuratObj)

GetAssayData(seuratObj, assay = "RNA", slot="counts")[1:5,1:5]
seuratObj[['RNA']]@counts[1:5,1:5]

##### Add to diagnostics #####
diagnostics[['dimAfterSeuratObj']]<-paste0(nrow(seuratObj)," genes - ",ncol(seuratObj)," cells")


seuratObj[["percent.mito"]] <- PercentageFeatureSet(object = seuratObj, pattern = "^mt-")
head(seuratObj@meta.data)

png(file=paste0(sampleFolder,"results/QC/4_vlnPlotSeurat.png"), width = 850, height = 642)
VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"))
dev.off()

##### Add orig ident
metaDataTable<-seuratObj@meta.data
metaDataTable$orig.ident<-as.character(metaDataTable$orig.ident)
metaDataTable[, which(colnames(metaDataTable)=="orig.ident")]<-sampleName

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
########## ANNOTATE
################################################################################
## Read in annotation Betsholtz paper
Dataset2_annotation<-read.xlsx("Betsholtz_new/Dissected_dural/Raw_data/dataset2_dissected_dural_annotation.xlsx")

rownames(Dataset2_annotation)<-Dataset2_annotation$Cells
all(rownames(Dataset2_annotation)==rownames(seuratObj@meta.data))

seuratObj@meta.data$Annotation_Betsholtz_paper<-Dataset2_annotation$annot

################################################################################
########## Save object
################################################################################

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_final_",sampleName,".rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_final_",sampleName,".rds"))


###################################################################################################################
###################################################################################################################

######################################################################################
#################### Dataset 4: GSE228882_Betsholtz_FACS #################### 
######################################################################################

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")

sampleName<-"GSE228882_Betsholtz_FACS"
sampleFolder<-"Betsholtz_new/FACS_sorted_cells/"

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

rawDataSparse <- read.table("Betsholtz_new/FACS_sorted_cells/Raw_data/dataset4_FACS_SS2_raw_counts.txt.gz", header = T)

rawDataRNA<-as.matrix(rawDataSparse)
dim(rawDataRNA) #30888   234

### Remove some variables
rm(rawDataSparse) 

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

################################################################################
########## CREATE SEURAT OBJECT
################################################################################

##### Create object #####
seuratObj <- CreateSeuratObject(counts = rawDataRNA, project = "seuratObj")

dim(seuratObj)

GetAssayData(seuratObj, assay = "RNA", slot="counts")[1:5,1:5]
seuratObj[['RNA']]@counts[1:5,1:5]

##### Add to diagnostics #####
diagnostics[['dimAfterSeuratObj']]<-paste0(nrow(seuratObj)," genes - ",ncol(seuratObj)," cells")


seuratObj[["percent.mito"]] <- PercentageFeatureSet(object = seuratObj, pattern = "^mt-")
head(seuratObj@meta.data)

png(file=paste0(sampleFolder,"results/QC/4_vlnPlotSeurat.png"), width = 850, height = 642)
VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"))
dev.off()

##### Add orig ident
metaDataTable<-seuratObj@meta.data
metaDataTable$orig.ident<-as.character(metaDataTable$orig.ident)
metaDataTable[, which(colnames(metaDataTable)=="orig.ident")]<-sampleName

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
########## ANNOTATE
################################################################################
## Read in annotation Betsholtz paper
Dataset4_annotation<-read.xlsx("Betsholtz_new/FACS_sorted_cells/Raw_data/dataset4_FACS_SS2_annotation.xlsx")

## Add X to start of cellnames starting with a number -> Can't start with number in seuratObj!!
Dataset4_annotation[59:154,"Cells"]<-paste0("X",Dataset4_annotation[59:154,"Cells"])
rownames(Dataset4_annotation)<-Dataset4_annotation$Cells
all(rownames(Dataset4_annotation)==rownames(seuratObj@meta.data))

seuratObj@meta.data$Annotation_Betsholtz_paper<-Dataset4_annotation$annot

################################################################################
########## Save object
################################################################################

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_final_",sampleName,".rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_final_",sampleName,".rds"))
