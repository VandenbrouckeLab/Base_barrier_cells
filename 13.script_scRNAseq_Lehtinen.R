## Lehtinen paper prep: scRNAseq embryonal data

# options(repos = c(getOption("repos"), BiocManager::repositories()))
# getOption("repos")
# packrat::get_opts()
# packrat::status()
# packrat::snapshot()

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

source('/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/RAW_DATA/script_functions_COVID.R') #KEVIN

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/FB_datasets/")

sampleName<-"Lehtinen_CP"
sampleFolder<-"Lehtinen_CP_scRNAseq/"

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
rawData_1 <- as.matrix(Read10X("~/VIB/DATA/Roos/Daan 1/FB_datasets/Lehtinen_CP_scRNAseq/Cell_atlas/SCP1365/expression/Counts/", gene.column = 1)) # Only 1 column for features (not first ensemble names!!)
rawData_2 <- as.matrix(Read10X("~/VIB/DATA/Roos/Daan 1/FB_datasets/Lehtinen_CP_scRNAseq/Cell_atlas/SCP1365/expression/Data/", gene.column = 1)) # Only 1 column for features (not first ensemble names!!)

## Create seuratobject from count data without removing any cells
seuratObj <- CreateSeuratObject(counts = rawData_1, project = "Lehtinen_CP", min.cells = 0, min.features = 0) # No filter
seuratObj

## Read in metadata
metaData<-read.table("~/VIB/DATA/Roos/Daan 1/FB_datasets/Lehtinen_CP_scRNAseq/Cell_atlas/SCP1365/metadata/scEmbryo_cell_metadata.txt", header = T, sep = "\t")
metaData<-metaData[-1,]
colnames(rawData_1)
metaData$NAME
intersect(colnames(rawData_1),metaData$NAME)
nrow(metaData[which(metaData$cell_type__ontology_label == "mesenchymal cell"),])

## Read in tSNE coordinates
UMAP_mesenchymal<-read.table("~/VIB/DATA/Roos/Daan 1/FB_datasets/Lehtinen_CP_scRNAseq/Cell_atlas/SCP1365/cluster/mesenchymal_tSNE_coordinates.txt", header = T, sep = "\t")
UMAP_all<-read.table("~/VIB/DATA/Roos/Daan 1/FB_datasets/Lehtinen_CP_scRNAseq/Cell_atlas/SCP1365/cluster/all_cells_tSNE_coordinates.txt", header = T, sep = "\t")

# Remove type row
UMAP_mesenchymal<-UMAP_mesenchymal[-1,]
UMAP_all<-UMAP_all[-1,]

## Check overlap
intersect(UMAP_mesenchymal$NAME,colnames(seuratObj))
intersect(UMAP_all$NAME,colnames(seuratObj))

## Info paper
# 15,620 single cells (scRNA-seq)
# 3894 Mesenchyml cells in tSNE coordinates
# 7810 single cells with metadata
# 1949 Mesenchymal cells annotated in metadata

## Issue that metadata doesn't contain all the cells!!!
## Subset based on tSNE coordinates of just mesenchymal cells!!
seuratObj_subset<-subset(seuratObj, cells = UMAP_mesenchymal$NAME)

## Prepro and add coordinates to mesenchymal object
rownames(UMAP_mesenchymal)<-UMAP_mesenchymal$NAME
UMAP_mesenchymal<-UMAP_mesenchymal[-1]
colnames(UMAP_mesenchymal)<-c("tSNE_1","tSNE_2") #Add key to column names!
UMAP_mesenchymal<-as.matrix(UMAP_mesenchymal) #Needs to be a matrix
class(UMAP_mesenchymal)<-"numeric" #Works!!!! Otherwise as character in seuratobject!!
seuratObj_subset[["tSNE"]]<-CreateDimReducObject(embeddings = UMAP_mesenchymal, key = "tSNE", assay="RNA")
seuratObj_subset@meta.data$Annotation<-"Mesenchymal cells" # Add annotation

## Prepro and add coordinates to full object
rownames(UMAP_all)<-UMAP_all$NAME
UMAP_all<-UMAP_all[-1]
colnames(UMAP_all)<-c("tSNE_1","tSNE_2")
UMAP_all<-as.matrix(UMAP_all)
class(UMAP_all)<-"numeric" #Works!!!!
seuratObj[["tSNE"]]<-CreateDimReducObject(embeddings = UMAP_all, key = "tSNE", assay="RNA")
# seuratObj@meta.data$Annotation<-"Mesenchymal cells"

## Create and save tSNEs as check!!
D1<- DimPlot(seuratObj_subset, reduction = "tSNE", label = T, group.by = "Annotation", label.size = 8)
D1.5<- DimPlot(seuratObj_subset, reduction = "tSNE", label = T, label.size = 8)
D2<- DimPlot(seuratObj, reduction = "tSNE", label = T,  label.size = 3)

ggsave(grid.arrange(D1, D1.5, ncol=2), file=paste0(sampleFolder,"results/QC/1_Mesenchymal_tSNE_",sampleName,".png"), width = 20, height=8)
ggsave(D2, file=paste0(sampleFolder,"results/QC/1_All_cells_tSNE_",sampleName,".png"), width = 12)


## Lognormalize data separately!! Don't use their data! Same workflow as our data!!
# Test<-CreateAssayObject(data =  rawData_2, min.cells = 0, min.features = 0)
seuratObj_subset <- NormalizeData(object = seuratObj_subset, normalization.method = "LogNormalize", scale.factor = 10000)

###############################################################################################

## Get extra metadata !!
## Convert to sce
sce<-SingleCellExperiment(list(counts=seuratObj_subset@assays$RNA@counts))

###############################################################################################
## Analyze sce object

##### Get mitochondrial genes #####
is.mito <- grepl("^MT-", rownames(sce), ignore.case = TRUE)
sum(is.mito)
##0
rownames(sce)[is.mito]

##### Calculate QC metrics #####
sce<- addPerCellQC(sce)
dim(colData(sce))

### List samples
listLabels<-"Lehtinen_CP"


##### Create metaData_extra matrix (used for downstream analysis) #####
metaData_extra<-data.frame("staticNr"=colnames(sce),"orig.ident"=listLabels[[1]], "nGene"=sce$detected,"nUMI"=sce$sum,
                     stringsAsFactors = F)
rownames(metaData_extra)<-metaData_extra$staticNr
metaData_extra$staticNr<-1

table(metaData_extra$orig.ident)

## Add to metadata3
seuratObj_subset@meta.data<-cbind(seuratObj_subset@meta.data,metaData_extra)

################################################################################
########## Save object
################################################################################

##### Save object
saveRDS(seuratObj_subset, file=paste0(sampleFolder,"Robjects/seuratObj_mesenchymal_",sampleName,".rds"))
