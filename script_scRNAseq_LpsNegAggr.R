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
# devtools::install_github('dambi/DisneyTool', host="github.ugent.be/api/v3", auth_token = 'e5ca75c8c2f815aa7f1195cb0b6b6a3190064707')
#cutils::download.file("https://github.ugent.be/api/v3/repos/dambi/DisneyTool/tarball/master", destfile = "test.zip", method = "curl")
library('DisneyTools')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/")

sampleName<-"LpsNegAggr"
sampleFolder<-paste0(sampleName,"/")

##add some subfolders
dir.create(paste0(sampleFolder,"results"))
dir.create(paste0(sampleFolder,"results/QC"))
dir.create(paste0(sampleFolder,"Robjects"))


########################################
##### Some variables
########################################

### Read from the file "aggregation.csv"
aggrFile<-read.csv(file=paste0(sampleFolder,"aggregation_csv.csv"), stringsAsFactors = F)
if(ncol(aggrFile)==3){
  listLabels<-as.list(paste0(aggrFile[,1],"-",aggrFile[,3]))
}else{
  listLabels<-as.list(aggrFile[,1])
}
listLabels
# #or manual
# listLabels<-list('CS35','CS36')

### General variables
diagnostics<-list()

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')


################################################################################
########## LOAD DATA
################################################################################
rawDataSparse <- Read10X(paste0(sampleFolder,"filtered_gene_bc_matrices_mex/mm10/"))
dim(rawDataSparse)

rawData<-as.matrix(rawDataSparse)
dim(rawData)

# rawData<-as.matrix(rawDataSparse$`Gene Expression`)
# dim(rawData)
# 
# rawDataADT<-as.matrix(rawDataSparse$`Antibody Capture`)
# dim(rawDataADT)

### Remove some variables
rm(rawDataSparse)

########################################
##### Characteristics of data
########################################
rawData[1:5,1:5]

nrZeros<-sum(rawData==0)/(nrow(rawData)*ncol(rawData))*100
nrZeros
## 94.747

##### In each cell: how many genes are expressed (count > 0) #####
cellCounts<-apply(rawData,2,function (x){sum(x>0)})
length(cellCounts[cellCounts<200])
## 148

##### For each gene: in how many cells is it expressed? (count > 0) #####
geneCounts<-apply(rawData,1,function (x){sum(x>0)})
length(geneCounts[geneCounts<3])
## 11279

##### Add to diagnostics #####
diagnostics[['dimRawData']]<-paste0(nrow(rawData)," genes - ",ncol(rawData)," cells")
diagnostics[['nrGenes']]<-nrow(rawData)
diagnostics[['nrCells']]<-ncol(rawData)
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
########## QC: CELLS
################################################################################

########################################
########## Calculate metrics
########################################

##### Create object #####
library("scater")
sce<-SingleCellExperiment(list(counts=rawData))
dim(sce)
diagnostics[['dimSce']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")

##### Get spike inns #####
is.spike <- grepl("^ercc-", rownames(sce))
sum(is.spike)
##0

##### Get mitochondrial genes #####
is.mito <- grepl("^mt-", rownames(sce))
sum(is.mito)
##13
rownames(sce)[is.mito]

##### Calculate QC metrics #####
### => pData(sce) is created
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))
dim(colData(sce))
# colnames(colData(sce))

##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("staticNr"=colnames(rawData),"orig.ident"=listLabels[[1]],
                     "nGene"=sce$total_features_by_counts,"nUMI"=sce$total_counts,"percent.mito"=sce$pct_counts_Mt, 
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

########################################
########## Get outliers
########################################
nmad_low_feature<-4
nmad_high_feature<-4

nmad_low_UMI<-4
nmad_high_UMI<-4

nmad_high_mito<-5

##### Aim: remove cells with low library sizes, low numbers of expressed features and with high mitochondrial proportions
##same as nGene in Seurat pipeline
feature.drop.low <- isOutlier(sce$total_features_by_counts, nmads=nmad_low_feature, type="lower", log=TRUE)
sum(feature.drop.low) ## 164

feature.drop.high <- isOutlier(sce$total_features_by_counts, nmads=nmad_high_feature, type="higher", log=TRUE)
sum(feature.drop.high) ## 0

feature.drop<-as.logical(feature.drop.low + feature.drop.high)
sum(feature.drop) ## 164

##same as UMI in Seurat pipeline
libsize.drop.low <- isOutlier(sce$total_counts, nmads=nmad_low_UMI, type="lower", log=TRUE)
sum(libsize.drop.low) ## 0

libsize.drop.high <- isOutlier(sce$total_counts, nmads=nmad_high_UMI, type="higher", log=TRUE)
sum(libsize.drop.high) ## 0

libsize.drop<-as.logical(libsize.drop.low+libsize.drop.high)
sum(libsize.drop) ## 0

##% mitochondrial genes
mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=nmad_high_mito, type="higher")
sum(mito.drop) ## 337

##### add to metaData matrix #####
metaData$nGene.drop=feature.drop
metaData$nUMI.drop=libsize.drop
metaData$mito.drop=mito.drop
metaData$final.drop=feature.drop | libsize.drop | mito.drop

##### Add to diagnostics #####
diagnostics[['nmad.low.feature']]<-nmad_low_feature
diagnostics[['nmad.high.feature']]<-nmad_high_feature
diagnostics[['nmad.low.libsize']]<-nmad_low_UMI
diagnostics[['nmad.high.libsize']]<-nmad_high_UMI
diagnostics[['nmad.high.mito']]<-nmad_high_mito

diagnostics[['feature.drop.low']]<-sum(feature.drop.low)
diagnostics[['feature.drop.high']]<-sum(feature.drop.high)
diagnostics[['feature.drop']]<-sum(feature.drop)
diagnostics[['libsize.drop.low']]<-sum(libsize.drop.low)
diagnostics[['libsize.drop.high']]<-sum(libsize.drop.high)
diagnostics[['libsize.drop']]<-sum(libsize.drop)
diagnostics[['mito.drop']]<-sum(mito.drop)

########################################
########## Create histogram + barplot
########################################
palette(c("#00BFC4","#F8766D","#7CAE00","#C77CFF"))
###00BFC4=cyan
###F8766D=red
###7CAE00=green
###C77CFF=purple

toPlot<-metaData
savePlots<-TRUE

##nGene
if(savePlots==TRUE){png(file=paste0(sampleFolder,"results/QC/1a_nGene.png"), width=850)}
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nGene),]
hist(tmp$nGene, breaks=30)
theColors<-as.factor(tmp$nGene.drop)
barplot(tmp$nGene, col=theColors, border=theColors)
if(savePlots==TRUE){dev.off()}

##nUMI
if(savePlots==TRUE){png(file=paste0(sampleFolder,"results/QC/1b_nUMI.png"), width=850)}
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nUMI),]
hist(tmp$nUMI, breaks=30)
theColors<-as.factor(tmp$nUMI.drop)
barplot(tmp$nUMI, col=theColors, border=theColors)
if(savePlots==TRUE){dev.off()}

##percent.mito
if(savePlots==TRUE){png(file=paste0(sampleFolder,"results/QC/1c_percMito.png"), width=850)}
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$percent.mito),]
hist(tmp$percent.mito, breaks=30)
theColors<-as.factor(tmp$mito.drop)
barplot(tmp$percent.mito, col=theColors, border=theColors)
if(savePlots==TRUE){dev.off()}

########################################
########## Create violinPlots
########################################

### Before filtering
toPlot<-metaData
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))
drawVlnPlot_split(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering_splitted.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))


drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering_nGene.png"), 
            colsToColor = c('nGene.drop','nGene.drop','nGene.drop'))
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering_nUMI.png"), 
            colsToColor = c('nUMI.drop','nUMI.drop','nUMI.drop'))
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering_mito.png"), 
            colsToColor = c('mito.drop','mito.drop','mito.drop'))


### After filtering
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2b_afterFiltering.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))
drawVlnPlot_split(toPlot, fileName = paste0(sampleFolder,"results/QC/2b_afterFiltering_splitted.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))


########################################
########## Remove outliers
########################################

sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
dim(sce)

### Number of cells removed
nrow(metaData)-ncol(sce)
## 340 (large overlap between mito and nfeature)

##### Add to diagnostics #####
diagnostics[['firstRemove']]<-nrow(metaData)-ncol(sce)
diagnostics[['dimFirstRemove']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")


########################################
########## Create PCA
########################################
library('mvoutlier')

# ##(check via raw code of the function runPCA)
# varsToUse <- c("pct_counts_in_top_100_features",
#                "total_features_by_counts", "pct_counts_feature_control",
#                "total_features_by_counts_feature_control", "log10_total_counts_endogenous",
#                "log10_total_counts_feature_control")
# setdiff(varsToUse, colnames(colData(sce)))
# exprs_to_plot <- scale(colData(sce)[,varsToUse], scale = T)
# x.mad = apply(exprs_to_plot, 2, mad)
# x.mad[x.mad==0]
# varsToUse<-setdiff(varsToUse, names(x.mad[x.mad==0]))
# sceNew<-runPCA(sce,use_coldata=T, detect_outliers=T, selected_variables = varsToUse)

##### Detect bad cells #####
sceNew<-runPCA(sce,use_coldata=T, detect_outliers=T)
table(sceNew$outlier)
#2 cells dropped as outliers

outs<-colnames(sceNew)[sceNew$outlier]
### Add to metaData
metaData$pca.drop<-metaData$final.drop
metaData[outs,which(colnames(metaData)=="pca.drop")]<-TRUE

##### Color bad cells on PCA plot #####
colorDF<-as.data.frame(cbind(colnames(sceNew),"1"), stringsAsFactors=F)
rownames(colorDF)<-colorDF[,1]
colorDF[outs,2]<-"2"
colorDF[,2]<-as.factor(colorDF[,2])
tmp<-colorDF[,2,drop=F]

png(file=paste0(sampleFolder,"results/QC/3a_pca.png"),  width = 850, height = 642)
plotReducedDim(sceNew, use_dimred = "PCA_coldata", colour_by='outlier',shape_by='outlier') + labs(title="PCA with outliers colored")
dev.off()


#### Add to metaData table ####
pca.drop<-metaData[colnames(sce),"pca.drop"]
sum(pca.drop)

##### Create violinplots ####
##Before
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/3b_beforePcaFiltering_2.png"), colsToColor = c(rep('pca.drop',3)))

##After
toPlot<-metaData[! metaData$pca.drop,]
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/3c_afterPcaFiltering.png"), colsToColor = c(rep('pca.drop',3)))


##### Add to diagnostics #####
diagnostics[['pcaRemove']]<-0
diagnostics[['pcaRemove']]<-sum(pca.drop)
diagnostics[['totalRemove']]<-nrow(metaData)-ncol(sce)

##### Remove outlier cells #####
sce <- sce[,!(pca.drop)] 
dim(sce)
diagnostics[['dimAfterPCA']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")

### Remove some variables
rm(sceNew)

################################################################################
########## FINALIZE QC
################################################################################

dim(sce)
saveRDS(sce, file=paste0(sampleFolder,"Robjects/sce.rds"))
# sce <- readRDS(file=paste0(sampleFolder,"Robjects/sce.rds"))

rawDataFiltered<-rawData[rownames(sce),colnames(sce)]
dim(rawDataFiltered)
# 27998 11383
diagnostics[['dimBeforeSeuratObj']]<-paste0(nrow(rawDataFiltered)," genes - ",ncol(rawDataFiltered)," cells")

### Remove some variables
rm(sce)
rm(rawData)


################################################################################
########## CREATE SEURAT OBJECT
################################################################################

##### Create object #####
seuratObj <- CreateSeuratObject(counts = rawDataFiltered, project = "seuratObj", min.cells = 3, min.features = 200)

### Explore object
names(seuratObj)
# "RNA"
seuratObj@active.assay
# "RNA"

dim(seuratObj[['RNA']])
dim(seuratObj)
#16708 11383

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
VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"), nCol = 3)
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

### Get normalised values
GetAssayData(seuratObj, assay = "RNA", slot="data")[1:5,1:4]
seuratObj[['RNA']]@data[1:5,1:5]


##### Check per group #####
metaDataTable<-seuratObj@meta.data
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@data),2,function(x){sum(x>0)})

drawVlnPlotSeurat_split(metaDataTable, paste0(sampleFolder,"results/QC/5_afterNorm_splitted.png"))

################################################################################
########## GET HVG
################################################################################

seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(seuratObj))

### Get more info about HVGs (mean, dispersion and dispersion scaled)
head(HVFInfo(seuratObj))

### Plot variable features with and without labels
top10 <- head(x = VariableFeatures(object = seuratObj), 10)
plot1 <- VariableFeaturePlot(object = seuratObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(file=paste0(sampleFolder,"results/QC/6_hvg.png"), width = 850, height = 642)
CombinePlots(plots = list(plot1, plot2))
dev.off()


##### Add to diagnostics #####
diagnostics[['varGenes']]<-length(VariableFeatures(seuratObj))


################################################################################
########## SCALE DATA
################################################################################
# Apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
# Scaling data => Shifts the expression of each gene, so that the mean expression across cells is 0
# Scaling data => Scales the expression of each gene, so that the variance across cells is 1 (so that highly-expressed genes do not dominate)

seuratObj <- ScaleData(object = seuratObj, features = rownames(seuratObj))

### Get scaled values
GetAssayData(seuratObj, assay = "RNA", slot="scale.data")[1:5,1:4]
seuratObj[['RNA']]@scale.data[1:5,1:5]

##### Check per group #####
head(metaDataTable)
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@scale.data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@scale.data),2,function(x){sum(x>0)})

drawVlnPlotSeurat_split(metaDataTable, paste0(sampleFolder,"results/QC/7_afterScale_splitted.png"))

################################################################################
########## PCA
################################################################################
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), npcs = 50, ndims.print = 1:5, nfeatures.print = 10)

# print(x = seuratObj[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(object = seuratObj, dims = 1:2, reduction = "pca")

names(seuratObj)
# "RNA" "pca"
seuratObj[['pca']]@cell.embeddings[1:5,1:5]

########################################
########## PCA PLOT
########################################
pdf(file=paste0(sampleFolder,"results/QC/8a_PCA.pdf"), width = 10)
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,3))
dev.off()


########################################
########## PCA PLOT 3D
########################################
library("rgl")

expTable<-seuratObj[['RNA']]@data
matrixPCAtmp<-expTable[VariableFeatures(seuratObj),]

### Prepare PCA-plot
pca<-seuratObj[['pca']]@cell.embeddings
matrixPCA<-cbind(pca[,1],pca[,2],pca[,3])

PoV <- seuratObj[['pca']]@stdev^2/sum(seuratObj[['pca']]@stdev^2)
summary(pca[,c(1,2,3)])

### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",pch=2, type="s",radius=0.5, legend=TRUE, xlab=paste0("pc1 (",round(PoV[1]*100,2),"%)"), 
                ylab=paste0("pc2 (",round(PoV[2]*100,2),"%)"), zlab=paste0("pc3 (",round(PoV[3]*100,2),"%)"))


rgl.viewpoint(0, 5)
rgl.snapshot(paste0(sampleFolder,"results/QC/8b_PCA_view1.png"))
rgl.viewpoint(35, 5)
rgl.snapshot(paste0(sampleFolder,"results/QC/8b_PCA_view2.png"))



########################################
########## HEATMAP OF PCs
########################################

### Create heatmap of PC 1-40
pdf(file=paste0(sampleFolder,"results/QC/9a_selectPC.pdf"))
PCHeatmap(seuratObj, dims = 1:12, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 13:24, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 25:36, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 37:40, cells = 500, balanced = TRUE)
dev.off()

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

# ### Run Jackstraw and create plot
# seuratObj <- JackStraw(object = seuratObj, num.replicate = 100)
# seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
# 
# pdf(file=paste0(sampleFolder,"results/QC/9c_jackStrawPlot.pdf")) 
# JackStrawPlot(object = seuratObj, dims = 1:20)
# dev.off()

### Create PCElbowplot
png(file=paste0(sampleFolder,"results/QC/9b_selectPC.png"), width = 850, height = 642)
ElbowPlot(object = seuratObj, ndims = 40)
dev.off()


################################################################################
########## CLUSTER THE CELLS
################################################################################
dimsToTry<-c(10,15,20,25)
resToUse<-0.8

### Final
dimsToTry<-c(25)
resToUse<-0.8
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse

library(reticulate)
use_condaenv(condaenv = "snowflakes", conda = "~/anaconda3/bin/conda", required = TRUE)
reticulate::py_install(packages ='umap-learn') #Restart R session and reload packages https://github.com/satijalab/seurat/issues/1760

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse)
  tsnePlot<-DimPlot(seuratObj, reduction = "tsne", label=T, label.size = 8, pt.size = 2)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "tsne", label=F, group.by="orig.ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0(sampleFolder,"results/QC/10a_test_tSNE_sliced_",min(dimsToUse),"-",max(dimsToUse),".png"), height = 8, width = 20)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n.neighbors = 30) #, umap.method = 'umap-learn', metric = 'correlation'
  umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(sampleFolder,"results/QC/10b_test_UMAP_sliced_",min(dimsToUse),"-",max(dimsToUse),".png"), height = 8, width = 20)
  
}

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,".rds"))

##### Save sliced object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_sliced_",sampleName,".rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_sliced_",sampleName,".rds"))

##### Read objects
seuratObj<-readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_sliced_",sampleName,".rds"))
diagnostics<-readRDS(file=paste0(sampleFolder,"Robjects/diagnostics_sliced_",sampleName,".rds"))

##### Create new clusters
c("CPE","CPE","CPE","CPE","MF",'CPE','EC','CPE',"FB1",'MF',"VAC",
  "Xist_CPE",'NK','Doublets','Doublets_2','FB2',"Mito","Microglia-like",
  "DC1","Doublets_10","Doublets_8", "MC_Unclear_12", "DC2_16","Baso","Neuron")

##new cluster for part of cl10
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 < -4.5)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 10))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 19)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

##new cluster for part of cl8
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

tsneSlice<-tsneTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., tSNE_1 > -34)
wantedCells<-intersect(tsneSlice$cell, WhichCells(seuratObj, idents = 8))
colorSomeCells(clusterMatrix, tsneTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 20)
DimPlot(seuratObj, reduction = "tsne", label = T, label.size = 8)

##new cluster for part of cl12
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

tsneSlice<-tsneTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., tSNE_1 > -5, tSNE_2 > 18)
wantedCells<-intersect(tsneSlice$cell, WhichCells(seuratObj, idents = 12))
colorSomeCells(clusterMatrix, tsneTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 21)
DimPlot(seuratObj, reduction = "tsne", label = T, label.size = 8)

##new cluster for part of cl16
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

tsneSlice<-tsneTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., tSNE_1 > -20)
wantedCells<-intersect(tsneSlice$cell, WhichCells(seuratObj, idents = 16))
colorSomeCells(clusterMatrix, tsneTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 22)
DimPlot(seuratObj, reduction = "tsne", label = T, label.size = 8)

##new cluster for part of cl21
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

tsneSlice<-tsneTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., tSNE_1 > -1, tSNE_2 < 20)
wantedCells<-intersect(tsneSlice$cell, WhichCells(seuratObj, idents = 21))
colorSomeCells(clusterMatrix, tsneTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 23)
DimPlot(seuratObj, reduction = "tsne", label = T, label.size = 8)

##new cluster for part of cl11
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

tsneSlice<-tsneTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., tSNE_1 > -16, tSNE_1 < -12,)
wantedCells<-intersect(tsneSlice$cell, WhichCells(seuratObj, idents = 11))
colorSomeCells(clusterMatrix, tsneTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 24)
DimPlot(seuratObj, reduction = "tsne", label = T, label.size = 8)


##### Test gene
FeaturePlot(object = seuratObj, features = c("Mki67","Ccnb2","Cdk1","Stmn1"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98')

FeaturePlot(object = seuratObj, features = "Mki67", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98')
FeaturePlot(object = seuratObj, features = "Xcr1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98')

################################################################################
########## PLOTS
################################################################################
clusterMatrix<-seuratObj@meta.data
# logTable<-as.matrix(seuratObj[['RNA']]@data)
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

########## UMI plot ##########
p1<-drawUMI_mitoPlot(tsneTable, 'tsne', clusterMatrix, 'nCount_RNA',"Tsne")
p2<-drawUMI_mitoPlot(umapTable, 'umap', clusterMatrix, 'nCount_RNA',"Umap")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results/QC/11a_UMI.png"), width = 20)

########## mito.genes plot ##########
p1<-drawUMI_mitoPlot(tsneTable, 'tsne', clusterMatrix, 'percent.mito',"Tsne")
p2<-drawUMI_mitoPlot(umapTable, 'umap', clusterMatrix, 'percent.mito',"Umap")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results/QC/11b_percMito.png"), width = 20)


########## PCA plot ##########
pdf(file=paste0(sampleFolder,"results/QC/13a_PCA.pdf"), width=10)
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,3))
dev.off()

pdf(file=paste0(sampleFolder,"results/QC/13b_PCA_split.pdf"), width=10)
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,2), group.by = "orig.ident")
DimPlot(object = seuratObj, reduction = "pca", dims = c(2,3), group.by = "orig.ident")
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,3), group.by = "orig.ident")
dev.off()


###### Barplot aggregate #######
# create a dataset
Sample <- seuratObj@meta.data$orig.ident
cluster <- seuratObj@active.ident
Aggr <- rep(sampleName,length(cluster)) 

cluster<- factor(cluster,levels(cluster)[c(7:25,6,5,4,3,2,1)]) #reorder levels
levels(cluster) <-c(rep("Epithelial cells",4),"Macrophages","Epithelial cells","Endothelial cells",
                    "Epithelial cells","Fibroblasts Type 1",'Macrophages',"Vascular associated cells",
                    "Xist+ Epithelial cells",'NK cells','Doublets','Doublets_2','Fibroblasts Type 2',
                    "cDC2","Microglia-like Macrophages","cDC1","Doublets_3","Doublets_4", "Monocytes and Unclear", 
                    "Mitotic Cells","Basophils","Neuronal Cells")
data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

# Stacked
library("ggthemes")
png(file=paste0(sampleFolder,"results/QC/14_SampleDistribution_ggplot2_annotated.png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"results/QC/14_SampleDistribution_ggplot2_annotated_2.png"), width = 2000, height = 1800, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() + 
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png(file=paste0(sampleFolder,"results/QC/14_SampleDistribution_ggplot2_annotated_3.png"), width = 2000, height = 1800, res = 300)
ggplot(data2, aes(fill=cluster, y=Freq, x=Aggr)) + theme_bw() + 
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

### Create annotated UMAP
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident
seuratObj@meta.data$annotated_clusters<- factor(seuratObj@meta.data$annotated_clusters,levels(seuratObj@meta.data$annotated_clusters)[c(7:25,6,5,4,3,2,1)]) #reorder levels
levels(seuratObj@meta.data$annotated_clusters) <- c(rep("Epithelial cells",4),"Macrophages","Epithelial cells","Endothelial cells",
                                                    "Epithelial cells","Fibroblasts Type 1",'Macrophages',"Vascular associated cells",
                                                    "Xist+ Epithelial cells",'NK cells','Doublets','Doublets_2','Fibroblasts Type 2',
                                                    "cDC2","Microglia-like Macrophages","cDC1","Doublets_3","Doublets_4", "Monocytes and Unclear", 
                                                    "Mitotic Cells","Basophils","Neuronal Cells")
U1 <- DimPlot(seuratObj, reduction = "umap", label = T, repel = T, label.size = 4, group.by="annotated_clusters")
ggsave(U1, file=paste0(sampleFolder,"results/QC/15_Annotated_UMAP_sliced.png"), width =15, height= 12, dpi = "retina")

T1 <- DimPlot(seuratObj, reduction = "tsne", label = T, repel = T, label.size = 4, group.by="annotated_clusters")
ggsave(T1, file=paste0(sampleFolder,"results/QC/15_Annotated_TSNE_sliced.png"), width =15, height= 12, dpi = "retina")

### Subset UMAP
seuratObj_subset<-seuratObj
Idents(seuratObj_subset)<-seuratObj_subset@meta.data$annotated_clusters
seuratObj_subset<-subset(seuratObj_subset, idents = c("Epithelial cells","Macrophages","Endothelial cells",
                                                      "Fibroblasts Type 1","Vascular associated cells",
                                                      "Xist+ Epithelial cells",'NK cells','Fibroblasts Type 2',
                                                      "cDC2","Microglia-like Macrophages","cDC1","Monocytes and Unclear", 
                                                      "Mitotic Cells","Basophils","Neuronal Cells"))
seuratObj_subset@meta.data$annotated_clusters<-factor(seuratObj_subset@meta.data$annotated_clusters, levels=c("Epithelial cells","Xist+ Epithelial cells","Macrophages","Microglia-like Macrophages",
                                                                                                              "Endothelial cells","Vascular associated cells","Fibroblasts Type 1",'Fibroblasts Type 2',
                                                                                                              'NK cells',"cDC2","cDC1","Monocytes and Unclear","Basophils",
                                                                                                              "Mitotic Cells","Neuronal Cells"))
Idents(seuratObj_subset)<-seuratObj_subset@meta.data$annotated_clusters

U2<-DimPlot(seuratObj_subset, reduction = "umap", label = T, repel = T, label.size = 5, pt.size = 1)
ggsave(U2, file=paste0(sampleFolder,"results/QC/15_Annotated_UMAP_subset.png"), width =16, height= 12, dpi = "retina")

T2<-DimPlot(seuratObj_subset, reduction = "tsne", label = T, repel = T, label.size = 5, pt.size = 1)
ggsave(T2, file=paste0(sampleFolder,"results/QC/15_Annotated_TSNE_subset.png"), width =16, height= 12, dpi = "retina")


## New annotation January 2022 for Daan FB paper
seuratObj_subset@meta.data$annotated_clusters2<-seuratObj_subset@meta.data$annotated_clusters
levels(seuratObj_subset@meta.data$annotated_clusters2)<-c("Epithelial cells","Epithelial cells","Macrophages","Other Immune cells",
                                                          "Endothelial cells","Vascular associated cells","Stromal Fibroblasts",'Base Fibroblasts',
                                                          rep("Other Immune cells",6),"Neuronal Cells")
Idents(seuratObj_subset)<-seuratObj_subset@meta.data$annotated_clusters2

U3<-DimPlot(seuratObj_subset, reduction = "umap", label = T, repel = T, label.size = 5, pt.size = 1, group.by = "annotated_clusters2")
pdf(file=paste0(sampleFolder,"results/QC/17_Annotated_UMAP_subset_v2_",sampleName,".pdf"), width =16, height= 12)
U3
dev.off()

##############################################

##### Save paper object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_paper_",sampleName,".rds"))

##### Read paper object
seuratObj<-readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_paper_",sampleName,".rds"))

##############################################

# Look at ventricle distribution
D1<-DimPlot(seuratObj_subset, reduction = "umap", label = T, repel = T, label.size = 5, pt.size = 1, split.by = "orig.ident")
D2<-DimPlot(seuratObj_subset, reduction = "umap", label = F, pt.size = 1, group.by = "orig.ident")

# Stacked barplot
Sample <- seuratObj_subset@meta.data$orig.ident
cluster <- seuratObj_subset@active.ident
Aggr <- rep(sampleName,length(cluster)) 

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

library("ggthemes")
S1<-ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split2<-data[which(data$Sample=="LpsNegFourVentr"),]
data_split3<-data[which(data$Sample=="LpsNegLatVentr"),]
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})
data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-round((x/sum(data_split3$Freq))*100,2)})

data_new<-rbind(data_split2,data_split3)

S2<-ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")

# Combo figure
pdf(file=paste0(sampleFolder,"results/QC/17_Annotated_UMAP_v2_split_by_ventricle_",sampleName,".pdf"), width =16, height= 12)
D1
D2
S1
S2
dev.off()

pdf(file=paste0(sampleFolder,"results/QC/17_Annotated_UMAP_ventricle_annotation_",sampleName,".pdf"), width =16, height= 12)
D2
dev.off()

pdf(file=paste0(sampleFolder,"results/QC/17_SampleDistribution_ggplot2_annotated_v2_",sampleName,"_adjusted.pdf"), width =8, height= 6)
S2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #Remove grid lines
dev.off()

# Featureplots Dpep1 and Igfbp6
F1<-FeaturePlot(object = seuratObj_subset, features = "Dpep1", cols = c("yellow", "red"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=F)

F2<-FeaturePlot(object = seuratObj_subset, features = "Igfbp6", cols = c("yellow", "red"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=F)

pdf(file=paste0(sampleFolder,"results/QC/17_Featureplots_FB_markers_",sampleName,".pdf"), width =13, height= 12)
F1
F2
dev.off()


## Feature plot 04/08/21
F1<-FeaturePlot(object = seuratObj_subset, features = "Folr1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98')

ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plot_Folr1_LpsNegAggr.png"), height = 10, width = 15, dpi = "retina")

## Feature plot 01/09/21 Chirantan
F1<-FeaturePlot(object = seuratObj_subset, features = "Ttr", cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T)

ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plot_Ttr_LpsNegAggr.png"), height = 10, width = 15, dpi = "retina")

D1<-DimPlot(seuratObj_subset, reduction = "umap", label = T, repel = T, label.size = 5, pt.size = 1,
        group.by = "orig.ident")

ggsave(D1, file=paste0(sampleFolder,"results/QC/Feature_plot_Ttr_accompanying_dimplot.png"), height = 10, width = 15, dpi = "retina")

## Feature plot 07/12/21 Smpd3
F1<-FeaturePlot(object = seuratObj_subset, features = "Smpd3", cols = c("yellow", "red"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, label = T, repel = T, label.size = 3)

D1<-DotPlot(seuratObj_subset, features = "Smpd3")

ggsave(F1+D1, file=paste0(sampleFolder,"results/QC/Plots_Smpd3_",sampleName,".png"), height = 15, width = 15, dpi = "retina")

## Feature plot 11/02/22 Pigr (No expr)
FeaturePlot(object = seuratObj_subset, features = "Pigr", cols = c("yellow", "red"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, label = T, repel = T, label.size = 3)

## Feature plot 06/06/22 Ace2 (No expr)
F1<-FeaturePlot(object = seuratObj_subset, features = "Ace2", cols = c("yellow", "red"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, label = T, repel = T, label.size = 3)

pdf(file=paste0(sampleFolder,"results/QC/17_Featureplot_Ace2_",sampleName,".pdf"), width =15, height= 10)
F1
dev.off()

## Feature plot Laminins 05/12/22 

Features_Lam<-c("Lamc2","Lamc1","Lamb3","Lamc3","Lama5","Lama2","Lama4","Lamb2","Lamb1","Lama1","Lama3","Col1a1")

F1<-FeaturePlot(object = seuratObj_subset, features = "Ace2", cols = c("yellow", "red"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, label = T, repel = T, label.size = 3)

pdf(file=paste0(sampleFolder,"results/QC/17_Featureplot_Laminins_",sampleName,".pdf"), width =15, height= 10)
for (feature in Features_Lam) {
  F1<-FeaturePlot(object = seuratObj_subset, features = feature, cols = c("yellow", "red"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, label = T, repel = T, label.size = 3)
  print(F1)
}
dev.off()

Colors_dotplot<-c("#071AE5","#F50635") #030720

D1<-DotPlot(seuratObj_subset, features = Features_Lam, cols = Colors_dotplot) + RotatedAxis()

pdf(file = paste0(sampleFolder,"results/QC/17_Dotplot_Laminin_genes_",sampleName,".pdf"), width = 10, height = 12)
D1
dev.off()

########## Split plot advanced: per sample ##########
##### Preparation
nrOtherSamples<-length(listLabels)-1

##### Create plots
allPlots<-list()
for(reductionType in c('umap','tsne')){
  
  listPlots<-list()
  for(i in 1:length(listLabels)){
    print(paste0('Working on ',listLabels[[i]]))
    myOrder<-c(unlist(listLabels[-i]), listLabels[[i]])
    tmp<-clusterMatrix[order(match(clusterMatrix$orig.ident, myOrder)),]
    
    myOrder2<-c(listLabels[[i]], unlist(listLabels[-i]))
    highlight<-list("color"=c(listColors[[i]], rep("lightgray",nrOtherSamples)), "alpha"=c(1, rep(0.3,nrOtherSamples)))
    highlight<-lapply(highlight, setNames, nm = myOrder2)
    print(highlight)
    
    pSplit<-drawTSNEplot_colorSample(highlight, tmp, columnName="orig.ident", reductionType=reductionType)
    listPlots[[i]]<-pSplit
  }
  allPlots[[reductionType]]<-listPlots
}


##Write
pdf(paste0(sampleFolder,'results/QC/12a_splitAdvanced.pdf'), width=20)
for(j in 1:length(allPlots)){
  listPlots<-allPlots[[j]]
  
  for(i in seq(from=1,to=length(listPlots),by=2)){
    x<-i+1
    if(x<=length(listPlots)){
      grid.arrange(listPlots[[i]], listPlots[[x]] ,ncol=2)
    }else{
      grid.arrange(listPlots[[i]],ncol=2)
    }
    
  }
}
dev.off()

### Remove some variables
#rm(logTable)
rm(tsneTable)
rm(umapTable)

################################################################################
########## ANNOTATION
################################################################################

pdf(file=paste0(sampleFolder,"results/QC/15_annotation.pdf"), width = 15)
grid.arrange(umapPlot, umapPlotSplit, ncol=2)
FeaturePlot(object = seuratObj, features = c("Mki67","Ccnb2","Cdk1","Stmn1"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98')
dev.off()


clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

wantedCells<-WhichCells(seuratObj, idents = 12)
colorSomeCells(clusterMatrix, umapTable, wantedCells)

################################################################################
########## DOUBLET FINDER
################################################################################
# devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
library('DoubletFinder')
library('fields')
library('modes')

runDoubletFinder <- function(object, PCs, minPCT = 1, maxPCT = 10, pN = 0.25){

  if(!require('DoubletFinder')){
    devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
  }
  require("dplyr")
  require("stringr")

  nSamples <- str_split_fixed(rownames(object@meta.data), "-", 2)[,2] %>% unique() %>% length()

  DFPredictions <- c()

  for(i in 1:nSamples) {

    cellsToUse <- grep(paste0("-", i, "$"), rownames(object@meta.data), value = T)
    object_sub <- SubsetData(object, cells = cellsToUse)
    nDoubletsMin <- (length(colnames(object)) * (minPCT/100) ) %>% round()
    nDoubletsMax <- (length(colnames(object)) * (maxPCT/100) ) %>% round()

    findpK <- paramSweep_v3(object_sub, PCs = PCs) %>%
      summarizeSweep() %>%
      find.pK()

    maxScore <- findpK %>% pull('BCmetric') %>% which.max()
    pKValue <- findpK[maxScore, 'pK'] %>% as.character() %>% as.numeric()

    object_sub <- doubletFinder_v3(seu = object_sub, pN = pN, pK = pKValue, nExp = nDoubletsMin, reuse.pANN = F, PCs = PCs)
    object_sub <- doubletFinder_v3(seu = object_sub, pN = pN, pK = pKValue, nExp = nDoubletsMax, reuse.pANN = paste0("pANN_", pN, "_", pKValue, "_", nDoubletsMin), PCs = PCs)

    object_sub@meta.data$DFPrediction <- object_sub@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_", nDoubletsMax)]
    object_sub@meta.data$DFPrediction[ object_sub@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_", nDoubletsMin)] == 'Doublet'] <- "High Confidence"
    object_sub@meta.data$DFPrediction <- gsub('Doublet', "Low Confidence", object_sub@meta.data$DFPrediction)


    newData <- object_sub@meta.data$DFPrediction
    names(newData) <- rownames(object_sub@meta.data)

    DFPredictions <- c(DFPredictions, newData)
  }

  object@meta.data$DFPrediction <- DFPredictions[rownames(object@meta.data)]

  return(object)
}

seuratObj <- runDoubletFinder(seuratObj, 1:25) #<NUMBER OF PCs USED IN TSNE/UMAP/CLUSTERING>
DimPlot(object = seuratObj, reduction = 'tsne', group.by = 'DFPrediction', cols = c("red", "yellow", "#C9C9C9"))
DimPlot(object = seuratObj, reduction = 'umap', group.by = 'DFPrediction', cols = c("red", "yellow", "#C9C9C9"))

################################################################################
########## GET DE GENES
################################################################################

dir.create(paste0(sampleFolder,"results/QC/Feature_plots"))

c("CPE","CPE","CPE","CPE","MF",'CPE','EC','CPE',"FB1",'MF',"VAC",
  "Xist_CPE",'NK','Doublets','Doublets_2','FB2',"Mito","Microglia-like",
  "DC1","Doublets_10","Doublets_8", "MC_12", "DC2_16","Baso","Neuron")

##### Epithelial marker -> large cluster (0+1+2+3+5+7) + Xist 11 <-> doublets (13!!)
FeaturePlot(object = seuratObj, features = "Otx2", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Ttr", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1 <- FeaturePlot(object = seuratObj, features = c("Otx2", "Ttr"), cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Epithelial_Cells.png"), height = 10, width = 20, dpi = "retina")

##### Endothelial marker -> cluster 6 (10?? SLice??) <-> 14 doublets!!
FeaturePlot(object = seuratObj, features = "Pecam1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Flt1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = c("Plvap"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1<-FeaturePlot(object = seuratObj, features = c("Pecam1", "Flt1", "Plvap"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Endothelial_Cells.png"), height = 10, width = 20, dpi = "retina")

##### Vascular associated marker -> cluster 10
FeaturePlot(object = seuratObj, features = "Pdgfrb", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Mylk", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = c("Myh11"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = c("Tagln"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1<-FeaturePlot(object = seuratObj, features = c("Pdgfrb", "Mylk", "Myh11", "Tagln"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Vasc_Assoc_Cells.png"), height = 10, width = 20, dpi = "retina")

##### Fibroblast marker -> cluster 8 (1) + 15 (2) <-> doublets 10 left cluster + doublets tsne
FeaturePlot(object = seuratObj, features = "Dcn", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Col1a1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1<-FeaturePlot(object = seuratObj, features = c("Dcn", "Col1a1"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Fibroblasts.png"), height = 10, width = 20, dpi = "retina")

##### Macrophage marker -> cluster 4+9 + Doublets: 14
FeaturePlot(object = seuratObj, features = "Adgre1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Csf1r", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Fcgr1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1<-FeaturePlot(object = seuratObj, features = c("Adgre1", "Csf1r","Fcgr1"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Macrophages.png"), height = 10, width = 20, dpi = "retina")

##### T-cell marker -> 
F1<-FeaturePlot(object = seuratObj, features = c("Cd3d","Cd4","Il7r","Cd8a","Cd8b","Ccr7","Id3"), cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_T_Cells.png"), height = 10, width = 20, dpi = "retina")


##### B-cell marker -> 
F1<-FeaturePlot(object = seuratObj, features = c("Ms4a1","Cd79a","Cd19","Ebf1"), cols = c("grey", "blue"),
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_B_Cells.png"), height = 10, width = 20, dpi = "retina")


##### Microglia marker -> cluster 17 (macrophage-microglia like cells)
F1 <- FeaturePlot(object = seuratObj, features = "P2ry12", cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Microglia like.png"), height = 10, width = 10, dpi = "retina")


##### NK cell marker -> cluster 12
FeaturePlot(object = seuratObj, features = "Klrb1c", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Gzmb", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1<- FeaturePlot(object = seuratObj, features = c("Klrb1c", "Gzmb"), cols = c("grey", "blue"), 
                 reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_NK_Cells.png"), height = 10, width = 20, dpi = "retina")

##### cDC marker -> cluster 16?? (Cd209a and Ccr7+) + 18 (Xcr1+) 
FeaturePlot(object = seuratObj, features = "Cd209a", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Ccr7", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Xcr1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1<-FeaturePlot(object = seuratObj, features = c("Cd209a", "Ccr7","Xcr1"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_cDC.png"), height = 20, width = 20, dpi = 250)


##### Neutrophil marker -> 
FeaturePlot(object = seuratObj, features = "S100a8", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Ngp", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Retnlg", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1<-FeaturePlot(object = seuratObj, features = c("S100a8", "Ngp", "Retnlg"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_neutrophils.png"), height = 20, width = 20, dpi = 250)

##### Basophil marker -> ??
FeaturePlot(object = seuratObj, features = "Fcer1a", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Fcer1g", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = c("Fcer1a", "Fcer1g"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

##### Monocyte marker -> top part of cluster 12
FeaturePlot(object = seuratObj, features = "Plac8", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Ly6c2", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Ccr2", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1<-FeaturePlot(object = seuratObj, features = c("Plac8", "Ly6c2", "Ccr2"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_monocytes.png"), height = 20, width = 20, dpi = 250)

##### Mast cell marker -> cluster 16
F1<-FeaturePlot(object = seuratObj, features = "Mcemp1", cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_mast_cells.png"), height = 10, width = 10, dpi = 250)


##### Neuronal cell marker -> unknown
FeaturePlot(object = seuratObj, features = "Tubb3", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

##### Mitotic cell marker -> cluster 16 right
FeaturePlot(object = seuratObj, features = "Birc5", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Mki67", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F1<-FeaturePlot(object = seuratObj, features = c("Birc5", "Mki67"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_mitotic_cells.png"), height = 10, width = 20, dpi = 250)

##### Glial marker -> ??
FeaturePlot(object = seuratObj, features = "Slc1a3", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Fabp7", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
# FeaturePlot(object = seuratObj, features = "Olig1", cols = c("grey", "blue"), 
#             reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = c("Slc1a3", "Fabp7"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)


##### Dissociation effect -> especially immune cells 
FeaturePlot(object = seuratObj, features = "Fos", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Junb", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Atf3", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Dusp1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Ccl4", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

########################################
##### all clusters vs all clusters
########################################

### Find markers for every cluster compared to all remaining cells, report only the positive ones
allMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allMarkers$cluster)
saveRDS(allMarkers, file=paste0(sampleFolder,"Robjects/allMarkers_",sampleName,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,".rds"))
diagnostics[['markersPerCluster']]<-paste0(table(allMarkers$cluster)," markers for cluster ",rownames(table(allMarkers$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_",sampleName,".rds"))


### Create list with markers
totalNrClusters<-max(as.numeric(names(table(allMarkers$cluster))))
totalNrClustersPlusOne<-totalNrClusters+1
markersList<-list()

for(i in 1:totalNrClustersPlusOne){
  clusterNr<-i-1
  
  tmp<-allMarkers[allMarkers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  markersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(markersList)<-paste0("cluster",0:totalNrClusters)

### Write to Excel
library('openxlsx')
write.xlsx(markersList, file =paste0(sampleFolder, "results/allClusters_",sampleName,".xlsx"))

FeaturePlot(object = seuratObj, features = c("Gp9"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98')


########################################
##### some clusters vs some clusters
########################################

group1_cl9.14.4 <- FindMarkers(seuratObj, ident.1 = c(4,9,14), min.pct = 0.10, 
                              min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = TRUE)

group2_cl15.8 <- FindMarkers(seuratObj, ident.1 = c(8,15), min.pct = 0.10, 
                              min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = TRUE)

group3_cl6.10.13 <- FindMarkers(seuratObj, ident.1 = c(6,10,13), min.pct = 0.10, 
                                            min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = TRUE)

group4_cl12.16.17.18 <- FindMarkers(seuratObj, ident.1 = c(12,16,17,18), min.pct = 0.10, 
                                       min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = TRUE)

group5_cl0.1.2.3.5.7 <- FindMarkers(seuratObj, ident.1 = c(0,1,2,3,5,7), min.pct = 0.10, 
                              min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = TRUE)


##### Create list
listDEgenesGroups<-tibble::lst(group1_cl9.14.4, group2_cl15.8, group3_cl6.10.13,
                               group4_cl12.16.17.18, group5_cl0.1.2.3.5.7)

##Add geneSymbol in column (for the export)
listDEgenesGroups<-lapply(listDEgenesGroups,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenesGroups<-lapply(listDEgenesGroups, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDEgenesGroups<-lapply(listDEgenesGroups, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
listDEgenesGroups<-lapply(listDEgenesGroups,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDEgenesGroups,file=paste0(sampleFolder,"Robjects/groupClusters_",sampleName,".rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDEgenesGroups, paste0(sampleFolder,"results/groupClusters_",sampleName,".xlsx"))


##Test gene
FeaturePlot(object = seuratObj, features = c("Mki67"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98')


########################################
##### different markers
########################################
########## 1. CREATE NEW CLUSTERS ##########
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

### Create new clusters: split on source
seuratObjNew<-seuratObj
seuratObjNew@meta.data$newClustersTmp<-seuratObjNew@meta.data$annotated_clusters
seuratObjNew@meta.data$newClusters<-paste0(seuratObjNew@meta.data$newClustersTmp,"_",seuratObjNew@meta.data$orig.ident)
head(seuratObjNew@meta.data)

### Use the new clusters
Idents(seuratObjNew)<-seuratObjNew@meta.data$newClusters
DimPlot(seuratObjNew, reduction = "umap", label = T, label.size = 3)

########## 2. GET MARKERS ##########
getDEgenes<-function(ident1, ident2){
  markersDiff <- FindMarkers(seuratObjNew, ident.1 = ident1, ident.2 = ident2, 
                             min.pct = 0.10) # No min diff pct!! , min.diff.pct = 0.15 or logFC logfc.threshold = 0.30,
  markersDiff<-markersDiff[markersDiff$p_val_adj < 0.01,]
  markersDiff<-markersDiff[order(markersDiff$avg_logFC, decreasing = T),]
  
  markersDiff$geneSymbol<-rownames(markersDiff)
  markersDiff$pct.1<-markersDiff$pct.1+0.001
  markersDiff$pct.2<-markersDiff$pct.2+0.001
  
  markersDiff<-rbind(markersDiff[markersDiff$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/pct.2*avg_logFC),
                     markersDiff[markersDiff$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/pct.1*avg_logFC))
  markersDiff<-markersDiff[order(markersDiff$score, decreasing = T),]
  return(markersDiff)
}


##### Get diff markers #####
Epithelial_4VnegvsLVneg<-getDEgenes(c('Epithelial cells_LpsNegFourVentr'),
                                    c('Epithelial cells_LpsNegLatVentr'))
Epithelial_4VnegvsLVneg<-Epithelial_4VnegvsLVneg[order(Epithelial_4VnegvsLVneg$avg_logFC,decreasing = T),]
head(Epithelial_4VnegvsLVneg)
dim(Epithelial_4VnegvsLVneg)

Xist_Epithelial_4VnegvsLVneg<-getDEgenes(c('Xist+ Epithelial cells_LpsNegFourVentr'),
                                         c('Xist+ Epithelial cells_LpsNegLatVentr'))
Xist_Epithelial_4VnegvsLVneg<-Xist_Epithelial_4VnegvsLVneg[order(Xist_Epithelial_4VnegvsLVneg$avg_logFC,decreasing = T),]
head(Xist_Epithelial_4VnegvsLVneg)
dim(Xist_Epithelial_4VnegvsLVneg)

Endothelial_4VnegvsLVneg<-getDEgenes(c('Endothelial cells_LpsNegFourVentr'),
                                     c('Endothelial cells_LpsNegLatVentr'))
Endothelial_4VnegvsLVneg<-Endothelial_4VnegvsLVneg[order(Endothelial_4VnegvsLVneg$avg_logFC,decreasing = T),]
head(Endothelial_4VnegvsLVneg)
dim(Endothelial_4VnegvsLVneg)

Macrophages_4VnegvsLVneg<-getDEgenes(c('Macrophages_LpsNegFourVentr'),
                                     c('Macrophages_LpsNegLatVentr'))
Macrophages_4VnegvsLVneg<-Macrophages_4VnegvsLVneg[order(Macrophages_4VnegvsLVneg$avg_logFC,decreasing = T),]
head(Macrophages_4VnegvsLVneg)
dim(Macrophages_4VnegvsLVneg)

Fibroblasts1_4VnegvsLVneg<-getDEgenes(c('Fibroblasts Type 1_LpsNegFourVentr'),
                                     c('Fibroblasts Type 1_LpsNegLatVentr'))
Fibroblasts1_4VnegvsLVneg<-Fibroblasts1_4VnegvsLVneg[order(Fibroblasts1_4VnegvsLVneg$avg_logFC,decreasing = T),]
head(Fibroblasts1_4VnegvsLVneg)
dim(Fibroblasts1_4VnegvsLVneg)

## No diff markers
# Fibroblasts2_4VnegvsLVneg<-getDEgenes(c('FB2_LpsNegFourVentr'),
#                                       c('FB2_LpsNegLatVentr'))
# Fibroblasts2_4VnegvsLVneg<-Fibroblasts2_4VnegvsLVneg[order(Fibroblasts2_4VnegvsLVneg$avg_logFC,decreasing = T),]
# head(Fibroblasts2_4VnegvsLVneg)
# dim(Fibroblasts2_4VnegvsLVneg)

# VAC_4VnegvsLVneg<-getDEgenes(c('Vascular associated cells_LpsNegFourVentr'),
#                              c('Vascular associated cells_LpsNegLatVentr'))
# VAC_4VnegvsLVneg<-VAC_4VnegvsLVneg[order(VAC_4VnegvsLVneg$avg_logFC,decreasing = T),]
# head(VAC_4VnegvsLVneg)
# dim(VAC_4VnegvsLVneg)

# Microglia_Macrophage_4VnegvsLVneg<-getDEgenes(c('Microglia-like Macrophages_LpsNegFourVentr'),
#                                               c('Microglia-like Macrophages_LpsNegLatVentr'))
# Microglia_Macrophage_4VnegvsLVneg<-Microglia_Macrophage_4VnegvsLVneg[order(Microglia_Macrophage_4VnegvsLVneg$avg_logFC,decreasing = T),]
# head(Microglia_Macrophage_4VnegvsLVneg)
# dim(Microglia_Macrophage_4VnegvsLVneg)

# NK_Cells_4VnegvsLVneg<-getDEgenes(c('NK cells_LpsNegFourVentr'),
#                                   c('NK cells_LpsNegLatVentr'))
# NK_Cells_4VnegvsLVneg<-NK_Cells_4VnegvsLVneg[order(NK_Cells_4VnegvsLVneg$avg_logFC,decreasing = T),]
# head(NK_Cells_4VnegvsLVneg)
# dim(NK_Cells_4VnegvsLVneg)
# 
# DC1_Cells_4VnegvsLVneg<-getDEgenes(c('cDC1_LpsNegFourVentr'),
#                                   c('cDC1_LpsNegLatVentr'))
# DC1_Cells_4VnegvsLVneg<-DC1_Cells_4VnegvsLVneg[order(DC1_Cells_4VnegvsLVneg$avg_logFC,decreasing = T),]
# head(DC1_Cells_4VnegvsLVneg)
# dim(DC1_Cells_4VnegvsLVneg)
# 
# DC2_Cells_4VnegvsLVneg<-getDEgenes(c('cDC2_LpsNegFourVentr'),
#                                    c('cDC2_LpsNegLatVentr'))
# DC2_Cells_4VnegvsLVneg<-DC2_Cells_4VnegvsLVneg[order(DC2_Cells_4VnegvsLVneg$avg_logFC,decreasing = T),]
# head(DC2_Cells_4VnegvsLVneg)
# dim(DC2_Cells_4VnegvsLVneg)

# Mitotic_Cells_4VnegvsLVneg<-getDEgenes(c('Mitotic Cells_LpsNegFourVentr'),
#                                      c('Mitotic Cells_LpsNegLatVentr'))
# Mitotic_Cells_4VnegvsLVneg<-Mitotic_Cells_4VnegvsLVneg[order(Mitotic_Cells_4VnegvsLVneg$avg_logFC,decreasing = T),]
# head(Mitotic_Cells_4VnegvsLVneg)
# dim(Mitotic_Cells_4VnegvsLVneg)

##add to list
listDiffMarkers<-tibble::lst(Epithelial_4VnegvsLVneg, Macrophages_4VnegvsLVneg, Endothelial_4VnegvsLVneg, Fibroblasts1_4VnegvsLVneg)

lapply(listDiffMarkers, dim)

saveRDS(listDiffMarkers,file=paste0(sampleFolder,"Robjects/markersDiffSamples.rds"))

## Read DiffMarkers
listDiffMarkers<- readRDS(file=paste0(sampleFolder,"Robjects/markersDiffSamples.rds"))

### Write to Excel
library('openxlsx')
write.xlsx(listDiffMarkers, file = paste0(sampleFolder,"results/QC/summaryDiffMarkers.xlsx"))


### DE genes plot
Gene_vector_Up <- vector()
Gene_vector_Down <- vector()

for (x in 1:length(listDiffMarkers)) {
  test <- data.frame(listDiffMarkers[x])
  Gene_vector_Up <- c(Gene_vector_Up, nrow(test[which(test[2] > 0),]))
  Gene_vector_Down <- c(Gene_vector_Down, nrow(test[which(test[2] < 0),]))
}
Gene_vector <- c(Gene_vector_Up, Gene_vector_Down)
Cells<-as.factor(c('Epithelial Cells', "Macrophages", "Endothelial Cells", "Fibroblasts Type 1"))
Genes <- data.frame(c=rep(Cells, 2), g=Gene_vector, Regulation=c(rep("UP", length(Cells)), rep("DOWN",length(Cells))))

g1 <- ggplot(data = Genes, 
             mapping = aes(x = c, 
                           y = ifelse(test = Regulation == "DOWN",  yes = -g, no = g), 
                           fill = Regulation)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x=c, y = ifelse(test = Regulation == "DOWN",  yes = -g, no = g), label = Genes$g, 
                hjust = ifelse(test = Regulation == "DOWN",  yes = 1, no = 0))) +
  scale_y_continuous(breaks=seq(-500,500,100),labels=abs(seq(-500,500,100))) +
  labs(y = "Gene counts", x="Cell populations", title = "DE genes between 4V Lps Neg and LV Lps Neg") + 
  coord_flip()

ggsave(g1, file=paste0(sampleFolder,"results/QC/16_DE_gene_plot.png"), width = 20, dpi = "retina")

#New plot Daan
library(RColorBrewer)
library(ggrepel)
Colorset<-brewer.pal(12,"Set3")
Titles<-as.factor(c('Epithelial Cells', "Macrophages", "Endothelial Cells", "Fibroblasts Type 1"))
myList=list()
for (i in 1:length(listDiffMarkers)) {
  testData<-listDiffMarkers[[i]]
  testData<-testData[,c(2,6:7)]
  testData[,"Jitter"]<-runif(length(testData$geneSymbol)) #Add random numbers for x axis uniform distribution (Jitter to see dots)
  testData<- rbind(testData[testData$score>0,] %>% dplyr::mutate(.,log_score=log10(score)),
                   testData[testData$score<0,] %>% dplyr::mutate(.,log_score=-log10(-score)))
  rownames(testData)<-testData$geneSymbol
  testData<-testData[-2]
  
  p <- ggplot(testData, aes(x=Jitter, y=avg_logFC)) + 
    geom_point(color=Colorset[i]) +
    coord_cartesian(ylim = c(-3, 6)) +
    geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-0.29,ymax=0.29),alpha=0.01, fill="lightblue")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(title=as.character(Titles[i]))
  
  myList[[i]]<-p
}

grid.arrange(grobs = myList, ncol = length(Titles), nrow = 1)



for (i in 1:length(listDiffMarkers)) {
  testData<-listDiffMarkers[[i]]
  testData<-testData[,c(2,6:7)]
  testData[,"Jitter"]<-runif(length(testData$geneSymbol)) #Add random numbers for x axis uniform distribution (Jitter to see dots)
  testData<- rbind(testData[testData$score>0,] %>% dplyr::mutate(.,log_score=log10(score)),
                   testData[testData$score<0,] %>% dplyr::mutate(.,log_score=-log10(-score)))
  
  p <- ggplot(testData, aes(x=Jitter, y=avg_logFC, label = geneSymbol)) + 
    geom_point(color=Colorset[i]) +
    coord_cartesian(ylim = c(-2, 2)) +
    # geom_hline(yintercept = 0, color = "lightblue", size = 10) +
    geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-0.29,ymax=0.29),alpha=0.01,fill="lightblue")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(title=as.character(Titles[i])) +
    ### geom_text_repel
    # only label players with PTS > 25 or < 18
    # align text vertically with nudge_y and allow the labels to
    # move horizontally with direction = "x"
    geom_label_repel(data          = rbind(subset(testData, avg_logFC > 3), subset(testData, avg_logFC < -2)),
                     # nudge_y       = 6 - subset(testData, avg_logFC > 3)$avg_logFC,
                     size          = 4,
                     box.padding   = 1.5,
                     point.padding = 0.5,
                     force         = 100,
                     segment.size  = 0.2,
                     segment.color = "grey50")
  
  myList[[i]]<-p
}

ggsave(grid.arrange(grobs = myList, ncol = length(Titles), nrow = 1), file=paste0(sampleFolder,"results/QC/16_DE_gene_plot_fancy.png"), height= 8, width = 12, dpi = 300)



drawMultipleFeaturePlot<-function(groupName, markers, ncol, reductionType){
  myList=list()
  for(i in 1:length(markers)){
    p<-drawFeaturePlot(groupName,markers[i], reductionType)
    myList[[2]]<-p2
  }
  return(grid.arrange(grobs = myList, ncol = 2, nrow = 1))
}


################################################################################
########## PRINT DIAGNOSTICS
################################################################################
printDiagnostics(sampleName = sampleName, fileName = paste0(sampleFolder,'results/diagnostics_',sampleName,'.txt'), diagnostics = diagnostics)



################################################################################
########## EXPORT FOR TSNE TOOL
################################################################################
exportShiny<-DisneyTools::expShiny(seuratObj, conditions = unlist(listLabels), assay = "RNA")
saveRDS(exportShiny, file=paste0(sampleFolder,"results/QC/neededFilesOnlineTool/SO_",sampleName,"_040320.rds"))
annotation <- exportShiny@meta.data$annotated_clusters %>% as.factor() 
saveRDS(annotation, file=paste0(sampleFolder,"results/QC/neededFilesOnlineTool/SO_",sampleName,"_annotation_040320.rds"))

################################################################################
########## EXPORT FOR SCATTER TOOL
################################################################################
dir.create(paste0(sampleFolder,"results/neededFilesOnlineTool"))

########## UMAP plot ########## 
png(file=paste0(sampleFolder,"results/neededFilesOnlineTool/UMAP_",sampleName,".png"), width = 850, height = 620)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
dev.off()

png(file=paste0(sampleFolder,"results/neededFilesOnlineTool/UMAP_",sampleName,"_split.png"), width = 850, height = 620)
DimPlot(seuratObj, reduction = "umap", label = F, group.by="orig.ident")
dev.off()

########## Feature plot ##########
myColPal <- colorRampPalette(c('gray','blue'))

logTable<-as.matrix(seuratObj[['RNA']]@data)
logTableResult<-logTable
for(i in 1:nrow(logTable)){
  tmp<-myColPal(100)[as.numeric(cut(logTable[i,],breaks = 100))]
  logTableResult[i,]<-tmp
}

### Create data frame and reorder columns
logTableResult<-as.data.frame(logTableResult,stringsAsFactors = FALSE)
logTableResult$Gene<-rownames(logTableResult)
maxNrColumn<-ncol(logTableResult)
maxNrColumnMinusOne<-maxNrColumn-1
logTableFinal<-logTableResult[,c(maxNrColumn, 1:maxNrColumnMinusOne)]

tsneTableFinal<-umapTable
colnames(tsneTableFinal)<-c("tSNE_1","tSNE_2")
tsneTableFinal$Sample<-rownames(umapTable)
tsneTableFinal<-tsneTableFinal[,c(3,1:2)]

write.table(tsneTableFinal,file=paste0(sampleFolder,"results/neededFilesOnlineTool/tsneTable_",sampleName,".csv"),sep=",",row.names = FALSE)
write.table(logTableFinal,file=paste0(sampleFolder,"results/neededFilesOnlineTool/logTable_",sampleName,".csv"),sep=",",row.names = FALSE)


allGenes<-as.matrix(rownames(logTable))
write.table(allGenes,file=paste0(sampleFolder,"results/neededFilesOnlineTool/allGenes_",sampleName,".txt"),sep=",",
            row.names = FALSE, col.names = FALSE)



