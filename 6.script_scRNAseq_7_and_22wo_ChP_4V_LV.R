## Script for processing 7- and 22-wo ChP 4V and LV samples from our lab
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

sampleName<-"Urvb_datasets"
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

### Change colnames
colnames(rawData_1)<-paste0(colnames(rawData_1),"-1")
colnames(rawData_2)<-paste0(colnames(rawData_2),"-2")
colnames(rawData_3)<-paste0(colnames(rawData_3),"-3")
colnames(rawData_4)<-paste0(colnames(rawData_4),"-4")

### Check before merge
cbind(head(rownames(rawData_1)),head(rownames(rawData_2)),head(rownames(rawData_3)),head(rownames(rawData_4)))
cbind(tail(rownames(rawData_1)),tail(rownames(rawData_2)),tail(rownames(rawData_3)),tail(rownames(rawData_4)))

dim(rawData_1)
dim(rawData_2)
dim(rawData_3)
dim(rawData_4)

sum(ncol(rawData_1),ncol(rawData_2),ncol(rawData_3),ncol(rawData_4))
paste0("Nr of cells: ",c(ncol(rawData_1),ncol(rawData_2),ncol(rawData_3),ncol(rawData_4)))

### Do merge
rawDataRNA<-cbind(rawData_1, rawData_2,rawData_3,rawData_4)
dim(rawDataRNA)
diagnostics[['dimBeforeSeuratObj']]<-paste0(nrow(rawDataRNA)," genes - ",ncol(rawDataRNA)," cells")

### Remove some variables
rm(rawData_1)
rm(rawData_2)
rm(rawData_3)
rm(rawData_4)

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

#### Get RBC genes ####
# working with a list detected in 3 different samples
#rbc.genes <- c("ADIPOR1", "ALAS2", "BNIP3L", "HBA1", "HBB", "HBD", "HEMGN", "SNCA")
rbc.genes <- c("ADIPOR1", "ALAS2", "ATP5E", "BAG1", "BCL2L1", "BNIP3L",
               "BPGM", "BTF3", "CA1", "DCAF12", "EPB42", "FBXO7", "FKBP8",
               "FTL", "GSPT1", "GUK1", "GYPC", "HBA1","HBB", "HBD", "HEMGN",
               "MPP1", "MYL6", "NCOA4", "OAZ1", "PFDN5", "RNF10", "RPL12",
               "RPL21", "RPL27A", "RPL30", "RPL31", "RPL32", "RPL38", "RPL41",
               "RPL7", "RPLP1", "RPLP2", "RPS11", "RPS12", "RPS13", "RPS14",
               "RPS24", "SELENBP1", "SERF2", "SLC25A37", "SLC25A39", "SNCA",
               "TPT1", "UBA52", "UBB", "YBX3")
is.rbc <- rownames(sce) %in% rbc.genes

#### Get COVID genes
covid.genes <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
is.covid <- rownames(sce) %in% covid.genes

##### Calculate QC metrics #####
### => pData(sce) is created
sce<- addPerCellQC(sce, subsets=list(Mito=is.mito, RBC=is.rbc, COVID=is.covid))
dim(colData(sce))
# colnames(colData(sce))

### List samples
listLabels<-c("RVD1_LpsNegFour","RVD2_LpsNegLat","RVD5_Y4V","RVD6_YLV")

##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("staticNr"=colnames(rawDataRNA),"orig.ident"=listLabels[[1]], "nGene"=sce$detected,"nUMI"=sce$sum,
                     "percent.mito"=sce$subsets_Mito_percent, "percent.rbc"=sce$subsets_RBC_percent,
                     "percent.COVID"=sce$subsets_COVID_percent,
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

##############################
########## FILTERING
##############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ < MANUAL BIT?
# start with 5 as more lenient thresholding or just do the low thresholding
# and rely on doubletfinder for the high ranges?
nmad_low_feature<-3
nmad_high_feature<-3

nmad_low_UMI<-3
nmad_high_UMI<-3

nmad_high_mito<-3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ > MANUAL BIT?

##### Aim: remove cells with low library sizes, low numbers of expressed features and with high mitochondrial proportions
##same as nGene in Seurat pipeline
feature.drop.low <- isOutlier(sce$detected, nmads=nmad_low_feature, type="lower", log=TRUE,batch = batch)
sum(feature.drop.low)

feature.drop.high <- isOutlier(sce$detected, nmads=nmad_high_feature, type="higher", log=TRUE,batch = batch)
sum(feature.drop.high)

feature.drop<-as.logical(feature.drop.low + feature.drop.high,batch = batch)
sum(feature.drop)

##same as UMI in Seurat pipeline
libsize.drop.low <- isOutlier(sce$sum, nmads=nmad_low_UMI, type="lower", log=TRUE,batch = batch)
sum(libsize.drop.low)

libsize.drop.high <- isOutlier(sce$sum, nmads=nmad_high_UMI, type="higher", log=TRUE,batch = batch)
sum(libsize.drop.high)

libsize.drop<-as.logical(libsize.drop.low+libsize.drop.high)
sum(libsize.drop)

##% mitochondrial genes
mito.drop.high <- isOutlier(sce$subsets_Mito_percent, nmads=nmad_high_mito, type="higher",batch = batch)
sum(mito.drop.high)

mito.drop<-as.logical(mito.drop.high)
sum(mito.drop)

##### add to metaData matrix #####
metaData$nGene.drop=feature.drop
metaData$nUMI.drop=libsize.drop
metaData$mito.drop=mito.drop
metaData$final.drop=feature.drop | libsize.drop | mito.drop
metaData$batch <- batch

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
diagnostics[['mito.drop.high']]<-sum(mito.drop.high)
diagnostics[['mito.drop']]<-sum(mito.drop)

########################################
########## PLOTS
########################################
palette(c("#00BFC4","#F8766D","#7CAE00","#C77CFF"))
###F8766D=red
###00BFC4=cyan
###7CAE00=green
###C77CFF=purple

toPlot<-metaData

##nGene
png(file=paste0(output.dir,"results/QC/01a_nGene.png"), width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nGene),]
hist(tmp$nGene, breaks=30)
theColors<-as.factor(tmp$nGene.drop)
barplot(tmp$nGene, col=theColors, border=theColors)
dev.off()

##nUMI
png(file=paste0(output.dir,"results/QC/01a_nUMI.png"), width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nUMI),]
hist(tmp$nUMI, breaks=30)
theColors<-as.factor(tmp$nUMI.drop)
barplot(tmp$nUMI, col=theColors, border=theColors)
dev.off()

##percent.mito
png(file=paste0(output.dir,"results/QC/01a_percMito.png"), width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$percent.mito),]
hist(tmp$percent.mito, breaks=30)
theColors<-as.factor(tmp$mito.drop)
barplot(tmp$percent.mito, col=theColors, border=theColors)
dev.off()

# calculating some limits for the plots
# ul.mito <- median(metaData$percent.mito)+nmad_high_mito*mad(metaData$percent.mito,na.rm = TRUE)
# ll.cd <- median(metaData$nUMI)-nmad_low_UMI*mad(metaData$nUMI,na.rm = TRUE)
# ul.cd <- median(metaData$nUMI)+nmad_high_UMI*mad(metaData$nUMI,na.rm = TRUE)
# ll.ng <- median(metaData$nGene)-nmad_low_feature*mad(metaData$nGene,na.rm = TRUE)
# ul.ng <- median(metaData$nGene)+nmad_high_feature*mad(metaData$nGene,na.rm = TRUE)

ul.mito <- attr(mito.drop.high,'thresholds')[2]
ll.cd <- attr(libsize.drop.low,'thresholds')[1]
ul.cd <- attr(libsize.drop.high,'thresholds')[2]
ll.ng <- attr(feature.drop.low,'thresholds')[1]
ul.ng <- attr(feature.drop.high,'thresholds')[2]

png(file=paste0(output.dir,"results/QC/01b_histogram_countdepth.png"), width=850)
ggplot(metaData, aes(x = nUMI)) +
  geom_histogram(binwidth=100) +
  xlab("Count depth (total UMI count)") +
  ylab("Frequency") +
  geom_vline(xintercept=ul.cd,
             color = "red", size=1) +
  geom_vline(xintercept=ll.cd,
             color = "red", size=1) +
  facet_grid(.~batch) +
  theme_bw() # histogram of count depth
dev.off()

png(file=paste0(output.dir,"results/QC/01b_histogram_genes.png"), width=850)
ggplot(metaData, aes(x = nGene)) +
  geom_histogram(binwidth=20) +
  xlab("Number of Genes") +
  ylab("Frequency") +
  geom_vline(xintercept=ul.ng,
             color = "red", size=1) +
  geom_vline(xintercept=ll.ng,
             color = "red", size=1) +
  facet_grid(.~batch) +
  theme_bw() # histogram of nr of genes
dev.off()

png(file=paste0(output.dir,"results/QC/01b_histogram_percentageMito.png"), width=850)
ggplot(metaData, aes(x = percent.mito)) +
  geom_histogram(binwidth=0.1) +
  xlab("% Mitochondrial counts") +
  ylab("Frequency") +
  geom_vline(xintercept=ul.mito,
             color = "red", size=1) +
  facet_grid(.~batch) +
  theme_bw() # histogram of % mito
dev.off()

png(file=paste0(output.dir,"results/QC/01c_scatterplot_filtering.png"), width=850)
ggplot(metaData, aes(x = nUMI,y=nGene,colour=percent.mito)) +
  geom_point(size=0.5) +
  scale_color_gradient2(midpoint=ul.mito, low="black", mid="white",
                        high="red", space ="Lab" )+
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Genes") +
  geom_vline(xintercept=ul.cd,
             color = "red", size=1) +
  geom_vline(xintercept=ll.cd,
             color = "red", size=1) +
  geom_hline(yintercept=ul.ng,
             color = "red", size=1) +
  geom_hline(yintercept=ll.ng,
             color = "red", size=1) +
  geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
  facet_grid(.~batch) +
  theme_bw()
dev.off()

########################################
########## Create violinPlots
########################################

toPlot<-metaData
drawVlnPlot_out(toPlot, fileName = paste0(output.dir,"results/QC/02a_Outliers_nGene.png"), colsToColor = c('nGene.drop','nGene.drop','nGene.drop','nGene.drop','nGene.drop'))
drawVlnPlot_out(toPlot, fileName = paste0(output.dir,"results/QC/02a_Outliers_nUMI.png"), colsToColor = c('nUMI.drop','nUMI.drop','nUMI.drop','nUMI.drop','nUMI.drop'))
drawVlnPlot_out(toPlot, fileName = paste0(output.dir,"results/QC/02a_Outliers_percMito.png"), colsToColor = c('mito.drop','mito.drop','mito.drop','mito.drop','mito.drop'))

GeneUmiMitodrops<-sum(metaData$final.drop)
toPlot<-metaData

#filename can be added to write the output to png file
drawVlnPlot_color(toPlot,paste0(output.dir,"results/QC/02a_Outliers_colors_base.png"))
drawVlnPlot_color_nGene(toPlot,paste0(output.dir,"results/QC/02a_Outliers_colors_nGene.png"))
drawVlnPlot_color_nGene_mito(toPlot,paste0(output.dir,"results/QC/02a_Outliers_colors_nGene&mito.png"))
drawVlnPlot_color_mito_extra(toPlot,paste0(output.dir,"results/QC/02a_Outliers_colors_Extra_mito.png"))
drawVlnPlot_color_nGene_extra(toPlot,paste0(output.dir,"results/QC/02a_Outliers_colors_Extra_nGene.png"))
drawVlnPlot_color_nUMI_extra(toPlot,paste0(output.dir,"results/QC/02a_Outliers_colors_Extra_nUMI.png"))

### Before filtering
toPlot<-metaData
drawVlnPlot(toPlot, fileName = paste0(output.dir,"results/QC/02b_beforeFiltering.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop',"percent.rbc","percent.COVID"))

### After filtering
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot(toPlot, fileName = paste0(output.dir,"results/QC/02c_afterFiltering.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop',"percent.rbc","percent.COVID"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ < MANUAL BIT
# CHECK OUTPUT BEFORE REMOVING OUTLIERS!!


########## Remove outliers
########################################
sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
# dim(sce)
### Number of cells removed
# nrow(metaData)-ncol(sce)
##### Add to diagnostics #####
diagnostics[['firstRemove']]<-nrow(metaData)-ncol(sce)
diagnostics[['dimFirstRemove']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ > MANUAL BIT

########################################
########## Create PCA
########################################

##(check via raw code of the function runPCA)
varsToUse <- c("sum","detected", "subsets_Mito_percent")
setdiff(colnames(colData(sce)),varsToUse)
exprs_to_plot <- scale(colData(sce)[,varsToUse], scale = T)
x.mad <- apply(exprs_to_plot, 2, mad)
x.mad[x.mad==0]
varsToUse<-setdiff(varsToUse, names(x.mad[x.mad==0]))
sceNew<-runColDataPCA(sce, outliers=T, variables = varsToUse)

##### Detect bad cells #####
#sceNew<-runPCA(sce,use_coldata=T, detect_outliers=T)
# table(sceNew$outlier)

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

png(file=paste0(output.dir,"results/QC/03a_PCA.png"),  width = 850, height = 642)
plotReducedDim(sceNew, dimred = "PCA_coldata", colour_by='outlier',shape_by='outlier') + labs(title="PCA with outliers colored")
dev.off()


#### Add to metaData table ####
pca.drop<-metaData[colnames(sce),"pca.drop"]
# sum(pca.drop)

##### Create violinplots ####
##Before
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot_out(toPlot, fileName = paste0(output.dir,"results/QC/03b_beforePcaFiltering.png"), colsToColor = c(rep('pca.drop',5)))

##After
toPlot<-metaData[! metaData$pca.drop,]
drawVlnPlot_out(toPlot, fileName = paste0(output.dir,"results/QC/03c_afterPcaFiltering.png"), colsToColor = c(rep('pca.drop',5)))


# ##### Add to diagnostics #####
# diagnostics[['pcaRemove']]<-0
# diagnostics[['pcaRemove']]<-sum(pca.drop)
# diagnostics[['totalRemove']]<-nrow(metaData)-ncol(sce)
# #
# # ##### Remove outlier cells #####
# sce <- sce[,!(pca.drop)]
# dim(sce)
# diagnostics[['dimAfterPCA']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")
#
# ### Remove some variables
rm(sceNew)



################################################################################
########## FINALIZE QC
################################################################################

dim(sce)
# saveRDS(sce, file=paste0(sampleFolder,"Robjects/sce.rds"))
# sce <- readRDS(file=paste0(sampleFolder,"Robjects/sce.rds"))

rawDataFiltered<-rawDataRNA[rownames(sce),colnames(sce)]
dim(rawDataFiltered)
# 27998  28646
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
#15937  28646

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

### Get normalised values
# GetAssayData(seuratObj, assay = "RNA", slot="data")[1:5,1:4]
# seuratObj[['RNA']]@data[1:5,1:5]

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
