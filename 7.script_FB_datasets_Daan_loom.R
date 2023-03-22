## Work with loom files single cell atlas

# devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

library(loomR)
library(Matrix)
library(Seurat)
library(scran)
library(scater)

setwd("/home/clintdn/VIB/DATA/Roos/Daan 1/FB_datasets/Atlas_Zeisel/")

source('/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/RAW_DATA/script_functions_COVID.R') #KEVIN
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')

# ## Vascular taxonomy level 1
# lfile <- connect(filename = "l6_r1_vascular_cells.loom", mode = "r+")
# MM<-lfile$matrix[, ]
# 
# cells<-lfile$col.attrs$CellID[]
# genes<-lfile$row.attrs$Gene[]
# 
# rownames(MM) = cells
# colnames(MM) = genes
# 
# mat<-t(MM)
# 
# lfile$close_all()
# 
# ## Vascular taxonomy level 2
# lfile2 <- connect(filename = "l6_r2_vascular_cells.loom", mode = "r+")
# MM2<-lfile2$matrix[, ]
# 
# cells2<-lfile2$col.attrs$CellID[]
# genes2<-lfile2$row.attrs$Gene[]
# 
# rownames(MM2) = cells2
# colnames(MM2) = genes2
# 
# mat2<-t(MM2)
# 
# lfile2$close_all()

## Vascular taxonomy level 3
lfile3 <- connect(filename = "l6_r3_vascular_cells.loom", mode = "r+")
MM3<-lfile3$matrix[, ]

cells3<-lfile3$col.attrs$CellID[]
genes3<-lfile3$row.attrs$Gene[]

rownames(MM3) = cells3
colnames(MM3) = genes3

# Pull metadata from the column attributes: issue with duplicated cells!!
# metadata3 <- lfile3$get.attribute.df(MARGIN = 2)
# names(lfile3$col.attrs)
clusters3<-lfile3$col.attrs$ClusterName[]
n_genes<-lfile3$col.attrs$`_NGenes`[]
Tot_genes<-lfile3$col.attrs$`_Total`[]
Mito_ribo_ratio<-lfile3$col.attrs$MitoRiboRatio[]

metadata3<-data.frame(clusters3[!duplicated(cells3)])
rownames(metadata3)<-cells3[!duplicated(cells3)]

metadata3$nGene<-n_genes[!duplicated(cells3)]
metadata3$TotGene<-Tot_genes[!duplicated(cells3)]
metadata3$Mito_ribo_ratio<-Mito_ribo_ratio[!duplicated(cells3)]

mat3<-t(MM3)

lfile3$close_all()

# ## All vascular loom files have exactly the same dimensions and values!!
# ## Work with taxonomy level 3!!
# 
# rm(cells)
# rm(genes)
# rm(mat)
# rm(MM)
# 
# rm(cells2)
# rm(genes2)
# rm(mat2)
# rm(MM2)
# 
# gc()

################################

# ## Astroependymal taxonomy level 3
# lfile4 <- connect(filename = "l6_r3_astroependymal_cells.loom", mode = "r+")
# MM4<-lfile4$matrix[, ]
# 
# cells4<-lfile4$col.attrs$CellID[]
# genes4<-lfile4$row.attrs$Gene[]
# 
# # Pull metadata from the column attributes
# # metadata4 <- lfile4$get.attribute.df(MARGIN = 2 )
# clusters4<-lfile4$col.attrs$ClusterName[]
# 
# rownames(MM4) = cells4
# colnames(MM4) = genes4
# 
# mat4<-t(MM4)
# 
# lfile4$close_all()


## Ependymal taxonomy level 4
lfile5 <- connect(filename = "l6_r4_ependymal_cells.loom", mode = "r+")
MM5<-lfile5$matrix[, ]

cells5<-lfile5$col.attrs$CellID[]
genes5<-lfile5$row.attrs$Gene[]

rownames(MM5) = cells5
colnames(MM5) = genes5

# Pull metadata from the column attributes
metadata5 <- lfile5$get.attribute.df(MARGIN = 2 )

mat5<-t(MM5)

lfile5$close_all()

# ## Only focus on ependymal. Issue with astroependymal. Can't retrieve metadata (duplicated)
# rm(cells4)
# rm(genes4)
# rm(mat4)
# rm(MM4)
# 
# gc()

#################################

## Check matrices
mat3[1:5,1:5]
mat5[1:5,1:5]

## Convert to sce
sce3<-SingleCellExperiment(list(counts=mat3))
sce5<-SingleCellExperiment(list(counts=mat5))

## Normalisation
set.seed(123)
gc()
q.clust3<- quickCluster(sce3)
q.clust5<- quickCluster(sce5)

# table(q.clust)
sce3 <- scran::computeSumFactors(sce3, cluster=q.clust3)
sce5 <- scran::computeSumFactors(sce5, cluster=q.clust5)
# if you get a warning of negative size factors, too many low quality cells
sce3 <- scater::logNormCounts(sce3)
sce5 <- scater::logNormCounts(sce5)

## Remove duplicated genes and cells for conversion to seurat
sce3<-sce3[!duplicated(rownames(sce3)),!duplicated(colnames(sce3))]
sce5<-sce5[!duplicated(rownames(sce5)),!duplicated(colnames(sce5))]

## Convert to seurat
experiment3 <- "Vascular_cells_atlas"
experiment5 <- "Ependymal_cells_atlas"

seuratObj3 <- as.Seurat(sce3, assay = "RNA", counts = "counts", data = "logcounts", project = experiment3) #Take along the normalized data!
seuratObj5 <- as.Seurat(sce5, assay = "RNA", counts = "counts", data = "logcounts", project = experiment5) #Take along the normalized data!

dim(seuratObj3)
#27933  12144

dim(seuratObj5)
#27933  1257

## Add metadata
seuratObj3@meta.data<-metadata3
seuratObj5@meta.data<-metadata5

##### Save object
dir.create("Robjects")

saveRDS(seuratObj3, file=paste0("Robjects/seuratObj_",experiment3,".rds"))
saveRDS(seuratObj5, file=paste0("Robjects/seuratObj_",experiment5,".rds"))

##################################################################
##################################################################
##################################################################

## Get extra metadata !!
## Convert to sce
sce3<-SingleCellExperiment(list(counts=mat3))
sce5<-SingleCellExperiment(list(counts=mat5))

## Remove duplicated genes and cells for conversion to seurat
sce3<-sce3[!duplicated(rownames(sce3)),!duplicated(colnames(sce3))]
sce5<-sce5[!duplicated(rownames(sce5)),!duplicated(colnames(sce5))]

###############################################################################################
## Analyze sce object 3

##### Get mitochondrial genes #####
is.mito <- grepl("^MT-", rownames(sce3), ignore.case = TRUE)
sum(is.mito)
##13
rownames(sce3)[is.mito]

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
is.rbc <- rownames(sce3) %in% rbc.genes

#### Get COVID genes
covid.genes <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
is.covid <- rownames(sce3) %in% covid.genes

##### Calculate QC metrics #####
### => pData(sce3) is created
sce3<- addPerCellQC(sce3, subsets=list(Mito=is.mito, RBC=is.rbc, COVID=is.covid))
dim(colData(sce3))
# colnames(colData(sce3))

### List samples
listLabels<-"Vascular_cells_atlas"


##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("staticNr"=colnames(sce3),"orig.ident"=listLabels[[1]], "nGene"=sce3$detected,"nUMI"=sce3$sum,
                     "percent.mito"=sce3$subsets_Mito_percent, "percent.rbc"=sce3$subsets_RBC_percent,
                     "percent.COVID"=sce3$subsets_COVID_percent,
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1

table(metaData$orig.ident)

## Add to metadata3
metadata3_New<-cbind(metadata3,metaData)

## Add to seuratobject3
seuratObj3@meta.data<-metadata3_New

##########################################################################################################

## Analyze sce object 5

##### Get mitochondrial genes #####
is.mito <- grepl("^MT-", rownames(sce5), ignore.case = TRUE)
sum(is.mito)
##13
rownames(sce5)[is.mito]

#### Get RBC genes ####
is.rbc <- rownames(sce5) %in% rbc.genes

#### Get COVID genes
is.covid <- rownames(sce5) %in% covid.genes

##### Calculate QC metrics #####
### => pData(sce5) is created
sce5<- addPerCellQC(sce5, subsets=list(Mito=is.mito, RBC=is.rbc, COVID=is.covid))
dim(colData(sce5))
# colnames(colData(sce5))

### List samples
listLabels<-"Ependymal_cells_atlas"

##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("staticNr"=colnames(sce5),"orig.ident"=listLabels[[1]], "nGene"=sce5$detected,"nUMI"=sce5$sum,
                     "percent.mito"=sce5$subsets_Mito_percent, "percent.rbc"=sce5$subsets_RBC_percent,
                     "percent.COVID"=sce5$subsets_COVID_percent,
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1

table(metaData$orig.ident)

## Add to metadata3
metadata5_New<-cbind(metadata5,metaData)

## Add to seuratobject3
seuratObj5@meta.data<-metadata5_New


################################################################################
########## Save object
################################################################################

##### Save object
saveRDS(seuratObj3, file=paste0("Robjects/seuratObj_",experiment3,".rds"))
saveRDS(seuratObj5, file=paste0("Robjects/seuratObj_",experiment5,".rds"))

