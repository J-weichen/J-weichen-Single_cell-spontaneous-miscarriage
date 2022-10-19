rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
options(stringsAsFactors = FALSE)
#Load packages
# To export/visualize in http://scope.aertslab.org
dyn.load('/mnt/data/chenwei/software/hdf/hdf5/lib/libhdf5.so.200')
dyn.load('/mnt/data/chenwei/software/hdf/hdf5/lib/libhdf5_hl.so.200.0.1')
#devtools::install_github("hhoeflin/hdf5r")
#devtools::install_github("aertslab/SCopeLoomR")
#devtools::install_github("aertslab/SCENIC")
library(hdf5r)
library(SCopeLoomR)
library(Seurat)
library(SCENIC)
library(arrow)

#
#library(DoubletFinder)
#library(monocle)
#step1 ：read Seurat object and load initial coldata and matrix

target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
head(target@meta.data)
target$final_major_subgroup_brief_Treat <- paste(target$final_major_subgroup_brief, target$Treat, sep = "_")
unique(target$final_major_subgroup_brief_Treat)
DefaultAssay(object = target)#"SCT"
DefaultAssay(target) <- "SCT"

ALL_coldata<-target@meta.data
exprMat0<-GetAssayData(object =target,slot = "counts")
cellInfo <- ALL_coldata

#step2 SCENIC init for filter gene  cisTarget database
#org="hgnc";dbs <- defaultDbNames[[org]]
dbs<-c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather","hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
scenicOptions <- initializeScenic(
  org="hgnc", # 物种名 or mgi, or dmel
  dbDir="/mnt/data/chenwei/10X_placenta/10x_data/191125-merge_nonormalization_add/SCENIC/human_hg38", # RcisTarget databases location
  dbs=dbs, # file name of motif database
  datasetTitle="SCENIC_on_jianhai_Cell_Atlas", # choose a name for your analysis
  nCores=10
)

#step3 : Exclude genes missing from database:
## following script was used to perform filter action usually using the following sentence
##genesKept <- geneFiltering(exprMat=dge, scenicOptions=scenicOptions, minCountsPerGene = 1,minSamples = 1)
## we do it because of the need of matrix as input for geneFiltering,but the RAM always show no enough for a matrix.

# Calculate stats
str(exprMat0);str(as.matrix(exprMat0))
nCountsPerGene <- rowSums(as.matrix(exprMat0), na.rm = T)
nCellsPerGene <- rowSums(as.matrix(exprMat0)>0, na.rm = T)
head(nCountsPerGene);head(nCellsPerGene)

## First filter:: # minCountsPerGene <- 3*.01*ncol(exprMat)
minCountsPerGene<-1
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene >  minCountsPerGene)]
# Second filter:: # minSamples <- ncol(exprMat)*.01
minSamples<-1
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells)## 25084

library(RcisTarget)
dbFilePath1 <- getDatabases(scenicOptions)[[1]]
motifRankings1 <- importRankings(dbFilePath1)
genesInDatabase1 <- colnames(getRanking(motifRankings1))

dbFilePath2 <- getDatabases(scenicOptions)[[2]]
motifRankings2 <- importRankings(dbFilePath2)
genesInDatabase2 <- colnames(getRanking(motifRankings2))
genesInDatabase<-unique(c(genesInDatabase1,genesInDatabase2));length(genesInDatabase)##27091

genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)#17447
genesKept <- genesLeft_minCells_inDatabases
exprMat_filtered <- exprMat0[genesKept, ]
dim(exprMat_filtered)##17447 49988

##The previous bug you encountered is now fixed in SCopeLoomR version 0.3.1 so no need anymore to convert the sparse matrix with as.matrix.
##https://github.com/aertslab/SCopeLoomR/issues/3
#step4  tramslated expression matrix into loom
source("/mnt/data/chenwei/gongchen/0.script/add_cellAnnotation.R")
loom <- build_loom("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/s1_ALLcell_all_genes.loom", dgem=exprMat_filtered)
loom <- add_cellAnnotation(loom, cellInfo)
close_loom(loom)