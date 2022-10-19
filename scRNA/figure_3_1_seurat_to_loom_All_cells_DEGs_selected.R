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

#step3 :filiter matrix for target genes or DEGs
merge_data0<-read.table(file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/DEG_01q_default_merge.txt",header =T,sep="\t") 
head(merge_data0)
#merge_data_rm0<-merge_data0[which(!(merge_data0$cluster %in% c("Ery_b","Mix_Ery_b","Mks_b","Ery_p","Mix_Ery_p"))),]
length(unique(merge_data0$gene))# 392
exprMat0<-as.data.frame(exprMat0)
exprMat<-exprMat0[unique(as.character(merge_data0$gene)),]
exprMat0<-as.matrix(exprMat);dim(exprMat0)#2859 49988

#step4. Kept genes in cisTarget database ##基因过滤/选择
##按每个基因的reads总数进行过滤。该filter旨在去除最可能是噪音的基因。默认情况下，它（minCountsPerGene）保留所有样品中至少带有6个UMI reads的基因（例如，如果在1％的细胞中以3的值表达，则基因将具有的总数）。
##通过基因的细胞数来实现过滤（例如 UMI > 0 ，或log 2（TPM）> 1 ）。默认情况下(minSamples)，保留下来的基因能在至少1％的细胞中检测得到。
##最后，只保留RcisTarget数据库中可用的基因。

genesKept <- geneFiltering(exprMat=exprMat0, scenicOptions=scenicOptions, minCountsPerGene = 1,minSamples = 1)
#Number of genes left after applying the following filters (sequential):
#2858	genes with counts per gene > 1
#2858	genes detected in more than 1 cells
#2775	genes available in RcisTarget database
length(genesKept)# 2775
setdiff(unique(merge_data0$gene), genesKept)
exprMat_filtered <- exprMat0[genesKept, ]
dim(exprMat_filtered)## 370 53232

##The previous bug you encountered is now fixed in SCopeLoomR version 0.3.1 so no need anymore to convert the sparse matrix with as.matrix.
##https://github.com/aertslab/SCopeLoomR/issues/3
#step4  tramslated expression matrix into loom
source("/mnt/data/chenwei/gongchen/0.script/add_cellAnnotation.R")
loom <- build_loom("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/s1_ALLcell_DEGs_selected.loom", dgem=exprMat_filtered)
loom <- add_cellAnnotation(loom, cellInfo)
close_loom(loom)
