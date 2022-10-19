###Comparison analysis of multiple datasets using CellChat
#ref:https://www.jianshu.com/p/03f9c2f3ee4f
##https://www.jianshu.com/p/b3d26ac51c5a 
##important note concept
##Q1: "Interaction weights/strength" calculated:
###A1: It is the communication probability defined in our paper. The matrix is cellchat@net$prob or cellchat@netP$prob
##Q2: Interaction strength and probability values #182 map or not map to PPI 
#######https://github.com/sqjin/CellChat/issues/182
##ou can take a look at the updated tutorial, where you may want to change the methods for computing average expression if some important signaling genes only express in a small portion of cells.

#26 April, 2021
rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

#Create a directory to save figures
data.dir <- '/mnt/data/chenwei/gongchen/7.cellchat/CTRL'
dir.create(data.dir)
setwd(data.dir)

##Load the required libraries
library(reticulate)
use_python('/mnt/data/chenwei/software/miniconda2_new/envs/scenic/bin/',required = T)
library(CellChat)
library(scales)
library(ggsci)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(NMF)
library(VennDiagram)

future::plan("multiprocess", workers =20) # do parallel

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB_human_interaction_meta<- data.frame(CellChatDB$interaction)[,c("pathway_name","interaction_name","ligand","receptor")]
pathways_all <- unique(CellChatDB_human_interaction_meta$pathway_name)
pathyway_gene<-list();pathyway_CPIpair<-list();pathyway_identify<-c()
for (i in 1:length(pathways_all)) {
  pathways_name<-pathways_all[i]
  #pathways_name<-"BMP"
  pairs_name<-CellChatDB_human_interaction_meta[which(CellChatDB_human_interaction_meta$pathway_name ==pathways_name ),]$interaction_name 
  pairs_genes<-unique(as.character(unlist(lapply(strsplit(pairs_name,"_"), function(x) x))))
  pathyway_CPIpair<-c(pathyway_CPIpair,list(pairs_name))
  pathyway_gene<-c(pathyway_gene,list(pairs_genes))
  pathyway_identify<-c(pathyway_identify,pathways_name)
}
length(pathyway_CPIpair);length(pathyway_gene);length(pathyway_identify)
names(pathyway_CPIpair)<-pathyway_identify
names(pathyway_gene)<-pathyway_identify
saveRDS(pathyway_CPIpair, file = "/mnt/data/chenwei/gongchen/7.cellchat/cellchat_pathyway_CPIpair.rds")
saveRDS(pathyway_gene, file = "/mnt/data/chenwei/gongchen/7.cellchat/cellchat_pathyway_gene.rds")

#library(Seurat)
##plan(strategy = "multicore", workers = 42)
##options(future.globals.maxSize = 1500 * 1024^12)
#loading final seurat object
##target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
##subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts")
##target_final2<-subset(x = target_final, subset = re_annotation_TFac %in%  subgroup_order0)
##target_final2$re_annotation_TFac<- factor(target_final2$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
##DefaultAssay(target_final2) <- "RNA"
##target_final2 <- NormalizeData(target_final2, verbose = FALSE)
##data.input <- GetAssayData(target_final2, assay = "RNA", slot = "data") # normalized data matrix
##meta = target_final2@meta.data # a dataframe with rownames containing cell mata data
##saveRDS(data.input, file = "/mnt/data/chenwei/gongchen/7.cellchat/data.input.rds")
##saveRDS(meta, file = "/mnt/data/chenwei/gongchen/7.cellchat/meta.rds")

data.input <- readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/data.input.rds")
meta <- readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/meta.rds")

##for CTRL
cell.CTRL = rownames(meta)[meta$Treat == "CTRL"] # extract the cell names from disease data
# Prepare input data for CelChat analysis
ctrl_data.input = data.input[, cell.CTRL]
ctrl_meta = meta[cell.CTRL, ]
#Create a CellChat object
ctrl_cellchat <- createCellChat(object = ctrl_data.input, meta = ctrl_meta, group.by = "re_annotation_TFac")


#CellChatDB_human_interaction_meta[which(),]
# use a subset of CellChatDB for cell-cell communication analysis
####CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
ctrl_cellchat@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
ctrl_cellchat <- subsetData(ctrl_cellchat) # This step is necessary even if using the whole database
ctrl_cellchat <- identifyOverExpressedGenes(ctrl_cellchat)#####time consume
ctrl_cellchat <- identifyOverExpressedInteractions(ctrl_cellchat)
# project gene expression data onto PPI network (optional)
ctrl_cellchat <- projectData(ctrl_cellchat, PPI.human)
#saveRDS(ctrl_cellchat, file = "/mnt/data/chenwei/gongchen/7.cellchat/ctrl_cellchat.rds")
#ctrl_cellchat<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/ctrl_cellchat.rds")

#Part II: Inference of cell-cell communication network
##Compute the communication probability and infer cellular communication network
ctrl_cellchat <- computeCommunProb(ctrl_cellchat)#####time consume

##Filter out the cell-cell communication if there are only few number of cells in certain cell groups
ctrl_cellchat <- filterCommunication(ctrl_cellchat, min.cells = 10)

##Extract the inferred cellular communication network as a data frame
##df.net <- subsetCommunication(cellchat)
##df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#Infer the cell-cell communication at a signaling pathway level
ctrl_cellchat <- computeCommunProbPathway(ctrl_cellchat)
#Calculate the aggregated cell-cell communication network
##We can calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. 
#USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
ctrl_cellchat <- aggregateNet(ctrl_cellchat)

### Systems analysis of cell-cell communication network
#1）It can determine major signaling sources and targets as well as mediators and influencers within a given signaling network using centrality measures from network analysis
#2）It can predict key incoming and outgoing signals for specific cell types as well as coordinated responses among different cell types by leveraging pattern recognition approaches.
#3）It can group signaling pathways by defining similarity measures and performing manifold learning from both functional and topological perspectives.
#4）It can delineate conserved and context-specific signaling pathways by joint manifold learning of multiple networks
##Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
# Compute the network centrality scores
ctrl_cellchat <- netAnalysis_computeCentrality(ctrl_cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

##Manifold and classification learning analysis of signaling networks
#CellChat 能够量化所有重要信号通路之间的相似性，然后根据其CellChat 网络的相似性对其进行分组。分组可以基于功能或结构相似性进行。
##功能相似性分析要求两个数据集之间的细胞群组成相同。
##结构相似性：结构相似性用于比较其信号网络结构，而不考虑发送器和接收器的相似性。
##Identify signaling groups based on their functional similarity
ctrl_cellchat <- computeNetSimilarity(ctrl_cellchat, type = "functional")
ctrl_cellchat <- netEmbedding(ctrl_cellchat, type = "functional")
ctrl_cellchat <- netClustering(ctrl_cellchat, type = "functional")

#Identify signaling groups based on structure similarity
ctrl_cellchat <- computeNetSimilarity(ctrl_cellchat, type = "structural")
ctrl_cellchat <- netEmbedding(ctrl_cellchat, type = "structural")
ctrl_cellchat <- netClustering(ctrl_cellchat, type = "structural")

#Part V: Save the CellChat object
saveRDS(ctrl_cellchat, file = "/mnt/data/chenwei/gongchen/7.cellchat/CTRL/ctrl_cellchat.rds")

###for Abortion
# Prepare input data for CelChat analysis
cell.Abortion = rownames(meta)[meta$Treat == "Abortion"] # extract the cell names from disease data
Abortion_data.input = data.input[, cell.Abortion]
Abortion_meta = meta[cell.Abortion, ]
#Create a CellChat object
Abortion_cellchat <- createCellChat(object = Abortion_data.input, meta = Abortion_meta, group.by = "re_annotation_TFac")
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
####CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
Abortion_cellchat@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
Abortion_cellchat <- subsetData(Abortion_cellchat) # This step is necessary even if using the whole database
Abortion_cellchat <- identifyOverExpressedGenes(Abortion_cellchat)#####time consume
Abortion_cellchat <- identifyOverExpressedInteractions(Abortion_cellchat)
# project gene expression data onto PPI network (optional)
Abortion_cellchat <- projectData(Abortion_cellchat, PPI.human)
#saveRDS(Abortion_cellchat, file = "/mnt/data/chenwei/gongchen/7.cellchat/Abortion/Abortion_cellchat")
Abortion_cellchat<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/Abortion/Abortion_cellchat")
##Compute the communication probability and infer cellular communication network
Abortion_cellchat <- computeCommunProb(Abortion_cellchat)#####time consume

##Filter out the cell-cell communication if there are only few number of cells in certain cell groups
Abortion_cellchat <- filterCommunication(Abortion_cellchat, min.cells = 10)
Abortion_cellchat <- computeCommunProbPathway(Abortion_cellchat)
#Calculate the aggregated cell-cell communication network
Abortion_cellchat <- aggregateNet(Abortion_cellchat)

### Systems analysis of cell-cell communication network
# Compute the network centrality scores
Abortion_cellchat <- netAnalysis_computeCentrality(Abortion_cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#Abortion_cellchat<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/Abortion/Abortion_cellchat.rds")
#source("/mnt/data/chenwei/gongchen/7.cellchat/functional.R")
##Manifold and classification learning analysis of signaling networks
##Identify signaling groups based on their functional similarity
Abortion_cellchat <- computeNetSimilarity(Abortion_cellchat, type = "functional")
Abortion_cellchat <- netEmbedding(Abortion_cellchat, type = "functional")
#UserWarning: A few of your vertices were disconnected from the manifold.  This shouldn't cause problems.
## f"A few of your vertices were disconnected from the manifold.
source("/mnt/data/chenwei/gongchen/7.cellchat/netcluster_function.R")
Abortion_cellchat <- netClustering(Abortion_cellchat, type = "functional")
#Error in do_one(nmeth) : NA/NaN/Inf in foreign function call (arg 1)

netVisual_embedding(Abortion_cellchat, type = "functional", label.size = 3.5)

#Identify signaling groups based on structure similarity
Abortion_cellchat <- computeNetSimilarity(Abortion_cellchat, type = "structural")
Abortion_cellchat <- netEmbedding(Abortion_cellchat, type = "structural")
Abortion_cellchat <- netClustering(Abortion_cellchat, type = "structural")
netVisual_embedding(Abortion_cellchat, type = "structural", label.size = 3.5)

#Part V: Save the CellChat object
saveRDS(Abortion_cellchat, file = "/mnt/data/chenwei/gongchen/7.cellchat/Abortion/Abortion_cellchat.rds")

###for combined all cells
# Prepare input data for CelChat analysis
cell.combine = rownames(meta) # extract the cell names from disease data
combine_data.input = data.input
combine_meta = meta[cell.combine, ]
#Create a CellChat object
combine_cellchat <- createCellChat(object = combine_data.input, meta = combine_meta, group.by = "re_annotation_TFac")
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
####CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
combine_cellchat@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
combine_cellchat <- subsetData(combine_cellchat) # This step is necessary even if using the whole database
combine_cellchat <- identifyOverExpressedGenes(combine_cellchat)#####time consume
combine_cellchat <- identifyOverExpressedInteractions(combine_cellchat)
# project gene expression data onto PPI network (optional)
combine_cellchat <- projectData(combine_cellchat, PPI.human)
#saveRDS(combine_cellchat, file = "/mnt/data/chenwei/gongchen/7.cellchat/combine/combine_cellchat.rds")
combine_cellchat<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/combine/combine_cellchat.rds")
##Compute the communication probability and infer cellular communication network
combine_cellchat <- computeCommunProb(combine_cellchat)#####time consume

##Filter out the cell-cell communication if there are only few number of cells in certain cell groups
combine_cellchat <- filterCommunication(combine_cellchat, min.cells = 10)
combine_cellchat <- computeCommunProbPathway(combine_cellchat)
#Calculate the aggregated cell-cell communication network
combine_cellchat <- aggregateNet(combine_cellchat)

### Systems analysis of cell-cell communication network
# Compute the network centrality scores
combine_cellchat <- netAnalysis_computeCentrality(combine_cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#combine_cellchat<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/combine/combine_cellchat.rds")
##Manifold and classification learning analysis of signaling networks
##Identify signaling groups based on their functional similarity
combine_cellchat <- computeNetSimilarity(combine_cellchat, type = "functional")
combine_cellchat <- netEmbedding(combine_cellchat, type = "functional")
## f"A few of your vertices were disconnected from the manifold.
combine_cellchat <- netClustering(combine_cellchat, type = "functional")
netVisual_embedding(combine_cellchat, type = "functional", label.size = 3.5)

#Identify signaling groups based on structure similarity
combine_cellchat <- computeNetSimilarity(combine_cellchat, type = "structural")
combine_cellchat <- netEmbedding(combine_cellchat, type = "structural")
combine_cellchat <- netClustering(combine_cellchat, type = "structural")
netVisual_embedding(combine_cellchat, type = "structural", label.size = 3.5)

#Part V: Save the CellChat object
saveRDS(combine_cellchat, file = "/mnt/data/chenwei/gongchen/7.cellchat/combine/combine_cellchat.rds")