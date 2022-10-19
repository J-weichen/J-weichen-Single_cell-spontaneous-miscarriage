##ref:https://www.jianshu.com/p/49a0a0b50987
##https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

##Load the required libraries
library(reticulate)
use_python('/mnt/data/chenwei/software/miniconda2_new/envs/scenic/bin/',required = T)
library(CellChat)
library(patchwork)
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
library(ggsci)
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

#Create a directory to save figures
data.dir <- '/mnt/data/chenwei/gongchen/7.cellchat/Merge'
#dir.create(data.dir)
setwd(data.dir)

#Load CellChat object of each dataset and then merge together
cellchat.CTRL <- readRDS("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/ctrl_cellchat.rds")
cellchat.Abortion <- readRDS( "/mnt/data/chenwei/gongchen/7.cellchat/Abortion/Abortion_cellchat.rds")
object.list <- list(CTRL = cellchat.CTRL, Abortion = cellchat.Abortion)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

##Part II: Identify the conserved and context-specific signaling pathways
##CellChat performs joint manifold learning and classification of the inferred communication networks based on their functional and topological similarity. 
##NB: Such analysis is applicable to more than two datasets.
#Identify signaling groups based on their functional similarity
## NB: Functional similarity analysis is not applicable to multiple datsets with different cell type composition.
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

##Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")##time consume
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_Dim_plot_functional_similarity.pdf",width=20,height=20)
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 3)
###Error in base::colSums(x, na.rm = na.rm, dims = dims, ...) :  'x' must be an array of at least two dimensions

# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

##Part V: Compare the signaling gene expression distribution between different datasets
##We can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CTRL", "Abortion")) # set factor level
#Save the merged CellChat object
saveRDS(cellchat, file = "/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_cellchat.rds")
#sessionInfo()

#
cellchat<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_cellchat.rds")


##Part I: Predict general principles of cell-cell communication
##
##1)Whether the cell-cell communication is enhanced or not
##2)The interaction between which cell types is significantly changed
##3)How the major sources and targets change from one condition to another

##The differential network analysis only works for pairwise datasets. 
##If there are more datasets for comparison, we can directly show the number of interactions or interaction strength between any two cell populations in each dataset.
##To better control the node size and edge weights of the inferred networks across different datasets, 
##we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_cellchat_network_circle_number_for_abortion_CTRL_split__plot.pdf",width=19,height=9)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F,#color.use=my_morandi_colors,
                  vertex.weight = as.numeric(table(object.list[[i]]@idents)), vertex.size.max = 10,
                   edge.weight.max = weight.max[2], edge.width.max = 12, arrow.width = 0.1,arrow.size = 0.1, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_cellchat_network_circle_strength_for_USM_CTRL_split_plot.pdf",width=19,height=9)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F,#color.use=my_morandi_colors,
                   vertex.weight = as.numeric(table(object.list[[i]]@idents)), vertex.size.max = 10,
                   edge.weight.max = weight.max[2], edge.width.max = 12, arrow.width = 0.1,arrow.size = 0.1, title.name = paste0("strength of interactions - ", names(object.list)[i]))
}
dev.off()


#Compare the total number of interactions and interaction strength
#CellChat compares the the total number of interactions and interaction strength of the inferred cell-cell communication networks from different biological conditions.
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_barplot_total_number_strength_of_interactions_cellchat_network.pdf",width=6,height=6)
gg1 + gg2
dev.off()

#Compare the number of interactions and interaction strength among different cell populations
##circle plot, where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.

groupSize_all <- as.numeric(table(cellchat@idents$joint))
groupSize_USM <- as.numeric(table(cellchat@idents$Abortion))
groupSize_CTRL <- as.numeric(table(cellchat@idents$CTRL))
groupSize_subtract<-groupSize_USM-groupSize_CTRL
color_sbstract <-ifelse(groupSize_subtract>0,"red","blue")

groupSize_USM2 <- as.numeric(table(cellchat@idents$Abortion))/sum(as.numeric(table(cellchat@idents$Abortion)))*100
groupSize_CTRL2 <- as.numeric(table(cellchat@idents$CTRL))/sum(as.numeric(table(cellchat@idents$CTRL)))*100
groupSize_subtract2<-groupSize_USM2-groupSize_CTRL2

pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_cellchat_network_circle_change_number_strength_plot_all_cell_number.pdf",width=13,height=6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat,vertex.weight = groupSize_all,vertex.size.max = 20,arrow.width = 0.1,arrow.size = 0.1,   weight.scale = T)
netVisual_diffInteraction(cellchat,vertex.weight = groupSize_all,vertex.size.max = 20,arrow.width = 0.1,arrow.size = 0.1,   weight.scale = T, measure = "weight")
dev.off()

pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_cellchat_network_circle_change_number_strength_plot_cell_abs_number_change.pdf",width=13,height=6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat,color.use =color_sbstract,vertex.weight = abs(groupSize_subtract),vertex.size.max = 20,arrow.width = 0.1,arrow.size = 0.1,   weight.scale = T)
netVisual_diffInteraction(cellchat,color.use =color_sbstract,vertex.weight = abs(groupSize_subtract),vertex.size.max = 20,arrow.width = 0.1,arrow.size = 0.1,   weight.scale = T, measure = "weight")
dev.off()

pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_cellchat_network_circle_change_number_strength_plot_ratio_change.pdf",width=13,height=6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat,color.use =color_sbstract,vertex.weight = abs(groupSize_subtract2),vertex.size.max = 20,arrow.width = 0.1,arrow.size = 0.1,   weight.scale = T)
netVisual_diffInteraction(cellchat,color.use =color_sbstract,vertex.weight = abs(groupSize_subtract2),vertex.size.max = 20,arrow.width = 0.1,arrow.size = 0.1,   weight.scale = T, measure = "weight")
dev.off()

##differential number of interactions or interaction strength in a greater details using a heatmap. 
##The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). 
##The right colored bar plot represents the sum of row of values (outgoing signaling). 
##In the colorbar, red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_change_CPI_number_strength_Heatmap_plot.pdf",width=12,height=8)
gg1 + gg2
dev.off()
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_change_CPI_number_strength_Heatmap_plot_number.pdf",width=8,height=8)
gg1
dev.off()
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_change_CPI_number_strength_Heatmap_plot_strength.pdf",width=8,height=8)
gg2
dev.off()
##Differential number of interactions or interaction strength among different cell types
#To simplify,we can aggregate the cell-cell communication based on the defined cell groups. 
#Here we categorize the cell populations into three cell types, and then re-merge the list of CellChat object.
group.cellType <- c(rep("Trophoblast", 8), rep("Endothelial_Cell",1), rep("Myeloid_Cell", 2), rep("Stromal_cells", 2), rep("Natural_killer_Cell",1),rep("T_Cell",1)) # grouping cell clusters into fibroblast, DC and TC cells
group.cellType <- factor(group.cellType, levels = c("Trophoblast", "Endothelial_Cell", "Myeloid_Cell","Stromal_cells","Natural_killer_Cell","T_Cell"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
##We then can show the number of interactions or interaction strength between any two cell types in each dataset.

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_major_celltype_cellchat_network_circle_number_for_abortion_CTRL_split_plot.pdf",width=13,height=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

##Simialrly, we can also show the differential number of interactions or interaction strength between any two cell types using circle plot. 
##Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_major_celltype_cellchat_network_circle_changed_number_strenth_for_abortion_CTRL_split_plot.pdf",width=13,height=6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
dev.off()

##Compare the major sources and targets in 2D space allows ready identification of the cell populations with significant changes in sending or receiving signals between different datasets.
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
## Dot size is proportional to the number of inferred links (both outgoing and incoming) associated with each cell group. 
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_out_income_strength_all_signaling_pathways_2D_plot.pdf",width=16,height=8)
patchwork::wrap_plots(plots = gg)
dev.off()


######From the scatter plot, we can see that Inflam.DC and cDC1 emerge as one of the major source and targets in LS compared to NL. Fibroblast populations also become the major sources in LS.
#Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS. 
## Identify signaling changes associated with one cell group
#gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "EVTs_3", signaling.exclude = "MIF")
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "EVTs_3")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "EVTs_2")
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_specific_signaling_changes_of_EVTs_3_EVTs_2.pdf",width=16,height=8)
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()

cell_names<-levels(cellchat@idents$joint)
gg <- list()
for (i in 1:length(cell_names)) {
  gg[[i]] <- netAnalysis_signalingChanges_scatter(cellchat, idents.use =cell_names[i])
}
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_out_income_strength_all_signaling_pathways_2D_plot_for_all_cell.pdf",width=50,height=45)
patchwork::wrap_plots(plots = gg,ncol=4)
dev.off()


# Access all the signaling pathways showing significant communications
pathyway_CPIpair<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/cellchat_pathyway_CPIpair.rds")
CellChatDB <- CellChatDB.human 
CellChatDB_human_interaction_meta<- data.frame(CellChatDB$interaction)[,c("pathway_name","interaction_name","ligand","receptor")]
pathways_all <- unique(CellChatDB_human_interaction_meta$pathway_name)

pathways.show.all <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
length(pathways.show.all)#113
head(CellChatDB_human_interaction_meta)
#all_genes<-cellchat@DB$geneInfo$Symbo
#all_genes<-unique(c(cellchat@var.features$CTRL$features,cellchat@var.features$Abortion$features))
for (i in 1:length(pathways.show.all)) {
  # i=3
  violin_plot<-plotGeneExpression(cellchat, signaling =  pathways.show.all[i],split.by = "datasets", colors.ggplot = T,color.use = ppCor[c(6,10)])
  gene_number<-length(violin_plot[]$patches$plots)
  plot_gg <- list()
  for (j in 1:gene_number) {
   # j=4
    plot_gg[[j]]<- violin_plot[[j]]+ stat_compare_means(method = "wilcox.test",data=violin_plot[[j]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =max(violin_plot[[j]]$data[,1])-0.25)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
  }
  #gene_names<-unique(as.character(unlist(lapply(strsplit(as.character(pathyway_CPIpair[[pathways.show.all[i]]]),"_"), function(x) x))))
  gene_list_merge<-patchwork::wrap_plots(plots = plot_gg,ncol=1)
  ggsave(filename=paste0("/mnt/data/chenwei/gongchen/7.cellchat/Merge/gene_expression/Merge_",pathways.show.all[i], "_gene_expression_violin_plot.pdf"), plot=gene_list_merge, width = 12, height =gene_number*1.5,limitsize = FALSE)
}
###plot delected pathway
pathways_final<-c("WNT","BMP","SEMA7","LEP")

for (i in 1:length(pathways_final)) {
  # i=1
  violin_plot<-plotGeneExpression(cellchat, signaling =  pathways_final[i],split.by = "datasets", colors.ggplot = T,color.use = ppCor[c(6,10)])
  gene_number<-length(violin_plot[]$patches$plots)
  plot_gg <- list()
  for (j in 1:gene_number) {
    # j=4
    plot_gg[[j]]<- violin_plot[[j]]+ stat_compare_means(method = "wilcox.test",data=violin_plot[[j]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =max(violin_plot[[j]]$data[,1])-0.25)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
  }
  #gene_names<-unique(as.character(unlist(lapply(strsplit(as.character(pathyway_CPIpair[[pathways.show.all[i]]]),"_"), function(x) x))))
  gene_list_merge<-patchwork::wrap_plots(plots = plot_gg,ncol=1)
  ggsave(filename=paste0("/mnt/data/chenwei/gongchen/7.cellchat/Merge/gene_expression/final_selected_Merge_",pathways_final[i], "_gene_expression_violin_plot.pdf"), plot=gene_list_merge, width = 12, height =gene_number*1.5,limitsize = FALSE)
}



#Compute and visualize the pathway distance in the learned joint manifold
##可以根据信号网络在共享双维空间中的欧几里德距离来识别差异较大（或更少）的信号网络。
##更大的距离意味着两个数据集之间的通信网络在功能或结构相似性方面存在更大的差异。
#NB：我们只计算两个数据集之间重叠信号通路的距离。此处未考虑仅在一个数据集中标识的信号通路。
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_conserved_and_context_specific_signaling_pathways.pdf",width=8,height=16)
rankSimilarity(cellchat, type = "functional")
dev.off()

##Identify and visualize the conserved and context-specific signaling pathways

##By comparing the information flow/interaction strengh of each signaling pathway, 
##we can identify signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv) increase, 
##by change their information flow at one condition as compared to another condition.
#which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).

#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_barplot_conserved_and_context_specific_signaling_pathways.pdf",width=12,height=16)
gg1 + gg2
dev.off()

rankNet(cellchat, mode = "comparison", stacked = F, do.stat = F)

#Compare outgoing (or incoming) signaling associated with each cell population
##The above analysis summarize the information from the outgoing and incoming signaling together. 
##We can also compare the outgoing (or incoming) signaling pattern between two datasets, allowing to identify signaling pathways/ligand-receptors that exhibit different signaling patterns.
#We can combine all the identified signaling pathways from different datasets and thus compare them side by side
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 26)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height =26)
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_heatmap_outgoing_all_signaling_pathways_abortion_CTRL.pdf",width=13,height=28)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 26, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 26, color.heatmap = "GnBu")
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_heatmap_incoming_all_signaling_pathways_abortion_CTRL.pdf",width=13,height=28)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 26, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 26, color.heatmap = "OrRd")
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_heatmap_all_signaling_pathways_abortion_CTRL.pdf",width=13,height=28)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabities
###We can compare the communication probabilities mediated by ligand-receptor pairs from some cell groups to other cell groups. 
##This can be done by setting comparison in the function netVisual_bubble.
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)

pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_out_EVT2_EVT3_interaction_pair_dot_plot.pdf",width=10,height=28)
netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8:15),comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE)
netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:7,9:15),comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE)
dev.off()
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_in_EVT2_EVT3_interaction_pair_dot_plot.pdf",width=10,height=28)
netVisual_bubble(cellchat, sources.use =c(1:6,8:15), targets.use = 7,comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE)
netVisual_bubble(cellchat, sources.use = c(1:7,9:15), targets.use = 8,comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE)
dev.off()

##Moreover, we can identify the upgulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset 
##compared to the other dataset. This can be done by specifying max.dataset and min.dataset in the function netVisual_bubble. 
##The increased signaling means these signaling have higher communication probability (strength) in one dataset compared to the other dataset.
#> Comparing communications on a merged object
gg1 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8:15),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Abortion", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8:15),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Abortion", angle.x = 45, remove.isolate = T)
gg3 <- netVisual_bubble(cellchat, sources.use =c(1:6,8:15), targets.use = 7,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Abortion", angle.x = 45, remove.isolate = T)
gg4 <- netVisual_bubble(cellchat, sources.use =c(1:6,8:15), targets.use = 7,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Abortion", angle.x = 45, remove.isolate = T)


pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_diff_EVT2_interaction_pair_dot_plot.pdf",width=15,height=22)
gg1 + gg2
gg3 + gg4
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8:15),  comparison = c(1, 2), max.dataset = 2, signaling = c("WNT","BMP","SEMA7","LEP"), title.name = "Increased signaling in Abortion", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8:15),  comparison = c(1, 2), max.dataset = 1,signaling = c("WNT","BMP","SEMA7","LEP"), title.name = "Decreased signaling in Abortion", angle.x = 45, remove.isolate = T)
gg3 <- netVisual_bubble(cellchat, sources.use =c(1:6,8:15), targets.use = 7,  comparison = c(1, 2), max.dataset = 2,signaling = c("WNT","BMP","SEMA7","LEP"), title.name = "Increased signaling in Abortion", angle.x = 45, remove.isolate = T)
gg4 <- netVisual_bubble(cellchat, sources.use =c(1:6,8:15), targets.use = 7,  comparison = c(1, 2), max.dataset = 1,signaling = c("WNT","BMP","SEMA7","LEP"), title.name = "Decreased signaling in Abortion", angle.x = 45, remove.isolate = T)

pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_diff_EVT2_interaction_pair_dot_plot_selected_pathway.pdf",width=15,height=22)
gg1 + gg2
gg3 + gg4
dev.off()

hh1 <- netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:7,9:15),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Abortion", angle.x = 45, remove.isolate = T)
hh2 <- netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:7,9:15),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Abortion", angle.x = 45, remove.isolate = T)
hh3 <- netVisual_bubble(cellchat, sources.use =c(1:7,9:15), targets.use = 8,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Abortion", angle.x = 45, remove.isolate = T)
hh4 <- netVisual_bubble(cellchat, sources.use =c(1:7,9:15), targets.use = 8,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Abortion", angle.x = 45, remove.isolate = T)

pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_diff_EVT3_interaction_pair_dot_plot.pdf",width=15,height=26)
hh1 + hh2
hh3 + hh4
dev.off()

hh1 <- netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:7,9:15),  comparison = c(1, 2), max.dataset = 2,signaling = c("WNT","BMP","SEMA7","LEP"),  title.name = "Increased signaling in Abortion", angle.x = 45, remove.isolate = T)
hh2 <- netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:7,9:15),  comparison = c(1, 2), max.dataset = 1,signaling = c("WNT","BMP","SEMA7","LEP"),  title.name = "Decreased signaling in Abortion", angle.x = 45, remove.isolate = T)
hh3 <- netVisual_bubble(cellchat, sources.use =c(1:7,9:15), targets.use = 8,  comparison = c(1, 2), max.dataset = 2,signaling = c("WNT","BMP","SEMA7","LEP"),  title.name = "Increased signaling in Abortion", angle.x = 45, remove.isolate = T)
hh4 <- netVisual_bubble(cellchat, sources.use =c(1:7,9:15), targets.use = 8,  comparison = c(1, 2), max.dataset = 1,signaling = c("WNT","BMP","SEMA7","LEP"), title.name = "Decreased signaling in Abortion", angle.x = 45, remove.isolate = T)

pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_diff_EVT3_interaction_pair_dot_plot_selected_pathway.pdf",width=15,height=26)
hh1 + hh2
hh3 + hh4
dev.off()

#NB：气泡图中显示的配体受体对可以通过signaling.LS Increased = gg1$data 访问。

##Identify dysfunctional signaling by using differential expression analysis
##The above method for identifying the upgulated and down-regulated signaling is perfomed by comparing the communication probability between two datasets for each L-R pair and each pair of cell groups. 
##Alternative, we can identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential gene expression analysis. 
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Abortion"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in USM
net.up1 <- subsetCommunication(cellchat, net = net, datasets = "Abortion",ligand.logFC = 0.1, receptor.logFC = NULL)
net.up2 <- subsetCommunication(cellchat, net = net, datasets = "Abortion",ligand.logFC = NULL, receptor.logFC = 0.1)
net.up<-distinct(rbind(net.up1,net.up2));dim(net.up)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down1 <- subsetCommunication(cellchat, net = net, datasets = "CTRL",ligand.logFC = -0.1, receptor.logFC = NULL)
net.down2 <- subsetCommunication(cellchat, net = net, datasets = "CTRL",ligand.logFC = NULL, receptor.logFC = -0.1)
net.down<-distinct(rbind(net.down1,net.down2));dim(net.down)

##Since the signaling genes in the net.up and net.down might be complex with multi-subunits, we can do further deconvolution to obtain the individual signaling genes.
gene.up1 <- extractGeneSubsetFromPair(net.up1, cellchat)
gene.up2 <- extractGeneSubsetFromPair(net.up2, cellchat)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
length(unique(c(gene.up1,gene.up2)));length(unique(gene.up1));length(unique(gene.up2))#330 283  288

gene.down1 <- extractGeneSubsetFromPair(net.down1, cellchat)
gene.down2 <- extractGeneSubsetFromPair(net.down2, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

length(unique(c(gene.down1,gene.down2)));length(unique(gene.down1));length(unique(gene.down2))#320 255  250

##We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,sources.use = 7, targets.use = c(1:6,8:15),comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,sources.use = 7,targets.use = c(1:6,8:15),comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg3 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,sources.use =c(1:6,8:15), targets.use = 7,  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg4 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,sources.use =c(1:6,8:15), targets.use = 7,comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_DEG_diff_EVT2_interaction_pair_dot_plot.pdf",width=15,height=22)
gg1 + gg2
gg3 + gg4
dev.off()

hh1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 8, targets.use = c(1:7,9:15),comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
hh2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 8, targets.use = c(1:7,9:15), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
hh3 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,sources.use =c(1:7,9:15), targets.use = 8,  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
hh4 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,sources.use =c(1:7,9:15), targets.use = 8,comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Merge_DEG_diff_EVT3_interaction_pair_dot_plot.pdf",width=15,height=26)
hh1 + hh2
hh3 + hh4
dev.off()

##待细看
##Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]],sources.use = 7, targets.use = c(1:6,8:15),slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]],sources.use = 7, targets.use = c(1:6,8:15),slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


##Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
#在所有可视化图中，边缘颜色与发送者源一致，边缘权重与交互强度成正比。较厚的边缘线表示信号更强
groupSize_USM <- as.numeric(table(cellchat@idents$Abortion))
groupSize_CTRL <- as.numeric(table(cellchat@idents$CTRL))
groupSize_list<-list(CTRL=groupSize_CTRL,USM=groupSize_USM)

###plot delected pathway

pathways_final<-c("WNT","BMP","SEMA7")#,"LEP"
for (i in 1:length(pathways_final)) {
 #i=2
  pathways.show <- pathways_final[[i]] # c("WNT","BMP","SEMA7","LEP"),
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
  
  pdf(paste0("/mnt/data/chenwei/gongchen/7.cellchat/Merge/gene_expression/final_selected_Merge_",pathways_final[i], "_circle_plot.pdf"),width=17,height=8)
  par(mfrow = c(1,2), xpd=TRUE)
  for (j in 1:length(object.list)) {
 # j=2
  netVisual_aggregate(object.list[[j]],vertex.weight=groupSize_list[[j]], signaling = pathways_final[i], layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways_final[i], names(object.list)[j]))
  }
  dev.off()
}
pdf(paste0("/mnt/data/chenwei/gongchen/7.cellchat/Merge/gene_expression/final_selected_Merge_LEP_circle_plot.pdf"),width=8,height=8)
netVisual_aggregate(object.list[[2]],vertex.weight=groupSize_list[[2]], signaling = "LEP", layout = "circle", edge.width.max = 10, signaling.name = paste("LEP", names(object.list)[2]))
dev.off()
#> Do heatmap based on a single object 
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Chord diagram two 
group.cellType <- c(rep("Trophoblast", 8), rep("Endothelial_Cell",1), rep("Myeloid_Cell", 2), rep("Stromal_cells", 2), rep("Natural_killer_Cell",1),rep("T_Cell",1)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

##netVisual_chord_gene is used for visualizing the cell-cell communication mediated by mutiple ligand-receptors or signaling pathways (where each sector in the chord diagram is a ligand, receptor or signaling pathway.)
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 4, targets.use = c(5:8), lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}

# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(8,10),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
}

# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}

