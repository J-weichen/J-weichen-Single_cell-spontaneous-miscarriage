rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(scales)
library(ggsci)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(NMF)
library(VennDiagram)

##可视化
#library(Seurat)
future::plan("multiprocess", workers =40) # do parallel
##plan(strategy = "multicore", workers = 42)
#for ctrl 
cellchat<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/CTRL/ctrl_cellchat.rds")
data.dir <- '/mnt/data/chenwei/gongchen/7.cellchat/CTRL'
setwd(data.dir)

#Part III: Visualization of cell-cell communication network
##For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_cellchat_network_circle_plot.pdf",width=13,height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F,arrow.width = 0.2,arrow.size = 0.2, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F,arrow.width = 0.2,arrow.size = 0.2, title.name = "Interaction weights/strength")
dev.off()

##Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
##Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
mat <- cellchat@net$weight
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_cellchat_network_circle_plot_split.pdf",width=25,height=16)
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize,  margin = 0.1,  edge.curved = 0.3, weight.scale = T, edge.weight.max = max(mat),arrow.width = 0.25,arrow.size = 0.01,title.name = rownames(mat)[i])
}
dev.off()

##Automatically save the plots of the all inferred network for quick exploration
pathyway_CPIpair<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/cellchat_pathyway_CPIpair.rds")

CellChatDB <- CellChatDB.human 
CellChatDB_human_interaction_meta<- data.frame(CellChatDB$interaction)[,c("pathway_name","interaction_name","ligand","receptor")]
pathways_all <- unique(CellChatDB_human_interaction_meta$pathway_name)

pathways_name<-"SEMA3"
#pairs_name<-CellChatDB_human_interaction_meta[which(CellChatDB_human_interaction_meta$pathway_name ==pathways_name ),]$interaction_name 
pathyway_CPIpair[[pathways_name]]
pathyway_gene[[pathways_name]]

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(10,11,12,13) # a numeric vector. 
for (i in 1:length(pathways.show.all)) {
  # i=3
  #length(pathyway_CPIpair[[pathways.show.all[i]]])
  
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i],out.format="pdf", vertex.receiver = vertex.receiver, layout = "hierarchy")
  netVisual(cellchat, signaling = pathways.show.all[i],out.format="pdf",layout = "circle")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(data.dir,"/CTRL_",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = length(pathyway_CPIpair[[pathways.show.all[i]]])*0.3,limitsize = FALSE)
}
#head(cellchat@DB$interaction)
#unique(cellchat@DB$interaction$pathway_name)

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
#不同层次的细胞通信可视化：可以使用netVisual_aggregate可视化信号通路的推断通信网络
#使用netVisual_individual可视化与该信号通路相关的单个L-R对的推断通信网络。

###Hierarchy plot
###This hierarchical plot consist of two components: the left portion shows autocrine and paracrine signaling to certain cell groups of interest (i.e, the defined vertex.receiver), 
###and the right portion shows autocrine and paracrine signaling to the remaining cell groups in the dataset. 
##Here we take input of one signaling pathway as an example. All the signaling pathways showing significant communications can be accessed by 
cellchat@netP$pathways#105
levels(cellchat@meta$re_annotation_TFac)
###for single pathway 
pathways.show <- c("SEMA3") #define one special signaling pathway
# Circle plot
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_SEMA3_path_network_circle_plot.pdf",width=8,height=8)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

# Heatmap
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_SEMA3_path_network_Heatmap_plot.pdf",width=10,height=10)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
# Chord diagram
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_SEMA3_path_network_Chord_plot.pdf",width=8,height=8)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
##CellChat has an independent function netVisual_chord_cell to flexibly visualize the signaling network
# Chord diagram2
group.cellType <- c(rep("Trophoblast", 8), rep("Endothelial_Cell",1), rep("Myeloid_Cell", 2), rep("Stromal_cells", 2), rep("Natural_killer_Cell",1),rep("T_Cell",1)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
dev.off()

#hierarchy
levels(cellchat@meta$re_annotation_TFac)
vertex.receiver = c(10,11,12,13) # a numeric vector. 
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_SEMA3_path_network_hierarchy_plot.pdf",width=12,height=6)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")
dev.off()

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_SEMA3_path_contribution_of_each_ligand-receptor_pair_plot.pdf",width=6,height=6)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

##visualize the cell-cell communication mediated by a single ligand-receptor pair. We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Circle plot
pdf(paste0("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_SEMA3_path_individualpair_",LR.show,"_circle_plot.pdf"),width=8,height=8)
netVisual_individual(cellchat, signaling = pathways.show,pairLR.use = LR.show,  layout = "circle")
dev.off()

# Chord diagram
pdf(paste0("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_SEMA3_path_individualpair_",LR.show,"_Chord_plot.pdf"),width=8,height=8)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

#hierarchy
levels(cellchat@meta$re_annotation_TFac)
vertex.receiver = c(10,11,12,13) # a numeric vector. 
pdf(paste0("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_SEMA3_path_individualpair_",LR.show,"_hierarchy_plot.pdf"),width=12,height=6)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout = "hierarchy")
dev.off()

levels(cellchat@meta$re_annotation_TFac)

##Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
#Bubble plot show all the significant interactions (L-R pairs) from some cell groups to other cell groups.
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
pdf(paste0("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_outcomming_EVT2_EVT3_interaction_pair_dot_plot.pdf"),width=10,height=30)
netVisual_bubble(cellchat, sources.use = c(7,8), targets.use = c(1:6,8:15), remove.isolate = TRUE)
dev.off()
pdf(paste0("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_intcomming_EVT2_EVT3_interaction_pair_dot_plot.pdf"),width=10,height=30)
netVisual_bubble(cellchat, sources.use = c(1:6,8:15), targets.use = c(7,8), remove.isolate = TRUE)
dev.off()

#Chord diagram

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_netP_EVT2_EVT3_interaction_pair_circle_plot.pdf",width=25,height=25)
netVisual_chord_gene(cellchat, sources.use =c(7,8), targets.use =  c(1:6,8:15), slot.name = "netP", legend.pos.x = 10)
netVisual_chord_gene(cellchat, sources.use =c(1:6,8:15), targets.use = c(7,8), slot.name = "netP", legend.pos.x = 10)
dev.off()

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("FGF","CCL"), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:11), pairLR.use = pairLR.use, remove.isolate = TRUE)

#Chord diagram
par(mfrow=c(1,1))
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y =20)
# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)


pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_sender_receiver_mediator_influencer_for_each_pathway_circle_plot.pdf",width=5,height=10)
netAnalysis_signalingRole_network(cellchat, width = 10, height = 5, font.size = 10)
dev.off()
for (i in 1:length(pathways.show.all)) {
  # i=3
  violin_plot<-plotGeneExpression(cellchat, signaling =  pathways.show.all[i])
  ggsave(filename=paste0(data.dir,"/pathway_gene/CTRL_",pathways.show.all[i], "_gene_expression_violin_plot.pdf"), plot=violin_plot, width = 10, height = length(pathyway_CPIpair[[pathways.show.all[i]]])*0.3,limitsize = FALSE)
  }

##Plot the signaling gene expression distribution using violin/dot plot
#We can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.
plotGeneExpression(cellchat, signaling = "CXCL")

###Part IV: Systems analysis of cell-cell communication network
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show <- c("CXCL") #define one special signaling pathway
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 12, height = 3.5, font.size = 10)

##Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_out_income_all_signaling_pathways_2D_plot.pdf",width=8,height=8)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()
##Identify signals contributing most to outgoing or incoming signaling of certain cell groups
## We can also answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 10,height =25)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 10,height =25)
ht1 + ht2

pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_heatmap_out_income_all_signaling_pathways_plot.pdf",width=21,height=25)
ht1 + ht2
dev.off()

grid.newpage()

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))
ht

##Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
##除了探索单个通路的详细通信外，一个重要问题是多个细胞组和信号通路如何协调功能。
#CellChat 采用模式识别方法识别全局通信模式。
#As the number of patterns increases, there might be redundant patterns, making it difficult to interpret the communication patterns. We chose five patterns as default. Generally, it is biologically meaningful with the number of patterns greater than 2. 
##Cophenetic and Silhouette. Both metrics measure the stability for a particular number of patterns based on a hierarchical clustering of the consensus matrix. For a range of the number of patterns, a suitable number of patterns is the one at which Cophenetic and Silhouette values begin to drop suddenly.

##Identify and visualize outgoing communication pattern of secreting cells
#we constructed a dot plot in which the dot size is proportion to the contribution score to show association between cell group and their enriched signaling pathways. 
#USERS can also decrease the parameter cutoff to show more enriched signaling pathways associated each cell group.

##run selectK to infer the number of patterns.
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_lineplot_for_pattern_number.pdf",width=5,height=5)
selectK(cellchat, pattern = "outgoing")##tine consume
selectK(cellchat, pattern = "incoming")
dev.off()
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 5.
nPatterns = 3
grid.newpage()
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_3_outgoing_pattern_identification.pdf",width=15,height=25)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 15,height =25)
dev.off()
# river plot
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_alluvial_plot_3_outgoing_pattern.pdf",width=21,height=20)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
# dot plot
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_dot_plot_3_outgoing_pattern.pdf",width=20,height=12)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

grid.newpage()
#Cophenetic values begin to drop when the number of incoming patterns is 3.
nPatterns = 3
grid.newpage()
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_4_incoming_pattern_identification.pdf",width=17,height=25)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, width = 17,height =25)
dev.off()
# river plot
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_alluvial_plot_4_incoming_pattern.pdf",width=21,height=20)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
# dot plot
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_dot_plot_4_incoming_pattern.pdf",width=20,height=12)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

# Visualization  functional similarity  in 2D-space
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_Dim_plot_4_functional_similarity.pdf",width=20,height=20)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()

# Visualization  structure similarity in 2D-space
pdf("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/CTRL_Dim_plot_4_structural_similarity.pdf",width=20,height=20)
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

#sessionInfo()
