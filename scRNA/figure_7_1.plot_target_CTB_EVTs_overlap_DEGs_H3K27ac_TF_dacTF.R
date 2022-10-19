rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

options(stringsAsFactors = FALSE)
cell_name="ALLcell_all_genes"

TF_gene_ALL <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".Cyto_Network.txt"), header=T)
head(TF_gene_ALL)
TF_gene_weight <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s2_",cell_name,".adj.tsv"), header=T)  
head(TF_gene_weight);dim(TF_gene_weight)

##read single cell DEGs
#get overlapped gene between developped genes and abortion DEGs 
merge_data0 <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge_new.txt", header=T)
head(merge_data0)
cell_subtype<-c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3")
USM_CTBs_EVTs_DEGs_data1<-merge_data0[which(merge_data0$cluster %in% cell_subtype),];dim(USM_CTBs_EVTs_DEGs_data1)#1976    6

##filter TF_regulon with only DEGs in CTBs and EVTs
USM_CTBs_EVTs_DEGs<-unique(as.character(USM_CTBs_EVTs_DEGs_data1$gene))
TF_gene_filter<-TF_gene_ALL[which(TF_gene_ALL$target_final_gene %in% USM_CTBs_EVTs_DEGs),];dim(TF_gene_filter)#16510     2
head(TF_gene_filter)
head(USM_CTBs_EVTs_DEGs_data1)

###plot network for all target TFs and DEGs
##for CTBs
TFs_name<-c("TEF","TFAP2A","TFAP2C")
USM_CTBs_DEGs_data1<-USM_CTBs_EVTs_DEGs_data1[which(USM_CTBs_EVTs_DEGs_data1$cluster %in% c("CTBs_1","CTBs_2")),]
USM_CTBs_DEGs<-unique(as.character(USM_CTBs_DEGs_data1$gene))
TF_gene_filter_CTBs<-TF_gene_ALL[which(TF_gene_ALL$target_final_gene %in% USM_CTBs_DEGs),];dim(TF_gene_filter_CTBs)#5119    2

dim(USM_CTBs_DEGs_data1)# 598   6
TFs_name<-c("TEF","TFAP2A","TFAP2C")
#TFs_name<-c("TEF","TFAP2A")
target_TF_filter_CTBs<-TF_gene_filter_CTBs[which(TF_gene_filter_CTBs$TFs %in% TFs_name),]
targets_trend<-USM_CTBs_DEGs_data1[which(USM_CTBs_DEGs_data1$gene %in% unique(c(unique(target_TF_filter_CTBs$target_final_gene),TFs_name))),]

####select target_final regulon and corresponding weight
TF_gene_weight1<-TF_gene_weight[which(TF_gene_weight$TF %in% TFs_name),]
head(TF_gene_weight1);dim(TF_gene_weight1)
TF_gene_weight1$pair<-paste0(TF_gene_weight1$TF,"_",TF_gene_weight1$target)
nrow(TF_gene_weight1)#216822
head(TF_gene_weight1)

target_TF_gene_final<-target_TF_filter_CTBs
target_TF_gene_final$pair<-paste0(target_TF_gene_final$TFs,"_",target_TF_gene_final$target_final_gene)
head(target_TF_gene_final)
target_TF_gene_final_1<-TF_gene_weight1[which(TF_gene_weight1$pair %in% target_TF_gene_final$pair),]
target_TF_gene_final_1
nrow(target_TF_gene_final_1);nrow(target_TF_gene_final_1)#20
range(target_TF_gene_final_1$importance) ##0.4312386 61.1091759
table(target_TF_gene_final_1$TF);length(table(target_TF_gene_final_1$target))
#TEF TFAP2A TFAP2C 
#1      6     13 
##get gene pair link number
TF_gene_final_stat<-data.frame(table(c(target_TF_gene_final_1$TF,target_TF_gene_final_1$target)))
colnames(TF_gene_final_stat)<-c("gene","count")
TF_gene_final_stat<-TF_gene_final_stat[order(TF_gene_final_stat$count,decreasing = T),]
TF_gene_final_stat$type<-ifelse(TF_gene_final_stat$gene %in% unique(target_TF_gene_final_1$TF),"TFs","targets")
head(TF_gene_final_stat)
table(TF_gene_final_stat$type)
TF_gene_final_stat$node<-as.character(TF_gene_final_stat$gene)

##get Trend and cell number of DEGs
head(targets_trend)
targets_trend2<-targets_trend[which(!(targets_trend$gene %in% TFs_name)),]
targets_split <-dcast(targets_trend2, gene ~ trend)
targets_split$Cell_sum<-targets_split$Down+targets_split$Up
targets_split$Up_percent<-targets_split$Up/targets_split$Cell_sum*100
targets_split$DEG_type<- ifelse(targets_split$Down> 0  & targets_split$Up > 0, "Contrary",ifelse(targets_split$Down == 0, "Up_only","Down_only"))
table(targets_split$DEG_type)
#Down_only   Up_only 
# 5        11 

head(targets_split)

##merge for node information
node_info_merge<-merge(TF_gene_final_stat,targets_split,by="gene",all=T)

algorithm_plot<-"nicely";regulon_DEG<-"CTBs_TEF_TFAP2A_TFAP2C"
node_info_merge[which(is.na(node_info_merge$DEG_type)),]$DEG_type<-"TFs"
node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$Cell_sum<-node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$count
node_info_merge[which(node_info_merge$DEG_type=="TFs"),]
node_info_merge$color<-node_info_merge$DEG_type
node_info_merge$color<-ifelse(node_info_merge$color =="Up_only","indianred3",ifelse(node_info_merge$color =="Down_only","steelblue",
                                                                                    ifelse(node_info_merge$color =="Contrary","lightgreen","purple")))

nodes<-node_info_merge
edges<-target_TF_gene_final_1[,c("TF","target","importance")]

net<-graph_from_data_frame(d=edges,vertices = nodes,directed = T)
node_color<-node_info_merge$color

set.seed(19921010)
layout <- create_layout(net, layout ="nicely")
head(layout)
CTBs_all_network<-ggraph(layout)+
  geom_edge_fan(aes(edge_width=importance),color="lightblue") + #,show.legend=FALSE
  geom_node_point(aes(size=Cell_sum,color=as.factor(DEG_type)))+
  geom_node_text(aes(label=node),size=3)+
  #geom_node_text(aes(filter= DEG_type =="TFs",label=node),size=3)+
  #geom_node_text(aes(filter= (DEG_type =="TFs" | Cell_sum>1),label=node),size=3)+
  scale_color_manual(limits = as.factor(layout$DEG_type), values =node_color)+ 
  #  scale_edge_color_continuous(low = "cyan",high = "red")+
  scale_size_area(max_size = 20)+guides(size="none")+
  scale_edge_width(range = c(0.5, 3))+ 
  guides(color=guide_legend(order=3))+  theme_graph()+ 
  ggtitle(paste(regulon_DEG,":",algorithm_plot)) + 
  theme(plot.title = element_text(hjust = 0.5))
CTBs_all_network
#pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/CTBs_TEF_TFAP2A_TFAP2C_all_network_plot.pdf",width = 6,height=6)
#ggsave(CTBs_all_network,"/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/CTBs_TEF_TFAP2A_TFAP2C_all_network_plot.pdf", width=6, height =6)

##for EVTs
TFs_name<-c("ARNT2","SMAD3","HSF1","ELF3")
USM_EVTs_DEGs_data1<-USM_CTBs_EVTs_DEGs_data1[which(USM_CTBs_EVTs_DEGs_data1$cluster %in% c("EVTs_1","EVTs_2","EVTs_3")),]
USM_EVTs_DEGs<-unique(as.character(USM_EVTs_DEGs_data1$gene))
TF_gene_filter_EVTs<-TF_gene_ALL[which(TF_gene_ALL$target_final_gene %in% USM_EVTs_DEGs),];dim(TF_gene_filter)#14377     2
target_TF_filter_EVTs<-TF_gene_filter_EVTs[which(TF_gene_filter_EVTs$TFs %in% TFs_name),]
targets_trend<-USM_EVTs_DEGs_data1[which(USM_EVTs_DEGs_data1$gene %in% unique(c(unique(target_TF_filter_EVTs$target_final_gene),TFs_name))),]

####select target_final regulon and corresponding weight
TF_gene_weight1<-TF_gene_weight[which(TF_gene_weight$TF %in% TFs_name),]
head(TF_gene_weight1);dim(TF_gene_weight1)
TF_gene_weight1$pair<-paste0(TF_gene_weight1$TF,"_",TF_gene_weight1$target)
nrow(TF_gene_weight1)#216822
head(TF_gene_weight1)

target_TF_gene_final<-target_TF_filter_EVTs
target_TF_gene_final$pair<-paste0(target_TF_gene_final$TFs,"_",target_TF_gene_final$target_final_gene)
head(target_TF_gene_final)
target_TF_gene_final_1<-TF_gene_weight1[which(TF_gene_weight1$pair %in% target_TF_gene_final$pair),]
target_TF_gene_final_1
nrow(target_TF_gene_final_1);nrow(target_TF_gene_final_1)#22 22
range(target_TF_gene_final_1$importance) ##3.008961 32.770553
table(target_TF_gene_final_1$TF);length(table(target_TF_gene_final_1$target))

##get gene pair link number
TF_gene_final_stat<-data.frame(table(c(target_TF_gene_final_1$TF,target_TF_gene_final_1$target)))
colnames(TF_gene_final_stat)<-c("gene","count")
TF_gene_final_stat<-TF_gene_final_stat[order(TF_gene_final_stat$count,decreasing = T),]
TF_gene_final_stat$type<-ifelse(TF_gene_final_stat$gene %in% unique(target_TF_gene_final_1$TF),"TFs","targets")
head(TF_gene_final_stat)
table(TF_gene_final_stat$type)
TF_gene_final_stat$node<-as.character(TF_gene_final_stat$gene)

##get Trend and cell number of DEGs
head(targets_trend)
targets_trend2<-targets_trend[which(!(targets_trend$gene %in% TFs_name)),]
targets_split <-dcast(targets_trend2, gene ~ trend)
targets_split$Cell_sum<-targets_split$Down+targets_split$Up
targets_split$Up_percent<-targets_split$Up/targets_split$Cell_sum*100
targets_split$DEG_type<- ifelse(targets_split$Down> 0  & targets_split$Up > 0, "Contrary",ifelse(targets_split$Down == 0, "Up_only","Down_only"))
table(targets_split$DEG_type)
#Contrary Down_only   Up_only 
#  4       182        96 

head(targets_split)

##merge for node information
node_info_merge<-merge(TF_gene_final_stat,targets_split,by="gene",all=T)
algorithm_plot<-"nicely";regulon_DEG<-"EVTs_HSF1_ELF3_SMAD3_ARNT2"
node_info_merge[which(is.na(node_info_merge$DEG_type)),]$DEG_type<-"TFs"
node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$Cell_sum<-node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$count
node_info_merge[which(node_info_merge$DEG_type=="TFs"),]
node_info_merge$color<-node_info_merge$DEG_type
node_info_merge$color<-ifelse(node_info_merge$color =="Up_only","indianred3",ifelse(node_info_merge$color =="Down_only","steelblue",
                                                                                    ifelse(node_info_merge$color =="Contrary","lightgreen","purple")))

nodes<-node_info_merge
edges<-target_TF_gene_final_1[,c("TF","target","importance")]

net<-graph_from_data_frame(d=edges,vertices = nodes,directed = T)
node_color<-node_info_merge$color

set.seed(19921010)
layout <- create_layout(net, layout ="nicely")
head(layout)

EVTs_all_network<-ggraph(layout)+
  geom_edge_fan(aes(edge_width=importance),color="lightblue") + #,show.legend=FALSE
  geom_node_point(aes(size=Cell_sum,color=as.factor(DEG_type)))+
  #geom_node_text(aes(filter= DEG_type =="TFs",label=node),size=3)+
  geom_node_text(aes(filter= (DEG_type =="TFs" | count>1),label=node),size=3)+
  scale_color_manual(limits = as.factor(layout$DEG_type), values =node_color)+ 
  #  scale_edge_color_continuous(low = "cyan",high = "red")+
  scale_size_area(max_size = 20)+guides(size="none")+
  scale_edge_width(range = c(0.5, 3))+ 
  guides(color=guide_legend(order=3))+  theme_graph()+ 
  ggtitle(paste(regulon_DEG,":",algorithm_plot)) + 
  theme(plot.title = element_text(hjust = 0.5))
EVTs_all_network
#pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/EVTs_HSF1_ELF3_SMAD3_ARNT2_all_network_plot.pdf",width = 6,height=6)
#ggsave(EVTs_all_network,"/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/EVTs_HSF1_ELF3_SMAD3_ARNT2_all_network_plot.pdf", width=6, height =6)


##obtain DEGs in specific TFs
##for EVTs
##repress
TFs_name<-c("HSF1","ELF3")
#TFs_name<-c("SMAD3");TFs_name<-c("ARNT2")
#TFs_name<-c("ARNT2","SMAD3","HSF1","ELF3")
USM_EVTs_DEGs_data1<-USM_CTBs_EVTs_DEGs_data1[which(USM_CTBs_EVTs_DEGs_data1$cluster %in% c("EVTs_1","EVTs_2","EVTs_3")),]
USM_EVTs_DEGs<-unique(as.character(USM_EVTs_DEGs_data1$gene))
TF_gene_filter_EVTs<-TF_gene_ALL[which(TF_gene_ALL$target_final_gene %in% USM_EVTs_DEGs),];dim(TF_gene_filter)#16510     2
target_TF_filter_EVTs<-TF_gene_filter_EVTs[which(TF_gene_filter_EVTs$TFs %in% TFs_name),]
targets_trend<-USM_EVTs_DEGs_data1[which(USM_EVTs_DEGs_data1$gene %in% unique(c(unique(target_TF_filter_EVTs$target_final_gene),TFs_name))),]

####select target_final regulon and corresponding weight
TF_gene_weight1<-TF_gene_weight[which(TF_gene_weight$TF %in% TFs_name),]
head(TF_gene_weight1);dim(TF_gene_weight1)
TF_gene_weight1$pair<-paste0(TF_gene_weight1$TF,"_",TF_gene_weight1$target)
nrow(TF_gene_weight1)#16916
head(TF_gene_weight1)

target_TF_gene_final<-target_TF_filter_EVTs
target_TF_gene_final$pair<-paste0(target_TF_gene_final$TFs,"_",target_TF_gene_final$target_final_gene)
head(target_TF_gene_final)
target_TF_gene_final_1<-TF_gene_weight1[which(TF_gene_weight1$pair %in% target_TF_gene_final$pair),]
target_TF_gene_final_1
nrow(target_TF_gene_final_1);nrow(target_TF_gene_final_1)#22 22
range(target_TF_gene_final_1$importance) ##3.008961 32.770553
table(target_TF_gene_final_1$TF);length(table(target_TF_gene_final_1$target))

##get gene pair link number
TF_gene_final_stat<-data.frame(table(c(target_TF_gene_final_1$TF,target_TF_gene_final_1$target)))
colnames(TF_gene_final_stat)<-c("gene","count")
TF_gene_final_stat<-TF_gene_final_stat[order(TF_gene_final_stat$count,decreasing = T),]
TF_gene_final_stat$type<-ifelse(TF_gene_final_stat$gene %in% unique(target_TF_gene_final_1$TFs),"TFs","targets")
head(TF_gene_final_stat)
table(TF_gene_final_stat$type)
TF_gene_final_stat$node<-as.character(TF_gene_final_stat$gene)

##get Trend and cell number of DEGs
head(targets_trend)
targets_trend2<-targets_trend[which(!(targets_trend$gene %in% TFs_name)),]
targets_split <-dcast(targets_trend2, gene ~ trend)
targets_split$Cell_sum<-targets_split$Down+targets_split$Up
targets_split$Up_percent<-targets_split$Up/targets_split$Cell_sum*100
targets_split$DEG_type<- ifelse(targets_split$Down> 0  & targets_split$Up > 0, "Contrary",ifelse(targets_split$Down == 0, "Up_only","Down_only"))
table(targets_split$DEG_type)
#Contrary Down_only   Up_only 
#   4       180        76 

head(targets_split)

##merge for node information
algorithm_plot<-"nicely";regulon_DEG<-"EVTs_repress_HSF1_ELF3"
node_info_merge<-merge(TF_gene_final_stat,targets_split,by="gene",all=T)
node_info_merge[which(is.na(node_info_merge$DEG_type)),]$DEG_type<-"TFs"
node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$Cell_sum<-node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$count
node_info_merge[which(node_info_merge$DEG_type=="TFs"),]
node_info_merge$color<-node_info_merge$DEG_type
node_info_merge$color<-ifelse(node_info_merge$color =="Up_only","indianred3",ifelse(node_info_merge$color =="Down_only","steelblue",ifelse(node_info_merge$color =="Contrary","lightgreen","purple")))

nodes<-node_info_merge
edges<-target_TF_gene_final_1[,c("TF","target","importance")]

net<-graph_from_data_frame(d=edges,vertices = nodes,directed = T)
node_color<-node_info_merge$color

set.seed(921010)
layout <- create_layout(net, layout ="nicely")
head(layout)

EVTs_repress_network<-ggraph(layout)+
  geom_edge_fan(aes(edge_width=importance),color="lightblue") + #,show.legend=FALSE
  geom_node_point(aes(size=Cell_sum,color=as.factor(DEG_type)))+
  geom_node_text(aes(filter= (DEG_type =="TFs"| count>1),label=node),size=3)+
  scale_color_manual(limits = as.factor(layout$DEG_type), values =node_color)+ 
  #  scale_edge_color_continuous(low = "cyan",high = "red")+
  scale_size_area(max_size = 20)+guides(size="none")+
  scale_edge_width(range = c(0.5, 3))+ 
  guides(color=guide_legend(order=3))+  theme_graph()+ 
  ggtitle(paste(regulon_DEG,":",algorithm_plot)) + 
  theme(plot.title = element_text(hjust = 0.5))

EVTs_repress_network
#ggsave(EVTs_repress_network,"/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/EVTs_repress_HSF1_ELF3_network_plot.pdf", width=6, height =6)

##Active
TFs_name<-c("ARNT2","SMAD3")

#TFs_name<-c("HSF1","ELF3")
#TFs_name<-c("SMAD3");TFs_name<-c("ARNT2")
#TFs_name<-c("ARNT2","SMAD3","HSF1","ELF3")
USM_EVTs_DEGs_data1<-USM_CTBs_EVTs_DEGs_data1[which(USM_CTBs_EVTs_DEGs_data1$cluster %in% c("EVTs_1","EVTs_2","EVTs_3")),]
USM_EVTs_DEGs<-unique(as.character(USM_EVTs_DEGs_data1$gene))
TF_gene_filter_EVTs<-TF_gene_ALL[which(TF_gene_ALL$target_final_gene %in% USM_EVTs_DEGs),];dim(TF_gene_filter)#14377     2
target_TF_filter_EVTs<-TF_gene_filter_EVTs[which(TF_gene_filter_EVTs$TFs %in% TFs_name),]
targets_trend<-USM_EVTs_DEGs_data1[which(USM_EVTs_DEGs_data1$gene %in% unique(c(unique(target_TF_filter_EVTs$target_final_gene),TFs_name))),]

####select target_final regulon and corresponding weight
TF_gene_weight1<-TF_gene_weight[which(TF_gene_weight$TF %in% TFs_name),]
head(TF_gene_weight1);dim(TF_gene_weight1)
TF_gene_weight1$pair<-paste0(TF_gene_weight1$TF,"_",TF_gene_weight1$target)
nrow(TF_gene_weight1)#216822
head(TF_gene_weight1)

target_TF_gene_final<-target_TF_filter_EVTs
target_TF_gene_final$pair<-paste0(target_TF_gene_final$TFs,"_",target_TF_gene_final$target_final_gene)
head(target_TF_gene_final)
target_TF_gene_final_1<-TF_gene_weight1[which(TF_gene_weight1$pair %in% target_TF_gene_final$pair),]
target_TF_gene_final_1
nrow(target_TF_gene_final_1);nrow(target_TF_gene_final_1)#22 22
range(target_TF_gene_final_1$importance) ##3.008961 32.770553
table(target_TF_gene_final_1$TF);length(table(target_TF_gene_final_1$target))

##get gene pair link number
TF_gene_final_stat<-data.frame(table(c(target_TF_gene_final_1$TF,target_TF_gene_final_1$target)))
colnames(TF_gene_final_stat)<-c("gene","count")
TF_gene_final_stat<-TF_gene_final_stat[order(TF_gene_final_stat$count,decreasing = T),]
TF_gene_final_stat$type<-ifelse(TF_gene_final_stat$gene %in% unique(target_TF_gene_final_1$TFs),"TFs","targets")
head(TF_gene_final_stat)
table(TF_gene_final_stat$type)
TF_gene_final_stat$node<-as.character(TF_gene_final_stat$gene)

##get Trend and cell number of DEGs
head(targets_trend)
targets_trend2<-targets_trend[which(!(targets_trend$gene %in% TFs_name)),]
targets_split <-dcast(targets_trend2, gene ~ trend)
targets_split$Cell_sum<-targets_split$Down+targets_split$Up
targets_split$Up_percent<-targets_split$Up/targets_split$Cell_sum*100
targets_split$DEG_type<- ifelse(targets_split$Down> 0  & targets_split$Up > 0, "Contrary",ifelse(targets_split$Down == 0, "Up_only","Down_only"))
table(targets_split$DEG_type)
#Down_only   Up_only 
#   3        22 

head(targets_split)

##merge for node information
node_info_merge<-merge(TF_gene_final_stat,targets_split,by="gene",all=T)
node_info_merge[which(is.na(node_info_merge$DEG_type)),]$DEG_type<-"TFs"
node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$Cell_sum<-node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$count

node_info_merge$color<-node_info_merge$DEG_type
node_info_merge$color<-ifelse(node_info_merge$color =="Up_only","indianred3",ifelse(node_info_merge$color =="Down_only","steelblue","orange") )

nodes<-node_info_merge
edges<-target_TF_gene_final_1[,c("TF","target","importance")]

net<-graph_from_data_frame(d=edges,vertices = nodes,directed = T)
node_color<-node_info_merge$color

algorithm_plot<-"nicely";regulon_DEG<-"EVTs_active_SMAD3_ARNT2"

set.seed(921010)
layout <- create_layout(net, layout ="nicely")
head(layout)

EVTs_active_network<-ggraph(layout)+
  #  geom_edge_link(aes(edge_width=importance,edge_color=importance))+#,edge_color=color
  geom_edge_fan(aes(edge_width=importance),color="lightblue") + #,show.legend=FALSE
  #geom_node_point(aes(size=count,color=as.factor(DEG_type)))+
  geom_node_point(aes(size=Cell_sum,color=as.factor(DEG_type)))+
  geom_node_text(aes(label= node),size=3)+
  scale_color_manual(limits = as.factor(layout$DEG_type), values =node_color)+ 
  #  scale_edge_color_continuous(low = "cyan",high = "red")+
  scale_size_area(max_size = 20)+guides(size="none")+
  scale_edge_width(range = c(0.5, 3))+ 
  guides(color=guide_legend(order=3))+  theme_graph()+ 
  ggtitle(paste(regulon_DEG,":",algorithm_plot)) + 
  theme(plot.title = element_text(hjust = 0.5))
#labs(tag = "Figure4b") + theme_bw()#+ xlim(-0.8,0.8)+ylim(-0.8,0.8)
EVTs_active_network
#ggsave(EVTs_active_network,"/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/EVTs_SMAD3_ARNT2_active_network_plot.pdf", width=6, height =6)

##for single TFs

##Active
TFs_name<-c("SMAD3")
#TFs_name<-c("ARNT2","SMAD3","HSF1","ELF3")
USM_EVTs_DEGs_data1<-USM_CTBs_EVTs_DEGs_data1[which(USM_CTBs_EVTs_DEGs_data1$cluster %in% c("EVTs_1","EVTs_2","EVTs_3")),]
USM_EVTs_DEGs<-unique(as.character(USM_EVTs_DEGs_data1$gene))
TF_gene_filter_EVTs<-TF_gene_ALL[which(TF_gene_ALL$target_final_gene %in% USM_EVTs_DEGs),];dim(TF_gene_filter)#14377     2
target_TF_filter_EVTs<-TF_gene_filter_EVTs[which(TF_gene_filter_EVTs$TFs %in% TFs_name),]
targets_trend<-USM_EVTs_DEGs_data1[which(USM_EVTs_DEGs_data1$gene %in% unique(c(unique(target_TF_filter_EVTs$target_final_gene),TFs_name))),]

####select target_final regulon and corresponding weight
TF_gene_weight1<-TF_gene_weight[which(TF_gene_weight$TF %in% TFs_name),]
head(TF_gene_weight1);dim(TF_gene_weight1)
TF_gene_weight1$pair<-paste0(TF_gene_weight1$TF,"_",TF_gene_weight1$target)
nrow(TF_gene_weight1)#216822
head(TF_gene_weight1)

target_TF_gene_final<-target_TF_filter_EVTs
target_TF_gene_final$pair<-paste0(target_TF_gene_final$TFs,"_",target_TF_gene_final$target_final_gene)
head(target_TF_gene_final)
target_TF_gene_final_1<-TF_gene_weight1[which(TF_gene_weight1$pair %in% target_TF_gene_final$pair),]
target_TF_gene_final_1
nrow(target_TF_gene_final_1);nrow(target_TF_gene_final_1)#22 22
range(target_TF_gene_final_1$importance) ##3.008961 32.770553
table(target_TF_gene_final_1$TF);length(table(target_TF_gene_final_1$target))

##get gene pair link number
TF_gene_final_stat<-data.frame(table(c(target_TF_gene_final_1$TF,target_TF_gene_final_1$target)))
colnames(TF_gene_final_stat)<-c("gene","count")
TF_gene_final_stat<-TF_gene_final_stat[order(TF_gene_final_stat$count,decreasing = T),]
TF_gene_final_stat$type<-ifelse(TF_gene_final_stat$gene %in% unique(target_TF_gene_final_1$TFs),"TFs","targets")
head(TF_gene_final_stat)
table(TF_gene_final_stat$type)
TF_gene_final_stat$node<-as.character(TF_gene_final_stat$gene)

##get Trend and cell number of DEGs
head(targets_trend)
targets_trend2<-targets_trend[which(!(targets_trend$gene %in% TFs_name)),]
targets_split <-dcast(targets_trend2, gene ~ trend)
targets_split$Cell_sum<-targets_split$Down+targets_split$Up
targets_split$Up_percent<-targets_split$Up/targets_split$Cell_sum*100
targets_split$DEG_type<- ifelse(targets_split$Down> 0  & targets_split$Up > 0, "Contrary",ifelse(targets_split$Down == 0, "Up_only","Down_only"))
table(targets_split$DEG_type)
#Down_only   Up_only 
#   3        22 

head(targets_split)

##merge for node information
node_info_merge<-merge(TF_gene_final_stat,targets_split,by="gene",all=T)
node_info_merge[which(is.na(node_info_merge$DEG_type)),]$DEG_type<-"TFs"
node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$Cell_sum<-node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$count

node_info_merge$color<-node_info_merge$DEG_type
node_info_merge$color<-ifelse(node_info_merge$color =="Up_only","indianred3",ifelse(node_info_merge$color =="Down_only","steelblue","orange") )

nodes<-node_info_merge
edges<-target_TF_gene_final_1[,c("TF","target","importance")]

net<-graph_from_data_frame(d=edges,vertices = nodes,directed = T)
node_color<-node_info_merge$color
algorithm_plot<-"tree_circular";regulon_DEG<-"EVTs_SMAD3"

set.seed(19921010)
layout <- create_layout(net, layout = 'tree', circular = TRUE)
head(layout)

Single_EVTs_active_network<-ggraph(layout)+
  #  geom_edge_link(aes(edge_width=importance,edge_color=importance))+#,edge_color=color
  geom_edge_fan(aes(edge_width=importance),color="lightblue") + #,show.legend=FALSE
  #geom_node_point(aes(size=count,color=as.factor(DEG_type)))+
  geom_node_point(aes(size=Cell_sum,color=as.factor(DEG_type)))+
  geom_node_text(aes(label= node),size=3)+
  scale_color_manual(limits = as.factor(layout$DEG_type), values =node_color)+ 
  #  scale_edge_color_continuous(low = "cyan",high = "red")+
  scale_size_area(max_size = 20)+guides(size="none")+
  scale_edge_width(range = c(0.5, 3))+ 
  guides(color=guide_legend(order=3))+  theme_graph()+ 
  ggtitle(paste(regulon_DEG,":",algorithm_plot)) + 
  theme(plot.title = element_text(hjust = 0.5))
#labs(tag = "Figure4b") + theme_bw()#+ xlim(-0.8,0.8)+ylim(-0.8,0.8)
Single_EVTs_active_network
#ggsave(Single_EVTs_active_network,"/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/single_EVTs_SMAD3_active_network_plot.pdf", width=6, height =6)

##以下备用
#graph_layouts <- c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl', 'lgl')
#ggraph(net)+ geom_edge_link()+ geom_node_point()
#ggraph(net, 'tree') +  geom_edge_diagonal()
#ggraph(net, 'treemap', weight = count) + geom_node_tile(aes(fill = count), size = 0.25)
layout <- create_layout(net, layout = 'centrality',cent=count)
layout <- create_layout(net, layout = 'igraph', algorithm = 'nicely')

#ggraph(net, layout = 'kk') + 
layout <- create_layout(net, layout = 'circle')
layout <- create_layout(net, layout = 'star')
layout <- create_layout(net, layout = 'stress')
layout <- create_layout(net, layout = 'eigen')
layout <- create_layout(net, layout = 'kk')
layout <- create_layout(net, layout = 'centrality',cent=count)
layout <- create_layout(net, layout = 'linear')
layout <- create_layout(net, layout = 'linear', circular = TRUE)
layout <- create_layout(net, layout = 'igraph', algorithm = 'nicely')
layout <- create_layout(net, layout = 'centrality',cent=count,circular = TRUE)
layout <- create_layout(net, layout = 'dendrogram')
layout <- create_layout(net, layout = 'dendrogram', circular = TRUE)
layout <- create_layout(net, layout = 'tree', circular = TRUE)
layout <- create_layout(net, layout ="nicely")


##展示形式：https://cran.r-project.org/web/packages/ggraph/vignettes/Layouts.html

set.seed(1)
ggraph(net, 'circlepack', weight = count) + 
  geom_edge_link() +  geom_node_point(aes(colour = type)) +coord_fixed()
ggraph(net, 'circlepack', weight = count) + 
  geom_node_circle(aes(fill = count), size = 0.25, n = 50) + coord_fixed()


lapply(c('stress', 'fr', 'lgl', 'graphopt'), function(layout) {
  
  ggraph(net, layout = layout) + 
    geom_edge_link(aes(colour = factor(type)), show.legend = FALSE) +
    geom_node_point() + 
    labs(caption = paste0('Layout: ', layout))
})  

#ggplot2有一个分片的功能，可以根据不同的属性取值绘制子图。ggraph也“继承”了这样一个功能，使用facet_nodes()，我们对4个小型社区进行分片绘制，并添加节点标签（图3）：

graph_gt%>% filter(group %in% 5:8)%>% ggraph(layout = 'kk') + geom_edge_fan(color="lightblue") + 
  geom_node_point(aes(size = deg),fill="orange",shape=21)+
  geom_node_text(aes(label=name),size=3,repel=TRUE)+ 
  guides(size=F)+   facet_nodes(~group) + 
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white') 