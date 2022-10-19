rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
dyn.load('/mnt/data/chenwei/software/share_road/libudunits2/libudunits2.so.0')
dyn.load('/mnt/data/chenwei/software/share_road/libproj/libproj.so.0')
dyn.load('/mnt/data/chenwei/software/share_road/libgdal.so.20')
dyn.load('/mnt/data/chenwei/software/share_road/libgeos-3.4.2.so')
dyn.load('/mnt/data/chenwei/software/share_road/libgeos_c.so.1')

library(scales)
library(ggsci)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dendextend)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)

pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

###expression file reading
data_plot0_USM<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/USM_All_Trophoblast_CTBs2EVTs_k27ac_RNA_overlapped_gene_expression.txt", sep = '\t', row.names = 1, header = T)
data_plot0_CTRL<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/CTRL_All_Trophoblast_CTBs2EVTs_k27ac_RNA_overlapped_gene_expression.txt", sep = '\t', row.names = 1, header = T)

plot_matrix_selected<-data_plot0_USM[!colnames(data_plot0_USM) %in% c("Treat","re_annotation_TFac","pseudotime")]
plot_matrix_selected2<-as.data.frame(t(plot_matrix_selected))
colData_USM<-data_plot0_USM[,c("Treat","re_annotation_TFac","pseudotime")]
plot_matrix_selected3<-plot_matrix_selected2[,which(colSums(plot_matrix_selected2)>0)]
data_expr_USM <-plot_matrix_selected3[which(rowSums(plot_matrix_selected3,na.rm = T)>0),]


plot_matrix_selected<-data_plot0_CTRL[!colnames(data_plot0_CTRL) %in% c("Treat","re_annotation_TFac","pseudotime")]
plot_matrix_selected2<-as.data.frame(t(plot_matrix_selected))
colData_CTRL<-data_plot0_CTRL[,c("Treat","re_annotation_TFac","pseudotime")]
plot_matrix_selected3 <-plot_matrix_selected2[,which(colSums(plot_matrix_selected2)>0)]
data_expr_CTRL<-plot_matrix_selected3[which(rowSums(plot_matrix_selected3,na.rm = T)>0),]

#conflict_prefer("cbind", "spam")
head(colData_USM);head(colData_CTRL)
colData_all<-rbind(colData_CTRL,colData_USM)
data_expr_USM[1:4,1:4];data_expr_CTRL[1:4,1:4]
data_expr_USM<-data_expr_USM[rownames(data_expr_CTRL),]
data_expr_all<-cbind(data_expr_CTRL,data_expr_USM)


dim(colData_all);dim(data_expr_all)#9420
#pheatmap(data_expr_all, cluster_rows=T, cluster_cols =F,scale="row",color = colorRampPalette(colors = c("blue","white","red"))(100))

data_expr2<-data_expr_all
data_expr2[data_expr2==0]<-NA
data_expr2[1:4,1:4]
sum(is.na(as.matrix(dist(data_expr2))))
giveNAs = which(is.na(as.matrix(dist(data_expr2))),arr.ind=TRUE)
#          row col

heat_plot1<-pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row",show_rownames = T,show_colnames = F,annotation_col=colData_all,color = colorRampPalette(colors = c("blue","white","red"))(100))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/ALL_CRB2EVT_RNA_K27ac_expression_pattern_no_NA_heatmaps.pdf",width = 25,height =10)
print(heat_plot1)
dev.off()
heat_plot2<-pheatmap(data_expr2, cluster_rows=T, cluster_cols =F,scale="row",show_rownames = T,show_colnames = F,annotation_col=colData_all,color = colorRampPalette(colors = c("blue","white","red"))(100))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/ALL_CRB2EVT_RNA_K27ac_expression_pattern_NA_heatmaps.pdf",width = 25,height =10)
print(heat_plot2)
dev.off()

###plot using complexheatmap
data_trans<-data_expr2
head(colData_all)
all(rownames(colData_all) %in% colnames(data_trans))
unique(colData_all$Treat)
table(colData_all$re_annotation_TFac)

anno_cell<-as.character(colData_all$re_annotation_TFac)
seudotime = as.numeric(colData_all$pseudotime)
ann_colors = list(re_annotation_TFac=c(CTBs_1=my_morandi_colors[1],CTBs_2=my_morandi_colors[2],EVTs_1=my_morandi_colors[6],
                                       EVTs_2=my_morandi_colors[7],EVTs_3=my_morandi_colors[8],other =my_morandi_colors[4]))

col_fun= colorRamp2(seq(from=0,to=24,length=length(seudotime)),colorRampPalette(colors = c("yellow", "red", "purple4"))(length(seudotime)))

column_ha = HeatmapAnnotation(pseudotime =seudotime,re_annotation_TFac = anno_cell,
                              col = list(pseudotime=col_fun,re_annotation_TFac = ann_colors),
                              na_col = "black")
#anno_cell2<-c(rep("ZGA_pre",24),rep("ZGA_after",36),rep("Differention",30))

plot_matrix <- t(apply(data_expr2,1,function(x){(x-mean(x,,na.rm = T))/sd(x,na.rm = T)}))
plot_matrix[1:4,1:4]
dend = as.dendrogram(hclust(dist(plot_matrix),method="ward.D2"))
##try different cluster_number
for ( d_number in 8:16){
  
  d_num<-d_number
  #d_num<-9
  print(d_num)
  dend = color_branches(dend, k =d_num)
  htkm <- Heatmap(plot_matrix,
                  #name= "z-score", 
                  #border = TRUE,
                  # top_annotation = column_ha,
                  na_col = "grey",
                  #column_split = factor(anno_cell2, levels = c("ZGA_pre","ZGA_after","Differention")),
                  # column_gap = unit(c(4), "mm"),
                  cluster_column_slices = FALSE,
                  column_title = "Gene expression in Human USM_CRB2EVT for RNA_K27ac_overlapped genes", 
                  col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                  show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                  row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,#cluster_rows= FALSE,
                  cluster_rows = dend, 
                  row_split = d_num,row_gap = unit(rep(2, d_num-1), "mm"),
                  row_dend_width = unit(4, "cm"))
  print(htkm)
  
  pdf(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Cluster_",d_num,"_ALL_CRB2EVT_RNA_K27ac_expression_pattern_NA_heatmaps.pdf"),height=25,width=10)
  print(htkm)
  dev.off()
  
  ###step three :: plot trend line for gene expression of each cluster
  #exact genes names in each cluster
  conflict_prefer("rbind", "spam")
  clusterlist = row_order(htkm)
  names(clusterlist)<-1:d_num
  htkm_module <- lapply(names(clusterlist), function(i){
    out <- data.frame(GeneID = rownames(plot_matrix[clusterlist[[i]],]),Cluster = paste0("cluster", i),stringsAsFactors = FALSE)
    return(out)}) %>%   do.call(rbind, .)
  
  head(htkm_module);colnames(htkm_module)<-c("row_genes","Cluster")
  table(htkm_module$Cluster)
  #cluster1 cluster2 cluster3 cluster4 
  #  158       80       83      179
  head(plot_matrix)
  write.table(htkm_module, file = paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/gene_Cluster_",d_num,"_ALL_CRB2EVT_RNA_K27ac.txt"), quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
  
  plot_matrix[1:4,1:4]
  data_plot<-data.frame(plot_matrix[htkm_module$row_genes,])
  data_plot$type<-as.character(htkm_module$Cluster)
  cell_names<-colnames(data_plot)
  
  conflict_prefer("melt", "reshape2")
  data_plot[1:4,1:4]
  data_preplot <- melt(data_plot,variable.name="cell",value.name = "expression",id.vars = "type")
  head(data_preplot)
  data_preplot2<-aggregate(expression ~ type+cell, data = data_preplot, mean)#median
  data_preplot2$cell<-gsub("[.]","-",data_preplot2$cell)
  head(data_preplot2);head(colData_all)
  colData_all$cell<-rownames(colData_all)
  colData_all$cell_order<-1:nrow(colData_all)
  
  data_preplot3=merge(data_preplot2,colData_all,by="cell")
  range(data_preplot3$expression)#-4.512997  5.022262
  data_preplot3<-data_preplot3[order(data_preplot3$re_annotation_TFac),]
  #data_preplot3<-data_preplot3[order(data_preplot3$re_annotation_TFac),]
  #data_preplot3$re_annotation_TFac<-factor(data_preplot3$re_annotation_TFac,levels = unique(data_preplot3$re_annotation_TFac))
  
  head(data_preplot3)
  
  library(moRandi)
  x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
  my_morandi_colors<-morandi_diy(my_colors = x)
  show_col(my_morandi_colors)
  annote_color<-my_morandi_colors[c(1:2,6:8,4)]
  
  data_plot<-data_preplot3[,c("expression","Treat","re_annotation_TFac","pseudotime","type")]
  colnames(data_plot)<-c("expression","Treat","cell_names","pseudotime","type")
  data_plot<-data_plot[order(data_plot$pseudotime,decreasing = F),]
  head(data_plot)
  method_type<-"loess"
  
  trend_plot22<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_rug(aes(colour = cell_names),sides="b") +
    #  scale_colour_manual(values=annote_color)+
    # geom_point(aes(colour = Treat),size=1,shape=19,alpha = 0.5)+ 
    scale_colour_manual(values=c("red",my_morandi_colors[1:2],"blue",my_morandi_colors[6:8],my_morandi_colors[4]))+
    # scale_y_continuous(limits=c(0,3),breaks=c(0,0.3,0.6,0.9,1.2,1.5,1.8))+
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5), panel.grid=element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0("all:",method_type))
  trend_plot3<-trend_plot22+ facet_wrap(~type,scales="free_y",ncol = 4)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),method=method_type,se=TRUE)
  #gam glm loess
  trend_plot3 
  ggsave(trend_plot3,file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Cluster_",d_num,"_CTBs2EVTs_CTRL_USM_trend_line_loess.pdf"),width = 20, height =20)
}


###using USM CTRL toghether cluster
d_num<-16
select_RNA_gene_module<-read.table(file =  paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/gene_Cluster_",d_num,"_ALL_CRB2EVT_RNA_K27ac.txt"), sep = '\t', row.names = 1, header = T)
head(select_RNA_gene_module)
select_RNA_gene_module$module<-select_RNA_gene_module$Cluster
table(select_RNA_gene_module$Cluster)
select_RNA_gene_module[which(select_RNA_gene_module$module %in% c("cluster5","cluster6","cluster7")),]$module<-"cluster5"
select_RNA_gene_module[which(select_RNA_gene_module$module %in% c("cluster8", "cluster9","cluster10","cluster11")),]$module<-"cluster6"
select_RNA_gene_module[which(select_RNA_gene_module$module %in% c("cluster12", "cluster13","cluster14", "cluster15", "cluster16")),]$module<-"cluster7"
select_RNA_gene_module$module<-factor(select_RNA_gene_module$module,levels = paste0("cluster",1:7))
select_RNA_gene_module<-select_RNA_gene_module[order(select_RNA_gene_module$module,decreasing = F),]
write.table(select_RNA_gene_module, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_rearrange_gene_Cluster_ALL_CRB2EVT_RNA_K27ac.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

###plot modules heatmaps for rearrange
RNA_all_gene_module<-read.table(file ="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_rearrange_gene_Cluster_ALL_CRB2EVT_RNA_K27ac.txt", sep = '\t', row.names = 1, header = T)
dim(RNA_all_gene_module)
RNA_all_gene_module$module<-factor(RNA_all_gene_module$module,levels = paste0("cluster",1:7))
RNA_all_gene_module<-RNA_all_gene_module[order(RNA_all_gene_module$module,decreasing = F),]
table(RNA_all_gene_module$module)
#cluster1 cluster2 cluster3 cluster4 cluster5 cluster6 cluster7 
#    8       31        9       51       51       96      255 

gap_number0<-as.numeric(table(RNA_all_gene_module$module))
gap_number<-gap_number0[-length(gap_number0)]
head(RNA_all_gene_module)
tail(RNA_all_gene_module)

plot_matrix2<-plot_matrix[as.character(RNA_all_gene_module$row_genes),]
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE, na_col = "grey",
                # row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(c("Early_time","Mid_low","Later_time","Mid_high"),gap_number0),
                row_gap = unit(rep(2, length(gap_number)), "mm"), 
                row_split =rep( paste0("cluster",1:7),gap_number0),
                cluster_column_slices = FALSE,column_title = "K27ac RNA overlapped Gene Expression in direction of CTBs2EVTs", 
                #   col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/rearrange_Cluster_ALL_CRB2EVT_RNA_K27ac_expression_pattern_NA_heatmaps_6clusters.pdf",height=25,width=10)
print(htkm)
dev.off()

###step three :: plot trend line for gene expression of each cluster
#exact genes names in each cluster
plot_matrix[1:4,1:4]
data_plot<-data.frame(plot_matrix[RNA_all_gene_module$row_genes,])
data_plot$type<-as.character(RNA_all_gene_module$module)
cell_names<-colnames(data_plot)

#conflict_prefer("melt", "reshape2")
data_plot[1:4,1:4]
data_preplot <- melt(data_plot,variable.name="cell",value.name = "expression",id.vars = "type")
head(data_preplot)

#data_preplot2<-aggregate(expression ~ type+cell, data = data_preplot, mean)#median
data_preplot2<-aggregate(expression ~ type+cell, data = data_preplot, median)#mean
data_preplot2$cell<-gsub("[.]","-",data_preplot2$cell)
head(data_preplot2);head(colData_all)
colData_all$cell<-rownames(colData_all)
colData_all$cell_order<-1:nrow(colData_all)

data_preplot3=merge(data_preplot2,colData_all,by="cell")
range(data_preplot3$expression)#-5.152833  5.022262
data_preplot3<-data_preplot3[order(data_preplot3$re_annotation_TFac),]
#data_preplot3<-data_preplot3[order(data_preplot3$re_annotation_TFac),]
#data_preplot3$re_annotation_TFac<-factor(data_preplot3$re_annotation_TFac,levels = unique(data_preplot3$re_annotation_TFac))

head(data_preplot3)

library(moRandi)
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
annote_color<-my_morandi_colors[c(1:2,6:8,4)]

data_plot<-data_preplot3[,c("expression","Treat","re_annotation_TFac","pseudotime","type")]
colnames(data_plot)<-c("expression","Treat","cell_names","pseudotime","type")
data_plot<-data_plot[order(data_plot$pseudotime,decreasing = F),]
head(data_plot)
method_type<-"loess"
data_plot$type<-factor(data_plot$type,levels = paste0("cluster",1:7))

trend_plot22<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
  geom_rug(aes(colour = cell_names),sides="b") +
  #  scale_colour_manual(values=annote_color)+
  # geom_point(aes(colour = Treat),size=1,shape=19,alpha = 0.5)+ 
  scale_colour_manual(values=c("red",my_morandi_colors[1:2],"blue",my_morandi_colors[6:8],my_morandi_colors[4]))+
  # scale_y_continuous(limits=c(0,3),breaks=c(0,0.3,0.6,0.9,1.2,1.5,1.8))+
  theme_bw()+ theme(plot.title = element_text(hjust = 0.5), panel.grid=element_blank())+
  labs(x = "pseudotime", y = "median of LogNormalized count in cluster", title =paste0("all:",method_type))
trend_plot3<-trend_plot22+ facet_wrap(~type,scales="free_y",ncol = 4)+
  geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),method=method_type,se=TRUE)
#gam glm loess
trend_plot3 
ggsave(trend_plot3,file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/rearrange_Cluster_ALL_CTBs2EVTs_CTRL_USM_trend_line_loess_6_cluster.pdf",width = 20, height =12)


###合并两组学变化趋势table
## reading gene in diff clusters
k27_USM_gene_module<-read.csv(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/k27_CTB2EVT_usm_gene_4cluster.csv", sep = ',', row.names = 1, header = F)
k27_USM_gene_module$V2<-paste0("cluster",k27_USM_gene_module$V2)
colnames(k27_USM_gene_module)<-c("k27_USM_cluster")

k27_CTRL_gene_module<-read.csv(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/k27_CTB2EVT_ct_gene_4cluster.csv", sep = ',', row.names = 1, header = F)
k27_CTRL_gene_module$V2<-paste0("cluster",k27_CTRL_gene_module$V2)
colnames(k27_CTRL_gene_module)<-c("k27_CTRL_cluster")

RNA_all_gene_module<-read.table(file ="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_rearrange_gene_Cluster_ALL_CRB2EVT_RNA_K27ac.txt", sep = '\t', row.names = 1, header = T)
head(RNA_all_gene_module)

rownames(RNA_all_gene_module)<-RNA_all_gene_module$row_genes
RNA_all_gene_module<-RNA_all_gene_module[,-1]
colnames(RNA_all_gene_module)<-c("RNA_cluster","RNA_module")

dim(k27_USM_gene_module);dim(k27_CTRL_gene_module);dim(RNA_all_gene_module)
k27_gene_module<-merge(k27_USM_gene_module,k27_CTRL_gene_module,by=0)
head(k27_gene_module)
rownames(k27_gene_module)<-k27_gene_module$Row.names
k27_gene_module<-k27_gene_module[,-1]
All_gene_module<-merge(k27_gene_module,RNA_all_gene_module,by=0)
rownames(All_gene_module)<-All_gene_module$Row.names
All_gene_module<-All_gene_module[,-1]
head(All_gene_module)


All_gene_module$k27_USM_cluster<-ifelse(All_gene_module$k27_USM_cluster =="cluster1","Up_Down_Up",
                                        ifelse(All_gene_module$k27_USM_cluster =="cluster2","Up_along",
                                               ifelse(All_gene_module$k27_USM_cluster =="cluster3","Down_Up_Down","Down_along")))
All_gene_module$k27_CTRL_cluster<-ifelse(All_gene_module$k27_CTRL_cluster =="cluster1","Up_along",
                                         ifelse(All_gene_module$k27_CTRL_cluster =="cluster2","Up_pre_along",
                                                ifelse(All_gene_module$k27_CTRL_cluster =="cluster3","Down_Up_Down","Down_along")))
gene_signature0<-c("KRT19","MMP15","SERPINE1","MALAT1","TGFB1")
All_gene_module[gene_signature0,]

##rename RNA CLUSTER
#All_gene_module$RNA_cluster  <-ifelse(All_gene_module$RNA_cluster   =="cluster1","min_Up_Down_Up",
#                                         ifelse(All_gene_module$RNA_cluster   =="cluster2","Down_along",
#                                                ifelse(All_gene_module$RNA_cluster   =="cluster3","Down_Up","Up_Down_Up")))




##总结类群数目
All_gene_module$all_class<-paste(All_gene_module$k27_USM_cluster,All_gene_module$k27_CTRL_cluster,All_gene_module$RNA_module,sep="-")
table(All_gene_module$all_class)
write.table(All_gene_module, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_trend_annotation_for_gene_overlap_CRB2EVT_RNA_K27ac.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##reading abortion related genes
Abortion_related_genes<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Abortion_regulated_gene_list_new.txt",header =T,sep="\t",check.names = F)
head(Abortion_related_genes);dim(Abortion_related_genes)
mark_gene <- as.character(Abortion_related_genes$gene_names);length(mark_gene)
gene_pos <- which(rownames(All_gene_module) %in% mark_gene)
rownames(All_gene_module)[gene_pos]
# "CCR7"     "MALAT1"   "SERPINE1" "STAT3"    "TGFB1"
gene_signature0<-rownames(All_gene_module)[gene_pos]

##绘制单基因假时序
All_gene_module<-read.table(file ="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_trend_annotation_for_gene_overlap_CRB2EVT_RNA_K27ac.txt", sep = '\t', row.names = 1, header = T)
gene_signature1<-c("ACTN4","CEBPB","COL4A2","EFNA1","MYH9","ARHGDIB","ASCL2","KRT19","ORC6","SNAI1","SOCS3","SPINT2","STK4","TPD52L1","ACAN","ADAM12","ADAM19","ADAM9","ANGPT4")
gene_signature2<-c("ADAM12","AKAP12","AURKA","AURKB","CD74","COL4A2","DNMT1","DUSP1","EDF1","EEF1E1","EFNA1","FLT1","GADD45B","GMNN","GTSE1","HMGB1","HSP90AB1","IQGAP1","ITGA5","KLF6","NPM1","NRN1","NUSAP1","PTMS","PTPRF","QSOX1","RPL29","RRM1","SERF2","SPINT1","TCF7L2","TOP2A","WWC3")
gene_signature<-unique(c(gene_signature0,gene_signature1,gene_signature2))
table(All_gene_module$RNA_module)
#cluster1 cluster2 cluster3 cluster4 cluster5 cluster6 cluster7 
#    8       31        9       51       51       94      253 

#gene_signature<-c("ARHGDIB", "ASCL2","KRT19","PGF","TPD52L1")
length(gene_signature)#54
All_gene_module[gene_signature,]

library(moRandi)
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
annote_color<-my_morandi_colors[c(1:2,6:8,4)]

data_plot0<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/All_Trophoblast_CTBs2EVTs_k27ac_RNA_overlapped_gene_expression.txt", sep = '\t', row.names = 1, header = T)
head(data_plot0)
data_preplot <- melt(data_plot0,id.vars=c("Treat","re_annotation_TFac","pseudotime"),variable.name="gene",value.name = "expression")
head(data_preplot)

data_plot_used<-data_preplot[which(data_preplot$gene %in% gene_signature),]
head(data_plot_used)

colnames(data_plot_used)<-c("Treat","cell_names","pseudotime","gene","expression")
data_plot_used<-data_plot_used[order(data_plot_used$pseudotime,decreasing = F),]

data_plot_used<-data_plot_used[which(data_plot_used$expression>0),]

method_type<-"loess"

trend_plot22<-ggplot(data_plot_used,aes(x= pseudotime,y=expression,colour=Treat))+
  geom_rug(aes(colour = cell_names),sides="b") +
  #  scale_colour_manual(values=annote_color)+
  # geom_point(aes(colour = Treat),size=0.5,shape=19,alpha = 0.2)+ 
  scale_colour_manual(values=c("red",my_morandi_colors[1:2],"blue",my_morandi_colors[6:8],my_morandi_colors[4]))+
  # scale_y_continuous(limits=c(0,3),breaks=c(0,0.3,0.6,0.9,1.2,1.5,1.8))+
  theme_bw()+ theme(plot.title = element_text(hjust = 0.5), panel.grid=element_blank())+
  labs(x = "pseudotime", y = "LogNormalized count", title =paste0("all:",method_type))
trend_plot3<-trend_plot22+ facet_wrap(~gene,scales="free_y",ncol = 7)+
  geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),method=method_type,se=TRUE)
#gam glm loess
trend_plot3
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/seperative_gene_for_each_module_CTBs2EVTs_CTRL_USM_trend_line_loess.pdf",trend_plot3,width=35, height=40)
#ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Abortion_seperative_gene_for_each_module_CTBs2EVTs_CTRL_USM_trend_line_loess_no_point.pdf",trend_plot3,width=26, height=5)

#ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/SERPINE1_for_CTBs2EVTs_CTRL_USM_trend_line_loess_point.pdf",trend_plot3,width=5, height=5)


##重新绘制带注释的热图
head(colData_all);plot_matrix[1:4,1:4];dim(plot_matrix)
colData_all$Treat<-as.character(colData_all$Treat)
colData_all[which(colData_all$Treat =="Abortion"),]$Treat<-"USM"

##get gaps number for row
###plot modules heatmaps for rearrange
RNA_all_gene_module<-read.table(file ="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_rearrange_gene_Cluster_ALL_CRB2EVT_RNA_K27ac.txt", sep = '\t', row.names = 1, header = T)
dim(RNA_all_gene_module)
RNA_all_gene_module$row_genes<-as.character(RNA_all_gene_module$row_genes)
#RNA_all_gene_module$module<-factor(RNA_all_gene_module$module,levels = paste0("cluster",1:11))
RNA_all_gene_module<-RNA_all_gene_module[order(RNA_all_gene_module$module,decreasing = F),]
table(RNA_all_gene_module$module)
#cluster1 cluster2 cluster3 cluster4 cluster5 cluster6 cluster7 
#    8       31        9       51       51       96      255 

groups<-names(table(RNA_all_gene_module$module))
gap_number0<-as.numeric(table(RNA_all_gene_module$module))
gap_number<-gap_number0[-length(gap_number0)]

head(RNA_all_gene_module)
##building expression matrix
plot_matrix2<-plot_matrix[RNA_all_gene_module$row_genes,];dim(plot_matrix2)#501 9420
plot_matrix2[1:4,1:4]

#set annotation for each cells  in column
head(colData_all);dim(colData_all)
colData_all[,c("re_annotation_TFac","Treat")]
colData_all$re_annotation_TFac<-as.character(colData_all$re_annotation_TFac)
table(colData_all$re_annotation_TFac)
#CTBs_1 CTBs_2 EVTs_1 EVTs_2 EVTs_3  other 
# 5076    568   1567   1649    442    118 
table(colData_all$Treat)
#CTRL  USM 
#8070 1350 

column_colors = list(group=c(USM=my_morandi_colors[1],CTRL=my_morandi_colors[21]),
                     cell_types=c(CTBs_1=my_morandi_colors[1],CTBs_2=my_morandi_colors[2],EVTs_1=my_morandi_colors[6],EVTs_2=my_morandi_colors[7],EVTs_3=my_morandi_colors[8],other =my_morandi_colors[4]))

column_ha = HeatmapAnnotation( group=as.character(colData_all$Treat),
                               cell_types = as.character(colData_all$re_annotation_TFac),
                               col = column_colors)

#set annotation for each genes in row
DEGs_class<-RNA_all_gene_module$module
ann_colors = list(DEGs_class=c(cluster1=ppCor[1], cluster2=ppCor[2], cluster3=ppCor[3], cluster4=ppCor[4], 
                               cluster5=ppCor[5], cluster6=ppCor[6], cluster7=ppCor[7], cluster8=ppCor[8]))

##reading abortion related genes
Abortion_related_genes<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Abortion_regulated_gene_list_new.txt",header =T,sep="\t",check.names = F)
head(Abortion_related_genes);dim(Abortion_related_genes)
mark_gene <- as.character(Abortion_related_genes$gene_names);length(mark_gene)
##or
gene_signature1<-c("ACTN4","CEBPB","COL4A2","EFNA1","MYH9","ARHGDIB","ASCL2","KRT19","ORC6","SNAI1","SOCS3","SPINT2","STK4","TPD52L1","ACAN","ADAM12","ADAM19","ADAM9","ANGPT4")
gene_signature2<-c("ADAM12","AKAP12","AURKA","AURKB","CD74","COL4A2","DNMT1","DUSP1","EDF1","EEF1E1","EFNA1","FLT1","GADD45B","GMNN","GTSE1","HMGB1","HSP90AB1","IQGAP1","ITGA5","KLF6","NPM1","NRN1","NUSAP1","PTMS","PTPRF","QSOX1","RPL29","RRM1","SERF2","SPINT1","TCF7L2","TOP2A","WWC3")
gene_signature<-c(gene_signature1,gene_signature2)

mark_gene_used<-unique(c(mark_gene,gene_signature));length(mark_gene_used)#173

gene_pos <- which(rownames(plot_matrix2) %in% mark_gene_used)
mark_gene2<-rownames(plot_matrix2[gene_pos,])
length(mark_gene2)#54

head(RNA_all_gene_module)
RNA_all_gene_module[which(RNA_all_gene_module$row_genes %in% mark_gene2),]


row_mark <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene2))
row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)

##plot rearrange heatmap
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,top_annotation = column_ha, 
                right_annotation =row_mark,left_annotation =row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                column_split =rep(c("CTRL","USM"),c("8070","1350")),
                cluster_column_slices = FALSE,column_title = "Overlap genes Developmental Gene Expression in direction of CTBs2EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                # col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_Abortion_and_GC_selected_genes_in_K27ac_RNA_overlapped)CTBs2EVTs_rearrange_module_annotation_add2.pdf",height=15,width=12)
print(htkm)
dev.off()

##
#GO BP enrichment
allcluster_BP_GO=data.frame();allcluster_BP_GO_whole=data.frame()

for (cluster_name in levels(gene_module$type)) {
  #cluster_name<-"Mid_low"
  print(paste0("cluster number: ",cluster_name))
  small_gene_group=unique(as.character(gene_module[which(gene_module$type == cluster_name),]$id))
  gene_change=bitr(small_gene_group, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  BP_GO <- enrichGO(gene= unique(gene_change$ENTREZID),OrgDb= org.Hs.eg.db,keyType= 'ENTREZID',ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff= 0.1,readable= TRUE)
  BP_GO2 <- simplify(BP_GO, cutoff=0.7, by="p.adjust", select_fun=min)
  go_res=BP_GO2@result
  go_res_whole=BP_GO@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=cluster_name;allcluster_BP_GO=rbind(allcluster_BP_GO,go_res)}
  if (dim(go_res_whole)[1] != 0) {
    go_res_whole$cluster=cluster_name;allcluster_BP_GO_whole=rbind(allcluster_BP_GO_whole,go_res_whole)}
}
head(allcluster_BP_GO[,c("ID","Description","qvalue","cluster")])
table(allcluster_BP_GO$cluster)
# Early Early_Mid     Later       Mid Mid_Later 
#   9        97       141        51       184 
write.table(allcluster_BP_GO, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_for_gene_module_along_CTBs2EVTs_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(allcluster_BP_GO_whole, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_whole_for_module_along_CTBs2EVTs_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

