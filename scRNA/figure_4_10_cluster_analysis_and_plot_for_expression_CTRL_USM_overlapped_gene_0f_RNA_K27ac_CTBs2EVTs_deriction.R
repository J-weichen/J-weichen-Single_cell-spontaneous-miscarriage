rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
dyn.load('/mnt/data/chenwei/software/share_road/libudunits2/libudunits2.so.0')
dyn.load('/mnt/data/chenwei/software/share_road/libproj/libproj.so.0')
dyn.load('/mnt/data/chenwei/software/share_road/libgdal.so.20')
dyn.load('/mnt/data/chenwei/software/share_road/libgeos-3.4.2.so')
dyn.load('/mnt/data/chenwei/software/share_road/libgeos_c.so.1')

library(Seurat)
library(monocle3)
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
dim(data_plot0_USM);dim(data_plot0_CTRL)# 1350  8070 
##for USM
plot_matrix_selected<-data_plot0_USM[!colnames(data_plot0_USM) %in% c("Treat","re_annotation_TFac","pseudotime")]
plot_matrix_selected2<-as.data.frame(t(plot_matrix_selected))
colData<-data_plot0_USM[,c("Treat","re_annotation_TFac","pseudotime")]
colData[1:4,]
plot_matrix_selected2[1:4,1:4]
plot_matrix_selected3 <-plot_matrix_selected2[,which(colSums(plot_matrix_selected2)>0)]

data_expr<-plot_matrix_selected3[which(rowSums(plot_matrix_selected3,na.rm = T)>0),]
dim(data_expr);dim(plot_matrix_selected2);dim(plot_matrix_selected3)

pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row")
pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row",color = colorRampPalette(colors = c("blue","white","red"))(100))

data_expr2<-data_expr
data_expr2[data_expr2==0]<-NA
data_expr2[1:4,1:4]
sum(is.na(as.matrix(dist(data_expr2))))#2
giveNAs = which(is.na(as.matrix(dist(data_expr2))),arr.ind=TRUE)
#          row col
#KCNQ1.AS1 272  83
#EPCAM      83 272
data_expr2[c(83,272),1:4]
##We get the rows out and start checking what to remove:
tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
  sum(is.na(as.matrix(dist(data_expr2[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]

data_expr2 = data_expr2[-as.numeric(rmv),]

pheatmap(data_expr2, cluster_rows=T, cluster_cols =F,na_col = "black",scale="row",show_rownames = T,show_colnames = F,annotation_col=colData,color = colorRampPalette(colors = c("blue","white","red"))(100))

heat_plot1<-pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row",show_rownames = T,show_colnames = F,annotation_col=colData,color = colorRampPalette(colors = c("blue","white","red"))(100))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/USM_CRB2EVT_RNA_K27ac_expression_pattern_no_NA_heatmaps.pdf",width = 20,height =8)
print(heat_plot1)
dev.off()
heat_plot2<-pheatmap(data_expr2, cluster_rows=T, cluster_cols =F,scale="row",show_rownames = T,show_colnames = F,annotation_col=colData,color = colorRampPalette(colors = c("blue","white","red"))(100))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/USM_CRB2EVT_RNA_K27ac_expression_pattern_NA_heatmaps.pdf",width = 20,height =8)
print(heat_plot2)
dev.off()

###plot using complexheatmap
data_trans<-data_expr2
head(colData)
all(rownames(colData) %in% colnames(data_trans))
unique(colData$Treat)
table(colData$re_annotation_TFac)

anno_cell<-as.character(colData$re_annotation_TFac)
seudotime = as.numeric(colData$pseudotime)
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
d_num<-4;dend = color_branches(dend, k =d_num)
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

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Cluster_USM_CRB2EVT_RNA_K27ac_expression_pattern_NA_heatmaps.pdf",height=20,width=8)
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

write.table(htkm_module, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/gene_Cluster_USM_CRB2EVT_RNA_K27ac.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

head(plot_matrix)
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
head(data_preplot2);head(colData)
colData$cell<-rownames(colData)
colData$cell_order<-1:nrow(colData)

data_preplot3=merge(data_preplot2,colData,by="cell")
range(data_preplot3$expression)#-1.612694  2.250575
data_preplot3<-data_preplot3[order(data_preplot3$re_annotation_TFac),]
#data_preplot3<-data_preplot3[order(data_preplot3$re_annotation_TFac),]
#data_preplot3$re_annotation_TFac<-factor(data_preplot3$re_annotation_TFac,levels = unique(data_preplot3$re_annotation_TFac))

head(data_preplot3)

library(moRandi)
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
annote_color<-my_morandi_colors[c(1:2,6:8,4)]

data_plot<-data_preplot3
method_type<-"loess"

trend_plot<-ggplot(data_plot,aes(x= cell_order,y=expression))+ 
  geom_point(aes(colour = re_annotation_TFac),size=1,shape=19,alpha = 0.5)+
  geom_rug(aes(color = re_annotation_TFac),sides="b") +  scale_colour_manual(values=annote_color)+#ylim(-2,7)+
  geom_smooth(method=method_type,se=TRUE)+ #gam glm loess
  theme_bw()+ theme(plot.title = element_text(hjust = 0.5),axis.text=element_blank(),panel.grid=element_blank(),axis.ticks = element_blank())+
  labs(x = "pseudotime", y = "scaled gene expression(LogNormalized count)", title =method_type)

trend_plot1<-trend_plot+ facet_wrap(~type,scales="free_y",ncol = 2)
trend_plot2<-trend_plot+ facet_wrap(~type,scales="free_y",ncol = 1)

ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/expression_trends_for_different_cluster_USM_CRB2EVT_RNA_K27ac_expression_pattern_NA_1.pdf",trend_plot1,width =7, height =6,limitsize = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/expression_trends_for_different_cluster_USM_CRB2EVT_RNA_K27ac_expression_pattern_NA_2.pdf",trend_plot2,width =7, height =18,limitsize = FALSE)

##for CTRL
library(moRandi)
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
annote_color<-my_morandi_colors[c(1:2,6:8,4)]

plot_matrix_selected<-data_plot0_CTRL[!colnames(data_plot0_CTRL) %in% c("Treat","re_annotation_TFac","pseudotime")]
plot_matrix_selected2<-as.data.frame(t(plot_matrix_selected))
colData<-data_plot0_CTRL[,c("Treat","re_annotation_TFac","pseudotime")]
colData[1:4,]
plot_matrix_selected2[1:4,1:4]
plot_matrix_selected3 <-plot_matrix_selected2[,which(colSums(plot_matrix_selected2)>0)]

data_expr<-plot_matrix_selected3[which(rowSums(plot_matrix_selected3,na.rm = T)>0),]
dim(data_expr);dim(plot_matrix_selected2);dim(plot_matrix_selected3)#8070

pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row")
pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row",color = colorRampPalette(colors = c("blue","white","red"))(100))

data_expr2<-data_expr
data_expr2[data_expr2==0]<-NA
data_expr2[1:4,1:4];dim(data_expr2)
sum(is.na(as.matrix(dist(data_expr2))))#2
giveNAs = which(is.na(as.matrix(dist(data_expr2))),arr.ind=TRUE)
#          row col

##We get the rows out and start checking what to remove:
##no remove

pheatmap(data_expr2, cluster_rows=T, cluster_cols =F,na_col = "black",scale="row",show_rownames = T,show_colnames = F,annotation_col=colData,color = colorRampPalette(colors = c("blue","white","red"))(100))

heat_plot1<-pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row",show_rownames = T,show_colnames = F,annotation_col=colData,color = colorRampPalette(colors = c("blue","white","red"))(100))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/CTRL_CRB2EVT_RNA_K27ac_expression_pattern_no_NA_heatmaps.pdf",width = 20,height =8)
print(heat_plot1)
dev.off()
heat_plot2<-pheatmap(data_expr2, cluster_rows=T, cluster_cols =F,scale="row",show_rownames = T,show_colnames = F,annotation_col=colData,color = colorRampPalette(colors = c("blue","white","red"))(100))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/CTRL_CRB2EVT_RNA_K27ac_expression_pattern_NA_heatmaps.pdf",width = 20,height =8)
print(heat_plot2)
dev.off()

###plot using complexheatmap
data_trans<-data_expr2
head(colData)
all(rownames(colData) %in% colnames(data_trans))
unique(colData$Treat)
table(colData$re_annotation_TFac)


library(circlize)


anno_cell<-as.character(colData$re_annotation_TFac)
seudotime = as.numeric(colData$pseudotime)
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
d_num<-4;dend = color_branches(dend, k =d_num)
htkm <- Heatmap(plot_matrix,
                #name= "z-score", 
                #border = TRUE,
                # top_annotation = column_ha,
                na_col = "grey",
                #column_split = factor(anno_cell2, levels = c("ZGA_pre","ZGA_after","Differention")),
                # column_gap = unit(c(4), "mm"),
                cluster_column_slices = FALSE,
                column_title = "Gene expression in Human CTRL_CRB2EVT for RNA_K27ac_overlapped genes", 
                col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,#cluster_rows= FALSE,
                cluster_rows = dend, 
                row_split = d_num,row_gap = unit(rep(2, d_num-1), "mm"),
                row_dend_width = unit(4, "cm"))
print(htkm)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Cluster_CTRL_CRB2EVT_RNA_K27ac_expression_pattern_NA_heatmaps.pdf",height=20,width=8)
print(htkm)
dev.off()

###step three :: plot trend line for gene expression of each cluster
#exact genes names in each cluster
clusterlist = row_order(htkm)
names(clusterlist)<-1:d_num
htkm_module <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(plot_matrix[clusterlist[[i]],]),Cluster = paste0("cluster", i),stringsAsFactors = FALSE)
  return(out)}) %>%   do.call(rbind, .)

head(htkm_module);colnames(htkm_module)<-c("row_genes","Cluster")
table(htkm_module$Cluster)
#cluster1 cluster2 cluster3 cluster4 
# 123       93      132      153 
write.table(htkm_module, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/gene_Cluster_CTRL_CRB2EVT_RNA_K27ac.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

head(plot_matrix)
plot_matrix[1:4,1:4]
data_plot<-data.frame(plot_matrix[htkm_module$row_genes,])
data_plot$type<-as.character(htkm_module$Cluster)
cell_names<-colnames(data_plot)

data_plot[1:4,1:4]
data_preplot <- melt(data_plot,variable.name="cell",value.name = "expression",id.vars = "type")
head(data_preplot)
data_preplot2<-aggregate(expression ~ type+cell, data = data_preplot, mean)#median
data_preplot2$cell<-gsub("[.]","-",data_preplot2$cell)
head(data_preplot2);head(colData)
colData$cell<-rownames(colData)
colData$cell_order<-1:nrow(colData)

data_preplot3=merge(data_preplot2,colData,by="cell")
range(data_preplot3$expression)#-1.612694  2.250575
data_preplot3<-data_preplot3[order(data_preplot3$re_annotation_TFac),]
#data_preplot3<-data_preplot3[order(data_preplot3$re_annotation_TFac),]
#data_preplot3$re_annotation_TFac<-factor(data_preplot3$re_annotation_TFac,levels = unique(data_preplot3$re_annotation_TFac))
head(data_preplot3)
data_plot<-data_preplot3


method_type<-"loess"
trend_plot<-ggplot(data_plot,aes(x= cell_order,y=expression))+ 
  geom_point(aes(colour = re_annotation_TFac),size=1,shape=19,alpha = 0.5)+
  geom_rug(aes(color = re_annotation_TFac),sides="b") +  scale_colour_manual(values=annote_color)+#ylim(-2,7)+
  geom_smooth(method=method_type,se=TRUE)+ #gam glm loess
  theme_bw()+ theme(plot.title = element_text(hjust = 0.5),axis.text=element_blank(),panel.grid=element_blank(),axis.ticks = element_blank())+
  labs(x = "pseudotime", y = "scaled gene expression(LogNormalized count)", title =method_type)

trend_plot1<-trend_plot+ facet_wrap(~type,scales="free_y",ncol = 2)
trend_plot2<-trend_plot+ facet_wrap(~type,scales="free_y",ncol = 1)

ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/expression_trends_for_different_cluster_CTRL_CRB2EVT_RNA_K27ac_expression_pattern_NA_1.pdf",trend_plot1,width =7, height =6,limitsize = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/expression_trends_for_different_cluster_CTRL_CRB2EVT_RNA_K27ac_expression_pattern_NA_2.pdf",trend_plot2,width =7, height =18,limitsize = FALSE)
