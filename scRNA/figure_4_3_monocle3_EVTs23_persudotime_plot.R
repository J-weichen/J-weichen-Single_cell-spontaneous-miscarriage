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
library(tidyverse)
library(patchwork)
library(scales)
library(ggsci)
library(ggplot2)
library(harmony)
library(ggbeeswarm)
library(ggridges)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dendextend)
library(pheatmap)

##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)

pal <- pal_npg("nrc", alpha=0.5)(9)
nejm<-pal_nejm("default",alpha =0.5)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)
show_col(pal)

cds_branch<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_EVT_selected_branch_monocle3_object.rds")
pseu_time<-data.frame(pseudotime(cds_branch))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)
#pseu_time$number<-ifelse(rownames(pseu_time) %in% rownames(colData(cds_branch)),1,-1)
pseu_time$cell_name<-rownames(pseu_time)
pseu_time<-pseu_time[order(pseu_time$pseudotime,decreasing = F),]
head(pseu_time);tail(pseu_time)
#reading development genes
gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_gene_module_along_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
head(gene_module)
gene_module$id<-as.character(gene_module$id)
gene_module<-gene_module[order(gene_module$module,decreasing = F),]
#get metadata
branch_meta<-as.data.frame(colData(cds_branch))
head(branch_meta);dim(branch_meta)
#get expression matrix
pt.matrix <- exprs(cds_branch)[match(gene_module$id,rownames(rowData(cds_branch))),]
pt.matrix[1:4,1:4];dim(pt.matrix)#1296 8232
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#1296 8232

##reorder branch cells
branch_meta<-branch_meta[colnames(pt.matrix),]
head(branch_meta);dim(branch_meta)

##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix[1:4,1:4];dim(pt.matrix)#1296 8232

#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

#get overlapped gene between developped genes and abortion DEGs 
merge_data0 <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge_new.txt", header=T)
head(merge_data0)
subgroup_order2<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3")
merge_data1<-distinct(merge_data0[which(merge_data0$cluster %in% subgroup_order2),c("gene","cluster_trend")])
merge_data1$cluster_trend<-factor(merge_data1$cluster_trend,levels =c(paste0(subgroup_order2,"_Up"),paste0(subgroup_order2,"_Down")))
Abortion_EVTs_2_3_DEGs<-unique(as.character(merge_data1[which(merge_data1$cluster_trend %in% c("EVTs_2_Up","EVTs_3_Up","EVTs_2_Down","EVTs_3_Down")),]$gene))
overlap_genes<-intersect(Abortion_EVTs_2_3_DEGs,unique(gene_module$id))

##rearange genes by different classes
row_anno2 <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/DEGs_EVT2_EVT3_overlapped_development_genes_along_sleceted_EVT_branchment_pseudotime.txt", header=T)
head(row_anno2)
##get gaps number for row
groups<-names(table(row_anno2$Cluster))
gap_number0<-as.numeric(table(row_anno2$Cluster))
gap_number<-gap_number0[-length(gap_number0)]

##building expression matrix
plot_matrix2<-plot_matrix[overlap_genes,];dim(plot_matrix2)#577 8232
plot_matrix2<-plot_matrix2[as.character(row_anno2$row_genes),]
plot_matrix2[1:4,1:4]

#set annotation for each cells  in column
head(branch_meta);dim(branch_meta)
branch_meta[,c("re_annotation_TFac2","Treat")]
branch_meta$re_annotation_TFac2<-as.character(branch_meta$re_annotation_TFac2)
branch_meta$re_annotation_TFac<-as.character(branch_meta$re_annotation_TFac)
branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac == "EVTs_2","EVTs_2",ifelse(branch_meta$re_annotation_TFac == "EVTs_3","EVTs_3","others"))

column_colors = list(group=c(Abortion=my_morandi_colors[1],CTRL=my_morandi_colors[21]),
  cell_types=c(EVTs_2=my_morandi_colors[3],EVTs_3=my_morandi_colors[9],others="darkgray"))
column_ha = HeatmapAnnotation( group=as.character(branch_meta$Treat),
                               cell_types = as.character(branch_meta$re_annotation_TFac2),
                               col = column_colors)

#set annotation for each genes in row
DEGs_class<-row_anno2$trend
ann_colors = list(DEGs_class=c(Both_Down=ppCor[1],Both_Up=ppCor[2],Down_Up=ppCor[3],EVTs_2_Down=ppCor[9],EVTs_2_Up=ppCor[5],EVTs_3_Down=ppCor[6],EVTs_3_Up=ppCor[7],Up_Down=ppCor[8]))

##reading abortion related genes
Abortion_related_genes<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Abortion_regulated_gene_list_new.txt",header =T,sep="\t",check.names = F)
head(Abortion_related_genes);dim(Abortion_related_genes)
mark_gene <- as.character(Abortion_related_genes$gene_names);length(mark_gene)
gene_pos <- which(rownames(plot_matrix2) %in% mark_gene)
mark_gene2<-rownames(plot_matrix2[gene_pos,])

row_mark <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene2))
row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)

#EVTs_2_Up<-row_anno2$EVTs_2_Up;EVTs_3_Up<-row_anno2$EVTs_3_Up;EVTs_2_Down<-row_anno2$EVTs_2_Down;EVTs_3_Down<-row_anno2$EVTs_3_Down
#ann_colors = list(EVTs_2_Up=c(EVTs_2_Up=ppCor[1],No="darkgray"),EVTs_3_Up=c(EVTs_3_Up=ppCor[2],No="darkgray"),
#                  EVTs_2_Down=c(EVTs_2_Down=ppCor[3],No="darkgray"),EVTs_3_Down=c(EVTs_3_Down=ppCor[5],No="darkgray"),
#                  DEGs_class=c(Both_Down=ppCor[1],Both_Up=ppCor[2],Down_Up=ppCor[3],EVTs_2_Down=ppCor[4],EVTs_2_Up=ppCor[5],EVTs_3_Down=ppCor[6],EVTs_3_Up=ppCor[7],Up_Down=ppCor[8]))
#row_ha = rowAnnotation(EVTs_2_Up =EVTs_2_Up,EVTs_3_Up =EVTs_3_Up,EVTs_2_Down =EVTs_2_Down,EVTs_3_Down =EVTs_3_Down,DEGs_class=DEGs_class,col = ann_colors)

row_anno2[which(row_anno2$row_genes %in% mark_gene2),]
##plot rearrange heatmap
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,top_annotation = column_ha, 
                right_annotation =row_mark,left_annotation =row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of EVTs_2 -> EVTs_3", 
               # col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "PuOr"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_overlap_DEGs_selected_direction_EVTs_rearrange_module_annotation_add2.pdf",height=8,width=7)
print(htkm)
dev.off()

##split matrix based on treat group
#get expression matrix
pt.matrix <- exprs(cds_branch)[match(gene_module$id,rownames(rowData(cds_branch))),]
pt.matrix[1:4,1:4];dim(pt.matrix)#1296 8232
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#1296 8232

CTRL_meta<-branch_meta[which(branch_meta$Treat == "CTRL"),]
pt.matrix_ctrl<-pt.matrix[,rownames(CTRL_meta)]

##recalculation for smooth expression along psudotime
pt.matrix_ctrl <- t(apply(pt.matrix_ctrl,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix_ctrl <- t(apply(pt.matrix_ctrl,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix_ctrl[1:4,1:4];dim(pt.matrix_ctrl)#1296 8232
#plot gene trend in pheatmap
plot_matrix<-pt.matrix_ctrl
colnames(plot_matrix) = 1:ncol(plot_matrix)

#set annotation for each cells  in column
column_colors = list(cell_types=c(EVTs_2=my_morandi_colors[3],EVTs_3=my_morandi_colors[9],others="darkgray"))
column_ha = HeatmapAnnotation(cell_types = as.character(CTRL_meta$re_annotation_TFac2),col = column_colors)

#set annotation for each genes in row
DEGs_class<-row_anno2$trend
ann_colors = list(DEGs_class=c(Both_Down=ppCor[1],Both_Up=ppCor[2],Down_Up=ppCor[3],EVTs_2_Down=ppCor[9],EVTs_2_Up=ppCor[5],EVTs_3_Down=ppCor[6],EVTs_3_Up=ppCor[7],Up_Down=ppCor[8]))
row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)
##building expression matrix
plot_matrix2<-plot_matrix[overlap_genes,];dim(plot_matrix2)#577 8232
plot_matrix2<-plot_matrix2[as.character(row_anno2$row_genes),]

##plot rearrange heatmap
htkm_ctrl <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,top_annotation = column_ha, 
                     right_annotation = row_ha,#left_annotation = ha2,
                     row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                     cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of EVTs_2 -> EVTs_3", 
                     col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                     show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                     row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm_ctrl)


##For Abortion
Abortion_meta<-branch_meta[which(branch_meta$Treat == "Abortion"),]
pt.matrix_abortion<-pt.matrix[,rownames(Abortion_meta)]

##recalculation for smooth expression along psudotime
pt.matrix_abortion <- t(apply(pt.matrix_abortion,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix_abortion <- t(apply(pt.matrix_abortion,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix_abortion[1:4,1:4];dim(pt.matrix_abortion)#1296 8232
#plot gene trend in pheatmap
plot_matrix<-pt.matrix_abortion
colnames(plot_matrix) = 1:ncol(plot_matrix)

#set annotation for each cells  in column
column_colors = list(cell_types=c(EVTs_2=my_morandi_colors[3],EVTs_3=my_morandi_colors[9],others="darkgray"))
column_ha = HeatmapAnnotation(cell_types = as.character(Abortion_meta$re_annotation_TFac2),col = column_colors)

#set annotation for each genes in row
DEGs_class<-row_anno2$trend
ann_colors = list(DEGs_class=c(Both_Down=ppCor[1],Both_Up=ppCor[2],Down_Up=ppCor[3],EVTs_2_Down=ppCor[9],EVTs_2_Up=ppCor[5],EVTs_3_Down=ppCor[6],EVTs_3_Up=ppCor[7],Up_Down=ppCor[8]))
row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)

##building expression matrix
plot_matrix2<-plot_matrix[overlap_genes,];dim(plot_matrix2)#577 4182
plot_matrix2<-plot_matrix2[as.character(row_anno2$row_genes),]

##plot rearrange heatmap
#colorRamp2，出自circlize包，可以根据你指定的几个颜色，生成一组渐变色。
htkm_abortion <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,top_annotation = column_ha, 
                         right_annotation = row_ha,#left_annotation = ha2,
                         row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                         cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of EVTs_2 -> EVTs_3", 
                         col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                         show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                         row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm_abortion)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/CTRL_abortion_pheatmap_overlap_DEGs_selected_direction_EVTs_rearrange_module_annotation_add2.pdf",height=8,width=6)
print(htkm_ctrl)
print(htkm_abortion)
dev.off()

##step 2
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

##selected celltypes 
cell_name<-c("EVTs_2","EVTs_3")
selected_cell<-subset(x = target_final, subset = re_annotation_TFac %in% cell_name)
subgroup_order<-c("EVTs_2","EVTs_3")
selected_cell$re_annotation_TFac<- factor(selected_cell$re_annotation_TFac,levels = rev(subgroup_order))
selected_cell$re_annotation_TFac_Treat<- paste0(selected_cell$re_annotation_TFac,"_",selected_cell$Treat)
table(selected_cell$re_annotation_TFac_Treat)
##EVTs_2_Abortion     EVTs_2_CTRL EVTs_3_Abortion     EVTs_3_CTRL 
##        911            5794            4626             120 

##construction matrix of averageExpression of EVTs
DefaultAssay(selected_cell) <- "RNA"
selected_cell <- NormalizeData(selected_cell, verbose = FALSE)
Idents(object = selected_cell) <- 're_annotation_TFac_Treat'
head(selected_cell@meta.data)

##
maker_gene<-unique(as.character(row_anno2$row_genes))
AverageExp<-AverageExpression(selected_cell,assays="RNA",slot = "data",features=maker_gene)

typeof(AverageExp);head(AverageExp$RNA)
range(AverageExp$RNA)# 1.571088e-03 3.254462e+02
AverageExp_data<-log2(AverageExp$RNA+1)
range(AverageExp_data);head(AverageExp_data)
colnames(AverageExp_data)
subgroup_order2<-c("EVTs_2_CTRL","EVTs_3_CTRL","EVTs_2_Abortion","EVTs_3_Abortion")
AverageExp_data<-AverageExp_data[as.character(row_anno2$row_genes),subgroup_order2]
head(AverageExp_data)

#set annotation for each cells  in column
column_colors = list(group=c(Abortion=my_morandi_colors[1],CTRL=my_morandi_colors[21]),
                     cell_types=c(EVTs_2_CTRL=my_morandi_colors[3],EVTs_3_CTRL=my_morandi_colors[9],EVTs_2_Abortion=my_morandi_colors[3],EVTs_3_Abortion=my_morandi_colors[9]))
column_ha = HeatmapAnnotation(group=c("CTRL","CTRL","Abortion","Abortion"),cell_types = c("EVTs_2_CTRL","EVTs_3_CTRL","EVTs_2_Abortion","EVTs_3_Abortion"),col = column_colors)

##building expression matrix
AverageExp_data_scal<-scale(t(AverageExp_data), center=T,scale=T)
range(AverageExp_data_scal)#-1.499623  1.498822
color_set = colorRamp2(seq(min(AverageExp_data_scal), max(AverageExp_data_scal), length = 3), c("steelblue3", "#EEEEEE", "indianred2"), space = "RGB")
AverageExp_data_scal[1:4,1:4]

#mark_gene <- as.character(Abortion_related_genes$gene_names);length(mark_gene)
#gene_pos <- which(colnames(AverageExp_data_scal) %in% mark_gene)
#mark_gene2<-colnames(AverageExp_data_scal[,gene_pos])
row_mark <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene2))
row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)

##plot rearrange heatmap
htkm <- Heatmap(t(AverageExp_data_scal),name= "z-score",border = TRUE,top_annotation = column_ha, 
                left_annotation = row_ha,right_annotation = row_mark,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                cluster_column_slices = FALSE,column_title = "Expression mean for selected genes in direction of EVTs_2 -> EVTs_3", 
                col= color_set,
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/heatmap_forExpression_mean_overlap_DEGs_selected_direction_EVTs_rearrange_module_annotation_add2.pdf",height=8,width=7)
print(htkm)
dev.off()

##step3
##plot expression smooth line along persudo time for each DEG class in different clusters
pt.matrix0 <- exprs(cds_branch)[match(gene_module$id,rownames(rowData(cds_branch))),]
pt.matrix_pre <- t(apply(pt.matrix0,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix_pre[1:4,1:4]
head(row_anno2)
table(row_anno2$type)

data_plot<-data.frame(pt.matrix_pre[row_anno2$row_genes,])
#data_plot<-data.frame(pt.matrix0[row_anno2$row_genes,])
cell_names<-colnames(data_plot)
data_plot<-data_plot[row_anno2$row_genes,]
data_plot$type<-as.character(row_anno2$type)
#data_plot$gene<-as.character(row_anno2$type)

data_plot[1:4,1:4]
data_preplot <- melt(data_plot,variable.name="cell",value.name = "expression",id.vars = "type")
head(data_preplot)
data_preplot2<-aggregate(expression ~ type+cell, data = data_preplot, mean)
data_preplot2$cell<-gsub("[.]","-",data_preplot2$cell)
head(data_preplot2)

metadata<-as.data.frame(colData(cds_branch))
pseudotime <- pseudotime(cds_branch, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(metadata)]
metadata$pseudotime <- pseudotime
metadata2 <- metadata[,c("Treat","re_annotation_TFac","pseudotime")]
table(metadata2$re_annotation_TFac)
#metadata2$re_annotation_TFac<-factor(metadata2$re_annotation_TFac,levels = c("EVTs_2","EVTs_3"))
metadata2$cell <-rownames(metadata2)
str(metadata2)
data_preplot3=merge(data_preplot2,metadata2,by="cell")
head(data_preplot3)
#data_preplot4<-data_preplot3[which(data_preplot3$re_annotation_TFac %in%  c("EVTs_2","EVTs_3")),]
#data_preplot4$re_annotation_TFac<-factor(data_preplot4$re_annotation_TFac,levels = c("EVTs_2","EVTs_3"))
data_preplot3$re_annotation_TFac<-factor(data_preplot3$re_annotation_TFac,levels = c("EVTs_2","EVTs_3"))

range(data_preplot3$expression)#-0.8526407 61.1569577

trend_plot<-ggplot(data_preplot3,aes(x= pseudotime,y=expression,colour = Treat))+ 
  geom_smooth(method = "loess") + 
  geom_rug(aes(color = re_annotation_TFac),sides="b") +
  scale_colour_manual(values=my_morandi_colors[c(1,21,3,9)])+ylim(-1,1)+
  theme_bw()+ theme(plot.title = element_text(hjust = 0.5),
                    axis.text=element_blank(),panel.grid=element_blank(),axis.ticks = element_blank())+
  xlab("pseudotime") + ylab("gene expression") + facet_wrap(~type,scales="free_y",ncol = 2)
trend_plot
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/expression_trends_along_pseudotime_of_selected_EVT_branch_for_cluster.png",trend_plot,width =10, height =30,limitsize = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/expression_trends_along_pseudotime_of_selected_EVT_branch_for_cluster.pdf",trend_plot,width =10, height =30,limitsize = FALSE)
