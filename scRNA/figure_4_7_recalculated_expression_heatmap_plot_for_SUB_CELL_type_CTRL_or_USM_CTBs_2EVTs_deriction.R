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

###analysis for CTBs2EVTs branch
cds_branch<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2EVTs_selected_branch_monocle3_object.rds")
plot_cells(cds_branch, color_cells_by = "re_annotation_TFac", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)

#基因表达趋势图
Track_genes<-c("MKI67","SERPINE1","TIMP3")
plot_genes_in_pseudotime(cds_branch[Track_genes,], color_cells_by="re_annotation_TFac", min_expr=0.5, ncol =1)

###reananlysis for choosen cells
USM_cds_subset = cds_branch[,colData(cds_branch)$Treat == "Abortion" ]
USM_cds_subset2 = USM_cds_subset[,colData(USM_cds_subset)$re_annotation_TFac %in% c("EVTs_1","EVTs_2","EVTs_3")]
USM_anno_plot0<-plot_cells(USM_cds_subset2, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
USM_cds_subset<-USM_cds_subset2

CTRL_cds_subset = cds_branch[,colData(cds_branch)$Treat == "CTRL"]
CTRL_cds_subset2 = CTRL_cds_subset[,colData(CTRL_cds_subset)$re_annotation_TFac %in% c("EVTs_1","EVTs_2")]
CTRL_anno_plot0<-plot_cells(CTRL_cds_subset2, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
CTRL_cds_subset<-CTRL_cds_subset2

##information substract for CTRL
meta_info<-as.data.frame(colData(CTRL_cds_subset))
pseu_time<-data.frame(pseudotime(CTRL_cds_subset))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)
#pseu_time$number<-ifelse(rownames(pseu_time) %in% rownames(colData(cds_branch)),1,-1)
pseu_time$cell_name<-rownames(pseu_time)
pseu_time<-pseu_time[order(pseu_time$pseudotime,decreasing = F),]
head(pseu_time);tail(pseu_time);dim(pseu_time)
dim(meta_info)
meta_info2<-meta_info[rownames(pseu_time),]
meta_info3<-merge(meta_info2,pseu_time,by=0)
head(meta_info3)
write.table(meta_info3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_pseudotime_info_for_cells_in_CTRL_EVT1_to_EVT2_branch.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##绘制变化动态
gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/developped_genes_rearranged_along_CTBs2EVTs_branchment_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
head(gene_module)
head(gene_module)
gene_module$module<-factor(gene_module$module,levels = c(8,7,1,6,5,2,3,4,9))
gene_module$type<-factor(gene_module$type,levels = c("Early","Early_Mid","Mid","Mid_Later","N_Later"))
gene_module<-gene_module[order(gene_module$type,gene_module$module,decreasing = F),]
gene_module$id<-as.character(gene_module$id)
head(gene_module)
table(gene_module$type)

#get metadata
branch_meta<-as.data.frame(colData(CTRL_cds_subset))
head(branch_meta);dim(branch_meta)

#get expression matrix
pt.matrix <- exprs(CTRL_cds_subset)[match(gene_module$id,rownames(rowData(CTRL_cds_subset))),]
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 7339
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 7339
pt.matrix<-pt.matrix[as.character(gene_module$id),]

##reorder branch cells
branch_meta<-branch_meta[colnames(pt.matrix),]
head(branch_meta);dim(branch_meta)#  7339   67

##target data for yaxi
Yaxi_gene_matrix1<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix1[,1:4];dim(Yaxi_gene_matrix1)# 3 7339
##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix_yaxi<-pt.matrix
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix[1:4,1:4];dim(pt.matrix)# 1174 7339

##target genes for yaxi
Yaxi_gene_matrix2<-pt.matrix_yaxi[which(rownames(pt.matrix_yaxi) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix2[,1:4];dim(Yaxi_gene_matrix2)#3 7339

Yaxi_gene_matrix3<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix3[,1:4];dim(Yaxi_gene_matrix3)#3 7339

write.table(Yaxi_gene_matrix1, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/CTRL_EVT1_to_EVT2_raw_selected_genes_expression_for_Yaxi1.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix2, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/CTRL_EVT1_to_EVT2_smooth_spline_selected_genes_expression_for_Yaxi2.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/CTRL_EVT1_to_EVT2_smooth_spline_scale_selected_genes_expression_for_Yaxi3.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##expression merge for coldata

#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

##get gaps number for row
groups<-names(table(gene_module$type))
gap_number0<-as.numeric(table(gene_module$type))
gap_number<-gap_number0[-length(gap_number0)]

##building expression matrix
plot_matrix2<-plot_matrix[gene_module$id,];dim(plot_matrix2)#1174 7339
plot_matrix2[1:4,1:4]

#set annotation for each cells  in column
head(branch_meta);dim(branch_meta)
branch_meta$re_annotation_TFac2<-as.character(branch_meta$re_annotation_TFac2)
branch_meta$re_annotation_TFac<-as.character(branch_meta$re_annotation_TFac)
table(branch_meta$re_annotation_TFac)
#EVTs_1 EVTs_2 
# 1683   5656 
#branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_2","STBs_1","STBs_2","STBs_3"),"others",branch_meta$re_annotation_TFac)
#branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("STBs_1","STBs_2","STBs_3"),"others",
#                                        ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_1","CTBs_2"),"CTBs","EVTs"))
table(branch_meta$re_annotation_TFac2)
#EVTs_1 EVTs_2 
# 1683   5656  

column_colors = list(#group=c(Abortion=my_morandi_colors[1],CTRL=my_morandi_colors[21]),
  cell_types=c(EVTs_1=my_morandi_colors[10],EVTs_2=my_morandi_colors[15],others="gray"))
column_ha = HeatmapAnnotation( #group=as.character(branch_meta$Treat),
  cell_types = as.character(branch_meta$re_annotation_TFac2),
  col = column_colors)

#set annotation for each genes in row
DEGs_class<-gene_module$type
ann_colors = list(DEGs_class=c(Early=ppCor[1],Early_Mid=ppCor[2],Mid=ppCor[3],Mid_Later=ppCor[9],N_Later=ppCor[5]))

##reading abortion related genes
Abortion_related_genes<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Abortion_regulated_gene_list_new.txt",header =T,sep="\t",check.names = F)
head(Abortion_related_genes);dim(Abortion_related_genes)
mark_gene <- as.character(Abortion_related_genes$gene_names);length(mark_gene)
gene_pos <- which(rownames(plot_matrix2) %in% mark_gene)
mark_gene2<-rownames(plot_matrix2[gene_pos,])

row_mark <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene2))
row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)

##plot rearrange heatmap
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,top_annotation = column_ha, 
                right_annotation =row_mark,left_annotation =row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of CTBs2EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_CTRL_EVT1_to_EVT2_rearrange_module_annotation_add2.pdf",height=8,width=7)
print(htkm)
dev.off()


##for USM
##information substract for USM

meta_info<-as.data.frame(colData(USM_cds_subset))
pseu_time<-data.frame(pseudotime(USM_cds_subset))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)
#pseu_time$number<-ifelse(rownames(pseu_time) %in% rownames(colData(cds_branch)),1,-1)
pseu_time$cell_name<-rownames(pseu_time)
pseu_time<-pseu_time[order(pseu_time$pseudotime,decreasing = F),]
head(pseu_time);tail(pseu_time);dim(pseu_time)
dim(meta_info)
meta_info2<-meta_info[rownames(pseu_time),]
meta_info3<-merge(meta_info2,pseu_time,by=0)
write.table(meta_info3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_pseudotime_info_for_cells_in_USM_EVT1EVT3_branch.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##绘制变化动态
gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/developped_genes_rearranged_along_CTBs2EVTs_branchment_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
head(gene_module)
head(gene_module)
gene_module$module<-factor(gene_module$module,levels = c(8,7,1,6,5,2,3,4,9))
gene_module$type<-factor(gene_module$type,levels = c("Early","Early_Mid","Mid","Mid_Later","N_Later"))
gene_module<-gene_module[order(gene_module$type,gene_module$module,decreasing = F),]
gene_module$id<-as.character(gene_module$id)
head(gene_module)
table(gene_module$type)

#get metadata
branch_meta<-as.data.frame(colData(USM_cds_subset))
head(branch_meta);dim(branch_meta)

#get expression matrix
pt.matrix <- exprs(USM_cds_subset)[match(gene_module$id,rownames(rowData(USM_cds_subset))),]
pt.matrix[1:4,1:4];dim(pt.matrix)# 1174 4954
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 4954
pt.matrix<-pt.matrix[as.character(gene_module$id),]

##reorder branch cells
branch_meta<-branch_meta[colnames(pt.matrix),]
head(branch_meta);dim(branch_meta)#  4954   67

##target data for yaxi
Yaxi_gene_matrix1<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix1[,1:4];dim(Yaxi_gene_matrix1)# 3 4954
##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix_yaxi<-pt.matrix
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix[1:4,1:4];dim(pt.matrix)#  3 4954

##target genes for yaxi
Yaxi_gene_matrix2<-pt.matrix_yaxi[which(rownames(pt.matrix_yaxi) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix2[,1:4];dim(Yaxi_gene_matrix2)# 3 4954

Yaxi_gene_matrix3<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix3[,1:4];dim(Yaxi_gene_matrix3)# 3 4954

write.table(Yaxi_gene_matrix1, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/USM_EVT1EVT3_raw_selected_genes_expression_for_Yaxi1.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix2, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/USM_EVT1EVT3_smooth_spline_selected_genes_expression_for_Yaxi2.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/USM_EVT1EVT3_smooth_spline_scale_selected_genes_expression_for_Yaxi3.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##expression merge for coldata

#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

##get gaps number for row
groups<-names(table(gene_module$type))
gap_number0<-as.numeric(table(gene_module$type))
gap_number<-gap_number0[-length(gap_number0)]

##building expression matrix
plot_matrix2<-plot_matrix[gene_module$id,];dim(plot_matrix2)# 1174 4954
plot_matrix2[1:4,1:4]

#set annotation for each cells  in column
head(branch_meta);dim(branch_meta)
branch_meta$re_annotation_TFac2<-as.character(branch_meta$re_annotation_TFac2)
branch_meta$re_annotation_TFac<-as.character(branch_meta$re_annotation_TFac)
table(branch_meta$re_annotation_TFac)
#EVTs_1 EVTs_2 EVTs_3 
# 79    796   4079 
#branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_2","STBs_1","STBs_2","STBs_3"),"others",branch_meta$re_annotation_TFac)
#branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("STBs_1","STBs_2","STBs_3"),"others",
#                                        ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_1","CTBs_2"),"CTBs","EVTs"))
table(branch_meta$re_annotation_TFac2)
#EVTs_1 EVTs_2 EVTs_3 
# 79    796   4079

column_colors = list(#group=c(Abortion=my_morandi_colors[1],CTRL=my_morandi_colors[21]),
  cell_types=c(EVTs_1=my_morandi_colors[10],EVTs_2=my_morandi_colors[15],EVTs_3=my_morandi_colors[2]))
column_ha = HeatmapAnnotation( #group=as.character(branch_meta$Treat),
  cell_types = as.character(branch_meta$re_annotation_TFac2),
  col = column_colors)

#set annotation for each genes in row
DEGs_class<-gene_module$type
ann_colors = list(DEGs_class=c(Early=ppCor[1],Early_Mid=ppCor[2],Mid=ppCor[3],Mid_Later=ppCor[9],N_Later=ppCor[5]))

##reading abortion related genes
Abortion_related_genes<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Abortion_regulated_gene_list_new.txt",header =T,sep="\t",check.names = F)
head(Abortion_related_genes);dim(Abortion_related_genes)
mark_gene <- as.character(Abortion_related_genes$gene_names);length(mark_gene)
gene_pos <- which(rownames(plot_matrix2) %in% mark_gene)
mark_gene2<-rownames(plot_matrix2[gene_pos,])

row_mark <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene2))
row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)

##plot rearrange heatmap
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,top_annotation = column_ha, 
                right_annotation =row_mark,left_annotation =row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of CTBs2EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_USM_EVT1_to_EVT3_rearrange_module_annotation_add2.pdf",height=8,width=7)
print(htkm)
dev.off()
