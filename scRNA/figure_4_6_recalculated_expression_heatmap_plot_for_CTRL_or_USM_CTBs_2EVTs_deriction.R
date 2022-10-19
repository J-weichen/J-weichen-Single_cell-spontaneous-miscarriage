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

###reananlysis for choosen cells
USM_cds_subset = cds_branch[,colData(cds_branch)$Treat == "Abortion"]
cds_subset <- USM_cds_subset
#cds_subset <- choose_cells(cds_second)
cds_subset <- cluster_cells(cds_subset)
cds_subset <- learn_graph(cds_subset)
cds_subset <- order_cells(cds_subset) 
subset_anno_plot0<-plot_cells(cds_subset, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
subset_anno_plot1<-plot_cells(cds_subset, color_cells_by = "Treat", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
subset_anno_plot2<-plot_cells(cds_subset, color_cells_by = "Phase", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
subset_anno_plot3<-plot_cells(cds_subset, color_cells_by = "re_annotation_TFac", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
subset_anno_plot4<-plot_cells(cds_subset, color_cells_by = "sample_code", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
subset_anno_plot5<-plot_cells(cds_subset, color_cells_by = "sample", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
#ggarrange(plot_j_gene_TRA_sig2, plot_j_gene_TRB_sig2,heights = c(0.7, 0.3),widths=c(1,1),ncol = 1, nrow = 2)
subset_anno_plot<-CombinePlots(plots = list(subset_anno_plot0,subset_anno_plot1,subset_anno_plot2,subset_anno_plot3,subset_anno_plot4,subset_anno_plot5),ncol = 3,legend = NULL)
USM_anno_plot<-subset_anno_plot
USM_cds_subset<-cds_subset

##for CTRL object
CTRL_cds_subset = cds_branch[,colData(cds_branch)$Treat == "CTRL"]
cds_subset <- CTRL_cds_subset
cds_subset <- cluster_cells(cds_subset)
cds_subset <- learn_graph(cds_subset)
cds_subset <- order_cells(cds_subset) 
subset_anno_plot0<-plot_cells(cds_subset, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
subset_anno_plot1<-plot_cells(cds_subset, color_cells_by = "Treat", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
subset_anno_plot2<-plot_cells(cds_subset, color_cells_by = "Phase", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
subset_anno_plot3<-plot_cells(cds_subset, color_cells_by = "re_annotation_TFac", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
subset_anno_plot4<-plot_cells(cds_subset, color_cells_by = "sample_code", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
subset_anno_plot5<-plot_cells(cds_subset, color_cells_by = "sample", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
#ggarrange(plot_j_gene_TRA_sig2, plot_j_gene_TRB_sig2,heights = c(0.7, 0.3),widths=c(1,1),ncol = 1, nrow = 2)
subset_anno_plot<-CombinePlots(plots = list(subset_anno_plot0,subset_anno_plot1,subset_anno_plot2,subset_anno_plot3,subset_anno_plot4,subset_anno_plot5),ncol = 3,legend = NULL)
subset_anno_plot
CTRL_anno_plot<-subset_anno_plot
CTRL_cds_subset<-cds_subset

combined_anno_plot<-USM_anno_plot/CTRL_anno_plot
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Split_CTB2EVTs_monocle3_pseudotime.pdf",combined_anno_plot,width = 15, height =20)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Split_CTB2EVTs_monocle3_pseudotime.png",combined_anno_plot,width = 15, height =20)

saveRDS(CTRL_cds_subset, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/subset_CTRL_CTB2EVTs_monocle3_pseudotime_object_reanalysis.rds")
saveRDS(USM_cds_subset, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/subset_USM_CTB2EVTs_monocle3_pseudotime_object_reanalysis.rds")


##substract the expression and plot heatmap  for CTRL or USM
###recalculated for expression ternd along pesudotime of target genes 
CTRL_CTBs2EVTs_branch_plot<-plot_cells(CTRL_cds_subset,  reduction_method="UMAP",color_cells_by = "pseudotime", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/recalculated_CTRL_CTBs2EVTs_branch_umap.pdf",CTRL_CTBs2EVTs_branch_plot,width = 8, height =7)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/recalculatedCTRL_CTBs2EVTs_branch_umap.png",CTRL_CTBs2EVTs_branch_plot,width = 8, height =7)

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
write.table(meta_info3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_pseudotime_info_for_cells_in_CTRL_CTBs2EVTs_branch.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

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
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 12710
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 12710
pt.matrix<-pt.matrix[as.character(gene_module$id),]

##reorder branch cells
branch_meta<-branch_meta[colnames(pt.matrix),]
head(branch_meta);dim(branch_meta)# 18330    67

##target data for yaxi
Yaxi_gene_matrix1<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix1[,1:4];dim(Yaxi_gene_matrix1)#1174 18330
##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix_yaxi<-pt.matrix
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix[1:4,1:4];dim(pt.matrix)# 1174 18330

##target genes for yaxi
Yaxi_gene_matrix2<-pt.matrix_yaxi[which(rownames(pt.matrix_yaxi) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix2[,1:4];dim(Yaxi_gene_matrix2)#1174 18330

Yaxi_gene_matrix3<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix3[,1:4];dim(Yaxi_gene_matrix3)#1174 18330

write.table(Yaxi_gene_matrix1, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/CTRL_CTBs2EVTs_raw_selected_genes_expression_for_Yaxi1.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix2, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/CTRL_CTBs2EVTs_smooth_spline_selected_genes_expression_for_Yaxi2.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/CTRL_CTBs2EVTs_smooth_spline_scale_selected_genes_expression_for_Yaxi3.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##expression merge for coldata

#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

##get gaps number for row
groups<-names(table(gene_module$type))
gap_number0<-as.numeric(table(gene_module$type))
gap_number<-gap_number0[-length(gap_number0)]

##building expression matrix
plot_matrix2<-plot_matrix[gene_module$id,];dim(plot_matrix2)#577 8232
plot_matrix2[1:4,1:4]

#set annotation for each cells  in column
head(branch_meta);dim(branch_meta)
branch_meta$re_annotation_TFac2<-as.character(branch_meta$re_annotation_TFac2)
branch_meta$re_annotation_TFac<-as.character(branch_meta$re_annotation_TFac)
table(branch_meta$re_annotation_TFac)
#branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_2","STBs_1","STBs_2","STBs_3"),"others",branch_meta$re_annotation_TFac)
branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("STBs_1","STBs_2","STBs_3"),"others",
                                        ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_1","CTBs_2"),"CTBs","EVTs"))
table(branch_meta$re_annotation_TFac2)
#CTBs   EVTs others 
#5647  12410    273 

column_colors = list(#group=c(Abortion=my_morandi_colors[1],CTRL=my_morandi_colors[21]),
  cell_types=c(CTBs=my_morandi_colors[10],EVTs=my_morandi_colors[15],others="gray"))
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
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_CTRL_CTBs2EVTs_rearrange_module_annotation_add2.pdf",height=8,width=7)
print(htkm)
dev.off()


##for USM
###recalculated for expression ternd along pesudotime of target genes 
USM_CTBs2EVTs_branch_plot<-plot_cells(USM_cds_subset,  reduction_method="UMAP",color_cells_by = "pseudotime", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/recalculated_USM_CTBs2EVTs_branch_umap.pdf",USM_CTBs2EVTs_branch_plot,width = 8, height =7)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/recalculatedUSM_CTBs2EVTs_branch_umap.png",USM_CTBs2EVTs_branch_plot,width = 8, height =7)

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
write.table(meta_info3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_pseudotime_info_for_cells_in_USM_CTBs2EVTs_branch.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

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
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 12710
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 12710
pt.matrix<-pt.matrix[as.character(gene_module$id),]

##reorder branch cells
branch_meta<-branch_meta[colnames(pt.matrix),]
head(branch_meta);dim(branch_meta)# 18330    67

##target data for yaxi
Yaxi_gene_matrix1<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix1[,1:4];dim(Yaxi_gene_matrix1)#1174 18330
##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix_yaxi<-pt.matrix
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix[1:4,1:4];dim(pt.matrix)# 1174 18330

##target genes for yaxi
Yaxi_gene_matrix2<-pt.matrix_yaxi[which(rownames(pt.matrix_yaxi) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix2[,1:4];dim(Yaxi_gene_matrix2)#1174 18330

Yaxi_gene_matrix3<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix3[,1:4];dim(Yaxi_gene_matrix3)#1174 18330

write.table(Yaxi_gene_matrix1, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/USM_CTBs2EVTs_raw_selected_genes_expression_for_Yaxi1.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix2, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/USM_CTBs2EVTs_smooth_spline_selected_genes_expression_for_Yaxi2.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/USM_CTBs2EVTs_smooth_spline_scale_selected_genes_expression_for_Yaxi3.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##expression merge for coldata

#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

##get gaps number for row
groups<-names(table(gene_module$type))
gap_number0<-as.numeric(table(gene_module$type))
gap_number<-gap_number0[-length(gap_number0)]

##building expression matrix
plot_matrix2<-plot_matrix[gene_module$id,];dim(plot_matrix2)# 1174 5620
plot_matrix2[1:4,1:4]

#set annotation for each cells  in column
head(branch_meta);dim(branch_meta)
branch_meta$re_annotation_TFac2<-as.character(branch_meta$re_annotation_TFac2)
branch_meta$re_annotation_TFac<-as.character(branch_meta$re_annotation_TFac)
table(branch_meta$re_annotation_TFac)
#branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_2","STBs_1","STBs_2","STBs_3"),"others",branch_meta$re_annotation_TFac)
branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("STBs_1","STBs_2","STBs_3"),"others",
                                        ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_1","CTBs_2"),"CTBs","EVTs"))
table(branch_meta$re_annotation_TFac2)
#CTBs   EVTs others 
#505   4954    161 

column_colors = list(#group=c(Abortion=my_morandi_colors[1],CTRL=my_morandi_colors[21]),
  cell_types=c(CTBs=my_morandi_colors[10],EVTs=my_morandi_colors[15],others="gray"))
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
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_USM_CTBs2EVTs_rearrange_module_annotation_add2.pdf",height=8,width=7)
print(htkm)
dev.off()
