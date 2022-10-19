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

#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
#major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
#subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
#target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
#target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)
##calculation 
#DefaultAssay(target_final) <- "RNA"
# Normalize RNA data for visualization purposes
#target_final <- NormalizeData(target_final, verbose = FALSE)
##for maternal and fetal gene special argument
#gene_signature<-c("HLA-DRB1","HLA-DRA","HLA-DR","FOLR2","CD206")
#plot_maker <- FeaturePlot(object = target_final, features = gene_signature,cols= c("grey", "purple"),ncol=3)
#plot_maker

###analysis for CTBs2EVTs branch
cds_branch<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2EVTs_selected_branch_monocle3_object.rds")
plot_cells(cds_branch, color_cells_by = "re_annotation_TFac", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)

##information substract monocle3
pseu_time<-data.frame(pseudotime(cds_branch))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)

##提取目标轨迹seurat3对象
Troph_target<-subset(x = target_final, cells =rownames(pseu_time))
DimPlot(Troph_target, group.by = "re_annotation_TFac",label = F,cols = ppCor)

##calculation 
DefaultAssay(Troph_target) <- "RNA"
# Normalize RNA data for visualization purposes
Troph_target <- NormalizeData(Troph_target, verbose = FALSE)
#express_data <- as.matrix(GetAssayData(Troph_target, slot = "data"))
conflict_prefer("which", "Matrix")
## reading gene
gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_CTBs2EVTs_gene_module_along_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
gene_module$id<-as.character(gene_module$id)
Trans_K27ac<-read.table(file ="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Tans_k27_CTB2EVT_gene.txt", sep = ' ', row.names = 1, header = TRUE)
head(Trans_K27ac)
gene_module_2<-gene_module[which(gene_module$id %in% as.character(Trans_K27ac$x)),]

##提取表达矩阵
exprs <- data.frame(FetchData(object = Troph_target,slot = "data",vars = c(gene_module_2$id,"Treat","re_annotation_TFac")))
exprs[1:4,1:4]

##添加假时序信息
data_plot0<-merge(exprs,pseu_time,by=0)
Troph_group<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3")
data_plot0$re_annotation_TFac<- factor(data_plot0$re_annotation_TFac,levels = Troph_group)
data_plot0<-data_plot0[order(data_plot0$pseudotime,decreasing = F),]
data_plot0[1:4,]


##绘制时序拟合曲线
head(data_plot0)
data_plot0$re_annotation_TFac<-as.character(data_plot0$re_annotation_TFac)
data_plot0[which(!(data_plot0$re_annotation_TFac %in% c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3"))),]$re_annotation_TFac<-"other"
table(data_plot0$re_annotation_TFac)
data_plot0$re_annotation_TFac<-factor(data_plot0$re_annotation_TFac,levels = c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3","other"))
rownames(data_plot0)<-data_plot0$Row.names
data_plot0<-data_plot0[,-1]
data_plot_raw<-data_plot0

data_plot0<-data_plot0[which(data_plot0$pseudotime<24),]
data_plot0_CTRL<-data_plot0[which(data_plot0$Treat =="CTRL"),]
data_plot0_USM<-data_plot0[which(data_plot0$Treat =="Abortion"),]

data_plot0
write.table(data_plot_raw, file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/All_no_cut_pseudotime_Trophoblast_CTBs2EVTs_k27ac_RNA_overlapped_gene_expression.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(data_plot0, file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/All_Trophoblast_CTBs2EVTs_k27ac_RNA_overlapped_gene_expression.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(data_plot0_CTRL, file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/CTRL_All_Trophoblast_CTBs2EVTs_k27ac_RNA_overlapped_gene_expression.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(data_plot0_USM, file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/USM_All_Trophoblast_CTBs2EVTs_k27ac_RNA_overlapped_gene_expression.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
