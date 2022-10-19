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
library(dendextend)
library(ComplexHeatmap)
library(scales)
library(ggsci)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:46))]

##step 2
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)
target_final$re_annotation_TFac_Treat<- paste0(target_final$re_annotation_TFac,"_",target_final$Treat)

DefaultAssay(target_final) <- "RNA"
target_final <- NormalizeData(target_final, verbose = FALSE)
Idents(object = target_final) <- 're_annotation_TFac_Treat'
head(target_final@meta.data)

#get overlapped gene between developped genes and abortion DEGs 
merge_data0 <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge_new.txt", header=T)
head(merge_data0)

merge_data0$cluster<- factor(merge_data0$cluster,levels=unique(merge_data0$cluster))
merge_data0$trend <- factor(merge_data0$trend,levels=c("Up","Down"))
merge_data0<-merge_data0[order(merge_data0$cluster,merge_data0$trend,decreasing = F),]

##construction matrix of averageExpression 
maker_gene<-unique(as.character(merge_data0$gene))
AverageExp<-AverageExpression(target_final,assays="RNA",slot = "data",features=maker_gene)

typeof(AverageExp);head(AverageExp$RNA)
range(AverageExp$RNA)# 0.00 2770.25
AverageExp_data<-log2(AverageExp$RNA+1)
head(AverageExp_data);range(AverageExp_data)# 0.00000 11.43632
colnames(AverageExp_data)
#AverageExp_data2<-AverageExp_data[maker_gene,c(paste0(levels(merge_data0$cluster),"_CTRL"),paste0(levels(merge_data0$cluster),"_Abortion"))]

subgroup_order1<-c("CTBs_1","CTBs_1","CTBs_2","CTBs_2","STBs_1","STBs_1","STBs_2","STBs_2","STBs_3","STBs_3","EVTs_1","EVTs_1",
                   "EVTs_2","EVTs_2","EVTs_3","EVTs_3","Endo","Endo","STCs","STCs","FBs","FBs","Mycs","Mycs","HCs","HCs","NKs","NKs","Ts","Ts","Epi","Epi","Bs","Bs")
subgroup_order2<-c("CTBs_1_Abortion","CTBs_1_CTRL","CTBs_2_Abortion","CTBs_2_CTRL","STBs_1_Abortion","STBs_1_CTRL",
"STBs_2_Abortion","STBs_2_CTRL","STBs_3_Abortion","STBs_3_CTRL","EVTs_1_Abortion","EVTs_1_CTRL",
"EVTs_2_Abortion","EVTs_2_CTRL","EVTs_3_Abortion","EVTs_3_CTRL","Endo_Abortion","Endo_CTRL","STCs_Abortion","STCs_CTRL",
"FBs_Abortion","FBs_CTRL","Mycs_Abortion","Mycs_CTRL","HCs_Abortion","HCs_CTRL","NKs_Abortion","NKs_CTRL",
"Ts_Abortion","Ts_CTRL","Epi_Abortion","Epi_CTRL","Bs_Abortion","Bs_CTRL")
AverageExp_data2<-AverageExp_data[maker_gene,subgroup_order2]
head(AverageExp_data2)
#set annotation for each cells  in column
column_colors = list(cell_types=c(CTBs_1=ppCor_all[1],CTBs_2=ppCor_all[2],STBs_1=ppCor_all[3],STBs_2=ppCor_all[4],STBs_3=ppCor_all[5],
                       EVTs_1=ppCor_all[6],EVTs_2=ppCor_all[7],EVTs_3=ppCor_all[8],Endo=ppCor_all[9],STCs=ppCor_all[10],FBs=ppCor_all[11],
                       Mycs=ppCor_all[12],HCs=ppCor_all[13],NKs=ppCor_all[14],Ts=ppCor_all[15],Epi=ppCor_all[16],Bs=ppCor_all[17]),group=c(Abortion=ppCor[1],CTRL=ppCor[2]))
column_ha = HeatmapAnnotation(cell_types = subgroup_order1,group=rep(c("CTRL","Abortion"),17),col = column_colors)

##building expression matrix
AverageExp_data_scal<-scale(t(AverageExp_data2), center=T,scale=T)
range(AverageExp_data_scal)#-3.939897  5.659453
#color_set = colorRamp2(seq(min(AverageExp_data_scal), max(AverageExp_data_scal), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
color_set = colorRamp2(c(-5, 0, 5), c("navy", "white", "red"))

AverageExp_data_scal[1:4,1:4]

#mark_gene <- as.character(Abortion_related_genes$gene_names);length(mark_gene)
#gene_pos <- which(colnames(AverageExp_data_scal) %in% mark_gene)
#mark_gene2<-colnames(AverageExp_data_scal[,gene_pos])
#row_mark <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene2))
#row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)

##plot rearrange heatmap
htkm <- Heatmap(t(AverageExp_data_scal),name= "z-score",border = TRUE,top_annotation = column_ha, 
               # left_annotation = row_ha,right_annotation = row_mark,
                row_gap = unit(rep(2, length(subgroup_order1)/2), "mm"), row_split =rep(groups,length(subgroup_order1)/2-1),
                cluster_column_slices = FALSE,column_title = "Expression mean for all Abortion related DEGs", 
                col= color_set,
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/heatmap_for_Expression_mean_overlap_all_DEGs.pdf",height=8,width=8)
print(htkm)
dev.off()
