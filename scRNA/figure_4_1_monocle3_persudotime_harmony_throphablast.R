#Monocle2 Pseudotime analysis
#ref:https://nbisweden.github.io/workshop-scRNAseq/oldlabs/monocle_analysis.html
#ref:https://cloud.tencent.com/developer/article/1692225
#ref: https://www.sohu.com/a/327133920_278730
#ref：https://www.jianshu.com/p/5d6fd4561bc0
#ref:https://zhuanlan.zhihu.com/p/378365295
#integrated abject
#ref；https://github.com/satijalab/seurat/issues/1658
#https://github.com/cole-trapnell-lab/monocle3/issues/148
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

pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#reading final seurat object
#target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object.rds")
#major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
#subgroup_order0<-c("CTBs_1","CTBs_2","CTBs_3","CTBs_4","STBs_1","STBs_2","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
#target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
#target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)
#metadata_new<-target_final@meta.data
#metadata_new$re_annotation_TFac<-as.character(metadata_new$re_annotation_TFac)
#metadata_new[which(metadata_new$re_annotation_TFac %in% c("CTBs_2","CTBs_3")),]$re_annotation_TFac <- "CTBs_2"
#metadata_new[which(metadata_new$re_annotation_TFac == "STBs_2"),]$re_annotation_TFac <- "STBs_3"
#metadata_new[which(metadata_new$re_annotation_TFac =="STBs_1"),]$re_annotation_TFac <- "STBs_2"
#metadata_new[which(metadata_new$re_annotation_TFac =="CTBs_4"),]$re_annotation_TFac <- "STBs_1"
#target_final@meta.data<-metadata_new
#saveRDS(target_final,file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")

#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

cell_name<-"Trophoblast_KRT7"
selected_cell<-subset(x = target_final, subset = re_anno_TFac_major == cell_name)
table(selected_cell$sample_code)
#CTRL_1   CTRL_2   CTRL_3   CTRL_4   CTRL_5 Arrest_1 Arrest_2 Arrest_3 Arrest_4 
# 10970     6153     4139     2282     1063     1217     5128     1751      409 
subgroup_order1<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3")
selected_cell$re_annotation_TFac<-factor(target_final$re_annotation_TFac,levels=subgroup_order1,ordered=TRUE)
table(selected_cell$re_annotation_TFac)
#CTBs_1 CTBs_2 STBs_1 STBs_2 STBs_3 EVTs_1 EVTs_2 EVTs_3 
#7837   8746    372   1955    903   1848   6705   4746

## 子集去除批次效应
DefaultAssay(selected_cell) <- "RNA"
scRNAsub <- SCTransform(selected_cell) %>% RunPCA() 
anno_plot<-DimPlot(scRNAsub, reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle("selected trophoblast Intergration object")+ NoLegend()
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Trophoblast_KRT7_selected_annoplot.pdf", anno_plot,width=10, height=10)
saveRDS(scRNAsub, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_selected_seurat_object.rds")
#scRNAsub<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_selected_seurat_object.rds")

#set.seed(19921010)
pc.num=1:30
scRNA_list<-list();round_number<-c()
for (run in 2:15) {
 # run<-2#test line
  print(paste0("round: ",run))
  scRNAsub_temp<- RunHarmony(scRNAsub, group.by.vars = "sample_code",assay.use = "SCT", max.iter.harmony = run)  
  scRNAsub_temp <- RunUMAP(scRNAsub_temp, reduction="harmony", dims=pc.num) %>%
    FindNeighbors(reduction="harmony", dims=pc.num) %>% FindClusters(resolution=0.8) 
  scRNA_list<-c(scRNA_list,list(scRNAsub_temp))
  round_number<-c(round_number,paste0("round_",run))
  DimPlot(scRNAsub_temp, reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",run,"round"))+ NoLegend()
}
names(scRNA_list)<-round_number
saveRDS(scRNA_list, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_different_iter_harmony_seurat_object.rds")

hplot2<-DimPlot(scRNA_list[[1]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",2,"round"))+ NoLegend()
hplot3<-DimPlot(scRNA_list[[2]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",3,"round"))+ NoLegend()
hplot4<-DimPlot(scRNA_list[[3]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",4,"round"))+ NoLegend()
hplot5<-DimPlot(scRNA_list[[4]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",5,"round"))+ NoLegend()
hplot6<-DimPlot(scRNA_list[[5]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",6,"round"))+ NoLegend()
hplot7<-DimPlot(scRNA_list[[6]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",7,"round"))+ NoLegend()
hplot8<-DimPlot(scRNA_list[[7]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",8,"round"))+ NoLegend()
hplot9<-DimPlot(scRNA_list[[8]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",9,"round"))+ NoLegend()
hplot10<-DimPlot(scRNA_list[[9]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",10,"round"))+ NoLegend()
hplot11<-DimPlot(scRNA_list[[10]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",11,"round"))+ NoLegend()
hplot12<-DimPlot(scRNA_list[[11]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",12,"round"))+ NoLegend()
hplot13<-DimPlot(scRNA_list[[12]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",13,"round"))+ NoLegend()
hplot14<-DimPlot(scRNA_list[[13]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",14,"round"))+ NoLegend()
hplot15<-DimPlot(scRNA_list[[14]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle(paste0("Intergration by harmony: ",15,"round"))+ NoLegend()

hplot_merge<-CombinePlots(plots = list(hplot2,hplot3,hplot4,hplot5,hplot6,hplot7, hplot8,hplot9, hplot10,hplot11,hplot12, hplot13,hplot14, hplot15),ncol = 4,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Trophoblast_KRT7_subcyte_harmony_different_cycle.png", hplot_merge,width=20, height=20)

#final used cycle:5
#scRNAsub_final<- scRNA_list[[4]]
#DimPlot(scRNA_list[[4]], reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle("Intergration by harmony")+ NoLegend()
pc.num=1:30
scRNAsub_final<- RunHarmony(scRNAsub, group.by.vars = "sample_code",assay.use = "SCT", max.iter.harmony = 10)  
scRNAsub_final <- RunUMAP(scRNAsub_final, reduction="harmony", dims=pc.num) %>% FindNeighbors(reduction="harmony", dims=pc.num) %>% FindClusters(resolution=0.8) 
harmony_plot<-DimPlot(scRNAsub_final, reduction = "umap", group.by = "re_annotation_TFac", label = T,cols =ppCor) +ggtitle("Intergration by harmony")+ NoLegend()
harmony_plot
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Trophoblast_KRT7_subcyte_harmony_intergration.pdf", harmony_plot,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Trophoblast_KRT7_subcyte_harmony_intergration.png", harmony_plot,width=8, height=8)
saveRDS(scRNAsub_final, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_harmony_seurat_object_final.rds")
#scRNAsub_final<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_harmony_seurat_object_final.rds")
scRNAsub_final$log10_nCount_RNA<-log10(scRNAsub_final$nCount_RNA)
plot_feature<-FeaturePlot(object = scRNAsub_final, features =c("nFeature_RNA","log10_nCount_RNA","percent.mt","percent.rps"),cols= c("grey", "purple"))
plot_feature
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Trophoblast_KRT7_subcyte_harmony_four_feature_plot.png", plot_feature,width=10, height=10)

#perform monocle3 
mono_object<-selected_cell
##step2 ：Construct CellDataSet object 10X的数据使用UMI count矩阵 
count_raw <- GetAssayData(mono_object, assay = 'RNA', slot = 'counts')
cell_metadata <- mono_object@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(count_raw))
rownames(gene_annotation) <- rownames(count_raw)

#Construct monocle cds
cds <- new_cell_data_set(count_raw,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
count_raw[1:4,1:4]
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)#可至100 # 这里设置这么多主成分，可能是因为细胞数太多(4万多个)，成分太少不足以代表整体
plot_pc_variance_explained(cds)#检查一下使用的主成分（PCs）是否能够抓取最主要的基因表达变化信息 # 很像Seurat的ElbowPlot()

#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA",reduction_method = c("UMAP"))
cds_raw<-cds
saveRDS(cds_raw, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_raw_monocle3_pseudotime_object.rds")
#cds_raw<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_raw_monocle3_pseudotime_object.rds")

# 关于这个函数：多种继续降维的算法：如"UMAP", "tSNE", "PCA" and "LSI"，不过函数默认使用UMAP
plot_monocle_reduce_raw  <- plot_cells(cds_raw, reduction_method="UMAP", color_cells_by="re_annotation_TFac")+scale_colour_manual(values=ppCor)+ ggtitle('re_annotation_TFac:cds.umap')
plot_monocle_reduce_raw2 <- plot_cells(cds_raw, reduction_method="UMAP", color_cells_by="sample_code")+scale_colour_manual(values=ppCor)+ ggtitle('sample_code:cds.umap')

table(as.character(colData(cds)$"re_annotation_TFac"))
#CTBs_1 CTBs_2 EVTs_1 EVTs_2 EVTs_3 STBs_1 STBs_2 STBs_3 
#7837   8746   1848   6705   4746    372   1955    903 
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Trophoblast_KRT7_monocle_reduce_raw_re_annotation_TFac.pdf", plot_monocle_reduce_raw,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Trophoblast_KRT7_monocle_reduce_raw_sample_code.pdf", plot_monocle_reduce_raw2,width=10, height=10)

##Transfer Seurat UMAP embeddings
cds@reducedDims$UMAP <- scRNAsub_final@reductions$umap@cell.embeddings
plot_cells(cds, reduction_method="UMAP", color_cells_by="re_annotation_TFac")+scale_colour_manual(values=ppCor)
plot_cells(cds, reduction_method="UMAP", color_cells_by="Phase")+scale_colour_manual(values=ppCor)
prolifer_maker<-FeaturePlot(object = scRNAsub_final, features = c("PCNA","MKI67"),cols= c("grey", "purple"))
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Trophoblast_KRT7_monocle3_reduce_raw_prolifer_maker.pdf", prolifer_maker,width=14, height=7)

ciliated_genes <- c("PCNA","MKI67")
plot_cells(cds,genes=ciliated_genes,label_cell_groups=FALSE,show_trajectory_graph=FALSE)

## Monocle3聚类分区
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) +ggtitle("label by partitionID")
p02 <-wrap_plots(p1, p2)
p02
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Trophoblast_KRT7_monocle3_reduce_raw_cluster_find.pdf", p02,width=12, height=6)

## 识别轨迹
cds <- learn_graph(cds)
p <-plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE,label_branch_points = FALSE)
p

plot_cells(cds, reduction_method="UMAP", color_cells_by="Phase")+scale_colour_manual(values=ppCor)
plot_cells(cds, reduction_method="UMAP", color_cells_by="re_annotation_TFac")+scale_colour_manual(values=ppCor)
plot_cells(cds, reduction_method="UMAP", color_cells_by="sample_code")+scale_colour_manual(values=ppCor)
plot_cells(cds, reduction_method="UMAP", color_cells_by="sample")+scale_colour_manual(values=ppCor)

##细胞按拟时排序
cds <- order_cells(cds) 
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)

#如果存在bug，使用辅助线选择root细胞
##p + geom_vline(xintercept = seq(3,5,0.5)) + geom_hline(yintercept = seq(-8,-6,0.5))
##embed <- data.frame(Embeddings(scRNAsub, reduction = "umap"))
##range(embed$UMAP_1);range(embed$UMAP_2)
##embed2 <- subset(embed, UMAP_1 > 3 & UMAP_1 < 4 & UMAP_2 > -7.5 & UMAP_2 < -6.5)
##root.cell <- rownames(embed2)
##cds <- order_cells(cds, root_cells = root.cell)
##plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
##利用函数识别与初始细胞类群最接近的设置为根
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="CTBs_1"){
  cell_ids <- which(colData(cds)[, "re_annotation_TFac"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
saveRDS(cds, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_harmony_monocle3_pseudotime_object_1208.rds")

#saveRDS(cds, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_harmony_monocle3_pseudotime_object.rds")


###formal cds object used
#cds_used2<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_harmony_monocle3_pseudotime_object.rds")
cds_used<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_harmony_monocle3_pseudotime_object_final_used.rds")
plot_cells(cds_used, reduction_method="UMAP", color_cells_by="re_annotation_TFac")+scale_colour_manual(values=ppCor)
#colData(cds_used)<-colData(cds)
plot_cells(cds_used, reduction_method="UMAP", color_cells_by="re_annotation_TFac")+scale_colour_manual(values=ppCor)

all_anno_plot0<-plot_cells(cds_used, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
all_anno_plot1<-plot_cells(cds_used, color_cells_by = "Treat", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
all_anno_plot2<-plot_cells(cds_used, color_cells_by = "Phase", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
all_anno_plot3<-plot_cells(cds_used, color_cells_by = "re_annotation_TFac", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
all_anno_plot4<-plot_cells(cds_used, color_cells_by = "sample_code", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
all_anno_plot5<-plot_cells(cds_used, color_cells_by = "sample", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
#ggarrange(plot_j_gene_TRA_sig2, plot_j_gene_TRB_sig2,heights = c(0.7, 0.3),widths=c(1,1),ncol = 1, nrow = 2)
all_anno_plot<-CombinePlots(plots = list(all_anno_plot0,all_anno_plot1,all_anno_plot2,all_anno_plot3,all_anno_plot4,all_anno_plot5),ncol = 3,legend = NULL)
all_anno_plot_split1<-all_anno_plot0+facet_wrap(~ Treat)
all_anno_plot_split2<-all_anno_plot1+facet_wrap(~ Treat)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/all_Trophoblast_monocle3_pseudotime_split.pdf"),all_anno_plot_split1,width = 13, height =6)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/all_Trophoblast_monocle3_Treat_split.pdf"),all_anno_plot_split2,width = 13, height =6)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/all_Trophoblast_monocle3_pseudotime_six.pdf"),all_anno_plot,width = 20, height =12)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/all_Trophoblast_monocle3_pseudotime_six.png"),all_anno_plot,width = 20, height =12)


##时序基因分析
#step 1 selected cells wanna annalysis
cds_final <- choose_cells(cds_used)
p <-plot_cells(cds_final, label_groups_by_cluster = FALSE, label_leaves = FALSE,label_branch_points = FALSE)
p
final_anno_plot0<-plot_cells(cds_final, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
final_anno_plot1<-plot_cells(cds_final, color_cells_by = "Treat", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
final_anno_plot2<-plot_cells(cds_final, color_cells_by = "Phase", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
final_anno_plot3<-plot_cells(cds_final, color_cells_by = "re_annotation_TFac", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
final_anno_plot4<-plot_cells(cds_final, color_cells_by = "sample_code", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
final_anno_plot5<-plot_cells(cds_final, color_cells_by = "sample", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
#ggarrange(plot_j_gene_TRA_sig2, plot_j_gene_TRB_sig2,heights = c(0.7, 0.3),widths=c(1,1),ncol = 1, nrow = 2)
final_anno_plot<-CombinePlots(plots = list(final_anno_plot0,final_anno_plot1,final_anno_plot2,final_anno_plot3,final_anno_plot4,final_anno_plot5),ncol = 3,legend = NULL)
final_anno_plot_split1<-final_anno_plot0+facet_wrap(~ Treat)
final_anno_plot_split2<-final_anno_plot1+facet_wrap(~ Treat)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_pseudotime_split.pdf"),final_anno_plot_split1,width = 13, height =6)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_Treat_split.pdf"),final_anno_plot_split2,width = 13, height =6)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_pseudotime_six.pdf"),final_anno_plot,width = 20, height =12)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_pseudotime_six.png"),final_anno_plot,width = 20, height =12)

###reananlysis for choosen cells
cds_subset <- cds_final
cds_subset <- choose_cells(cds_second)
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
subset_anno_plot<-CombinePlots(plots = list(subset_anno_plot0,subset_anno_plot1,subset_anno_plot2,subset_anno_plot3,final_anno_plot4,subset_anno_plot5),ncol = 3,legend = NULL)
subset_anno_plot_split1<-subset_anno_plot0+facet_wrap(~ Treat)
subset_anno_plot_split2<-subset_anno_plot1+facet_wrap(~ Treat)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/subset_Trophoblast_monocle3_pseudotime_split.pdf"),subset_anno_plot_split1,width = 13, height =6)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/subset_Trophoblast_monocle3_Treat_split.pdf"),subset_anno_plot_split2,width = 13, height =6)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/subset_Trophoblast_monocle3_pseudotime_six.pdf"),subset_anno_plot,width = 20, height =12)
ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/subset_Trophoblast_monocle3_pseudotime_six.png"),subset_anno_plot,width = 20, height =12)
saveRDS(cds_subset, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/subset_Trophoblast_harmony_monocle3_pseudotime_object_reanalysis.rds")

#monocle3差异分析
##寻找拟时轨迹相关基因：轨迹上位置相似的细胞是否有相关的表达：
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
  #空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
Track_genes <- graph_test(cds_final, neighbor_graph="principal_graph", cores=30)
head(Track_genes)
write.table(Track_genes, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_gene_cor_along_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##寻找共表达模块
###no selection
genelist <- pull(Track_genes, gene_short_name) %>% as.character()##不做筛选
gene_module_unselect <- find_gene_modules(cds_final[genelist,], resolution=1e-2, cores = 30)
##selection for Track_genes
pr_deg_ids <- row.names(subset(Track_genes, q_value < 0.05 & morans_I > 0.25))#进一步筛选
length(unique(pr_deg_ids));#length(unique(genelist))### 15904 ## 24099
gene_module_df <- find_gene_modules(cds_final[pr_deg_ids,], resolution=1e-2, cores = 30)## 将轨迹可变基因收集到模块中：
head(gene_module_df);dim(gene_module_df) #1296    5
# 1   2   3   4   5   6   7   8   9  10  11  12 
#244 183 168 156 104  98  92  74  64  62  30  21
write.table(gene_module_df, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_gene_module_along_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
#基因表达趋势图
plot_genes_in_pseudotime(cds_final[Track_genes_sig,],color_cells_by = "re_annotation_TFac", min_expr=0.5, ncol = 2)
#FeaturePlot图
plot_cells(cds_final, genes=Track_genes_sig, show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)

##moldule expression in different cell types
gene_module <- gene_module_df
cell_group <- tibble::tibble(cell=row.names(colData(cds_final)), cell_group=colData(cds_final)$re_annotation_TFac)
agg_mat <- aggregate_gene_expression(cds_final, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
head(agg_mat)
pheatmap::pheatmap(agg_mat, scale="column",clustering_method="ward.D2")
module_heatmap<-pheatmap::pheatmap(agg_mat, scale="column", cluster_cols=F,clustering_method="ward.D2")

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/heatmap_for_average_expression_of_psudotime_module_ALL_in_cellsubtype_of_Trophoblast_KRT7.pdf",width =8,height = 12)
print(module_heatmap)
dev.off()
##We can also pass gene_module_df to plot_cells() as we did when we compared clusters in the L2 data above.
all_module_expression<-plot_cells(cds_final, genes=gene_module_df ,label_cell_groups=FALSE,show_trajectory_graph=FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_all_module_expression.png",all_module_expression,width = 20, height =20,limitsize = FALSE)
all_module_expression2<-all_module_expression + scale_color_gradient2(low = "blue", mid = "gray", high =  "red")
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_all_module_expression2.png",all_module_expression2,width = 20, height =20,limitsize = FALSE)


plot_cells(cds_final, genes=gene_module_df %>% filter(module %in% c(27, 10, 7, 30)),label_cell_groups=FALSE,show_trajectory_graph=FALSE)
###save object
saveRDS(cds_final,file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_final_used.rds")

pseudotime <- pseudotime(cds_final, reduction_method = 'UMAP')
mono_object2<-subset(x = mono_object, cells = names(pseudotime))
pseudotime <- pseudotime[rownames(mono_object2@meta.data)]
mono_object2$pseudotime <- pseudotime
p <-FeaturePlot(mono_object2, reduction = "umap", features = "pseudotime")+ggtitle("pseudotime for original Umap for selected Trophoblast")
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_pseudotime_in_original_Umap.pdf",p,width = 8, height =8,limitsize = FALSE)
saveRDS(mono_object2,file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_seurat_object_monocle3_pseudotime_add.rds")

scRNAsub_final2<-subset(x = scRNAsub_final, cells = names(pseudotime))
pseudotime <- pseudotime[rownames(scRNAsub_final2@meta.data)]
scRNAsub_final2$pseudotime <- pseudotime
p = FeaturePlot(scRNAsub_final2, reduction = "umap", features = "pseudotime")+ggtitle("pseudotime for rearrange Umap for selected Trophoblast")
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_pseudotime_in_rearrange_Umap.pdf",p,width = 8, height =8,limitsize = FALSE)

saveRDS(scRNAsub_final2,file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_seurat_object_rearrange_monocle3_pseudotime_add.rds")

################
#We can compare the inferred pseudotime to the known sampling timepoints.
meta_plot<-mono_object2@meta.data
meta_plot$re_annotation_TFac_Treat<-paste0(meta_plot$re_annotation_TFac,"_",meta_plot$Treat)
cell_distribution<-ggplot(meta_plot,aes(x = pseudotime,y = re_annotation_TFac, colour = re_annotation_TFac)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = ppCor) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("cell subtypes") +
  ggtitle("Cells ordered by monocle3 pseudotime")
cell_distribution1<-cell_distribution+facet_wrap(~Treat,ncol = 2)
cell_distribution1
cell_distribution2<-ggplot(meta_plot,aes(x = pseudotime,y = re_annotation_TFac_Treat, colour = Treat)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = ppCor[c(3,5)]) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("cell subtypes") +
  ggtitle("Cells ordered by monocle3 pseudotime")
cell_distribution2

dis_plot1<-ggplot(meta_plot,aes(x=pseudotime,fill=Treat,alpha = 1/10))+geom_density(position="identity")+scale_fill_manual(values = ppCor[c(3,5)]) + theme_classic()
#ggplot(meta_plot2,aes(x=Pseudotime,fill=Treat,alpha = 1/10)) +geom_density(position="stack")+scale_fill_manual(values = ppCor[c(3,5)]) + theme_classic()
dis_plot2<-ggplot(meta_plot,aes(x=pseudotime,fill=Treat,alpha = 1/10))+geom_density()+scale_fill_manual(values = ppCor[c(3,5)]) + theme_classic()+facet_wrap(~re_annotation_TFac, ncol = 3)
dis_plot3<-ggplot(meta_plot,aes(x=pseudotime,colour=Treat))+geom_density()+scale_color_manual(values = ppCor[c(3,5)]) + theme_classic() +facet_wrap(~re_annotation_TFac, ncol = 3)
#ggplot(meta_plot,aes(x=pseudotime,fill=re_annotation_TFac,alpha = 1/10,colour=Treat))+
# geom_density(position="stack")+scale_color_manual(values = ppCor[c(3,5)]) + scale_fill_manual(values = ppCor) +theme_classic()

ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_pseudotime_in_dotplot_treatsplit.pdf",cell_distribution1,width = 8, height =8,limitsize = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_pseudotime_in_dotplot.pdf",cell_distribution2,width = 6, height =8,limitsize = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_pseudotime_in_whole_density.pdf",dis_plot1,width = 6, height =6,limitsize = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_pseudotime_in_each_celltype_density.pdf",dis_plot2,width = 9, height =9,limitsize = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_pseudotime_in_each_celltype_density2.pdf",dis_plot3,width = 9, height =9,limitsize = FALSE)


#############未运行

#Subset cells by branch

##Working with 3D trajectories
cds_3d <- reduce_dimension(cds_used, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
cds_3d_plot_obj
saveRDS(cds_3d,file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_raw_monocle3_pseudotime_object_3D.rds")

marker_genes<-c("CD93","CD34","CLDN5","PECAM1")

##绘制基因的假时序表达热图 因为未具体划分分化方向,因此意义不大
##https://www.jianshu.com/p/a8c9ff0fe4a8
##这是另外的实现方法
genes <- row.names(subset(Track_genes, q_value <= 0.05 & morans_I > 0.25))
pt.matrix <- exprs(cds_final)[match(genes,rownames(rowData(cds_final))),order(pseudotime(cds_final))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
pt.matrix[1:4,1:4]
#K means with 6 groups
htkm <- Heatmap(pt.matrix,name= "z-score",
  col= colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
  km = 6,row_title_rot= 0,cluster_rows= TRUE,cluster_row_slices= FALSE,cluster_columns= FALSE)

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(pt.matrix,name= "z-score",
  col= colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
  row_title_rot= 0,cluster_rows = TRUE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
print(hthc)

##选择特定类群中特定基因的表达趋势图
##You can select a path with choose_cells() or by subsetting the cell data set by cluster, cell type, 
###or other annotation that's restricted to the path. 
##Let's pick one such path, the AFD cells:
AFD_genes <- c("gcy-8", "dac-1", "oig-8");
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,colData(cds)$cell.type %in% c("AFD")]
plot_genes_in_pseudotime(AFD_lineage_cds,color_cells_by="embryo.time.bin",min_expr=0.5)
##关于内皮部分的特定基因表达
plot_cells(cds_final, genes=c("CD93","CD34","CLDN5","PECAM1"), show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)
plot_genes_in_pseudotime(cds_final[c("CD93","CD34","CLDN5","PECAM1"),], color_cells_by="re_annotation_TFac", min_expr=0.5, ncol = 2)
plot_genes_in_pseudotime(cds_final[c("CDH5"),], color_cells_by="re_annotation_TFac", min_expr=0.5, ncol = 2)
plot_cells(cds_final, genes=c("EDNRB","PDPN","TIE1","ACKR1", "PCDH17","FAM167B","AQP1","HP","IFI27"), show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)
plot_cells(cds_final, genes=c("IFITM1"), show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)
plot_genes_in_pseudotime(cds_final[c("CDH5","IFITM1"),], color_cells_by="re_annotation_TFac", min_expr=0.5, ncol = 2)

 
##精细化分析与特定分支或者是分化路径相关的基因表达模式
###Analyzing branches in single-cell trajectories
###Analyzing the genes that are regulated around trajectory branch points often provides insights into the genetic circuits that control cell fate decisions. 
###Monocle can help you drill into a branch point that corresponds to a fate decision in your system.
#Doing so is as simple as selecting the cells (and branch point) of interest with
cds_subset <- choose_cells(cds)
#And then calling graph_test() on the subset. 
#This will identify genes with interesting patterns of expression that fall only within the region of the trajectory you selected, giving you a more refined and relevant set of genes.
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
#Grouping these genes into modules can reveal fate specific genes or those that are activate immediate prior to or following the branch point:
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)
#We will organize the modules by their similarity (using hclust) over the trajectory so it's a little easier to see which ones come on before others:
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module,levels = row.names(agg_mat)[module_dendro$order])
plot_cells(cds_subset,genes=gene_module_df,label_cell_groups=FALSE,show_trajectory_graph=FALSE)
str(cds)
rowData(cds)
colData(cds)
colnames(colData(cds))


