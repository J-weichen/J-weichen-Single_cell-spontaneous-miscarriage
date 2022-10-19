rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(Seurat)
library(future)

library(data.table)
library(scales)
library(dplyr)
library(psych)

library(ggplot2)
library(cowplot)
library(ggsci)
library(grid)
library(gridExtra)

library(pheatmap)
library(ComplexHeatmap)

grid.newpage()

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
#n=5;barplot(rep(1,n), col=colorRampPalette(colors = c('red', 'white'))( n ))
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
show_col(ppCor_all2)
##plot cells
Cells_col<-colorRampPalette(colors = rev(Cells_col_raw))(80)
length(unique(Cells_col))
barplot(rep(1,80), col=Cells_col)
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(80)
barplot(rep(1,80), col=Cells_col)
#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
plan(strategy = "multicore", workers = 42)
options(future.globals.maxSize = 1500 * 1024^12)

#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
plan(strategy = "multicore", workers = 42)
options(future.globals.maxSize = 1500 * 1024^12)

target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
target_final<-subset(x = target,subset = final_major_subgroup== "Troph_mix_group",invert = TRUE)
major_order<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN",
               "Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order<-c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2","MyCs","HCs","STCs","FBs","Endo","Epi", "NKs","Ts","Bs","Ery","Masts")
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief, levels=subgroup_order,ordered=TRUE)
target_final$final_major_group_brief<- factor(target_final$final_major_group_brief, levels=major_order,ordered=TRUE)

#plot  group
DimPlot(target_final, group.by = "final_major_subgroup_brief",cols =Cells_col)
DimPlot(target_final, group.by = "final_major_group_brief",label = F,cols = ppCor_all2)
DimPlot(object = target_final, group.by = "final_major_group_brief", split.by = "Treat",reduction = "umap",cols = ppCor_all2,ncol = 2)

subgroup_order<-c("CTBs_1","CTBs_0","CTBs_3","CTBs_2","EVTs_2","EVTs_3","EVTs_1","STBs_1","STBs_2",
                  "Epi","MyCs","HCs","Endo","STCs","FBs","NKs","Ts","Bs","Masts","Ery")
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief,levels = rev(subgroup_order))

cell_name<-"All_final_major_subgroup"
Allcluster.markers_pos_005<-read.table(file = paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_one_result/inter_celltype_DEGs/",cell_name,"_default0.25_pos_001_genes.xls",sep=""),header =T,sep="\t")
head(Allcluster.markers_pos_005)
table(Allcluster.markers_pos_005$cluster)
Allcluster.markers_pos_005$cluster<-factor(Allcluster.markers_pos_005$cluster,levels = rev(subgroup_order))

#identify the relationship among different cell types
Allcluster.markers_pos_005<-Allcluster.markers_pos_005[order(Allcluster.markers_pos_005$cluster,decreasing = T),]
plot_maker<-Allcluster.markers_pos_005 %>% group_by(cluster)
top50 <- plot_maker %>% top_n(n = 50, wt = avg_logFC)
top10 <- plot_maker %>% top_n(n = 10, wt = avg_logFC)
top2 <- plot_maker %>% top_n(n = 2, wt = avg_logFC)
top1 <- plot_maker %>% top_n(n = 1, wt = avg_logFC)

target_final
DefaultAssay(target_final) <- "RNA"
target_final <- NormalizeData(target_final, verbose = FALSE)
head(target_final@meta.data)
table(target_final$final_major_subgroup_brief)
names(table(target_final$final_major_subgroup))

#有5个可视化方法，分别是：坐标映射图，峰峦图,小提琴图，，气泡图，热图。
#后三类也可参考https://www.cnblogs.com/TOP-Bio/p/14471750.html
plot_umap_maker <- FeaturePlot(object = target_final, features = as.character(top1$gene),cols= c("grey", "red"),ncol=5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Umap_final_submaker_top_one.png", plot_umap_maker,width=25, height=20)

plot_Ridge_maker <-RidgePlot(target_final,group.by = "final_major_subgroup_brief", features = as.character(top1$gene), ncol = 5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Ridge_final_submaker_top_one.png", plot_Ridge_maker,width=25, height=20)

#2.经典气泡图
Idents(object = target_final) <- "final_major_subgroup_brief"
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief,levels = rev(subgroup_order))
maker_gene<-unique(as.character(top2$gene))
Plotdot1 <-DotPlot(target_final, features = maker_gene,group.by ="final_major_subgroup_brief",cols = c("blue","red")) + RotatedAxis()
Plotdot2 <-DotPlot(target_final, features = maker_gene,group.by ="final_major_subgroup_brief",split.by = "Treat",cols =  c("blue","red"),assay = "RNA") + RotatedAxis()
Plotdot3 <-DotPlot(target_final, features = maker_gene,split.by = "Treat",cols = ppCor,assay = "RNA") + RotatedAxis()
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Plotdot_final_top2_submaker1.pdf", Plotdot1,width=10, height=8)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Plotdot_final_top2_submaker2.pdf", Plotdot2,width=10, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Plotdot_final_top2_submaker3.pdf", Plotdot3,width=10, height=15)

##3. heatmap for averageExpression
maker_gene<-unique(as.character(top2$gene))
AverageExp<-AverageExpression(target_final,assays="RNA",slot = "data",features=maker_gene)
typeof(AverageExp);head(AverageExp$RNA)
range(AverageExp$RNA)#0.000 1528.887
AverageExp_data<-log2(AverageExp$RNA+1)
range(AverageExp_data);head(AverageExp_data)
#AverageExp_data[AverageExp_data == "Inf"]<- 2000
colnames(AverageExp_data)
AverageExp_data<-AverageExp_data[,subgroup_order]

p1<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
             show_rownames = T,show_colnames = T,border=TRUE,border_color ="grey", 
             scale ="row", color = colorRampPalette(c("navy","white","firebrick3"))(60),
             main ="log2(mean of expression level of maker + 1):scale by row")
pdf("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_average_heatmap1.pdf",width =8,height = 8)
print(p1)
dev.off()

#similarity of all cell types
maker_gene<-unique(as.character(plot_maker$gene))
AverageExp<-AverageExpression(target_final,assays="RNA",slot = "data",features=maker_gene)
typeof(AverageExp);
head(AverageExp$RNA)
range(AverageExp$RNA)#0.000 1528.887
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap(coorda$r)

##4.整体细胞类群表达水平绘制
target_1K<-subset(target_final, downsample = 1000)
target_1K@assays$RNA@scale.data <- scale(target_1K@assays$RNA@data, scale = TRUE)
#target_1K <- ScaleData(object = target_1K, features = rownames(target_1K))
PH<-DoHeatmap(object = target_1K, features =maker_gene,draw.lines = TRUE,label = F,
              disp.min = -2.5, disp.max = 2.5,assay = "RNA",slot = "scale.data", size = 3, angle = -50, hjust=0.8,
              group.by ="final_major_subgroup_brief",group.bar.height = 0.1,group.colors=ppCor[1:length(levels(target_1K$final_major_subgroup_brief))]) 
PH2<-PH+scale_fill_gradientn(colors = c("lightblue", "white", "red"))+theme(axis.text.y = element_text(size = 10))
PH2
ggsave(PH2,file="/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/intersubgroup_genes_all_heatmap1.png",width = 15,height = 10)
ggsave(PH,file="/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/intersubgroup_genes_all_heatmap2.png",width = 15,height = 10)

#手动提取
levels(target_final)
DefaultAssay(target_final)<-"SCT"
#mat0 <- GetAssayData(target, slot = "counts")
#mat0 <- log2(mat0 + 1)
#其他SCT的数据说明：https://www.jianshu.com/p/e639cc257d51 
mat0<- GetAssayData(target_final, slot = "data")

#FC_10<-Allcluster.markers[Allcluster.markers$p_val_adj<0.01 & Allcluster.markers$avg_logFC>1,]%>% group_by(cluster)
#FC_1.7<-Allcluster.markers[Allcluster.markers$p_val_adj<0.01 & Allcluster.markers$avg_logFC>0.25,]%>% group_by(cluster)
#top10 <- FC_1.7 %>% top_n(n = 10, wt = avg_logFC)

#get genes and cluster information
#gene_features <-top10$gene
FC_1.5<-plot_maker[plot_maker$p_val_adj<0.01 & plot_maker$avg_logFC>0.5849625,]%>% group_by(cluster)
top50_FC_1.5 <- FC_1.5 %>% top_n(n = 50, wt = avg_logFC)
table(top50_FC_1.5$cluster)
names(table(top50_FC_1.5$cluster))
gene_features <- unique(top50_FC_1.5$gene);length(gene_features)
cluster_info <- target_final$final_major_subgroup_brief

#sort the expression matrix 
mat1 <- as.matrix(mat0[as.character(gene_features), names(cluster_info)])
mat1[1:6,1:6];dim(mat1)
name_mean1<-as.data.frame(t(unlist(apply(mat1,1,mean))))
name_mean1<-as.data.frame(rbind(c(1:length(rownames(mat1))),name_mean1))
name_mean1<-as.data.frame(rbind(colnames(name_mean1),name_mean1))
rownames(name_mean1)<-c("Gname","row","mean")
dim(name_mean1);str(name_mean1)
name_mean1<-as.data.frame(t(name_mean1))
name_mean1$row<-as.numeric(as.character(name_mean1$row))
name_mean1$mean<-as.numeric(as.character(name_mean1$mean))
name_mean1$Gname<-as.character(name_mean1$Gname)
dim(name_mean1)
head(name_mean1)
#去重
name_mean1<-name_mean1[order(name_mean1$Gname,name_mean1$mean,decreasing = T),]
name_mean1<-name_mean1[!duplicated(name_mean1$Gname),]
name_mean1<-name_mean1[order(name_mean1$row),]
length(unique(gene_features))#2105
#name_mean1$row
head(name_mean1)
gene_features_uniq<-name_mean1$Gname
gene_features_dup<-gene_features[which(duplicated(gene_features))]
length(gene_features);length(gene_features_uniq);length(gene_features_dup)
#重新去unique数据框
mat2 <- mat1[gene_features_uniq,]
#手动scale
mat3<-scale(t(mat2), center=T,scale=T)
range(mat3)#-4.798357 162.184118
mat4<-t(mat3)
mat4[mat4>=2.5]= 2.5
mat4[mat4 < (-2.5)]= -2.5 #小于负数时，加括号！
#set color for cell types
col <- ppCor_all[1:length(levels(cluster_info))]
#colourCount<-length(unique(target$sub_cluster))
names(col) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col),labels = levels(cluster_info), labels_gp = gpar(cex = 0.5, col = "black")))

#important genes notes
mark_gene0 <- c("VIM","HLA-B","CDH1","EGFR","HLA-G","MMP2","CSH1","CGA")
#mark_gene0<-as.character(top2$gene)
gene_pos <- which(rownames(mat4) %in% mark_gene0)
mark_gene1<-rownames(mat4)[gene_pos]
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = mark_gene1))

#legend adjust
col_fun1  <- circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red"))
col_fun2  <- circlize::colorRamp2(c(-2,0,2), c("purple", "black", "yellow"))
col_fun<-col_fun1
#show_col(col_fun(seq(-3, 3)))

plot_HP<-Heatmap(mat4,cluster_rows = FALSE,cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 #show_row_names = TRUE,
                 column_split = cluster_info,
                 top_annotation = top_anno,
                 #  right_annotation = row_anno,
                 column_title = NULL,
                 col = col_fun,
                 heatmap_legend_param = list(legend_height = unit(4, "cm"),grid_width = unit(0.4, "cm"),title = "log2(count+1)",title_position = "leftcenter-rot",at = c(-2, 0, 2), labels = c("low", "median", "high") ))
pdf("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/submaker_top50_FC_1.5_heatmap1.pdf",width = 15,height = 8)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("major_submaker_heatmap_SCT_data"))
dev.off()
