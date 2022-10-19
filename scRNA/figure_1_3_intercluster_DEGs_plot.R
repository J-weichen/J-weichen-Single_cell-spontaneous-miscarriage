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
library(viridis)
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

#########################
###reading data of  Integrate for Loxy_2 and Loxy_1
########################
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

DefaultAssay(target_final) <- "RNA"
# Normalize RNA data for visualization purposes
target_final <- NormalizeData(target_final, verbose = FALSE)

##################
###plot intercluster DEGs
##################
cell_name<-"All_final_subgroup"
Allcluster.markers_pos_001<-read.table(file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_default0.25_pos_001_genes.xls",sep=""),header =T,sep="\t")
head(Allcluster.markers_pos_001)
table(Allcluster.markers_pos_001$cluster)
#   Bs CTBs_1 CTBs_2   Endo    Epi    Ery EVTs_1 EVTs_2 EVTs_3    FBs    HCs  Masts   Mycs    NKs STBs_1 STBs_2 STBs_3   STCs     Ts 
#   276    861    598   1680    867    258    613    753    874   1157   1009    894    968   1412   1088    559    745    915    656 

cell_order<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
Allcluster.markers_pos_001$cluster<-factor(Allcluster.markers_pos_001$cluster,levels = rev(cell_order))

#identify the relationship among different cell types
Allcluster.markers_pos_001<-Allcluster.markers_pos_001[order(Allcluster.markers_pos_001$cluster,decreasing = T),]
plot_maker<-Allcluster.markers_pos_001 %>% group_by(cluster)
top50 <- plot_maker %>% top_n(n = 50, wt = avg_logFC)
top10 <- plot_maker %>% top_n(n = 10, wt = avg_logFC)
top5 <- plot_maker %>% top_n(n = 5, wt = avg_logFC)
top2 <- plot_maker %>% top_n(n = 2, wt = avg_logFC)
top1 <- plot_maker %>% top_n(n = 1, wt = avg_logFC)

target_final
#DefaultAssay(target_final) <- "RNA"
#target_final <- NormalizeData(target_final, verbose = FALSE)
head(target_final@meta.data)
table(target_final$re_annotation_TFac)
#CTBs_1 CTBs_2 STBs_1 STBs_2 STBs_3 EVTs_1 EVTs_2 EVTs_3    Epi   Endo    HCs   Mycs   STCs    FBs    NKs     Ts     Bs  Masts    Ery 
#  7837   8746    372   1955    903   1848   6705   4746     61    133   1902   4574   1020   2859    625    420     37    113    944 
names(table(target_final$re_annotation_TFac))

#有5个可视化方法，分别是：坐标映射图，峰峦图,小提琴图，，气泡图，热图。
#后三类也可参考https://www.cnblogs.com/TOP-Bio/p/14471750.html
plot_umap_maker <- FeaturePlot(object = target_final, features = as.character(top1$gene),cols= c("grey", "red"),ncol=5)
plot_umap_maker
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/final_intercluster_submaker_top_one.png", plot_umap_maker,width=20, height=16)

plot_Ridge_maker <-RidgePlot(target_final,group.by = "re_annotation_TFac", features = as.character(top1$gene),cols = ppCor_all,ncol = 5)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/final_Ridge_submaker_top_one.pdf", plot_Ridge_maker,width=20, height=16)

#2.经典气泡图
Idents(object = target_final) <- "re_annotation_TFac"
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels = rev(cell_order))
maker_gene<-unique(as.character(top2$gene))
Plotdot1 <-DotPlot(target_final, features = maker_gene,group.by ="re_annotation_TFac",cols = c("blue","red"),assay = "RNA") + RotatedAxis()
Plotdot2 <-DotPlot(target_final, features = maker_gene,group.by ="re_annotation_TFac",split.by = "Treat",cols =  c("blue","red"),assay = "RNA") + RotatedAxis()
Plotdot3 <-DotPlot(target_final, features = maker_gene,split.by = "Treat",cols = ppCor,assay = "RNA") + RotatedAxis()
Plotdot1;Plotdot2;Plotdot3
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker1.pdf", Plotdot1,width=12, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker2.pdf", Plotdot2,width=12, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker3.pdf", Plotdot3,width=12, height=8)

Plotdot12<-Plotdot1+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ coord_flip()+
  theme(plot.title = element_text(size=12,colour = "black",face = "bold"),
        axis.title.x = element_text(size=12,colour = "black",face = "bold"),
        axis.title.y = element_text(size=12,colour = "black",face = "bold"),
        axis.text.x = element_text(size=12,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=12,colour = "black"))
Plotdot13<-Plotdot1+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2) +
  scale_colour_viridis(option="magma") +  coord_flip()+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(plot.title = element_text(size=12,colour = "black",face = "bold"),
        axis.title.x = element_text(size=12,colour = "black",face = "bold"),
        axis.title.y = element_text(size=12,colour = "black",face = "bold"),
        axis.text.x = element_text(size=12,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=12,colour = "black"))
Plotdot02<-Plotdot1+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ 
  theme(axis.text.x = element_text(size=12,angle=90,hjust=1, vjust=0.5))
Plotdot12;Plotdot13;Plotdot02
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker4.pdf", Plotdot12,width=6, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker5.pdf", Plotdot13,width=6, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker6.pdf", Plotdot02,width=8, height=6)

target_final$celltype_sample<- paste(target_final$re_annotation_TFac,target_final$Treat,sep="-")
table(target_final$celltype_sample)
names(table(target_final$celltype_sample))
cell_order2<-c("CTBs_1-Abortion","CTBs_1-CTRL","CTBs_2-Abortion","CTBs_2-CTRL",
               "STBs_1-Abortion","STBs_1-CTRL","STBs_2-Abortion","STBs_2-CTRL","STBs_3-Abortion","STBs_3-CTRL",
               "EVTs_1-Abortion","EVTs_1-CTRL","EVTs_2-Abortion","EVTs_2-CTRL","EVTs_3-Abortion","EVTs_3-CTRL",
               "Epi-Abortion","Epi-CTRL","Endo-Abortion","Endo-CTRL",
               "HCs-Abortion","HCs-CTRL","Mycs-Abortion","Mycs-CTRL",
               "STCs-Abortion","STCs-CTRL","FBs-Abortion","FBs-CTRL","NKs-Abortion","NKs-CTRL",
               "Ts-Abortion","Ts-CTRL","Bs-Abortion","Bs-CTRL",
               "Masts-Abortion","Masts-CTRL","Ery-Abortion","Ery-CTRL")
target_final$celltype_sample<- factor(target_final$celltype_sample,levels = rev(cell_order2))
Idents(object = target_final) <- "celltype_sample"
maker_gene<-unique(as.character(top2$gene))
Plotdot4 <-DotPlot(target_final, features = maker_gene,group.by ="celltype_sample",cols = c("blue","red"),assay = "RNA") + RotatedAxis()
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker_sample_split_2.pdf", Plotdot4,width=10, height=10)

Plotdot42<-Plotdot4+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ coord_flip()+
  theme(plot.title = element_text(size=12,colour = "black",face = "bold"),
        axis.title.x = element_text(size=12,colour = "black",face = "bold"),
        axis.title.y = element_text(size=12,colour = "black",face = "bold"),
        axis.text.x = element_text(size=12,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=12,colour = "black"))
Plotdot43<-Plotdot4+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2) +
  scale_colour_viridis(option="magma") +  coord_flip()+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(plot.title = element_text(size=12,colour = "black",face = "bold"),
        axis.title.x = element_text(size=12,colour = "black",face = "bold"),
        axis.title.y = element_text(size=12,colour = "black",face = "bold"),
        axis.text.x = element_text(size=12,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=12,colour = "black"))
Plotdot04<-Plotdot4+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ 
  theme(axis.text.x = element_text(size=12,angle=90,hjust=1, vjust=0.5))
Plotdot42;Plotdot43;Plotdot04

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker_sample_split_3.pdf", Plotdot42,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker_sample_split_4.pdf", Plotdot43,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Plotdot_final_top2_submaker_sample_split_5.pdf", Plotdot04,width=10, height=10)

##3. heatmap for averageExpression
Idents(object = target_final) <- "re_annotation_TFac"
maker_gene<-unique(as.character(top2$gene))
AverageExp<-AverageExpression(target_final,assays="RNA",slot = "data",features=maker_gene)
typeof(AverageExp);head(AverageExp$RNA)
range(AverageExp$RNA)#0.000 2537.101
AverageExp_data<-log2(AverageExp$RNA+1)
range(AverageExp_data);head(AverageExp_data)
#AverageExp_data[AverageExp_data == "Inf"]<- 2000
colnames(AverageExp_data)
AverageExp_data<-AverageExp_data[,cell_order]

p1<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
             show_rownames = T,show_colnames = T,border=TRUE,border_color ="grey", 
             scale ="row", color = colorRampPalette(c("navy","white","firebrick3"))(60),
             main ="log2(mean of expression level of maker + 1):scale by row")
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/TOP2_submaker_average_heatmap_sample_no_split.pdf",width =8,height = 8)
print(p1)
dev.off()

##for sample split
Idents(object = target_final) <- "celltype_sample"
maker_gene<-unique(as.character(top2$gene))
AverageExp<-AverageExpression(target_final,assays="RNA",slot = "data",features=maker_gene)
typeof(AverageExp);head(AverageExp$RNA)
range(AverageExp$RNA)#0.00 2770.25
AverageExp_data<-log2(AverageExp$RNA+1)
range(AverageExp_data);head(AverageExp_data)# 0.00000 11.43632
#AverageExp_data[AverageExp_data == "Inf"]<- 2000
colnames(AverageExp_data)
AverageExp_data<-AverageExp_data[,cell_order2]

p2<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
             show_rownames = T,show_colnames = T,border=TRUE,border_color ="grey", 
             scale ="row", color = colorRampPalette(c("navy","white","firebrick3"))(60),
             main ="log2(mean of expression level of maker + 1):scale by row")
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/TOP2_submaker_average_heatmap_sample_split.pdf",width =8,height = 8)
print(p2)
dev.off()

#based on intercluster DEGs for similarity of all cell types
Idents(object = target_final) <- "re_annotation_TFac"
maker_gene<-unique(as.character(plot_maker$gene))
AverageExp<-AverageExpression(target_final,assays="RNA",slot = "data",features=maker_gene)
typeof(AverageExp);
head(AverageExp$RNA)
range(AverageExp$RNA)#0.000 2537.101
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
p3<-pheatmap(coorda$r)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/Intercluster_DEGs_based_sample_spearman_similarity_heatmap.pdf",width =6,height = 6)
print(p3)
dev.off()

##4.整体细胞类群表达水平绘制
Idents(object = target_final) <- "re_annotation_TFac"
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels = cell_order)

levels(target_final)
DefaultAssay(target_final)<-"RNA"
#DefaultAssay(target_final)<-"SCT"
#mat0 <- GetAssayData(target, slot = "counts")
#mat0 <- log2(mat0 + 1)
#其他SCT的数据说明：https://www.jianshu.com/p/e639cc257d51 
mat0<- GetAssayData(target_final, slot = "data")
conflict_prefer("which", "Matrix");conflict_prefer("rbind", "spam")
#FC_10<-Allcluster.markers[Allcluster.markers$p_val_adj<0.01 & Allcluster.markers$avg_logFC>1,]%>% group_by(cluster)
#FC_1.7<-Allcluster.markers[Allcluster.markers$p_val_adj<0.01 & Allcluster.markers$avg_logFC>0.25,]%>% group_by(cluster)
#top10 <- FC_1.7 %>% top_n(n = 10, wt = avg_logFC)

#get genes and cluster information
#gene_features <-top10$gene
head(plot_maker)
FC_1.5<-plot_maker[plot_maker$p_val_adj<0.01 & plot_maker$avg_logFC>0.5849625,]%>% group_by(cluster)
top50_FC_1.5 <- FC_1.5 %>% top_n(n = 50, wt = avg_logFC)
table(top50_FC_1.5$cluster)
names(table(top50_FC_1.5$cluster))
gene_features <- unique(top50_FC_1.5$gene);length(gene_features)

gene_features <- unique(plot_maker$gene);length(gene_features)
cluster_info <- target_final$re_annotation_TFac

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
range(mat3)#  -8.043521 143.459217
mat4<-t(mat3)
mat4[mat4>=2.5]= 2.5
mat4[mat4 < (-2.5)]= -2.5 #小于负数时，加括号！

dim(mat4)#7258 45800
#set color for cell types
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
####col <- ppCor_all[1:length(levels(cluster_info))]
col <-my_morandi_colors[1:19]
#colourCount<-length(unique(target$sub_cluster))
names(col) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col),labels = levels(cluster_info), labels_gp = gpar(cex = 0.5, col = "black")))

#important genes notes
#mark_gene0 <- c("POU5F1","NANOG","SOX2", "KLF17", "GATA6","PDGFRA","GATA4","FN1", "GATA3","GATA2","KRT7","TFAP2C","KRT8","KRT18","PDGFA","LRP2")
mark_gene0<-c("PCNA","PAGE4","CDH1","ERVFRD-1","CYP19A1","EGFR","MKI67","HLA-G","PAPPA2",
              "PAEP","PECAM1","LYVE1","AIF1","DCN","DLK1","NCAM1","CD3D","CD79A","MS4A2","PF4","HBA1")

#mark_gene0<-as.character(top2$gene)
gene_pos <- which(rownames(mat4) %in% mark_gene0)
gene_pos
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
                 right_annotation = row_anno,
                 column_title = NULL,
                 col = col_fun,
                 heatmap_legend_param = list(legend_height = unit(4, "cm"),grid_width = unit(0.4, "cm"),
                                             title = "log2(count+1)",title_position = "leftcenter-rot",at = c(-2, 0, 2),
                                             labels = c("low", "median", "high") ))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/final_all_19_intercuster_all_DEGs_heatmap.pdf",width = 10,height =8)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("major_submaker_heatmap_RNA_data"))
dev.off()

png("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/final_all_19_intercuster_all_DEGs_heatmap.png",height=1000,width=1000)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("major_submaker_heatmap_RNA_data"))
dev.off()

###for top 50
gene_features <- unique(top50$gene);length(gene_features)
cluster_info <- target_final$re_annotation_TFac
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
range(mat3)#-5.304552 143.459217
mat4<-t(mat3)
mat4[mat4>=2.5]= 2.5
mat4[mat4 < (-2.5)]= -2.5 #小于负数时，加括号！
#set color for cell types
####col <- ppCor_all[1:length(levels(cluster_info))]
col <-my_morandi_colors[1:19]
#colourCount<-length(unique(target$sub_cluster))
names(col) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col),labels = levels(cluster_info), labels_gp = gpar(cex = 0.5, col = "black")))

#important genes notes
#mark_gene0 <- c("POU5F1","NANOG","SOX2", "KLF17", "GATA6","PDGFRA","GATA4","FN1", "GATA3","GATA2","KRT7","TFAP2C","KRT8","KRT18","PDGFA","LRP2")
mark_gene0<-c("PCNA","PAGE4","CDH1","ERVFRD-1","CYP19A1","EGFR","MKI67","HLA-G","PAPPA2",
              "PAEP","PECAM1","LYVE1","AIF1","DCN","DLK1","NCAM1","CD3D","CD79A","MS4A2","PF4","HBA1")

#mark_gene0<-as.character(top2$gene)
gene_pos <- which(rownames(mat4) %in% mark_gene0)
gene_pos
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
                 right_annotation = row_anno,
                 column_title = NULL,
                 col = col_fun,
                 heatmap_legend_param = list(legend_height = unit(4, "cm"),grid_width = unit(0.4, "cm"),
                                             title = "log2(count+1)",title_position = "leftcenter-rot",at = c(-2, 0, 2),
                                             labels = c("low", "median", "high") ))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/final_all_19_intercuster_top50_DEGs_heatmap.pdf",width = 8,height =8)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = "Four_cluster_top50_heatmap_RNA_data")
dev.off()

png("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/final_all_19_intercuster_top50_DEGs_heatmap.png",width =1000,height =1000)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = "Four_cluster_top50_heatmap_RNA_data")
dev.off()
###以下未跑
##more strict limitation in gene select
#get genes and cluster information
#gene_features <-top10$gene
head(plot_maker)
FC_1.5<-plot_maker[plot_maker$p_val_adj<0.01 & plot_maker$avg_logFC>0.5849625,]%>% group_by(cluster)
top50_FC_1.5 <- FC_1.5 %>% top_n(n = 50, wt = avg_logFC)
table(top50_FC_1.5$cluster)
names(table(top50_FC_1.5$cluster))
gene_features <- unique(top50_FC_1.5$gene);length(gene_features)

gene_features <- unique(plot_maker$gene);length(gene_features)
cluster_info <- target_final$celltype

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
####col <- ppCor_all[1:length(levels(cluster_info))]
col <-c(ppCor_all[c(7,37,40)],"grey")
#colourCount<-length(unique(target$sub_cluster))
names(col) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col),labels = levels(cluster_info), labels_gp = gpar(cex = 0.5, col = "black")))

#important genes notes
mark_gene0 <- c("POU5F1","NANOG","SOX2", "KLF17", "GATA6","PDGFRA","GATA4","FN1", "GATA3","GATA2","KRT7","TFAP2C","KRT8","KRT18","PDGFA","LRP2")
#mark_gene0<-as.character(top2$gene)
gene_pos <- which(rownames(mat4) %in% mark_gene0)
gene_pos
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
                 right_annotation = row_anno,
                 column_title = NULL,
                 col = col_fun,
                 heatmap_legend_param = list(legend_height = unit(4, "cm"),grid_width = unit(0.4, "cm"),
                                             title = "log2(count+1)",title_position = "leftcenter-rot",at = c(-2, 0, 2),
                                             labels = c("low", "median", "high") ))

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/final_four_intercuster_all_DEGs_heatmap.pdf",width = 8,height = 6)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("major_submaker_heatmap_RNA_data"))
dev.off()
