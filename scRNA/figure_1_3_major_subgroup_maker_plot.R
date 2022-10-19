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

#DefaultAssay(object = vitDcca) <- "RNA"
#vitDcca <- NormalizeData(object =vitDcca, normalization.method = "LogNormalize", scale.factor = 10000)
#all.genes <- rownames(x = vitDcca)
#vitDcca <- ScaleData(object = vitDcca, features = all.genes)
target_final
DefaultAssay(target_final) <- "RNA"
target_final <- NormalizeData(target_final, verbose = FALSE)

head(target_final@meta.data)
table(target_final$final_major_subgroup_brief)
names(table(target_final$final_major_subgroup))

subgroup_order<-c("CTBs_1","CTBs_0","CTBs_3","CTBs_2","EVTs_2","EVTs_3","EVTs_1","STBs_1","STBs_2",
                  "Epi","MyCs","HCs","Endo","STCs","FBs","NKs","Ts","Bs","Masts","Ery")
#major_submaker<-c("VIM","HLA-B","KRT7","PERP","MKI67","PCNA","HLA-G","PAPPA2","PAGE4","CDH1","EGFR","CGA","CSH2","ERVFRD-1",
#               "AIF1","LYVE1","DCN","DLK1","PECAM1","EPCAM","PAEP","CD3D","IL7R","CD8A","NCAM1","NKG7","CD79A","HBA1","MS4A2")
major_submaker<-c("KRT7","PAGE4","PCNA","MKI67","EGFR","HLA-G","PAPPA2","CYP19A1","CSH2","ERVFRD-1",
                  "VIM","PAEP","AIF1","LYVE1","PECAM1","DCN","DLK1","NCAM1","CD3D","CD79A","MS4A2","HBA1")
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief,levels = rev(subgroup_order))

#有5个可视化方法，分别是：坐标映射图，峰峦图,小提琴图，，气泡图，热图。
#后三类也可参考https://www.cnblogs.com/TOP-Bio/p/14471750.html
plot_umap_maker <- FeaturePlot(object = target_final, features = major_submaker,cols= c("grey", "red"),ncol=5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Umap_final_major_submaker.png", plot_umap_maker,width=30, height=30)

plot_Ridge_maker <-RidgePlot(target_final,group.by = "final_major_subgroup_brief", features = major_submaker, ncol = 5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Ridge_final_major_submaker.png", plot_Ridge_maker,width=30, height=30)

#maker plot
#小提琴图
#ref:https://www.rdocumentation.org/packages/Seurat/versions/3.1.0/topics/FetchData
#CR maker
maker_Vln_plots <- VlnPlot(target_final, features = rev(major_submaker),group.by = "final_major_subgroup_brief",cols = ppCor_all2,pt.size = 0, combine = FALSE)
maker_Vln_plots2<-CombinePlots(plots = maker_Vln_plots, ncol =1)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Violin_final_major_submaker.png", maker_Vln_plots2,width=10, height=45)

exprs <- data.frame(FetchData(object = target_final,slot = "data",vars = c(major_submaker,"Treat","final_major_subgroup_brief")))
data_plot <- melt(exprs)
data_plot$final_major_subgroup_brief<- factor(data_plot$final_major_subgroup_brief,levels = rev(subgroup_order))
maker_plot_line<-ggplot(data_plot,aes(x = final_major_subgroup_brief,y=value,fill = final_major_subgroup_brief)) + 
  geom_violin(scale = "width") +  #geom_jitter(size = 0.00)+
  coord_flip()+ facet_grid(. ~ variable,scale='free') +
  theme_bw()+scale_fill_manual(values=ppCor_all2)+
  labs(y = "Expression level", x = "main_Cell_type", title = "main_cluster marker")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x.bottom = element_blank(),axis.ticks.y.left = element_blank())+
  theme(axis.title = element_text(size=10),axis.text.x = element_text(face="italic",size=4), axis.text.y = element_text(face="italic",size=6),legend.position = "right")+
  theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))

maker_plot_up<-ggplot(data_plot,aes(x = final_major_subgroup_brief,y=value,fill = final_major_subgroup_brief)) + 
  geom_violin(scale = "width") + #geom_jitter(size = 0.00)+
  facet_grid(variable ~ .,scale='free') +  theme_bw()+scale_fill_manual(values=ppCor_all2)+
  labs(x = "Expression level", y = "main_Cell_type", title = "main_cluster marker")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y.left = element_blank(),axis.ticks.x.bottom = element_blank())+
  theme(axis.title = element_text(size=15),axis.text.x = element_text(face="italic",size=15,angle = 45,vjust =0.5,hjust=1), axis.text.y = element_text(face="italic",size=15),legend.position = "right")+
  theme(strip.text=element_text(colour = "red", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))
maker_plot_up
maker_plot_line
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Violin_final_major_submaker1.png", maker_plot_up,width=12, height=18)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Violin_final_major_submaker2.png", maker_plot_line,width=20, height=12)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Violin_final_major_submaker1.pdf", maker_plot_up,width=12, height=18)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Violin_final_major_submaker2.pdf", maker_plot_line,width=18, height=10)

#有人认为应该按照seurat中一样的参数，使用scale ="width",同时添加随机误差，但是Y叔有一篇文章对使用scale = "width"持反对态度,本人也持反对态度
#ref1:https://www.jianshu.com/p/be67a9511bef
#ref2:https://mp.weixin.qq.com/s/tSOR68y7F0CCkAV6pUoV8Q
#Aver_plot<-top10
#2.经典气泡图
Idents(object = target_final) <- "final_major_subgroup_brief"
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief,levels = rev(subgroup_order))

Plotdot1 <-DotPlot(target_final, features = major_submaker,group.by ="final_major_subgroup_brief",cols = c("blue","red")) + RotatedAxis()
Plotdot2 <-DotPlot(target_final, features = major_submaker,group.by ="final_major_subgroup_brief",split.by = "Treat",cols =  c("blue","red"),assay = "RNA") + RotatedAxis()
Plotdot3 <-DotPlot(target_final, features = major_submaker,split.by = "Treat",cols = ppCor,assay = "RNA") + RotatedAxis()
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Plotdot_final_major_submaker1.pdf", Plotdot1,width=10, height=8)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Plotdot_final_major_submaker2.pdf", Plotdot2,width=10, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Plotdot_final_major_submaker3.pdf", Plotdot3,width=10, height=15)

##3. heatmap for averageExpression
Aver_plot<-major_submaker
AverageExp<-AverageExpression(target_final,assays="RNA",slot = "data",features=major_submaker)
typeof(AverageExp);
head(AverageExp$RNA)
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

p2<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
         show_rownames = T,show_colnames = T,border=TRUE,border_color ="grey", 
         scale ="row", color = colorRampPalette(c("navy","blue","white","orange","firebrick3"))(60),
         main ="log2(mean of expression level of maker + 1):scale by row")
pdf("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_average_heatmap1.pdf",width =8,height = 8)
print(p1)
dev.off()
pdf("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_average_heatmap2.pdf",width =8,height = 8)
print(p2)
dev.off()
dev.new();dev.new()

#similarity of all cell types
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap(coorda$r)

##4.整体细胞类群表达水平绘制
names(target_final@assays)
DefaultAssay(target_final)<-"SCT"
scaled_data <- target_final@assays$SCT@scale.data;dim(scaled_data)# 4947 45800
makergenes<-major_submaker
gene_panel_intersect <- intersect(makergenes, rownames(scaled_data))
#NOTE:  downsampling is crucial to get the heatmap plotted.
#number of cell must less than 30K
#NOTE2: default 3000 genes in scale.data is calculated automatically during SCTransform where
# So, some of the genes in my panel may be absent from these 3000.
#NOTE3:if DoHeatmap found No requested Features found in the scale.data slot for the RNA assay.
#then  DefaultAssay(x)<-"RNA"  When you normalize data by SCTransform, the data and scale.data slots are empty in the RNA assay. 
#You need to change the default assay to RNA and run NormalizeData and ScaleData to generate data for those two slots.
##target@assays$RNA@scale.data <- scale(target@assays$RNA@data, scale = TRUE)
##target <- ScaleData(object = target, features = rownames(target))
##DoHeatmap(target, features = makergenes,size = 0.5,slot = "scale.data") + NoLegend()

target_1K<-subset(target_final, downsample = 1000)
scaled_data_1K <- target_1K@assays$SCT@scale.data
range(scaled_data_1K)

target_1K$final_major_subgroup_brief <- factor(x =target_1K$final_major_subgroup_brief, levels = subgroup_order)
#heatmap for all cells
PH<-DoHeatmap(object = target_1K, features =gene_panel_intersect,draw.lines = TRUE,label = F,
              disp.min = -2.5, disp.max = 2.5,assay = "SCT",size = 3, angle = -50, hjust=0.8,
              group.by ="final_major_subgroup_brief",group.bar.height = 0.1,group.colors=ppCor[1:length(levels(target_1K$final_major_subgroup_brief))]) 
PH2<-PH+scale_fill_gradientn(colors = c("blue", "white", "red"))+theme(axis.text.y = element_text(size = 10))
ggsave(PH2,file="/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_heatmap.png",width = 12,height = 10)
ggsave(PH,file="/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_heatmap0.png",width = 12,height = 10)

target_1K@assays$RNA@scale.data <- scale(target_1K@assays$RNA@data, scale = TRUE)
#target_1K <- ScaleData(object = target_1K, features = rownames(target_1K))
PH<-DoHeatmap(object = target_1K, features =major_submaker,draw.lines = TRUE,label = F,
              disp.min = -2.5, disp.max = 2.5,assay = "RNA",slot = "scale.data", size = 3, angle = -50, hjust=0.8,
              group.by ="final_major_subgroup_brief",group.bar.height = 0.1,group.colors=ppCor[1:length(levels(target_1K$final_major_subgroup_brief))]) 
PH2<-PH+scale_fill_gradientn(colors = c("lightblue", "white", "red"))+theme(axis.text.y = element_text(size = 10))
PH2
ggsave(PH2,file="/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_heatmap1.png",width = 15,height = 10)
ggsave(PH,file="/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_heatmap2.png",width = 15,height = 10)

#手动提取数据绘制所有细胞热图
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
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief,levels = subgroup_order)
gene_features <- major_submaker
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
range(mat3)#-1.715613 63.262483
mat4<-t(mat3)
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1.5 #小于负数时，加括号！
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
                 #show_row_names = FALSE,
                 show_row_names = TRUE,
                 column_split = cluster_info,
                 top_annotation = top_anno,
                 #  right_annotation = row_anno,
                 column_title = NULL,
                 col = col_fun,
                 heatmap_legend_param = list(legend_height = unit(4, "cm"),grid_width = unit(0.4, "cm"),title = "log2(count+1)",title_position = "leftcenter-rot",at = c(-2, 0, 2), labels = c("low", "median", "high") ))
pdf("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_heatmap1.pdf",width = 15,height = 8)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("major_submaker_heatmap_SCT_data"))
dev.off()
pdf("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_heatmap1.png",width = 15,height = 8)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("major_submaker_heatmap_SCT_data"))
dev.off()

##换用pheatmap
#pheatmap画，在布局上可以自由发挥

colanno=target_final@meta.data[,c("sample_code","final_major_subgroup_brief")]
colanno$final_major_subgroup_brief=factor(colanno$final_major_subgroup_brief,levels = subgroup_order)
colanno$sample_code<-as.character(colanno$sample_code)
colanno=colanno%>%arrange(final_major_subgroup_brief)
#rownames(colanno)=colanno$CB;colanno$CB=NULL
#先对细胞进行排序，按照celltype的顺序，然后对基因排序
#rowanno=markerdf1
#
#提取scale矩阵的行列时，按照上面的顺序
mark_gene<-major_submaker

celltype_names<-c(rep("no_Trophoblast",2),rep("Trophoblast_KRT7",2),rep("Myeloid_Cell_AIF1",2),rep("Stromal_cells_DCN",2),
                  "Endothelial_Cell_PECAM1", rep("Epithelial_Cell_EPCAM_PAEP",2),"Natural_killer_Cell_NCAM1",
                  rep("T_Cell_CD3D",3),"B_Cell_CD79A","Erythrocyte_HBA1", "Mast_Cell_MS4A2")
rowanno<-data.frame(gene=mark_gene,major_cluster=celltype_names)
rowanno$major_cluster=factor(rowanno$major_cluster,levels = c("no_Trophoblast",subgroup_order))
rowanno=rowanno%>%arrange(major_cluster)

mat2 <- mat1[major_submaker,]
#手动scale
mat3<-scale(t(mat2), center=T,scale=T)
mat4<-t(mat3)
range(mat4)
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1.5 #小于负数时，加括号！
#下面就是绘图代码了，我加了分界线，使其看上去更有区分度
mat4[1:4,1:4]

p1<-pheatmap(mat4[,rownames(colanno)],cluster_rows = F,cluster_cols = F,
         show_colnames = F,annotation_col = colanno,
         #gaps_row=as.numeric(cumsum(table(rowanno$major_cluster))[-11]),
         gaps_col=as.numeric(cumsum(table(colanno$final_major_subgroup_brief))[-20]))
#,filename="heatmap.2.pdf",width=11,height = 7
pdf("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_pheatmap2.pdf",width = 12,height = 8)
print(p1)
dev.off()
pdf("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_pheatmap2.png",width = 12,height = 8)
print(p1)
dev.off()
##==========================
#old and young split
Idents(object = target_final) <- 'Treat'
levels(target_final)
#DoHeatmap(object = target_sub, features = as.character(Aver_plot$gene),draw.lines = TRUE,group.colors=ppCor_all[1:length(levels(Aver_plot$cluster))]) 
#The following features were omitted as they were not found in the scale.data slot for the SCT assay
cells.1 <- WhichCells(target_final, idents ='Abortion')
cells.2 <- WhichCells(target_final,  idents = 'CTRL')
Idents(object = target_final) <- 'final_major_subgroup_brief'
makergenes<-major_submaker
hm1 <- DoHeatmap(target_final,features =makergenes, cells = cells.1, draw.lines = TRUE,label = F,group.colors=ppCor_all[1:length(levels(target_final$final_major_subgroup_brief))])
hm2 <- DoHeatmap(target_final,features = makergenes, cells = cells.2, draw.lines = TRUE,label = F,group.colors=ppCor_all[1:length(levels(target_final$final_major_subgroup_brief))])
hm1<-hm1+theme(axis.text.y = element_blank())
hm2<-hm2+theme(axis.text.y = element_blank())
hm22<-CombinePlots(plots = list(hm1, hm2),ncol = 2)
ggsave(hm22,file="/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_submaker_treat_split_heatmap.png",width = 22,height = 10)

