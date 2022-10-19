rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(data.table)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(scales)
library(ggsci)
library(data.table)
library(grid)
library(gridExtra)
library(cowplot)
library(psych)
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

#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
library(future)
plan(strategy = "multicore", workers = 42)
options(future.globals.maxSize = 1500 * 1024^12)

target<-readRDS(file ="/mnt/data/chenwei/jianghai/3.seurat_result/Cell_0518.rds")
head(target@meta.data)
#plot  group
DimPlot(target, group.by = "major_cluster",label = F,cols = ppCor_all2)

DimPlot(target, group.by = "merge_cluster",label = F,cols = ppCor_all2)
DimPlot(object = target, group.by = "merge_cluster", split.by = "Treat",reduction = "umap",cols = ppCor_all2,ncol = 2)
DimPlot(target, group.by = "Tissue",label = F,cols = ppCor_all[c(3,5)])

#DefaultAssay(object = vitDcca) <- "RNA"
#vitDcca <- NormalizeData(object =vitDcca, normalization.method = "LogNormalize", scale.factor = 10000)
#all.genes <- rownames(x = vitDcca)
#vitDcca <- ScaleData(object = vitDcca, features = all.genes)

#有5个可视化方法，分别是：坐标映射图，峰峦图,小提琴图，，气泡图，热图。
#后三类也可参考https://www.cnblogs.com/TOP-Bio/p/14471750.html
FeaturePlot(target, features = c("MS4A1", "CD79A"))
RidgePlot(target, features = c("MS4A1", "CD79A"), ncol = 1)

#maker plot
#小提琴图
#ref:https://www.rdocumentation.org/packages/Seurat/versions/3.1.0/topics/FetchData
#CR maker
maker_plots <- VlnPlot(target, features = rev(c("VIM","KRT7")),group.by = "major_cluster",cols = ppCor_all2,pt.size = 0, combine = FALSE)
CombinePlots(plots = maker_plots, ncol =1)
classic_markers.to.plot <- c("VIM","KRT7","CD68","AIF1","CD3D","NKG7","FGFBP2","CD79A","MS4A1","HBA1","PPBP")
#GPB maker
classic_markers.to.plot <- c("CD79A","CD3D","CD8A","SELL","CCR7","S100A4","FOXP3","ZBTB16",
                              "KLRB1","FGFBP2","FCGR3A","NCAM1","ENTPD1","CD160","KIT","CD14",
                              "LYZ","FCER1A","TPSAB1")


exprs <- data.frame(FetchData(object = target, vars = c(classic_markers.to.plot,"Treat","major_cluster")))
data_plot <- melt(exprs)
data_plot$major_cluster<- factor(data_plot$major_cluster,levels = c("Trophoblast_KRT7","myeloid_Cell_AIF1_pos","T_Cell_CD3D_pos","NK_Cell_GNLY_pos_CD8A_neg",
                                                         "B_Cell_CD79A_pos","Erythrocyte_HBA1_pos","Megakaryocytes_PPBP_pos","Mix_cell_Erythrocyte_HBA1_pos"))
maker_plot_line<-ggplot(data_plot,aes(x = major_cluster,y=value,fill = major_cluster)) + 
  geom_violin(scale = "width") +  #geom_jitter(size = 0.00)+
  coord_flip()+ facet_grid(. ~ variable,scale='free') +
  theme_bw()+scale_fill_manual(values=ppCor)+
  labs(y = "Expression level", x = "main_Cell_type", title = "main_cluster marker")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x.bottom = element_blank(),axis.ticks.y.left = element_blank())+
  theme(axis.title = element_text(size=10),axis.text.x = element_text(face="italic",size=4), axis.text.y = element_text(face="italic",size=6),legend.position = "right")+
  theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))
maker_plot_up<-ggplot(data_plot,aes(x = major_cluster,y=value,fill = major_cluster)) + 
  geom_violin(scale = "width") + #geom_jitter(size = 0.00)+
  facet_grid(variable ~ .,scale='free') +  theme_bw()+scale_fill_manual(values=ppCor)+
  labs(x = "Expression level", y = "main_Cell_type", title = "main_cluster marker")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y.left = element_blank(),axis.ticks.x.bottom = element_blank())+
  theme(axis.title = element_text(size=15),axis.text.x = element_text(face="italic",size=15), axis.text.y = element_text(face="italic",size=15),legend.position = "right")+
  theme(strip.text=element_text(colour = "red", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))
maker_plot_up
maker_plot_line

#有人认为应该按照seurat中一样的参数，使用scale ="width",同时添加随机误差，但是Y叔有一篇文章对使用scale = "width"持反对态度,本人也持反对态度
#ref1:https://www.jianshu.com/p/be67a9511bef
#ref2:https://mp.weixin.qq.com/s/tSOR68y7F0CCkAV6pUoV8Q
#Aver_plot<-top10
Idents(object = target) <- "major_cluster"
#2.经典气泡图
DotPlot(target, features = unique(classic_markers.to.plot),cols = c("blue","red")) + RotatedAxis()
DotPlot(target, features = unique(classic_markers.to.plot),split.by = "Treat",cols = ppCor,assay = "SCT") + RotatedAxis()

##3. heatmap for averageExpression
Aver_plot<-classic_markers.to.plot
AverageExp<-AverageExpression(target,features=unique(Aver_plot))
typeof(AverageExp);
head(AverageExp$RNA)
range(AverageExp$RNA)
AverageExp_data<-log2(AverageExp$RNA+1)
range(AverageExp_data);head(AverageExp_data)
AverageExp_data[AverageExp_data == "Inf"]<- 2000
colnames(AverageExp_data)
AverageExp_data<-AverageExp_data[,c("Trophoblast_KRT7","Myeloid_Cell_AIF1_pos","T_Cell_CD3D_pos","NK_Cell_CD3D_neg",
                                    "B_Cell_CD79A_pos","Erythrocyte_HBB_pos","Megakaryocytes_PPBP_pos")]

pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
         show_rownames = T,show_colnames = T,border=TRUE,border_color ="grey", 
         scale ="row", color = colorRampPalette(c("navy","white","firebrick3"))(60),
         main ="log2(mean of expression level of maker + 1):scale by row")

pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
         show_rownames = T,show_colnames = T,border=TRUE,border_color ="grey", 
         scale ="row", color = colorRampPalette(c("navy","blue","white","orange","firebrick3"))(60),
         main ="log2(mean of expression level of maker + 1):scale by row")

#similarity of all cell types
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap(coorda$r)

##4.整体细胞类群表达水平绘制
names(target@assays)
DefaultAssay(target)<-"SCT"
scaled_data <- target@assays$SCT@scale.data;dim(scaled_data)# 2981 53232
makergenes<-unique(Aver_plot)
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

target_1K<-subset(target, downsample = 1000)
scaled_data_1K <- target_1K@assays$SCT@scale.data
range(scaled_data_1K)

celltype_order<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1_pos","T_Cell_CD3D_pos","NK_Cell_CD3D_neg",
                  "B_Cell_CD79A_pos","Erythrocyte_HBB_pos","Megakaryocytes_PPBP_pos")

target_1K$major_cluster <- factor(x =target_1K$major_cluster, levels = celltype_order)
#Idents(object = target_1K) <- 'major_cluster'
#heatmap for all cells
PH<-DoHeatmap(object = target_1K, features =unique(Aver_plot),draw.lines = TRUE,label = F,
              disp.min = -2.5, disp.max = 2.5,assay = "SCT",size = 3, angle = -50, hjust=0.8,
              group.by ="major_cluster",group.bar.height = 0.1,group.colors=ppCor[1:length(levels(target_25K$major_cluster))]) 
PH2<-PH+scale_fill_gradientn(colors = c("blue", "white", "red"))+theme(axis.text.y = element_text(size = 10))+NoLegend()
ggsave(PH2,file="/mnt/data/chenwei/jianghai/3.seurat_result/2.DEGs_call/heatmap.png",width = 12,height = 10)
ggsave(PH,file="/mnt/data/chenwei/jianghai/3.seurat_result/2.DEGs_call/heatmap0.png",width = 12,height = 10)

#手动提取数据绘制所有细胞热图
levels(target)
DefaultAssay(target)<-"SCT"
#mat0 <- GetAssayData(target, slot = "counts")
#mat0 <- log2(mat0 + 1)
#其他SCT的数据说明：https://www.jianshu.com/p/e639cc257d51 
mat0<- GetAssayData(target, slot = "data")
#FC_10<-Allcluster.markers[Allcluster.markers$p_val_adj<0.01 & Allcluster.markers$avg_logFC>1,]%>% group_by(cluster)
#FC_1.7<-Allcluster.markers[Allcluster.markers$p_val_adj<0.01 & Allcluster.markers$avg_logFC>0.25,]%>% group_by(cluster)
#top10 <- FC_1.7 %>% top_n(n = 10, wt = avg_logFC)

#get genes and cluster information
#gene_features <-top10$gene
gene_features <- classic_markers.to.plot#top10
cluster_info <- target$major_cluster

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
range(mat3)
mat4<-t(mat3)
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
                 #right_annotation = row_anno,
                 column_title = NULL,
                 col = col_fun,
                 heatmap_legend_param = list(legend_height = unit(4, "cm"),grid_width = unit(0.4, "cm"),title = "log2(count+1)",title_position = "leftcenter-rot",at = c(-2, 0, 2), labels = c("low", "median", "high") ))
pdf("/mnt/data/chenwei/jianghai/3.seurat_result/2.DEGs_call/figure2_",cell_name,"_",def_range,"_pos_005_genes.pdf"))
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("DE in ",cell_name,"_",def_range," sub_cluster_brief"))
dev.off()


##换用pheatmap
#pheatmap画，在布局上可以自由发挥

colanno=target@meta.data[,c("major_cluster","sample")]
colanno$major_cluster=factor(colanno$major_cluster,levels = celltype_order)
colanno=colanno%>%arrange(major_cluster)
#rownames(colanno)=colanno$CB;colanno$CB=NULL
#先对细胞进行排序，按照celltype的顺序，然后对基因排序
#rowanno=markerdf1
#
#提取scale矩阵的行列时，按照上面的顺序
mark_gene<-c("KRT7","CD68","AIF1","CD3D","NKG7","FGFBP2","CD79A","MS4A1","HBA1","PPBP") 
celltype_names<-c("Trophoblast_KRT7",rep("Myeloid_Cell_AIF1_pos",2),"T_Cell_CD3D_pos",
                  rep("NK_Cell_CD3D_neg",2),rep("B_Cell_CD79A_pos",2),"Erythrocyte_HBB_pos",   
                  "Megakaryocytes_PPBP_pos")
rowanno<-data.frame(gene=mark_gene,major_cluster=celltype_names)
rowanno$major_cluster=factor(rowanno$major_cluster,levels = celltype_order)
rowanno=rowanno%>%arrange(major_cluster)

mat4=target[["SCT"]]@scale.data[as.character(rowanno$gene),rownames(colanno)]
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1.5 #小于负数时，加括号！
#下面就是绘图代码了，我加了分界线，使其看上去更有区分度
mat4[1:4,1:4]
pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno,
         gaps_row=as.numeric(cumsum(table(rowanno$major_cluster))[-7]),
         gaps_col=as.numeric(cumsum(table(colanno$major_cluster))[-7])
         #,filename="heatmap.2.pdf",width=11,height = 7
)

##==========================
#old and young split
Idents(object = target) <- 'Treat'
levels(target)
#DoHeatmap(object = target_sub, features = as.character(Aver_plot$gene),draw.lines = TRUE,group.colors=ppCor_all[1:length(levels(Aver_plot$cluster))]) 
#The following features were omitted as they were not found in the scale.data slot for the SCT assay
cells.1 <- WhichCells(target, idents ='PPH')
cells.2 <- WhichCells(target,  idents = 'CTRL')
Idents(object = target) <- 'major_cluster'
makergenes<-unique(Aver_plot)
hm1 <- DoHeatmap(target,features =makergenes, cells = cells.1, draw.lines = TRUE,label = F,group.colors=ppCor_all[1:length(levels(target$major_cluster))])
hm2 <- DoHeatmap(target,features = makergenes, cells = cells.2, draw.lines = TRUE,label = F,group.colors=ppCor_all[1:length(levels(target$major_cluster))])
hm1<-hm1+theme(axis.text.y = element_blank())
hm2<-hm2+theme(axis.text.y = element_blank())
hm22<-CombinePlots(plots = list(hm1, hm2),ncol = 2)
ggsave(hm22,file="/mnt/data/chenwei/jianghai/3.seurat_result/2.DEGs_call/treat_split_heatmap.png",width = 22,height = 10)
