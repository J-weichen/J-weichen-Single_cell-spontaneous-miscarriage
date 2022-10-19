rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
dyn.load('/mnt/data/chenwei/software/hdf/hdf5/lib/libhdf5.so.200')
dyn.load('/mnt/data/chenwei/software/hdf/hdf5/lib/libhdf5_hl.so.200')

library(hdf5r)
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
library(grid)
library(gridExtra)
library(dplyr)
library(reshape2)
library(psych)
library(pheatmap)
library(ComplexHeatmap)

##set colour
x <- c(30,4,1,2,3,20,
       26,29,37,41,6,
       7,8,51,39,42,
       56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)

#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
subgroup_order<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order,ordered=TRUE)

#extract metadata information
metadata_Cell<-target_final@meta.data
final_umap<-data.frame(Embeddings(target_final[["umap"]]))
#final_umap$raw_id <- rownames(final_umap)
head(metadata_Cell);head(final_umap)
metadata_Cell_Umap_add <- merge(metadata_Cell,final_umap,by=0)
head(metadata_Cell_Umap_add)
rownames(metadata_Cell_Umap_add)<-metadata_Cell_Umap_add$Row.names
write.table(metadata_Cell_Umap_add, file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/final_metadata_of_scRNA_seq.txt",quote=F, row.names=T, col.names=T,sep="\t") 

#major figure 1b plot 
Umap_anno_plot<-DimPlot(target_final,group.by = "re_annotation_TFac",label = T,cols = my_morandi_colors)
ggsave(Umap_anno_plot,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1b_Umap_anno_plot.png",width = 11,height = 10)
ggsave(Umap_anno_plot,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1b_Umap_anno_plot.pdf",width = 11,height = 10)
Umap_anno_plot<-DimPlot(target_final,group.by = "re_annotation_TFac",label =F,cols = my_morandi_colors)
Umap_anno_plot2<-Umap_anno_plot+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
ggsave(Umap_anno_plot2,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1b_Umap_anno_plot2.png",width = 10.5,height = 10)
Umap_anno_plot<-DimPlot(target_final,group.by = "re_annotation_TFac",label =F,cols = my_morandi_colors)+NoLegend()
Umap_anno_plot3<-Umap_anno_plot+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
Umap_anno_plot4<-Umap_anno_plot+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), 
                                                            panel.grid.minor = element_blank(),
                                                            plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                      axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                      axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
ggsave(Umap_anno_plot3,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1b_Umap_anno_plot3.png",width = 8,height = 8)
ggsave(Umap_anno_plot4,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1b_Umap_anno_plot4.png",width = 8,height = 8)


#plot2 heatmaps for major makers
major_submaker<-c("PCNA","PAGE4","CDH1","ERVFRD-1","CYP19A1","EGFR","MKI67","HLA-G","PAPPA2",
                  "PAEP","PECAM1","LYVE1","AIF1","DCN","DLK1","NCAM1","CD3D","CD79A","MS4A2","HBA1")
#手动提取数据绘制所有细胞热图
levels(target_final)
DefaultAssay(target_final)<-"SCT"
FeaturePlot(object = target_final, features = c("PAPPA2"),split.by = "sample",cols= c("grey", "purple"))
FeaturePlot(object = target_final, features = c("HLA-G"),split.by = "sample",cols= c("grey", "red"))
FeaturePlot(object = target_final, features = c("CDH1","EGFR"),cols= c("grey", "red"))
FeaturePlot(object = target_final, features = c("KIT","MS4A2"),cols= c("grey", "red"))
proliferate_markers<-c("MKI67","TOP2A","TK1","PCNA")
FeaturePlot(object = target_final, features =proliferate_markers,cols= c("grey", "red"))
FeaturePlot(object = target_final, features =proliferate_markers,split.by = "Treat",cols= c("grey", "red"))

mat0<- GetAssayData(target_final, slot = "data")
gene_features <- major_submaker
#cluster_info <- target_final$re_annotation_TFac
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
range(mat3)# -1.122253 63.262483
mat4<-t(mat3)
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1.5 #小于负数时，加括号！
#set color for cell types
col <- ppCor_all[1:length(levels(cluster_info))]
#colourCount<-length(unique(target$sub_cluster))
names(col) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col),labels = levels(cluster_info), labels_gp = gpar(cex = 0.5, col = "black")))

#legend adjust
col_fun1  <- circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red"))
col_fun2  <- circlize::colorRamp2(c(-2,0,2), c("purple", "black", "yellow"))
col_fun<-col_fun1
#show_col(col_fun(seq(-3, 3)))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/major_submaker_heatmap.pdf",width = 15,height = 8)
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

plot_HP
dev.off()

#绘制Umap coordination
DefaultAssay(target_final) <- "RNA"
target_final <- NormalizeData(target_final, verbose = FALSE)
head(target_final@meta.data)
table(target_final$re_annotation_TFac)
names(table(target_final$re_annotation_TFac))


major_submaker0<-c("VIM","HLA-B","KRT7","PAGE4","PCNA","MKI67","EGFR","CDH1","CYP19A1","ERVFRD-1","HLA-G","PAPPA2","PAEP","PECAM1","DCN","DLK1","AIF1","LYVE1","NCAM1","CD3D","CD79A","MS4A2","PF4","HBA1")
#有5个可视化方法，分别是：坐标映射图，峰峦图,小提琴图，，气泡图，热图。
#后三类也可参考https://www.cnblogs.com/TOP-Bio/p/14471750.html
##supplement figure 1a plot 
plot_umap_maker <- FeaturePlot(object = target_final, features = major_submaker0,cols= c("grey", "red"),ncol=5)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/S1a_Umap_final_major_submaker.png", plot_umap_maker,width=30, height=30)
plot_umap_maker <- FeaturePlot(object = target_final, features = major_submaker0,cols= c("grey", "red"),ncol=4)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/S1a_Umap_final_major_submaker2.png", plot_umap_maker,width=24, height=26)

##hightlight postive cells
target_gene<-major_submaker0
#target_gene<-Three_linege_maker
length(target_gene)#24
target_gene3<-target_gene[which(target_gene %in% rownames(target_final))]
length(target_gene3)#24

target_expression <- GetAssayData(object = target_final,assay = "RNA", slot = "data")[target_gene3,]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))
metadata_Cell_Umap_add2<-metadata_Cell_Umap_add[rownames(target_expression),]
merge_data2 <- merge(metadata_Cell_Umap_add2, target_expression,by=0,sort=FALSE)
head(merge_data2)
#绘制基因活性的 umap distribution
gene_list<-list();gene_list2<-list()
for ( genename in target_gene3){
  # genename <-"PAEP"
  print(genename)
  #for gene expression
  express_merge_data3<-merge_data2[order(merge_data2[,genename],decreasing = T),]
  express_merge_target<- express_merge_data3[,c("UMAP_1","UMAP_2",genename)]
  colnames(express_merge_target)<-c("UMAP_1","UMAP_2","Expression")
  head(express_merge_target)
  express_merge_target$selected_cluster<-"other"
  express_merge_target[which(express_merge_target$Expression>0),]$selected_cluster<-"Pos"
  express_merge_target$selected_cluster <-factor(express_merge_target$selected_cluster,levels=c("Pos","other"))
  pos_cell_umap<-ggplot(data=express_merge_target[order(express_merge_target$selected_cluster,decreasing = T),],
                        mapping=aes(x=UMAP_1,y=UMAP_2,color = Expression))+
    geom_point(stat= "identity",size=1,alpha=0.8,show.legend = TRUE)+
    scale_color_gradient2(low="#4393C3", mid="pink", high="red", midpoint=max(express_merge_target$Expression)/2)+#"#156077","#f46f20"
    labs(title =paste0("Expression of ",genename))+theme_bw()+
    theme(panel.border = element_rect(colour="grey",fill=NA),
          panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
          axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=15),legend.position = "right")
  pos_cell_umap
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/umap_expression_",genename,"_pos_value.png"),pos_cell_umap,width=8, height=6)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/umap_expression_",genename,"_pos_value.pdf"),pos_cell_umap,width=8, height=6)
  gene_list<-c(gene_list,list(pos_cell_umap))
  
  pos_cell_umap2<-pos_cell_umap+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/umap_expression_",genename,"_pos_value_notext.png"),pos_cell_umap2,width=8, height=6)
  gene_list2<-c(gene_list2,list(pos_cell_umap2))
}

gene_list_merge<-grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],
                              gene_list[[5]], gene_list[[6]],gene_list[[7]], gene_list[[8]],
                              gene_list[[9]], gene_list[[10]],gene_list[[11]], gene_list[[12]],
                              gene_list[[13]], gene_list[[14]],gene_list[[15]], gene_list[[16]],
                              gene_list[[17]], gene_list[[18]],gene_list[[19]], gene_list[[20]],
                              gene_list[[21]], gene_list[[22]],gene_list[[23]], gene_list[[24]],
                              ncol=4)
gene_list_merge2<-grid.arrange(gene_list2[[1]], gene_list2[[2]],gene_list2[[3]],gene_list2[[4]],
                               gene_list2[[5]], gene_list2[[6]],gene_list2[[7]],gene_list2[[8]],
                               gene_list2[[9]], gene_list2[[10]],gene_list2[[11]],gene_list2[[12]],
                               gene_list2[[13]], gene_list2[[14]],gene_list2[[15]],gene_list2[[16]],
                               gene_list2[[17]], gene_list2[[18]],gene_list2[[19]], gene_list2[[20]],
                               gene_list2[[21]], gene_list2[[22]],gene_list2[[23]], gene_list2[[24]],
                               ncol=4)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1a_final_target_major_umap_expression2.png",gene_list_merge,width=24, height=26)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1a_final_target_major_umap_expression_notext2.png",gene_list_merge2,width=24, height=26)

plot_Ridge_maker <-RidgePlot(target_final,group.by = "re_annotation_TFac", features = major_submaker0, ncol = 5)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Ridge_final_major_submaker.png", plot_Ridge_maker,width=40, height=40)

#plot3 
#ref:https://www.rdocumentation.org/packages/Seurat/versions/3.1.0/topics/FetchData
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","Mycs","HCs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)

major_submaker1  <- c("KRT7","PAGE4","PCNA","CDH1","EGFR","ERVFRD-1","CYP19A1","MKI67","HLA-G","PAPPA2","PAEP","PECAM1","AIF1","LYVE1","DCN","DLK1","NCAM1","CD3D","CD79A","MS4A2","VIM","HLA-B","HBA1")

maker_Vln_plots  <- VlnPlot(target_final, features = rev(major_submaker1),group.by = "re_annotation_TFac",cols = my_morandi_colors,pt.size = 0, combine = FALSE)+NoLegend()
maker_Vln_plots2 <- CombinePlots(plots = maker_Vln_plots, ncol =1)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Violin_final_major_submaker.png", maker_Vln_plots2,width=10, height=45)

##supplement figure 1b plot 
exprs <- data.frame(FetchData(object = target_final,slot = "data",vars = c(major_submaker1,"Treat","re_annotation_TFac")))
data_plot <- melt(exprs)
data_plot$re_annotation_TFac<- factor(data_plot$re_annotation_TFac,levels = subgroup_order0)
maker_plot_up<-ggplot(data_plot,aes(x = re_annotation_TFac,y=value,fill = re_annotation_TFac)) + 
  geom_violin(scale = "width") + #geom_jitter(size = 0.00)+
  facet_grid(variable ~ .,scale='free') +  theme_bw()+scale_fill_manual(values=my_morandi_colors[1:19])+
  labs(x = "Expression level", y = "main_Cell_type", title = "main_cluster marker")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y.left = element_blank(),axis.ticks.x.bottom = element_blank())+
  theme(axis.title = element_text(size=15),axis.text.x = element_text(face="bold",size=15,angle = 90,vjust =0.5,hjust=1), axis.text.y = element_text(face="italic",size=15),legend.position = "right")+
  theme(strip.text=element_text(colour = "red", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))
maker_plot_up
data_plot$re_annotation_TFac<- factor(data_plot$re_annotation_TFac,levels = rev(subgroup_order0))
maker_plot_line<-ggplot(data_plot,aes(x = re_annotation_TFac,y=value,fill = re_annotation_TFac)) + 
  geom_violin(scale = "width") +  #geom_jitter(size = 0.00)+
  coord_flip()+ facet_grid(. ~ variable,scale='free') +
  theme_bw()+scale_fill_manual(values=rev(my_morandi_colors[1:19]))+
  labs(y = "Expression level", x = "main_Cell_type", title = "main_cluster marker")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x.bottom = element_blank(),axis.ticks.y.left = element_blank())+
  theme(axis.title = element_text(size=10),axis.text.x = element_text(face="italic",size=4), axis.text.y = element_text(face="bold",size=6),legend.position = "right")+
  theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))
maker_plot_line

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1b_Violin_final_major_submaker_up.pdf", maker_plot_up,width=12, height=18)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1b_Violin_final_major_submaker_line.pdf", maker_plot_line,width=20, height=12)

#有人认为应该按照seurat中一样的参数，使用scale ="width",同时添加随机误差，但是Y叔有一篇文章对使用scale = "width"持反对态度,本人也持反对态度
#ref1:https://www.jianshu.com/p/be67a9511bef
#ref2:https://mp.weixin.qq.com/s/tSOR68y7F0CCkAV6pUoV8Q
#Aver_plot<-top10

#2.经典气泡图
Idents(object = target_final) <- "re_annotation_TFac"
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels = rev(subgroup_order0))

Plotdot1 <-DotPlot(target_final, features = major_submaker1,group.by ="re_annotation_TFac",cols = c("blue","red")) + RotatedAxis()
Plotdot2 <-DotPlot(target_final, features = major_submaker1,group.by ="re_annotation_TFac",split.by = "Treat",cols =  c("blue","red"),assay = "RNA") + RotatedAxis()
Plotdot3 <-DotPlot(target_final, features = major_submaker1,split.by = "Treat",cols = my_morandi_colors[c(21,1)],assay = "RNA") + RotatedAxis()
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Plotdot_final_major_submaker1.pdf", Plotdot1,width=10, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Plotdot_final_major_submaker2.pdf", Plotdot2,width=10, height=15)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Plotdot_final_major_submaker3.pdf", Plotdot3,width=10, height=15)


##DotPlot
target_gene<-major_submaker1
Plotdot2 <-DotPlot(target_final, features = target_gene,group.by ="re_annotation_TFac",cols =  c("blue","red"),assay = "SCT") #+ RotatedAxis()
Plotdot2+  coord_flip()
Plotdot32<-Plotdot2+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2) +
  scale_colour_viridis(option="magma") +  coord_flip()+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(plot.title = element_text(size=8,colour = "black",face = "bold"),
        axis.title.x = element_text(size=8,colour = "black",face = "bold"),
        axis.title.y = element_text(size=8,colour = "black",face = "bold"),
        axis.text.x = element_text(size=8,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=8,colour = "black"))
Plotdot32
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Plotdot_final_major_submaker4.pdf",Plotdot32,width=6, height=6)


##3. heatmap for averageExpression
Aver_plot<-major_submaker1
AverageExp<-AverageExpression(target_final,assays="RNA",slot = "data",features=Aver_plot)
typeof(AverageExp);
head(AverageExp$RNA)
range(AverageExp$RNA)#0.000 1775.358
AverageExp_data<-log2(AverageExp$RNA+1)
range(AverageExp_data);head(AverageExp_data)
#AverageExp_data[AverageExp_data == "Inf"]<- 2000
colnames(AverageExp_data)
AverageExp_data<-AverageExp_data[,subgroup_order0]

p1<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
             show_rownames = T,show_colnames = T,border=TRUE,border_color ="grey", 
             scale ="row", color = colorRampPalette(c("navy","white","firebrick3"))(60),
             main ="log2(mean of expression level of maker + 1):scale by row")

p2<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
             show_rownames = T,show_colnames = T,border=TRUE,border_color ="grey", 
             scale ="row", color = colorRampPalette(c("navy","blue","white","orange","firebrick3"))(60),
             main ="log2(mean of expression level of maker + 1):scale by row")
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/major_submaker_average_heatmap1.pdf",width =8,height = 8)
print(p1)
dev.off()
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/major_submaker_average_heatmap2.pdf",width =8,height = 8)
print(p2)
dev.off()

#similarity of all cell types
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
p3<-pheatmap(coorda$r)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/cell_type_cluster_by_major_submaker_average_heatmap.pdf",width =8,height = 8)
print(p3)
dev.off()

##4.整体细胞类群表达水平绘制
names(target_final@assays)
DefaultAssay(target_final)<-"SCT"
scaled_data <- target_final@assays$SCT@scale.data;dim(scaled_data)# 4947 45800

makergenes<-major_submaker1
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

target_1K$re_annotation_TFac <- factor(x =target_1K$re_annotation_TFac, levels = subgroup_order0)
#heatmap for all cells
PH<-DoHeatmap(object = target_1K, features =gene_panel_intersect,draw.lines = TRUE,label = F,
              disp.min = -2.5, disp.max = 2.5,assay = "SCT",size = 3, angle = -50, hjust=0.8,
              group.by ="re_annotation_TFac",group.bar.height = 0.1,group.colors=my_morandi_colors[1:length(levels(target_1K$re_annotation_TFac))]) 
PH2<-PH+scale_fill_gradientn(colors = c("blue", "white", "red"))+theme(axis.text.y = element_text(size = 10))
ggsave(PH2,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/major_submaker_1k_SCT_heatmap.png",width = 12,height = 10)
ggsave(PH,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/major_submaker_1k_SCT_heatmap0.png",width = 12,height = 10)

target_1K@assays$RNA@scale.data <- scale(target_1K@assays$RNA@data, scale = TRUE)
#target_1K <- ScaleData(object = target_1K, features = rownames(target_1K))
PH<-DoHeatmap(object = target_1K, features =gene_panel_intersect,draw.lines = TRUE,label = F,
              disp.min = -2.5, disp.max = 2.5,assay = "RNA",slot = "scale.data", size = 3, angle = -50, hjust=0.8,
              group.by ="re_annotation_TFac",group.bar.height = 0.1,group.colors=my_morandi_colors[1:length(levels(target_1K$re_annotation_TFac))]) 
PH2<-PH+scale_fill_gradientn(colors = c("lightblue", "white", "red"))+theme(axis.text.y = element_text(size = 10))
PH2
ggsave(PH2,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/major_submaker_1k_RNA_heatmap1.png",width = 12,height = 10)
ggsave(PH,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/major_submaker_1k_RNA_heatmap2.png",width = 12,height = 10)

##==========================
#old and young split
Idents(object = target_final) <- 'Treat'
levels(target_final)
target_final$re_annotation_TFac <- factor(x =target_final$re_annotation_TFac, levels = subgroup_order0)
cells.1 <- WhichCells(target_final, idents ='Abortion')
cells.2 <- WhichCells(target_final,  idents = 'CTRL')
Idents(object = target_final) <- 're_annotation_TFac'
makergenes<-major_submaker1
hm1<-DoHeatmap(object = target_final, features =makergenes,cells = cells.1,draw.lines = TRUE,label = F,
              disp.min = -2.5, disp.max = 2.5,assay = "RNA",slot = "data", size = 3, angle = -50, hjust=0.8,
              group.by ="re_annotation_TFac",group.bar.height = 0.1,group.colors=my_morandi_colors[1:length(levels(target_final$re_annotation_TFac))]) 

hm1<-hm1+scale_fill_gradientn(colors = c("lightblue", "white", "red"))+theme(axis.text.y = element_text(size = 10))
hm1
hm2<-DoHeatmap(object = target_final, features =makergenes,cells = cells.2,draw.lines = TRUE,label = F,
               disp.min = -2.5, disp.max = 2.5,assay = "RNA",slot = "data", size = 3, angle = -50, hjust=0.8,
               group.by ="re_annotation_TFac",group.bar.height = 0.1,group.colors=my_morandi_colors[1:length(levels(target_final$re_annotation_TFac))]) 

hm2<-hm2+scale_fill_gradientn(colors = c("lightblue", "white", "red"))+theme(axis.text.y = element_text(size = 10))
hm2
hm22<-CombinePlots(plots = list(hm1, hm2),ncol = 2)
ggsave(hm22,file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/major_submaker_RNA_treat_split_heatmap.png",width = 30,height = 10)
