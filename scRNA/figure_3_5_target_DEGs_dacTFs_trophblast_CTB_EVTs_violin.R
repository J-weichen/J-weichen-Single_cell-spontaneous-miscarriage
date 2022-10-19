##提取某条KEGG通路中的所有基因
##ref：https://zhuanlan.zhihu.com/p/434383719
##更详细解析：https://www.shenxt.info/post/2019-11-14-kegg/ 
#kegg
#1, hippo pathway:hsa04390
#2, TGF-b pathway: hsa04350
#3, Wnt pathway: hsa04310
rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(data.table)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(grid)
library(gridExtra)
library(viridis)
library(ggpubr)
library(purrr)
library(cowplot)
grid.newpage()
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:46))]
#reading expression data
#########################
###KS_panc_fl_data
########################
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)
Umap_coordinate<-data.frame(Embeddings(target_final[["umap"]]))
target_final <- AddMetaData(target_final, Umap_coordinate)
DimPlot(target_final, group.by = "re_annotation_TFac",label = F,cols = ppCor_all2)

##substract Trophoblast
Troph_group<-c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3")#"STBs_1","STBs_2","STBs_3",
Troph_target<-subset(x = target_final, subset = re_annotation_TFac %in% Troph_group)
Troph_target$re_annotation_TFac<- factor(Troph_target$re_annotation_TFac, levels=Troph_group,ordered=TRUE)

##reannotation for CTBs and EVTs
Troph_target$merge_cluster <- as.character(Troph_target$re_annotation_TFac)
sub_cell<-subset(x = Troph_target,subset = re_annotation_TFac %in% c("CTBs_1","CTBs_2"))
Troph_target$merge_cluster[Cells(sub_cell)] <- "CTBs"
sub_cell<-subset(x = Troph_target,subset = re_annotation_TFac %in% c("EVTs_1","EVTs_2","EVTs_3"))
Troph_target$merge_cluster[Cells(sub_cell)] <- "EVTs"
DimPlot(Troph_target, group.by = "merge_cluster",label = F,cols = ppCor_all2)
Troph_target$merge_cluster<- factor(Troph_target$merge_cluster, levels=c("CTBs","EVTs"),ordered=TRUE)

##calculation 
DefaultAssay(Troph_target) <- "RNA"
# Normalize RNA data for visualization purposes
Troph_target <- NormalizeData(Troph_target, verbose = FALSE)
express_data <- as.matrix(GetAssayData(Troph_target, slot = "data"))

##for merge CTBs and EVTs
Group_tag<-"merge_USM_TF_select_0531"
gene_signature<-c("TP53","GATA3","STAT3")
gene_signature<-c("NR5A1","SOX4","STAT3")
Group_tag<-"merge_USM_TF_select_0705"
gene_signature<-c("STAT5B","STAT5A","STAT3","CAVIN1")

Group_tag<-"merge_USM_TF_select_0713"
gene_signature<-c("ELF3","GATA3","STAT3")

Group_tag<-"merge_USM_TF_select_0713"
gene_signature<-c("ARNT2","ATF3","SMAD3","TEAD1","ELF3","GATA2","HSF1","SNAI1","TFAP2A","TEF")

Group_tag<-"merge_USM_TF_select_0713"
gene_signature<-c("ARNT2","ATF3","SMAD3","TEAD1","ELF3","GATA2","HSF1","SNAI1","TFAP2A","TEF")

##expression matrix for yaxi for figure four
target_exprs <- data.frame(FetchData(object = Troph_target,slot = "data",vars = c("Treat","merge_cluster",gene_signature)))
target_exprs$Treat<-as.character(target_exprs$Treat)
target_exprs[which(target_exprs$Treat =="Abortion"),]$Treat<-"USM"
target_exprs$group<-paste(target_exprs$merge_cluster,target_exprs$Treat,sep="_")
target_exprs<-target_exprs[,-c(1:2)]
target_exprs[1:4,];dim(target_exprs)#29882    13
write.table(target_exprs,file ="/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/normalized_gene_expression_CTBs_EVTs_TFs_regulon_activity.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

Idents(object = Troph_target) <- "merge_cluster"
VlnPlot(Troph_target, features =gene_signature, pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 5)
VlnPlot(Troph_target, features =gene_signature[2], pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
VlnPlot(Troph_target, features =gene_signature, pt.size = 0.01, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
target_violin<-VlnPlot(Troph_target, features =gene_signature, pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot.pdf"), target_violin,width=5, height =length(gene_signature)*3)

plot1<-VlnPlot(Troph_target, features = gene_signature[2], y.max =4.5,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot1_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =4,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[4], y.max =4,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot2_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =3.5,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[9], y.max =3.1,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot3_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =3,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[5], y.max =4,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot4_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =3.5,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[6], y.max =3.5,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot5_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =3,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[7], y.max =2.8,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot6_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =2.6,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[8], y.max =2.6,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot7_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =2.5,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))

gene_list_merge<-grid.arrange(plot1_1,plot2_1,plot3_1,plot4_1,plot5_1,plot6_1,plot7_1,ncol=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot_stat.pdf"),gene_list_merge,width=12, height=6)

##############

plot1<-VlnPlot(Troph_target, features = gene_signature[1], y.max =2.5,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot1_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =2,label = "p.signif",hide.ns = F, label.x = 2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[2], y.max =4,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot1_2<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[3], y.max =3,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))

gene_list_merge<-grid.arrange(plot1_1,plot1_2,plot1_3,ncol=1)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot_stat.pdf"),gene_list_merge,width=5, height=length(gene_signature)*3)

##plot dot
Troph_target$cell_class<-paste0(Troph_target$merge_cluster,"_",Troph_target$Treat)

##正式绘图
Plotdot <-DotPlot(Troph_target, features = gene_signature,group.by ="merge_cluster",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()
Plotdot2 <-DotPlot(Troph_target, features = gene_signature,group.by ="cell_class",cols =  c("blue","red"),assay = "RNA") #+ RotatedAxis()
Plotdot2;Plotdot
Plotdot31<-Plotdot2+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ coord_flip()+
  theme(plot.title = element_text(size=8,colour = "black",face = "bold"),
        axis.title.x = element_text(size=8,colour = "black",face = "bold"),
        axis.title.y = element_text(size=8,colour = "black",face = "bold"),
        axis.text.x = element_text(size=8,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=8,colour = "black"))

Plotdot32<-Plotdot2+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2) +
  scale_colour_viridis(option="magma") +  coord_flip()+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(plot.title = element_text(size=8,colour = "black",face = "bold"),
        axis.title.x = element_text(size=8,colour = "black",face = "bold"),
        axis.title.y = element_text(size=8,colour = "black",face = "bold"),
        axis.text.x = element_text(size=8,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=8,colour = "black"))

Plotdot02<-Plotdot2+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ 
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1, vjust=0.5))

ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_genes_dotplot_all_celltype_split.pdf"),Plotdot31,width=4, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_genes_dotplot_all_celltype_split2.pdf"),Plotdot32,width=4, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_genes_dotplot_all_celltype_split3.pdf"),Plotdot02,width=4, height=4)



###for celltypes
Idents(object = Troph_target) <- "re_annotation_TFac"

#panc_fl <- subset(x = KS_raw, idents = c("Common Progenitor","Leydig Progenitor","Leydig Precursor #1","Leydig Precursor #2","Differented Leydig Cells", "Peritubular Myoid Cells", "Sertoli Cells #1", "Sertoli Cells #2", "Sertoli Progenitor"))
umapplot<-DimPlot(Troph_target, group.by = "re_annotation_TFac", label = T, repel = T,cols = ppCor_all2)
umapplot

Group_tag<-"Six_USM_predict_TFs"
gene_signature<-c("ZEB1", "ZBTB2","ZBTB6","TBX18","SOX17","PLAGL1","HIC1")
FeaturePlot(object = target_final, features = gene_signature,cols= c("grey", "red"),split.by = "Treat",pt.size = 1)
FeaturePlot(object = Troph_target, features = gene_signature,cols= c("grey", "red"),split.by = "Treat",pt.size = 1)

FeaturePlot(object = target_final, features = "ZEB2",cols= c("grey", "red"),split.by = "Treat",pt.size = 1)


Group_tag<-"Three_USM_SNP_overlapped"
gene_signature<-c("MMP9","SERPINE1","CD320")
gene_signature<-c("CD320","VEGFA","SERPINE1")
Group_tag<-"USM_TF_select_0531"
#gene_signature<-c("NR5A1","SOX4","STAT3")
gene_signature<-c("TP53","GATA3","STAT3")

Group_tag<-"USM_interaction_gene"
gene_signature<-c("BMPR2","ACVR2B","ACVR2A","ACVR1","BMP4","WNT7A","FZD5","WNT2","LRP6","ITGA1", "SEMA7A","LEP","LEPR")

#gene_signature<-c("VEGFA","MTHFR")
#gene_signature<-c("NR5A1","SOX4","STAT3","SERPINE1","CD320")
#gene_signature<-c("SMCHD1","TNFRSF1B","TMEM176A")
VlnPlot(Troph_target, features =gene_signature[1], pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
target_violin<-VlnPlot(Troph_target, features =gene_signature, pt.size = 0.01, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot.pdf"), target_violin,width=8, height =length(gene_signature)*3)

plot1<-VlnPlot(Troph_target, features = gene_signature[1], y.max =2,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot1_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =2,label = "p.signif",hide.ns = F, label.x = 2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[2], y.max =4,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot1_2<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[3], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))

gene_list_merge<-grid.arrange(plot1_1,plot1_2,plot1_3,ncol=1)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot_stat.pdf"),gene_list_merge,width=4, height=length(gene_signature)*3)

##plot dot
Troph_target$cell_class<-paste0(Troph_target$re_annotation_TFac,"_",Troph_target$Treat)

##正式绘图
Plotdot <-DotPlot(Troph_target, features = gene_signature,group.by ="re_annotation_TFac",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()
Plotdot2 <-DotPlot(Troph_target, features = gene_signature,group.by ="cell_class",cols =  c("blue","red"),assay = "RNA") #+ RotatedAxis()
Plotdot2;Plotdot
Plotdot31<-Plotdot2+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ coord_flip()+
  theme(plot.title = element_text(size=8,colour = "black",face = "bold"),
        axis.title.x = element_text(size=8,colour = "black",face = "bold"),
        axis.title.y = element_text(size=8,colour = "black",face = "bold"),
        axis.text.x = element_text(size=8,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=8,colour = "black"))

Plotdot32<-Plotdot2+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2) +
  scale_colour_viridis(option="magma") +  coord_flip()+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(plot.title = element_text(size=8,colour = "black",face = "bold"),
        axis.title.x = element_text(size=8,colour = "black",face = "bold"),
        axis.title.y = element_text(size=8,colour = "black",face = "bold"),
        axis.text.x = element_text(size=8,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=8,colour = "black"))

Plotdot02<-Plotdot2+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ 
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1, vjust=0.5))
Plotdot02




##supplement figure 1b plot 
exprs <- data.frame(FetchData(object = Troph_target,slot = "data",vars = c(gene_signature,"Treat","re_annotation_TFac")))
data_plot <- melt(exprs)
data_plot$re_annotation_TFac<- factor(data_plot$re_annotation_TFac,levels = Troph_group)
head(data_plot)
for ( genename in gene_signature){
  # genename <-"ZEB1"
  print(genename)
  merge_data_target<- data_plot[which(data_plot$variable==genename),]
  maker_plot<-ggplot(merge_data_target,aes(x = Treat,y=value,fill = Treat)) + geom_violin(scale = "width") + 
    ylim(0,max(merge_data_target$value)+0.5)+
    facet_grid(re_annotation_TFac ~ .,scales = "free_y") +  theme_bw()+
    scale_fill_manual(values=ppCor[c(10,6)])+  stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
    geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +
    stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))+
    labs(x = "Group", y = "LogNormalized count", title =genename)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x.bottom = element_blank(),legend.position = "right")
  maker_plot
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Up_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=4, height=12)
}

##for expression value max adjust
#gene_signature<-c("CD320","VEGFA","SERPINE1")
##For CD320
genename <-"CD320"
print(genename)
#  VlnPlot(Troph_target, features =genename, pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
merge_data_target<- data_plot[which(data_plot$variable==genename),]
head(merge_data_target)
merge_data_target_original<-merge_data_target
max_value<-1
merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot1<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  # geom_violin(scale = "area",trim = TRUE) +  geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  # ylim(0,max(merge_data_target$value)+0.1)+
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot1
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)

##For VEGFA
genename <-"VEGFA"
print(genename)
#  VlnPlot(Troph_target, features =genename, pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
merge_data_target<- data_plot[which(data_plot$variable==genename),]
head(merge_data_target)
merge_data_target_original<-merge_data_target
max_value<-0.5
merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot2<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  # geom_violin(scale = "area",trim = TRUE) +  geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  # ylim(0,max(merge_data_target$value)+0.1)+
  scale_y_continuous(limits = c(0,0.5),breaks = c(0,0.1,0.2,0.3,0.4,0.5))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot2
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)

##For VEGFA
genename <-"SERPINE1"
print(genename)
#  VlnPlot(Troph_target, features =genename, pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
merge_data_target<- data_plot[which(data_plot$variable==genename),]
head(merge_data_target)
merge_data_target_original<-merge_data_target
max_value<-5
merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot3<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  # geom_violin(scale = "area",trim = TRUE) +  geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  # ylim(0,max(merge_data_target$value)+0.1)+
  scale_y_continuous(limits = c(0,5),breaks = c(0,1,2,3,4,5))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot3
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)

gene_list_merge<-grid.arrange(maker_plot1,maker_plot2,maker_plot3,ncol=1)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_CD320_VEGFA_SERPINE1_USM_SNP_overlapped_max_value_adjust_violin_target_gene_plot_stat.pdf",gene_list_merge,width=5, height=10)



##plot in all cell types

DefaultAssay(target_final) <- "RNA"
# Normalize RNA data for visualization purposes
target_final <- NormalizeData(target_final, verbose = FALSE)
express_data <- as.matrix(GetAssayData(target_final, slot = "data"))

Group_tag<-"USM_interaction_gene"
gene_signature<-c("BMPR2","ACVR2B","ACVR2A","ACVR1","BMP4","WNT7A","FZD5","WNT2","LRP6","ITGA1", "SEMA7A","LEP","LEPR")
target_violin<-VlnPlot(target_final, features =gene_signature, pt.size = 0.01, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot.pdf"), target_violin,width=8, height =length(gene_signature)*3)
target_violin<-VlnPlot(target_final, features =gene_signature, pt.size = 0.01, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 2)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot2.pdf"), target_violin,width=13, height =length(gene_signature)*3)
target_violin<-VlnPlot(target_final, features =gene_signature, pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 2)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot3.pdf"), target_violin,width=13, height =length(gene_signature)*3)


#y.max =5, 
plot1<-VlnPlot(target_final, features = gene_signature, group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot1_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1[[1]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.1)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_2<-plot1[[2]]+ stat_compare_means(method = "wilcox.test",data=plot1[[2]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_3<-plot1[[3]]+ stat_compare_means(method = "wilcox.test",data=plot1[[3]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =1.7)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_4<-plot1[[4]]+ stat_compare_means(method = "wilcox.test",data=plot1[[4]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_5<-plot1[[5]]+ stat_compare_means(method = "wilcox.test",data=plot1[[5]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_6<-plot1[[6]]+ stat_compare_means(method = "wilcox.test",data=plot1[[6]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_7<-plot1[[7]]+ stat_compare_means(method = "wilcox.test",data=plot1[[7]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_8<-plot1[[8]]+ stat_compare_means(method = "wilcox.test",data=plot1[[8]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_9<-plot1[[9]]+ stat_compare_means(method = "wilcox.test",data=plot1[[9]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.8)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_10<-plot1[[10]]+ stat_compare_means(method = "wilcox.test",data=plot1[[10]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_11<-plot1[[11]]+ stat_compare_means(method = "wilcox.test",data=plot1[[11]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_12<-plot1[[12]]+ stat_compare_means(method = "wilcox.test",data=plot1[[12]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =5.7)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_13<-plot1[[13]]+ stat_compare_means(method = "wilcox.test",data=plot1[[13]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
gene_list_merge<-grid.arrange(plot1_1,plot1_2,plot1_3,plot1_4,plot1_5,plot1_6,plot1_7,
                              plot1_8,plot1_9,plot1_10,plot1_11,plot1_12,plot1_13,ncol=2)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot_stat2.pdf"),gene_list_merge,width=14, height=16)


plot1<-VlnPlot(target_final, features = gene_signature[1], y.max =0.5,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot1_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =0.45,label = "p.signif",hide.ns = F, label.x = 2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[2], y.max =5,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot1_2<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =4.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[3], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[3], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[4], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[5], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[6], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[7], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[8], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[9], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[10], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1<-VlnPlot(target_final, features = gene_signature[11], y.max =3,  group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
plot1_3<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))

gene_list_merge<-grid.arrange(plot1_1,plot1_2,plot1_3,ncol=1)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_violin_target_gene_plot_stat.pdf"),gene_list_merge,width=4, height=length(gene_signature)*3)


#extract metadata information
metadata_Cell<-target_final@meta.data
final_umap<-data.frame(Embeddings(target_final[["umap"]]))
#final_umap$raw_id <- rownames(final_umap)
head(metadata_Cell);head(final_umap)
metadata_Cell_Umap_add <- merge(metadata_Cell,final_umap,by=0)
head(metadata_Cell_Umap_add)
rownames(metadata_Cell_Umap_add)<-metadata_Cell_Umap_add$Row.names

##plot for target genes
target_final$cell_class<-paste0(target_final$cell_type,"_",target_final$Sex)
table(target_final$cell_class)

##get genes in target 
target_gene<-c("GLI1","PTCH1","SF1","INSL3","CYP17A1","IGF2")
tag<-"qingmeng_genes"

##正式绘图
Plotdot <-DotPlot(target_final, features = target_gene,group.by ="cell_type",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()
Plotdot2 <-DotPlot(target_final, features = target_gene,group.by ="cell_class",cols =  c("blue","red"),assay = "RNA") #+ RotatedAxis()
Plotdot2;Plotdot
Plotdot31<-Plotdot2+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ coord_flip()+
  theme(plot.title = element_text(size=8,colour = "black",face = "bold"),
        axis.title.x = element_text(size=8,colour = "black",face = "bold"),
        axis.title.y = element_text(size=8,colour = "black",face = "bold"),
        axis.text.x = element_text(size=8,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=8,colour = "black"))

Plotdot32<-Plotdot2+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2) +
  scale_colour_viridis(option="magma") +  coord_flip()+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(plot.title = element_text(size=8,colour = "black",face = "bold"),
        axis.title.x = element_text(size=8,colour = "black",face = "bold"),
        axis.title.y = element_text(size=8,colour = "black",face = "bold"),
        axis.text.x = element_text(size=8,angle=90,hjust=1, vjust=0.5),
        axis.text.y = element_text(size=8,colour = "black"))
Plotdot01<-Plotdot+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ 
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1, vjust=0.5))
Plotdot02<-Plotdot2+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ 
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1, vjust=0.5))

ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_",tag,"_genes_dotplot_all_celltype_split.pdf"),Plotdot31,width=8, height=5)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_",tag,"_genes_dotplot_all_celltype_split_2.pdf"),Plotdot32,width=8, height=5)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_",tag,"_genes_dotplot_all_celltype_1.pdf"),Plotdot01,width=7, height=6)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_",tag,"_genes_dotplot_all_celltype_split_3.pdf"),Plotdot02,width=7, height=8)

##plot violin with p value
target_gene2<-target_gene[which(target_gene %in% rownames(target_final))]
target_expression <- GetAssayData(object = target_final,assay = "RNA", slot = "data")[target_gene2,]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))

metadata_Cell_Umap_add2<-metadata_Cell_Umap_add[rownames(target_expression),]
merge_data2 <- merge(metadata_Cell_Umap_add2, target_expression,by=0,sort=FALSE)
head(merge_data2)

##绘制两比较组间的表达量
gene_list<-list()
for ( genename in rev(target_gene2)){
  # genename <-"PTCH1"
  print(genename)
  merge_data_target<- merge_data2[,c("cell_type","Sex",genename)]
  colnames(merge_data_target)<-c("cell_type","Sex","Expression_level")
  plot_Cell_exp<- ggplot(merge_data_target,aes(x = cell_type,y=Expression_level,fill = Sex)) + 
    geom_violin(scale = "width") +  #geom_jitter(size = 0.00)+coord_flip()+ facet_grid(. ~ variable,scale='free') +
    theme_bw()+scale_fill_manual(values=ppCor)+
    labs(y = "Expression level", x = "main_Cell_type", title = genename)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.y.left = element_blank())+
    theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),legend.position = "right")+
    theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))+
    stat_compare_means(method = "wilcox.test",data=merge_data_target,aes(x=cell_type,y = Expression_level,group = Sex,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((merge_data_target$Expression_level)))
  plot_Cell_exp2<- plot_Cell_exp+ geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=4, col="orange", position = position_dodge(0.9))
  ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_stat_barplot_expression_",tag,"_",genename,".pdf"),plot_Cell_exp2,width=12, height=7)
  gene_list<-c(gene_list,list(plot_Cell_exp2))
}

gene_list_merge<-grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],gene_list[[5]], gene_list[[6]],ncol=3)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/merge_six_",Group_tag,"_stat_barplot_expression_",tag,"_",genename,".pdf"),gene_list_merge,width=33, height=12)


#绘制基因活性的 umap distribution
gene_list<-list()
for ( genename in target_gene2){
  # genename <-"PTCH1"
  print(genename)
  #for gene expression
  express_merge_data3<-merge_data2[order(merge_data2[,genename],decreasing = T),]
  express_merge_target<- express_merge_data3[,c("UMAP_1","UMAP_2","Sex",genename)]
  colnames(express_merge_target)<-c("UMAP_1","UMAP_2","Sex","Expression")
  head(express_merge_target)
  express_merge_target$selected_cluster<-"other"
  express_merge_target[which(express_merge_target$Expression>0),]$selected_cluster<-"Pos"
  express_merge_target$selected_cluster <-factor(express_merge_target$selected_cluster,levels=c("Pos","other"))
  
  pos_cell_umap<-ggplot(data=express_merge_target[order(express_merge_target$selected_cluster,decreasing = T),],
                        mapping=aes(x=UMAP_1,y=UMAP_2,color = Expression))+
    geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
    scale_color_gradient(low="darkgrey",high = "red")+
    labs(title =paste0("Expression of ",genename))+theme_bw()+
    theme(panel.border = element_rect(colour="grey",fill=NA),
          panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
          axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=15),legend.position = "right")+facet_grid(.~Sex)
  pos_cell_umap
  ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_umap_",Group_tag,"_umap_expression_",tag,"_",genename,"_pos_value.png"),pos_cell_umap,width=10, height=6)
  ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_umap_",Group_tag,"_umap_expression_",tag,"_",genename,"_pos_value.pdf"),pos_cell_umap,width=10, height=6)
  gene_list<-c(gene_list,list(pos_cell_umap))
}
gene_list_merge<-grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],gene_list[[5]], gene_list[[6]],ncol=1)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/merge_six_",Group_tag,"_umap_expression_",tag,"_merge.png"),gene_list_merge,width=8, height=25)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/merge_six_",Group_tag,"_umap_expression_",tag,"_merge.pdf"),gene_list_merge,width=8, height=25)

gene_list_merge<-grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],gene_list[[5]], gene_list[[6]],ncol=2)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/merge_six_",Group_tag,"_umap_expression_",tag,"_merge2.png"),gene_list_merge,width=16, height=13)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/merge_six_",Group_tag,"_umap_expression_",tag,"_merge2.pdf"),gene_list_merge,width=16, height=13)
