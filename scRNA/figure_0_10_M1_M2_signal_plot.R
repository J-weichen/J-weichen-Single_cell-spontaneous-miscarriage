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
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
target_new<-subset(x = target_final, subset = re_annotation_TFac %in% c("Mycs","HCs"))
target_new$major_cluster_Treat <- paste(target_new$re_annotation_TFac, target_new$Treat, sep = "_")

DefaultAssay(target_new) <- "RNA"
target_new <- SCTransform(object = target_new,verbose = FALSE) 
#DefaultAssay(target_new) <- "RNA"
#target_new <- NormalizeData(target_new, verbose = FALSE)

#read Y special maker ##maker come from "Transcriptional Profiling of the Human Monocyte-to-Macrophage Differentiation and Polarization: New Molecules and Patterns of Gene Expression"
Mac_special_makers <- read.csv("/mnt/data/chenwei/jianghai/0.script/M1_M2_signaling.csv",header = T)
head(Mac_special_makers)
M1_signal0<-unique(as.character(Mac_special_makers[which(Mac_special_makers$Signal == "M1_signal"),]$Gene_Symbol))
M2_signal0<-unique(as.character(Mac_special_makers[which(Mac_special_makers$Signal == "M2_signal"),]$Gene_Symbol))

#M1_signal <-intersect(M1_signal0,rownames(target_new[["RNA"]]))
#M2_signal <-intersect(M2_signal0,rownames(target_new[["RNA"]]))
#target_new$M1_signalsum<-colSums(x = GetAssayData(object = target_new,assay = "RNA", slot = "counts")[M1_signal, , drop = FALSE])
#target_new$M2_signalsum<-colSums(x = GetAssayData(object = target_new,assay = "RNA", slot = "counts")[M2_signal, , drop = FALSE])
M1_signal <-intersect(M1_signal0,rownames(target_new[["SCT"]]))
M2_signal <-intersect(M2_signal0,rownames(target_new[["SCT"]]))
target_new$M1_signalsum<-colSums(x = GetAssayData(object = target_new,assay = "SCT", slot = "counts")[M1_signal, , drop = FALSE])
target_new$M2_signalsum<-colSums(x = GetAssayData(object = target_new,assay = "SCT", slot = "counts")[M2_signal, , drop = FALSE])

#violin ploy in manuscript
exprs <- data.frame(FetchData(object = target_new,slot = "data",vars = c("M1_signalsum","M2_signalsum","re_annotation_TFac","Treat")))
data_plot <- melt(exprs)
data_plot$log2_value<-log2(data_plot$value)
data_plot$major_cluster_Treat <- paste(data_plot$re_annotation_TFac, data_plot$Treat, sep = "_")

head(data_plot)

plot_MyCs_signal1 <-ggviolin(data_plot, x = "variable",y="log2_value",fill = "variable", palette = ppCor[c(3,4)], add = "boxplot", add.params = list(fill="white"))+
  labs(y = "Expression level log2_normalized_count", x = "re_annotation_TFac")+ylim(0,max(data_plot$log2_value)+1)+
  facet_wrap(~major_cluster_Treat,scales="free",nrow=2)+
  stat_compare_means(aes(group = variable,method = "wilcox.test", label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.5,label.y =10)

plot_Cell_exp<-ggplot(data_plot,aes(x = major_cluster_Treat,y=log2_value,fill = variable)) + 
  #geom_violin(scale ="area")+# "width"
  geom_violin(scale ="width")+# "width"
  theme_bw()+scale_fill_manual(values=ppCor[c(3,4)])+
  labs(y = "Expression level log2_normalized_count", x = "main_Cell_type")+ylim(0,max(data_plot$log2_value)+1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x.bottom = element_blank())+
  theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
        axis.title = element_text(size=15),axis.text.x = element_text(face="bold",size=8,angle = 0,vjust =1,hjust=1), 
        axis.text.y = element_text(face="bold",size=15),legend.position = "right")+
  theme(strip.text=element_text(colour = "black", face="bold", size=7), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))
plot_MyCs_signal2 <-plot_Cell_exp+  stat_compare_means(aes(group = variable,method = "wilcox.test", label=paste0(..method..,"\n",..p.signif..)), label.x = 1.5,label.y =10)
plot_MyCs_signal3<-plot_MyCs_signal2+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=3, col="purple", position = position_dodge(0.9))

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1_M2_global_signal_Violin_expression_in_global_MyCs_HCs.pdf",plot_MyCs_signal1,width=12, height=12)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1_M2_global_signal_Violin_expression_in_global_MyCs_HCs_together.pdf",plot_MyCs_signal3,width=12, height=6)

target_new$M1_M2_polar<-target_new$M1_signalsum/target_new$M2_signalsum
exprs <- data.frame(FetchData(object = target_new,slot = "data",vars = c("M1_M2_polar","re_annotation_TFac","Treat")))
data_plot <- melt(exprs)
data_plot$major_cluster_Treat <- paste(data_plot$re_annotation_TFac, data_plot$Treat, sep = "_")
head(data_plot)
data_plot$major_cluster_Treat <- factor(data_plot$major_cluster_Treat,levels =c("HCs_Abortion","HCs_CTRL","Mycs_Abortion","Mycs_CTRL"))

comparison_list= list(c("HCs_Abortion", "HCs_CTRL"), c("Mycs_Abortion", "Mycs_CTRL"))
                   
plot_MyCs_signal4<-ggviolin(data_plot, x = "major_cluster_Treat",y="value",fill = "Treat", palette = ppCor,outlier.colour=NA)+#, add = "boxplot", add.params = list(fill="white")
  labs(y = "Foldchange_M1_M2", x = "main_Cell_type")+
  #  ylim(0,max(data_plot$value)+1)+
  ylim(0,5)+
  stat_compare_means(comparisons =comparison_list,aes(group = Treat,method = "wilcox.test", label=paste0(..method..,"\n",..p.signif..)), label.x = 1.5,label.y =3.5)
plot_MyCs_HCs_signal5<-plot_MyCs_signal4+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4,outlier.colour=NA) +
  stat_summary(fun=mean, geom="point", shape=18, size=3, col="orange", position = position_dodge(0.9))
plot_MyCs_HCs_signal5
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1_M2_polar_signal_Violin_expression_in_global_MyCs_HCs_together.pdf",plot_MyCs_HCs_signal5,width=8, height=7)

##compare for each genes
target_new<-subset(x = target_final, subset = re_annotation_TFac %in% c("Mycs","HCs"))
target_new$major_cluster_Treat <- paste(target_new$re_annotation_TFac, target_new$Treat, sep = "_")
table(target_new$major_cluster_Treat)
#HCs_Abortion      HCs_CTRL Mycs_Abortion     Mycs_CTRL 
#      1727           175          4174           400 
DefaultAssay(target_new) <- "RNA"
target_new <- NormalizeData(target_new, verbose = FALSE)

#reading final seurat object
express_data <- as.matrix(GetAssayData(target_new, slot = "data"))
Idents(object = target_new) <- "re_annotation_TFac"

#extract metadata information
metadata_Cell<-target_new@meta.data
final_umap<-data.frame(Embeddings(target_new[["umap"]]))
#final_umap$raw_id <- rownames(final_umap)
head(metadata_Cell);head(final_umap)
metadata_Cell_Umap_add <- merge(metadata_Cell,final_umap,by=0)
head(metadata_Cell_Umap_add)
rownames(metadata_Cell_Umap_add)<-metadata_Cell_Umap_add$Row.names

##plot for target genes
target_new$cell_class<-paste0(target_new$re_annotation_TFac,"_",target_new$Treat)
table(target_new$cell_class)
#HCs_Abortion      HCs_CTRL Mycs_Abortion     Mycs_CTRL 
#     1727           175          4174           400 
length(c(M1_signal,M2_signal))# 80
##get genes in target 
target_gene<-rev(c(M1_signal,M2_signal))
tag<-"Macrophage_HCs_polar_genes"
##正式绘图
Plotdot2 <-DotPlot(target_new, features = target_gene,group.by ="cell_class",cols =  c("blue","red"),assay = "RNA") #+ RotatedAxis()
Plotdot2
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

ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/",tag,"_genes_dotplot_all_celltype_split.pdf"),Plotdot31,width=5, height=10)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/",tag,"_genes_dotplot_all_celltype_split_2.pdf"),Plotdot32,width=5, height=10)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/",tag,"_genes_dotplot_all_celltype_split_3.pdf"),Plotdot02,width=15, height=5)

##plot violin with p value
target_gene2<-target_gene[which(target_gene %in% rownames(target_new))]
target_expression <- GetAssayData(object = target_new,assay = "RNA", slot = "data")[target_gene2,]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))

metadata_Cell_Umap_add2<-metadata_Cell_Umap_add[rownames(target_expression),]
merge_data2 <- merge(metadata_Cell_Umap_add2, target_expression,by=0,sort=FALSE)
head(merge_data2)

##绘制两比较组间的表达量
group_tag<-"M1_signal"
for ( genename in M1_signal){
  # genename <-"ADK" 
  print(genename)
  merge_data_target<- merge_data2[,c("re_annotation_TFac","Treat",genename)]
  colnames(merge_data_target)<-c("cell_type","Treat","Expression_level")
  plot_Cell_exp<- ggplot(merge_data_target,aes(x = cell_type,y=Expression_level,fill = Treat)) + 
    geom_violin(scale = "width") +ylim(0,max(merge_data_target$Expression_level)+0.5)+ #geom_jitter(size = 0.00)+coord_flip()+ facet_grid(. ~ variable,scale='free') +
    theme_bw()+scale_fill_manual(values=ppCor)+
    labs(y = "Expression level", x = "main_Cell_type", title = genename)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.y.left = element_blank())+
    theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),legend.position = "right")+
    theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))+
    stat_compare_means(method = "wilcox.test",data=merge_data_target,aes(x=cell_type,y = Expression_level,group = Treat,label= ..p.signif..),hide.ns = TRUE, label.x = 1.5,label.y =max((merge_data_target$Expression_level)))
  plot_Cell_exp2<- plot_Cell_exp+ geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=4, col="orange", position = position_dodge(0.9))
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/macrophage/stat_barplot_expression_",group_tag,"_",tag,"_",genename,".pdf"),plot_Cell_exp2,width=8, height=6)
}
group_tag<-"M2_signal"
for ( genename in M2_signal){
  # genename <-"ADK" 
  print(genename)
  merge_data_target<- merge_data2[,c("re_annotation_TFac","Treat",genename)]
  colnames(merge_data_target)<-c("cell_type","Treat","Expression_level")
  plot_Cell_exp<- ggplot(merge_data_target,aes(x = cell_type,y=Expression_level,fill = Treat)) + 
    geom_violin(scale = "width") +ylim(0,max(merge_data_target$Expression_level)+0.5)+ #geom_jitter(size = 0.00)+coord_flip()+ facet_grid(. ~ variable,scale='free') +
    theme_bw()+scale_fill_manual(values=ppCor)+
    labs(y = "Expression level", x = "main_Cell_type", title = genename)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.y.left = element_blank())+
    theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),legend.position = "right")+
    theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))+
    stat_compare_means(method = "wilcox.test",data=merge_data_target,aes(x=cell_type,y = Expression_level,group = Treat,label= ..p.signif..),hide.ns = TRUE, label.x = 1.5,label.y =max((merge_data_target$Expression_level)))
  plot_Cell_exp2<- plot_Cell_exp+ geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=4, col="orange", position = position_dodge(0.9))
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/macrophage/stat_barplot_expression_",group_tag,"_",tag,"_",genename,".pdf"),plot_Cell_exp2,width=6, height=6)
}

##以下未跑

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


##for selected genes
Group_tag<-"Macrophage"
length(Mac_special_makers)
plot1<-VlnPlot(target_final, features = gene_signature, group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], y.max =4, combine =F) 
plot1_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1[[1]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_2<-plot1[[2]]+ stat_compare_means(method = "wilcox.test",data=plot1[[2]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_3<-plot1[[3]]+ stat_compare_means(method = "wilcox.test",data=plot1[[3]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_4<-plot1[[4]]+ stat_compare_means(method = "wilcox.test",data=plot1[[4]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_5<-plot1[[5]]+ stat_compare_means(method = "wilcox.test",data=plot1[[5]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
plot1_6<-plot1[[6]]+ stat_compare_means(method = "wilcox.test",data=plot1[[6]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =3.5)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
gene_list_merge<-grid.arrange(plot1_1,plot1_2,ncol=1)
gene_list_merge<-grid.arrange(plot1_1,plot1_2,plot1_3,plot1_4,plot1_5,ncol=1)
ggsave("/mnt/data/chenwei/qinmen_BR/keshi/merge_selecect_four_regulated_genes_qinmeng_stat_violin_plot_expression_seurat.pdf",gene_list_merge,width=20, height=12)

gene_list_merge<-grid.arrange(plot1_1,plot1_2,plot1_3,plot1_4,plot1_5, plot1_6,ncol=3)
ggsave("/mnt/data/chenwei/qinmen_BR/keshi/merge_six_genes_qinmeng_stat_violin_plot_expression_seurat.pdf",gene_list_merge,width=20, height=12)
ggsave("/mnt/data/chenwei/qinmen_BR/keshi/merge_six_genes_qinmeng_stat_violin_plot_expression_seurat.png",gene_list_merge,width=20, height=12)

