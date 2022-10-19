rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(Seurat)
library(scales)
library(ggsci)
library(ggplot2)
library(reshape2)
library(ggpubr)

library("ggbeeswarm")
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
#set colors
#step1 ：read Seurat object and load initial coldata and matrix
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
target_new<-subset(x = target_final, subset = re_annotation_TFac == "Mycs")
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
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1_M2_global_signal_Violin_expression_in_global_MyCs.pdf",plot_MyCs_signal1,width=8, height=12)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1_M2_global_signal_Violin_expression_in_global_MyCs_together.pdf",plot_MyCs_signal3,width=10, height=6)


##four group camparesion
data_plot$major_cluster_Treat_type <- paste(data_plot$major_cluster_Treat, data_plot$variable, sep = "_")

head(data_plot)
plot_Cell_exp<-ggplot(data_plot,aes(x = major_cluster_Treat_type,y=log2_value,fill = variable)) + 
  geom_violin(scale ="width")+# "width"
  theme_bw()+scale_fill_manual(values=ppCor[c(3,4)])+
  labs(y = "Expression level log2_normalized_count", x = "main_Cell_type")+
  scale_y_continuous(limits=c(0,16),breaks = seq(0,16,2),minor_breaks = seq(0, 16,0.5),expand = c(0,0.5))+
  theme(panel.grid.minor.x  = element_blank(), 
        panel.grid.minor.y =element_line(linetype = 2,colour="grey",size = 0.2),
        axis.ticks.x.bottom = element_blank())+
  theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
        axis.title = element_text(size=15),axis.text.x = element_text(face="bold",size=8,angle = 45,vjust =1,hjust=1), 
        axis.text.y = element_text(face="bold",size=15),legend.position = "right")+
  theme(strip.text=element_text(colour = "black", face="bold", size=7), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))
comparison_list= list(c("Mycs_Abortion_M1_signalsum", "Mycs_Abortion_M2_signalsum"),
                      c("Mycs_CTRL_M1_signalsum", "Mycs_CTRL_M2_signalsum"), 
                      c("Mycs_Abortion_M1_signalsum", "Mycs_CTRL_M1_signalsum"),
                      c("Mycs_Abortion_M2_signalsum", "Mycs_CTRL_M2_signalsum"))
plot_MyCs_signal2 <-plot_Cell_exp+  stat_compare_means(comparisons =comparison_list,aes(group = variable,method = "wilcox.test", label=paste0(..method..,"\n",..p.signif..)))
plot_MyCs_signal3<-plot_MyCs_signal2+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=3, col="purple", position = position_dodge(0.9))
plot_MyCs_signal3

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1_M2_global_signal_Violin_expression_in_global_MyCs_together_multi_stat.pdf",plot_MyCs_signal3,width=8, height=8)

##for M1/M2 ratio
target_new$M1_M2_polar<-target_new$M1_signalsum/target_new$M2_signalsum
exprs <- data.frame(FetchData(object = target_new,slot = "data",vars = c("M1_M2_polar","re_annotation_TFac","Treat")))
data_plot <- melt(exprs)
data_plot$major_cluster_Treat <- paste(data_plot$re_annotation_TFac, data_plot$Treat, sep = "_")
head(data_plot)
data_plot$major_cluster_Treat <- factor(data_plot$major_cluster_Treat,levels =c("Mycs_Abortion","Mycs_CTRL"))

plot_MyCs_signal4<-ggviolin(data_plot, x = "major_cluster_Treat",y="value",fill = "Treat", palette = ppCor,outlier.colour=NA)+#, add = "boxplot", add.params = list(fill="white")
  labs(y = "Foldchange_M1_M2", x = "main_Cell_type")+
#  ylim(0,max(data_plot$value)+1)+
  ylim(0,5)+
  stat_compare_means(aes(group = Treat,method = "wilcox.test", label=paste0(..method..,"\n",..p.signif..)), label.x = 1.5,label.y =3.5)
plot_MyCs_signal5<-plot_MyCs_signal4+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4,outlier.colour=NA) +
  stat_summary(fun=mean, geom="point", shape=18, size=3, col="orange", position = position_dodge(0.9))
plot_MyCs_signal5
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/M1_M2_polar_signal_Violin_expression_in_global_MyCs_together.pdf",plot_MyCs_signal5,width=6, height=6)

##绘制两比较组间的表达量
target_new<-subset(x = target_final, subset = re_annotation_TFac =="Mycs")
table(target_new$Treat)
#Abortion     CTRL 
#  4174      400 
target_new$major_cluster_Treat <- paste(target_new$re_annotation_TFac, target_new$Treat, sep = "_")

#extract metadata information
metadata_Cell<-target_new@meta.data
final_umap<-data.frame(Embeddings(target_new[["umap"]]))
#final_umap$raw_id <- rownames(final_umap)
head(metadata_Cell);head(final_umap)
metadata_Cell_Umap_add <- merge(metadata_Cell,final_umap,by=0)
head(metadata_Cell_Umap_add)
rownames(metadata_Cell_Umap_add)<-metadata_Cell_Umap_add$Row.names
#normalization
DefaultAssay(target_new) <- "RNA"
target_new <- NormalizeData(target_new, verbose = FALSE)

##plot for target genes
target_new$cell_class<-paste0(target_new$re_annotation_TFac,"_",target_new$Treat)
table(target_new$cell_class)
#Mycs_Abortion     Mycs_CTRL 
#     4174           400 
length(c(M1_signal,M2_signal))# 80
##get genes in target 
target_gene<-rev(c(M1_signal,M2_signal))
tag<-"Macrophage_polar_genes"

##plot violin with p value
target_gene2<-target_gene[which(target_gene %in% rownames(target_new))]
target_expression <- GetAssayData(object = target_new,assay = "RNA", slot = "data")[target_gene2,]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))

metadata_Cell_Umap_add2<-metadata_Cell_Umap_add[rownames(target_expression),]
merge_data2 <- merge(metadata_Cell_Umap_add2, target_expression,by=0,sort=FALSE)
head(merge_data2)
group_tag<-"M1_signal"
for ( genename in M1_signal){
  # genename <-"ADK" 
  print(genename)
  merge_data_target<- merge_data2[,c("re_annotation_TFac","Treat",genename)]
  colnames(merge_data_target)<-c("cell_type","Treat","Expression_level")
  plot_Cell_exp<- ggplot(merge_data_target,aes(x = Treat,y=Expression_level,fill = Treat)) + 
    geom_violin(scale = "width") +ylim(0,max(merge_data_target$Expression_level)+0.5)+ #geom_jitter(size = 0.00)+coord_flip()+ facet_grid(. ~ variable,scale='free') +
    theme_bw()+scale_fill_manual(values=ppCor)+
    labs(y = "Expression level", x = "main_Cell_type", title = genename)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.y.left = element_blank())+
    theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),legend.position = "right")+
    theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))+
    stat_compare_means(method = "wilcox.test",data=merge_data_target,aes(x=Treat,y = Expression_level,group = Treat,label= ..p.signif..),hide.ns = TRUE, label.x = 1.5,label.y =max((merge_data_target$Expression_level)))
  plot_Cell_exp2<- plot_Cell_exp+ geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=4, col="orange", position = position_dodge(0.9))
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/macrophage/stat_barplot_expression_",group_tag,"_",tag,"_",genename,".pdf"),plot_Cell_exp2,width=6, height=6)
}
group_tag<-"M2_signal"
for ( genename in M2_signal){
  # genename <-"ADK" 
  print(genename)
  merge_data_target<- merge_data2[,c("re_annotation_TFac","Treat",genename)]
  colnames(merge_data_target)<-c("cell_type","Treat","Expression_level")
  plot_Cell_exp<- ggplot(merge_data_target,aes(x = Treat,y=Expression_level,fill = Treat)) + 
    geom_violin(scale = "width") +ylim(0,max(merge_data_target$Expression_level)+0.5)+ #geom_jitter(size = 0.00)+coord_flip()+ facet_grid(. ~ variable,scale='free') +
    theme_bw()+scale_fill_manual(values=ppCor)+
    labs(y = "Expression level", x = "main_Cell_type", title = genename)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.y.left = element_blank())+
    theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),legend.position = "right")+
    theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))+
    stat_compare_means(method = "wilcox.test",data=merge_data_target,aes(x=Treat,y = Expression_level,group = Treat,label= ..p.signif..),hide.ns = TRUE, label.x = 1.5,label.y =max((merge_data_target$Expression_level)))
  plot_Cell_exp2<- plot_Cell_exp+ geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=4, col="orange", position = position_dodge(0.9))
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/macrophage/stat_barplot_expression_",group_tag,"_",tag,"_",genename,".pdf"),plot_Cell_exp2,width=6, height=6)
}


FeaturePlot(object = target_new, features = c("M1_signalsum","M2_signalsum","M1_M2_polar"),split.by = "Treat", cols= c("grey", "purple"))
