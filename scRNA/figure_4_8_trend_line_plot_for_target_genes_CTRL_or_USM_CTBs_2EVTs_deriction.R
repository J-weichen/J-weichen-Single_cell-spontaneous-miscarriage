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


#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

###analysis for CTBs2EVTs branch
cds_branch<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2EVTs_selected_branch_monocle3_object.rds")
plot_cells(cds_branch, color_cells_by = "re_annotation_TFac", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)

##information substract monocle3
pseu_time<-data.frame(pseudotime(cds_branch))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)

##提取目标轨迹seurat3对象
Troph_target<-subset(x = target_final, cells =rownames(pseu_time))
DimPlot(Troph_target, group.by = "re_annotation_TFac",label = F,cols = ppCor_all2)
##calculation 
DefaultAssay(Troph_target) <- "RNA"
# Normalize RNA data for visualization purposes
Troph_target <- NormalizeData(Troph_target, verbose = FALSE)
#express_data <- as.matrix(GetAssayData(Troph_target, slot = "data"))
##提取表达矩阵
gene_signature<-c("MKI67","SERPINE1","TIMP3","STMN1","TIMP2","TIMP1","HLA-G","THBD","MALAT1","TGFB1","LAMA4")
gene_signature<-c("PAPPA2","CDH1","KRT19","KRT7","CYP19A1","EGFR","ERVFRD-1","RPL39","TIMP4","MMP11","MMP15","MMP17")

##final maker used
gene_signature<-c("TIMP3","SERPINE1","CDH1","STMN1","THBD")

plot_maker <- FeaturePlot(object = Troph_target, features = gene_signature,cols= c("grey", "purple"),ncol=3)
plot_maker
exprs <- data.frame(FetchData(object = Troph_target,slot = "data",vars = c(gene_signature,"Treat","re_annotation_TFac")))
exprs[1:4,]

##添加假时序信息
data_plot0<-merge(exprs,pseu_time,by=0)
Troph_group<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3")
data_plot0$re_annotation_TFac<- factor(data_plot0$re_annotation_TFac,levels = Troph_group)
data_plot0<-data_plot0[order(data_plot0$pseudotime,decreasing = F),]
data_plot0[1:4,]


##绘制时序拟合曲线
head(data_plot0)
data_plot0$re_annotation_TFac<-as.character(data_plot0$re_annotation_TFac)
data_plot0[which(!(data_plot0$re_annotation_TFac %in% c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3"))),]$re_annotation_TFac<-"other"
table(data_plot0$re_annotation_TFac)
data_plot0$re_annotation_TFac<-factor(data_plot0$re_annotation_TFac,levels = c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3","other"))

data_plot0<-data_plot_raw
data_plot0<-data_plot0[which(data_plot0$pseudotime<24),]

gene_list<-list()
for ( feature in gene_signature){
  #feature<-"TIMP3"
  print(feature)
  genename<-feature;method_type="loess"
  data_plot<-data_plot0[,c(genename,"Treat","re_annotation_TFac","pseudotime")]
  colnames(data_plot)<-c("expression","Treat","cell_names","pseudotime")
  trend_plot01<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_point(aes(colour=cell_names),shape=19, size = 1,alpha = 0.5) + 
    #  geom_point(colour="green", shape=19, size = 1, aes(fill = factor(Treat))) + 
    #  scale_fill_manual(values=c("blue", "cyan4"))+
    geom_rug(aes(colour = cell_names),sides="b") +  
    scale_colour_manual(values=ppCor)+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess#https://ggplot2.tidyverse.org/reference/geom_smooth.html
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  
  trend_plot01
  
  genename<-feature;method_type="gam"
  trend_plot02<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_point(aes(colour=cell_names),shape=19, size = 1,alpha = 0.5) + 
    #  geom_point(colour="green", shape=19, size = 1, aes(fill = factor(Treat))) + 
    #  scale_fill_manual(values=c("blue", "cyan4"))+
    geom_rug(aes(colour = cell_names),sides="b") +  
    scale_colour_manual(values=ppCor)+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  trend_plot02
  
  genename<-feature;method_type="loess"
  data_plot11<-data_plot0[,c(genename,"Treat","re_annotation_TFac","pseudotime")]
  colnames(data_plot)<-c("expression","Treat","cell_names","pseudotime")
  trend_plot11<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_point(aes(colour=Treat),shape=19, size = 1,alpha = 0.5) + 
    #  geom_point(colour="green", shape=19, size = 1, aes(fill = factor(Treat))) + 
    #  scale_fill_manual(values=c("blue", "cyan4"))+
    geom_rug(aes(colour = Treat),sides="b") +  
    scale_colour_manual(values=ppCor)+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess#https://ggplot2.tidyverse.org/reference/geom_smooth.html
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  
  trend_plot11
  
  genename<-feature;method_type="gam"
  trend_plot12<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_point(aes(colour=Treat),shape=19, size = 1,alpha = 0.5) + 
    #  geom_point(colour="green", shape=19, size = 1, aes(fill = factor(Treat))) + 
    #  scale_fill_manual(values=c("blue", "cyan4"))+
    geom_rug(aes(colour = Treat),sides="b") +  
    scale_colour_manual(values=ppCor)+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  trend_plot12
  
  ###for cell only in pos expression
  genename<-feature;method_type="loess"
  data_plot<-data_plot0[,c(genename,"Treat","re_annotation_TFac","pseudotime")]
  colnames(data_plot)<-c("expression","Treat","cell_names","pseudotime")
  head(data_plot)
  data_plot<-data_plot[which(data_plot$expression>0),]
  dim(data_plot0);dim(data_plot)#18330    #7709 
  trend_plot21<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_point(aes(colour=cell_names),shape=19, size = 1,alpha = 0.5) + 
    #  geom_point(colour="green", shape=19, size = 1, aes(fill = factor(Treat))) + 
    #  scale_fill_manual(values=c("blue", "cyan4"))+
    geom_rug(aes(colour = cell_names),sides="b") +  
    scale_colour_manual(values=ppCor)+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess#https://ggplot2.tidyverse.org/reference/geom_smooth.html
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  
  trend_plot21
  
  genename<-feature;method_type="gam"
  trend_plot22<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_point(aes(colour=cell_names),shape=19, size = 1,alpha = 0.5) + 
    #  geom_point(colour="green", shape=19, size = 1, aes(fill = factor(Treat))) + 
    #  scale_fill_manual(values=c("blue", "cyan4"))+
    geom_rug(aes(colour = cell_names),sides="b") +  
    scale_colour_manual(values=ppCor)+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  trend_plot22
  
  genename<-feature;method_type="loess"
  trend_plot31<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_point(aes(colour=Treat),shape=19, size = 1,alpha = 0.5) + 
    #  geom_point(colour="green", shape=19, size = 1, aes(fill = factor(Treat))) + 
    #  scale_fill_manual(values=c("blue", "cyan4"))+
    geom_rug(aes(colour = Treat),sides="b") +  
    scale_colour_manual(values=ppCor)+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess#https://ggplot2.tidyverse.org/reference/geom_smooth.html
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  
  trend_plot31
  
  genename<-feature;method_type="gam"
  trend_plot32<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_point(aes(colour=Treat),shape=19, size = 1,alpha = 0.5) + 
    #  geom_point(colour="green", shape=19, size = 1, aes(fill = factor(Treat))) + 
    #  scale_fill_manual(values=c("blue", "cyan4"))+
    geom_rug(aes(colour = Treat),sides="b") +  
    scale_colour_manual(values=ppCor)+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  trend_plot32
  
  plot_Cell<-ggarrange(trend_plot01,trend_plot02,trend_plot11,trend_plot12,trend_plot21,trend_plot22,trend_plot31, trend_plot32,ncol=4, nrow = 2)
  gene_list<-c(gene_list,list(plot_Cell))
}
gene_list[[1]]

####formal plot  绘制时序拟合曲线
library(moRandi)
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)

head(data_plot0)
data_plot0$re_annotation_TFac<-as.character(data_plot0$re_annotation_TFac)
data_plot0[which(!(data_plot0$re_annotation_TFac %in% c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3"))),]$re_annotation_TFac<-"other"
table(data_plot0$re_annotation_TFac)
data_plot0$re_annotation_TFac<-factor(data_plot0$re_annotation_TFac,levels = c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3","other"))

data_plot0<-data_plot0[which(data_plot0$pseudotime<24),]
symmary_data<-melt(data_plot0, id=c("Treat","re_annotation_TFac","pseudotime","Row.names"),variable.name="gene",value.name="expression",na.rm = TRUE)
head(symmary_data)
symmary_data<-symmary_data[which(symmary_data$expression>0),]
table(symmary_data$gene,symmary_data$Treat)

for ( feature in gene_signature){
  #feature<-"THBD"# 6394
  print(feature)
  ###for cell only in pos expression
  genename<-feature;method_type="loess"
  data_plot<-data_plot0[,c(genename,"Treat","re_annotation_TFac","pseudotime")]
  colnames(data_plot)<-c("expression","Treat","cell_names","pseudotime")
  head(data_plot)
  data_plot<-data_plot[which(data_plot$expression>0),]
  dim(data_plot0);dim(data_plot)# #TIMP3 7709  #SERPINE1:2330 #CDH1 2998 #STMN1 8991 #THBD 1494 
  trend_plot21<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
     geom_rug(aes(colour = cell_names),sides="b") +# ylim(0,1.8)+
    scale_colour_manual(values=c("red",my_morandi_colors[1:2],"blue",my_morandi_colors[6:8],my_morandi_colors[4]))+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess#https://ggplot2.tidyverse.org/reference/geom_smooth.html
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  
  trend_plot21
  
  genename<-feature;method_type="gam"
  trend_plot22<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_rug(aes(colour = cell_names),sides="b") +  
    scale_colour_manual(values=c("red",my_morandi_colors[1:2],"blue",my_morandi_colors[6:8],my_morandi_colors[4]))+#ylim(-2,7)+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,se=TRUE)+ #gam glm loess
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  trend_plot22
  ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/CTBs2EVTs_CTRL_USM_trend_line_loess_",genename,".pdf"),trend_plot21,width = 5, height =5)
  ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/CTBs2EVTs_CTRL_USM_trend_line_gam_",genename,".pdf"),trend_plot22,width = 5, height =5)
}

#         Abortion CTRL
#TIMP3        1144 5250
#SERPINE1      804 1526
#CDH1          293 2705
#STMN1        1099 7892
#THBD          149 1345

##final used file
for ( feature in gene_signature){
  #feature<-"CDH1"# 6394
  print(feature)
  genename<-feature;method_type="gam"
  data_plot<-data_plot0[,c(genename,"Treat","re_annotation_TFac","pseudotime")]
  colnames(data_plot)<-c("expression","Treat","cell_names","pseudotime")
  head(data_plot)
  data_plot<-data_plot[which(data_plot$expression>0),]
  trend_plot22<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_rug(aes(colour = cell_names),sides="b") +
    scale_colour_manual(values=c("red",my_morandi_colors[1:2],"blue",my_morandi_colors[6:8],my_morandi_colors[4]))+
   # scale_y_continuous(limits=c(0,3),breaks=c(0,0.3,0.6,0.9,1.2,1.5,1.8))+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,formula = y ~ s(x, bs = "cs"),se=TRUE)+ #gam glm loess
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5),#axis.text=element_blank(),
                      panel.grid=element_blank(),axis.ticks.x = element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  trend_plot22
  ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/CTBs2EVTs_CTRL_USM_trend_line_gam_",genename,".pdf"),trend_plot22,width = 5, height =5)
}

for ( feature in gene_signature){
  #feature<-"CDH1"# 6394
  print(feature)
  genename<-feature;method_type="loess"
  data_plot<-data_plot0[,c(genename,"Treat","re_annotation_TFac","pseudotime")]
  colnames(data_plot)<-c("expression","Treat","cell_names","pseudotime")
  head(data_plot)
  data_plot<-data_plot[which(data_plot$expression>0),]
  trend_plot22<-ggplot(data_plot,aes(x= pseudotime,y=expression,colour=Treat))+
    geom_rug(aes(colour = cell_names),sides="b") +
    scale_colour_manual(values=c("red",my_morandi_colors[1:2],"blue",my_morandi_colors[6:8],my_morandi_colors[4]))+
    # scale_y_continuous(limits=c(0,3),breaks=c(0,0.3,0.6,0.9,1.2,1.5,1.8))+
    geom_smooth(mapping = aes(x= pseudotime,y=expression,colour=Treat),data = data_plot,method=method_type,formula = y ~ s(x, bs = "cs"),se=TRUE)+ #gam glm loess
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5), panel.grid=element_blank())+
    labs(x = "pseudotime", y = "LogNormalized count", title =paste0(genename,":",method_type))
  trend_plot22
  ggsave(file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/CTBs2EVTs_CTRL_USM_trend_line_loess_",genename,".pdf"),trend_plot22,width = 5, height =5)
}
