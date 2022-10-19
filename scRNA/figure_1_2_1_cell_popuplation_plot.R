rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(Seurat)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(grid)
library(gridExtra)
library(cowplot)
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:46))]
show_col(ppCor_all2)

#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

plot_Cell0<-DimPlot(target_final, group.by = "re_anno_TFac_major",cols =ppCor_all2[c(1:9,11)])#+ NoLegend()
plot_Cell1<-DimPlot(target_final, group.by = "re_annotation_TFac",cols =ppCor_all2)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Umap_anno_plot_label.pdf",width = 10,height = 10)
LabelClusters(plot = plot_Cell0, id = 're_anno_TFac_major',segment.color ="grey",colour = "black",fontface = "bold",size = 4)+ NoLegend()
LabelClusters(plot = plot_Cell1, id = 're_annotation_TFac',segment.color ="grey",colour = "black",fontface = "bold",size = 4)+ NoLegend()
dev.off()

#绘图
##plot cell number sperated by Abortion and CTRL
cell_prop<-as.data.frame(table(target_final$re_annotation_TFac, target_final$Treat))
cell_prop_raw1<-cell_prop[which(cell_prop$Freq>0),]

cell_prop_raw2 <- dcast(cell_prop_raw1,Var1 ~ Var2, value.var = "Freq")
cell_prop_raw2[is.na(cell_prop_raw2)] <- 0
cell_prop_raw2$sum<-cell_prop_raw2$Abortion+cell_prop_raw2$CTRL

colnames(cell_prop_raw1)<-c( "Cell_type","Group","Freq")
head(cell_prop_raw1);dim(cell_prop_raw1)##34 
cell_prop_raw1$Cell_type<-factor(cell_prop_raw1$Cell_type,levels=subgroup_order0)

plot_cell_prop_num_decrease<-ggplot(data=cell_prop_raw1, mapping=aes(x= Cell_type,y=Freq,fill=Group))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+scale_fill_manual(name="Types",values=my_morandi_colors[c(1,21)])+
  theme_classic()+labs(x="Cell types",y="Cell number",title="The Cell number for each cells")+
  scale_y_continuous(breaks = seq(0,15000,by=1000))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_cell_prop_num_decrease2<-plot_cell_prop_num_decrease+geom_text(aes(label=Freq),size=3,position=position_stack(vjust=0.5))#3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
plot_cell_prop_num_decrease2
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Barplot_cell_number_for_each_sig_HGt_sub_celltype_group_highlight1.pdf", plot_cell_prop_num_decrease2,width=8, height=8)

##order by total cell number
head(cell_prop_raw2)
cell_name<-as.character(cell_prop_raw2[order(cell_prop_raw2$sum,decreasing = T),"Var1"])
plot_cell_prop_num<-cell_prop_raw1[which(cell_prop_raw1$Cell_type %in% cell_name),]
plot_cell_prop_num$Cell_type<-factor(plot_cell_prop_num$Cell_type,levels=cell_name)
plot_cell_prop_num_decrease<-ggplot(data=plot_cell_prop_num, mapping=aes(x= Cell_type,y=Freq,fill=Group))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+scale_fill_manual(name="Types",values=my_morandi_colors[c(1,21)])+
  theme_classic()+labs(x="Cell types",y="Cell number",title="The Cell number for each cells")+
  scale_y_continuous(breaks = seq(0,15000,by=1000))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 8,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_cell_prop_num_decrease2<-plot_cell_prop_num_decrease+geom_text(aes(label=Freq),size=4,position=position_stack(vjust=0.5))#3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
plot_cell_prop_num_decrease2
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2c_Barplot_cell_number_for_each_sig_HGt_sub_celltype_group_highlight1.pdf", plot_cell_prop_num_decrease2,width=8, height=8)

##plot supplement figure 2a
P1<-DimPlot(object = target_final, reduction = "umap", group.by = "re_anno_TFac_major",cols =my_morandi_colors[1:10])
P2<-DimPlot(object = target_final, reduction = "umap", group.by = "re_annotation_TFac",cols =my_morandi_colors)#+guides(fill = guide_legend(ncols =2))
P3<-DimPlot(object = target_final, reduction = "umap",group.by = "Phase",cols = my_morandi_colors[c(8,3,14)])
P4<-DimPlot(object = target_final, reduction = "umap",  group.by = "Treat",cols =my_morandi_colors[c(1,21)])
P5<-DimPlot(object = target_final, reduction = "umap",group.by = "sample",cols = my_morandi_colors)
P6<-DimPlot(object = target_final, reduction = "umap",group.by = "sample_code",cols = my_morandi_colors)
type_plot<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/all_cell_type_umap.png",type_plot,width=22, height=12)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/all_cell_type_umap_re_anno_TFac_major.png",P1,width=10, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/all_cell_type_umap_re_annotation_TFac.png",P2,width=9, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/all_cell_type_umap_Phase.png",P3,width=9, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/all_cell_type_umap_Treat.png",P4,width=9, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/all_cell_type_umap_sample.png",P5,width=9, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/all_cell_type_umap_sample_code.png",P6,width=9, height=8)

##plot supplement figure 2a
type_plot<-CombinePlots(plots = list(P1,P3,P4,P6),ncol =2,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2a_all_cell_type_umap.png",type_plot,width=22, height=12)
P12<-P1+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
P32<-P3+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
P42<-P4+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
P62<-P6+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2a_major_cell_type_umap.png",P12,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2a_Phase_umap.png",P32,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2a_Treat_umap.png",P42,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2a_sample_code_umap.png",P62,width=8, height=8)

#cell distribution by each type split
P21<-DimPlot(object = target_final, reduction = "umap", group.by = "Phase",split.by="Phase",cols =  my_morandi_colors[c(8,3,14)],ncol=3)
P22<-DimPlot(object = target_final, reduction = "umap",group.by = "Treat",split.by="Treat",cols = my_morandi_colors[c(1,21)],ncol=2)
P23<-DimPlot(object = target_final, reduction = "umap",group.by = "sample_code",split.by="sample_code",cols = my_morandi_colors,ncol=3)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Umap_Phase_split_cell_type_umap.png", P21,width=25, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Umap_Treat_split_cell_type_umap.png", P22,width=17, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Umap_sample_split_cell_type_umap.png", P23,width=25, height=25)

P212<-P21+NoLegend()+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
P222<-P22+NoLegend()+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
P232<-P23+NoLegend()+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Umap_Phase_split_cell_type_umap_notype.png", P212,width=24, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Umap_Treat_split_cell_type_umap_notype.png", P222,width=16, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Umap_sample_split_cell_type_umap_notype.png", P232,width=12, height=12)


P31<-DimPlot(object = target_final, group.by = "re_annotation_TFac",cols =my_morandi_colors,split.by="Treat",ncol=2)
P32<-DimPlot(object = target_final, group.by = "re_anno_TFac_major",cols =my_morandi_colors,split.by="Treat",ncol=2)#ppCor_all2[c(1:9,11)]
P33<-DimPlot(object = target_final,  group.by = "Phase",cols = my_morandi_colors[c(8,3,14)],split.by="Treat",ncol=2)
P34<-DimPlot(object = target_final,  group.by = "sample_code",cols = my_morandi_colors,split.by="Treat",ncol=2)
combine_plot<-P31/P32/P33/P34
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Umap_different_stage_treat_split_cell_type.png", combine_plot,width=17, height=32)

P312<-P31+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                      plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                      axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                      legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Umap_cell_subtypes_anno_umap_notype.png", P312,width=12, height=6)


#propetion bar sperated without Ery and mast
#target_final<-subset(x = target_final,subset = re_anno_TFac_major %in% c("Erythrocyte_HBA1","Mast_Cell_MS4A2"),invert = TRUE)
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

#计算每组细胞数目
table(target_final$re_anno_TFac_major)
table(target_final$re_annotation_TFac,target_final$sample_code)
as.data.frame(table(target_final$re_annotation_TFac,target_final$Treat))
as.data.frame(table(target_final$Treat,target_final$re_annotation_TFac))

major_subgroup_pieplot<-ggplot(as.data.frame(table(target_final$re_annotation_TFac)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors)+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank()) 
major_group_pieplot<-ggplot(as.data.frame(table(target_final$re_anno_TFac_major)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors[1:10])+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  
sample_pieplot<-ggplot(as.data.frame(table(target_final$sample_code)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors)+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  
Phase_pieplot<-ggplot(as.data.frame(table(target_final$Phase)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors[c(8,3,14)])+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  
group_pieplot<-ggplot(as.data.frame(table(target_final$Treat)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors[c(1,21)])+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Peiplot_Type_ratio_sub_celltype.pdf", major_subgroup_pieplot,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Peiplot_Type_ratio_major_celltype.pdf", major_group_pieplot,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Peiplot_Type_ratio_sample.pdf", sample_pieplot,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Peiplot_Type_ratio_phase.pdf", Phase_pieplot,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Peiplot_Type_ratio_group.pdf", group_pieplot,width=8, height=8)

#各细胞类群中各样本占比
cell.prop1<-as.data.frame(prop.table(table(target_final$re_annotation_TFac, target_final$sample_code)))
colnames(cell.prop1)<-c("Cell_type","Cell_origin","proportion")
Ratio_1<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Cell_origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors)+guides(fill=guide_legend(title=NULL))
Ratio_1
##各细胞类群中不同状态占比
cell.prop1<-as.data.frame(prop.table(table(target_final$re_annotation_TFac, target_final$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
Ratio_2<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Group))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors[c(1,21)])+guides(fill=guide_legend(title=NULL))
Ratio_2
#各细胞类群中phase比例
cell.prop1<-as.data.frame(prop.table(table(target_final$re_annotation_TFac, target_final$Phase)))
colnames(cell.prop1)<-c("Cell_type","Phase","proportion")
Ratio_3<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Phase))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors[c(8,3,14)])+guides(fill=guide_legend(title=NULL))
Ratio_3
cell_sample_combineplot<-CombinePlots(plots = list(Ratio_1,Ratio_2,Ratio_3),ncol = 1,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2cd_Barplot_sample_ratio_for_each_subcell_type.pdf", cell_sample_combineplot,width=10, height=10)

#各细胞大类群中各样本占比
cell.prop1<-as.data.frame(prop.table(table(target_final$re_anno_TFac_major, target_final$sample_code)))
colnames(cell.prop1)<-c("Cell_type","Cell_origin","proportion")
Ratio_1<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Cell_origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=30,hjust=1, vjust=1,size=6), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors[1:10])+guides(fill=guide_legend(title=NULL))
Ratio_1
##各细胞大类群中不同状态占比
cell.prop1<-as.data.frame(prop.table(table(target_final$re_anno_TFac_major, target_final$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
Ratio_2<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Group))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=30,hjust=1, vjust=1,size=6), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors[c(1,21)])+guides(fill=guide_legend(title=NULL))
Ratio_2

#各细胞大类群中年龄比例
cell.prop1<-as.data.frame(prop.table(table(target_final$re_anno_TFac_major, target_final$Phase)))
colnames(cell.prop1)<-c("Cell_type","Phase","proportion")
Ratio_3<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Phase))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=30,hjust=1, vjust=1,size=6), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors[c(8,3,14)])+guides(fill=guide_legend(title=NULL))
Ratio_3
cell_sample_combineplot<-CombinePlots(plots = list(Ratio_1,Ratio_2,Ratio_3),ncol = 1,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Barplot_sample_ratio_for_each_major_cell_type.pdf", cell_sample_combineplot,width=10, height=10)

#propetion bar sperated without Ery and mast
target_final_noEry_mast<-subset(x = target_final,subset = re_anno_TFac_major %in% c("Erythrocyte_HBA1","Mast_Cell_MS4A2"),invert = TRUE)
major_order1<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A")
subgroup_order1<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs")
target_final_noEry_mast$re_annotation_TFac<- factor(target_final_noEry_mast$re_annotation_TFac, levels=subgroup_order1,ordered=TRUE)
target_final_noEry_mast$re_anno_TFac_major<- factor(target_final_noEry_mast$re_anno_TFac_major,levels=major_order1,ordered=TRUE)

#计算去除红细胞后计算每组细胞数目
table(target_final_noEry_mast$re_anno_TFac_major)
table(target_final_noEry_mast$re_annotation_TFac,target_final_noEry_mast$sample_code)
as.data.frame(table(target_final_noEry_mast$re_annotation_TFac,target_final_noEry_mast$Treat))
as.data.frame(table(target_final_noEry_mast$Treat,target_final_noEry_mast$re_annotation_TFac))

#pieplot for each type and add frequence
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry_mast$re_annotation_TFac, target_final_noEry_mast$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
head(cell.prop1)
pieplot_ring1<-ggplot(cell.prop1,aes(x=Group,y=proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+coord_polar("y",start=1) + 
  theme_minimal()+ labs(title = "USM versus CTRL")+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=my_morandi_colors)#+facet_wrap(~Disease_group)
pieplot_ring1

cell.prop1<-as.data.frame(prop.table(table(target_final_noEry_mast$re_anno_TFac_major, target_final_noEry_mast$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
head(cell.prop1)
pieplot_ring2<-ggplot(cell.prop1,aes(x=Group,y=proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+coord_polar("y",start=1) + 
  theme_minimal()+ labs(title = "USM vs CTRL")+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=my_morandi_colors[1:10])#+facet_wrap(~Disease_group)
pieplot_ring2

cell.prop1<-as.data.frame(prop.table(table(target_final_noEry_mast$Phase, target_final_noEry_mast$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
head(cell.prop1)
pieplot_ring3<-ggplot(cell.prop1,aes(x=Group,y=proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+coord_polar("y",start=1) + 
  theme_minimal()+ labs(title = "USM vs CTRL")+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=my_morandi_colors[c(8,3,14)])#+facet_wrap(~Disease_group)
pieplot_ring3

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2f_Peiplot_Type_ratio_treat_split_subgroup_type.pdf", pieplot_ring1,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2f_Peiplot_Type_ratio_treat_split_major_type.pdf", pieplot_ring2,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF2f_Peiplot_Type_ratio_treat_split_phase.pdf", pieplot_ring3,width=8, height=8)

major_subgroup_pieplot<-ggplot(as.data.frame(table(target_final_noEry_mast$re_annotation_TFac)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors)+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank()) 
major_group_pieplot<-ggplot(as.data.frame(table(target_final_noEry_mast$re_anno_TFac_major)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors[1:10])+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  
sample_pieplot<-ggplot(as.data.frame(table(target_final_noEry_mast$sample_code)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors)+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  
Phase_pieplot<-ggplot(as.data.frame(table(target_final_noEry_mast$Phase)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors[c(8,3,14)])+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  
group_pieplot<-ggplot(as.data.frame(table(target_final_noEry_mast$Treat)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=my_morandi_colors[c(1,21)])+
  theme(panel.grid=element_blank(),axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/no_Ery_Peiplot_Type_ratio_sub_celltype.pdf", major_subgroup_pieplot,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/no_Ery_Peiplot_Type_ratio_major_celltype.pdf", major_group_pieplot,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/no_Ery_Peiplot_Type_ratio_sample.pdf", sample_pieplot,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/no_Ery_Peiplot_Type_ratio_phase.pdf", Phase_pieplot,width=8, height=8)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/no_Ery_Peiplot_Type_ratio_group.pdf", group_pieplot,width=8, height=8)

##细胞占比的超几何检验
cell_prop<-as.data.frame(table(target_final_noEry_mast$re_annotation_TFac, target_final_noEry_mast$Treat))
cell_prop_raw1<-cell_prop[which(cell_prop$Freq>0),]
cell_prop_raw2 <- dcast(cell_prop_raw1,Var1  ~  Var2   , value.var = "Freq")
cell_prop_raw2[is.na(cell_prop_raw2)] <- 0
Abortion_sum<-sum(cell_prop[which(cell_prop$Var2=="Abortion"),]$Freq)
CTRL_sum<-sum(cell_prop[which(cell_prop$Var2=="CTRL"),]$Freq)
cell_prop_raw2$Abortion_other<-Abortion_sum-cell_prop_raw2$Abortion
cell_prop_raw2$CTRL_other<-CTRL_sum-cell_prop_raw2$CTRL
head(cell_prop_raw2)
#卡方检验和fisher精确性检验细胞比例变动：https://www.jianshu.com/p/cf2c8bdcbeae
Element_percentage<-c();cnane_region<-c()
for ( row_code in c(1: nrow(cell_prop_raw2))){
  # row_code<-13
  Cell_name<-as.character(cell_prop_raw2[row_code,]$Var1)
  print(Cell_name)
  como<-cell_prop_raw2[row_code,c("CTRL","CTRL_other","Abortion","Abortion_other")]
  como2 <- as.table(cbind(c(como$CTRL,como$Abortion),c(como$CTRL_other,como$Abortion_other)))
  dimnames(como2) <- list(c('CTRL','Abortion'),c('Cell','Other_cell'))
  como3<-chisq.test(como2) 
  p_HG <- ifelse(min(como3$expected)>=5, round(chisq.test(como2,correct = F)$p.value,10),round(fisher.test(como2)$p.value,10))
  Element_percentage<-c(Element_percentage,p_HG)
  cnane_region<-c(cnane_region,Cell_name)
}

length(Element_percentage);length(cnane_region)#17
pvalue_data<-data.frame(Cell_type=cnane_region,p_HG=Element_percentage)
pvalue_data$q_adj_HG<-p.adjust(pvalue_data$p_HG,"BH")
colnames(cell_prop_raw2)<-c("Cell_type","Abortion","CTRL","Abortion_other","CTRL_other")
cell_prop_raw3<-merge(cell_prop_raw2,pvalue_data);head(cell_prop_raw3)
length(which(cell_prop_raw3$p_HG <=0.05));length(which(cell_prop_raw3$q_adj_HG<=0.05))
head(cell_prop_raw3)

#绘图
##plot cell number sperated by Abortion and CTRL
head(cell_prop_raw1);head(cell_prop_raw3)
cell_prop_raw3$sum<-cell_prop_raw3$Abortion+cell_prop_raw3$CTRL
colnames(cell_prop_raw1)<-c( "Cell_type","Group","Freq")
head(cell_prop_raw1);dim(cell_prop_raw1)##34 
cell_name<-as.character(cell_prop_raw3[order(cell_prop_raw3$sum,decreasing = T),"Cell_type"])

plot_cell_prop_num<-cell_prop_raw1[which(cell_prop_raw1$Cell_type %in% cell_name),]
plot_cell_prop_num$Cell_type<-factor(plot_cell_prop_num$Cell_type,levels=cell_name)

plot_cell_prop_num_decrease<-ggplot(data=plot_cell_prop_num, mapping=aes(x= Cell_type,y=Freq,fill=Group))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+scale_fill_manual(name="Types",values=ppCor_all[1:2])+
  theme_classic()+labs(x="Cell types",y="Cell number",title="The Cell number for each cells")+
  scale_y_continuous(breaks = seq(0,15000,by=1000))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 8,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_cell_prop_num_decrease
plot_cell_prop_num_decrease2<-plot_cell_prop_num_decrease+geom_text(aes(label=Freq),position=position_stack(vjust=0.5))#3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/no_Ery_Barplot_cell_number_for_each_sig_HGt_sub_celltype_group_highlight1.pdf", plot_cell_prop_num_decrease2,width=8, height=8)

##Barplot showing proportion Seperated by Abortion or CTRL 
plot_cell_prop_num2<-merge(cell_prop_raw1,cell_prop_raw3[, c("Cell_type","sum","p_HG","q_adj_HG")])
head(plot_cell_prop_num2)
plot_cell_prop_num2$Cell_type<-factor(plot_cell_prop_num2$Cell_type,levels=rev(subgroup_order1))
plot_cell_prop_num2$proportion<-ifelse(grepl("CTRL",plot_cell_prop_num2$Group),plot_cell_prop_num2$Freq*100/CTRL_sum,plot_cell_prop_num2$Freq*-100/Abortion_sum)
head(plot_cell_prop_num2);range(plot_cell_prop_num2$proportion)

plot_all_cell_pro <- ggplot(plot_cell_prop_num2, aes(Cell_type, proportion, fill = Group)) +
  geom_col(position = position_dodge(width = 0), width = 1.6, size = 0.1, colour = 'black') + #柱形图绘制
  geom_text(aes(label=abs(round(proportion,2)),y=proportion,x= Cell_type),size=3,vjust=0.5)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(),
        axis.text.x = element_text(vjust =1, hjust = 1,angle = 0)) +   #调整背景
  labs(x = '', y = 'Cell types Proportion(%)',title="Significant different Cell types(Hypergeometric test)") +   #设置坐标轴标签
  geom_hline(yintercept = 0, size = 0.5) + ylim(-30, 30)+
  scale_fill_manual(name="Age_group",values=my_morandi_colors[c(21,1)],breaks=c("CTRL", "Abortion"),labels=c("CTRL", "USM"))+
  coord_flip()#+ NoLegend()   scale_fill_brewer(palette = 'Accent')+coord_flip() 
plot_all_cell_pro 
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/MF1g_Barplot_cell_number_for_each_sub_celltype_group_split.pdf", plot_all_cell_pro,width=6, height=8)

##plot for significant different cell types
plot_cell_prop_sig<-cell_prop_raw3[which(cell_prop_raw3$p_HG<=0.05),]
plot_cell_prop_sig<-plot_cell_prop_sig[order(plot_cell_prop_sig$sum,decreasing = T),]
Cell_sig_name<-as.character(plot_cell_prop_sig$Cell_type)

plot_cell_prop_sig_num<-cell_prop_raw1[which(cell_prop_raw1$Cell_type %in% Cell_sig_name),]
plot_cell_prop_sig_num2<-merge(plot_cell_prop_sig_num,plot_cell_prop_sig[, c("Cell_type","sum","p_HG","q_adj_HG")])
plot_cell_prop_sig_num2$Cell_type<-factor(plot_cell_prop_sig_num2$Cell_type,levels=rev(Cell_sig_name))
plot_cell_prop_sig_num2$proportion<-ifelse(grepl("CTRL",plot_cell_prop_sig_num2$Group),plot_cell_prop_sig_num2$Freq*100/CTRL_sum,plot_cell_prop_sig_num2$Freq*-100/Abortion_sum)
head(plot_cell_prop_sig_num2);range(plot_cell_prop_sig_num2$proportion)

plot_all_cell_sig <- ggplot(plot_cell_prop_sig_num2, aes(Cell_type, proportion, fill = Group)) +
  geom_col(position = position_dodge(width = 0), width = 1.6, size = 0.1, colour = 'black') + #柱形图绘制
  geom_text(aes(label=abs(round(proportion,2)),y=proportion,x= Cell_type),size=3,vjust=0.5)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(),
        axis.text.x = element_text(vjust =1, hjust = 1,angle = 0)) +   #调整背景
  labs(x = '', y = 'Cell types Proportion(%)',title="Significant different Cell types(Hypergeometric test)") +   #设置坐标轴标签
  geom_hline(yintercept = 0, size = 0.5) + ylim(-30, 30)+
  scale_fill_manual(name="Age_group",values=my_morandi_colors[c(21,1)],breaks=c("CTRL", "Abortion"),labels=c("CTRL", "Abortion"))+
  coord_flip()#+ NoLegend()   scale_fill_brewer(palette = 'Accent')+coord_flip() 

plot_cell_prop_sig_pro<-plot_all_cell_sig+geom_text(aes(x= Cell_type, y=30, label = round(p_HG,3)),vjust = 0.5,color="black", size=3 ) 
plot_cell_prop_sig_pro
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Barplot_two_direction_cell_ratio_for_sig_HGt_sub_celltype.pdf", plot_cell_prop_sig_pro,width=8, height=6)

###Boxplo_cell_ratio wilcox test stastics for two group 
table_sample<-data.frame(table(target_final_noEry_mast$re_annotation_TFac,target_final_noEry_mast$sample_code))
table_sample2 <- dcast(table_sample, Var2 ~ Var1 , value.var = "Freq")
head(table_sample2)
table_sample3<-as.data.frame(t(apply(table_sample2[,2:ncol(table_sample2)],1,function(x) round(prop.table(x)*100,2))))
table_sample3$sample<-as.character(table_sample2$Var2)
table_sample3$Group<-as.character(table_sample3$sample)
table_sample3$Group<-ifelse(table_sample3$sample %in% paste0("CTRL_",1:5),"CTRL","Abortion")
table_sample3$Group<-factor(table_sample3$Group,levels=c("Abortion","CTRL"))
head(table_sample2)
head(table_sample3)
table_sample4 <- melt(table_sample3,variable.name="Cell_type",value.name = "cell_ratio",id.vars = c("sample","Group") )
head(table_sample4)
my_comparisons <- list(c("Abortion","CTRL"))
ratio_plot_subgroup<-ggboxplot(table_sample4, x = "Group", y = "cell_ratio",color = "Group", palette = my_morandi_colors[c(1,21)],add = "jitter", facet.by = c("Cell_type"))+#ylim(c(0,90))+
  theme(legend.position ="right",legend.direction = "vertical")+
  labs(x="Abortion group",y="Cell ratio",title="Cell ratio of diffetent cell types")+guides(fill=guide_legend(title="Abortion group"))
ratio_plot_wilcoxtest_1<-ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "p.format")+facet_wrap(~Cell_type, scales = "free",ncol =6 )
ratio_plot_wilcoxtest_2<-ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "..p.signif..")+facet_wrap(~Cell_type, scales = "free",ncol =6 )
ratio_plot_wilcoxtest_3<-ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "..p.signif..")+ facet_wrap(~Cell_type, scales = "free",ncol =6 )+ylim(0,max(table_sample4$cell_ratio)+0.5)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Boxplo_cell_ratio_for_wilcoxtest_sub_celltype1.pdf", ratio_plot_wilcoxtest_1,width=18, height=10)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Boxplo_cell_ratio_for_wilcoxtest_sub_celltype2.pdf", ratio_plot_wilcoxtest_2,width=18, height=10)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Boxplo_cell_ratio_for_wilcoxtest_sub_celltype3.pdf", ratio_plot_wilcoxtest_3,width=18, height=10)

###for major group
table_sample<-data.frame(table(target_final_noEry_mast$re_anno_TFac_major,target_final_noEry_mast$sample_code))
table_sample2 <- dcast(table_sample, Var2 ~ Var1 , value.var = "Freq")

table_sample3<-as.data.frame(t(apply(table_sample2[,2:ncol(table_sample2)],1,function(x) round(prop.table(x)*100,2))))
table_sample3$sample<-as.character(table_sample2$Var2)
table_sample3$Group<-as.character(table_sample3$sample)
table_sample3$Group<-ifelse(table_sample3$sample %in% paste0("CTRL_",1:5),"CTRL","Abortion")
table_sample3$Group<-factor(table_sample3$Group,levels=c("Abortion","CTRL"))
head(table_sample2)
head(table_sample3)
table_sample4 <- melt(table_sample3,variable.name="Cell_type",value.name = "cell_ratio",id.vars = c("sample","Group") )
head(table_sample4)
my_comparisons <- list(c("Abortion","CTRL"))
ratio_plot_subgroup<-ggboxplot(table_sample4, x = "Group", y = "cell_ratio",color = "Group", palette = my_morandi_colors[c(1,21)],add = "jitter", facet.by = c("Cell_type"))+#ylim(c(0,90))+
  theme(legend.position ="right",legend.direction = "vertical")+
  labs(x="Abortion group",y="Cell ratio",title="Cell ratio of diffetent cell types")+guides(fill=guide_legend(title="Abortion group"))

ratio_plot_wilcoxtest1<-ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "p.format")+facet_wrap(~Cell_type, scales = "free",ncol =4)
ratio_plot_wilcoxtest2<-ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "..p.signif..")+facet_wrap(~Cell_type, scales = "free",ncol =4)
ratio_plot_wilcoxtest3<-ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "..p.signif..")+ facet_wrap(~Cell_type, scales = "free",ncol =4)+ylim(0,max(table_sample4$cell_ratio)+0.5)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Boxplo_cell_ratio_for_wilcoxtest_major_celltype1.pdf", ratio_plot_wilcoxtest1,width=18, height=10)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Boxplo_cell_ratio_for_wilcoxtest_major_celltype2.pdf", ratio_plot_wilcoxtest2,width=18, height=10)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/Boxplo_cell_ratio_for_wilcoxtest_major_celltype3.pdf", ratio_plot_wilcoxtest3,width=18, height=10)