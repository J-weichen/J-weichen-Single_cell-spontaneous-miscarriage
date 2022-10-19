rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

require("Matrix")
library(data.table)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(data.table)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpubr)

#library(ggpubr)
grid.newpage()

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
show_col(pal)
pal
#[1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF"
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
show_col(nejm)
nejm
#[1] "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"
#decide color
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

#n=5;barplot(rep(1,n), col=colorRampPalette(colors = c('red', 'white'))( n ))
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
show_col(Cells_col_raw)

##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
show_col(pal1)
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
show_col(pal2)
pal3<- pal_aaas("default",alpha=1)(10)
show_col(pal3)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
show_col(pal4)
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
show_col(pal5)
pal5
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor_all)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]

##plot cells
Cells_col<-colorRampPalette(colors = rev(Cells_col_raw))(80)
length(unique(Cells_col))
barplot(rep(1,80), col=Cells_col)
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(80)

#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
library(future)
plan(strategy = "multicore", workers = 40)
options(future.globals.maxSize = 1500 * 1024^12)

#plot subtype for each major group
Cell_submian_list<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_major_recal_selected_population_list.rds")
names(Cell_submian_list)
#[1] "Endothelial_Cell_PECAM1"   "Erythrocyte_HBA1"          "Myeloid_Cell_AIF1"         "Natural_killer_Cell_NCAM1" "Stromal_cells_DCN"        
#[6] "T_Cell_CD3D"               "Troph_mix_group"           "Trophoblast_KRT7"         
P1<-DimPlot(object = Cell_submian_list[["Erythrocyte_HBA1"]], reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Erythrocyte_HBA1')
P2<-DimPlot(object =Cell_submian_list[["Endothelial_Cell_PECAM1"]], reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Endothelial_Cell_PECAM1')
P3<-DimPlot(object =Cell_submian_list[["Stromal_cells_DCN"]], reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Stromal_cells_DCN')
P4<-DimPlot(object =Cell_submian_list[["Trophoblast_KRT7"]], reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Trophoblast_KRT7')
P5<-DimPlot(object =Cell_submian_list[["Myeloid_Cell_AIF1"]], reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Myeloid_cells_AIF1')
P6<-DimPlot(object =Cell_submian_list[["Natural_killer_Cell_NCAM1"]], reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Natural_killer_Cell_NCAM1')
P7<-DimPlot(object =Cell_submian_list[["T_Cell_CD3D"]], reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('T_Cell_CD3D')
P8<-DimPlot(object =Cell_submian_list[["Troph_mix_group"]], reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Troph_mix_group')

type_plot_each_cells<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6,P7,P8),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/figure_one_result/prefinal_cluster_brief_type_for_each_major_cell.png", type_plot_each_cells,width=30, height=30)


##plot for all cells
target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")

DimPlot(target, group.by = "final_major_subgroup_brief",cols =Cells_col)
DimPlot(target, group.by = "final_major_subgroup_brief",cols =rev(Cells_col))

umap_pn<-DimPlot(object = target, group.by ="seurat_cluster_nature",label = TRUE,cols = ppCor_all2)  +NoLegend()
umap_p1<-DimPlot(object = target, reduction = "umap", group.by = "final_major_group_brief",cols =ppCor_all2,label =F) +ggtitle('final_major_group_brief') 
umap_p2<-DimPlot(object = target, reduction = "umap", group.by = "final_major_subgroup",cols =ppCor_all2,label =T) +ggtitle('prefinal_cluster_brief') + NoLegend() 
umap_p3<-DimPlot(object = target, reduction = "umap", group.by = "final_major_subgroup_brief",cols =ppCor_all2,label =T) +ggtitle('raw_cluster_brif') + NoLegend() 
umap_pn+umap_p1+umap_p2+umap_p3

target_final<-subset(x = target,subset = final_major_subgroup== "Troph_mix_group",invert = TRUE)
dim(target_final@meta.data)# 45800    56
table(target_final$final_major_subgroup_brief)
#Bs CTBs_0 CTBs_1 CTBs_2 CTBs_3   Endo    Epi    Ery EVTs_1 EVTs_2 EVTs_3    FBs    HCs  Masts   MyCs    NKs STBs_1 STBs_2   STCs     Ts 
# 37   2144   3283  11650    913    131     71   1136   7058   1858   3422   2712   1974     19   4413    664    963   1900   1033    419 
names(table(target_final$final_major_subgroup_brief))
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief, 
                                                 levels=c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2",
                                                          "MyCs","HCs","STCs","FBs","Endo","Epi",    
                                                          "NKs","Ts","Bs","Ery","Masts"),ordered=TRUE)
table(target_final$final_major_group_brief)
#B_Cell_CD79A    Endothelial_Cell_PECAM1 Epithelial_Cell_EPCAM_PAEP           Erythrocyte_HBA1            Mast_Cell_MS4A2 
#   37                        131                         71                       1136                         19 
#Myeloid_Cell_AIF1  Natural_killer_Cell_NCAM1          Stromal_cells_DCN                T_Cell_CD3D           Trophoblast_KRT7 
#  6387                        664                       3745                        419                      33191 

target_final$final_major_group_brief<- factor(target_final$final_major_group_brief, 
                                                 levels=c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN",
                                                          "Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP",
                                                          "Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A",
                                                          "Erythrocyte_HBA1","Mast_Cell_MS4A2"),ordered=TRUE)


umap_class_plot<-DimPlot(target_final, group.by = "final_major_subgroup_brief", reduction = "umap", cols =ppCor_all2)
tsne_class_plot<-DimPlot(target_final, group.by = "final_major_subgroup_brief", reduction = "tsne", cols =ppCor_all2)
pca_class_plot<-DimPlot(target_final, group.by = "final_major_subgroup_brief", reduction = "pca", cols =ppCor_all2)
CombinePlots(plots = list(umap_class_plot,tsne_class_plot,pca_class_plot),ncol = 3,legend = NULL)
umap_class_plot2<-DimPlot(target_final, group.by = "final_major_subgroup_brief", reduction = "umap", cols =ppCor_all2,label = TRUE)+ NoLegend()
tsne_class_plot2<-DimPlot(target_final, group.by = "final_major_subgroup_brief", reduction = "tsne", cols =ppCor_all2,label = TRUE)+ NoLegend()
pca_class_plot2<-DimPlot(target_final, group.by = "final_major_subgroup_brief", reduction = "pca", cols =ppCor_all2)
CombinePlots(plots = list(umap_class_plot2,tsne_class_plot2,pca_class_plot2),ncol = 3)

plot_Cell0<-DimPlot(target_final, group.by = "final_major_group_brief",cols =ppCor_all2[c(1:9,11)])#+ NoLegend()
plot_Cell0
LabelClusters(plot = plot_Cell0, id = 'final_major_group_brief',segment.color ="grey",colour = "black",fontface = "bAbortion",size = 4)+ NoLegend()
plot_Cell<-DimPlot(target_final, group.by = "final_major_subgroup_brief",cols =ppCor_all2)
LabelClusters(plot = plot_Cell, id = 'final_major_subgroup_brief',segment.color ="grey",colour = "black",fontface = "bAbortion",size = 4)+ NoLegend()
plot_Cell2<-DimPlot(target_final, group.by = "final_major_subgroup",cols =ppCor_all2)#+ NoLegend()
LabelClusters(plot = plot_Cell2, id = 'final_major_subgroup',segment.color ="grey",colour = "black",fontface = "bAbortion",size = 4)+ NoLegend()

P1<-DimPlot(object = target_final, reduction = "umap", group.by = "final_major_group_brief",cols =ppCor_all2[c(1:9,11)])
P2<-DimPlot(object = target_final, reduction = "umap", group.by = "final_major_subgroup_brief",cols =ppCor_all2)#+guides(fill = guide_legend(ncols =2))
P3<-DimPlot(object = target_final, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P4<-DimPlot(object = target_final, reduction = "umap",  group.by = "Treat",cols =ppCor_all2)
P5<-DimPlot(object = target_final, reduction = "umap",group.by = "sample",cols = ppCor_all2)
P6<-DimPlot(object = target_final, reduction = "umap",group.by = "sample_code",cols = ppCor_all2)
type_plot<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/figure_one_result/all_cell_type_umap.png"), type_plot,width=22, height=12)

#highlight each cell groups
##for major group
raw_population<-names(table(target_final$final_major_group_brief)) 
Cell_highlight_list<-list();i=1
for (population in raw_population) {
  print(as.character(population))
  i=i+1
  high_plot<-DimPlot(target_final, cols.highlight = ppCor[i], cols = "grey",
        cells.highlight = Cells(subset(x=target_final, subset = final_major_group_brief == population)))+ ggtitle(population)+NoLegend()
  Cell_highlight_list<-c(Cell_highlight_list,list(high_plot))
  }
CombinePlots(plots =Cell_highlight_list,ncol = 5,legend = NULL)
#for subgroups
raw_population<-names(table(target_final$final_major_subgroup_brief)) 
Cell_highlight_list<-list();i=1
for (population in raw_population) {
  print(as.character(population))
  i=i+1
  high_plot<-DimPlot(target_final, cols.highlight = ppCor_all2[i], cols = "grey",
                     cells.highlight = Cells(subset(x=target_final, subset = final_major_subgroup_brief == population)))+ ggtitle(population)+NoLegend()
  Cell_highlight_list<-c(Cell_highlight_list,list(high_plot))
}
CombinePlots(plots =Cell_highlight_list,ncol = 5,legend = NULL)

#cell distribution by each type split
P21<-DimPlot(object = target_final, reduction = "umap", group.by = "Phase",split.by="Phase",cols = ppCor[c(3:5)],ncol=3)
#P42<-DimPlot(object = target_final, reduction = "umap",group.by = "gender",split.by="gender",cols = ppCor[c(3:4)],ncol=2)
P22<-DimPlot(object = target_final, reduction = "umap",group.by = "Treat",split.by="Treat",cols = ppCor_all2[c(5,10)],ncol=2)
P23<-DimPlot(object = target_final, reduction = "umap",group.by = "sample_code",split.by="sample_code",cols = ppCor,ncol=3)

P31<-DimPlot(object = target_final, group.by = "final_major_group_brief",cols =ppCor_all2[c(1:9,11)],split.by="Treat",ncol=2)
P32<-DimPlot(object = target_final, group.by = "final_major_subgroup_brief",cols =ppCor_all2,split.by="Treat",ncol=2)
P33<-DimPlot(object = target_final,  group.by = "Phase",cols = ppCor[c(3:5)],split.by="Treat",ncol=2)
P31/P32/P33


#propetion bar sperated without Ery 
target_final_noEry<-subset(x = target_final,subset = final_major_group_brief== "Erythrocyte_HBA1",invert = TRUE)
target_final_noEry$final_major_group_brief<-factor(target_final_noEry$final_major_group_brief,levels = c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN",
                                                                                             "Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP",
                                                                                             "Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A",
                                                                                            "Mast_Cell_MS4A2"))
target_final_noEry$final_major_subgroup_brief<- factor(target_final_noEry$final_major_subgroup_brief, 
                                                 levels=c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2",
                                                          "MyCs","HCs","STCs","FBs","Endo","Epi",    
                                                          "NKs","Ts","Bs","Masts"),ordered=TRUE)

#计算去除红细胞后计算每组细胞数目
table(target_final_noEry$final_major_group_brief)
as.data.frame(table(target_final_noEry$final_major_subgroup_brief,target_final_noEry$Treat))
as.data.frame(table(target_final_noEry$Treat,target_final_noEry$final_major_subgroup_brief))


#pieplot for each type and add frequence
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_subgroup_brief, target_final_noEry$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
head(cell.prop1)
pieplot_ring1<-ggplot(cell.prop1,aes(x=Group,y=proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+coord_polar("y",start=1) + 
  theme_minimal()+ labs(title = "Abortion vs CTRL")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all2)#+facet_wrap(~Disease_group)
pieplot_ring1
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_group_brief, target_final_noEry$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
head(cell.prop1)
pieplot_ring2<-ggplot(cell.prop1,aes(x=Group,y=proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+coord_polar("y",start=1) + 
  theme_minimal()+ labs(title = "Abortion vs CTRL")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all2)#+facet_wrap(~Disease_group)
pieplot_ring2

major_subgroup_pieplot<-ggplot(as.data.frame(table(target_final_noEry$final_major_subgroup_brief)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=ppCor_all2)+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank()) 
major_group_pieplot<-ggplot(as.data.frame(table(target_final_noEry$final_major_group_brief)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=ppCor_all2)+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  
sample_pieplot<-ggplot(as.data.frame(table(target_final_noEry$sample_code)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=ppCor_all2)+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  
Phase_pieplot<-ggplot(as.data.frame(table(target_final_noEry$Phase)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=ppCor_all2)+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  
group_pieplot<-ggplot(as.data.frame(table(target_final_noEry$Treat)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ coord_polar("y",start=1) + 
  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ scale_fill_manual(values=ppCor_all2)+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())  

#propetion bar sperated by Abortion and CTRL 
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_subgroup_brief, target_final_noEry$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
Ratio_01<-ggplot(cell.prop1,aes(Group,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+ theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=ppCor_all2)+coord_flip() 
Ratio_01

cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_subgroup_brief, target_final_noEry$sample_code)))
colnames(cell.prop1)<-c("Cell_type","sample","proportion")
Ratio_02<-ggplot(cell.prop1,aes(sample,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+ theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=ppCor_all2)+coord_flip() 
Ratio_02

cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_group_brief, target_final_noEry$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
Ratio_03<-ggplot(cell.prop1,aes(Group,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+ theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=ppCor_all2)+coord_flip() 
Ratio_03

cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_group_brief, target_final_noEry$sample_code)))
colnames(cell.prop1)<-c("Cell_type","sample","proportion")
Ratio_04<-ggplot(cell.prop1,aes(sample,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+ theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=ppCor_all2)+coord_flip() 
Ratio_04

#各细胞类群中各样本占比
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_subgroup_brief, target_final_noEry$sample_code)))
colnames(cell.prop1)<-c("Cell_type","Cell_origin","proportion")
Ratio_1<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Cell_origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=ppCor)+guides(fill=guide_legend(title=NULL))
Ratio_1
##各细胞类群中不同状态占比
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_subgroup_brief, target_final_noEry$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
Ratio_2<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Group))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=ppCor[c(1:4)])+guides(fill=guide_legend(title=NULL))
Ratio_2
#各细胞类群中gender占比
##cell.prop1<-as.data.frame(prop.table(table(target_final$final_major_subgroup_brief, target_final$gender)))
##colnames(cell.prop1)<-c("Cell_type","gender","proportion")
##Ratio_3<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=gender))+
##  geom_bar(stat="identity",position="fill")+ggtitle("")+theme_bw()+
##  theme(axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
##  scale_fill_manual(values=ppCor[c(3:4)])+ guides(fill=guide_legend(title=NULL))
#各细胞类群中年龄比例
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_subgroup_brief, target_final_noEry$Phase)))
colnames(cell.prop1)<-c("Cell_type","Phase","proportion")
Ratio_3<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Phase))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=ppCor[c(1:4)])+guides(fill=guide_legend(title=NULL))
CombinePlots(plots = list(Ratio_1,Ratio_2,Ratio_3),ncol = 1,legend = NULL)

#各细胞类群中各样本占比
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_group_brief, target_final_noEry$sample_code)))
colnames(cell.prop1)<-c("Cell_type","Cell_origin","proportion")
Ratio_1<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Cell_origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=30,hjust=1, vjust=1,size=6), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=ppCor)+guides(fill=guide_legend(title=NULL))
Ratio_1
##各细胞类群中不同状态占比
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_group_brief, target_final_noEry$Treat)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
Ratio_2<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Group))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=30,hjust=1, vjust=1,size=6), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=ppCor[c(1:4)])+guides(fill=guide_legend(title=NULL))
Ratio_2
#各细胞类群中gender占比
##cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_group_brief, target_final_noEry$gender)))
##colnames(cell.prop1)<-c("Cell_type","gender","proportion")
##Ratio_3<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=gender))+
##  geom_bar(stat="identity",position="fill")+ggtitle("")+theme_bw()+
##  theme(axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
##  scale_fill_manual(values=ppCor[c(3:4)])+ guides(fill=guide_legend(title=NULL))
#各细胞类群中年龄比例
cell.prop1<-as.data.frame(prop.table(table(target_final_noEry$final_major_group_brief, target_final_noEry$Phase)))
colnames(cell.prop1)<-c("Cell_type","Phase","proportion")
Ratio_3<-ggplot(cell.prop1,aes(Cell_type,proportion,fill=Phase))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=30,hjust=1, vjust=1,size=6), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=ppCor[c(1:4)])+guides(fill=guide_legend(title=NULL))
CombinePlots(plots = list(Ratio_1,Ratio_2,Ratio_3),ncol = 1,legend = NULL)

##细胞占比的超几何检验
cell_prop<-as.data.frame(table(target_final_noEry$final_major_subgroup_brief, target_final_noEry$Treat))
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

length(Element_percentage);length(cnane_region)
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
head(cell_prop_raw1);dim(cell_prop_raw1)##49 
cell_name<-as.character(cell_prop_raw3[order(cell_prop_raw3$sum,decreasing = T),"Cell_type"])

plot_cell_prop_num<-cell_prop_raw1[which(cell_prop_raw1$Cell_type %in% cell_name),]
plot_cell_prop_num$Cell_type<-factor(plot_cell_prop_num$Cell_type,levels=cell_name)

plot_cell_prop_num_decrease<-ggplot(data=plot_cell_prop_num, mapping=aes(x= Cell_type,y=Freq,fill=Group))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+scale_fill_manual(name="Types",values=ppCor_all[1:2])+
  theme_classic()+labs(x="Cell types",y="Cell number",title="The Cell number for each cells")+
  scale_y_continuous(breaks = seq(0,15000,by=1000))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 8,colour = 'black',vjust=1,hjust=1,angle = 60),legend.title = element_text(size = 9))
plot_cell_prop_num_decrease
plot_cell_prop_num_decrease2<-plot_cell_prop_num_decrease+geom_text(aes(label=Freq),position=position_stack(vjust=0.5))#3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部

##plot for significant different cell types
plot_cell_prop_sig<-cell_prop_raw3[which(cell_prop_raw3$p_HG<=0.05),]
plot_cell_prop_sig<-plot_cell_prop_sig[order(plot_cell_prop_sig$sum,decreasing = T),]
Cell_sig_name<-as.character(plot_cell_prop_sig$Cell_type)

plot_cell_prop_sig_num<-cell_prop_raw1[which(cell_prop_raw1$Cell_type %in% Cell_sig_name),]
plot_cell_prop_sig_num2<-merge(plot_cell_prop_sig_num,plot_cell_prop_sig[, c("Cell_type","sum","p_HG","q_adj_HG")])
plot_cell_prop_sig_num2$Cell_type<-factor(plot_cell_prop_sig_num2$Cell_type,levels=Cell_sig_name)

plot_cell_prop_sig_num3<-ggplot(data=plot_cell_prop_sig_num2, mapping=aes(x= Cell_type,y=Freq,fill=Group))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+scale_fill_manual(name="Types",values=ppCor_all[1:2])+
  theme_classic()+labs(x="Cell types",y="Cell number",title="The Cell number for each cells(Hypergeometric test)")+
  scale_y_continuous(breaks = seq(0,15000,by=1000))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 8,colour = 'black',vjust=1,hjust=1,angle = 60),legend.title = element_text(size = 9))
plot_cell_prop_sig_num4<-plot_cell_prop_sig_num3+geom_text(aes(x= Cell_type, y=sum, label = round(p_HG,3)),vjust = -1,color="black", size=3 )
plot_cell_prop_sig_num5<-plot_cell_prop_sig_num3+geom_text(aes(x= Cell_type, y=sum, label = sum),vjust = -1,color="black", size=3 )

##Barplot showing proportion Seperated by Abortion or CTRL 
plot_cell_prop_sig_num2$Cell_type<-factor(plot_cell_prop_sig_num2$Cell_type,levels=rev(Cell_sig_name))
plot_cell_prop_sig_num2$proportion<-ifelse(grepl("CTRL",plot_cell_prop_sig_num2$Group),plot_cell_prop_sig_num2$Freq*100/CTRL_sum,plot_cell_prop_sig_num2$Freq*-100/Abortion_sum)
head(plot_cell_prop_sig_num2);range(plot_cell_prop_sig_num2$proportion)

plot_all_cell_sig <- ggplot(plot_cell_prop_sig_num2, aes(Cell_type, proportion, fill = Group)) +
  geom_col(position = position_dodge(width = 0), width = 1.6, size = 0.1, colour = 'black') + #柱形图绘制
  geom_text(aes(label=abs(round(proportion,2)),y=proportion,x= Cell_type),size=3,vjust=0.5)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(),
        axis.text.x = element_text(vjust =1, hjust = 1,angle = 0)) +   #调整背景
  labs(x = '', y = 'Cell types Proportion(%)',title="Significant different Cell types(Hypergeometric test)") +   #设置坐标轴标签
  geom_hline(yintercept = 0, size = 0.5) + ylim(-40, 40)+
  scale_fill_manual(name="Age_group",values=ppCor_all[2:1],breaks=c("CTRL", "Abortion"),labels=c("CTRL", "Abortion"))+
  coord_flip()#+ NoLegend()   scale_fill_brewer(palette = 'Accent')+coord_flip() 

plot_cell_prop_sig_pro<-plot_all_cell_sig+geom_text(aes(x= Cell_type, y=40, label = round(p_HG,3)),vjust = 0.5,color="black", size=3 ) 

#plot_cell_prop_sig_num5+ coord_flip()+ scale_y_reverse()
plot_cell_prop_sig_pro
plot_cell_prop_sig_num4
plot_cell_prop_sig_num5
plot_cell_prop_num_decrease2


###T-test stastics for two group 
table_sample<-data.frame(table(target_final_noEry$final_major_subgroup_brief,target_final_noEry$sample_code))
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
ratio_plot_subgroup<-ggboxplot(table_sample4, x = "Group", y = "cell_ratio",color = "Group", palette = ppCor,add = "jitter", facet.by = c("Cell_type"))+#ylim(c(0,90))+
  theme(legend.position ="right",legend.direction = "vertical")+
  labs(x="Abortion group",y="Cell ratio",title="Cell ratio of diffetent cell types")+guides(fill=guide_legend(title="Abortion group"))
ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "p.format")+facet_wrap(~Cell_type, scales = "free",ncol =5 )
ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "..p.signif..")+facet_wrap(~Cell_type, scales = "free",ncol =5 )
ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "..p.signif..")+
  facet_wrap(~Cell_type, scales = "free",ncol =5 )+ylim(0,max(table_sample4$cell_ratio)+0.5)
###for major group
table_sample<-data.frame(table(target_final_noEry$final_major_group_brief,target_final_noEry$sample_code))
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
ratio_plot_subgroup<-ggboxplot(table_sample4, x = "Group", y = "cell_ratio",color = "Group", palette = ppCor,add = "jitter", facet.by = c("Cell_type"))+#ylim(c(0,90))+
  theme(legend.position ="right",legend.direction = "vertical")+
  labs(x="Abortion group",y="Cell ratio",title="Cell ratio of diffetent cell types")+guides(fill=guide_legend(title="Abortion group"))
ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "p.format")+facet_wrap(~Cell_type, scales = "free",ncol =5 )
ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "..p.signif..")+facet_wrap(~Cell_type, scales = "free",ncol =5 )
ratio_plot_subgroup+ stat_compare_means(method = "wilcox.test",label.x=1.5,label = "..p.signif..")+
  facet_wrap(~Cell_type, scales = "free",ncol =5 )+ylim(0,max(table_sample4$cell_ratio)+0.5)


unique(target_final$merge_cluster)
table(target_final$merge_cluster)

#"Trophoblast_mix_EECs_KRT7" :: 5837
#"B_Cell_CD79A" ::5407 
#"T_NK_raw_Erythrocyte_Mast_Cell": 80209 
#"Myeloid_Cell_AIF1":37951                            
#"Decidual_stromal_cells_DKK1_ACTA2"::11439           
#"Fetal_Stromal_cells_Endothelial_cells_DLK1" :: 5316
#"Lymphatic_Endothelial_Cells_PECAM1_LYVE1_pos" ::1824
#"Vascular_Endothelial_Cells_PECAM1_LYVE1_neg" ::720 
#"Erythrocyte_HBA1" :652                                                           
#"Megakaryocytes_PPBP"::214          
#"Epithelial_Cell_EPCAM_PAEP":: 129 

#For major cell population
FeaturePlot(object = target_final, features = c("HLA-B","KRT7","PERP","PAEP","DCN","DLK1","PECAM1","LYVE1","AIF1","CD3D","NKG7","CD79A","MS4A2","CD34","HBA1","PPBP"),cols= c("grey", "red"),ncol=4)

#For no_immunine Cells 
FeaturePlot(object = target_final, features = c("VIM","CDH1","EGFR","HLA-G","MMP11","CGA","CYP19A1","CSH2","ERVFRD-1",
                                                "DCN","DLK1","THY1","COL1A1","LAbortion2","TIMP1",
                                                "MYH11","RGS5","PRL","IGFBP1","APOD","COL6A2",
                                                "PECAM1","LYVE1","EGFL7","PPBP","PF4","EPCAM","PAEP"),cols= c("grey", "red"),ncol=6)

#For immunine Cells 
FeaturePlot(object = target_final, features = c("AIF1","CD14","S100A8","FCN1","FCGR3A","CLEC9A","CD1C",
                                                "CSF1R","CD163","CD209","CD69","LYVE1","LYZ","APOE","MS4A3",
                                                "MS4A1","CD79A","CD79B","JCHAIN",
                                                "CD3D","CD8A","IL7R","NKG7","KLRB1","FGFBP2","S100A4","HBB","KIT"),cols= c("grey", "red"),ncol=6)



#target_final<-readRDS(file = "/home/chenwei/10x_data/191125-merge_nonormalization_add/Cell_0512.rds")
#plot for main nine population
target_final$main_nine_population<-as.character(target_final$merge_cluster)
Idents(object = target_final) <- "final_major_subgroup_brief"
levels(Idents(target_final))
sub_cell<-subset(x = target_final,idents=c("Ery"))
table(sub_cell@meta.data$Tissue)
target_final$main_nine_population[Cells(sub_cell)] <- "Erythrocyte_HBA1"
sub_cell<-subset(x = target_final,idents=c("D_LECs"))
table(sub_cell@meta.data$Tissue)
target_final$main_nine_population[Cells(sub_cell)] <- "Decidual_lymphatic_endothelial_Cell_PECAM1_LYVE1"

Idents(object = target_final) <- "main_nine_population"
levels(Idents(target_final))
target_final$main_nine_population[Cells(subset(x = target_final,idents=c("Trophoblast_mix_EECs_KRT7")))] <- "Trophoblast_Cell_PERP"
target_final$main_nine_population[Cells(subset(x = target_final,idents=c("T_NK_raw_Erythrocyte_Mast_Cell")))] <- "T_NK_Cell_CD3D_NKG7"
target_final$main_nine_population[Cells(subset(x = target_final,idents=c("Endometrial_Epithelial_Cell_EPCAM_PAEP")))] <- "Epithelial_Cell_PAEP"
target_final$main_nine_population[Cells(subset(x = target_final,idents=c("Endothelial_Cell_PECAM1")))] <- "Vascular_Endothelial_Cell_PECAM1"
target_final$main_nine_population[Cells(subset(x = target_final,idents=c("Fetal_Stromal_cells_Endothelial_cells_DLK1")))] <- "Fetal_Stromal_cells_DLK1"
target_final$main_nine_population[Cells(subset(x = target_final,idents=c("Decidual_stromal_cells_DKK1_ACTA2")))] <- "Decidual_stromal_cells_DCN_DKK1_ACTA2"

names(table(target_final$main_nine_population))
target_final$main_nine_population

target_final$main_nine_population<- factor(target_final$main_nine_population, 
                                           levels=c("T_NK_Cell_CD3D_NKG7","Myeloid_Cell_AIF1","Decidual_stromal_cells_DCN_DKK1_ACTA2",
                                                    "Trophoblast_Cell_PERP", "B_Cell_CD79A","Fetal_Stromal_cells_DLK1",
                                                    "Decidual_lymphatic_endothelial_Cell_PECAM1_LYVE1","Vascular_Endothelial_Cell_PECAM1",
                                                    "Erythrocyte_HBA1","Megakaryocytes_PPBP","Epithelial_Cell_PAEP"),ordered=TRUE)

DimPlot(object = target_final, reduction = "umap",group.by = "main_nine_population",cols = ppCor)
main_one<-DimPlot(object = target_final, reduction = "umap",group.by = "main_nine_population",cols = ppCor,label = TRUE)+ NoLegend() 
main_two<-FeaturePlot(object = target_final, features = c("HLA-B","KRT7","PERP","DCN","DKK1","ACTA2","DLK1","PECAM1","LYVE1","EPCAM","PAEP","AIF1","CD3D","NKG7","CD79A","PPBP","PF4","HBA1"),cols= c("grey", "red"),ncol=5)
main_one+main_two

#plot for fetal maternal mix population
target_final$Fetal_maternal_origin<-as.character(target_final$final_major_subgroup_brief)
Idents(object = target_final) <- "final_major_subgroup_brief"
Idents(target_final)
sub_cell<-subset(x = target_final,idents=c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2","STBs_3",
                                           "F_STCs","F_FBs_1","F_FBs_2","F_FBs_3","F_FBs_VECs","F_VECs","HBs","F_Nai_CD4_Ts"))
table(sub_cell@meta.data$Tissue)
target_final$Fetal_maternal_origin[Cells(sub_cell)] <- "pure_fetal"

sub_cell<-subset(x = target_final,idents=c("D_STCs","D_FBs_1","D_FBs_2","D_FBs_3","MFs","D_LECs","D_VECs","Epis","Mac_1","Mac_2",
                                           "nc_Mon_1","B_c_Mon","B_Bs","MALT_Bs","Aty_Bs","D_NKs","D_Cyt_Ts","D_Mem_CD4_Ts","B_NKs",
                                           "B_Cyt_Ts","B_NKTs","B_Mem_CD4_Ts", "B_Nai_CD4_Ts","Mast"))

table(sub_cell@meta.data$Tissue)
target_final$Fetal_maternal_origin[Cells(sub_cell)] <- "pure_maternal"

sub_cell<-subset(x = target_final,idents=c("TroB_VIM","F_AIF1_NK","F_AIF1_T","PVs","MKCs_1","MKCs_2","GMPs","DCs","Mac_3","nc_Mon_2","c_Mon_2","c_Mon_3","Trm_Bs","NB_NKs","Ery"))
table(sub_cell@meta.data$Tissue)
target_final$Fetal_maternal_origin[Cells(sub_cell)] <- "Mix"
table(target_final$Fetal_maternal_origin)
##  Mix    pure_fetal pure_maternal 
##21486         20801        113748 
