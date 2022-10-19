rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(Seurat)
require("Matrix")
library(data.table)
library(cowplot)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(data.table)
library(grid)
library(gridExtra)
library(cowplot)
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
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(80)
show_col(Cells_col)

#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
library(future)
plan(strategy = "multicore", workers = 8)
options(future.globals.maxSize = 1500 * 1024^12)

#step1 ：read Seurat object and load initial coldata and matrix
single_cell_placent_all<-readRDS(file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_all.rds")
DefaultAssay(single_cell_placent_all) <- "SCT"

target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
table(target_final$sample_code)
# CTRL_1   CTRL_2   CTRL_3   CTRL_4   CTRL_5 Arrest_1 Arrest_2 Arrest_3 Arrest_4 
# 11044     6520     4533     4249     3171     2905     5628     4170     3580 
mydata_list<-c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5")
#data_name<-c()
seurat_data_list<-list()
for ( sample_name in mydata_list){
  # sample_name<-"CTRL_1"#test line
  print(as.character(sample_name))
  selected_cell<-subset(x = target_final, subset = sample_code == sample_name)
  DefaultAssay(selected_cell) <- "RNA"    ## very important
  selected_cell <- SCTransform(object = selected_cell,verbose = FALSE) 
  seurat_data_list<-c(seurat_data_list,list(selected_cell))
}

seurat_data_list2<-c(single_cell_placent_all,seurat_data_list)
names(seurat_data_list2)<-c("nature",mydata_list)

dtlist <- seurat_data_list2
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 5000) #2976 
dtlist <- PrepSCTIntegration(object.list = dtlist,anchor.features = intfts)
anchors <- FindIntegrationAnchors(object.list = dtlist, normalization.method = "SCT", anchor.features = intfts)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated <- RunPCA(integrated)
ElbowPlot(object = integrated)
VizDimLoadings(object = integrated, dims = 1:2, reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:25)
integrated<- RunTSNE(integrated, dims = 1:25)

integrated <- FindNeighbors(object = integrated, dims = 1:25,verbose = FALSE)
integrated <- FindClusters(object = integrated,resolution = 0.6, verbose = FALSE)

######
plot0<-DimPlot(integrated, reduction = "pca")
plot1<-DimPlot(integrated, reduction = "umap")
plot2<-DimPlot(integrated, reduction = "tsne")
CombinePlots(plots = list(plot0,plot1,plot2),legend="top",ncol=3)

head(integrated@meta.data);tail(integrated@meta.data)
table(integrated@meta.data$orig.ident)
DimPlot(object = integrated, group.by ="annotation",label = TRUE) + NoLegend() 

meta_data_new<-integrated@meta.data
meta_data_new$data_source<-meta_data_new$annotation
meta_data_new[which(!(is.na(meta_data_new$annotation))),]$data_source<-"nature_data"
meta_data_new[which(is.na(meta_data_new$annotation)),]$data_source<-"gong_chen_data"
meta_data_new[which(is.na(meta_data_new$annotation)),]$annotation<-"gong_chen_data"
table(meta_data_new$data_source)
#gong_chen_data    nature_data 
#   29517          64734 
table(meta_data_new$annotation)
meta_data_new$data_source<- factor(meta_data_new$data_source,levels = c("gong_chen_data","nature_data"))
integrated@meta.data<-meta_data_new
saveRDS(integrated, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/final_CTRL_GC_seurat_object_CCA_with_nature_paper.rds")

##plot
integrated<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/final_CTRL_GC_seurat_object_CCA_with_nature_paper.rds")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
integrated$re_annotation_TFac<- factor(integrated$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
integrated$re_anno_TFac_major<- factor(integrated$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

DimPlot(object = integrated, group.by ="annotation",split.by ="data_source",label = TRUE,cols = ppCor_all2) 
DimPlot(object = integrated, group.by ="annotation",label = TRUE,cols = ppCor_all2)
DimPlot(object = integrated, group.by ="re_annotation_TFac",split.by ="data_source",label = TRUE,cols = ppCor_all2) 
DimPlot(object = integrated, group.by ="re_annotation_TFac",label = TRUE,cols = ppCor_all2) 

plot_annotation1<-DimPlot(object = integrated, group.by ="re_annotation_TFac",split.by ="Treat",label = TRUE,cols = ppCor_all2) 
plot_annotation2<-DimPlot(object = integrated, group.by ="re_annotation_TFac",label = TRUE,cols = ppCor_all2) 
plot_annotation3<-DimPlot(object = integrated, group.by ="annotation",split.by ="Treat",label = TRUE,cols = ppCor_all2) 
plot_annotation4<-DimPlot(object = integrated, group.by ="annotation",label = TRUE,cols = ppCor_all2) 

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_1.pdf",plot_annotation1,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_2.pdf",plot_annotation2,width=7, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_3.pdf",plot_annotation3,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_4.pdf",plot_annotation4,width=8, height=6)

meta_data_new<-integrated@meta.data
head(meta_data_new)
table(meta_data_new$location)
# Blood  Decidua Placenta 
# 10001    36186    18547 
table(meta_data_new$annotation)
table(meta_data_new$re_anno_TFac_major)
table(meta_data_new$re_annotation_TFac)

conflict_prefer("which", "Matrix")
meta_data_new$major_group_brief<-as.character(meta_data_new$re_anno_TFac_major)
meta_data_new[which(meta_data_new$annotation %in% c("DC1","DC2","dM1","dM2","dM3","HB","MO")),]$major_group_brief<- "Myeloid_Cell_AIF1"
meta_data_new[which(meta_data_new$annotation %in% c("dNK p","dNK1","dNK2","dNK3","NK CD16+","NK CD16-","Tcells","Plasma","ILC3","Granulocytes","Plasma")),]$major_group_brief<- "T_NK_B_Cell"
meta_data_new[which(meta_data_new$re_anno_TFac_major %in% c("B_Cell_CD79A","T_Cell_CD3D","Natural_killer_Cell_NCAM1","Mast_Cell_MS4A2")),]$major_group_brief<- "T_NK_B_Cell"

meta_data_new[which(meta_data_new$annotation %in% c("dP1","dP2","dS1","dS2","dS3")),]$major_group_brief<- "Stromal_cells_DCN"
meta_data_new[which(meta_data_new$re_annotation_TFac %in% c("STCs")),]$major_group_brief<- "Stromal_cells_DCN"
meta_data_new[which(meta_data_new$annotation %in% c("fFB1","fFB2")),]$major_group_brief<- "Stromal_cells_DCN"
meta_data_new[which(meta_data_new$re_annotation_TFac %in% c("FBs")),]$major_group_brief<- "Stromal_cells_DCN"

meta_data_new[which(meta_data_new$annotation %in% c("Endo (f)","Endo (m)","Endo L")),]$major_group_brief<- "Endothelial_Cell_PECAM1"
meta_data_new[which(meta_data_new$annotation %in% c("SCT","VCT","EVT")),]$major_group_brief<- "Trophoblast_KRT7"
meta_data_new[which(meta_data_new$annotation %in% c("Epi1","Epi2")),]$major_group_brief<- "Epithelial_Cell_EPCAM_PAEP"

table(meta_data_new$re_anno_TFac_major)
table(meta_data_new$major_group_brief)
table(meta_data_new$annotation)
integrated@meta.data<-meta_data_new

integrated2<-subset(x = integrated, subset = major_group_brief == "Erythrocyte_HBA1",invert = TRUE)
integrated2$major_group_brief<-factor(integrated2$major_group_brief,levels = c("Stromal_cells_DCN","Endothelial_Cell_PECAM1","Myeloid_Cell_AIF1","T_NK_B_Cell","Epithelial_Cell_EPCAM_PAEP","Trophoblast_KRT7"))
plot_annotation1<-DimPlot(object = integrated2, group.by ="major_group_brief",split.by ="data_source",cols = ppCor_all2[c(1:3,7,5:6)]) 

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_data_source_split_no_Ery.pdf",plot_annotation1,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_data_source_split_no_Ery.png",plot_annotation1,width=13, height=6)

##for cell ratio 
Metadate_plot<-integrated2@meta.data
Metadate_plot$class<-Metadate_plot$Treat
Metadate_plot[which(Metadate_plot$location =="Placenta"),]$class<- "Placenta"
Metadate_plot[which(Metadate_plot$location =="Decidua"),]$class<- "Decidua"
Metadate_plot[which(Metadate_plot$location =="Blood"),]$class<- "Blood"
Metadate_plot$major_group_brief<-factor(Metadate_plot$major_group_brief,levels = c("Stromal_cells_DCN","Endothelial_Cell_PECAM1","Myeloid_Cell_AIF1","T_NK_B_Cell","Epithelial_Cell_EPCAM_PAEP","Trophoblast_KRT7"))
Metadate_plot2<-Metadate_plot[which(Metadate_plot$class %in% c("Placenta","CTRL")),]

dim(Metadate_plot2)#47210    71
cell.prop1<-as.data.frame(prop.table(table(Metadate_plot2$major_group_brief, Metadate_plot2$class)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
cell.prop1$Group<-factor(cell.prop1$Group,levels = c("Placenta","CTRL"))

Ratio_01<-ggplot(cell.prop1,aes(Group,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+ theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=ppCor_all2[c(1:3,7,5:6)])+coord_flip() 
Ratio_01
pieplot_ring3<-ggplot(cell.prop1,aes(x=Group,y=proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+coord_polar("y",start=1) + 
  theme_minimal()+ labs(title = "scNature vs CTRL")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all2[c(1:3,7,5:6)])#+facet_wrap(~Disease_group)
pieplot_ring3

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_CTRL_scNATURE/major_Cell_ratio_for_CTRL_GC_seurat_object_sc_nature_placenta_intergrate_Group_bar.pdf",Ratio_01,width=8, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_CTRL_scNATURE/major_Cell_ratio_for_CTRL_GC_seurat_object_sc_nature_placenta_intergrate_Group_ring_pie.pdf",pieplot_ring3,width=8, height=6)

plot_annotation3<-DimPlot(object = integrated, group.by ="annotation",split.by ="data_source",label = TRUE,cols = ppCor_all2) 
plot_annotation4<-DimPlot(object = integrated, group.by ="annotation",label = TRUE,cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_split.pdf",plot_annotation3,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_split.png",plot_annotation3,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_merge.pdf",plot_annotation4,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annotation_merge.png",plot_annotation4,width=6, height=6)

plot_annotation10<-DimPlot(object = integrated, group.by ="location",split.by ="Treat",label = TRUE,cols = ppCor_all2) 
plot_annotation11<-DimPlot(object = integrated, group.by ="location",label = TRUE,cols = ppCor_all2) 
plot_annotation12<-DimPlot(object = integrated, group.by ="location",split.by ="location",cols = ppCor_all2) 
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annot_split_5.pdf",plot_annotation10,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annot_split_6.pdf",plot_annotation11,width=7, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/CTRL_scNATURE/CTRL_GC_sample_sc_nature_intergrate_Group_umap_annot_split_7.pdf",plot_annotation12,width=25, height=6)