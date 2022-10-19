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
library(moRandi)
library(data.table)
library(grid)
library(gridExtra)
library(cowplot)
##set colour
x <- c(30,4,1,2,3,20,
       26,29,37,41,6,
       7,8,51,39,42,
       56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
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
mydata_list<-c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","Arrest_1","Arrest_2","Arrest_3","Arrest_4")
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
saveRDS(integrated, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/final_GC_seurat_object_CCA_with_nature_paper.rds")

###plot figure
integrated<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/final_GC_seurat_object_CCA_with_nature_paper.rds")

#extract metadata information
metadata_Cell<-integrated@meta.data
final_umap<-data.frame(Embeddings(integrated[["umap"]]))
#final_umap$raw_id <- rownames(final_umap)
head(metadata_Cell);head(final_umap)
metadata_Cell_Umap_add <- merge(metadata_Cell,final_umap,by=0)
head(metadata_Cell_Umap_add)
rownames(metadata_Cell_Umap_add)<-metadata_Cell_Umap_add$Row.names

metadata_Cell_Umap_add$data_source<-metadata_Cell_Umap_add$annotation
metadata_Cell_Umap_add[which(!(is.na(metadata_Cell_Umap_add$annotation))),]$data_source<-"nature_data"
metadata_Cell_Umap_add[which(is.na(metadata_Cell_Umap_add$annotation)),]$data_source<-"gong_chen_data"
#metadata_Cell_Umap_add[which(is.na(metadata_Cell_Umap_add$annotation)),]$annotation<-"gong_chen_data"
metadata_Cell_Umap_add[which(is.na(metadata_Cell_Umap_add$Treat)),]$Treat<-"Nature"


table(metadata_Cell_Umap_add$data_source)
#gong_chen_data    nature_data 
#    45800           64734 
table(metadata_Cell_Umap_add$Treat)
#Abortion     CTRL   Nature 
# 16283    29517    64734 

table(metadata_Cell_Umap_add$annotation)
metadata_Cell_Umap_add$data_source<- factor(metadata_Cell_Umap_add$data_source,levels = c("gong_chen_data","nature_data"))
metadata_Cell_Umap_add$Treat<- factor(metadata_Cell_Umap_add$Treat,levels = c("Nature","CTRL","Abortion"))

head(metadata_Cell_Umap_add)
table(metadata_Cell_Umap_add$location)
# Blood  Decidua Placenta 
# 10001    36186    18547 
table(metadata_Cell_Umap_add$annotation)
table(metadata_Cell_Umap_add$re_anno_TFac_major)
table(metadata_Cell_Umap_add$re_annotation_TFac)

metadata_Cell_Umap_add$major_group_brief<-as.character(metadata_Cell_Umap_add$re_anno_TFac_major)
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("DC1","DC2","dM1","dM2","dM3","HB","MO")),]$major_group_brief<- "Myeloid_Cell_AIF1"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("dNK p","dNK1","dNK2","dNK3","NK CD16+","NK CD16-","Tcells","Plasma","ILC3","Granulocytes")),]$major_group_brief<- "T_NK_B_Cell"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$re_anno_TFac_major %in% c("B_Cell_CD79A","T_Cell_CD3D","Natural_killer_Cell_NCAM1","Mast_Cell_MS4A2")),]$major_group_brief<- "T_NK_B_Cell"

metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("dP1","dP2","dS1","dS2","dS3")),]$major_group_brief<- "Stromal_cells_DCN"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$re_annotation_TFac %in% c("STCs")),]$major_group_brief<- "Stromal_cells_DCN"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("fFB1","fFB2")),]$major_group_brief<- "Stromal_cells_DCN"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$re_annotation_TFac %in% c("FBs")),]$major_group_brief<- "Stromal_cells_DCN"

metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("Endo (f)","Endo (m)","Endo L")),]$major_group_brief<- "Endothelial_Cell_PECAM1"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("SCT","VCT","EVT")),]$major_group_brief<- "Trophoblast_KRT7"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("Epi1","Epi2")),]$major_group_brief<- "Epithelial_Cell_EPCAM_PAEP"
table(metadata_Cell_Umap_add$major_group_brief)

##reannotation for cell subtypes
metadata_Cell_Umap_add$reannotation<-metadata_Cell_Umap_add$annotation
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("MO","dM1","dM2","dM3","DC1","DC2","Granulocytes")),]$reannotation <- "MyCs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("dNK p","dNK1","dNK2","dNK3","Tcells", "ILC3","Plasma","NK CD16-","NK CD16+")),]$reannotation <- "TsBsNK"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("dP1","dP2","dS1","dS2","dS3")),]$reannotation <- "dSs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("Endo (f)","Endo (m)","Endo L")),]$reannotation <- "Endos"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("Epi1","Epi2")),]$reannotation <- "Epis"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("fFB1","fFB2")),]$reannotation <- "fFBs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("VCT")),]$reannotation <- "CTBs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("EVT")),]$reannotation <- "EVTs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$annotation %in% c("SCT")),]$reannotation <- "STBs"

table(metadata_Cell_Umap_add$reannotation)
names(table(metadata_Cell_Umap_add$reannotation))

metadata_Cell_Umap_add$formal_cell_types<-ifelse(metadata_Cell_Umap_add$data_source =="nature_data",metadata_Cell_Umap_add$reannotation,as.character(metadata_Cell_Umap_add$re_annotation_TFac))
metadata_Cell_Umap_add$formal_cell_types2<-metadata_Cell_Umap_add$formal_cell_types
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types %in% c("Epi","Epis")),]$formal_cell_types2<- "Epis"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types %in% c( "FBs","fFBs")),]$formal_cell_types2<- "FBs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types %in% c("dSs","STCs")),]$formal_cell_types2<- "STCs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types %in% c("Endo","Endos")),]$formal_cell_types2<- "Endos"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types %in% c("Mycs","MyCs")),]$formal_cell_types2<- "MyCs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types %in% c("NKs","Ts","Bs","TsBsNK")),]$formal_cell_types2<- "TsBsNK"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types %in% c("HB","HCs")),]$formal_cell_types2<- "HCs"

metadata_Cell_Umap_add$formal_cell_types3<-metadata_Cell_Umap_add$formal_cell_types2
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types2 %in% c("CTBs","CTBs_1","CTBs_2")),]$formal_cell_types3<- "CTBs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types2 %in% c("EVTs","EVTs_1","EVTs_2","EVTs_3")),]$formal_cell_types3<- "EVTs"
metadata_Cell_Umap_add[which(metadata_Cell_Umap_add$formal_cell_types2 %in% c("STBs","STBs_1","STBs_2","STBs_3")),]$formal_cell_types3<- "STBs"


##plot for major cell types
table(metadata_Cell_Umap_add$major_group_brief)
major_order<-c("Trophoblast_KRT7","Stromal_cells_DCN","Myeloid_Cell_AIF1","Endothelial_Cell_PECAM1","T_NK_B_Cell","Epithelial_Cell_EPCAM_PAEP","Erythrocyte_HBA1")
metadata_Cell_Umap_add$major_group_brief<- factor(metadata_Cell_Umap_add$major_group_brief, levels=major_order,ordered=TRUE)
pos_cell_umap1<-ggplot(data=metadata_Cell_Umap_add[order(metadata_Cell_Umap_add$major_group_brief ,decreasing = F),], 
                       mapping=aes(x=UMAP_1,y=UMAP_2,colour = major_group_brief ))+
  geom_point(stat= "identity",size=0.5,alpha=0.5,show.legend = TRUE)+
  #labs(title =paste0("Expression of ",genename))+
  scale_color_manual(values=my_morandi_colors[c(2,1,3,9,6,7,8)])+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                   axis.text.y = element_text(size=15),legend.position = "right")
pos_cell_umap12<-pos_cell_umap1+facet_wrap(~data_source)
pos_cell_umap13<-pos_cell_umap1+facet_wrap(~Treat)


cell_order2<-c("CTBs","STBs","EVTs","Epis","Endos","HCs","MyCs","STCs","FBs","TsBsNK","Ery","Masts")
metadata_Cell_Umap_add$formal_cell_types3<- factor(metadata_Cell_Umap_add$formal_cell_types3, levels=cell_order2,ordered=TRUE)

pos_cell_umap2<-ggplot(data=metadata_Cell_Umap_add[order(metadata_Cell_Umap_add$formal_cell_types3 ,decreasing = F),], 
                       mapping=aes(x=UMAP_1,y=UMAP_2,colour = formal_cell_types3 ))+
  geom_point(stat= "identity",size=0.5,alpha=0.5,show.legend = TRUE)+
  #labs(title =paste0("Expression of ",genename))+
  scale_color_manual(values=my_morandi_colors[c(1,3,6,9,10,11,12,13,14,15,18,19)])+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                   axis.text.y = element_text(size=15),legend.position = "right")
pos_cell_umap22<-pos_cell_umap2+facet_wrap(~data_source)
pos_cell_umap23<-pos_cell_umap2+facet_wrap(~Treat)

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_major_annotation_treat_split.pdf",pos_cell_umap13,width=16, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_major_annotation_treat_split.png",pos_cell_umap13,width=16, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_major_annotation_data_source_split.pdf",pos_cell_umap12,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_major_annotation_data_source_split.png",pos_cell_umap12,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_major_annotation_merge.pdf",pos_cell_umap1,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_major_annotation_merge.png",pos_cell_umap1,width=6, height=6)

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_subtype_annotation_treat_split.pdf",pos_cell_umap23,width=16, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_subtype_annotation_treat_split.png",pos_cell_umap23,width=16, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_subtype_annotation_data_source_split.pdf",pos_cell_umap22,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_subtype_annotation_data_source_split.png",pos_cell_umap22,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_subtype_annotation_merge.pdf",pos_cell_umap2,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_subtype_annotation_merge.png",pos_cell_umap2,width=6, height=6)



pos_cell_umap112<-pos_cell_umap1+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                             plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                             axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                             axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
pos_cell_umap122<-pos_cell_umap12+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                              plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                              axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                              axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
pos_cell_umap132<-pos_cell_umap13+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                              plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                              axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                              axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
pos_cell_umap212<-pos_cell_umap2+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                             plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                             axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                             axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
pos_cell_umap222<-pos_cell_umap22+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                              plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                              axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                              axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
pos_cell_umap232<-pos_cell_umap23+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                              plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                              axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                              axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())


ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_major_annotation_treat_split_null.png",pos_cell_umap132,width=16, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_major_annotation_data_source_split_null.png",pos_cell_umap122,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_major_annotation_merge_null.png",pos_cell_umap112,width=6, height=6)

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_subtype_annotation_treat_split_null.png",pos_cell_umap232,width=16, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_subtype_annotation_data_source_split_null.png",pos_cell_umap222,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/GC_sample_all_sc_nature_intergrate_Group_umap_subtype_annotation_merge_null.png",pos_cell_umap212,width=6, height=6)

##plot ratio
metadata_Cell_Umap_add2<-subset(x = metadata_Cell_Umap_add, subset = major_group_brief != "Erythrocyte_HBA1")
major_order2<-c("Trophoblast_KRT7","Stromal_cells_DCN","Myeloid_Cell_AIF1","Endothelial_Cell_PECAM1","T_NK_B_Cell","Epithelial_Cell_EPCAM_PAEP")
conflict_prefer("which", "Matrix")
metadata_Cell_Umap_add2$major_group_brief<-factor(metadata_Cell_Umap_add2$major_group_brief,levels =major_order2)
metadata_Cell_Umap_add2$class<-metadata_Cell_Umap_add2$Treat
metadata_Cell_Umap_add2[which(metadata_Cell_Umap_add2$location =="Placenta"),]$class<- "Placenta"
metadata_Cell_Umap_add2[which(metadata_Cell_Umap_add2$location =="Decidua"),]$class<- "Decidua"
metadata_Cell_Umap_add2[which(metadata_Cell_Umap_add2$location =="Blood"),]$class<- "Blood"
metadata_Cell_Umap_add2$major_group_brief<-factor(metadata_Cell_Umap_add2$major_group_brief,levels =major_order2)
metadata_Cell_Umap_add2<-metadata_Cell_Umap_add2[which(metadata_Cell_Umap_add2$class %in% c("Placenta","Abortion","CTRL")),]

dim(metadata_Cell_Umap_add2)#63403    78
cell.prop1<-as.data.frame(prop.table(table(metadata_Cell_Umap_add2$major_group_brief, metadata_Cell_Umap_add2$class)))
colnames(cell.prop1)<-c("Cell_type","Group","proportion")
Ratio_01<-ggplot(cell.prop1,aes(Group,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+ theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=my_morandi_colors[c(2,1,3,9,6,7)])+coord_flip() 
Ratio_01
pieplot_ring3<-ggplot(cell.prop1,aes(x=Group,y=proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+coord_polar("y",start=1) + 
  theme_minimal()+ labs(title = "USM vs CTRL")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  #scale_fill_manual(values=rev(ppCor_all2[c(1:3,7,5:6)]))
  scale_fill_manual(values=my_morandi_colors[c(2,1,3,9,6,7)])#+facet_wrap(~Disease_group)
pieplot_ring3

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1f_major_Cell_ratio_for_GC_seurat_object_sc_nature_placenta_intergrate_Group_bar.pdf",Ratio_01,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1f_major_Cell_ratio_for_GC_seurat_object_sc_nature_placenta_intergrate_Group_ring_pie.pdf",pieplot_ring3,width=6, height=6)
