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

##add umap for GEX seurat
GEX_object0<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
GEX_object <-subset(x = GEX_object0,subset = final_major_subgroup== "Troph_mix_group",invert = TRUE)
metadata<-GEX_object@meta.data
final_umap<-data.frame(Embeddings(GEX_object[["umap"]]))
head(metadata_Cell);head(final_umap)
GEX_metadata_Umap_add <- merge(metadata,final_umap,by=0)
rownames(GEX_metadata_Umap_add)<-GEX_metadata_Umap_add$Row.names
head(GEX_metadata_Umap_add)
GEX_Umap<-GEX_metadata_Umap_add[,c("UMAP_1","UMAP_2")]
colnames(GEX_Umap)<-c("GEX_UMAP_1","GEX_UMAP_2")
#extract metadata information
metadata_Cell<-target_final@meta.data
TFs_metadata_GEX_Umap_add <- merge(metadata_Cell,GEX_Umap,by=0)
head(TFs_metadata_GEX_Umap_add)
rownames(TFs_metadata_GEX_Umap_add)<-TFs_metadata_GEX_Umap_add$Row.names

write.table(TFs_metadata_GEX_Umap_add, file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/supplement_Seurat_add_GEX_Seurat_UMAP_distribution.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

## calculating the position of cluster labels
get_label_pos <- function(data, emb = "UMAP", group.by="ClusterID") {
  new.data <- data[, c(paste(emb, 1:2, sep = "_"), group.by)]
  colnames(new.data) <- c("x","y","cluster")
  clusters <- names(table(new.data$cluster))
  new.pos <- lapply(clusters, function(i) {
    tmp.data = subset(new.data, cluster == i)
    data.frame(
      x = median(tmp.data$x),
      y = median(tmp.data$y),
      label = i)
  })
  do.call(rbind, new.pos)
}


##plot for major cell types in GEX style
pos_cell_umap0<-ggplot(data=TFs_metadata_GEX_Umap_add, mapping=aes(x=GEX_UMAP_1,y=GEX_UMAP_2,colour = sample_code ))+
  geom_point(stat= "identity",size=1,alpha=1,show.legend = TRUE)+
 # scale_color_manual(values=my_morandi_colors[c(2,1,3,9,6,7,8,4,5,10)])+
  scale_color_manual(values=my_morandi_colors)+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                   axis.text.y = element_text(size=15),legend.position = "right")
pos_cell_umap0
pos_cell_umap03<-pos_cell_umap0+facet_wrap(~sample_code,ncol = 3)
pos_cell_umap03

table(TFs_metadata_GEX_Umap_add$re_anno_TFac_major)
major_order<-c("Trophoblast_KRT7","Epithelial_Cell_EPCAM_PAEP","Endothelial_Cell_PECAM1","Myeloid_Cell_AIF1",
               "Stromal_cells_DCN","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
TFs_metadata_GEX_Umap_add$re_anno_TFac_major<- factor(TFs_metadata_GEX_Umap_add$re_anno_TFac_major, levels=major_order,ordered=TRUE)
pos_cell_umap1<-ggplot(data=TFs_metadata_GEX_Umap_add[order(TFs_metadata_GEX_Umap_add$re_anno_TFac_major ,decreasing = F),], 
                       mapping=aes(x=GEX_UMAP_1,y=GEX_UMAP_2,colour = re_anno_TFac_major ))+
  geom_point(stat= "identity",size=0.5,alpha=0.5,show.legend = TRUE)+
  #labs(title =paste0("Expression of ",genename))+
  scale_color_manual(values=my_morandi_colors[c(2,1,3,9,6,7,8,4,5,10)])+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                   axis.text.y = element_text(size=15),legend.position = "right")
pos_cell_umap1
pos_cell_umap13<-pos_cell_umap1+facet_wrap(~Treat)

subgroup_order<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
TFs_metadata_GEX_Umap_add$re_annotation_TFac<- factor(TFs_metadata_GEX_Umap_add$re_annotation_TFac, levels=subgroup_order,ordered=TRUE)
TFs_metadata_GEX_Umap_add<-TFs_metadata_GEX_Umap_add[order(TFs_metadata_GEX_Umap_add$re_annotation_TFac ,decreasing = F),]
pos_cell_umap2<-ggplot(data=TFs_metadata_GEX_Umap_add, 
                       mapping=aes(x=GEX_UMAP_1,y=GEX_UMAP_2,colour = re_annotation_TFac ))+
  geom_point(stat= "identity",size=0.5,alpha=1,show.legend = TRUE)+
  scale_color_manual(values=my_morandi_colors)+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
pos_cell_umap21<-pos_cell_umap2+geom_label_repel(inherit.aes = F, data = get_label_pos(TFs_metadata_GEX_Umap_add, emb = "GEX_UMAP",group.by="re_annotation_TFac"), aes(x,y,label=label),size=4, 
                                                 box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))
pos_cell_umap23<-pos_cell_umap21+facet_wrap(~Treat)
pos_cell_umap24<-pos_cell_umap2+facet_wrap(~sample_code,ncol = 3)

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_subtype_annotation_family_split.png",pos_cell_umap24,width=16, height=15)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_subtype_annotation_family_split.pdf",pos_cell_umap24,width=16, height=15)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_subtype_annotation_treat_split.pdf",pos_cell_umap23,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_subtype_annotation_treat_split.png",pos_cell_umap23,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_subtype_annotation_merge.pdf",pos_cell_umap21,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_subtype_annotation_merge.png",pos_cell_umap21,width=6, height=6)

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_major_annotation_treat_split.pdf",pos_cell_umap13,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_major_annotation_treat_split.png",pos_cell_umap13,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_major_annotation_merge.pdf",pos_cell_umap1,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_major_annotation_merge.png",pos_cell_umap1,width=6, height=6)

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_family_merge.pdf",pos_cell_umap0,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_family_merge.png",pos_cell_umap0,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_family_split.png",pos_cell_umap03,width=16, height=15)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_family_split.pdf",pos_cell_umap03,width=16, height=15)


pos_cell_umap111<-pos_cell_umap1+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                             plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                             axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                             axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
pos_cell_umap132<-pos_cell_umap13+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                              plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                              axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                              axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())

pos_cell_umap112<-pos_cell_umap2+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                             plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                             axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                             axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())

pos_cell_umap122<-pos_cell_umap2+facet_wrap(~Treat)+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                                                plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                                                axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                                                axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())


pos_cell_umap211<-pos_cell_umap03+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                              plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                              axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),legend.title = element_blank(),
                                                              axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_subtype_annotation_treat_split_no_text.png",pos_cell_umap122,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_subtype_annotation_merge_no_text.png",pos_cell_umap112,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_major_annotation_treat_split_no_text.png",pos_cell_umap132,width=13, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_major_annotation_merge.png",pos_cell_umap111,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1c_GEX_seurat_umap_family_split_no_text.png",pos_cell_umap211,width=16, height=15)

##plot maker genes in GEX distribution
#绘制Umap coordination
DefaultAssay(GEX_object) <- "RNA"
GEX_object <- NormalizeData(GEX_object, verbose = FALSE)

major_submaker0<-c("VIM","HLA-B","KRT7","PAGE4","PCNA","MKI67","EGFR","CDH1","CYP19A1","ERVFRD-1","HLA-G","PAPPA2","PAEP","PECAM1","DCN","DLK1","AIF1","LYVE1","NCAM1","CD3D","CD79A","MS4A2","PF4","HBA1")
##supplement figure XXX plot 
plot_umap_maker <- FeaturePlot(object = GEX_object, features = major_submaker0,cols= c("grey", "red"),ncol=5)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/S1a_GEX_seurat_Umap_final_major_submaker.png", plot_umap_maker,width=30, height=30)
plot_umap_maker <- FeaturePlot(object = GEX_object, features = major_submaker0,cols= c("grey", "red"),ncol=4)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/S1a_GEX_seurat_Umap_final_major_submaker2.png", plot_umap_maker,width=24, height=26)

##hightlight postive cells
target_gene<-major_submaker0
#target_gene<-Three_linege_maker
length(target_gene)#24
target_gene3<-target_gene[which(target_gene %in% rownames(GEX_object))]
length(target_gene3)#24

target_expression <- GetAssayData(object = GEX_object,assay = "RNA", slot = "data")[target_gene3,]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))
TFs_metadata_GEX_Umap_add2<-TFs_metadata_GEX_Umap_add[rownames(target_expression),]
merge_data2 <- merge(TFs_metadata_GEX_Umap_add2, target_expression,by=0,sort=FALSE)
head(merge_data2)
#绘制基因活性的 umap distribution
gene_list<-list();gene_list2<-list()
for ( genename in target_gene3){
  # genename <-"PAEP"
  print(genename)
  #for gene expression
  express_merge_data3<-merge_data2[order(merge_data2[,genename],decreasing = T),]
  express_merge_target<- express_merge_data3[,c("GEX_UMAP_1","GEX_UMAP_2",genename)]
  colnames(express_merge_target)<-c("GEX_UMAP_1","GEX_UMAP_2","Expression")
  head(express_merge_target)
  express_merge_target$selected_cluster<-"other"
  express_merge_target[which(express_merge_target$Expression>0),]$selected_cluster<-"Pos"
  express_merge_target$selected_cluster <-factor(express_merge_target$selected_cluster,levels=c("Pos","other"))
  pos_cell_umap<-ggplot(data=express_merge_target[order(express_merge_target$selected_cluster,decreasing = T),],
                        mapping=aes(x=GEX_UMAP_1,y=GEX_UMAP_2,color = Expression))+
    geom_point(stat= "identity",size=0.2,alpha=0.8,show.legend = TRUE)+
    scale_color_gradient2(low="#4393C3", mid="pink", high="red", midpoint=max(express_merge_target$Expression)/2)+#"#156077","#f46f20"
    labs(title =paste0("Expression of ",genename))+theme_bw()+
    theme(panel.border = element_rect(colour="grey",fill=NA),
          panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
          axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=15),legend.position = "right")
  pos_cell_umap
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/GEX_umap_expression_",genename,"_pos_value.png"),pos_cell_umap,width=8, height=6)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/GEX_umap_expression_",genename,"_pos_value.pdf"),pos_cell_umap,width=8, height=6)
  gene_list<-c(gene_list,list(pos_cell_umap))
  
  pos_cell_umap2<-pos_cell_umap+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/single_gene_selected/GEX_umap_expression_",genename,"_pos_value_notext.png"),pos_cell_umap2,width=8, height=6)
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
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1a_GEX_seurat_target_major_umap_expression2.png",gene_list_merge,width=24, height=26)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/SF1a_GEX_seurat_target_major_umap_expression_notext2.png",gene_list_merge2,width=24, height=26)
