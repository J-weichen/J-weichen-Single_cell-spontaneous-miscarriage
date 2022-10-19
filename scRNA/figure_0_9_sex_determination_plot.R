rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
dyn.load('/mnt/data/chenwei/software/hdf/hdf5/lib/libhdf5.so.200')
dyn.load('/mnt/data/chenwei/software/hdf/hdf5/lib/libhdf5_hl.so.200')

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
subgroup_order<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order,ordered=TRUE)

#绘制Umap coordination
DefaultAssay(target_final) <- "RNA"
target_final <- NormalizeData(target_final, verbose = FALSE)

HBs_gene<- c("CD14","FCGR3A","FCGR3B","ITGAX","ACP3","CD1A","CD1B","CD1C","CD1D","CCR5","CD68","CR3","CXCR4","G6PD","IFNA1","IFNA3","IFNA8")
plot_HBs_maker <- FeaturePlot(object = target_final, features = HBs_gene,cols= c("grey", "purple"),ncol=4)
plot_HBs_maker
#  The following requested variables were not found: ACP3, CR3, IFNA3, IFNA8
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/HBs_maker_for_intergrated_nine_Cell.png", plot_HBs_maker,width=20, height=20)

#read Y special maker
Y_special_makers <- read.table("/mnt/data/chenwei/gongchen/0.script/Y_special_gene.txt",header = T)
Y_special_gene<-unique(as.character(Y_special_makers$Y_speical_gene))
gene_target<-intersect(Y_special_gene,rownames(target_final[["RNA"]]))
plot_Y_special_gene<-FeaturePlot(object = target_final, features = c(gene_target),cols= c("grey", "red"),ncol=5)
plot_Y_special_gene
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/Y_special_gene_umap.png",plot_Y_special_gene,width=25, height=20)

Sex_gene<- c('RBMY2FP','RBMY1B','TTTY15','KDM5D','RBMY1J','RBMY1F','RBMY1D','RBMY1E','TSPY4','TSPY2','TSPY8','TSPY1','RPS4Y1','EIF1AY','DDX3Y')
plot_YP_sex_gene<-FeaturePlot(object = target_final, features = c(Sex_gene),cols= c("grey", "purple"),ncol=3)
plot_YP_sex_gene
#  The following requested variables were not found: RBMY2FP, RBMY1B, RBMY1J, RBMY1D, TSPY4, TSPY2, TSPY8, TSPY1, RPS4Y1
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/YP_sex_gene_umap.png",plot_YP_sex_gene,width=15, height=10)

target_final$Y_genesum<-colSums(x = GetAssayData(object = target_final,assay = "RNA", slot = "counts")[gene_target, , drop = FALSE])
plot_EIF1AY_XIST_gene<-FeaturePlot(object = target_final, features = c('EIF1AY',"XIST"),cols= c("grey", "purple"),ncol=2)
plot_EIF1AY_XIST_gene
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/EIF1AY_XIST_gene_umap.png",plot_EIF1AY_XIST_gene,width=20, height=11)

#FeaturePlot(object = target, features = c("XIST"),cols= c("grey", "purple"))
P0<-FeaturePlot(object = target_final, features = c("Y_genesum"),cols= c("grey", "blue"))
P1<-FeaturePlot(object = target_final, features = c("Y_genesum"),cols= c("grey", "blue"),max.cutoff = 5)
P2<-FeaturePlot(object = target_final, features = c('EIF1AY'),cols= c("grey", "red"),max.cutoff = 2)
P3<-FeaturePlot(object = target_final, features = c("XIST"),cols= c("grey", "purple"),max.cutoff = 3)
plot_three_gene<-CombinePlots(plots = list(P1,P2,P3),legend="top",ncol=2)
plot_three_gene
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/three_sex_maker_umap.png",plot_three_gene,width=30, height=11)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/Y_sum_no_limit_sex_maker_umap.png",P0,width=10, height=10)

#one:EIF1AY
target_final$selected_cluster <- rep("EIF1AY_neg",nrow(target_final@meta.data))
sub_cell<-subset(x = target_final,subset = EIF1AY >0) 
target_final$selected_cluster[Cells(sub_cell)] <- "EIF1AY_pos"
prop.table(table(target_final$selected_cluster))
#EIF1AY_neg EIF1AY_pos 
#0.7212009  0.2787991  

metadata_for_target<-target_final@meta.data
target_coordinate<-data.frame(Embeddings(target_final[["umap"]]))
target_coordinate$raw_id<-rownames(target_coordinate)
target_coordinate2 <- merge(metadata_for_target,target_coordinate,by="raw_id")
target_coordinate2$selected_cluster <-factor(target_coordinate2$selected_cluster,levels=c("EIF1AY_pos","EIF1AY_neg"))
head(target_coordinate2)

selected_cluster_umap1<-ggplot(data=target_coordinate2[order(target_coordinate2$selected_cluster,decreasing = T),], 
                               mapping=aes(x=UMAP_1,y=UMAP_2,colour = selected_cluster))+
  geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
  labs(title ="Distribution of EIF1AY postive Cell")+scale_color_manual(values=c(ppCor_all[2],"grey"))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                   axis.text.y = element_text(size=15),legend.position = "right")
selected_cluster_umap1#+facet_grid(.~Age_group,margins=TRUE)

#two XIST
target_final$selected_cluster <- rep("XIST_neg",nrow(target_final@meta.data))
sub_cell<-subset(x = target_final,subset = XIST >0) 
target_final$selected_cluster[Cells(sub_cell)] <- "XIST_pos"
prop.table(table(target_final$selected_cluster))
# XIST_neg  XIST_pos 
#0.5946507 0.4053493 

metadata_for_target<-target_final@meta.data
target_coordinate<-data.frame(Embeddings(target_final[["umap"]]))
target_coordinate$raw_id<-rownames(target_coordinate)
target_coordinate2 <- merge(metadata_for_target,target_coordinate,by="raw_id")
target_coordinate2$selected_cluster <-factor(target_coordinate2$selected_cluster,levels=c("XIST_pos","XIST_neg"))
head(target_coordinate2)

selected_cluster_umap2<-ggplot(data=target_coordinate2[order(target_coordinate2$selected_cluster,decreasing = T),], 
                               mapping=aes(x=UMAP_1,y=UMAP_2,colour = selected_cluster))+
  geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
  labs(title ="Distribution of XIST postive Cell")+scale_color_manual(values=c(ppCor_all[12],"grey"))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                   axis.text.y = element_text(size=15),legend.position = "right")
selected_cluster_umap2#+facet_grid(.~Age_group,margins=TRUE)


#three Y gene sum
target_final$selected_cluster <- rep("Y_gene_neg",nrow(target_final@meta.data))
sub_cell<-subset(x = target_final,Y_genesum>0) 
target_final$selected_cluster[Cells(sub_cell)] <- "Y_gene_pos"
prop.table(table(target_final$selected_cluster))
#Y_gene_neg Y_gene_pos 
#0.6201965  0.3798035 

metadata_for_target<-target_final@meta.data
target_coordinate<-data.frame(Embeddings(target_final[["umap"]]))
target_coordinate$raw_id<-rownames(target_coordinate)
target_coordinate2 <- merge(metadata_for_target,target_coordinate,by="raw_id")
target_coordinate2$selected_cluster <-factor(target_coordinate2$selected_cluster,levels=c("Y_gene_pos","Y_gene_neg"))
head(target_coordinate2)

selected_cluster_umap3<-ggplot(data=target_coordinate2[order(target_coordinate2$selected_cluster,decreasing = T),], 
                               mapping=aes(x=UMAP_1,y=UMAP_2,colour = selected_cluster))+
  geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
  labs(title ="Distribution of Y gene sum postive Cell")+scale_color_manual(values=c(ppCor_all[3],"grey"))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                   axis.text.y = element_text(size=15),legend.position = "right")
selected_cluster_umap3#+facet_grid(.~Age_group,margins=TRUE)

plot_sex_maker<-selected_cluster_umap1+selected_cluster_umap2+selected_cluster_umap3
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/pos_three_sex_maker_umap.png",plot_sex_maker,width=30, height=11)


##expression highlight
##hightlight postive cells
conflict_prefer("which", "Matrix")
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
subgroup_order<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order,ordered=TRUE)

#绘制Umap coordination
DefaultAssay(target_final) <- "RNA"
target_final <- NormalizeData(target_final, verbose = FALSE)
#read Y special maker
Y_special_makers <- read.table("/mnt/data/chenwei/gongchen/0.script/Y_special_gene.txt",header = T)
Y_special_gene<-unique(as.character(Y_special_makers$Y_speical_gene))
gene_target<-intersect(Y_special_gene,rownames(target_final[["RNA"]]))
target_final$Y_genesum<-colSums(x = GetAssayData(object = target_final,assay = "RNA", slot = "counts")[gene_target, , drop = FALSE])
target_final$Y_genesum_log2<-log2(target_final$Y_genesum+1)
range(target_final$Y_genesum)#0 54
range(target_final$Y_genesum_log2)#0.00000 5.78136

metadata_for_target<-target_final@meta.data
target_coordinate<-data.frame(Embeddings(target_final[["umap"]]))
target_coordinate$raw_id<-rownames(target_coordinate)
metadata_Cell_Umap_add <- merge(metadata_for_target,target_coordinate,by="raw_id")
head(metadata_Cell_Umap_add)
rownames(metadata_Cell_Umap_add)<-metadata_Cell_Umap_add$raw_id 

target_expression <- GetAssayData(object = target_final,assay = "RNA", slot = "data")[c('EIF1AY',"XIST"),]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))
dim(target_expression)

metadata_Cell_Umap_add2<-metadata_Cell_Umap_add[rownames(target_expression),]
merge_data2 <- merge(metadata_Cell_Umap_add2, target_expression,by=0,sort=FALSE)
dim(merge_data2)
head(merge_data2)
write.table(as.data.frame(merge_data2), "/mnt/data/chenwei/gongchen/manuscript/manu_table/final_annotation_for_scRNAseq_cell.txt",quote = F, row.names = T,sep = "\t")


#绘制基因活性的 umap distribution
target_gene3<-c("Y_genesum","Y_genesum_log2",'EIF1AY',"XIST")
gene_list<-list();gene_list2<-list()
for ( genename in target_gene3){
  # genename <-"Y_genesum_log2"
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
    geom_point(stat= "identity",size=0.5,alpha=0.5,show.legend = TRUE)+
    scale_color_gradient2(low="#4393C3", mid="pink", high="red", midpoint=max(express_merge_target$Expression)/2)+#"#156077","#f46f20"
    labs(title =paste0("Expression of ",genename))+theme_bw()+
    theme(panel.border = element_rect(colour="grey",fill=NA),
          panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
          axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=15),legend.position = "right")
  pos_cell_umap
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/umap_expression_",genename,"_pos_value.png"),pos_cell_umap,width=8, height=6)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/umap_expression_",genename,"_pos_value.pdf"),pos_cell_umap,width=8, height=6)
  gene_list<-c(gene_list,list(pos_cell_umap))
  
  pos_cell_umap2<-pos_cell_umap+theme(plot.title=element_blank(),axis.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank())
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/umap_expression_",genename,"_pos_value_notext.png"),pos_cell_umap2,width=8, height=6)
  gene_list2<-c(gene_list2,list(pos_cell_umap2))
}

gene_list_merge<-grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],ncol=2)
gene_list_merge2<-grid.arrange(gene_list2[[1]], gene_list2[[2]],gene_list2[[3]],gene_list2[[4]], ncol=2)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/SF1a_final_target_major_umap_expression2.png",gene_list_merge,width=12, height=10)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/SF1a_final_target_major_umap_expression_notext2.png",gene_list_merge2,width=12, height=10)

gene_list_merge<-grid.arrange(gene_list[[2]], gene_list[[4]],ncol=2)
gene_list_merge2<-grid.arrange(gene_list2[[2]],gene_list2[[4]], ncol=2)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/SF1a_final_target_major_umap_expression3.png",gene_list_merge,width=11, height=5)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/sex_determination/SF1a_final_target_major_umap_expression3_notext2.png",gene_list_merge2,width=11, height=5)
