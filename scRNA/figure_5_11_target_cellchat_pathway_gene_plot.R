rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(data.table)
library(Seurat)
library(CellChat)
library(patchwork)
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


##read all genes provided by pathway recorded in CellChatDB.human of CellChat
pathyway_CPIpair<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/cellchat_pathyway_CPIpair.rds")
CellChatDB <- CellChatDB.human 
CellChatDB_human_interaction_meta<- data.frame(CellChatDB$interaction)[,c("pathway_name","interaction_name","ligand","receptor")]
pathways_all <- unique(CellChatDB_human_interaction_meta$pathway_name)

pathyway_gene<-readRDS( file = "/mnt/data/chenwei/gongchen/7.cellchat/cellchat_pathyway_gene.rds")
pathway_name<-"WNT"
pathyway_gene[[pathway_name]]

List_tag<-paste0("USM_interaction_gene_of_",pathway_name)

CellChatDB_target<-CellChatDB_human_interaction_meta[which(CellChatDB_human_interaction_meta$pathway_name ==pathway_name),]
CellChatDB_target[which(CellChatDB_target$ligand =="WNT7A"),]
#############
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
##calculation 
DefaultAssay(target_final) <- "RNA"
# Normalize RNA data for visualization purposes
target_final <- NormalizeData(target_final, verbose = FALSE)

####get the detected genes
Detect_path_genes<-Reduce(intersect,list(pathyway_gene[[pathway_name]],rownames(target_final)))
Detect_path_genes<-Detect_path_genes[order(Detect_path_genes)]
length(Detect_path_genes)#25
gene_signature<-Detect_path_genes


##4.整体细胞类群表达水平绘制
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts")
target_final2<-subset(x = target_final, subset = re_annotation_TFac %in% subgroup_order0)
target_final2$re_annotation_TFac<- factor(target_final2$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)

Idents(object = target_final2) <- "re_annotation_TFac"
target_final2$re_annotation_TFac<- factor(target_final2$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)

levels(target_final2)
DefaultAssay(target_final2)<-"RNA"

#其他SCT的数据说明：https://www.jianshu.com/p/e639cc257d51 
mat0<- GetAssayData(target_final2, slot = "data")

#get genes and cluster information
gene_features <- unique(gene_signature);length(gene_signature)
cluster_info <- target_final2$re_annotation_TFac

#sort the expression matrix 
mat1 <- as.matrix(mat0[as.character(gene_features), names(cluster_info)])
mat1[1:6,1:6];dim(mat1)
name_mean1<-as.data.frame(t(unlist(apply(mat1,1,mean))))
name_mean1<-as.data.frame(rbind(c(1:length(rownames(mat1))),name_mean1))
name_mean1<-as.data.frame(rbind(colnames(name_mean1),name_mean1))
rownames(name_mean1)<-c("Gname","row","mean")
dim(name_mean1);str(name_mean1)
name_mean1<-as.data.frame(t(name_mean1))
name_mean1$row<-as.numeric(as.character(name_mean1$row))
name_mean1$mean<-as.numeric(as.character(name_mean1$mean))
name_mean1$Gname<-as.character(name_mean1$Gname)
dim(name_mean1)
head(name_mean1)
#去重
name_mean1<-name_mean1[order(name_mean1$Gname,name_mean1$mean,decreasing = T),]
name_mean1<-name_mean1[!duplicated(name_mean1$Gname),]
name_mean1<-name_mean1[order(name_mean1$row),]
length(unique(gene_features))#25
#name_mean1$row
head(name_mean1)
gene_features_uniq<-name_mean1$Gname
gene_features_dup<-gene_features[which(duplicated(gene_features))]
length(gene_features);length(gene_features_uniq);length(gene_features_dup)
#重新去unique数据框
mat2 <- mat1[gene_features_uniq,]
#手动scale
mat3<-scale(t(mat2), center=T,scale=T)
range(mat3)#   -0.4667643 173.2193591
mat4<-t(mat3)

#mat4[mat4>=2.5]= 2.5
#mat4[mat4 < (-2.5)]= -2.5 #小于负数时，加括号！

dim(mat4)# 25 44706
#set color for cell types
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60,43)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
####col <- ppCor_all[1:length(levels(cluster_info))]
col <-my_morandi_colors[1:16]
#colourCount<-length(unique(target$sub_cluster))
names(col) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col),labels = levels(cluster_info), labels_gp = gpar(cex = 0.5, col = "black")))

#important genes notes
#mark_gene0<-c("PCNA","PAGE4","CDH1","ERVFRD-1","CYP19A1","EGFR","MKI67","HLA-G","PAPPA2","PAEP","PECAM1","LYVE1","AIF1","DCN","DLK1","NCAM1","CD3D","CD79A","MS4A2","PF4","HBA1")
#gene_pos <- which(rownames(mat4) %in% mark_gene0)
#mark_gene1<-rownames(mat4)[gene_pos]
#row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = mark_gene1))

#legend adjust
col_fun1  <- circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red"))
col_fun2  <- circlize::colorRamp2(c(-2,0,2), c("purple", "black", "yellow"))
col_fun<-col_fun1
#show_col(col_fun(seq(-3, 3)))

plot_HP<-Heatmap(mat4,cluster_rows = FALSE,cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = TRUE,
                 #show_row_names = TRUE,
                 column_split = cluster_info,
                 top_annotation = top_anno,
                 # right_annotation = row_anno,
                 column_title = NULL,
                 col = col_fun,
                 heatmap_legend_param = list(legend_height = unit(4, "cm"),grid_width = unit(0.4, "cm"),
                                             title = "log2(count+1)",title_position = "leftcenter-rot",at = c(-2, 0, 2),
                                             labels = c("low", "median", "high") ))
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/final_all_19_intercuster_all_DEGs_heatmap.pdf",width = 10,height =8)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("Pathway_heatmap_RNA_data"))
dev.off()

png("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_1/intercluster_maker/final_all_19_intercuster_all_DEGs_heatmap.png",height=1000,width=1000)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("major_submaker_heatmap_RNA_data"))
dev.off()


##绘制分组点图
##plot dot
target_final$cell_class<-paste0(target_final$re_annotation_TFac,"_",target_final$Treat)

##正式绘图
Plotdot <-DotPlot(target_final, features = gene_signature,group.by ="re_annotation_TFac",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()
Plotdot2 <-DotPlot(target_final, features = gene_signature,group.by ="cell_class",cols =  c("blue","red"),assay = "RNA") #+ RotatedAxis()
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
Plotdot31
Plotdot32
Plotdot02
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_genes_dotplot_all_celltype_split.pdf"),Plotdot31,width=4, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_genes_dotplot_all_celltype_split2.pdf"),Plotdot32,width=4, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/",Group_tag,"_genes_dotplot_all_celltype_split3.pdf"),Plotdot02,width=4, height=4)


Idents(object = target_final) <- "re_annotation_TFac"
plots2 <- VlnPlot(target_final, features = rev(gene_signature), split.by = "Treat", group.by = "re_annotation_TFac", split.plot = TRUE,pt.size = 0, combine = FALSE)
plots3<-CombinePlots(plots = plots2, ncol =1)

VlnPlot(target_final, features = gene_signature,group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 

violin_plot<-VlnPlot(target_final, features = gene_signature, group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],combine =F) 
gene_number<-length(gene_signature)
plot_gg <- list()
for (j in 1:gene_number) {
  # j=4
  plot_gg[[j]]<- violin_plot[[j]]+ stat_compare_means(method = "wilcox.test",data=violin_plot[[j]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =max(violin_plot[[j]]$data[,1]-0.25))+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
}
gene_list_merge<-patchwork::wrap_plots(plots = plot_gg,ncol=1)
plot_gg


###subset target expression matrix
exprs <- data.frame(FetchData(object = target_final,slot = "data",vars = c(gene_signature,"Treat","re_annotation_TFac")))
#conflict_prefer("melt", "reshape2")
data_plot <- melt(exprs)
head(data_plot);range(data_plot$value)

###
#cell_group<-"merge_cluster";split_class<-"merge_cluster"
cell_group<-"re_annotation_TFac";split_class<-"Treat"

##remove zero expression
#conflict_prefer("which", "Matrix")
data_plot2<-data_plot[which(data_plot$value>0),]
#maker_plot<-ggplot(data_plot2,aes(x = merge_cluster,y=value,fill = merge_cluster)) +
maker_plot<-ggplot(data_plot2,aes(x = re_annotation_TFac,y=value,fill = Treat)) +  
  # geom_point(size = 1,position = 'jitter',alpha = 0.1)+  
  geom_violin(aes(fill = Treat),scale = "width",alpha = 0.7) +  
  facet_wrap(~variable,scales = "free_y", ncol = 5) +
  scale_fill_manual(values=ppCor[c(3,9)])+  theme_bw()+ 
  labs(x = "Group", y = "LogNormalized count", title ="SNPs correlated TFs")+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x.bottom = element_blank(),legend.position = "right")+
  stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))

maker_plot
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_sub_and_split_",Group_tag,"_violin_target_gene_plot.pdf"),maker_plot,width=12, height=12)


##For WNT7A
genename <-"FZD6"#WNT7A
print(genename)
merge_data_target<- data_plot2[which(data_plot2$variable==genename),]
head(merge_data_target);dim(merge_data_target)#17723     5
merge_data_target_original<-merge_data_target
range(merge_data_target_original$value)# 0.1479054 4.6395716
max_value<-2.5
merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot3<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  scale_y_continuous(limits = c(0,2.6),breaks = c(0,0.5,1,1.5,2,2.5,3))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot3
#ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_sub_and_split_max_value_adjust_",genename,"_violin_target_gene_plot.pdf"),maker_plot3,width=5, height=4)
