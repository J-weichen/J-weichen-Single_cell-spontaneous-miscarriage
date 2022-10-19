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
###KS_panc_fl_data
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

DefaultAssay(target_final) <- "RNA"
# Normalize RNA data for visualization purposes
target_final <- NormalizeData(target_final, verbose = FALSE)
express_data <- as.matrix(GetAssayData(target_final, slot = "data"))

Idents(object = target_final) <- "re_annotation_TFac"
##plot genes in selected final pathway
# Access all the signaling pathways showing significant communications
pathyway_CPIpair<-readRDS(file = "/mnt/data/chenwei/gongchen/7.cellchat/cellchat_pathyway_CPIpair.rds")
CellChatDB <- CellChatDB.human 
CellChatDB_human_interaction_meta<- data.frame(CellChatDB$interaction)[,c("pathway_name","interaction_name","ligand","receptor")]
pathways_all <- unique(CellChatDB_human_interaction_meta$pathway_name)

Group_tag<-"Final_selected_pathway_WNT"

gene_signature<-unique(as.character(unlist(lapply(strsplit(as.character(pathyway_CPIpair[[pathways.show.all[i]]]),"_"), function(x) x))))

violin_plot<-VlnPlot(target_final, features = gene_signature, group.by = "re_annotation_TFac",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], y.max =4, combine =F) 

gene_number<-length(violin_plot[]$patches$plots)
plot_gg <- list()
for (j in 1:gene_number) {
  # j=4
  plot_gg[[j]]<- violin_plot[[j]]+ stat_compare_means(method = "wilcox.test",data=violin_plot[[j]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =max(violin_plot[[j]]$data[,1]-0.25))+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))
}
#gene_names<-unique(as.character(unlist(lapply(strsplit(as.character(pathyway_CPIpair[[pathways.show.all[i]]]),"_"), function(x) x))))
gene_list_merge<-patchwork::wrap_plots(plots = plot_gg,ncol=1)
ggsave(filename=paste0("/mnt/data/chenwei/gongchen/7.cellchat/Merge/gene_expression/Merge_",pathways.show.all[i], "_gene_expression_violin_plot.pdf"), plot=gene_list_merge, width = 12, height =gene_number*1.5,limitsize = FALSE)



Group_tag<-"Abortion_SNP"
#panc_fl <- subset(x = KS_raw, idents = c("Common Progenitor","Leydig Progenitor","Leydig Precursor #1","Leydig Precursor #2","Differented Leydig Cells", "Peritubular Myoid Cells", "Sertoli Cells #1", "Sertoli Cells #2", "Sertoli Progenitor"))

umapplot<-DimPlot(target_final, group.by = "re_annotation_TFac", label = T, repel = T,cols = ppCor_all2)
#ggsave("/mnt/data/chenwei/qinmen_BR/keshi/Umap_cell_type_plot.pdf", umapplot,width=8, height=6)
#ggsave("/mnt/data/chenwei/qinmen_BR/keshi/Umap_cell_type_plot.png", umapplot,width=8, height=6)


##The normalized and log-transformed values are used for the violin plot. 
##VlnPlot doesn't perform any additional transformations on the data. It will just plot what you have stored in @data. 
##You can verify this for yourself if you want by pulling the data out manually and inspecting the values.
##The argument y.log changes only the display of the data (scaling of the y axis).
#data: The data slot (object@data) stores normalized and log-transformed single cell expression. 
#This maintains the relative abundance levels of all genes, and contains only zeros or positive values. 
#See ?NormalizeData for more information. This data is used for visualizations, such as violin and feature plots, most differential expression tests, finding high-variance genes, and as input to ScaleData (see below).
gene_signature<-c("VEGFA","MTHFR")
gene_signature<-c("NR5A1","SOX4","STAT3","SERPINE1","CD320")
gene_signature<-c("SMCHD1","TNFRSF1B","TMEM176A")
gene_signature<-c("CCR7","MALAT1","SERPINE1","STAT3","TGFB1")
target_violin<-VlnPlot(target_final, features =gene_signature, pt.size = 0.01, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 1)
ggsave("/mnt/data/chenwei/qinmen_BR/keshi/sixviolin_target_gene_plot.pdf", sixviolin,width=16, height=11)

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


#extract metadata information
metadata_Cell<-target_final@meta.data
final_umap<-data.frame(Embeddings(target_final[["umap"]]))
#final_umap$raw_id <- rownames(final_umap)
head(metadata_Cell);head(final_umap)
metadata_Cell_Umap_add <- merge(metadata_Cell,final_umap,by=0)
head(metadata_Cell_Umap_add)
rownames(metadata_Cell_Umap_add)<-metadata_Cell_Umap_add$Row.names

##plot for target genes
target_final$cell_class<-paste0(target_final$cell_type,"_",target_final$Sex)
table(target_final$cell_class)

##get genes in target 
target_gene<-c("GLI1","PTCH1","SF1","INSL3","CYP17A1","IGF2")
tag<-"qingmeng_genes"

##正式绘图
Plotdot <-DotPlot(target_final, features = target_gene,group.by ="cell_type",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()
Plotdot2 <-DotPlot(target_final, features = target_gene,group.by ="cell_class",cols =  c("blue","red"),assay = "RNA") #+ RotatedAxis()
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
Plotdot01<-Plotdot+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ 
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1, vjust=0.5))
Plotdot02<-Plotdot2+scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+ 
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1, vjust=0.5))

ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_",tag,"_genes_dotplot_all_celltype_split.pdf"),Plotdot31,width=8, height=5)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_",tag,"_genes_dotplot_all_celltype_split_2.pdf"),Plotdot32,width=8, height=5)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_",tag,"_genes_dotplot_all_celltype_1.pdf"),Plotdot01,width=7, height=6)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_",tag,"_genes_dotplot_all_celltype_split_3.pdf"),Plotdot02,width=7, height=8)

##plot violin with p value
target_gene2<-target_gene[which(target_gene %in% rownames(target_final))]
target_expression <- GetAssayData(object = target_final,assay = "RNA", slot = "data")[target_gene2,]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))

metadata_Cell_Umap_add2<-metadata_Cell_Umap_add[rownames(target_expression),]
merge_data2 <- merge(metadata_Cell_Umap_add2, target_expression,by=0,sort=FALSE)
head(merge_data2)

##绘制两比较组间的表达量
gene_list<-list()
for ( genename in rev(target_gene2)){
  # genename <-"PTCH1"
  print(genename)
  merge_data_target<- merge_data2[,c("cell_type","Sex",genename)]
  colnames(merge_data_target)<-c("cell_type","Sex","Expression_level")
  plot_Cell_exp<- ggplot(merge_data_target,aes(x = cell_type,y=Expression_level,fill = Sex)) + 
    geom_violin(scale = "width") +  #geom_jitter(size = 0.00)+coord_flip()+ facet_grid(. ~ variable,scale='free') +
    theme_bw()+scale_fill_manual(values=ppCor)+
    labs(y = "Expression level", x = "main_Cell_type", title = genename)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.y.left = element_blank())+
    theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),legend.position = "right")+
    theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))+
    stat_compare_means(method = "wilcox.test",data=merge_data_target,aes(x=cell_type,y = Expression_level,group = Sex,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((merge_data_target$Expression_level)))
  plot_Cell_exp2<- plot_Cell_exp+ geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=4, col="orange", position = position_dodge(0.9))
  ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/six_",Group_tag,"_stat_barplot_expression_",tag,"_",genename,".pdf"),plot_Cell_exp2,width=12, height=7)
  gene_list<-c(gene_list,list(plot_Cell_exp2))
}

gene_list_merge<-grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],gene_list[[5]], gene_list[[6]],ncol=3)
ggsave(paste0("/mnt/data/chenwei/qinmen_BR/keshi/merge_six_",Group_tag,"_stat_barplot_expression_",tag,"_",genename,".pdf"),gene_list_merge,width=33, height=12)


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
