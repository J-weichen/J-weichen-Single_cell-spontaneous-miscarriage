rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
options(stringsAsFactors = FALSE)
## Load packages
##dyn.load('/usr/local/hdf5/lib/libhdf5_hl.so.100')
library(Seurat)
library(tidyverse)
library(reshape2)
library(dendextend)
library(ggsci)
library(scales)

##sample information plot in generation
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
#annalysis about acTFs
#read Seurat object and load initial cAbortionata and matrix
#target_final<-readRDS(file = "/home/chenwei/10x_data/191125-merge_nonormalization_add/Cell_0512.rds")
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
                     
target_group="re_annotation_TFac"
unique(target_final[[target_group]][,target_group])
metadata_for_all_cells<-target_final@meta.data
#metadata_for_all_cells2 <- merge(metadata_for_all_cells,Umap_coordinate,by=0,sort=FALSE)

## 读入RAS矩阵
cell_name="ALLcell_all_genes"
rasMat<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".rasMat.rds"))
rasMat<-rasMat[rownames(metadata_for_all_cells),]
dim(rasMat)
merge_data2 <- merge(metadata_for_all_cells, rasMat,by=0,sort=FALSE)
merge_data2$re_annotation_TFac_Treat <-  paste(merge_data2$re_annotation_TFac, merge_data2$Treat, sep="_")
#读入细胞信息
##数据集合并

#确定目标基因
gene_target<-c("ETV5","TEAD4","NFE2L3","TFAP2C","NFYB","TBPL1","HEY1","BHLHE40","SOX4","MYCN","ELF3","ELF1","TFAP2A","NFIL3")
Tag<-"common_TFs_DEGs"

gene_target<-c("HBZ","HBE1","HBG2","GAPDH","HBA1","XAGE3","COL1A1","HBA2","LY6E","PKM","APOE","APOC1","CSH2","CTSD","SPP1","TFPI2","CSH1","FCER1G","HLA-A","HLA-DRA")
Tag<-"Top10_DEGs"

gene_target<-c("POU2AF1","PAX4","XBP1","BHLHE40","ETV5","EZH2","YY2","ZBTB2","CEBPB","FOS","ZNF444","EGR1", "KLF4","SOX21","SP2","HIF1A","HOXC5","ELF1","ETV4","FOSB" )
Tag<-"Top10_TFs"

#conflict_prefer("melt", "reshape2");conflict_prefer("which", "Matrix")


gene_target<-c("MTHFR","VEGFA")
Tag<-"SNP_overlap_gene"

length(gene_target)
#for gene expression
express_data2<-data.frame(t(express_data[gene_target,]))
express_merge_data2 <- merge(metadata_for_all_cells, express_data2,by=0,sort=FALSE)
gene_target2<-gsub("-",".",gene_target)
express_merge_data3 <- express_merge_data2[,c("re_annotation_TFac","Treat",gene_target2)]
express_merge_data4<-melt(express_merge_data3, id=c("re_annotation_TFac","Treat"))
colnames(express_merge_data4)<-c("Cell_type","Treat","gene","Expression")
head(express_merge_data4)

#for TFs regulon activity
gene_target[which(gene_target %in% colnames(rasMat))]
colnames(merge_data2)
merge_data3 <- merge_data2[,c("re_annotation_TFac","Treat",gene_target[which(gene_target %in% colnames(rasMat))])]
head(merge_data3)
merge_data4<-melt(merge_data3, id=c("re_annotation_TFac","Treat"))
colnames(merge_data4)<-c("Cell_type","Treat","gene","Activity")
head(merge_data4)
merge_data4$gene<-gsub("-",".",merge_data4$gene)

#绘制基因活性的 umap distribution
for ( genename in gene_target2){
  # genename <-"MTHFR"
  print(genename)
  
  #for gene expression
  #for TFs activity
  express_merge_data3<-express_merge_data2[order(express_merge_data2[,genename],decreasing = T),]
  express_merge_target<- express_merge_data3[,c("UMAP_1","UMAP_2","Treat",genename)]
  colnames(express_merge_target)<-c("UMAP_1","UMAP_2","Treat","Expression")
  
  head(express_merge_target)
  express_merge_target$selected_cluster<-"other"
  express_merge_target[which(express_merge_target$Expression>0),]$selected_cluster<-"Pos"
  express_merge_target$selected_cluster <-factor(express_merge_target$selected_cluster,levels=c("Pos","other"))
  pos_cell_umap1<-ggplot(data=express_merge_target[order(express_merge_target$selected_cluster,decreasing = T),], 
                         mapping=aes(x=UMAP_1,y=UMAP_2,colour = selected_cluster))+
    geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
    labs(title =paste0("Expression of ",genename))+scale_color_manual(values=c(ppCor_all[3],"grey"))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                     panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                     axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                     axis.text.y = element_text(size=15),legend.position = "right")
  pos_cell_umap1+facet_grid(.~Treat)
  pos_cell_umap2<-ggplot(data=express_merge_target[order(express_merge_target$selected_cluster,decreasing = T),],
                         mapping=aes(x=UMAP_1,y=UMAP_2,color = Expression))+
    geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
    scale_color_gradient(low = "lightgrey",high = "red")+
    labs(title =paste0("Expression of ",genename))+theme_bw()+
    theme(panel.border = element_rect(colour="grey",fill=NA),
          panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
          axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=15),legend.position = "right")+facet_grid(.~Treat)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/",Tag,"_Expression_",genename,"_pos_highlight.pdf",sep=""),pos_cell_umap1,width=9, height=8)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/",Tag,"_Expression_",genename,"_value.pdf",sep=""),pos_cell_umap2,width=12, height=6)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/",Tag,"_Expression_",genename,"_pos_highlight.png",sep=""),pos_cell_umap1,width=9, height=8)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/",Tag,"_Expression_",genename,"_value.png",sep=""),pos_cell_umap2,width=12, height=6)
  
  #for TFs activity
  if (!(genename %in% colnames(merge_data2)))  next
  
  merge_data3<-merge_data2[order(merge_data2[,genename],decreasing = T),]
  merge_data_target<- merge_data3[,c("UMAP_1","UMAP_2","Treat",genename)]
  colnames(merge_data_target)<-c("UMAP_1","UMAP_2","Treat","Activity")
  
  head(merge_data_target)
  merge_data_target$selected_cluster<-"other"
  merge_data_target[which(merge_data_target$Activity>0),]$selected_cluster<-"Pos"
  merge_data_target$selected_cluster <-factor(merge_data_target$selected_cluster,levels=c("Pos","other"))
  pos_cell_umap1<-ggplot(data=merge_data_target[order(merge_data_target$selected_cluster,decreasing = T),], 
         mapping=aes(x=UMAP_1,y=UMAP_2,colour = selected_cluster))+
    geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
    labs(title =paste0("Activity of ",genename))+scale_color_manual(values=c(ppCor_all[2],"grey"))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                     panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                     axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                     axis.text.y = element_text(size=15),legend.position = "right")
  
  pos_cell_umap2<-ggplot(data=merge_data_target[order(merge_data_target$selected_cluster,decreasing = T),],
                        mapping=aes(x=UMAP_1,y=UMAP_2,color = Activity))+
    geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
    scale_color_gradient(low = "lightgrey",high = "purple")+
    labs(title =paste0("Activity of ",genename))+theme_bw()+
    theme(panel.border = element_rect(colour="grey",fill=NA),
          panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
          axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=15),legend.position = "right")+facet_grid(.~Treat)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/",Tag,"_activity_",genename,"_pos_highlight.pdf",sep=""),pos_cell_umap1,width=9, height=8)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/",Tag,"_activity_",genename,"_value.pdf",sep=""),pos_cell_umap2,width=12, height=6)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/",Tag,"_activity_",genename,"_pos_highlight.png",sep=""),pos_cell_umap1,width=9, height=8)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/",Tag,"_activity_",genename,"_value.png",sep=""),pos_cell_umap2,width=12, height=6)
}

#######################################
#以下未执行
#expression for single gene
gene_target2<-c("HIF1A")
genename<-"HIF1A"
#merge_data3_PBMCs <- merge_data3[which(merge_data3$Tissue == "PBMCs"),]
ggplot(merge_data4[which(merge_data4$gene == genename),], aes(x = Treat, y = Activity))+
  geom_jitter(alpha=0.5,size=0.5, aes(color = Treat ,fill = Treat))+
  geom_boxplot(alpha=0.5,aes(x = Treat, color = Treat,fill = Treat))+
  scale_color_manual(values=ppCor_all2)+  facet_wrap(~Cell_type ,nrow  = 4)+
  xlab("Disease State") +  ylab("Activity") + labs(title = paste0("Activity of ",genename," :: wilcox.test"))+
  theme_bw() +NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                               axis.text.x=element_text(size=15,angle=0,hjust=0.5, vjust=0.5),axis.text.y = element_text(size=15)) +
  theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))+
  stat_compare_means(aes(group = Treat,label=paste0(..p.format..,"\n",..p.signif..)),hide.ns = TRUE, label.x = 1.5,label.y =max(merge_data4[which(merge_data4$gene == genename),]$Activity))
ggplot(express_merge_data4[which(express_merge_data4$gene == genename),], aes(x = Treat, y = Expression))+
  geom_jitter(alpha=0.5,size=0.5, aes(color = Treat ,fill = Treat))+
  geom_boxplot(alpha=0.5,aes(x = Treat, color = Treat,fill = Treat))+
  scale_color_manual(values=ppCor_all2)+  facet_wrap(~Cell_type ,nrow  = 4)+
  xlab("Disease State") +  ylab("expression") + labs(title = paste0("expression of ",genename," :: wilcox.test"))+
  theme_bw() +NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                               axis.text.x=element_text(size=15,angle=0,hjust=0.5, vjust=0.5),axis.text.y = element_text(size=15)) +
  theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))+
  stat_compare_means(aes(group = Treat,label=paste0(..p.format..,"\n",..p.signif..)),hide.ns = TRUE, label.x = 1.5,label.y =max(express_merge_data4[which(express_merge_data4$gene == genename),]$Expression))

library(gghalves)
library(ggpubr)
library(ggsignif)
require(gridExtra)
#https://blog.csdn.net/weixin_43700050/article/details/107512448
#绘制所有细胞中活性值
gene_list<-list()
for ( genename in rev(gene_target)){
  # genename <-"KLF4"
  print(genename)
  merge_data_target<- merge_data2[,c("re_annotation_TFac","Treat",genename)]
  colnames(merge_data_target)<-c("Cell_type","Treat","Activity")
  plot_Cell_act<- ggplot()+  
    # geom_jitter(data=merge_data_target, aes(x=Cell_type,y = Activity),colour="black",alpha=0.5,shape=18,size = 0.5)+
    geom_half_violin(data=merge_data_target %>% filter(Treat =="old"),aes(x=Cell_type,y=Activity,fill=Treat) ,side = "l")+
    geom_half_violin(data=merge_data_target %>% filter(Treat =="young"),aes(x=Cell_type,y=Activity,fill=Treat) ,side = "r")+
    scale_color_manual(values=ppCor_all2)+
    xlab("Cell types") +  ylab("Activity") + labs(title = genename)+
    theme_bw()+NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                                axis.text.x=element_text(size=15,angle=90,hjust=1, vjust=0.5),axis.text.y = element_text(size=15)) +
    theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))
  
  plot_Cell_act2 <-plot_Cell_act+stat_compare_means(method = "wilcox.test",data=merge_data_target,aes(x=Cell_type,y = Activity,group = Treat,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((merge_data_target$Activity)))
  ggsave(paste0("C:/Users/xiaoji/Desktop/差异基因/DEGs_daTFs/barplot_activity_",genename,".pdf",sep=""),plot_Cell_act2,width=30, height=8)
  gene_list<-c(gene_list,list(plot_Cell_act2+ xlab("")+theme(axis.text.x=element_blank())))
}

grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],gene_list[[5]], gene_list[[6]],gene_list[[7]], gene_list[[8]],
             gene_list[[9]], gene_list[[10]],gene_list[[11]], gene_list[[12]],gene_list[[13]], gene_list[[14]],gene_list[[15]], gene_list[[16]],
             gene_list[[17]], gene_list[[18]],gene_list[[19]],  ncol=1)

#for target genes in sig cell types
D_ac_TFs_active_high_matrix <- read.table("D:/0.paper writing/0.paper writing/4.figure and table/中间文件存储/SCENIC/loom_result/ALLcell_all_genes_1215/Def_ac_TFs/merge_ALLcell_all_genes_1215_D_ac_TFs_active_high_matrix.txt", header=T)
DEGs_merge <- read.table("D:/0.paper writing/0.paper writing/4.figure and table/中间文件存储/DEG_Age_default_1215/DEG_005_default_merge_1215.txt", sep = "\t", header = T)

#for active TFs
for ( genename in gene_target[which(!(gene_target %in% c("JUN","IRF1","JUNB")))]){
  #  genename <- "IRF1"
  print(genename)
  cell_target <- as.character(D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$label == genename),]$cell_identity)
  data_target <-merge_data2[which(merge_data2$re_annotation_TFac %in% cell_target),]
  length(unique(as.character(data_target$re_annotation_TFac))) #8
  data_target2<-data_target[,c("re_annotation_TFac","re_annotation_TFac_Treat","Treat",genename)]
  colnames(data_target2)<-c("re_annotation_TFac","re_annotation_TFac_Treat","Treat","Activity")
  
  plot_Cell_pos<-ggplot(data_target2, aes(x = Treat, y = Activity ))+
    geom_jitter(alpha=0.5,size=0.5, aes(color = Treat ,fill = Treat))+
    geom_boxplot(alpha=0.5,aes(x = Treat, color = Treat,fill = Treat))+
    scale_color_manual(values=ppCor_all2)+  facet_wrap(~ re_annotation_TFac,nrow  = 1)+
    xlab("Age group") +  ylab("Activity") + labs(title = paste0("Activity of ",genename," :: wilcox.test"))+
    theme_bw() +NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                                 axis.text.x=element_text(size=15,angle=0,hjust=0.5, vjust=0.5),axis.text.y = element_text(size=15)) +
    theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))+
    stat_compare_means(aes(group = Treat,label=paste0(..p.format..,"\n",..p.signif..)),hide.ns = TRUE, label.x = 1.5,label.y =max(data_target2$Activity))
  ggsave(paste0("C:/Users/xiaoji/Desktop/差异基因/DEGs_daTFs/Activity_in_each_unselected_sig_type_",genename,".pdf",sep=""),plot_Cell_pos,width=length(unique(as.character(data_target$re_annotation_TFac)))*2, height=9)
}

#for active TFs DEGs in the same Cells

target_final<-readRDS(file = "D:/0.paper writing/0.paper writing/4.figure and table/中间文件存储/Cell_1215.rds")
target_expression <- GetAssayData(object = target_final,assay = "SCT", slot = "data")[gene_target,]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))
merge_data_DEGs <- merge(metadata_for_all_cells, target_expression,by=0,sort=FALSE)
merge_data_DEGs$re_annotation_TFac_Treat <-  paste(merge_data_DEGs$re_annotation_TFac, merge_data_DEGs$Treat, sep="_")

for ( genename in gene_target){
  #genename <-"IRF1"
  print(genename)
  cell_target1 <- as.character(D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$label == genename),]$cell_identity)
  cell_target2 <- as.character(DEGs_merge[which(DEGs_merge$gene == genename),]$cluster)
  cell_target <- intersect(cell_target1,cell_target2)
  if(length(cell_target)==0){next}
  data_target_TFs <-merge_data2[which(merge_data2$re_annotation_TFac %in% cell_target),]
  length(unique(as.character(data_target_TFs$re_annotation_TFac))) #8
  data_target_TFs2<-data_target_TFs[,c("re_annotation_TFac","re_annotation_TFac_Treat","Treat",genename)]
  colnames(data_target_TFs2)<-c("re_annotation_TFac","re_annotation_TFac_Treat","Treat","Activity")
  
  data_target_DEGs <-merge_data_DEGs[which(merge_data_DEGs$re_annotation_TFac %in% cell_target),]
  length(unique(as.character(data_target_DEGs$re_annotation_TFac))) #3
  data_target_DEGs2<-data_target_DEGs[,c("re_annotation_TFac","re_annotation_TFac_Treat","Treat",genename)]
  colnames(data_target_DEGs2)<-c("re_annotation_TFac","re_annotation_TFac_Treat","Treat","Expression")
  
  plot_Cell_TFs<-ggplot(data_target_TFs2, aes(x = Treat, y = Activity ))+
    geom_jitter(alpha=0.5,size=0.5, aes(color = Treat ,fill = Treat))+
    geom_boxplot(alpha=0.5,aes(x = Treat, color = Treat,fill = Treat))+
    scale_color_manual(values=ppCor_all2)+  facet_wrap(~ re_annotation_TFac,nrow  = 1)+
    xlab("Age group") +  ylab("Activity") + labs(title = paste0("Activity of ",genename," :: wilcox.test"))+
    theme_bw() +NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                                 axis.text.x=element_text(size=15,angle=0,hjust=0.5, vjust=0.5),axis.text.y = element_text(size=15)) +
    theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))+
    stat_compare_means(aes(group = Treat,label=paste0(..p.format..,"\n",..p.signif..)),hide.ns = TRUE, label.x = 1.5,label.y =max(data_target_TFs2$Activity))
  
  plot_Cell_DEGs<-ggplot(data_target_DEGs2, aes(x = Treat, y = Expression ))+
    geom_jitter(alpha=0.5,size=0.5, aes(color = Treat ,fill = Treat))+
    geom_boxplot(alpha=0.5,aes(x = Treat, color = Treat,fill = Treat))+
    scale_color_manual(values=ppCor_all2)+  facet_wrap(~ re_annotation_TFac,nrow  = 1)+
    xlab("Age group") +  ylab("Expression") + labs(title = paste0("Expression of ",genename," :: wilcox.test"))+
    theme_bw() +NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                                 axis.text.x=element_text(size=15,angle=0,hjust=0.5, vjust=0.5),axis.text.y = element_text(size=15)) +
    theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))+
    stat_compare_means(aes(group = Treat,label=paste0(..p.format..,"\n",..p.signif..)),hide.ns = TRUE, label.x = 1.5,label.y =max(data_target_DEGs2$Expression))
  ggsave(paste0("C:/Users/xiaoji/Desktop/差异基因/DEGs_daTFs/Activity_expression_in_each_selected_sig_type_",genename,".pdf",sep=""),plot_Cell_TFs/plot_Cell_DEGs,width=length(unique(as.character(data_target_DEGs2$re_annotation_TFac)))*2, height=18)
  
  plot_Cell_TFs2 <- ggdotplot(data_target_TFs2, x = "Treat",y = "Activity", combine = TRUE, color = "Treat", palette = ppCor[1:2],
                              fill = "white",binwidth = 0.00045,add = "mean_sd", add.params = list(size = 1,color="grey"),
                              xlab="Treat",ylab="Regulon Active Score")+labs(title=paste0("RAS for ",genename," :: wilcox.test")) + facet_wrap(~ re_annotation_TFac,nrow  = 1)+NoLegend()+
    stat_compare_means(aes(group = Treat,label=paste0(..p.format..,"\n",..p.signif..)),hide.ns = TRUE, label.x = 1.5,label.y =max(data_target_TFs2$Activity))
  
  plot_Cell_DEGs2<-ggdotplot(data_target_DEGs2, x = "Treat",y = "Expression", combine = TRUE, color = "Treat", palette = ppCor[1:2],
                             fill = "white",binwidth = 0.03,add = "mean_sd", add.params = list(size = 1,color="grey"),
                             xlab="Treat",ylab="Expression")+labs(title= paste0("Expression of ",genename," :: wilcox.test"))+ facet_wrap(~ re_annotation_TFac,nrow  = 1)+NoLegend()+
    stat_compare_means(aes(group = Treat,label=paste0(..p.format..,"\n",..p.signif..)),hide.ns = TRUE, label.x = 1.5,label.y =max(data_target_DEGs2$Expression))
  ggsave(paste0("C:/Users/xiaoji/Desktop/差异基因/DEGs_daTFs/dot_Activity_expression_in_each_selected_sig_type_",genename,".pdf",sep=""),plot_Cell_TFs2/plot_Cell_DEGs2,width=length(unique(as.character(data_target_DEGs2$re_annotation_TFac)))*2, height=18)
  
}


