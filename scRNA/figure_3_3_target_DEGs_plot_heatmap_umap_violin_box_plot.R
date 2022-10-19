rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(scales)
library(ggsci)
library(gghalves)
library(ggpubr)
library(ggsignif)
require(gridExtra)
library(VennDiagram)
library(reshape2)
library(Seurat)

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

conflict_prefer("filter", "dplyr")

#selected candiidated group for venn plot
#gene file 1: single cell DEGs
merge_data0<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge.txt",header =T,sep="\t",check.names = F)
head(merge_data0)
sc_Up_gene<-unique(as.character(merge_data0[which(merge_data0$avg_logFC > 0),]$gene))
length(sc_Up_gene)# 1528
sc_Down_gene<-unique(as.character(merge_data0[which(merge_data0$avg_logFC < 0),]$gene))
length(sc_Down_gene)# 1627

#gene file 2: DEGs from bulk
bulk_merge_data<-read.table(file="/mnt/data/chenwei/gongchen/2.map_result/count_file/4.count/Abortion_vs_CTRL.DEG_information_pvalue005_FC1.5.txt",header =T,sep="\t") 
head(bulk_merge_data);dim(bulk_merge_data)#5316    4
bulk_Up_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange > 0),]$ID))
length(bulk_Up_gene)# 2748
bulk_Down_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange < 0),]$ID))
length(bulk_Down_gene)# 2568

#gene file 3: single cell dacTFs
cell_name="ALLcell_all_genes"
D_ac_TFs_active_high_matrix <- read.table( paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Different_TFs/merge_",cell_name,"_D_ac_TFs_active_high_matrix.txt"), header=T)
head(D_ac_TFs_active_high_matrix);dim(D_ac_TFs_active_high_matrix)
#relationship from single cell genes and 
#FOR UP_gene VERSUS Down_gene
venn <-venn.diagram(list(sc_Up_DEGs=sc_Up_gene,sc_Down_DEGs= sc_Down_gene,
                         bulk_Up_DEGs=bulk_Up_gene,bulk_Down_DEGs=bulk_Down_gene),
                    alpha=c(0.5,0.5,0.5,0.5),lwd=1,lty=1, col="white",fill=ppCor[c(1,2,3,6)], 
                    cex = 1.5,cat.col=ppCor[c(1,2,3,6)], cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage(); grid.draw(venn)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/MF1i_Venn_plot_sc_DEGs_bulk_DEGs_for_Abortion.pdf",width =8,height = 6)
grid.draw(venn)
dev.off()

venn2<-venn.diagram(list(All_scDEGs=unique(merge_data0$gene),All_sig_DATF= unique(D_ac_TFs_active_high_matrix$lable)),
                                   alpha=c(0.7,0.7),
                                   lwd=1,lty=1,col="black" , fill=ppCor[c(5:6)], cex = 1.5, cat.col=ppCor[c(5:6)],
                                   cat.fontface=4, cat.cex = 1.5,main = paste0(cell_name,":p<0.01(active_high)"),
                                   main.cex = 1.5, #main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.newpage(); 
grid.draw(venn2)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Venn_plot_sc_DEGs_def_acTF_for_Abortion.pdf",width =6,height = 6)
 grid.draw(venn2)
dev.off()

Up_acTF<-as.character(unique(D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$change == "Up"),]$label))
Down_acTF<-as.character(unique(D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$change == "Down"),]$label))

#for common TFs
length(Up_acTF);length(Down_acTF)#290 210
length(Reduce(intersect,list(Up_acTF,Down_acTF)))#136

venn_acTF<-venn.diagram(list(Up_acTF=Up_acTF,Down_acTF= Down_acTF),
                                   alpha=c(0.7,0.7),lwd=1,lty=1,col="black" ,  fill=ppCor[c(1:2)], 
                                   cex = 1.5, cat.col=ppCor[c(1:2)],cat.fontface=4,  cat.cex = 1.5,    
                                   main = paste0(cell_name,":",":p<0.01(active_high)"),
                                   main.cex = 1.5,# main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.newpage();grid.draw(venn_acTF)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/SF4b_Venn_plot_up_down_def_acTF_for_Abortion.pdf",width =6,height = 6)
grid.draw(venn_acTF)
dev.off()

Reduce(intersect,list(sc_Up_gene,sc_Down_gene,Up_acTF,Down_acTF))
# "TFAP2A" "SOX4"   "ELF3"   "JUND"   "FOS"    "FOSB"   "HMGB3"  "EGR1"   "MAFF"   
Reduce(intersect,list(sc_Up_gene,Up_acTF,Down_acTF))
#[1] "TFAP2A" "SOX4"   "CEBPD"  "ELF3"   "ESRRA"  "JUND"   "KLF6"   "FOS"    "KLF4"   "FOSB"   "TEAD1"  "KLF5"   "TCF7L2" "ASCL2"  "STAT2"  "HMGB3" 
#[17] "GCM1"   "CEBPB"  "EGR1"   "ATF3"   "FOXP1"  "MAFF"   "NR3C1"  "HIF1A"  "REL" 
Reduce(intersect,list(sc_Up_gene,Up_acTF))#48
Reduce(intersect,list(sc_Down_gene,Up_acTF,Down_acTF))
# [1] "TFAP2C"  "BHLHE40" "SOX4"    "MYCN"    "ELF3"    "ELF1"    "TFAP2A"  "NFIL3"   "EGR1"    "FOS"     "FOSB"    "HMGB3"   "NR2F6"   "LMO2"   
#[15] "MAFF"    "JUND"  
Reduce(intersect,list(sc_Down_gene,Down_acTF))
#[1] "ETV5"    "TEAD4"   "NFE2L3"  "TFAP2C"  "NFYB"    "TBPL1"   "HEY1"    "BHLHE40" "SOX4"    "MYCN"    "ELF3"    "ELF1"    "TFAP2A"  "NFIL3"  
#[15] "EGR1"    "FOS"     "FOSB"    "HMGB3"   "NR2F6"   "LMO2"    "MAFF"    "JUND"   
setdiff(Reduce(intersect,list(sc_Up_gene,Up_acTF)),c(sc_Down_gene,Down_acTF))
#[1] "TFDP2"   "NR2F2"   "PIR"     "ZNF83"   "NFE2L2"  "CREM"    "IRF1"    "MAFG"    "MXD4"    "SMAD1"   "SOX18"   "TWIST2"  "PRDM1"   "MEIS1"  
#[15] "NR2F1"   "ZBTB7A"  "EBF1"    "BHLHE41" "BCL11A"  "FOXO3"  
setdiff(Reduce(intersect,list(sc_Down_gene,Down_acTF)),c(sc_Up_gene,Up_acTF))
#[1] "ETV5"   "TEAD4"  "NFE2L3" "NFYB"   "TBPL1"  "HEY1

venn_3<-venn.diagram(list(sc_Up_gene=sc_Up_gene,sc_Down_gene=sc_Down_gene,Up_DATF= Up_acTF,Down_DATF=Down_acTF ),
                                   alpha=c(0.5,0.5,0.5,0.5),
                                   lwd=1,lty=1, #lwd用于设定圆弧的宽度，lty用于设定圆弧的线型
                                   col="black" , #指定图形的圆周边缘颜色  transparent 透明   
                                   fill=ppCor[c(1,2,3,5)], #参数fill表示各个集合对应的圆的填充颜色,
                                   cex = 1.5,    #每个区域label名称的大小
                                   cat.col=ppCor[c(1,2,3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,  #字体格式
                                   cat.cex = 1.5,      #每个分类名称大小
                                   main = paste0(cell_name,"::DEG_005_default_and_DE_acTFq005_p<0.01(active_high)"),
                                   main.cex = 1.5,# main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.newpage();grid.draw(venn_3)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Venn_plot_sc_DEGs_def_acTF_for_Abortion2.pdf",width =8,height = 6)
grid.draw(venn_3)
dev.off()

bulk_sc_Up_DEGs<- setdiff(Reduce(intersect,list(bulk_Up_gene,sc_Up_gene)),c(sc_Down_gene,bulk_Down_gene))
bulk_sc_Down_DEGs<-setdiff(Reduce(intersect,list(sc_Down_gene,bulk_Down_gene)),c(bulk_Up_gene,sc_Up_gene))
setdiff(Reduce(intersect,list(bulk_sc_Up_DEGs,Up_acTF)),c(bulk_sc_Down_DEGs,Down_acTF))
#[1] "NFE2L2" "ZNF83"  "CREM"   "MAFG"   "IRF1"  

setdiff(Reduce(intersect,list(bulk_sc_Down_DEGs,Down_acTF)),c(bulk_sc_Up_DEGs,Up_acTF))
#[1] "TEAD4" "HEY1" 

venn_bulk_sc_DEGs_da_TFs<-venn.diagram(list(bulk_sc_Up_DEGs=bulk_sc_Up_DEGs,bulk_sc_Down_DEGs=bulk_sc_Down_DEGs,Up_DATF= Up_acTF,Down_DATF=Down_acTF ),
                                   alpha=c(0.5,0.5,0.5,0.5),lwd=1,lty=1,col="black" , fill=ppCor[c(1,2,3,5)], 
                                   cex = 1.5, cat.col=ppCor[c(1,2,3,5)],cat.fontface=4,cat.cex = 1.5,
                                   main = paste0(cell_name,"::DEG_005_default_and_DE_acTFq005_p<0.01(active_high)"),
                                   main.cex = 1.5,# main.fontface = 1.5, main.fontfamily =1.5,
                                   filename = NULL)
grid.newpage(); grid.draw(venn_bulk_sc_DEGs_da_TFs)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Venn_plot_bulk_sc_DEGs_def_acTF_for_Abortion.pdf",width =8,height = 6)
grid.draw(venn_bulk_sc_DEGs_da_TFs)
dev.off()
grid.newpage()


#merge DEGs imformation and TFs activity information
merge_data0$trend<-"UP"
merge_data0[which(merge_data0$avg_logFC <0),]$trend<-"Down"
merge_data0$gene_cell_trend<-paste(merge_data0$gene,merge_data0$cluster,merge_data0$trend,sep="_")

D_acTFs_data<-D_ac_TFs_active_high_matrix[,c("lable","cell_identity","change")]
D_acTFs_data$gene_cell_trend<-paste(D_acTFs_data$lable,D_acTFs_data$cell_identity,D_acTFs_data$change,sep="_")
merge(merge_data0,D_acTFs_data)

merge_data0$same_expression_activity<-"NO"
merge_data0[which(merge_data0$gene_cell_trend %in% D_acTFs_data$gene_cell_trend),]$same_expression_activity<-"Yes"
table(merge_data0$same_expression_activity)
# NO  Yes 
#5914   17 
unique(as.character(merge_data0[which(merge_data0$same_expression_activity == "Yes"),]$gene))
# [1] "ETV5"    "TEAD4"   "NFE2L3"  "TFAP2C"  "NFYB"    "TBPL1"   "HEY1"    "BHLHE40" "SOX4"    "MYCN"    "ELF3"    "ELF1"    "TFAP2A"  "NFIL3"  

unique(as.character(merge_data0[which(merge_data0$same_expression_activity == "Yes"),]$cluster))
##"STBs_1" "STBs_2" "STBs_3" "EVTs_1" "EVTs_2" "EVTs_3"
unique(as.character(merge_data0[which(merge_data0$same_expression_activity == "Yes"),]$gene_cell_trend))

#[1] "ETV5_STBs_1_Down"    "TEAD4_STBs_1_Down"   "NFE2L3_STBs_2_Down"  "TFAP2C_STBs_2_Down"  "NFYB_STBs_2_Down"    "TBPL1_STBs_3_Down"  
#[7] "NFE2L3_STBs_3_Down"  "HEY1_EVTs_1_Down"    "HEY1_EVTs_2_Down"    "BHLHE40_EVTs_2_Down" "SOX4_EVTs_2_Down"    "MYCN_EVTs_2_Down"   
#[13] "ELF3_EVTs_2_Down"    "ELF1_EVTs_2_Down"    "TFAP2A_EVTs_2_Down"  "NFIL3_EVTs_2_Down"   "SOX4_EVTs_3_Down" 

head(merge_data0)
write.table(as.data.frame(merge_data0), file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge_TFs_information_add.txt",row.names=T, col.names=T,sep="\t") 

#barplot2
cell_type_all0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","STCs","FBs","Mycs","HCs","NKs","Ts","Bs")
merge_data0<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge_TFs_information_add.txt",row.names=1, header =T,sep="\t") 
head(merge_data0)
merge_data0$trend_activity<-paste(merge_data0$trend,merge_data0$same_expression_activity,sep="_")
merge_data0$Num<-1
DEGs_merge_data2<-aggregate(merge_data0$Num, by=list(merge_data0$cluster,merge_data0$trend,merge_data0$same_expression_activity), FUN=sum)
colnames(DEGs_merge_data2)<-c("cell_identity","tendency","TFs_same","DEGs_num")
#DEGs_merge_data2$cell_identity<-as.character(DEGs_merge_data2$cell_identity)
DEGs_merge_data2$tendency<-factor(DEGs_merge_data2$tendency,levels=c("UP","Down"))
DEGs_merge_data2$TFs_same<-factor(DEGs_merge_data2$TFs_same,levels=c("NO","Yes"))
#DEGs_merge_data2<-DEGs_merge_data2[order(DEGs_merge_data2$TFs_same,decreasing = T),]
DEGs_merge_data2$cell_identity<-factor(DEGs_merge_data2$cell_identity,levels=cell_type_all0)

head(DEGs_merge_data2)
TF_DEG_plot<-ggplot(data=DEGs_merge_data2,mapping=aes(x= cell_identity,y= DEGs_num,group=tendency,fill=TFs_same))+
  geom_bar(stat="identity",width=0.8,position= 'dodge',color="grey")+
  scale_fill_manual(name="TFs same trend",values=ppCor[1:2])+
  geom_text(aes(label=DEGs_num,y=DEGs_num+2),size=3,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="DEGs Number",x="Cell types",title="The number of Abortion related DEGs in each Cells")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 10),axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")
TF_DEG_plot
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/barplot_same_TF_Up_Down_DEGs_number_order_in_all_celltypes.pdf",TF_DEG_plot,width=10, height=8)


#准备表达矩阵和细胞分布注释
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

DefaultAssay(target_final) <- "RNA"
target_final <- NormalizeData(target_final, verbose = FALSE)

#extract metadata information
metadata_Cell<-target_final@meta.data
final_umap<-data.frame(Embeddings(target_final[["umap"]]))
#final_umap$raw_id <- rownames(final_umap)
head(metadata_Cell);head(final_umap)
metadata_Cell_Umap_add <- merge(metadata_Cell,final_umap,by=0)
head(metadata_Cell_Umap_add)
rownames(metadata_Cell_Umap_add)<-metadata_Cell_Umap_add$Row.names

levels(target_final);Idents(target_final)
table(target_final$re_annotation_TFac)

gene_UP_sel<-setdiff(Reduce(intersect,list(bulk_sc_Up_DEGs,Up_acTF)),c(bulk_sc_Down_DEGs,Down_acTF))
gene_Down_sel<-setdiff(Reduce(intersect,list(bulk_sc_Down_DEGs,Down_acTF)),c(bulk_sc_Up_DEGs,Up_acTF))
gene_target<-c(gene_UP_sel,gene_Down_sel)
gene_target
#[1]  "NFE2L2" "ZNF83"  "CREM"   "MAFG"   "IRF1"   "TEAD4"  "HEY1" 
length(gene_target)#7
#umap plot for target genes
for ( gene_target2 in gene_target){
  genename<-gene_target2
  # genename <- "SERPINE1"
  print(genename)
  genename <-FeaturePlot(object = target_final, features = genename,cols= c("grey", "red"),split.by = "Treat",pt.size = 1)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/Umap_split_",gene_target2,".pdf",sep=""), genename,width=12, height=6)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/Umap_split_",gene_target2,".png",sep=""), genename,width=12, height=6)
  }


plots <- VlnPlot(target_final, features = gene_target, split.by = "Treat", group.by = "re_annotation_TFac", split.plot = TRUE,pt.size = 0.01, combine = FALSE)
plots1<-CombinePlots(plots = plots, ncol =1)
plots2 <- VlnPlot(target_final, features = rev(gene_target), split.by = "Treat", group.by = "re_annotation_TFac", split.plot = TRUE,pt.size = 0, combine = FALSE)
plots3<-CombinePlots(plots = plots2, ncol =1)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/sc_bulk_dacTF_gene_expression_violin_plot3.pdf", plots3,width=12, height=3*length(gene_target))
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/sc_bulk_dacTF_gene_expression_violin_plot2.pdf", plots1,width=12, height=3*length(gene_target))
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/sc_bulk_dacTF_gene_expression_violin_plot3.png", plots3,width=12, height=3*length(gene_target))
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/sc_bulk_dacTF_gene_expression_violin_plot2.png", plots1,width=12, height=3*length(gene_target))

#VlnPlot(object = target_final, features = gene_target, group.by = "sub_cluster_brief",cols = ppCor_all2, split.by = "Age_group", split.plot = TRUE,pt.size = 0.01, combine = TRUE,ncol=1)

#数据的提取
target_expression <- GetAssayData(object = target_final,assay = "RNA", slot = "data")[gene_target,]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))
merge_data2 <- merge(metadata_Cell_Umap_add, target_expression,by=0,sort=FALSE)
merge_data2$re_annotation_TFac_Treat <-  paste(merge_data2$re_annotation_TFac, merge_data2$Treat, sep="_")

#https://blog.csdn.net/weixin_43700050/article/details/107512448
#for target genes in all celltypes
gene_list<-list()
for ( genename in rev(gene_target)){
  # genename <-"SERPINE1"
  print(genename)
  merge_data_target<- merge_data2[,c("re_annotation_TFac","Treat",genename)]
  colnames(merge_data_target)<-c("Cell_type","Treat","Expression_level")
  plot_Cell_exp<- ggplot()+  
    # geom_jitter(data=merge_data_target, aes(x=Cell_type,y = Activity),colour="black",alpha=0.5,shape=18,size = 0.5)+
    geom_half_violin(data=merge_data_target[which(merge_data_target$Treat =="Abortion"),],aes(x=Cell_type,y=Expression_level,fill=Treat) ,side = "l")+
    geom_half_violin(data=merge_data_target[which(merge_data_target$Treat =="CTRL"),],aes(x=Cell_type,y=Expression_level,fill=Treat) ,side = "r")+
    scale_color_manual(values=ppCor)+
    xlab("Cell types") +  ylab("Expression_level") + labs(title = genename)+
    theme_bw()+NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                                axis.text.x=element_text(size=15,angle=90,hjust=1, vjust=0.5),axis.text.y = element_text(size=15)) +
    theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))
  
  plot_Cell_exp2 <-plot_Cell_exp+stat_compare_means(method = "wilcox.test",data=merge_data_target,aes(x=Cell_type,y = Expression_level,group = Treat,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((merge_data_target$Expression_level)))
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/Statistic_barplot_expression_in_all_celltypes_",genename,".pdf",sep=""),plot_Cell_exp2,width=30, height=8)
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/Statistic_barplot_expression_in_all_celltypes_",genename,".png",sep=""),plot_Cell_exp2,width=30, height=8)
  gene_list<-c(gene_list,list(plot_Cell_exp2+ xlab("")+theme(axis.text.x=element_blank())))
}
length(gene_target)

gene_list_merge<-grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],gene_list[[5]], gene_list[[6]],ncol=1)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/target_gene_TFs/Statistic_barplot_expression_in_all_celltypes_merge_all_genes.pdf",gene_list_merge,width=10, height=3*length(gene_target))

#for target genes in target celltypes
for ( gene_target2 in gene_target){
  #gene_target2 <-"SMPDL3A"
  genename <-gene_target2
  print(genename)
  #cell_target <-c("HCs","MyCs")
  cell_target <- as.character(merge_data0[which(merge_data0$gene == genename),]$cluster)
  data_target <-merge_data2[which(merge_data2$re_annotation_TFac%in% cell_target),]
  length(unique(as.character(data_target$re_annotation_TFac))) #8
  data_target2<-data_target[,c("re_annotation_TFac","re_annotation_TFac_Treat","Treat",genename)]
  colnames(data_target2)<-c("re_annotation_TFac","re_annotation_TFac_Treat","Treat","expression")
  plot_Cell_pos<-ggplot(data_target2, aes(x = Treat, y = expression ))+
    geom_jitter(alpha=0.5,size=0.5, aes(color = Treat ,fill = Treat))+
    geom_boxplot(alpha=0.5,aes(x = Treat, color = Treat,fill = Treat))+
    scale_color_manual(values=ppCor)+  facet_wrap(~ re_annotation_TFac,nrow  = 1)+
    xlab("Treat group") +  ylab("log adjust count") + labs(title = paste0("Expression level of ",genename," :: wilcox.test"))+
    theme_bw() +NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                                 axis.text.x=element_text(size=15,angle=0,hjust=0.5, vjust=0.5),axis.text.y = element_text(size=15)) +
    theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))+
    stat_compare_means(aes(group = Treat,label=paste0(..p.format..,"\n",..p.signif..)),hide.ns = TRUE, label.x = 1.5,label.y =max(data_target2$expression))
  ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/Statistic_barplot_expression_in_target_celltypes_",genename,".pdf",sep=""),plot_Cell_pos,width=2*length(cell_target), height=8)

  }

#heatmap for  target genes in target celltypes
DefaultAssay(target_final)<-"RNA"
target_final@assays$RNA@scale.data <- scale(target_final@assays$RNA@data, scale = TRUE)

#for target genes in target celltypes
for ( gene_target2 in gene_target){
  #gene_target2 <-"NFE2L2"
  genename <-gene_target2
  print(genename)
  #cell_target <-c("HCs","MyCs")
 cell_target <- as.character(merge_data0[which(merge_data0$gene == genename),]$cluster)
 data_target <-merge_data2[which(merge_data2$re_annotation_TFac %in% cell_target),]
 data_target$re_annotation_TFac<-factor(data_target$re_annotation_TFac,levels = cell_target)
 data_target <-data_target[order(data_target$re_annotation_TFac,data_target$Treat),]
 cell_treat<- as.character(unique(data_target$re_annotation_TFac_Treat))
 
 target_plot<-subset(x = target_final,subset = re_annotation_TFac %in% cell_target)

 target_plot$re_annotation_TFac_Treat <-  paste(target_plot$re_annotation_TFac, target_plot$Treat, sep="_")
 target_plot$re_annotation_TFac_Treat<-factor(target_plot$re_annotation_TFac_Treat,levels = cell_treat)
 Idents(object = target_plot) <- 're_annotation_TFac_Treat'

 PH<-DoHeatmap(target_plot,features =gene_target2, draw.lines = TRUE,label = F,group.colors=ppCor_all[1:length(levels(target_plot$re_annotation_TFac_Treat))])
 PH_all<-PH+scale_fill_gradientn(colors = c("lightblue", "orange", "red"))+theme(axis.text.y = element_text(size = 10))

 #selected 1K cell
 target_1K<-subset(target_plot, downsample = 1000)
 scaled_data_1K <- target_1K@assays$SCT@scale.data
 #heatmap for all cells
 PH<-DoHeatmap(object = target_1K, features =gene_target2,draw.lines = TRUE,label = F,disp.min = -2.5, disp.max = 2.5,assay = "RNA",size = 3, angle = -50, hjust=0.8,
              group.by ="re_annotation_TFac_Treat",group.bar.height = 0.1,group.colors=ppCor[1:length(levels(target_1K$re_annotation_TFac_Treat))]) 
 PH_1K<-PH+scale_fill_gradientn(colors = c("blue", "orange", "red"))+theme(axis.text.y = element_text(size = 10))
 ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/heatmap_expression_in_target_celltypes_",genename,".pdf",sep=""),PH_all,width=8, height=6)
 ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/target_gene/target_gene_TFs/heatmap_expression_in_target_celltypes_1K_",genename,".pdf",sep=""),PH_1K,width=8, height=6)
}

