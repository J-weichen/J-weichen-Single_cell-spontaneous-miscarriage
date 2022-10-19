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
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)
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
#selected candiidated group for venn plot
#gene file 1: single cell DEGs
merge_data0<-read.table(file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/DEG_01q_default_merge.txt",header =T,sep="\t") 
head(merge_data0)
merge_data_rm0<-merge_data0[which(!(merge_data0$cluster %in% c("Ery"))),]
sc_Up_gene<-unique(as.character(merge_data_rm0[which(merge_data_rm0$avg_logFC > 0),]$gene))
length(sc_Up_gene)# 1551
sc_Down_gene<-unique(as.character(merge_data_rm0[which(merge_data_rm0$avg_logFC < 0),]$gene))
length(sc_Down_gene)# 1601

#gene file 2: DEGs from bulk
bulk_merge_data<-read.table(file="/mnt/data/chenwei/gongchen/2.map_result/count_file/4.count/Abortion_vs_CTRL.DEG_information_pvalue005_FC1.5.txt",header =T,sep="\t") 
head(bulk_merge_data);dim(bulk_merge_data)#5316    4
bulk_Up_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange > 0),]$ID))
length(bulk_Up_gene)# 2748
bulk_Down_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange < 0),]$ID))
length(bulk_Down_gene)# 2568

#gene file 3: single cell dacTFs
cell_name="ALLcell_all_genes"
D_ac_TFs_active_high_matrix <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/merge_",cell_name,"_D_ac_TFs_active_high_matrix.txt"), header=T)
head(D_ac_TFs_active_high_matrix)

#relationship from single cell genes and 
#FOR UP_gene VERSUS Down_gene
venn <-venn.diagram(list(sc_Up_DEGs=sc_Up_gene,sc_Down_DEGs= sc_Down_gene,
                         bulk_Up_DEGs=bulk_Up_gene,bulk_Down_DEGs=bulk_Down_gene),
                    alpha=c(0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=ppCor[c(1,2,3,6)], 
                    cex = 1.5,cat.col=ppCor[c(1,2,3,6)], cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)

grid.newpage(); 
grid.draw(venn)
grid.newpage(); 


venn_Age_DEGs_genAge<-venn.diagram(list(All_DEGs=unique(merge_data_rm0$gene),All_sig_DATF= unique(D_ac_TFs_active_high_matrix$lable)),
                                   alpha=c(0.7,0.7),
                                   lwd=1,lty=1,col="black" , fill=ppCor[c(5:6)], cex = 1.5, cat.col=ppCor[c(5:6)],
                                   cat.fontface=4, cat.cex = 1.5,main = paste0(cell_name,":p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.newpage(); 
grid.draw(venn_Age_DEGs_genAge)
grid.newpage(); #清空画板，开始画新图

Up_acTF<-as.character(unique(D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$change == "Up"),]$label))
Down_acTF<-as.character(unique(D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$change == "Down"),]$label))

#for common TFs
length(Up_acTF);length(Down_acTF)#235 166
length(Reduce(intersect,list(Up_acTF,Down_acTF)))#87

grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(Up_acTF=Up_acTF,Down_acTF= Down_acTF),
                                   alpha=c(0.7,0.7),lwd=1,lty=1,col="black" ,  fill=ppCor[c(1:2)], 
                                   cex = 1.5, cat.col=ppCor[c(1:2)],cat.fontface=4,  cat.cex = 1.5,    
                                   main = paste0(cell_name,":",":p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)

Reduce(intersect,list(sc_Up_gene,sc_Down_gene,Up_acTF,Down_acTF))
#"FOS"    "JUND"   "HMGB3"  "FOSB"   "SOX4"   "TFDP2"  "GATA2"  "TFAP2A"
Reduce(intersect,list(sc_Up_gene,Up_acTF,Down_acTF))
#[1] "FOS"    "JUND"   "HMGB3"  "FOSB"   "CEBPB"  "ATF3"   "KLF6"   "SOX4"   "ASCL2"  "GATA3"  "STAT2" 
#[12] "TFDP2"  "FOXP1"  "REL"    "GATA2"  "SNAI1"  "TFAP2A" "ELF3"   "ESRRA"  "NR3C1"  "HIF1A"  "ELF1" 
Reduce(intersect,list(sc_Up_gene,Up_acTF))
Reduce(intersect,list(sc_Down_gene,Up_acTF,Down_acTF))
# [1] "MYCN"    "TFDP2"   "BHLHE40" "TFAP2A"  "UGP2"    "FOS"     "HMGB3"   "GATA2"   "NR2F6"  
#[10] "SOX4"    "JUND"    "BCLAF1"  "FOSB"   
Reduce(intersect,list(sc_Down_gene,Down_acTF))

setdiff(Reduce(intersect,list(sc_Up_gene,Up_acTF)),c(sc_Down_gene,Down_acTF))
setdiff(Reduce(intersect,list(sc_Down_gene,Down_acTF)),c(sc_Up_gene,Up_acTF))

grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(sc_Up_gene=sc_Up_gene,sc_Down_gene=sc_Down_gene,Up_DATF= Up_acTF,Down_DATF=Down_acTF ),
                                   alpha=c(0.9,0.9,0.9,0.9),
                                   lwd=1,lty=1, #lwd用于设定圆弧的宽度，lty用于设定圆弧的线型
                                   col="black" , #指定图形的圆周边缘颜色  transparent 透明   
                                   fill=ppCor[c(1,2,3,5)], #参数fill表示各个集合对应的圆的填充颜色,
                                   cex = 1.5,    #每个区域label名称的大小
                                   cat.col=ppCor[c(1,2,3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,  #字体格式
                                   cat.cex = 1.5,      #每个分类名称大小
                                   main = paste0(cell_name,"::DEG_005_default_and_DE_acTFq005_p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)
grid.newpage(); #清空画板，开始画新图


bulk_sc_Up_DEGs<- setdiff(Reduce(intersect,list(bulk_Up_gene,sc_Up_gene)),c(sc_Down_gene,bulk_Down_gene))
bulk_sc_Down_DEGs<-setdiff(Reduce(intersect,list(sc_Down_gene,bulk_Down_gene)),c(bulk_Up_gene,sc_Up_gene))
setdiff(Reduce(intersect,list(bulk_sc_Up_DEGs,Up_acTF)),c(bulk_sc_Down_DEGs,Down_acTF))
#[1] "XBP1"  "MAFF"  "SPI1"  "MAFK"  "CREM"  "CEBPA" "MAF"   "IRF1" 
setdiff(Reduce(intersect,list(bulk_sc_Down_DEGs,Down_acTF)),c(bulk_sc_Up_DEGs,Up_acTF))
#[1] "TEAD4" "HEY1"  "EZH2" 

venn_Age_DEGs_genAge<-venn.diagram(list(bulk_sc_Up_DEGs=bulk_sc_Up_DEGs,bulk_sc_Down_DEGs=bulk_sc_Down_DEGs,Up_DATF= Up_acTF,Down_DATF=Down_acTF ),
                                   alpha=c(0.9,0.9,0.9,0.9),lwd=1,lty=1,col="black" , fill=ppCor[c(1,2,3,5)], 
                                   cex = 1.5, cat.col=ppCor[c(1,2,3,5)],cat.fontface=4,cat.cex = 1.5,
                                   main = paste0(cell_name,"::DEG_005_default_and_DE_acTFq005_p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5,filename = NULL)
grid.newpage(); 
grid.draw(venn_Age_DEGs_genAge)
grid.newpage(); 


#准备表达矩阵和细胞分布注释
target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
target_final<-subset(x = target,subset = final_major_subgroup== "Troph_mix_group",invert = TRUE)
major_order<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN",
               "Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order<-c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2","MyCs","HCs","STCs","FBs","Endo","Epi", "NKs","Ts","Bs","Ery","Masts")
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief, levels=subgroup_order,ordered=TRUE)
target_final$final_major_group_brief<- factor(target_final$final_major_group_brief, levels=major_order,ordered=TRUE)

#extract metadata information
metadata_Cell<-target_final@meta.data
final_umap<-data.frame(Embeddings(target_final[["umap"]]))
final_umap$raw_id <- rownames(final_umap)
head(metadata_Cell);head(final_umap)
metadata_Cell_Umap_add <- merge(metadata_Cell,final_umap)
head(metadata_Cell_Umap_add)
rownames(metadata_Cell_Umap_add)<-metadata_Cell_Umap_add$raw_id

levels(target_final);Idents(target_final)
table(target_final$final_major_subgroup_brief)

gene_UP_sel<-setdiff(Reduce(intersect,list(bulk_sc_Up_DEGs,Up_acTF)),c(bulk_sc_Down_DEGs,Down_acTF))
gene_Down_sel<-setdiff(Reduce(intersect,list(bulk_sc_Down_DEGs,Down_acTF)),c(bulk_sc_Up_DEGs,Up_acTF))
gene_target<-c(gene_UP_sel,gene_Down_sel)
#gene_target<-c("MMP10","LEFTY2")
length(gene_target)
#umap plot for target genes
for ( gene_target2 in gene_target){
  genename<-gene_target2
  # genename <-"SPI1"
  print(genename)
  genename <-FeaturePlot(object = target_final, features = genename,cols= c("grey", "red"),split.by = "Treat",pt.size = 1)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.gene_expression/Umap_split_",gene_target2,".pdf",sep=""), genename,width=12, height=6)
}


plots <- VlnPlot(target_final, features = gene_target, split.by = "Treat", group.by = "final_major_subgroup_brief", split.plot = TRUE,pt.size = 0.01, combine = FALSE)
plots1<-CombinePlots(plots = plots, ncol =1)
plots2 <- VlnPlot(target_final, features = rev(gene_target), split.by = "Treat", group.by = "final_major_subgroup_brief", split.plot = TRUE,pt.size = 0, combine = FALSE)
plots3<-CombinePlots(plots = plots2, ncol =1)
ggsave("/mnt/data/chenwei/gongchen/3.gene_expression/sc_bulk_dacTF_violin_plot3.pdf", plots3,width=12, height=1.5*length(gene_target))
ggsave("/mnt/data/chenwei/gongchen/3.gene_expression/sc_bulk_dacTF_violin_plot2.pdf", plots1,width=12, height=1.5*length(gene_target))

#VlnPlot(object = target_final, features = gene_target, group.by = "sub_cluster_brief",cols = ppCor_all2, split.by = "Age_group", split.plot = TRUE,pt.size = 0.01, combine = TRUE,ncol=1)

#数据的提取
target_expression <- GetAssayData(object = target_final,assay = "SCT", slot = "data")[gene_target,]
target_expression <- as.data.frame(t(as.data.frame(target_expression)))
merge_data2 <- merge(metadata_Cell_Umap_add, target_expression,by=0,sort=FALSE)
merge_data2$final_major_subgroup_brief_Treat <-  paste(merge_data2$final_major_subgroup_brief, merge_data2$Treat, sep="_")

#https://blog.csdn.net/weixin_43700050/article/details/107512448
#for target genes in all celltypes
gene_list<-list()
for ( genename in rev(gene_target)){
  # genename <-"LEFTY2"  
  print(genename)
  merge_data_target<- merge_data2[,c("final_major_subgroup_brief","Treat",genename)]
  colnames(merge_data_target)<-c("Cell_type","Treat","Expression_level")
  plot_Cell_exp<- ggplot()+  
    # geom_jitter(data=merge_data_target, aes(x=Cell_type,y = Activity),colour="black",alpha=0.5,shape=18,size = 0.5)+
    geom_half_violin(data=merge_data_target %>% filter(Treat =="Abortion"),aes(x=Cell_type,y=Expression_level,fill=Treat) ,side = "l")+
    geom_half_violin(data=merge_data_target %>% filter(Treat =="CTRL"),aes(x=Cell_type,y=Expression_level,fill=Treat) ,side = "r")+
    scale_color_manual(values=ppCor)+
    xlab("Cell types") +  ylab("Expression_level") + labs(title = genename)+
    theme_bw()+NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                                axis.text.x=element_text(size=15,angle=90,hjust=1, vjust=0.5),axis.text.y = element_text(size=15)) +
    theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))
  
  plot_Cell_exp2 <-plot_Cell_exp+stat_compare_means(method = "wilcox.test",data=merge_data_target,aes(x=Cell_type,y = Expression_level,group = Treat,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((merge_data_target$Expression_level)))
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.gene_expression/Statistics/barplot_expression_in_all_celltypes_",genename,".pdf",sep=""),plot_Cell_exp2,width=30, height=8)
  gene_list<-c(gene_list,list(plot_Cell_exp2+ xlab("")+theme(axis.text.x=element_blank())))
}
length(gene_target)

gene_list_merge<-grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],gene_list[[5]], gene_list[[6]],gene_list[[7]], gene_list[[8]],
             gene_list[[9]], gene_list[[10]],gene_list[[11]],ncol=1)
ggsave("/mnt/data/chenwei/gongchen/3.gene_expression/Statistics/barplot_expression_in_all_celltypes_merge_all_genes.pdf",gene_list_merge,width=10, height=3*length(gene_target))

#for target genes in target celltypes
for ( gene_target2 in gene_target){
  #gene_target2 <-"SPI1"
  genename <-gene_target2
  print(genename)
  #cell_target <-c("HCs","MyCs")
  cell_target <- as.character(merge_data_rm0[which(merge_data_rm0$gene == genename),]$cluster)
  data_target <-merge_data2[which(merge_data2$final_major_subgroup_brief%in% cell_target),]
  length(unique(as.character(data_target$final_major_subgroup_brief))) #8
  data_target2<-data_target[,c("final_major_subgroup_brief","final_major_subgroup_brief_Treat","Treat",genename)]
  colnames(data_target2)<-c("final_major_subgroup_brief","final_major_subgroup_brief_Treat","Treat","expression")
  plot_Cell_pos<-ggplot(data_target2, aes(x = Treat, y = expression ))+
    geom_jitter(alpha=0.5,size=0.5, aes(color = Treat ,fill = Treat))+
    geom_boxplot(alpha=0.5,aes(x = Treat, color = Treat,fill = Treat))+
    scale_color_manual(values=ppCor)+  facet_wrap(~ final_major_subgroup_brief,nrow  = 1)+
    xlab("Treat group") +  ylab("log adjust count") + labs(title = paste0("Expression level of ",genename," :: wilcox.test"))+
    theme_bw() +NoLegend()+theme(plot.title = element_text(hjust=0.5,size=16,vjust=0.5,colour = "black",face = "bold"),
                                 axis.text.x=element_text(size=15,angle=0,hjust=0.5, vjust=0.5),axis.text.y = element_text(size=15)) +
    theme(axis.title.x = element_text(size=15,colour = "black",face = "bold"), axis.title.y = element_text(size=15,colour = "black",face = "bold"))+
    stat_compare_means(aes(group = Treat,label=paste0(..p.format..,"\n",..p.signif..)),hide.ns = TRUE, label.x = 1.5,label.y =max(data_target2$expression))
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.gene_expression/Statistics/barplot_expression_in_target_celltypes_",genename,".pdf",sep=""),plot_Cell_pos,width=2*length(cell_target), height=8)
}

#heatmap for  target genes in target celltypes
DefaultAssay(target_final)<-"SCT"
#for target genes in target celltypes
for ( gene_target2 in gene_target){
  #gene_target2 <-"SPI1"
  genename <-gene_target2
  print(genename)
  #cell_target <-c("HCs","MyCs")
 cell_target <- as.character(merge_data_rm0[which(merge_data_rm0$gene == genename),]$cluster)
 data_target <-merge_data2[which(merge_data2$final_major_subgroup_brief %in% cell_target),]
 data_target$final_major_subgroup_brief<-factor(data_target$final_major_subgroup_brief,levels = cell_target)
 data_target <-data_target[order(data_target$final_major_subgroup_brief,data_target$Treat),]
 cell_treat<- as.character(unique(data_target$final_major_subgroup_brief_Treat))
 
 target_plot<-subset(x = target_final,subset = final_major_subgroup_brief %in% cell_target)

 target_plot$final_major_subgroup_brief_Treat <-  paste(target_plot$final_major_subgroup_brief, target_plot$Treat, sep="_")
 target_plot$final_major_subgroup_brief_Treat<-factor(target_plot$final_major_subgroup_brief_Treat,levels = cell_treat)
 Idents(object = target_plot) <- 'final_major_subgroup_brief_Treat'

 PH<-DoHeatmap(target_plot,features =gene_target2, draw.lines = TRUE,label = F,group.colors=ppCor_all[1:length(levels(target_plot$final_major_subgroup_brief_Treat))])
 PH_all<-PH+scale_fill_gradientn(colors = c("lightblue", "orange", "red"))+theme(axis.text.y = element_text(size = 10))

 #selected 1K cell
 target_1K<-subset(target_plot, downsample = 1000)
 scaled_data_1K <- target_1K@assays$SCT@scale.data
 #heatmap for all cells
 PH<-DoHeatmap(object = target_1K, features =gene_target2,draw.lines = TRUE,label = F,disp.min = -2.5, disp.max = 2.5,assay = "SCT",size = 3, angle = -50, hjust=0.8,
              group.by ="final_major_subgroup_brief_Treat",group.bar.height = 0.1,group.colors=ppCor[1:length(levels(target_1K$final_major_subgroup_brief_Treat))]) 
 PH_1K<-PH+scale_fill_gradientn(colors = c("blue", "orange", "red"))+theme(axis.text.y = element_text(size = 10))
 ggsave(paste0("/mnt/data/chenwei/gongchen/3.gene_expression/Statistics/heatmap_expression_in_target_celltypes_",genename,".pdf",sep=""),PH_all,width=8, height=6)
 ggsave(paste0("/mnt/data/chenwei/gongchen/3.gene_expression/Statistics/heatmap_expression_in_target_celltypes_1K_",genename,".pdf",sep=""),PH_1K,width=8, height=6)
}
