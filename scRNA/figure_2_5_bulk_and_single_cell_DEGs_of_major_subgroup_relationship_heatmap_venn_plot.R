rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(scales)
library(ggsci)
library(VennDiagram)
library(reshape2)

##sample information plot in generation
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

#selected candiidated group for venn plot
merge_data0<-read.table(file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/DEG_01q_default_merge.txt",header =T,sep="\t") 
head(merge_data0)
merge_data_rm0<-merge_data0[which(!(merge_data0$cluster %in% c("Ery"))),]
sc_Up_gene<-unique(as.character(merge_data_rm0[which(merge_data_rm0$avg_logFC > 0),]$gene))
length(sc_Up_gene)# 1551
sc_Down_gene<-unique(as.character(merge_data_rm0[which(merge_data_rm0$avg_logFC < 0),]$gene))
length(sc_Down_gene)# 1601
#reading DEGs from bulk
bulk_merge_data<-read.table(file="/mnt/data/chenwei/gongchen/2.map_result/count_file/4.count/Abortion_vs_CTRL.DEG_information_pvalue005_FC1.5.txt",header =T,sep="\t") 
head(bulk_merge_data);dim(bulk_merge_data)#5316    4
bulk_Up_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange > 0),]$ID))
length(bulk_Up_gene)# 2748
bulk_Down_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange < 0),]$ID))
length(bulk_Down_gene)# 2568

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

#selected unique trend DEGs
Up_DEGs<-setdiff(Reduce(intersect,list(sc_Up_gene,bulk_Up_gene)),sc_Down_gene)
Down_DEGs<-setdiff(Reduce(intersect,list(sc_Down_gene,bulk_Down_gene)),sc_Up_gene)
length(Up_DEGs);length(Down_DEGs)
head(merge_data_rm0)

head(bulk_merge_data)
bulk_merge_data_out<-bulk_merge_data[which(bulk_merge_data$ID %in% c(Up_DEGs,Down_DEGs)),]
write.table(bulk_merge_data_out, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/bulk_DEGs_overlapped_DEG_01q_default_merge.txt",quote=F, row.names=F, col.names=T,sep="\t") 

merge_data1<-merge_data_rm0[which(merge_data_rm0$gene %in% c(Up_DEGs,Down_DEGs)),]
merge_data1$trend<-"Up"
merge_data1[which(merge_data1$gene %in% Down_DEGs),]$trend<-"Down"
merge_data1$trend<- factor(merge_data1$trend,levels=c("Up","Down"),ordered=TRUE)

table(merge_data1$cluster)
cell_type_all1<-c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2",
                  "STCs","FBs","Epi","Endo","MyCs","HCs","NKs","Ts","Bs")
merge_data1$cluster<- factor(merge_data1$cluster,levels=cell_type_all1,ordered=TRUE)
merge_data1<-merge_data1[order(merge_data1$trend,merge_data1$cluster,decreasing = F),]
gene_order<-unique(as.character(merge_data1$gene))

merge_data2<- dcast(merge_data1[,c("cluster","gene","avg_logFC")],gene~cluster)
merge_data2$gene<- factor(merge_data2$gene,levels=gene_order,ordered=TRUE)
merge_data2<-merge_data2[order(merge_data2$gene),]
head(merge_data2);dim(merge_data2)#2859   19
rownames(merge_data2)<-merge_data2$gene
merge_data2<-merge_data2[,-(1)]
merge_data2[is.na(merge_data2)] <- 0
merge_data2[merge_data2>0] <- 1
merge_data2[merge_data2<0] <- c(-1)
head(merge_data2);tail(merge_data2)
merge_data3<-t(merge_data2)
dim(merge_data3) ##18 713
range(merge_data3)# -1  1
merge_data3[1:6,1:6]

#Building anno 
cell_type_all2<-c(rep("Cytotrophoblasts_PAGE4",4),rep("Extravilous_trophoblast_HLA_G",3),rep("Syncytiotrophoblast_CGA", 2),
                  rep("Stromal_cells_DCN",2),rep("Epithelial_Cell_EPCAM",1),rep("Endothelial_Cell_PECAM1",1),
                  rep("Myeloid_Cell_AIF1",2),rep("Leukomonocyte",3))
cell_type_all3<-c(rep("Trophoblasts_KRT7",9),rep("Stromal_cells_DCN",2),rep("Epithelial_Cell_EPCAM",1),rep("Endothelial_Cell_PECAM1",1),
                  rep("Myeloid_Cell_AIF1",2),rep("Leukomonocyte",3))
annotation_cell_type <-data.frame(cluster=cell_type_all1,main_group=cell_type_all2,main_group2=cell_type_all3)
rownames(annotation_cell_type) = annotation_cell_type$cluster
head(annotation_cell_type)

anno_colors = list(
  main_group =c(Cytotrophoblasts_PAGE4=Cells_col[1],Extravilous_trophoblast_HLA_G=Cells_col[18],
                Syncytiotrophoblast_CGA=Cells_col[9],Stromal_cells_DCN=Cells_col[14],
                Epithelial_Cell_EPCAM=Cells_col[30],Endothelial_Cell_PECAM1=Cells_col[27],
                Myeloid_Cell_AIF1=Cells_col[22], Leukomonocyte=Cells_col[54]),
  main_group2 =c(Trophoblasts_KRT7=Cells_col[1],Stromal_cells_DCN=Cells_col[14],
                 Epithelial_Cell_EPCAM=Cells_col[30],Endothelial_Cell_PECAM1=Cells_col[27],
                 Myeloid_Cell_AIF1=Cells_col[22], Leukomonocyte=Cells_col[54]))
labels_col = c("")

pheatmap(merge_data3,cluster_rows=F, cluster_cols =F, 
         annotation_row=annotation_cell_type[,c("main_group","main_group2")],
         labels_col = labels_col,
         annotation_colors = anno_colors,
         #gaps_row = c(18),
         main = "Disease all single cell and bulk DEG_defualt(q=0.1)",
         legend_breaks = c(-0.6,0,0.6),legend_labels = c("Down","Unchange","Up"),
         color = colorRampPalette(colors = c("blue","gray","red"))(3))

pheatmap(merge_data3,cluster_rows=F, cluster_cols =T, 
         annotation_row=annotation_cell_type[,c("main_group","main_group2")],labels_col = labels_col,
         annotation_colors = anno_colors,
         main = "Disease all single cell and bulk DEG_defualt(q=0.1)",
         legend_breaks = c(-0.6,0,0.6),legend_labels = c("Down","Unchange","Up"),
         color = colorRampPalette(colors = c("blue","gray","red"))(3))

pheatmap(merge_data3,cluster_rows=T, cluster_cols =T, 
         annotation_row=annotation_cell_type[,c("main_group","main_group2")],labels_col = labels_col,
         annotation_colors = anno_colors,
         main = "Disease all single cell and bulk DEG_defualt(q=0.1)",
         legend_breaks = c(-0.6,0,0.6),legend_labels = c("Down","Unchange","Up"),
         color = colorRampPalette(colors = c("blue","gray","red"))(3))

#存储相关关系文件
head(merge_data1)
table(merge_data1$cluster)
write.table(merge_data1, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/bulk_DEGs_and_DEG_01q_default_merge.txt",quote=F, row.names=F, col.names=T,sep="\t") 
head(merge_data2)
write.table(merge_data2, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/relationship_for_Cellsubtypes_bulk_DEGs_and_DEG_01q_default_merge.txt",quote=F, row.names=F, col.names=T,sep="\t") 
