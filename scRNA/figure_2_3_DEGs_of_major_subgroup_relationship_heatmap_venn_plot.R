rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(reshape2)
library(gridExtra)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(ggsci)
library(psych)
library(pheatmap)
library(VennDiagram)
library(grid)
library(futile.logger)
##sample information plot in generation
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
show_col(pal)
pal
#[1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF"
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
show_col(nejm)
nejm
#[1] "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"
#decide color
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

#n=5;barplot(rep(1,n), col=colorRampPalette(colors = c('red', 'white'))( n ))
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
show_col(Cells_col_raw)
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
show_col(ppCor_all)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(80)
show_col(Cells_col)


#对差基因的不同组关系进行比较
#For default genes 
##for majar subtypes
#read DEGs between age list
files<- list.files("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/default_DEGs/Pos_genes_01q")
#list.files命令将input文件夹下所有文件名输入
data_list<-data_frame()
for ( f in files[1:length(files)]){
  # f=files[1]
  FF<-paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/default_DEGs/Pos_genes_01q/",f,sep="")
  column <- read.delim(FF,header = T,sep = "\t")[,c("gene","cluster","avg_logFC")]
  if (is.na(column[1,]$cluster)) {
    next
  }
  data_list<- rbind(data_list,column)
}
str(data_list)

files<- list.files("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/default_DEGs/Neg_genes_01q")
data_list2<-data_frame()
for ( f in files[1:length(files)]){
  #f<-files[1]
  FF<-paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/default_DEGs/Neg_genes_01q/",f,sep="")
  column <- read.delim(FF,header = T,sep = "\t")[,c("gene","cluster","avg_logFC")]
  if (is.na(column[1,]$cluster)) {
    next
  }
  data_list2<- rbind(data_list2,column)
}
str(data_list2)



cell_type<-unique(c(as.character(data_list$cluster),as.character(data_list2$cluster)))
cell_type
cell_type_all1<-c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2",
                  "STCs","FBs","Epi","Endo","MyCs","HCs","NKs","Ts","Bs","Ery")

head(data_list)
data_list$cluster<- factor(data_list$cluster,levels=cell_type_all1,ordered=TRUE)
data_list1<-data_list[order(data_list$cluster),]
data_list1$gene<-as.character(data_list1$gene)
dim(data_list1);str(data_list1);length(unique(data_list1$gene)) #up genes: 1876

data_list2$cluster<- factor(data_list2$cluster,levels=cell_type_all1,ordered=TRUE)
data_list3<-data_list2[order(data_list2$cluster),]
data_list3$gene<-as.character(data_list3$gene)
dim(data_list3);str(data_list3);length(unique(data_list3$gene)) #down genes:1821
data_list1[which(is.na(data_list1$gene)),];data_list3[which(is.na(data_list3$gene)),]
head(data_list1);head(data_list3)

merge_data<-as.data.frame(rbind(data_list1,data_list3))
gene_order<-unique(merge_data$gene)
#merge_data$gene<-as.factor(as.character(merge_data$gene))
dim(merge_data);str(merge_data);head(merge_data);length(unique(merge_data$gene))#3259
merge_data[which(is.na(merge_data$gene)),]

#remove Erythrocyte_mix and Megakaryocytes
merge_data1<-merge_data[which(!(merge_data$cluster %in% c("Ery"))),]
cell_type_all1
merge_data1$cluster<- factor(merge_data1$cluster,levels=cell_type_all1[-length(cell_type_all1)],ordered=TRUE)

merge_data2<- dcast(merge_data1[,c("cluster","gene","avg_logFC")],gene~cluster)
merge_data2$gene<- factor(merge_data2$gene,levels=gene_order,ordered=TRUE)
merge_data2<-merge_data2[order(merge_data2$gene),]
head(merge_data2);dim(merge_data2)#2859 
rownames(merge_data2)<-merge_data2$gene
merge_data2<-merge_data2[,-(1)]
merge_data2[is.na(merge_data2)] <- 0
merge_data2[merge_data2>0] <- 1
merge_data2[merge_data2<0] <- c(-1)
head(merge_data2);tail(merge_data2)
merge_data3<-t(merge_data2)
dim(merge_data3) ##18 2859
range(merge_data3)# -1  1
merge_data3[1:6,1:6]

#Building anno 
cell_type_all2<-c(rep("Cytotrophoblasts_PAGE4",4),rep("Extravilous_trophoblast_HLA_G",3),rep("Syncytiotrophoblast_CGA", 2),
                  rep("Stromal_cells_DCN",2),rep("Epithelial_Cell_EPCAM",1),rep("Endothelial_Cell_PECAM1",1),
                  rep("Myeloid_Cell_AIF1",2),rep("Leukomonocyte",3),rep("Erythrocyte_HBA1",1))
cell_type_all3<-c(rep("Trophoblasts_KRT7",9),rep("Stromal_cells_DCN",2),rep("Epithelial_Cell_EPCAM",1),rep("Endothelial_Cell_PECAM1",1),
                  rep("Myeloid_Cell_AIF1",2),rep("Leukomonocyte",3),rep("Erythrocyte_HBA1",1))
annotation_cell_type <-data.frame(cluster=cell_type_all1[-length(cell_type_all1)],main_group=cell_type_all2[-length(cell_type_all2)],main_group2=cell_type_all3[-length(cell_type_all3)])
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
         main = "Disease all DEG_defualt(q=0.1)",
         legend_breaks = c(-0.6,0,0.6),legend_labels = c("Down","Unchange","Up"),
         color = colorRampPalette(colors = c("blue","gray","red"))(3))

pheatmap(merge_data3,cluster_rows=F, cluster_cols =T, 
         annotation_row=annotation_cell_type[,c("main_group","main_group2")],labels_col = labels_col,
         annotation_colors = anno_colors,
         main = "Disease all DEG_defualt(q=0.1)",
         legend_breaks = c(-0.6,0,0.6),legend_labels = c("Down","Unchange","Up"),
         color = colorRampPalette(colors = c("blue","gray","red"))(3))

pheatmap(merge_data3,cluster_rows=T, cluster_cols =T, 
         annotation_row=annotation_cell_type[,c("main_group","main_group2")],labels_col = labels_col,
         annotation_colors = anno_colors,
         main = "Disease all DEG_defualt(q=0.1)",
         legend_breaks = c(-0.6,0,0.6),legend_labels = c("Down","Unchange","Up"),
         color = colorRampPalette(colors = c("blue","gray","red"))(3))

#存储相关关系文件
#merge_data0<-as.data.frame(merge_data)
#annotation_cell_type1<-annotation_cell_type
#annotation_cell_type1$cluster<-rownames(annotation_cell_type1)
#merge_data_save<-merge(merge_data0,annotation_cell_type1)
#dim(merge_data_save);dim(merge_data0)
#head(merge_data_save)
head(merge_data1)
table(merge_data1$cluster)
write.table(merge_data1, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/DEG_01q_default_merge.txt",quote=F, row.names=F, col.names=T,sep="\t") 
head(merge_data2)
write.table(merge_data2, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/relationship_for_Cellsubtypes_DEG_01q_default_merge.txt",quote=F, row.names=F, col.names=T,sep="\t") 

#selected candiidated group for venn plot
merge_data0<-read.table(file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/DEG_01q_default_merge.txt",header =T,sep="\t") 
head(merge_data0)
merge_data_rm0<-merge_data0[which(!(merge_data0$cluster %in% c("Ery"))),]

Vagina_up<-unique(as.character(merge_data_rm0[which(merge_data_rm0$avg_logFC > 0),]$gene))
length(Vagina_up)# 130
Vagina_down<-unique(as.character(merge_data_rm0[which(merge_data_rm0$avg_logFC < 0),]$gene))
length(Vagina_down)# 147

Up_gene<-unique(c(Vagina_up));Down_gene<-unique(c(Vagina_down))
#FOR UP_gene VERSUS Down_gene
venn <-venn.diagram(list(Up_gene=Up_gene,Down_gene= Down_gene),
                     alpha=c(0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","navy"), 
                     cex = 1.5,cat.col=c("red","navy"), cat.fontface=4,cat.cex = 1.5,    
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                     rotation.degree=90,filename = NULL)

grid.newpage(); 
grid.draw(venn)
grid.newpage(); 
