rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(scales)
library(ggsci)
library(ggplot2)
library(RColorBrewer)
library(circlize)

library(dendextend)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)


pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

##read gene expression matrix
target_exprs<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/normalized_gene_expression_CTBs_EVTs_TFs_regulon_activity.txt", sep = '\t', row.names = 1, header = TRUE)
target_exprs[1:4,];dim(target_exprs)#29882    11
table(target_exprs$group)
#CTBs_CTRL  CTBs_USM EVTs_CTRL  EVTs_USM 
##for CTBs
Exprs_CTBs_CTRL<-target_exprs[which(target_exprs$group == "CTBs_CTRL"),]
Exprs_CTBs_USM<-target_exprs[which(target_exprs$group == "CTBs_USM"),]
Exprs_CTBs<-rbind(Exprs_CTBs_CTRL,Exprs_CTBs_USM)
head(Exprs_CTBs)
gene_name<-c("TFAP2A","TEF")
Exprs_CTBs2<-Exprs_CTBs[,c(gene_name,"group")]
Exprs_CTBs2$median_value<-apply(Exprs_CTBs2[,gene_name],1,median)
Exprs_CTBs2<-Exprs_CTBs2[order(Exprs_CTBs2$group,Exprs_CTBs2$median_value,decreasing = F),]
head(Exprs_CTBs2)

table(Exprs_CTBs$group)
#CTBs_CTRL  CTBs_USM 
# 14503      2080 

Exprs_CTBs3<-as.data.frame(t(Exprs_CTBs2[,gene_name]))
Exprs_CTBs3[,1:4]

annotation_col<-data.frame(group = factor(Exprs_CTBs2$group))
rownames(annotation_col) = colnames(Exprs_CTBs3)
ann_colors = list(group=c(CTBs_CTRL=ppCor[2],CTBs_USM=ppCor[1]))

AverageExp_data<-Exprs_CTBs3
range(AverageExp_data)##0.000000 2.947309
#AverageExp_data[AverageExp_data == "Inf"]<- 2000
p1<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = T,show_colnames = F,
             annotation_col = annotation_col,#annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             #gaps_row = ,
             gaps_col =14503,cutree_col = 2,
             treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #scale ="row", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","orange","firebrick3"))(50),
             main ="Normalizated gene expression:no scale by row",angle_col ="90")

p2<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = T,show_colnames = F,
             annotation_col = annotation_col,#annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             #gaps_row = ,
             gaps_col =14503,cutree_col = 2,
             treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             scale ="row", border=FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(50),
             main ="Normalizated gene expression:scale by row",angle_col ="90")

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/noscaled_CTBs_TEF_TFAP2A_heatmap_plot.pdf",width = 7,height=3)
print(p1)
dev.off()
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/scaled_CTBs_TEF_TFAP2A_heatmap_plot.pdf",width = 7,height=3)
print(p2)
dev.off()
###scale by row
head(Exprs_CTBs2)
#get metadata
branch_meta<-annotation_col
head(branch_meta);dim(branch_meta)
#get expression matrix
pt.matrix <-Exprs_CTBs3
pt.matrix[1:4,1:4];dim(pt.matrix)#2 16583

##reorder branch cells
##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
dim(pt.matrix)# 2 16583

#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

##get gaps number for row
groups<-names(table(Exprs_CTBs$group))
gap_number0<-as.numeric(table(Exprs_CTBs$group))
gap_number<-gap_number0[-length(gap_number0)]

#set annotation for each cells  in column
head(branch_meta);dim(branch_meta)
table(branch_meta$group)
#CTBs_CTRL  CTBs_USM 
#14503        2080 

#set annotation for each genes in row
column_colors = list(group=c(CTBs_CTRL=ppCor[2],CTBs_USM=ppCor[1]))
column_ha = HeatmapAnnotation(group=as.character(Exprs_CTBs2$group),col = column_colors)
##plot rearrange heatmap
htkm <- Heatmap(plot_matrix,name= "z-score",border = TRUE,top_annotation = column_ha, 
                column_split = as.character(branch_meta$group),
                cluster_column_slices = FALSE,column_title = "Gene Expression in CTBs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/smooth_CTBs_TEF_TFAP2A_heatmap_plot.pdf",width = 7,height=3)
print(htkm)
dev.off()
##for EVTs
Exprs_EVTs_CTRL<-target_exprs[which(target_exprs$group == "EVTs_CTRL"),]
Exprs_EVTs_USM<-target_exprs[which(target_exprs$group == "EVTs_USM"),]
Exprs_EVTs<-rbind(Exprs_EVTs_CTRL,Exprs_EVTs_USM)
head(Exprs_EVTs)
gene_name<-c( "ARNT2","ATF3","SMAD3","TEAD1","ELF3","GATA2","HSF1","SNAI1","TFAP2A")
Exprs_EVTs2<-Exprs_EVTs[,c(gene_name,"group")]
Exprs_EVTs2$median_value<-apply(Exprs_EVTs2[,gene_name],1,median)
Exprs_EVTs2<-Exprs_EVTs2[order(Exprs_EVTs2$group,Exprs_EVTs2$median_value,decreasing = F),]
head(Exprs_EVTs2)

table(Exprs_EVTs$group)
#EVTs_CTRL  EVTs_USM 
# 7675      5624 

Exprs_EVTs3<-as.data.frame(t(Exprs_EVTs2[,gene_name]))
Exprs_EVTs3[,1:4]

annotation_col<-data.frame(group = factor(Exprs_EVTs2$group))
rownames(annotation_col) = colnames(Exprs_EVTs3)


annotation_row = data.frame(Class = factor(rep(c("UP_TFs","Down_TFs"), c(4,5))))
rownames(annotation_row) = rownames(Exprs_EVTs3)

ann_colors = list(group=c(EVTs_CTRL=ppCor[2],EVTs_USM=ppCor[1]),Class=c(UP_TFs=ppCor[3],Down_TFs=ppCor[4]))

AverageExp_data<-Exprs_EVTs3
range(AverageExp_data)##0.00000 3.82654
#AverageExp_data[AverageExp_data == "Inf"]<- 2000
p1<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = T,show_colnames = F,
             annotation_col = annotation_col,annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             gaps_row = 4,
             gaps_col =7675,cutree_col = 2,
             treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #scale ="row", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","orange","firebrick3"))(50),
             main ="Normalizated gene expression:no scale by row",angle_col ="90")


p2<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = T,show_colnames = F,
             annotation_col = annotation_col,annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             gaps_row = 4,
             gaps_col =7675,cutree_col = 2,
             treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             scale ="row", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(50),
             main ="Normalizated gene expression:scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/noscaled_EVTs_ATF3_TEAD1_GATA2_SNAI1_TFAP2A_HSF1_ELF3_SMAD3_ARNT2_heatmap_plot.pdf",width = 6,height=6)
print(p1)
dev.off()
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/scaled_EVTs_ATF3_TEAD1_GATA2_SNAI1_TFAP2A_HSF1_ELF3_SMAD3_ARNT2_heatmap_plot.pdf",width = 6,height=6)
print(p2)
dev.off()


###scale by row
head(Exprs_EVTs2)

#get metadata
branch_meta<-annotation_col
head(branch_meta);dim(branch_meta)
#get expression matrix
pt.matrix <-Exprs_EVTs3
pt.matrix[1:4,1:4];dim(pt.matrix)#316 18364

##reorder branch cells
##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix[1:4,1:4];dim(pt.matrix)# 1174 18330

#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

##get gaps number for row
groups<-names(table(Exprs_EVTs$group))
gap_number0<-as.numeric(table(Exprs_EVTs$group))
gap_number<-gap_number0[-length(gap_number0)]

#set annotation for each cells  in column
head(branch_meta);dim(branch_meta)
table(branch_meta$group)
#EVTs_CTRL  EVTs_USM 
#7675      5624 

#set annotation for each genes in row
column_colors = list(group=c(EVTs_CTRL=ppCor[2],EVTs_USM=ppCor[1]))
column_ha = HeatmapAnnotation(group=as.character(Exprs_EVTs2$group),col = column_colors)

row_colors = list(Class=c(UP_TFs=ppCor[3],Down_TFs=ppCor[4]))
row_ha = rowAnnotation(Class=rep(c("UP_TFs","Down_TFs"), c(4,5)),col = row_colors)

##plot rearrange heatmap
htkm <- Heatmap(plot_matrix,name= "z-score",border = TRUE,top_annotation = column_ha, 
                #right_annotation =row_mark,
                left_annotation =row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), 
                row_split =rep(c("UP_TFs","Down_TFs"), c(4,5)),
                column_split = as.character(branch_meta$group),
                cluster_column_slices = FALSE,column_title = "Gene Expression in EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/smooth_EVTs_ATF3_TEAD1_GATA2_SNAI1_TFAP2A_HSF1_ELF3_SMAD3_ARNT2_heatmap_plot.pdf",width = 7,height=3)
print(htkm)
dev.off()
