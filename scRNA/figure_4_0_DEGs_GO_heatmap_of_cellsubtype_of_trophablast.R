rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(Seurat)
library(future)

library(data.table)
library(scales)
library(dplyr)
library(psych)

library(ggplot2)
library(cowplot)
library(ggsci)
library(grid)
library(gridExtra)

library(pheatmap)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
grid.newpage()

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
#n=5;barplot(rep(1,n), col=colorRampPalette(colors = c('red', 'white'))( n ))
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
show_col(ppCor_all2)
##plot cells
Cells_col<-colorRampPalette(colors = rev(Cells_col_raw))(80)
length(unique(Cells_col))
barplot(rep(1,80), col=Cells_col)
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(80)
barplot(rep(1,80), col=Cells_col)
#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
plan(strategy = "multicore", workers = 42)
options(future.globals.maxSize = 1500 * 1024^12)


#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

cell_name<-"Trophoblast_KRT7"
selected_cell<-subset(x = target_final, subset = re_anno_TFac_major == cell_name)
subgroup_order<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3")
selected_cell$re_annotation_TFac<- factor(selected_cell$re_annotation_TFac,levels = rev(subgroup_order))
#selected_cell$re_annotation_TFac<- factor(selected_cell$re_annotation_TFac,levels = subgroup_order)

Allcluster.markers_pos_001<-read.table(file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_default0.25_pos_001_genes.xls",sep=""),header =T,sep="\t")
head(Allcluster.markers_pos_001)
table(Allcluster.markers_pos_001$cluster)
Allcluster.markers_pos_001$cluster<-factor(Allcluster.markers_pos_001$cluster,levels = rev(subgroup_order))

#identify the relationship among different cell types
Allcluster.markers_pos_001<-Allcluster.markers_pos_001[order(Allcluster.markers_pos_001$cluster,decreasing = T),]
plot_maker<-Allcluster.markers_pos_001 %>% group_by(cluster)
top50 <- plot_maker %>% top_n(n = 50, wt = avg_logFC)
top10 <- plot_maker %>% top_n(n = 10, wt = avg_logFC)
top2 <- plot_maker %>% top_n(n = 2, wt = avg_logFC)
top1 <- plot_maker %>% top_n(n = 1, wt = avg_logFC)


#GO BP enrichment
allcluster_BP_GO=data.frame();allcluster_BP_GO_whole=data.frame()
for (celltype in subgroup_order) {
  #celltype<-"CTBs_1"
  print(paste0("celltype : ",celltype))
  genenames=as.character(plot_maker[which(plot_maker$cluster==celltype),]$gene)
  gene_change=bitr(genenames, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  BP_GO <- enrichGO(gene= unique(gene_change$ENTREZID),OrgDb= org.Hs.eg.db,keyType= 'ENTREZID',ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff= 0.1,readable= TRUE)
  BP_GO2 <- simplify(BP_GO, cutoff=0.7, by="p.adjust", select_fun=min)
  go_res=BP_GO2@result
  go_res_whole=BP_GO@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=celltype;allcluster_BP_GO=rbind(allcluster_BP_GO,go_res)}
  if (dim(go_res_whole)[1] != 0) {
    go_res_whole$cluster=celltype;allcluster_BP_GO_whole=rbind(allcluster_BP_GO_whole,go_res_whole)}
}

head(allcluster_BP_GO[,c("ID","Description","qvalue","cluster")])
table(allcluster_BP_GO$cluster)
write.table(allcluster_BP_GO, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_simple_GO_BP_for_cellsubtype_of_Trophoblast_KRT7.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(allcluster_BP_GO_whole, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_whole_GO_BP_for_cellsubtype_of_Trophoblast_KRT7.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)


##plot cluster unique DEGs heatmap
DefaultAssay(selected_cell) <- "RNA"
selected_cell <- NormalizeData(selected_cell, verbose = FALSE)
Idents(object = selected_cell) <- 're_annotation_TFac'

head(selected_cell@meta.data)
table(selected_cell$re_annotation_TFac)
#EVTs_3 EVTs_2 EVTs_1 STBs_3 STBs_2 STBs_1 CTBs_2 CTBs_1 
# 4746   6705   1848    903   1955    372   8746   7837 
names(table(selected_cell$re_annotation_TFac))

##3. heatmap for averageExpression
#maker_gene<-unique(as.character(top2$gene))
maker_gene<-unique(as.character(plot_maker[which(plot_maker$avg_logFC>=log(2)),]$gene))
AverageExp<-AverageExpression(selected_cell,assays="RNA",slot = "data",features=maker_gene)

typeof(AverageExp);head(AverageExp$RNA)
range(AverageExp$RNA)#0.001988972 9.633643413
AverageExp_data<-log2(AverageExp$RNA+1)
range(AverageExp_data);head(AverageExp_data)
#AverageExp_data[AverageExp_data == "Inf"]<- 2000
colnames(AverageExp_data)
AverageExp_data<-AverageExp_data[,subgroup_order]
head(AverageExp_data)
p1<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
             show_rownames = F,show_colnames = T,border=TRUE,border_color ="grey", 
             scale ="row", color = colorRampPalette(c("navy","white","firebrick3"))(100),
             main ="log2(mean of expression level of DEGs with 2 FC + 1):scale by row")
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/heatmap_for_average_expression_of_submaker_FC2_for_cellsubtype_of_Trophoblast_KRT7.pdf",width =8,height = 12)
print(p1)
dev.off()

#for all DEGs
maker_gene<-unique(as.character(plot_maker$gene))
AverageExp<-AverageExpression(selected_cell,assays="RNA",slot = "data",features=maker_gene)
typeof(AverageExp);head(AverageExp$RNA)
range(AverageExp$RNA)#0.02774331 793.35683163
AverageExp_data<-log2(AverageExp$RNA+1)
range(AverageExp_data);head(AverageExp_data)
#AverageExp_data[AverageExp_data == "Inf"]<- 2000
colnames(AverageExp_data)
AverageExp_data<-AverageExp_data[,subgroup_order]

p1<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
             show_rownames = F,show_colnames = T,border=TRUE,border_color ="grey", 
             scale ="row", color = colorRampPalette(c("navy","white","firebrick3"))(60),
             main ="log2(mean of expression level of all DEGs + 1):scale by row")
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/heatmap_for_average_expression_of_submaker_ALL_for_cellsubtype_of_Trophoblast_KRT7.pdf",width =8,height = 12)
print(p1)
dev.off()

##类群特异性DEGs精选

#selected intercluster DEGs unique in each cell types and perform corresponding GO BP enrichment
allcluster_BP_GO=data.frame();allcluster_BP_GO_whole=data.frame();gene_uniq_list<-c()
for (celltype in subgroup_order) {
  #celltype<-"CTBs_1"
  print(paste0("celltype : ",celltype))
  genenames=as.character(plot_maker[which(plot_maker$cluster==celltype),]$gene)
  other_genes<-as.character(plot_maker[which(plot_maker$cluster !=celltype),]$gene)
  genename_uniq<-setdiff(genenames,other_genes)
  gene_uniq_list<-c(gene_uniq_list,genename_uniq)
  gene_change=bitr(genename_uniq, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  BP_GO <- enrichGO(gene= unique(gene_change$ENTREZID),OrgDb= org.Hs.eg.db,keyType= 'ENTREZID',ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff= 0.1,readable= TRUE)
  BP_GO2 <- simplify(BP_GO, cutoff=0.7, by="p.adjust", select_fun=min)
  go_res=BP_GO2@result
  go_res_whole=BP_GO@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=celltype;allcluster_BP_GO=rbind(allcluster_BP_GO,go_res)}
  if (dim(go_res_whole)[1] != 0) {
    go_res_whole$cluster=celltype;allcluster_BP_GO_whole=rbind(allcluster_BP_GO_whole,go_res_whole)}
}

head(allcluster_BP_GO[,c("ID","Description","qvalue","cluster")])
table(allcluster_BP_GO$cluster)
#CTBs_1 CTBs_2 EVTs_1 EVTs_2 EVTs_3 STBs_1 STBs_2 STBs_3 
# 116      8     20      6     72     60     10     18 
table(allcluster_BP_GO_whole$cluster)
#CTBs_1 CTBs_2 EVTs_1 EVTs_2 EVTs_3 STBs_1 STBs_2 STBs_3 
# 2522   1202   2079   1939   3780   4155   1884   2473 

write.table(allcluster_BP_GO, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Unique_gene_Enrichment_simple_GO_BP_for_cellsubtype_of_Trophoblast_KRT7.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(allcluster_BP_GO_whole, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Unique_gene_Enrichment_whole_GO_BP_for_cellsubtype_of_Trophoblast_KRT7.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

length(gene_uniq_list)
#1949
plot_maker2<-plot_maker[which(plot_maker$gene %in% gene_uniq_list),]
write.table(as.data.frame(plot_maker2), file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Unique_biomaker_identified_for_cellsubtype_of_Trophoblast_KRT7.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

DefaultAssay(selected_cell) <- "RNA"
Idents(object = selected_cell) <- 're_annotation_TFac'
##3. heatmap for averageExpression
#maker_gene<-unique(as.character(top2$gene))
maker_gene<-unique(gene_uniq_list)
AverageExp<-AverageExpression(selected_cell,assays="RNA",slot = "data",features=maker_gene)

typeof(AverageExp);head(AverageExp$RNA)
range(AverageExp$RNA)#0.0000 155.3957
AverageExp_data<-log2(AverageExp$RNA+1)
range(AverageExp_data);head(AverageExp_data)#0.000000 7.289057
colnames(AverageExp_data)
AverageExp_data<-AverageExp_data[,subgroup_order]
head(AverageExp_data)

table(plot_maker2$cluster)
#EVTs_3 EVTs_2 EVTs_1 STBs_3 STBs_2 STBs_1 CTBs_2 CTBs_1 
#340    132     92    172    125    568    119    401 
gap_number<-c(401,520,1088,1213,1385,1477,1609)
p1<-pheatmap(AverageExp_data, cluster_row =FALSE,cluster_col =FALSE,
             gaps_row =gap_number,
             show_rownames = F,show_colnames = T,border=TRUE,border_color ="grey", 
             scale ="row", color = colorRampPalette(c("navy","white","firebrick3"))(100),
             main ="log2(mean of expression level of DEGs+ 1):scale by row")
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/heatmap_for_average_expression_of_unique_submaker_for_cellsubtype_of_Trophoblast_KRT7.pdf",width =8,height = 12)
print(p1)
dev.off()




###待完善
#similarity of all cell types
maker_gene<-unique(as.character(plot_maker$gene))
AverageExp<-AverageExpression(selected_cell,assays="RNA",slot = "data",features=maker_gene)
typeof(AverageExp);
head(AverageExp$RNA)
range(AverageExp$RNA)#0.000 1528.887
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
p2<-pheatmap(coorda$r)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/correlation_heatmap_for_average_expression_of_submaker_ALL_for_cellsubtype_of_Trophoblast_KRT7.pdf",width =8,height = 8)
print(p2)
dev.off()

##4.整体细胞类群表达水平绘制
target_1K<-subset(selected_cell, downsample = 1000)
target_1K@assays$RNA@scale.data <- scale(target_1K@assays$RNA@data, scale = TRUE)
target_1K$re_annotation_TFac<-factor(target_1K$re_annotation_TFac,levels =subgroup_order)
#target_1K <- ScaleData(object = target_1K, features = rownames(target_1K))
PH<-DoHeatmap(object = target_1K, features =maker_gene,draw.lines = TRUE,label = F,
              disp.min = -2.5, disp.max = 2.5,assay = "RNA",slot = "scale.data", size = 3, angle = -50, hjust=0.8,
              group.by ="final_major_subgroup_brief",group.bar.height = 0.1,group.colors=ppCor[1:length(levels(target_1K$final_major_subgroup_brief))]) 
PH2<-PH+scale_fill_gradientn(colors = c("lightblue", "white", "red"))+theme(axis.text.y = element_text(size = 10))
PH2
ggsave(PH2,file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/heatmap_for_single_expression_of_submaker_ALL_for_cellsubtype_of_1K_Trophoblast_KRT7.png",width = 15,height = 10)
ggsave(PH,file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/heatmap_for_single_expression_of_submaker_ALL_for_cellsubtype_of_1K_Trophoblast_KRT7_2.png",width = 15,height = 10)

#手动提取
levels(target_final)
DefaultAssay(target_final)<-"SCT"
#mat0 <- GetAssayData(target, slot = "counts")
#mat0 <- log2(mat0 + 1)
#其他SCT的数据说明：https://www.jianshu.com/p/e639cc257d51 
mat0<- GetAssayData(target_final, slot = "data")

#FC_10<-Allcluster.markers[Allcluster.markers$p_val_adj<0.01 & Allcluster.markers$avg_logFC>1,]%>% group_by(cluster)
#FC_1.7<-Allcluster.markers[Allcluster.markers$p_val_adj<0.01 & Allcluster.markers$avg_logFC>0.25,]%>% group_by(cluster)
#top10 <- FC_1.7 %>% top_n(n = 10, wt = avg_logFC)

#get genes and cluster information
#gene_features <-top10$gene
FC_1.5<-plot_maker[plot_maker$p_val_adj<0.01 & plot_maker$avg_logFC>log(1.5),]%>% group_by(cluster)
top50_FC_1.5 <- FC_1.5 %>% top_n(n = 50, wt = avg_logFC)
table(top50_FC_1.5$cluster)
names(table(top50_FC_1.5$cluster))
gene_features <- unique(top50_FC_1.5$gene);length(gene_features)
cluster_info <- target_final$final_major_subgroup_brief

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
length(unique(gene_features))#2105
#name_mean1$row
head(name_mean1)
gene_features_uniq<-name_mean1$Gname
gene_features_dup<-gene_features[which(duplicated(gene_features))]
length(gene_features);length(gene_features_uniq);length(gene_features_dup)
#重新去unique数据框
mat2 <- mat1[gene_features_uniq,]
#手动scale
mat3<-scale(t(mat2), center=T,scale=T)
range(mat3)#-4.798357 162.184118
mat4<-t(mat3)
mat4[mat4>=2.5]= 2.5
mat4[mat4 < (-2.5)]= -2.5 #小于负数时，加括号！
#set color for cell types
col <- ppCor_all[1:length(levels(cluster_info))]
#colourCount<-length(unique(target$sub_cluster))
names(col) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col),labels = levels(cluster_info), labels_gp = gpar(cex = 0.5, col = "black")))

#important genes notes
mark_gene0 <- c("VIM","HLA-B","CDH1","EGFR","HLA-G","MMP2","CSH1","CGA")
#mark_gene0<-as.character(top2$gene)
gene_pos <- which(rownames(mat4) %in% mark_gene0)
mark_gene1<-rownames(mat4)[gene_pos]
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = mark_gene1))

#legend adjust
col_fun1  <- circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red"))
col_fun2  <- circlize::colorRamp2(c(-2,0,2), c("purple", "black", "yellow"))
col_fun<-col_fun1
#show_col(col_fun(seq(-3, 3)))

plot_HP<-Heatmap(mat4,cluster_rows = FALSE,cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 #show_row_names = TRUE,
                 column_split = cluster_info,
                 top_annotation = top_anno,
                 #  right_annotation = row_anno,
                 column_title = NULL,
                 col = col_fun,
                 heatmap_legend_param = list(legend_height = unit(4, "cm"),grid_width = unit(0.4, "cm"),title = "log2(count+1)",title_position = "leftcenter-rot",at = c(-2, 0, 2), labels = c("low", "median", "high") ))
pdf("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/submaker_top50_FC_1.5_heatmap1.pdf",width = 15,height = 8)
draw(plot_HP, row_title = "Genes", 
     #     padding = unit(c(2, 0, 2, 0), "mm"),
     adjust_annotation_extension = TRUE,  # default
     column_title = paste0("major_submaker_heatmap_SCT_data"))
dev.off()
