rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
dyn.load('/mnt/data/chenwei/software/share_road/libudunits2/libudunits2.so.0')
dyn.load('/mnt/data/chenwei/software/share_road/libproj/libproj.so.0')
dyn.load('/mnt/data/chenwei/software/share_road/libgdal.so.20')
dyn.load('/mnt/data/chenwei/software/share_road/libgeos-3.4.2.so')
dyn.load('/mnt/data/chenwei/software/share_road/libgeos_c.so.1')

library(Seurat)
library(monocle3)
library(scales)
library(ggsci)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dendextend)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)

pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)


##reading all 
cds_final<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_fielt_final_monocle3_object.rds")

###analysis for EVTs branch
##heatmap plot for persudotime related genes
cds_sub_CTBs2EVTs <- choose_graph_segments(cds_final)
cds_branch<-cds_final[,rownames(colData(cds_final)) %in% rownames(colData(cds_sub_CTBs2EVTs))]

CTBs2EVTs_branch_plot<-plot_cells(cds_branch,  reduction_method="UMAP",color_cells_by = "pseudotime", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2EVTs_branch_umap.pdf",CTBs2EVTs_branch_plot,width = 8, height =7)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2EVTs_branch_umap.png",CTBs2EVTs_branch_plot,width = 8, height =7)

plot_cells(cds_branch, color_cells_by = "re_annotation_TFac", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
saveRDS(cds_branch,file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2EVTs_selected_branch_monocle3_object.rds")

###analysis for CTBs2EVTs branch
cds_branch<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2EVTs_selected_branch_monocle3_object.rds")
pseu_time<-data.frame(pseudotime(cds_branch))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)
#pseu_time$number<-ifelse(rownames(pseu_time) %in% rownames(colData(cds_branch)),1,-1)
pseu_time$cell_name<-rownames(pseu_time)
pseu_time<-pseu_time[order(pseu_time$pseudotime,decreasing = F),]
head(pseu_time);tail(pseu_time)
write.table(pseu_time, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_pseudotime_info_for_cells_in_CTBs2EVTs_branch.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##重新寻找拟时轨迹相关基因
Track_genes <- graph_test(cds_branch, neighbor_graph="principal_graph", cores=30)
head(Track_genes)
write.table(Track_genes, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_CTBs2EVTs_branch_gene_cor_along_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
#基因表达趋势图
plot_genes_in_pseudotime(cds_branch[Track_genes_sig,],color_cells_by = "re_annotation_TFac", min_expr=0.5, ncol = 2)
#FeaturePlot图
plot_cells(cds_branch, genes=Track_genes_sig, show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)

##寻找共表达模块
###no selection
genelist <- pull(Track_genes, gene_short_name) %>% as.character()##不做筛选
gene_module_unselect <- find_gene_modules(cds_branch[genelist,], resolution=1e-2, cores = 30)
##selection for Track_genes
pr_deg_ids <- row.names(subset(Track_genes, q_value < 0.05 & morans_I > 0.25))#进一步筛选
length(unique(pr_deg_ids));#length(unique(genelist))### 15904 ## 24099
gene_module_df <- find_gene_modules(cds_branch[pr_deg_ids,], resolution=1e-2, cores = 30)## 将轨迹可变基因收集到模块中：
head(gene_module_df);dim(gene_module_df) #1174    5
table(gene_module_df$module)
# 1   2   3   4   5   6   7   8   9 
#283 245 190 135  93  86  70  48  24 
write.table(gene_module_df, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_CTBs2EVTs_gene_module_along_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

##moldule expression in different cell types
gene_module <- gene_module_df
cell_group <- tibble::tibble(cell=row.names(colData(cds_branch)), cell_group=colData(cds_branch)$re_annotation_TFac)
agg_mat <- aggregate_gene_expression(cds_branch, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
head(agg_mat)
agg_mat2<-agg_mat[,c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3")]
pheatmap::pheatmap(agg_mat2, scale="column",clustering_method="ward.D2")
module_heatmap<-pheatmap::pheatmap(agg_mat2, scale="column", cluster_cols=F,clustering_method="ward.D2")

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/heatmap_for_average_expression_of_psudotime_module_ALL_in_CTBs2EVTs.pdf",width =8,height = 12)
print(module_heatmap)
dev.off()
##根据聚类可以再归为5大类
##We can also pass gene_module_df to plot_cells() as we did when we compared clusters in the L2 data above.
all_module_expression<-plot_cells(cds_branch, genes=gene_module_df ,label_cell_groups=FALSE,show_trajectory_graph=FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_CTBs2EVTs_module_expression.png",all_module_expression,width = 20, height =20,limitsize = FALSE)
all_module_expression2<-all_module_expression + scale_color_gradient2(low = "blue", mid = "gray", high =  "red")
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_CTBs2EVTs_module_expression2.png",all_module_expression2,width = 20, height =20,limitsize = FALSE)


##绘制变化动态
gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_CTBs2EVTs_gene_module_along_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
head(gene_module)
gene_module$id<-as.character(gene_module$id)
gene_module<-gene_module[order(gene_module$module,decreasing = F),]
gap_number0<-as.numeric(table(gene_module$module))
gap_number<-gap_number0[-length(gap_number0)]
head(gene_module)

pt.matrix <- exprs(cds_branch)[match(gene_module$id,rownames(rowData(cds_branch))),]
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 18330
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 18330

##target data for yaxi
Yaxi_gene_matrix1<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix1[,1:4];dim(Yaxi_gene_matrix1)#1174 18330

##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix_yaxi<-pt.matrix
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

##target genes for yaxi
Yaxi_gene_matrix2<-pt.matrix_yaxi[which(rownames(pt.matrix_yaxi) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix2[,1:4];dim(Yaxi_gene_matrix2)#1174 18330

Yaxi_gene_matrix3<-pt.matrix[which(rownames(pt.matrix) %in% c("MKI67","SERPINE1","TIMP3")),]
Yaxi_gene_matrix3[,1:4];dim(Yaxi_gene_matrix3)#1174 18330

write.table(Yaxi_gene_matrix1, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/raw_selected_genes_expression_for_Yaxi1.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix2, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/smooth_spline_selected_genes_expression_for_Yaxi2.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(Yaxi_gene_matrix3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/smooth_spline_scale_selected_genes_expression_for_Yaxi3.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)


#rownames(pt.matrix) <- pr_DEGs_branch
pt.matrix[1:4,1:4];dim(pt.matrix)#1296 8232
#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

htkm <- Heatmap(plot_matrix,name= "z-score",#border = TRUE,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(paste0("Module_",1:length(gap_number0)),gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of CTBs2EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_direction_CTBs2EVTs_all_module.pdf",height=15,width=8)
print(htkm)
dev.off()

gene_module$module<-factor(gene_module$module,levels = c(8,7,1,6,5,2,3,4,9))
gene_module<-gene_module[order(gene_module$module,decreasing = F),]
gene_module$type<-"Later"
gene_module[which(gene_module$module %in% c(8)),]$type<-"Early"
gene_module[which(gene_module$module %in% c(7,1,6)),]$type<-"Early_Mid"
gene_module[which(gene_module$module %in% c(5,2)),]$type<-"Mid"
gene_module[which(gene_module$module %in% c(3)),]$type<-"Mid_Later"
gene_module[which(gene_module$module %in% c(4,9)),]$type<-"N_Later"

gene_module$type<-factor(gene_module$type,levels = c("Early","Early_Mid","Mid","Mid_Later","N_Later"))
#gene_module<-gene_module[order(gene_module$type,decreasing = F),]
gene_module$module<-factor(gene_module$module,levels = c(8,7,1,6,5,2,3,4,9))
gene_module<-gene_module[order(gene_module$module,decreasing = F),]
write.table(gene_module, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/developped_genes_rearranged_along_CTBs2EVTs_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)


#GO BP enrichment
allcluster_BP_GO=data.frame();allcluster_BP_GO_whole=data.frame()

for (cluster_name in levels(gene_module$type)) {
  #cluster_name<-"Mid_low"
  print(paste0("cluster number: ",cluster_name))
  small_gene_group=unique(as.character(gene_module[which(gene_module$type == cluster_name),]$id))
  gene_change=bitr(small_gene_group, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  BP_GO <- enrichGO(gene= unique(gene_change$ENTREZID),OrgDb= org.Hs.eg.db,keyType= 'ENTREZID',ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff= 0.1,readable= TRUE)
  BP_GO2 <- simplify(BP_GO, cutoff=0.7, by="p.adjust", select_fun=min)
  go_res=BP_GO2@result
  go_res_whole=BP_GO@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=cluster_name;allcluster_BP_GO=rbind(allcluster_BP_GO,go_res)}
  if (dim(go_res_whole)[1] != 0) {
    go_res_whole$cluster=cluster_name;allcluster_BP_GO_whole=rbind(allcluster_BP_GO_whole,go_res_whole)}
}
head(allcluster_BP_GO[,c("ID","Description","qvalue","cluster")])
table(allcluster_BP_GO$cluster)
# Early Early_Mid     Later       Mid Mid_Later 
#   9        97       141        51       184 
write.table(allcluster_BP_GO, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_for_gene_module_along_CTBs2EVTs_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(allcluster_BP_GO_whole, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_whole_for_module_along_CTBs2EVTs_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

###plot modules heatmaps for rearrange
gap_number0<-as.numeric(table(gene_module$type))
gap_number<-gap_number0[-length(gap_number0)]
head(gene_module)
tail(gene_module)

plot_matrix2<-plot_matrix[as.character(gene_module$id),]
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,
                # row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(c("Early_time","Mid_low","Later_time","Mid_high"),gap_number0),
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(c("Early","Early_Mid","Mid","Mid_Later","Later"),gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of CTBs2EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_direction_CTBs2EVTs_rearrange_module.pdf",height=15,width=8)
print(htkm)
dev.off()

plot_matrix2[1:4,1:4]

####for CTBs2STBs
##reading all 
cds_final<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_fielt_final_monocle3_object.rds")
###analysis for STBs branch
##heatmap plot for persudotime related genes
cds_sub_CTBs2STBs <- choose_graph_segments(cds_final)
cds_branch<-cds_final[,rownames(colData(cds_final)) %in% rownames(colData(cds_sub_CTBs2STBs))]

CTBs2STBs_branch_plot<-plot_cells(cds_branch,  reduction_method="UMAP",color_cells_by = "pseudotime", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2STBs_branch_umap.pdf",CTBs2STBs_branch_plot,width = 8, height =7)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2STBs_branch_umap.png",CTBs2STBs_branch_plot,width = 8, height =7)

plot_cells(cds_branch, color_cells_by = "re_annotation_TFac", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
saveRDS(cds_branch,file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2STBs_selected_branch_monocle3_object.rds")

###analysis for CTBs2STBs branch
cds_branch<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2STBs_selected_branch_monocle3_object.rds")
pseu_time<-data.frame(pseudotime(cds_branch))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)
#pseu_time$number<-ifelse(rownames(pseu_time) %in% rownames(colData(cds_branch)),1,-1)
pseu_time$cell_name<-rownames(pseu_time)
pseu_time<-pseu_time[order(pseu_time$pseudotime,decreasing = F),]
head(pseu_time);tail(pseu_time)

##重新寻找拟时轨迹相关基因
Track_genes <- graph_test(cds_branch, neighbor_graph="principal_graph", cores=30)
head(Track_genes)
write.table(Track_genes, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_CTBs2STBs_branch_gene_cor_along_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

#CTBs_1 CTBs_2 STBs_1 STBs_2 STBs_3 EVTs_1 EVTs_2 EVTs_3 
# 6706   8304    308   1904    739    220    101     82 
#Track_genes<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_CTBs2STBs_branch_gene_cor_along_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)

#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
#基因表达趋势图
plot_genes_in_pseudotime(cds_branch[Track_genes_sig,],color_cells_by ="re_annotation_TFac", min_expr=0.5, ncol = 2)
#FeaturePlot图
plot_cells(cds_branch, genes=Track_genes_sig, show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)

##寻找共表达模块
###no selection
genelist <- pull(Track_genes, gene_short_name) %>% as.character()##不做筛选
gene_module_unselect <- find_gene_modules(cds_branch[genelist,], resolution=1e-2, cores = 30)
##selection for Track_genes
pr_deg_ids <- row.names(subset(Track_genes, q_value < 0.05 & morans_I > 0.25))#进一步筛选
length(unique(pr_deg_ids));#length(unique(genelist))###316
gene_module_df <- find_gene_modules(cds_branch[pr_deg_ids,], resolution=1e-2, cores = 30)## 将轨迹可变基因收集到模块中：
head(gene_module_df);dim(gene_module_df)#316   5
table(gene_module_df$module)
# 1  2  3  4 
# 92 91 88 45 
write.table(gene_module_df, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_CTBs2STBs_gene_module_along_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)


##moldule expression in different cell types
gene_module <- gene_module_df
cell_group <- tibble::tibble(cell=row.names(colData(cds_branch)), cell_group=colData(cds_branch)$re_annotation_TFac)
agg_mat <- aggregate_gene_expression(cds_branch, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
head(agg_mat)
agg_mat2<-agg_mat[,c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3")]
pheatmap::pheatmap(agg_mat2, scale="column",clustering_method="ward.D2")
module_heatmap<-pheatmap::pheatmap(agg_mat2, scale="column", cluster_cols=F,clustering_method="ward.D2")

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/heatmap_for_average_expression_of_psudotime_module_ALL_in_CTBs2STBs.pdf",width =8,height = 12)
print(module_heatmap)
dev.off()
##根据聚类可以再归为5大类
##We can also pass gene_module_df to plot_cells() as we did when we compared clusters in the L2 data above.
all_module_expression<-plot_cells(cds_branch, genes=gene_module_df ,label_cell_groups=FALSE,show_trajectory_graph=FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_CTBs2STBs_module_expression.png",all_module_expression,width = 20, height =20,limitsize = FALSE)
all_module_expression2<-all_module_expression + scale_color_gradient2(low = "blue", mid = "gray", high =  "red")
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_CTBs2STBs_module_expression2.png",all_module_expression2,width = 20, height =20,limitsize = FALSE)


##绘制变化动态
gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_CTBs2STBs_gene_module_along_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
head(gene_module)
gene_module$id<-as.character(gene_module$id)
gene_module<-gene_module[order(gene_module$module,decreasing = F),]
gap_number0<-as.numeric(table(gene_module$module))
gap_number<-gap_number0[-length(gap_number0)]
head(gene_module)
# 
pt.matrix <- exprs(cds_branch)[match(gene_module$id,rownames(rowData(cds_branch))),]
pt.matrix[1:4,1:4];dim(pt.matrix)#316 18364
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#316 18364
##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
#rownames(pt.matrix) <- pr_DEGs_branch
pt.matrix[1:4,1:4];dim(pt.matrix)#316 18364
#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

htkm <- Heatmap(plot_matrix,name= "z-score",#border = TRUE,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(paste0("Module_",1:length(gap_number0)),gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of CTBs2STBs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_direction_CTBs2STBs_all_module.pdf",height=15,width=8)
print(htkm)
dev.off()
gene_module$module<-factor(gene_module$module,levels = c(2,3,1,4))
gene_module<-gene_module[order(gene_module$module,decreasing = F),]
gene_module$type<-"N_Later"
gene_module[which(gene_module$module %in% c(2)),]$type<-"Early"
gene_module[which(gene_module$module %in% c(3)),]$type<-"Early_Mid"
gene_module[which(gene_module$module %in% c(1)),]$type<-"Mid_Later"
gene_module[which(gene_module$module %in% c(4)),]$type<-"N_Later"

gene_module$type<-factor(gene_module$type,levels = c("Early","Early_Mid","Mid_Later","N_Later"))
#gene_module<-gene_module[order(gene_module$type,decreasing = F),]
gene_module$module<-factor(gene_module$module,levels = c(2,3,1,4))
gene_module<-gene_module[order(gene_module$module,decreasing = F),]
write.table(gene_module, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/developped_genes_rearranged_along_CTBs2STBs_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

#GO BP enrichment
allcluster_BP_GO=data.frame();allcluster_BP_GO_whole=data.frame()
for (cluster_name in levels(gene_module$type)) {
  #cluster_name<-"Mid_low"
  print(paste0("cluster number: ",cluster_name))
  small_gene_group=unique(as.character(gene_module[which(gene_module$type == cluster_name),]$id))
  gene_change=bitr(small_gene_group, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  BP_GO <- enrichGO(gene= unique(gene_change$ENTREZID),OrgDb= org.Hs.eg.db,keyType= 'ENTREZID',ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff= 0.1,readable= TRUE)
  BP_GO2 <- simplify(BP_GO, cutoff=0.7, by="p.adjust", select_fun=min)
  go_res=BP_GO2@result
  go_res_whole=BP_GO@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=cluster_name;allcluster_BP_GO=rbind(allcluster_BP_GO,go_res)}
  if (dim(go_res_whole)[1] != 0) {
    go_res_whole$cluster=cluster_name;allcluster_BP_GO_whole=rbind(allcluster_BP_GO_whole,go_res_whole)}
}
head(allcluster_BP_GO[,c("ID","Description","qvalue","cluster")])
table(allcluster_BP_GO$cluster)
#  Early Early_Mid Mid_Later   N_Later 
#    78        21        49        33 
write.table(allcluster_BP_GO, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_for_gene_module_along_CTBs2STBs_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(allcluster_BP_GO_whole, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_whole_for_module_along_CTBs2STBs_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

###plot modules heatmaps for rearrange
gap_number0<-as.numeric(table(gene_module$type))
gap_number<-gap_number0[-length(gap_number0)]
head(gene_module)
tail(gene_module)

plot_matrix2<-plot_matrix[as.character(gene_module$id),]
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,
                # row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(c("Early_time","Mid_low","Later_time","Mid_high"),gap_number0),
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(c("Early","Early_Mid","Mid_Later","N_Later"),gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of CTBs2STBs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_direction_CTBs2STBs_rearrange_module.pdf",height=15,width=8)
print(htkm)
dev.off()


###modification for heatmaps

cds_branch<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2EVTs_selected_branch_monocle3_object.rds")
pseu_time<-data.frame(pseudotime(cds_branch))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)
#pseu_time$number<-ifelse(rownames(pseu_time) %in% rownames(colData(cds_branch)),1,-1)
pseu_time$cell_name<-rownames(pseu_time)
pseu_time<-pseu_time[order(pseu_time$pseudotime,decreasing = F),]
head(pseu_time);tail(pseu_time)

#reading development genes
gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/developped_genes_rearranged_along_CTBs2EVTs_branchment_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
head(gene_module)
gene_module$module<-factor(gene_module$module,levels = c(8,7,1,6,5,2,3,4,9))
gene_module$type<-factor(gene_module$type,levels = c("Early","Early_Mid","Mid","Mid_Later","N_Later"))
gene_module<-gene_module[order(gene_module$type,gene_module$module,decreasing = F),]
gene_module$id<-as.character(gene_module$id)
head(gene_module)
table(gene_module$type)

#get metadata
branch_meta<-as.data.frame(colData(cds_branch))
head(branch_meta);dim(branch_meta)
#get expression matrix
pt.matrix <- exprs(cds_branch)[match(gene_module$id,rownames(rowData(cds_branch))),]
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 18330
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#1174 18330

##reorder branch cells
branch_meta<-branch_meta[colnames(pt.matrix),]
head(branch_meta);dim(branch_meta)# 18330    67

##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix[1:4,1:4];dim(pt.matrix)# 1174 18330

#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

##get gaps number for row
groups<-names(table(gene_module$type))
gap_number0<-as.numeric(table(gene_module$type))
gap_number<-gap_number0[-length(gap_number0)]

##building expression matrix
plot_matrix2<-plot_matrix[gene_module$id,];dim(plot_matrix2)#577 8232
plot_matrix2[1:4,1:4]

#set annotation for each cells  in column
head(branch_meta);dim(branch_meta)
branch_meta[,c("re_annotation_TFac2","Treat")]
branch_meta$re_annotation_TFac2<-as.character(branch_meta$re_annotation_TFac2)
branch_meta$re_annotation_TFac<-as.character(branch_meta$re_annotation_TFac)
table(branch_meta$re_annotation_TFac)
#branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_2","STBs_1","STBs_2","STBs_3"),"others",branch_meta$re_annotation_TFac)
branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("STBs_1","STBs_2","STBs_3"),"others",
                                        ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_1","CTBs_2"),"CTBs","EVTs"))
table(branch_meta$re_annotation_TFac2)
#CTBs   EVTs others 
#5647  12410    273 

column_colors = list(group=c(Abortion=my_morandi_colors[1],CTRL=my_morandi_colors[21]),
                     cell_types=c(CTBs=my_morandi_colors[10],EVTs=my_morandi_colors[15],others="gray"))
column_ha = HeatmapAnnotation( group=as.character(branch_meta$Treat),
                               cell_types = as.character(branch_meta$re_annotation_TFac2),
                               col = column_colors)

#set annotation for each genes in row
DEGs_class<-gene_module$type
ann_colors = list(DEGs_class=c(Early=ppCor[1],Early_Mid=ppCor[2],Mid=ppCor[3],Mid_Later=ppCor[9],N_Later=ppCor[5]))

##reading abortion related genes
Abortion_related_genes<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Abortion_regulated_gene_list_new.txt",header =T,sep="\t",check.names = F)
head(Abortion_related_genes);dim(Abortion_related_genes)
mark_gene <- as.character(Abortion_related_genes$gene_names);length(mark_gene)
gene_pos <- which(rownames(plot_matrix2) %in% mark_gene)
mark_gene2<-rownames(plot_matrix2[gene_pos,])

row_mark <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene2))
row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)

row_anno2[which(row_anno2$row_genes %in% mark_gene2),]
##plot rearrange heatmap
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,top_annotation = column_ha, 
                right_annotation =row_mark,left_annotation =row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of CTBs2EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_CTBs2EVTs_rearrange_module_annotation_add2.pdf",height=8,width=7)
print(htkm)
dev.off()

##for CTBs2STBs
cds_branch<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_CTBs2STBs_selected_branch_monocle3_object.rds")
pseu_time<-data.frame(pseudotime(cds_branch))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)
#pseu_time$number<-ifelse(rownames(pseu_time) %in% rownames(colData(cds_branch)),1,-1)
pseu_time$cell_name<-rownames(pseu_time)
pseu_time<-pseu_time[order(pseu_time$pseudotime,decreasing = F),]
head(pseu_time);tail(pseu_time)

#reading development genes
gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/developped_genes_rearranged_along_CTBs2STBs_branchment_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
head(gene_module)
table(gene_module$type)

gene_module$module<-factor(gene_module$module,levels =  c(2,3,1,4))
gene_module$type<-factor(gene_module$type,levels = c("Early","Early_Mid","Mid_Later","N_Later"))
gene_module<-gene_module[order(gene_module$type,gene_module$module,decreasing = F),]
gene_module$id<-as.character(gene_module$id)
head(gene_module)

#get metadata
branch_meta<-as.data.frame(colData(cds_branch))
head(branch_meta);dim(branch_meta)
#get expression matrix
pt.matrix <- exprs(cds_branch)[match(gene_module$id,rownames(rowData(cds_branch))),]
pt.matrix[1:4,1:4];dim(pt.matrix)#316 18364
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#316 18364

##reorder branch cells
branch_meta<-branch_meta[colnames(pt.matrix),]
head(branch_meta);dim(branch_meta)# 18330    67

##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix[1:4,1:4];dim(pt.matrix)# 1174 18330

#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

##get gaps number for row
groups<-names(table(gene_module$type))
gap_number0<-as.numeric(table(gene_module$type))
gap_number<-gap_number0[-length(gap_number0)]

##building expression matrix
plot_matrix2<-plot_matrix[gene_module$id,];dim(plot_matrix2)#316 18364
plot_matrix2[1:4,1:4]

#set annotation for each cells  in column
head(branch_meta);dim(branch_meta)
branch_meta[,c("re_annotation_TFac2","Treat")]
branch_meta$re_annotation_TFac2<-as.character(branch_meta$re_annotation_TFac2)
branch_meta$re_annotation_TFac<-as.character(branch_meta$re_annotation_TFac)
table(branch_meta$re_annotation_TFac)
#branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_2","STBs_1","STBs_2","STBs_3"),"others",branch_meta$re_annotation_TFac)
branch_meta$re_annotation_TFac2<-ifelse(branch_meta$re_annotation_TFac %in% c("EVTs_1","EVTs_2","EVTs_3"),"others",
                                        ifelse(branch_meta$re_annotation_TFac %in% c("CTBs_1","CTBs_2"),"CTBs","STBs"))
table(branch_meta$re_annotation_TFac2)
#CTBs   STBs others 
#5647  12410    273 

column_colors = list(group=c(Abortion=my_morandi_colors[1],CTRL=my_morandi_colors[21]),
                     cell_types=c(CTBs=my_morandi_colors[10],STBs=ppCor[5],others="gray"))
column_ha = HeatmapAnnotation( group=as.character(branch_meta$Treat),
                               cell_types = as.character(branch_meta$re_annotation_TFac2),
                               col = column_colors)

#set annotation for each genes in row
DEGs_class<-gene_module$type
ann_colors = list(DEGs_class=c(Early=ppCor[1],Early_Mid=ppCor[2],Mid_Later=ppCor[9],N_Later=ppCor[5]))

##reading abortion related genes
Abortion_related_genes<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Abortion_regulated_gene_list_new.txt",header =T,sep="\t",check.names = F)
head(Abortion_related_genes);dim(Abortion_related_genes)
mark_gene <- as.character(Abortion_related_genes$gene_names);length(mark_gene)
gene_pos <- which(rownames(plot_matrix2) %in% mark_gene)
mark_gene2<-rownames(plot_matrix2[gene_pos,])
mark_gene2
row_mark <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene2))
row_ha = rowAnnotation(DEGs_class=DEGs_class,col = ann_colors)

##plot rearrange heatmap
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,top_annotation = column_ha, 
                right_annotation =row_mark,left_annotation =row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of CTBs2STBs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 4),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_CTBs2STBs_rearrange_module_annotation_add2.pdf",height=8,width=7)
print(htkm)
dev.off()
