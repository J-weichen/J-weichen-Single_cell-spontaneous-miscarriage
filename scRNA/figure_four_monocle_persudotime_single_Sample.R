#Monocle2 Pseudotime analysis
#ref:https://nbisweden.github.io/workshop-scRNAseq/oldlabs/monocle_analysis.html
#ref:https://cloud.tencent.com/developer/article/1692225
#ref: https://www.sohu.com/a/327133920_278730

rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
#dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
#dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40.3.1')
library(Seurat)
library(monocle)
library(scales)
library(ggsci)
library(ggplot2)
library(reshape2)
library(ggbeeswarm)
library(ggridges)

#library(clusterProfiler)
#library(org.Hs.eg.db)

pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)


#step1 ：read Seurat object and load initial coldata and matrix
target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
target$final_major_subgroup_brief_Treat <- paste(target$final_major_subgroup_brief, target$Treat, sep = "_")
table(target$final_major_group_brief)

cell_name<-"Trophoblast_KRT7"
selected_cell0<-subset(x = target, subset = final_major_group_brief == cell_name)
table(selected_cell0$sample_code)
#CTRL_1   CTRL_2   CTRL_3   CTRL_4   CTRL_5 Arrest_1 Arrest_2 Arrest_3 Arrest_4 
# 10970     6138     4120     2309     1033     1212     5208     1776      425 
sample_names<- "Arrest_1"
selected_cell<-subset(x = selected_cell0, subset = sample_code == sample_names)
DefaultAssay(selected_cell) <- "RNA"    ## very important
selected_cell <- SCTransform(object = selected_cell,verbose = FALSE) 
top10 <- head(VariableFeatures(selected_cell), 10)
length(VariableFeatures(selected_cell))

selected_cell <- RunPCA(object = selected_cell,verbose = FALSE)
ElbowPlot(object = selected_cell)
selected_cell <- RunUMAP(object = selected_cell, dims = 1:25,verbose = FALSE)
selected_cell <- FindNeighbors(object = selected_cell, dims = 1:25,verbose = FALSE)
selected_cell <- FindClusters(object = selected_cell,resolution = 0.6, verbose = FALSE)
DimPlot(object = selected_cell, reduction = "umap", group.by = "final_major_subgroup_brief",cols =ppCor)

metadata_new<-selected_cell@meta.data
dim(metadata_new)
table(metadata_new$final_major_subgroup_brief)
names(table(metadata_new$final_major_subgroup_brief))
metadata_new$monocle_names0<-metadata_new$final_major_subgroup_brief
metadata_new[which(metadata_new$final_major_subgroup_brief %in% c("CTBs_0","CTBs_1","CTBs_2","CTBs_3")),]$monocle_names0 <- "CTBs"
metadata_new[which(metadata_new$final_major_subgroup_brief %in% c("EVTs_1","EVTs_2","EVTs_3")),]$monocle_names0 <- "EVTs"
metadata_new[which(metadata_new$final_major_subgroup_brief %in% c("STBs_1","STBs_2")),]$monocle_names0 <- "STBs"
head(metadata_new)
selected_cell@meta.data<-metadata_new
saveRDS(selected_cell, file = paste("/mnt/data/chenwei/gongchen/6.monocle2/",sample_names,"_Trophoblast_seurat_pbject.rds"))

table(selected_cell$monocle_names0)
table(selected_cell$final_major_subgroup_brief)

##call DEGs inter group
Idents(object = selected_cell) <- "final_major_subgroup_brief" #'merge_cluster'
Allcluster.markers_all <- FindAllMarkers(object = selected_cell,min.pct = 0.25)
Allcluster.markers_all_default0.25<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=0.25,]
Allcluster.markers_all_default0.25_q01<-Allcluster.markers_all_default0.25[Allcluster.markers_all_default0.25$p_val_adj<0.1,]
DEGs_seurat_call<-Allcluster.markers_all_default0.25_q01

##step2 ：Construct CellDataSet object 10X的数据使用UMI count矩阵 
count_raw <- as(as.matrix(selected_cell@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = selected_cell@meta.data)
fData <- data.frame(gene_short_name = row.names(count_raw), row.names = row.names(count_raw))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(count_raw,phenoData = pd,featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

##step3 Size factor and dispersion will be used to normalize data
HSMM<-monocle_cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

##step4 Filtering gene and cells 
#Genes that aren’t highly expressed enough will not be used for clustering, since they may not give meaningful signal and would only add noise.
#1）Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 3 );print(head(fData(HSMM)))
#2）找出至少在10个细胞里表达的基因
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
length(expressed_genes)#5603
print(head(pData(HSMM)))
#设置上下限：
upper_bound_HSMM <- 10^(mean(log10(pData(HSMM)$nCount_RNA))+ 2*sd(log10(pData(HSMM)$nCount_RNA)))
lower_bound_HSMM <- 10^(mean(log10(pData(HSMM)$nCount_RNA)) - 2*sd(log10(pData(HSMM)$nCount_RNA)))
#可视化
qplot(nCount_RNA,data = pData(HSMM), geom = "density") + geom_vline(xintercept = lower_bound_HSMM) + geom_vline(xintercept = upper_bound_HSMM)
##留下两条竖线中间的部分：
HSMM <- HSMM[,pData(HSMM)$nCount_RNA > lower_bound_HSMM & pData(HSMM)$nCount_RNA < upper_bound_HSMM]

#step 5 set genes list for ordering cell 
# Selection of ordering genes can be done in a variety of ways. E.g.:
  ##1) Differential genes
  ##2) High dispersion genes
  ##3) Top PC loading genes
  ##4) Unsupervised feature selection based on density peak clustering
  ##5) Semi-supervised ordering with known marker genes

#Typically, to order cells by progress, we want to reduce the number of genes analyzed. So we can select for a subset of genes that we believe are important in setting said ordering, such as overdispersed genes.

##1) variable gene using top dispersed genes for pseudotime ordering
disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
#unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)$gene_id
disp_HSMM <- setOrderingFilter(HSMM,ordering_genes)
plot_ordering_genes(disp_HSMM)
disp_HSMM <- reduceDimension(disp_HSMM, max_components = 2, method = 'DDRTree')
#HSMM <-reduceDimension(HSMM,max_components = 2, reduction_method = "DDRTree",norm_method = 'log',residualModelFormulaStr = "~num_genes_expressed") 
#residualModelFormulaStr减少其他因素的影响，比如不同样本、不同批次

#HSMM <- reduceDimension(HSMM, num_dim = 20,reduction_method = 'tSNE', verbose = T)
disp_HSMM <- orderCells(disp_HSMM)
plot_cell_trajectory(disp_HSMM, color_by = "final_major_subgroup_brief",cell_size = 0.75)+facet_wrap(~final_major_subgroup_brief, nrow =2) ###绘制clusster的分面图

##1) deg called by monocle2
#Genes likely to be informative of ordering of cells along the pseudotime trajectory will be selected for pseudotime inference.
diff_genes <- differentialGeneTest(HSMM, fullModelFormulaStr = "~ Cluster + final_major_subgroup_brief")# Use top 3000 differentially expressed genes
ordering_genes <- row.names(subset(diff_genes, qval < 1e-3))[order(diff_genes$qval)][1:3000]
diff_HSMM <- setOrderingFilter(HSMM, ordering_genes)
diff_HSMM <- reduceDimension(diff_HSMM, max_components = 2, method = 'DDRTree')#Here Monocle 2 will first project the data to 2 dimensions with DDRTree, and then do trajectory inference (orderCells).
diff_HSMM <- orderCells(diff_HSMM)

plot_cell_trajectory(diff_HSMM, color_by = "final_major_subgroup_brief",cell_size = 0.75)+facet_wrap(~final_major_subgroup_brief, nrow =2) ###绘制clusster的分面图

#hvg called by seurat：
hvg_genes<-VariableFeatures(selected_cell)
ordering_genes<-hvg_genes
hvg_HSMM <- setOrderingFilter(HSMM, ordering_genes)
hvg_HSMM <- reduceDimension(hvg_HSMM, max_components = 2, method = 'DDRTree')
#HSMM <- reduceDimension(HSMM, num_dim = 20,reduction_method = 'tSNE', verbose = T)
hvg_HSMM <- orderCells(hvg_HSMM)
plot_ordering_genes(hvg_HSMM)

#DEGs called by seurat：
ordering_genes_seurat<-DEGs_seurat_call
ordering_genes<-ordering_genes_seurat$gene
length(ordering_genes)
DEGs_HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(DEGs_HSMM)
DEGs_HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
#HSMM <- reduceDimension(HSMM, num_dim = 20,reduction_method = 'tSNE', verbose = T)
DEGs_HSMM <- orderCells(DEGs_HSMM)

##可进行另一类降维--选作
HSMM <- reduceDimension(HSMM, num_dim = 20,reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, method = "louvain")

plot_cell_clusters(HSMM, cell_size = 0.5) +theme(legend.position = "none") +labs(x = "tSNE1", y = "tSNE2")
plot_cell_clusters(HSMM, cell_size = 0.5, color_by = "final_major_subgroup_brief") +
  scale_colour_manual(values=ppCor)+
  labs(x = "tSNE1", y = "tSNE2") +theme(legend.position = "right") +guides(color = guide_legend(override.aes = list(size = 3)))
plot_cell_clusters(HSMM, cell_size = 0.5, color_by = "monocle_names0") +
  scale_colour_manual(values=ppCor)+
  labs(x = "tSNE1", y = "tSNE2") +theme(legend.position = "right") +guides(color = guide_legend(override.aes = list(size = 3)))

#selected special HSMM
#DEGs_HSMM hvg_HSMM diff_HSMM disp_HSMM
HSMM<-hvg_HSMM
#Step6 See what the trajectory looks like. for pre root identify
#Monocle does not know which the root state is for the tree, we need to define the root.
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="State")+scale_colour_manual(values=rev(ppCor))
plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")

plot_cell_trajectory(HSMM, color_by = "monocle_names0",cell_size = 0.75)+facet_wrap(~monocle_names0, ncol = 2) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief",cell_size = 0.75)+scale_colour_manual(values=ppCor)+facet_wrap(~final_major_subgroup_brief, ncol = 3) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)+facet_wrap(~Phase, ncol = 2)

table(pData(HSMM)$final_major_subgroup_brief,pData(HSMM)$State)
# we clearly see that State1 should be the root # so now we can reorder the cells  again
plot_cell_trajectory(HSMM, color_by="State")+scale_colour_manual(values=rev(ppCor))

HSMM <- orderCells(HSMM, root_state=2)
# plot with coloring by pseudotime
plot_cell_trajectory(HSMM, color_by="Pseudotime")


#We can compare the inferred pseudotime to the known sampling timepoints.
meta_plot<-selected_cell@meta.data;head(meta_plot);dim(meta_plot)
pData_plot<-pData(HSMM);head(pData_plot);dim(pData_plot)
meta_plot2<-merge(meta_plot,pData_plot[,c("Pseudotime","State")],by=0)
dim(meta_plot2)
ggplot(meta_plot2,aes(x = Pseudotime,y = final_major_subgroup_brief, colour = final_major_subgroup_brief)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = ppCor) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("cell subtypes") +
  ggtitle("Cells ordered by monocle2 pseudotime")

plotdf=pData(HSMM)
#monocle 山脊图绘制
ggplot(plotdf, aes(x=Pseudotime,y=final_major_subgroup_brief,fill=final_major_subgroup_brief))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
ggsave("tmp1.pdf",width = 13,height = 7,units = "cm")
#添加颜色
ggplot(plotdf, aes(x=Pseudotime,y=final_major_subgroup_brief,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+theme(panel.grid = element_blank())
ggsave("tmp2.pdf",width = 13,height = 7,units = "cm")
#monocle 坐标提取后绘图
plotdf2=as.data.frame(t(HSMM@reducedDimS))
colnames(plotdf2)=c("component1","component2")
plotdf2$Pseudotime=HSMM$Pseudotime

plotdf2%>%ggplot(aes(component1,component2,color=Pseudotime))+
  geom_point()+theme_minimal()
ggsave("tmp3.pdf",width = 13,height = 10,units = "cm")


##Find genes that change along pseudotime
#Now we run a new differential gene test which test for changes along pseudotime.
#找随拟时序变化的差异基因，以及不同state之间的差异基因。
#这两个都是monocle里面的概念，与seurat里面找的差异基因不同。
#nocle中differentialGeneTest()函数可以按条件进行差异分析，
#我们可以按一定的条件筛选基因后进行差异分析，全部基因都输入会耗费比较长的时间。建议使用cluster差异基因或高变基因输入函数计算。
my_pseudotime_de <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)")
my_pseudotime_de <-my_pseudotime_de[order(my_pseudotime_de$qval,decreasing = F),]
#time comsuming
write.table(my_pseudotime_de, file = "/mnt/data/chenwei/gongchen/6.monocle2/pseudotime_de.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
saveRDS(my_pseudotime_de, file = "/mnt/data/chenwei/gongchen/6.monocle2/my_pseudotime_de.rds")
##Plot some top genes along pseudotime
# take top 8 genes
plotgenes <- rownames(my_pseudotime_de[order(my_pseudotime_de$qval)[1:8],])
plot_genes_in_pseudotime(HSMM[plotgenes,], color_by="Stage",ncol=2)

#Clustering genes by pseudotemporal expression pattern
sig_pseudo_gene_names <- row.names(subset(my_pseudotime_de, qval < 0.01))# qval < 1e-5
length(sig_pseudo_gene_names)
plot_pseudotime_heatmap(HSMM[sig_pseudo_gene_names[1:50], ], num_clusters =3,show_rownames=TRUE)

#状态点差异
my_states_de <- differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~State")
my_states_de <- my_states_de[order(my_states_de$qval), ]
sig_states_gene_names <- row.names(subset(my_states_de, qval < 0.1))
#plot_pseudotime_heatmap(HSMM[sig_states_gene_names[1:50], ], num_clusters =3,show_rownames=TRUE)
write.table(my_states_de, file = "/mnt/data/chenwei/gongchen/6.monocle2/states_de.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
saveRDS(my_states_de, file = "/mnt/data/chenwei/gongchen/6.monocle2/my_states_de.rds")


#对第一分支点进行分析
BEAM_res <- BEAM(HSMM,branch_point = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
plot_genes_branched_heatmap(HSMM[row.names(BEAM_res)[1:50]], branch_point = 1, num_clusters = 3, use_gene_short_name=TRUE, show_rownames=TRUE)
my_branched_heatmap <- plot_genes_branched_heatmap (HSMM[row.names (subset (BEAM_res, qval < 1e-4)), ], branch_point = 1, num_clusters = 4,use_gene_short_name = TRUE,show_rownames = TRUE,return_heatmap = TRUE)
saveRDS(BEAM_res, file = "my_BEAM_res.rds")

#优化热图绘制 以及 分组基因提取
tmp1=plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 4, #这些基因被分成几个group
                                 cores = 1,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, #默认值
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #是否返回一些重要信息
)

pdf("branched_heatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()
#提取不同classes gene
gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)
#富集分析

allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene= unique(df_name$ENTREZID), OrgDb= org.Hs.eg.db,keyType = 'ENTREZID',ont = "BP",pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,qvalueCutoff  = 0.2,readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])

#对Monocle基因可视化 目标基因的专门展示
test_genes=c("TCF7","CX3CR1")
pdf("genes_branched_pseudotime.pdf",width = 9,height = 4)
plot_genes_branched_pseudotime(HSMM[test_genes,],branch_point = 1,color_by = "merge_cluster",cell_size=2,ncol = 2)
dev.off()

s.genes <- c("ITGB1","CCR7","KLRB1","GNLY")
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("pseudotime/genes_visual.png", plot = plotc, width = 8, height = 4.5)


#highlight 目标细胞
my_colour <- c('#006dbb','#1cac55','#e52622','#4db3e8', "#c46cac", '#006dbb','#1cac55')
p1 <- plot_cell_trajectory(CDS, color_by = "State", cell_size = 1.25) +
  scale_colour_manual(values=my_colour)