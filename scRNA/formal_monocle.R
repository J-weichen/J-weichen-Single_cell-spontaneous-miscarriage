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
library(clusterProfiler)
library(org.Hs.eg.db)
library(reshape2)

pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#note: Monocle2直接输入Seurat object只适用于Seuratv2.0版本的Seurat object。
#Monocle3无法将Seurat Object 转为cds,
# 自己手动构建celldataset,--De novo construct monocle v2 的 celldataset
#Although Monocle can be used with raw read counts, these are not directly proportional to expression values unless you normalize them by length, so some Monocle functions could produce nonsense results. If you don't have UMI counts, We recommend you load up FPKM or TPM values instead of raw read counts.

#step1 ：read Seurat object and load initial coldata and matrix
target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
target$final_major_subgroup_brief_Treat <- paste(target$final_major_subgroup_brief, target$Treat, sep = "_")
table(target$final_major_group_brief)
cell_name<-"Trophoblast_KRT7"
selected_cell0<-subset(x = target, subset = final_major_group_brief == cell_name)
selected_cell<-subset(x = selected_cell0, subset = sample_code == "CTRL_1")

DefaultAssay(selected_cell) <- "RNA"    ## very important
selected_cell <- SCTransform(object = selected_cell,verbose = FALSE) 
top10 <- head(VariableFeatures(selected_cell), 10)
length(VariableFeatures(selected_cell))

selected_cell <- RunPCA(object = selected_cell,verbose = FALSE)
ElbowPlot(object = selected_cell)
selected_cell <- RunUMAP(object = selected_cell, dims = 1:25,verbose = FALSE)
#selected_cell <- FindNeighbors(object = selected_cell, dims = 1:25,verbose = FALSE)
#selected_cell <- FindClusters(object = selected_cell,resolution = 0.6, verbose = FALSE)

DimPlot(object = selected_cell, reduction = "umap", group.by = "final_major_subgroup_brief",cols =ppCor)
#plot_var <- LabelPoints(plot = selected_cell, points = top10, repel = TRUE)
#plot_var
#selected CD4 T cells
#Tcells<-c("CD4_T_Cell_memory_CD8A_neg_S100A4_pos","CD4_T_Cell_naive_CD8A_neg_SELL_pos")
#Ts_blood_seurat<-subset(x = selected_cell, subset = merge_cluster %in% Tcells)
#Ts_blood_seurat$merge_cluster_treat<-paste0(Ts_blood_seurat$merge_cluster,"_",Ts_blood_seurat$Treat)

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


table(selected_cell$monocle_names0)
#CTBs EVTs STBs 
#8266 1671 1033
table(selected_cell$final_major_subgroup_brief)
#CTBs_0 CTBs_1 CTBs_2 CTBs_3 EVTs_1 EVTs_2 EVTs_3 STBs_1 STBs_2 
# 809   1405   5783    269    417    492    762    210    823 

### step0 :导入数据，创建对象，预处理

#1.data@assays$RNA@data
##存放 relative expression values （TPM, FPKM/RPKM）
#2.data@assays$RNA@counts
##存放 absolute transcript counts （TPM, FPKM/RPKM）

#prepare input three files(data, phenotype data, and feature data) from the SeuratObject
#expression matrix(numeric) gene * cellnames
#phenoData: cellnames and other metadata
#featureData: genenames and other metadata(gene_short_name:symbol was need)  
#length of cellnames and genenames  should be the sames and matched as that in matrix 
#for building newCellDataSet of Monocle2
##大数据集使用稀疏矩阵，节省内存，加快运算
#10X的数据使用UMI count矩阵
count_raw <- as(as.matrix(selected_cell@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = selected_cell@meta.data)
fData <- data.frame(gene_short_name = row.names(count_raw), row.names = row.names(count_raw))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(count_raw,phenoData = pd,featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

HSMM<-monocle_cds

#Size factor and dispersion will be used to normalize data and select genes for clustering.
## Monocle的标准化和scale 计算SizeFactor和Dispersions
###size factors可帮助我们对不同细胞的mRNA差异进行标准化（normalize）；
###dispersion值将有助于我们在下面的分析中进行差异表达分析。
##下面这两步只在前面选择了负二项分布的情况下使用。
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#观察到这里多了一个Size_Factor的列
head(pData(HSMM))
#disp_table <- dispersionTable(cds)

#数据清洗 根据自身需求进行

#1）Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))
#2）找出至少在10个细胞里表达的基因
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
length(expressed_genes)#1556
print(head(pData(HSMM)))

#设置上下限：
upper_bound_HSMM <- 10^(mean(log10(pData(HSMM)$nCount_RNA))+ 2*sd(log10(pData(HSMM)$nCount_RNA)))
lower_bound_HSMM <- 10^(mean(log10(pData(HSMM)$nCount_RNA)) - 2*sd(log10(pData(HSMM)$nCount_RNA)))
#可视化
qplot(nCount_RNA,data = pData(HSMM), geom = "density") + geom_vline(xintercept = lower_bound_HSMM) + geom_vline(xintercept = upper_bound_HSMM)
##留下两条竖线中间的部分：
HSMM <- HSMM[,pData(HSMM)$nCount_RNA > lower_bound_HSMM & pData(HSMM)$nCount_RNA < upper_bound_HSMM]

#数据标准化结果展示 #将表达矩阵中所有值进行log标准化
L <- log(exprs(HSMM[expressed_genes,]))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))##将每个基因都标准化，melt方便作图
head(na.omit(melted_dens_df))
###作图，查看标准化的基因表达值的分布
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") + ylab("Density")


##细胞分类（Classifying）;这一步一般在seurat里面已经做完
#CDC20 <- row.names(subset(fData(HSMM), gene_short_name == "CDC20"))
#GABPB2 <- row.names(subset(fData(HSMM), gene_short_name == "GABPB2"))
#cth <- newCellTypeHierarchy()
#cth <- addCellType(cth, "CDC20", classify_func = function(x) { x[CDC20,] >= 1 })
#cth <- addCellType(cth, "GABPB2", classify_func = function(x){ x[GABPB2,] >= 1 })
#HSMM <- classifyCells(HSMM, cth, 0.1)
#table(pData(HSMM)$CellType)
#pie <- ggplot(pData(HSMM),aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
#pie + coord_polar(theta = "y") +theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#确定用于排序的基因 这一步一般也在seurat中做完了
#Clustering cells without marker genes 
#我们可以根据平均表达水平筛选基因，我们还可以选择细胞间异常变异的基因。这些基因往往对细胞状态有很高的信息含量。
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
#红线表示单片基于这种关系对色散的期望。我们标记用于聚类的基因用黑点表示，其他的用灰点表示。
# Plots the percentage of variance explained by the each component based on PCA from the normalized expression data using the same procedure used in reduceDimension function.
# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F) # norm_method='log' time consuming
#step2 reduceDimension: time consuming 
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 10,reduction_method = 'tSNE', verbose = T)
#HSMM <- reduceDimension(HSMM, max_components= 2, num_dim = 3, norm_method = 'log', reduction_method = 'tSNE',verbose = T)

#聚类很少用 可用于减少无趣变量的影响;对细胞进行聚类，在Seurat中采用的是分辨率来确定cluster的数目。而monocle中可以直接指定聚类数目。主要指出的这里所聚类获得的cluster数目比我们赋值的要少一个。即当num_clusters=3时，你只获得了2个cluster
#HSMM <- clusterCells(HSMM, num_clusters = 2)
#plot_cell_clusters(HSMM, 1, 2, color = "CellType",markers = c("CDC20", "GABPB2"))
#HSMM <- clusterCells(HSMM, num_clusters = 5)
#plot_cell_clusters(HSMM, 1, 2)


#正式进行拟时间分析
#step one :define genes for order: 1)DEG or 2) maker genes 选择输入基因用于机器学习
#理想情况下，我们希望尽可能少地使用正在研究的系统生物学的先验知识。
#我们希望从数据中发现重要的排序基因,而不是依赖于文献和教科书，因为这可能会在排序中引入偏见。
#我们将从一种更简单的方法开始，但是我们通常推荐一种更复杂的方法，称为“dpFeature”。
##这个过程称为feature selection（特征选择），这些基因对轨迹的形状有着最重要的影响。
#我们应该要选择的是最能反映细胞状态的基因。

diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~final_major_subgroup_brief")
dim(diff_test_res)
#ordering_genes <- row.names (subset(diff_test_res, qval < 0.1)) ## 不要也写0.1,而是要写0.01。
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
ordering_genes_cal<-as.character(ordering_genes)

#second list gene 也可以使用seurat筛选的hvg：
ordering_genes_seurat<-VariableFeatures(selected_cell)

ordering_genes<-ordering_genes_seurat
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
#cds_ordering_genes <- VariableFeatures(dt) ##这里要注意，如果是subset的数据，要重新计算VariableFeatures
#cds <-setOrderingFilter(cds,ordering_genes = cds_ordering_genes)


#step2 reduceDimension: time consuming
HSMM <- reduceDimension(HSMM, max_components = 2,method="DDRTress")#time consuming
#HSMM <-reduceDimension(HSMM,max_components = 2, reduction_method = "DDRTree",norm_method = 'log',residualModelFormulaStr = "~num_genes_expressed") 
#residualModelFormulaStr减少其他因素的影响，比如不同样本、不同批次

#saveRDS("/mnt/data/chenwei/jianghai/0.script/Trajectory/test_datafile/CDS_DDRTree.rds", file="CDS_DDRTree.rds")

#step three:  order cells along the trajectory   按照伪时间排序
#CDS <- readRDS("CDS_DDRTree.rds")
HSMM <- orderCells(HSMM)
#这个过程称为manifold learning（流形学习）。Monocle利用轨迹来描述细胞如何从一个状态转换到另一个状态。得到的是一个树状图，树的根部或树干表示的是细胞最初的状态（如果有的话），随着细胞的变化就沿着分叉树枝发展。一个细胞的伪时间值（pseudotime value）为它的位置沿着树枝到根部的距离

#step four: plot
colnames(pData(HSMM))
head(pData(HSMM))
#plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief_Treat")
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief")
plot_cell_trajectory(HSMM, color_by = "monocle_names0")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "Phase")
plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "monocle_names0",cell_size = 0.75)+facet_wrap(~monocle_names0, ncol = 2) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief",cell_size = 0.75)+facet_wrap(~final_major_subgroup_brief, ncol = 3) ##根据细胞聚类的进行着色


plot_cell_trajectory(HSMM, color_by = "State",cell_size = 0.75) ####以state进行着色
plot_cell_trajectory(HSMM, color_by = "State",cell_size = 0.75)+facet_wrap(~State, nrow = 2) ##绘制state的分面图
plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 0.75) ##根据拟时间值着色
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief",cell_size = 0.75) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief",cell_size = 0.75)+facet_wrap(~final_major_subgroup_brief, nrow =2) ###绘制clusster的分面图
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief_Treat",cell_size = 0.75) +scale_color_manual(values=ppCor) ##如果有Seurat生的rds文件的话，按照seurat中分的群进行着色，如果不想用ggplot的默认色，可以提供颜色列表col list。
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief_Treat", cell_size = 0.75) + facet_wrap (~final_major_subgroup_brief_Treat, nrow =5) ###按照seurat中分的群绘制分面图。
plot_cell_trajectory(HSMM, color_by = "Treat",cell_size = 0.75) ###按照样本进行着色
plot_cell_trajectory(HSMM, color_by = "Treat",cell_size = 0.75)+facet_wrap(~Treat, ncol = 2)
plot_cell_trajectory(HSMM, color_by = "merge_cluster",cell_size = 0.75)+facet_wrap(~Treat, ncol = 2) ##根据细胞聚类的进行着色
plot_cell_trajectory(test,color_by = "celltype")

##绘制样本着色的分面图。
my_colour <- c('#006dbb','#1cac55','#e52622','#4db3e8', "#c46cac", '#006dbb','#1cac55')
p1 <- plot_cell_trajectory(CDS, color_by = "State", cell_size = 1.25) +
  scale_colour_manual(values=my_colour)
p1
#保存到pdf文件
ggsave(file="/mnt/data/chenwei/jianghai/0.script/Trajectory/Ts_cell_monocle_State.pdf",p1,width = 8, height =8)


#step five:set root state #定义root
#经常拟时序的方向或是根节点弄错了，还需要手动更改

#绝大多数的非成熟的造血干细胞表达E-Slam，所以用这个标准定义轨迹的“root”:
table(pData(cds)$State, pData(cds)$celltype)[,"ESLAM"]
# 1  2  3 
# 0 10  0
#以第3个state作为root 
cds <- orderCells(cds, root_state = 3)
#根据与root的距离定义拟时间
plot_cell_trajectory(cds, color_by = "Pseudotime")

#根据拟时间模型，把每一个branch做差异基因 
#找随拟时序变化的差异基因，以及不同state之间的差异基因。这两个都是monocle里面的概念，与seurat里面找的差异基因不同。
my_pseudotime_de <- differentialGeneTest(HSMM, fullModelFormulaStr = "~sm.ns(Pseudotime)")
my_pseudotime_de <-my_pseudotime_de[order(my_pseudotime_de$qval,decreasing = F),]
#time comsuming
sig_pseudo_gene_names <- row.names(subset(my_pseudotime_de, qval < 0.1))
plot_pseudotime_heatmap(HSMM[sig_pseudo_gene_names[1:50], ], num_clusters =3,show_rownames=TRUE)

my_states_de <- differentialGeneTest(test[expressed_genes,],fullModelFormulaStr = "~State")
my_states_de <- states_de[order(states_de$qval), ]
sig_states_gene_names <- row.names(subset(my_states_de, qval < 0.1))
#plot_pseudotime_heatmap(HSMM[sig_states_gene_names[1:50], ], num_clusters =3,show_rownames=TRUE)
write.table(my_pseudotime_de, file = "pseudotime_de.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(my_states_de, file = "states_de.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
saveRDS(my_pseudotime_de, file = "my_pseudotime_de.rds")
saveRDS(my_states_de, file = "my_states_de.rds")

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

#对目标基因的专门展示
test_genes=c("TCF7","CX3CR1")
pdf("genes_branched_pseudotime.pdf",width = 9,height = 4)
plot_genes_branched_pseudotime(HSMM[test_genes,],branch_point = 1,color_by = "merge_cluster",cell_size=2,ncol = 2)
dev.off()


