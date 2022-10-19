rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

require("Matrix")
library(data.table)
library(Seurat)
library(DropletUtils)
library(DoubletFinder)

library(cowplot)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(data.table)
library(grid)
library(gridExtra)
library(cowplot)
#library(ggpubr)
grid.newpage()

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
show_col(Cells_col_raw)

##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)
pal2<-pal_jama("default",alpha = 1)(7)
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)
pal5 <- pal_npg("nrc", alpha=0.5)(10)
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
show_col(ppCor_all)

#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
library(future)
plan(strategy = "multicore", workers = 8)
options(future.globals.maxSize = 1500 * 1024^12)

#reading quality_evaluation_matrix 
sample_qualify_matrix <- read.csv("/mnt/data/chenwei/gongchen/3.seurat_result/sample_qualify_evaluated_matrix.csv")
head(sample_qualify_matrix)
str(sample_qualify_matrix)
sample_qualify_matrix$sample_name<-as.character(sample_qualify_matrix$sample_name)

mydata_list<-c("N1_out","N2_out","N3_out","N4_out","N5_out","A1_out","A3_out","A4_out","A8_out")
#data_name<-c()
#seurat_data_list<-list();seurat_data_list2<-list()
for ( sample_name in mydata_list){
#read 10X Raw matrix for evaluation
#sample_name<-"A3_out"
sample<- unlist(lapply(strsplit(sample_name,"_"), function(x) x[1]))
print(as.character(sample))

sce <- read10xCounts(paste0("/mnt/data/chenwei/gongchen/2.map_result/",sample_name,"/outs/filtered_feature_bc_matrix/"))
#evaluation by other tools
#https://davetang.org/muse/2018/08/09/getting-started-with-cell-ranger/
br.out <- barcodeRanks(counts(sce))
br.out.df <- as.data.frame(br.out)
br.out.df$barcode <- colData(sce)$Barcode
br.out.df %>%  filter(rank <= 10) %>% arrange(rank)
x_knee <- br.out.df %>% filter(total > metadata(br.out)$knee) %>% arrange(total) %>% select(rank) %>% head(1)
x_inflection <- br.out.df %>% filter(total > metadata(br.out)$inflection) %>% arrange(total) %>% select(rank) %>% head(1)
padding <- length(br.out$rank) / 10

Umi_rank_plot<-ggplot(br.out.df, aes(x = rank, y = total)) +
  geom_point() + scale_x_log10() +scale_y_log10() +theme_bw() +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size = 14), title = element_text(size = 16)) +
  geom_hline(yintercept = metadata(br.out)$knee, linetype = 2, colour = "dodgerblue") +
  geom_hline(yintercept = metadata(br.out)$inflection, linetype = 2, colour = "forestgreen") +
  labs(x = "Rank", y = "Total", title =paste0("Barcode Rank vs Total UMI in ",sample)) +
  annotate("text", label = paste0("Knee (", x_knee$rank, ")"), x = x_knee$rank + padding, y = metadata(br.out)$knee, size = 5) +
  annotate("text", label = paste0("Inflection (", x_inflection$rank, ")"), x = x_inflection$rank + padding, y = metadata(br.out)$inflection, size = 5)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/Umi_rank_",sample,"_plot.pdf"), Umi_rank_plot,width=15, height=10)

#br.out <- barcodeRanks(my.counts)
# Making a other plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), legend=c("knee", "inflection"))

length(na.omit(br.out$rank))#9632
length(na.omit(br.out$fitted))#8228
metadata(br.out)$inflection#1298
metadata(br.out)$knee#4984

#read 10X data
dataset_for10X <- Read10X(data.dir = paste0("/mnt/data/chenwei/gongchen/2.map_result/",sample_name,"/outs/filtered_feature_bc_matrix/"))
Cell <- CreateSeuratObject(counts = dataset_for10X, project = "Gongchen_filter_by_hand", min.cells = 3, min.features = 200)
Cell #21366 features across 9504 samples

##metadata built
head(Cell@meta.data)
#                            orig.ident nCount_RNA nFeature_RNA
#AAACCCAAGGTTGTTC-1 Gongchen_filter_by_hand       4415         1686

#Pre-View the data quality for un-fielter Cells
gene_List <- rownames(Cell)
grep(pattern = "^MT-", x = rownames(Cell), value = TRUE)#13 MT genes
length(grep(pattern = "^RPS", x = rownames(Cell), value = TRUE))#46 RPS genes
#mitochondrial genes percent
DefaultAssay(object = Cell)#RNA
#Cell[["percent.mt"]] <- PercentageFeatureSet(object = Cell, pattern = "^MT-")
##There are a few confusion in here
##It seems that I have to pre-perform remove the no-SNG cell in my dataset,or I will get a percent more than 100

MT_features <-  grep(pattern =  "^MT-", x = rownames(x = Cell[["RNA"]]),value = TRUE)
Cell$nCount_MT<-colSums(x = GetAssayData(object = Cell,assay = "RNA", slot = "counts")[MT_features, , drop = FALSE])
meta.data_zero<-Cell@meta.data
meta.data_zero$percent.mt<-100*meta.data_zero$nCount_MT/meta.data_zero$nCount_RNA
Cell@meta.data<-meta.data_zero
#Cell$percent.mt<-100*Cell$nCount_MT/Cell$nCount_RNA
range(Cell$percent.mt)#0.00000 96.08602
length(which(Cell$nCount_MT > Cell$nCount_RNA))#0

#rRNA genes percent
#Cell[["percent.rps"]] <- PercentageFeatureSet(Cell, pattern = "^RPS")
Rps_features <-  grep(pattern =  "^RPS", x = rownames(x = Cell[["RNA"]]),value = TRUE)
Cell$nCount_Rps<-colSums(x = GetAssayData(object = Cell,assay = "RNA", slot = "counts")[Rps_features, , drop = FALSE])
#Cell$percent.rps<-100*Cell$nCount_Rps/Cell$nCount_RNA
meta.data_zero<-Cell@meta.data
meta.data_zero$percent.rps<-100*meta.data_zero$nCount_Rps/meta.data_zero$nCount_RNA
Cell@meta.data<-meta.data_zero
#number of genes detected per UMI: this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
# Add number of genes per UMI for each cell to metadata
meta.data_zero<-Cell@meta.data
meta.data_zero$log10GenesPerUMI <- log10(meta.data_zero$nFeature_RNA) / log10(meta.data_zero$nCount_RNA)
Cell@meta.data<-meta.data_zero
range(Cell$log10GenesPerUMI)#0.6094434 0.9636129
range(Cell$percent.rps)#0.4207574 28.4694153
length(which(Cell$nCount_Rps > Cell$nCount_RNA))
head(meta.data_zero)
#QC by hand
#https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
# Visualize the number UMIs/transcripts per cell
P1<-meta.data_zero %>% 
  ggplot(aes(x=nCount_RNA)) +geom_density(alpha = 0.5,color = "pink", lwd = 2,linetype = 2,fill="#69b3a2") + 
  scale_x_log10() + theme_classic() +ylab("Cell density") + geom_vline(xintercept = 500)
# Visualize the number fearture per cell
P2<-meta.data_zero %>% 
  ggplot(aes(x=nFeature_RNA)) +geom_density(alpha = 0.5,color = "pink", lwd = 2,linetype = 2,fill="blue") + 
  scale_x_log10() + theme_classic() +ylab("Cell density") + geom_vline(xintercept = 500)
P3<-meta.data_zero %>% 
  ggplot(aes(x=percent.mt)) +geom_density(alpha = 0.5,color = "pink", lwd = 2,linetype = 2,fill="orange") + 
  scale_x_log10() + theme_classic() +ylab("Cell density") + geom_vline(xintercept = 10)
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
P4<-meta.data_zero %>% 
  ggplot(aes(x=log10GenesPerUMI)) +geom_density(alpha = 0.5,color = "pink", lwd = 2,linetype = 2,fill="grey") + 
  scale_x_log10() + theme_classic() +ylab("Cell density") + geom_vline(xintercept = 0.8)

combined_4_plot<-P1+P2+P3+P4
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/raw_evaluation_feature_",sample,"_plot.pdf"), combined_4_plot,width=10, height=10)

combined_4_plot2<-VlnPlot(object = Cell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps"), ncol = 2,pt.size = 0.01,cols = ppCor)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/raw_evaluation_feature_",sample,"_plot2.pdf"), combined_4_plot2,width=10, height=10)

#VlnPlot(object = Cell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps"),ncol = 2,pt.size = 0,cols = ppCor,group.by ="Tissue",split.by = "Age_group")
#在V3中 两两指标的相关性
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
opar<-par(no.readonly = TRUE)
plot1 <- FeatureScatter(object = Cell, feature1 = "nCount_RNA", feature2 = "percent.mt",cols =ppCor[2])+ NoLegend()
plot2 <- FeatureScatter(object = Cell, feature1 = "nFeature_RNA", feature2 = "percent.mt",cols =ppCor[3])+ NoLegend()
plot3 <- FeatureScatter(object = Cell, feature1 = "nCount_RNA", feature2 = "percent.rps",cols =ppCor[4])+ NoLegend()
plot4 <- FeatureScatter(object = Cell, feature1 = "nFeature_RNA", feature2 = "percent.rps",cols =ppCor[5])+ NoLegend()
plot5 <- FeatureScatter(object = Cell, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols =ppCor[6])+ NoLegend()
Combine_five_plot<-CombinePlots(plots = list(plot1,plot2,plot3, plot4,plot5),ncol = 3)
par(opar)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/raw_evaluation_feature_",sample,"_plot3.pdf"), Combine_five_plot,width=15, height=10)

#selected Cells  ##也可以保存去除线粒体后的rds,效果一样
#limitation
top_high_UMI<-sample_qualify_matrix[which(sample_qualify_matrix$sample_name == sample),]$Cell_idebtify_100.   
meta.data_zero2<-meta.data_zero[order(meta.data_zero$nCount_RNA,decreasing = T),]
Cell_top_high_UMI<-rownames(meta.data_zero2)[1:top_high_UMI]
Cell_new <- subset(Cell, cells = Cell_top_high_UMI)
Cell_new
#21366 features across 3208 samples
opar<-par(no.readonly = TRUE)
plot1 <- FeatureScatter(object = Cell_new, feature1 = "nCount_RNA", feature2 = "percent.mt",cols =ppCor[2])+ NoLegend()
plot2 <- FeatureScatter(object = Cell_new, feature1 = "nFeature_RNA", feature2 = "percent.mt",cols =ppCor[3])+ NoLegend()
plot3 <- FeatureScatter(object = Cell_new, feature1 = "nCount_RNA", feature2 = "percent.rps",cols =ppCor[4])+ NoLegend()
plot4 <- FeatureScatter(object = Cell_new, feature1 = "nFeature_RNA", feature2 = "percent.rps",cols =ppCor[5])+ NoLegend()
plot5 <- FeatureScatter(object = Cell_new, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols =ppCor[6])+ NoLegend()
Combine_five_plot<- CombinePlots(plots = list(plot1,plot2,plot3, plot4,plot5),ncol = 3)
par(opar)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/filter_evaluation_feature_",sample,"_plot3.pdf"), Combine_five_plot,width=15, height=10)
P1<-Cell_new@meta.data %>% 
  ggplot(aes(x=nCount_RNA)) +geom_density(alpha = 0.5,color = "pink", lwd = 2,linetype = 2,fill="#69b3a2") + 
  scale_x_log10() + theme_classic() +ylab("Cell density") + geom_vline(xintercept = 500)
# Visualize the number fearture per cell
P2<-Cell_new@meta.data %>% 
  ggplot(aes(x=nFeature_RNA)) +geom_density(alpha = 0.5,color = "pink", lwd = 2,linetype = 2,fill="blue") + 
  scale_x_log10() + theme_classic() +ylab("Cell density") + geom_vline(xintercept = 500)
P3<-Cell_new@meta.data %>% 
  ggplot(aes(x=percent.mt)) +geom_density(alpha = 0.5,color = "pink", lwd = 2,linetype = 2,fill="orange") + 
  scale_x_log10() + theme_classic() +ylab("Cell density") + geom_vline(xintercept = 10)
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
P4<-Cell_new@meta.data %>% 
  ggplot(aes(x=log10GenesPerUMI)) +geom_density(alpha = 0.5,color = "pink", lwd = 2,linetype = 2,fill="grey") + 
  scale_x_log10() + theme_classic() +ylab("Cell density") + geom_vline(xintercept = 0.8)

combined_4_plot<-P1+P2+P3+P4
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/filter_evaluation_feature_",sample,"_plot.pdf"), combined_4_plot,width=10, height=10)

combined_4_plot2<-VlnPlot(object = Cell_new, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps"), ncol = 2,pt.size = 0.01,cols = ppCor)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/filter_evaluation_feature_",sample,"_plot2.pdf"), combined_4_plot2,width=10, height=10)

Cell_new_final<- subset(x = Cell_new, subset=  percent.mt < 20)
Cell_new_final #21366 features across 3180 samples 
#Cell_new2<- subset(x = Cell_new,  subset= (nCount_RNA >= 500) & (nCount_RNA < 70000)&
#                                  (nFeature_RNA >= 200) &  (nFeature_RNA <7000) & (percent.mt < 10))
#Cell_new2 # 8325  samples

##genes which are expressed in 10 or more cells.
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = Cell_new_final, slot = "counts")
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 3 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 3
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
Cell_new_final2 <- CreateSeuratObject(filtered_counts, meta.data = Cell_new_final@meta.data)
Cell_new_final2
#20555 features across 3180 samples within 1 assay 

#other strict filter in not forced selected cells
#quantile(meta.data_zero$nCount_RNA,0.95)
#nUMI > 500 #nGene > 200 #log10GenesPerUMI > 0.8 #mitoRatio < 0.2
noana_cell_rm_mf_mix <- subset(x = Cell, subset = percent.mt > 20)#  323 samples 
noana_cell_rm_mf_mix
noana_cell <- subset(x = Cell, subset = nFeature_RNA <500)# 943 samples 
noana_cell
noana_cell <- subset(x = Cell, subset = nFeature_RNA > 4000)# 1225 samples 
noana_cell
noana_cell <- subset(x = Cell, subset = nCount_RNA > 20000)#  1177 samples 
noana_cell
Cell2 <- subset(x = Cell,  subset=(nCount_RNA >= 500) & (nCount_RNA < 70000)&(nFeature_RNA >= 200) &  (nFeature_RNA <7000) & (percent.mt < 20))
Cell2##20651  samples
saveRDS(Cell2, file = paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_no_forced_selected_before_double_finder_Cell_raw.rds"))
saveRDS(Cell_new_final2, file = paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_selected_before_double_finder_Cell_raw.rds"))

#####多细胞检测###
#perform double or multiple cells evaluation by DoubletFinder

counts <- GetAssayData(Cell_new_final2, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(MT_features,Rps_features))),]
seurat_data2 <- subset(Cell_new_final2, features = rownames(counts))

seurat_data2 <- SCTransform(object = seurat_data2,verbose = FALSE) 
seurat_data2 <- RunPCA(object = seurat_data2,verbose = FALSE)
ElbowPlot(object = seurat_data2)
seurat_data2 <- RunUMAP(object = seurat_data2, dims = 1:25,verbose = FALSE)
seurat_data2<- RunTSNE(seurat_data2, dims = 1:25)
plot0<-DimPlot(seurat_data2, reduction = "pca")
plot1<-DimPlot(seurat_data2, reduction = "umap")
plot2<-DimPlot(seurat_data2, reduction = "tsne")
plot_merge_three<-CombinePlots(plots = list(plot0,plot1, plot2),legend="top",ncol = 3)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_no_filter_reduce_dim_plot.pdf"),plot_merge_three,width=12, height=4)

#Seurat采用的是graph-based聚类方法，k-means方法在V3中已经不存在了
seurat_data2 <- FindNeighbors(object = seurat_data2, dims = 1:25,verbose = FALSE)
##different resolution
seurat_data2 <- FindClusters(object = seurat_data2, verbose = FALSE) #default=0.8
#Cell <- FindClusters(object = Cell, resolution = 1.5,verbose = FALSE
seurat_data2 <- FindClusters(object = seurat_data2, resolution = 5,verbose = FALSE)
seurat_data2 <- FindClusters(object = seurat_data2, resolution = 3,verbose = FALSE)
seurat_data2 <- FindClusters(object = seurat_data2, resolution = 2.5,verbose = FALSE)
seurat_data2 <- FindClusters(object = seurat_data2, resolution = 2,verbose = FALSE)
seurat_data2 <- FindClusters(object = seurat_data2, resolution = 1.5,verbose = FALSE)
seurat_data2 <- FindClusters(object = seurat_data2, resolution = 1.2,verbose = FALSE)
seurat_data2 <- FindClusters(object = seurat_data2, resolution = 1,verbose = FALSE)
seurat_data2 <- FindClusters(object = seurat_data2, verbose = FALSE)
seurat_data2 <- FindClusters(object = seurat_data2,resolution = 0.6, verbose = FALSE)
seurat_data2 <- FindClusters(object = seurat_data2, resolution = 0.4,verbose = FALSE)

plot1<-DimPlot(seurat_data2, group.by ="SCT_snn_res.0.8",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.0.8') 
plot2<-DimPlot(seurat_data2, group.by ="SCT_snn_res.1",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.1') 
plot3<-DimPlot(seurat_data2, group.by ="SCT_snn_res.1.2",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.1.2') 
plot4<-DimPlot(seurat_data2, group.by ="SCT_snn_res.1.5",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.1.5') 
plot5<-DimPlot(seurat_data2, group.by ="SCT_snn_res.2",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.2') 
plot6<-DimPlot(seurat_data2, group.by ="SCT_snn_res.2.5",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.2.5') 
plot7<-DimPlot(seurat_data2, group.by ="SCT_snn_res.3",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.3') 
plot8<-DimPlot(seurat_data2, group.by ="SCT_snn_res.5",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.5') 
plot_merge_eight<-CombinePlots(plots = list(plot1, plot2,plot3,plot4,plot5, plot6,plot7,plot8),ncol = 4,legend = NULL)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_no_filter_diff_res.pdf"),plot_merge_eight,width=16, height=8)

plot9<-DimPlot(seurat_data2, group.by ="SCT_snn_res.0.4",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.0.4') 
plot10<-DimPlot(seurat_data2, group.by ="SCT_snn_res.0.6",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.0.6') 
plot_merge<-CombinePlots(plots = list(plot9, plot10),ncol = 2,legend = NULL)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_no_filter_diff_res2.pdf"),plot_merge,width=10, height=5)

Idents(object = seurat_data2) <- "SCT_snn_res.0.6"
seurat_data2$seurat_cluster <- as.character(Idents(seurat_data2))

# 找最佳 PK
#sct 取决于你走的 Seurat 流程是用 NormalizeData() + FindVariableFeatures() + ScaleData() 还是 SCTransform()
sweep.res.list_sub_cell <- paramSweep_v3(seurat_data2, PCs = 1:25, sct = T)
#head(sweep.res.list_sub_cell)
sweep.stats_sub_cell <- summarizeSweep(sweep.res.list_sub_cell, GT = FALSE)
#head(sweep.res.list_sub_cell)
bcmvn_sub_cell <- find.pK(sweep.stats_sub_cell)
mpK<-as.numeric(as.vector(bcmvn_sub_cell$pK[which.max(bcmvn_sub_cell$BCmetric)]))
# 找最佳 nExp
#commend by Cell paper data
##Exp_num<- (0.0007*ncol(sub_cell@assays$SCT@data)-0.127)/100
Exp_num <-ifelse(ncol(seurat_data2@assays$SCT@data)<1000,0.008,ifelse(ncol(seurat_data2@assays$SCT@data)<2000,0.016,
                                                                      ifelse(ncol(seurat_data2@assays$SCT@data)<3000,0.023,
                                                                             ifelse(ncol(seurat_data2@assays$SCT@data)<4000,0.031,
                                                                                    ifelse(ncol(seurat_data2@assays$SCT@data)<5000,0.039,
                                                                                           ifelse(ncol(seurat_data2@assays$SCT@data)<6000,0.046,
                                                                                                  ifelse(ncol(seurat_data2@assays$SCT@data)<7000,0.054,
                                                                                                         ifelse(ncol(seurat_data2@assays$SCT@data)<8000,0.061,
                                                                                                                ifelse(ncol(seurat_data2@assays$SCT@data)<9000,0.069,
                                                                                                                       ifelse(ncol(seurat_data2@assays$SCT@data)<10000,0.076,0.1))))))))))
annotations <- seurat_data2@meta.data$seurat_cluster
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)    
nExp_poi <- round(Exp_num*ncol(seurat_data2@assays$SCT@data))  
#nExp_poi <- round(0.075*ncol(sub_cell@assays$SCT@data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#final # 找 Doublet
# sub_cell <- doubletFinder_v3(sub_cell, PCs = 1:25, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
seurat_data2 <- doubletFinder_v3(seurat_data2, PCs = 1:25, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE,sct = T)
seurat_data2$DF_hi.lo<- seurat_data2@meta.data[,paste("DF.classifications_0.25",mpK,nExp_poi.adj,sep="_")]
seurat_data2$DF_hi.lo[which(seurat_data2@meta.data$DF_hi.lo == "Doublet" & seurat_data2@meta.data[,paste("DF.classifications_0.25",mpK,nExp_poi.adj,sep="_")] == "Singlet")] <- "Doublet_lo"
seurat_data2$DF_hi.lo[which(seurat_data2@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"

table(seurat_data2$seurat_cluster,seurat_data2$DF_hi.lo)
#further filter cells
#filted unqualified cells
#seurat_data3 <- subset(x = seurat_data2, subset = percent.mt < 20)
#seurat_data3 <- subset(x = seurat_data3, subset = nFeature_RNA < 5000)
#seurat_data3 <- subset(x = seurat_data3, subset = nCount_RNA < 50000)
#A1: 17977 features across 15819 samples (16068)

p2<- DimPlot(seurat_data2, group.by = "DF_hi.lo",label = F,cols = ppCor)
#p3<- DimPlot(seurat_data3, group.by = "DF_hi.lo",label = F,cols = ppCor)
#seurat_data4 <- subset(x = seurat_data3, subset = DF_hi.lo == "Singlet" )
seurat_data4 <- subset(x = seurat_data2, subset = DF_hi.lo == "Singlet" )
p4<- DimPlot(seurat_data4, group.by = "DF_hi.lo",label = F,cols = ppCor[2])
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_no_filter_Doublet_umap.pdf"), p2,width=6, height=6)
#ggsave(paste0("D:/6.gongchen/",f_name,"/remove_low_cell_Doublet_umap.pdf"), p3,width=6, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_remove_low_cell_and_Doublet_umap.pdf"), p4,width=6, height=6)
#seurat_data2<-0
DefaultAssay(seurat_data4) <- "RNA"    ## very important
#normalization
seurat_data4 <- SCTransform(seurat_data4, verbose = FALSE)

DefaultAssay(seurat_data2) <- "RNA"    ## very important
#normalization
seurat_data2 <- SCTransform(seurat_data2, verbose = FALSE)

m1<-FeaturePlot(object = seurat_data2, features = c("VIM","HLA-B","KRT7","PERP","EPCAM","PAEP","DCN","DLK1","PECAM1","LYVE1","AIF1","CD3D","CD8A","IL7R","NKG7","CD79A","MS4A2","HBA1","PF4","PPBP"),cols= c("grey", "red"),ncol=4)
m2<- FeaturePlot(object = seurat_data4, features = c("VIM","HLA-B","KRT7","PERP","EPCAM","PAEP","DCN","DLK1","PECAM1","LYVE1","AIF1","CD3D","CD8A","IL7R","NKG7","CD79A","MS4A2","HBA1","PF4","PPBP"),cols= c("grey", "red"),ncol=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_main_maker_before_filter_Doublet.png"), m1,width=20, height=20)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_main_maker_after_filter_Doublet.png"), m2,width=20, height=20)

saveRDS(seurat_data2, file = paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_cell_low_cell_remove_double_remain.rds"))
saveRDS(seurat_data4, file = paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_cell_low_cell_double_both_remove.rds"))
seurat_data2<-seurat_data4<-sweep.res.list_sub_cell<-0

}



#deeply cleandata for single cells
#remove of multiple cells
#data_anno<-readRDS(file = "D:/2.jianhai/CTRL1_out/raw_cell_double_finder_matadata_1123.rds")
#data_anno1<-as.data.frame(data_anno)
#Cell_meta_data<-merge(Cell@meta.data,data_anno1[,c("BARCODE_new","library","DF_hi.lo")],by="BARCODE_new")
#rownames(Cell_meta_data)<-Cell_meta_data$BARCODE_new
#head(Cell_meta_data);dim(Cell_meta_data)##157617     20
#Cell@meta.data<-Cell_meta_data

#table(Cell$DF_hi.lo)
#Doublet_hi    Singlet 
#     6664     150953 
# multiple rate:0.04207668
#Cell<-subset(x = Cell, subset = DF_hi.lo == "Singlet")
#Cell #150953



Cell<-readRDS(file = "/mnt/data/chenwei/huangnana/3.seurat_result/SIUGR_L1_forced_selected_before_double_finder_Cell_raw.rds")

## seurat downsteam normalization and scale 
#not linelogist for mt and rps
Cell<- SCTransform(object = Cell,verbose = FALSE)
##remove gene of mt and rps
MT_gene<-grep(pattern = "^MT-", x = rownames(Cell), value = TRUE)#13 MT genes
RPS_gene<-grep(pattern = "^RPS", x = rownames(Cell), value = TRUE)

counts <- GetAssayData(Cell, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(MT_gene,RPS_gene))),]
Cell_new <- subset(Cell, features = rownames(counts))

#Perform linear dimensional reduction
Cell_new <- RunPCA(object = Cell_new,verbose = FALSE)
DimPlot(object = Cell_new, reduction = "pca")
#Determine statistically significant principal components
ElbowPlot(object = Cell_new)
# Examine and visualize PCA results a few different ways
print(x = Cell_new[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = Cell_new, dims = 1:2, reduction = "pca")
DimHeatmap(object = Cell_new, dims = 1:10, cells = 500, balanced = TRUE)

### 特征值选取确定标准Determine statistically significant principal components
#Seurat clusters cells based on their PCA scores, #with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set. 
#Cell_new2 <- JackStraw(object = Cell_new, num.replicate = 100)
#Cell_new2 <- ScoreJackStraw(object = Cell_new2, dims = 1:20)
#JackStrawPlot(object = Cell_new2, dims = 1:22)
#Cell_new2<-0

#可视化降维，使用sctransform时可以选择更多的PC---具体的选择标准应该是什么
Cell_new <- RunUMAP(object = Cell_new, dims = 1:25,verbose = FALSE)
Cell_new<- RunTSNE(Cell_new, dims = 1:25)
plot0<-DimPlot(Cell_new, reduction = "pca")
plot1<-DimPlot(Cell_new, reduction = "umap")
plot2<-DimPlot(Cell_new, reduction = "tsne")
CombinePlots(plots = list(plot0,plot1, plot2),legend="top",ncol=3)

#细胞周期评估
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
length(s.genes);length(g2m.genes)#43； 54
Cell_new <- CellCycleScoring(Cell_new, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) #presumably running on SCT assay
##Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
##Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
Cell_new$CC.Difference <- Cell_new$S.Score-Cell_new$G2M.Score

Cell_score_plot1<-DimPlot(Cell_new,reduction = "umap",group.by= "Phase",split.by = "Phase",cols =ppCor)
Cell_score_plot2<-DimPlot(Cell_new,reduction = "pca",group.by= "Phase",split.by = "Phase",cols =ppCor)
Cell_score_plot1/Cell_score_plot2
DimPlot(Cell_new,reduction = "umap",group.by= "Phase",cols =ppCor)+ ggtitle('Cell cycle evaluation') 
DimPlot(Cell_new,reduction = "umap",group.by= "Phase",split.by= "Phase",cols =ppCor)+ ggtitle('Cell cycle evaluation') 

# Visualize the distribution of cell cycle markers across
RidgePlot(Cell_new, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
FeaturePlot(object = Cell_new, features = c("PCNA", "TOP2A", "MCM6", "MKI67"),cols= c("grey", "red"),ncol=2)
table(Cell_new$Phase)
#  G1   G2M     S 
#  5970 1362 2097  
#高变基因中的cell cycle 基因
library(VennDiagram)
VRF_genes<-VariableFeatures(Cell_new)
length(VRF_genes)#2981
cell_cycle_genes<-c(s.genes,g2m.genes)

grid.newpage(); #清空画板，开始画新图
grid.draw(venn.diagram(list(VRF_genes=VRF_genes,cell_cycle_genes=cell_cycle_genes),fill=c(ppCor[2],ppCor[1]), filename = NULL))
grid.newpage(); #清空画板，开始画新图


#关注状态直观绘制
DimPlot(object = Cell_new, reduction = "umap",label = TRUE) + NoLegend() + ggtitle('Cell_sctransform') 
##sample information plot in generation
#plots <- DimPlot(Cell_new, group.by = c("Family","Tissue","cell.type","Age_group","gender","Age_gender_group"), combine = FALSE,cols = ppCor)
#plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(ncol =1, override.aes = list(size = 5))))
#CombinePlots(plots)

# Cluster the cells
#Seurat采用的是graph-based聚类方法，k-means方法在V3中已经不存在了
Cell_new <- FindNeighbors(object = Cell_new, dims = 1:25,verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new, verbose = FALSE) #default=0.8
#Cell <- FindClusters(object = Cell, resolution = 1.5,verbose = FALSE

##different resolution
#Cell_new <- FindClusters( object = Cell_new,resolution = c(seq(.4,0.6,0.8,1.2,1.5,2,2.5,3,5)))

Cell_new <- FindClusters(object = Cell_new, resolution = 5,verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new, resolution = 3,verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new, resolution = 2.5,verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new, resolution = 2,verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new, resolution = 1.5,verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new, resolution = 1.2,verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new, resolution = 1,verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new, verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new,resolution = 0.6, verbose = FALSE)
Cell_new <- FindClusters(object = Cell_new, resolution = 0.4,verbose = FALSE)


plot1<-DimPlot(Cell_new, group.by ="SCT_snn_res.0.8",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.0.8') 
plot2<-DimPlot(Cell_new, group.by ="SCT_snn_res.1",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.1') 
plot3<-DimPlot(Cell_new, group.by ="SCT_snn_res.1.2",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.1.2') 
plot4<-DimPlot(Cell_new, group.by ="SCT_snn_res.1.5",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.1.5') 
plot5<-DimPlot(Cell_new, group.by ="SCT_snn_res.2",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.2') 
plot6<-DimPlot(Cell_new, group.by ="SCT_snn_res.2.5",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.2.5') 
plot7<-DimPlot(Cell_new, group.by ="SCT_snn_res.3",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.3') 
plot8<-DimPlot(Cell_new, group.by ="SCT_snn_res.5",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.5') 
CombinePlots(plots = list(plot1, plot2,plot3,plot4,plot5, plot6,plot7,plot8),ncol = 4,legend = NULL)
plot9<-DimPlot(Cell_new, group.by ="SCT_snn_res.0.4",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.0.4') 
plot10<-DimPlot(Cell_new, group.by ="SCT_snn_res.0.6",label = TRUE)+ NoLegend()+ ggtitle('SCT_snn_res.0.6') 
CombinePlots(plots = list(plot9,plot10),ncol = 2,legend = NULL)

length(unique(Cell_new@meta.data$SCT_snn_res.0.4));length(unique(Cell_new@meta.data$SCT_snn_res.0.6));#24 33
length(unique(Cell_new@meta.data$SCT_snn_res.0.8));length(unique(Cell_new@meta.data$SCT_snn_res.1));#40 43
length(unique(Cell_new@meta.data$SCT_snn_res.1.2));length(unique(Cell_new@meta.data$SCT_snn_res.1.5)); #48 52
length(unique(Cell_new@meta.data$SCT_snn_res.2));length(unique(Cell_new@meta.data$SCT_snn_res.2.5)); #62 70
length(unique(Cell_new@meta.data$SCT_snn_res.3));length(unique(Cell_new@meta.data$SCT_snn_res.5)) #76 107


#分群评估
#参考：https://www.jianshu.com/p/99cb6dc8de45  https://www.jianshu.com/p/a7a6b8b11e3c
library(clustree)
clustree(Cell_new@meta.data, prefix = "SCT_snn_res.")

library(ggalluvial)
library(tidyverse)
##绘图时间长
head(Cell_new@meta.data)
ggplot(data = Cell_new@meta.data,
       aes(axis1 = SCT_snn_res.0.4, axis2 = SCT_snn_res.0.6,axis3 = SCT_snn_res.0.8,axis4 = SCT_snn_res.1,axis5 = SCT_snn_res.1.2,axis6 = SCT_snn_res.1.5,axis7 = SCT_snn_res.2)) +
  scale_x_discrete(limits = c(paste0("SCT_snn_res.",seq(.4,1.5,.2))), expand = c(.01, .05)) +
  geom_alluvium(aes(fill = SCT_snn_res.1.5)) +
  geom_stratum() + geom_text(stat = "stratum", infer.label = TRUE) +
  #coord_polar()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("cell number in each cluster")

saveRDS(Cell_new, file = "/mnt/data/chenwei/huangnana/3.seurat_result/SIUGR_L1_object.rds")
DimPlot(object = Cell_new, reduction = "umap",label = TRUE) + NoLegend() + ggtitle('Cell_default0.8_sctransform') 

#合并注释大类群
#saveRDS(Cell_new, file = "D:/2.jianhai/CTRL1_out/Cell_1201.rds")
target<-readRDS(file = "/mnt/data/chenwei/huangnana/3.seurat_result/SIUGR_L1_object.rds")
#target<-Cell_new
#For major cell population
major_plot<-FeaturePlot(object = target, features = c("VIM","HLA-B","KRT7","PERP","PAEP","DCN","DLK1","PECAM1","LYVE1","AIF1","CD3D","NKG7","CD79A","MS4A2","CD34","HBA1","PPBP"),cols= c("grey", "red"),ncol=4)
ggsave(file=paste0("/mnt/data/chenwei/huangnana/3.seurat_result/",sample_name,"_umap_plot_for_major_maker.png"),major_plot,width = 6*4, height = 6*4)

#For no_immunine Cells 
no_immunine_plot<-FeaturePlot(object = target, features = c("VIM","CDH1","EGFR","HLA-G","MMP11","CGA","CYP19A1","CSH2","ERVFRD-1",
                                                            "DCN","DLK1","THY1","COL1A1","LAMA2","TIMP1",
                                                            "MYH11","RGS5","PRL","IGFBP1","APOD","COL6A2",
                                                            "PECAM1","LYVE1","EGFL7","PPBP","PF4","EPCAM","PAEP"),cols= c("grey", "red"),ncol=5)
ggsave(file=paste0("/mnt/data/chenwei/huangnana/3.seurat_result/",sample_name,"_umap_plot_for_no_immunine_maker.png"),no_immunine_plot,width = 6*5, height = 6*6)

#For immunine Cells 
immunine_plot<-FeaturePlot(object = target,features = c("MS4A1","CD79A","CD79B","JCHAIN","TCL1A","FCER2", 
                                                        "CD3E", "CD3D", "CD8A","IL7R","SELL","CCR7","CCL5","S100A4","GZMK", 
                                                        "NKG7", "GNLY","KLRB1","FGFBP2","XCL1",
                                                        "AIF1","LYZ","CST3","CD14","FCGR3A","MS4A7","S100A8","FCN1","CLEC9A","CD1C",
                                                        "CSF1R","CD163","CD209","CD69","LYVE1","LYZ","APOE","MS4A3",
                                                        "HBB","KIT"), pt.size = 0.2,cols= c("grey", "red"),ncol=6)
ggsave(file=paste0("/mnt/data/chenwei/huangnana/3.seurat_result/",sample_name,"_umap_plot_for_immunine_maker.png"),immunine_plot,width = 6*6, height = 6*6)

#used in my placental article
immunine_plot2<-FeaturePlot(object = target, features = c("AIF1","CD14","S100A8","FCN1","FCGR3A","CLEC9A","CD1C",
                                                          "CSF1R","CD163","CD209","CD69","LYVE1","LYZ","APOE","MS4A3",
                                                          "MS4A1","CD79A","CD79B","JCHAIN",
                                                          "CD3D","CD8A","IL7R","NKG7","KLRB1","FGFBP2","S100A4","HBB","KIT"),cols= c("grey", "red"),ncol=6)
ggsave(file=paste0("/mnt/data/chenwei/huangnana/3.seurat_result/",sample_name,"_umap_plot_for_immunine_maker2.png"),immunine_plot2,width = 6*6, height = 6*5)

#for TNK cells
FeaturePlot(object = target, features = c("CD3D","CD4","CD8A","KLRB1","NKG7","NCAM1","FGFBP2","FCGR3A","CXCR4","XCL1","PRF1","GNLY","GZMA","GZMB","GZMK"),cols= c("grey", "red"),ncol= 5)
FeaturePlot(object = target, features = c("MS4A2","PTPRC","CD69","CD3D","CD4","CD8A","IL7R","CCR7","SELL","FOXP3", "KLRF1",
                                          "FGFBP2","NKG7","GNLY","GZMB","GZMK","GZMA","NCAM1","TRDC","TRGC2"),cols= c("grey", "red"),ncol= 5)

FeaturePlot(object = target, features = c("PTPRC","HBA1","MS4A2","CD34","CD69","CD3D","CD8A","CD4","KLRB1","KLRF1","IL7R","CCR7","SELL","FOXP3",
                                          "NKG7","NCAM1", "FGFBP2","FCGR3A","GNLY","GZMB","GZMK","TRDV2","TRGV9","MKI67","TOP2A"),cols= c("grey", "red"),ncol= 5)

#immune maker 
FeaturePlot(object =  target, features = c("PTPRC"),cols= c("grey", "red"))
#tissue resident markers 
FeaturePlot(object =  target, features = c("CD69","ITGA1","CD9"),cols= c("grey", "red"),ncol=3)
#Mitotic or  proliferating subpopulation: MKI67, TOP2A ,TK1
FeaturePlot(object = target, features = c("CD3D","CD4","CD8A","KLRB1","KLRF1","MKI67","TOP2A","TK1"),cols= c("grey", "red"),ncol=4)
#T cell maker
FeaturePlot(object = target, features = c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B"),cols= c("grey", "red"),ncol=4)
FeaturePlot(object = target, features = c("TRAC","TRBC1","TRBC2","TRDC","TRGC1","TRGC2"),cols= c("grey", "red"),ncol=3)
FeaturePlot(object = target, features = c("TRAV1-2","TRDV1","TRDV2","TRDV3","TRGV3","TRGV4","TRGV9"),cols= c("grey", "red"),ncol=3)

#immature erythrocytic cells:CD34,GATA2
FeaturePlot(object = target,features = c("CD34","GATA2", "GATA1","TFRC","CD47","HBA1"),cols= c("grey", "red"))
#Mast_Cells :MS4A2
FeaturePlot(object = target, features = c("KIT","MS4A2","CPA3","IL7R"),cols= c("grey", "red"),ncol=2)

#Natural killer Cells:  KLRB1 也存在一定的极化现象
FeaturePlot(object = target, features = c("KLRD1","KLRB1","KLRF1"),cols= c("grey", "red"),ncol= 3)
FeaturePlot(object = target, features = c("CD3D","CD4","CD8A","KLRB1","NKG7","NCAM1","FGFBP2","FCGR3A","CXCR4","XCL1","PRF1","GNLY","GZMA","GZMB","GZMK",
                                          "ITGAE","TBX21","CX3CR1","IL7R","SELL"),cols= c("grey", "red"),ncol= 5)

#circulating NK (cNK) cell markers FCGR3A, CX3CR1, and the transcription factor T-bet (TBX21)
FeaturePlot(object = target, features = c("FCGR3A","CX3CR1","TBX21","IFNG","CCR1","CCR7","XCL2","ENTPD1","CYP26A1","B4GALNT1","ANXA1", "ITGB2","ITGAE"),cols= c("grey", "red"),ncol=5)

#ILC1-like CELL or CD56bright NK cells : IL7R+TCF7+SELL+KLRD1+NCAM1+ :  innate lymphocyte cell marker CD127 (also known as IL7R)
#unclassical circulating NK (ucNK)   lacked FCGR3A and FGFBP2 
FeaturePlot(object = target, features = c("IL7R","TCF7","SELL","KLRD1","NCAM1","CD160","ITGAE","TBX21"),cols= c("grey", "red"),ncol=4)

#blood γδT Cells :CD3D nagetive; KLRB1+NCAM1, CX3CR1+ TRDV2 and TRGV9 pos
FeaturePlot(object =  target, features = c("KLRG1","KLRD1","KLRB1","KLRF1","TRDV2","TRGV9"),cols= c("grey", "red"),ncol=2)

##Naïve CD4 CD3D+CD4+CCR7 SELL  naive marker genes: SELL, TCF7,CCR7 and ID3  #Memory CD4+ T: IL7R, S100A4, CCR10
FeaturePlot(object = target, features = c("CD3D","CD4","CD8A","IL7R","CCR7","ID3","SELL","TCF7","S100A4", "CCR10"),cols= c("grey", "red"),ncol= 4)

#Active CD4T: CD3D+CD4+LTBhiGATA3+; cTh1: CD3D+CD4+CXCR3+TBX21+STAT4+IFNG+   ex-cTh17: CD3D+CD4+RORA+IL23R+CCR6+STAT4+EOMES+IFNG+
FeaturePlot(object = target, features = c("LTB","GATA3","CXCR3","TBX21","STAT1","STAT4","IFNG","RORA","RORC","CCR6","IL23R","EOMES"),cols= c("grey", "red"),ncol= 3)

### Regulatory T cells(Treg): FOXP3, CTLA-4, TNFRSF18; IL2RA (CD25 +)
FeaturePlot(object = target, features = c("FOXP3","TNFRSF18","IL2RA","CTLA4","TNFRSF9","TIGIT","IL10","IL23R","IL17A","RORA","CCR6"),cols= c("grey", "red"),ncol=3)

FeaturePlot(object = target, features = c("CD3D","CD3E","CD8A","CD8B","CD4","FOXP3"),cols= c("grey", "red"),ncol=3)
FeaturePlot(object = target, features = c("FOXP3"),cols= c("grey", "red"))

#for cell type sleceted
target$selected_cluster <- rep("other_cells",nrow(target@meta.data))
DimPlot(object = target, group.by ="SCT_snn_res.1.2", label = TRUE) + NoLegend()
Idents(object = target) <- "SCT_snn_res.1.2"
select_name0 <- "CD4_Ts_other"
sub_cell0<-subset(x = target,idents=c("3","5"))
target$selected_cluster[Cells(sub_cell0)] <- select_name0

# Regulatory T cells(Treg): FOXP3, CTLA-4, TNFRSF18; IL2RA (CD25 +)
select_name1 <- "Tregs_FOXP3_pos"
sub_cell1<-subset(x = sub_cell0,FOXP3>0) 
sub_cell0$selected_cluster[Cells(sub_cell1)] <- select_name1
target$selected_cluster[Cells(sub_cell1)] <- select_name1

prop.table(table(sub_cell0$selected_cluster))
#  CD4_Ts    Tregs_FOXP3_pos 
#0.95865633      0.04134367  
prop.table(table(target$selected_cluster))
#CD4_Ts_other     other_cells   Tregs_FOXP3_pos 
#0.15087434      0.84261895      0.00650671 

metadata_for_target<-target@meta.data
target_coordinate<-data.frame(Embeddings(target[["umap"]]))
target_coordinate$BARCODE_new<-rownames(target_coordinate)
target_coordinate2 <- merge(metadata_for_target,target_coordinate)
target_coordinate2$selected_cluster <-factor(target_coordinate2$selected_cluster,levels=c(select_name1,select_name0,"other_cells"))
head(target_coordinate2)

selected_cluster_umap<-ggplot(data=target_coordinate2[order(target_coordinate2$selected_cluster,decreasing = T),], 
                              mapping=aes(x=UMAP_1,y=UMAP_2,colour = selected_cluster))+
  geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
  labs(title ="Distribution of Tregs")+scale_color_manual(values=c(ppCor_all[3],ppCor_all[5],"grey"))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                   axis.text.y = element_text(size=15),legend.position = "right")
selected_cluster_umap#+facet_grid(.~Age_group,margins=TRUE)

# Regulatory T cells(Treg): FOXP3, CTLA-4, TNFRSF18; IL2RA (CD25 +)
target$selected_cluster[Cells(sub_cell0)] <- select_name0
select_name1 <- "Tregs_CD25_pos"
sub_cell1<-subset(x = sub_cell0,IL2RA>0) 
sub_cell0$selected_cluster[Cells(sub_cell1)] <- select_name1
target$selected_cluster[Cells(sub_cell1)] <- select_name1

prop.table(table(sub_cell0$selected_cluster))
# CD4_Ts_other Tregs_CD25_pos 
# 0.91989664     0.08010336  
prop.table(table(target$selected_cluster))
#CD4_Ts_other    other_cells Tregs_CD25_pos 
# 0.14477430     0.84261895     0.01260675  

metadata_for_target<-target@meta.data
target_coordinate<-data.frame(Embeddings(target[["umap"]]))
target_coordinate$BARCODE_new<-rownames(target_coordinate)
target_coordinate2 <- merge(metadata_for_target,target_coordinate)
target_coordinate2$selected_cluster <-factor(target_coordinate2$selected_cluster,levels=c(select_name1,select_name0,"other_cells"))
head(target_coordinate2)

selected_cluster_umap<-ggplot(data=target_coordinate2[order(target_coordinate2$selected_cluster,decreasing = T),], 
                              mapping=aes(x=UMAP_1,y=UMAP_2,colour = selected_cluster))+
  geom_point(stat= "identity",size=1,alpha=0.5,show.legend = TRUE)+
  labs(title ="Distribution of Tregs")+scale_color_manual(values=c(ppCor_all[3],ppCor_all[5],"grey"))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
                   axis.text.y = element_text(size=15),legend.position = "right")
selected_cluster_umap#+facet_grid(.~Age_group,margins=TRUE)


#identify cell subtypes
DimPlot(object = target, group.by ="SCT_snn_res.0.8", label = TRUE) + NoLegend()# + ggtitle('Cell_default0.6_sctransform') 
Idents(object = target) <- "SCT_snn_res.0.6"
Idents(target)
target$merge_cluster <- as.character(Idents(target))

#B
sub_cell<-subset(x = target,idents=c("10","31"))
table(sub_cell@meta.data$Tissue);table(sub_cell@meta.data$cell.type)
target$merge_cluster[Cells(sub_cell)] <- "B_Cell_CD79A"
DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

#T_NK_raw_Erythrocyte_Mast_Cell
sub_cell<-subset(x = target,idents=c("0","2","3","4","5","6","9","11","18","20","27","32"))
table(sub_cell@meta.data$Tissue)
target$merge_cluster[Cells(sub_cell)] <- "T_NK_raw_Erythrocyte_Mast_Cell"
DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

#Erythrocyte
sub_cell<-subset(x = target,idents=c("29"))
table(sub_cell@meta.data$Tissue)
target$merge_cluster[Cells(sub_cell)] <- "Erythrocyte_HBA1"
DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

#Megakaryocytes
sub_cell<-subset(x = target,idents=c("30"))
table(sub_cell@meta.data$Tissue)
target$merge_cluster[Cells(sub_cell)] <- "Megakaryocytes_PPBP"
DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

# Myeloid cells 
sub_cell<-subset(x = target,idents=c("1","7","12","13","15","17","22","26"))
table(sub_cell@meta.data$Tissue)
target$merge_cluster[Cells(sub_cell)] <- "Myeloid_Cell_AIF1"
DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

# Trophoblast_mix_EECs  
sub_cell<-subset(x = target,idents=c("14","23"))
table(sub_cell@meta.data$Tissue);table(sub_cell@meta.data$cell.type)
sub_cell@meta.data[which(sub_cell@meta.data$Tissue=="PBMCs"),]#待议
sub_cell@meta.data[which(sub_cell@meta.data$cell.type=="maternal"),]#待议
target$merge_cluster[Cells(sub_cell)] <- "Trophoblast_mix_EECs_KRT7"
DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 
##Base PBMCs  Peri  Vill 
##2999     1  2198   769 
##fetal maternal 
##5837      130 

#Epithelial_Cell_EPCAM_PAEP
sub_cell2<-subset(x = sub_cell, subset = cell.type  =="fetal",invert = TRUE)
table(sub_cell2$Tissue) #存在一个PBMCs
#Base PBMCs  Peri  Vill 
#39     1    90     0 
target$merge_cluster[Cells(sub_cell2)] <- "Epithelial_Cell_EPCAM_PAEP"
DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

#Decidual Stromal cells_raw
sub_cell<-subset(x = target,idents=c("8","19","21"))
table(sub_cell@meta.data$Tissue)
sub_cell@meta.data[which(sub_cell@meta.data$Tissue=="PBMCs"),]#待议 14 PBMCs
target$merge_cluster[Cells(sub_cell)] <- "Decidual_stromal_cells_DKK1_ACTA2"

sub_cell2<-subset(x = sub_cell,subset = cell.type == "fetal")
target$merge_cluster[Cells(sub_cell2)] <- "Fetal_Stromal_cells_Endothelial_cells_DLK1"

DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

#Fetal Stromal cells
sub_cell<-subset(x = target,idents=c("16","25"))
table(sub_cell@meta.data$Tissue)
target$merge_cluster[Cells(sub_cell)] <- "Fetal_Stromal_cells_Endothelial_cells_DLK1"

sub_cell2<-subset(x = sub_cell,subset = cell.type == "maternal")
target$merge_cluster[Cells(sub_cell2)] <- "Decidual_stromal_cells_DKK1_ACTA2"

DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

#Endothelial cells
## vascular endothelial Cells(F&D VECs) &  lymphatic endothelial Cells (LECs)
sub_cell<-subset(x = target,idents=c("24"))
table(sub_cell@meta.data$Tissue)
target$merge_cluster[Cells(sub_cell)] <- "Lymphatic_Endothelial_Cells_PECAM1_LYVE1_pos" #1 Vill  待议拟归于VECs
DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

#sub_cell2<-subset(x = sub_cell,subset = Tissue == "Vill")
#sub_cell$merge_cluster[Cells(sub_cell2)] <- "Vill_Lymphatic_Endothelial_Cells_PECAM1_LYVE1_pos" 
#table(sub_cell@meta.data$merge_cluster)
#DimPlot(object = sub_cell,group.by = "merge_cluster",cols = ppCor) + NoLegend() 

sub_cell<-subset(x = target,idents=c("28"))
table(sub_cell@meta.data$Tissue)
target$merge_cluster[Cells(sub_cell)] <- "Vascular_Endothelial_Cells_PECAM1_LYVE1_neg" 
DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

#plot merge group
DimPlot(target, group.by = "merge_cluster",label = TRUE,cols = ppCor)
DimPlot(object = target, reduction = "umap", split.by = "merge_cluster",group.by = "merge_cluster",cols = ppCor,ncol = 4)
table(target$merge_cluster,target$cell.type)
#               
#                                             fetal maternal
#B_Cell_CD79A                                   500     4907 
#Decidual_stromal_cells_DKK1_ACTA2                0    11454
#Epithelial_Cell_EPCAM_PAEP                       0      130
#Erythrocyte_HBA1                               301       16
#Fetal_Stromal_cells_Endothelial_cells_DLK1    5316        0
#Lymphatic_Endothelial_Cells_PECAM1_LYVE1_pos     1     1823 #待议
#Megakaryocytes_PPBP                             20      194 
#Myeloid_Cell_AIF1                             9521    28430 
#T_NK_raw_Erythrocyte_Mast_Cell                5983    74561
#Trophoblast_mix_EECs_KRT7                     5837      0
#Vascular_Endothelial_Cells_PECAM1_LYVE1_neg    131      589

#subset for each Cells groups
table(target$merge_cluster)
#B_Cell_CD79A            Decidual_stromal_cells_DKK1_ACTA2                   Epithelial_Cell_EPCAM_PAEP 
#  5407                                        11454                                          130 
#Erythrocyte_HBA1   Fetal_Stromal_cells_Endothelial_cells_DLK1 Lymphatic_Endothelial_Cells_PECAM1_LYVE1_pos 
#  317                                         5316                                         1824 
#Megakaryocytes_PPBP                            Myeloid_Cell_AIF1               T_NK_raw_Erythrocyte_Mast_Cell 
#  214                                        37951                                        80544 
#Trophoblast_mix_EECs_KRT7  Vascular_Endothelial_Cells_PECAM1_LYVE1_neg 
#  5837                                          720 

saveRDS(target, file = "D:/2.jianhai/CTRL1_out/Cell_1201.rds")

target<-readRDS(file = "D:/2.jianhai/CTRL1_out/Cell_1201.rds")

#1.For Trophoblast_mix_EECs_KRT7
sub_cell<-subset(x = target, subset = merge_cluster == "Trophoblast_mix_EECs_KRT7")
#sub_cell2<-subset(x = sub_cell, subset = cell.type  =="fetal")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Trophoblast_mix_EECs_KRT7 <- sub_cell
table(Trophoblast_mix_EECs_KRT7$Tissue)
#Base PBMCs  Peri  Vill 
#2999     1  2198   769 

#2.For Epithelial_Cell_EPCAM_PAEP
#sub_cell2<-subset(x = sub_cell, subset = cell.type  =="fetal",invert = TRUE)
sub_cell<-subset(x = target, subset = merge_cluster == "Epithelial_Cell_EPCAM_PAEP")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Epithelial_Cell_EPCAM_PAEP <- sub_cell
table(Epithelial_Cell_EPCAM_PAEP$Tissue) #存在一个PBMCs
# Base PBMCs  Peri  Vill 
#  39     1    90     0 

#3.For B_Cell_CD79A
sub_cell<-subset(x = target, subset = merge_cluster == "B_Cell_CD79A")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
B_Cell_CD79A <- sub_cell
table(B_Cell_CD79A$Tissue)
#Base PBMCs  Peri  Vill 
#507  4514   169   217 

#4. For Myeloid_Cell_AIF1
sub_cell<-subset(x = target, subset = merge_cluster == "Myeloid_Cell_AIF1")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Myeloid_Cell_AIF1 <- sub_cell

#5. For T_NK_raw_Erythrocyte_Mast_Cell
sub_cell<-subset(x = target, subset = merge_cluster == "T_NK_raw_Erythrocyte_Mast_Cell")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
T_NK_raw_Erythrocyte_Mast_Cell <- sub_cell

#6. For Megakaryocytes_PPBP
sub_cell<-subset(x = target, subset = merge_cluster == "Megakaryocytes_PPBP")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Megakaryocytes_PPBP <- sub_cell

#7. For Decidual_stromal_cells_DKK1_ACTA2 
sub_cell<-subset(x = target, subset = merge_cluster == "Decidual_stromal_cells_DKK1_ACTA2")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Decidual_stromal_cells_DKK1_ACTA2 <- sub_cell

#8. For Fetal_Stromal_cells_Endothelial_cells_DLK1
sub_cell<-subset(x = target, subset = merge_cluster == "Fetal_Stromal_cells_Endothelial_cells_DLK1")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Fetal_Stromal_cells_Endothelial_cells_DLK1 <- sub_cell

#9. For Endothelial_Cell_PECAM1
sub_cell<-subset(x = target, subset = merge_cluster %in% c("Lymphatic_Endothelial_Cells_PECAM1_LYVE1_pos","Vascular_Endothelial_Cells_PECAM1_LYVE1_neg"))
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Endothelial_Cell_PECAM1 <- sub_cell

Cell_submian_list<-c(list(Trophoblast_mix_EECs_KRT7),list(Epithelial_Cell_EPCAM_PAEP),list(Fetal_Stromal_cells_Endothelial_cells_DLK1),list(Decidual_stromal_cells_DKK1_ACTA2),
                     list(Megakaryocytes_PPBP),list(T_NK_raw_Erythrocyte_Mast_Cell),list(Myeloid_Cell_AIF1),list(B_Cell_CD79A),list(Endothelial_Cell_PECAM1))
names(Cell_submian_list)<-c("Trophoblast_mix_EECs_KRT7","Epithelial_Cell_EPCAM_PAEP","Fetal_Stromal_cells_Endothelial_cells_DLK1","Decidual_stromal_cells_DKK1_ACTA2",
                            "Megakaryocytes_PPBP","T_NK_raw_Erythrocyte_Mast_Cell","Myeloid_Cell_AIF1","B_Cell_CD79A","Endothelial_Cell_PECAM1" )
saveRDS(Cell_submian_list, file = "D:/2.jianhai/CTRL1_out/Cell_submian_list_1201.rds")


## 查看每一类细胞的数目
sort(table(x = Idents(object = target)),decreasing = TRUE)
#查看某一cluster的具体细胞barcode
head(WhichCells(target,idents="2"))
# 提取某一cluster细胞查看。
##
head(subset(as.data.frame(target@active.ident),target@active.ident=="2"))
## or
head(subset(as.data.frame(target@meta.data),target@meta.data$merge_cluster=="Trophoblast_mix_EECs_KRT7"))
