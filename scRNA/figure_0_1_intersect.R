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

p2<- DimPlot(seurat_data2, group.by = "DF_hi.lo",label = F,cols = ppCor)
seurat_data4 <- subset(x = seurat_data2, subset = DF_hi.lo == "Singlet" )
p4<- DimPlot(seurat_data4, group.by = "DF_hi.lo",label = F,cols = ppCor[2])
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_no_filter_Doublet_umap.pdf"), p2,width=6, height=6)
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

###############
###prepare lists for intergration
###############
###one:for double_removed
mydata_list<-c("N1_out","N2_out","N3_out","N4_out","N5_out","A1_out","A3_out","A4_out","A8_out")
#data_name<-c()
seurat_data_list<-list()
for ( sample_name in mydata_list){
  sample<- unlist(lapply(strsplit(sample_name,"_"), function(x) x[1]))
  print(as.character(sample))
  seurat_data<-readRDS(file = paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_cell_low_cell_double_both_remove.rds"))
  seurat_data_list<-c(seurat_data_list,list(seurat_data))
}
data_name<-c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","Arrest_1","Arrest_2","Arrest_3","Arrest_4")
names(seurat_data_list)<-data_name
saveRDS(seurat_data_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713.rds")

dtlist <- seurat_data_list
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 3000) #2976 
dtlist <- PrepSCTIntegration(object.list = dtlist,anchor.features = intfts)
anchors <- FindIntegrationAnchors(object.list = dtlist, normalization.method = "SCT", anchor.features = intfts)
saveRDS(anchors, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713_anchor.rds")
seurat_data_list<-anchors<-dtlist<-intfts<-0

########################
###perform Integrate
########################
anchors<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713_anchor.rds")
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated <- RunPCA(integrated)
ElbowPlot(object = integrated)
VizDimLoadings(object = integrated, dims = 1:2, reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:25)
integrated<- RunTSNE(integrated, dims = 1:25)
saveRDS(integrated, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713_integrated_seurat.rds")

##Two:for double remained
mydata_list<-c("N1_out","N2_out","N3_out","N4_out","N5_out","A1_out","A3_out","A4_out","A8_out")
#data_name<-c()
seurat_data_list<-list()
for ( sample_name in mydata_list){
  sample<- unlist(lapply(strsplit(sample_name,"_"), function(x) x[1]))
  print(as.character(sample))
  seurat_data<-readRDS(file = paste0("/mnt/data/chenwei/gongchen/3.seurat_result/",sample_name,"/",sample,"_forced_cell_cell_low_cell_remove_double_remain.rds"))
  seurat_data_list<-c(seurat_data_list,list(seurat_data))
}
data_name<-c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","Arrest_1","Arrest_2","Arrest_3","Arrest_4")
names(seurat_data_list)<-data_name
saveRDS(seurat_data_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_remove_double_remain_sample_list_0713.rds")
dtlist <- seurat_data_list
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 5000) #2976 
dtlist <- PrepSCTIntegration(object.list = dtlist,anchor.features = intfts)
anchors <- FindIntegrationAnchors(object.list = dtlist, normalization.method = "SCT", anchor.features = intfts)
saveRDS(anchors, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_remove_double_remain_sample_list_0713_anchor.rds")
seurat_data_list<-anchors<-dtlist<-intfts<-0
########################
### Integrate
########################
anchors<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_remove_double_remain_sample_list_0713_anchor.rds")
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated <- RunPCA(integrated)
ElbowPlot(object = integrated)
VizDimLoadings(object = integrated, dims = 1:2, reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:25)
integrated<- RunTSNE(integrated, dims = 1:25)
saveRDS(integrated, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_remove_double_remain_sample_list_0713_integrated_seurat.rds")
doublet_remain_object<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_remove_double_remain_sample_list_0713_integrated_seurat.rds")
doublet_remain_meta<-doublet_remain_object@meta.data
doublet_remain_meta$gemgroup<-as.numeric(unlist(lapply(strsplit(rownames(doublet_remain_meta),"_"), function(x) x[2])))

table(doublet_remain_meta$DF_hi.lo)
#Doublet_hi    Singlet 
#    2779       49988 
table(doublet_remain_meta$gemgroup,doublet_remain_meta$DF_hi.lo)
#   Doublet_hi Singlet
#1       1056   11046
#2        318    6524
#3        159    4544
#4        158    4354
#5        107    3872
#6         88    3092
#7        643    8751
#8        151    4215
#9         99    3590
