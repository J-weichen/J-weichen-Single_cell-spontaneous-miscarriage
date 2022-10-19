#Monocle2 Pseudotime analysis
#ref:https://nbisweden.github.io/workshop-scRNAseq/oldlabs/monocle_analysis.html
#ref:https://cloud.tencent.com/developer/article/1692225
#ref: https://www.sohu.com/a/327133920_278730
#ref：https://www.jianshu.com/p/5d6fd4561bc0
#ref:https://zhuanlan.zhihu.com/p/378365295
#integrated abject
#ref；https://github.com/satijalab/seurat/issues/1658
#https://github.com/cole-trapnell-lab/monocle3/issues/148
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

library("ggbeeswarm")
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)


#step1 ：read Seurat object and load initial coldata and matrix
target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
target$final_major_subgroup_brief_Treat <- paste(target$final_major_subgroup_brief, target$Treat, sep = "_")
table(target$final_major_subgroup_brief)

cell_name<-"MyCs"
selected_cell0<-subset(x = target, subset = final_major_subgroup_brief == cell_name)
table(selected_cell0$sample_code)
#CTRL_1   CTRL_2   CTRL_3   CTRL_4   CTRL_5 Arrest_1 Arrest_2 Arrest_3 Arrest_4 
# 16        1       16      207      159      907      366     1874      867 
#CTRL_2由于细胞太少不进行分析
mydata_list<-c("CTRL_4","CTRL_5","Arrest_1","Arrest_2","Arrest_3","Arrest_4")
#data_name<-c()
seurat_data_list<-list()
for ( sample_name in mydata_list){
  # sample_name<-"CTRL_1"#test line
  print(as.character(sample_name))
  selected_cell<-subset(x = selected_cell0, subset = sample_code == sample_name)
  DefaultAssay(selected_cell) <- "RNA"    ## very important
  selected_cell <- SCTransform(object = selected_cell,verbose = FALSE) 
  seurat_data_list<-c(seurat_data_list,list(selected_cell))
}
names(seurat_data_list)<-mydata_list
saveRDS(seurat_data_list, file = "/mnt/data/chenwei/gongchen/6.monocle2/All_MYCs_Cell_seurat_sample_list.rds")
seurat_data_list<-readRDS(file = "/mnt/data/chenwei/gongchen/6.monocle2/All_MYCs_Cell_seurat_sample_list.rds")

#dtlist <- seurat_data_list[c("CTRL_4","CTRL_5")]
dtlist <- seurat_data_list
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 5000) #2976 
dtlist <- PrepSCTIntegration(object.list = dtlist,anchor.features = intfts)

#k.filter <- min(sapply(dtlist, ncol)) + 1
k.filter <- min(150, min(sapply(dtlist, ncol))) #default=200 ref:https://github.com/satijalab/seurat/issues/997

anchors <- FindIntegrationAnchors(object.list = dtlist, k.filter = k.filter, normalization.method = "SCT", anchor.features = intfts)
#saveRDS(anchors, file = "/mnt/data/chenwei/gongchen/6.monocle2/Myeloid_Cell_sample_list.anchor.rds")
#seurat_data_list<-anchors<-dtlist<-intfts<-0
########################
### Integrate
########################
#anchors<-readRDS(file = "/mnt/data/chenwei/gongchen/6.monocle2/Myeloid_Cell_sample_list.anchor.rds")
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated <- RunPCA(integrated)
ElbowPlot(object = integrated)
VizDimLoadings(object = integrated, dims = 1:2, reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:25)
integrated<- RunTSNE(integrated, dims = 1:25)
selected_cell <- FindNeighbors(object = integrated, dims = 1:25,verbose = FALSE)
selected_cell <- FindClusters(object = selected_cell,resolution = 0.6, verbose = FALSE)
DimPlot(object = selected_cell, reduction = "umap", group.by = "final_major_subgroup_brief",cols =ppCor)
DimPlot(object = selected_cell, reduction = "umap", group.by = "sample_code",cols =ppCor)
DimPlot(object = selected_cell, reduction = "umap", group.by = "sample_code",split.by = "sample_code",cols =ppCor)


metadata_new<-selected_cell@meta.data
dim(metadata_new)# 4380   54
table(metadata_new$prefinal_cluster_brief)
names(table(metadata_new$prefinal_cluster_brief))
metadata_new$monocle_names0<-metadata_new$prefinal_cluster_brief
#metadata_new[which(metadata_new$final_major_subgroup_brief %in% c("CTBs_0","CTBs_1","CTBs_2","CTBs_3")),]$monocle_names0 <- "CTBs"
head(metadata_new)
selected_cell@meta.data<-metadata_new
saveRDS(selected_cell, file = "/mnt/data/chenwei/gongchen/6.monocle2/All_Myeloid_Cell_seurat_pbject.rds")

table(selected_cell$monocle_names0)
#cMon  DCs  Mac 
#136   69 4175 
table(selected_cell$final_major_subgroup_brief)

##call DEGs inter group
#Idents(object = selected_cell) <- "final_major_subgroup_brief" #'merge_cluster'
#Allcluster.markers_all <- FindAllMarkers(object = selected_cell,min.pct = 0.25)
#Allcluster.markers_all_default0.25<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=0.25,]
#Allcluster.markers_all_default0.25_q01<-Allcluster.markers_all_default0.25[Allcluster.markers_all_default0.25$p_val_adj<0.1,]
#DEGs_seurat_call<-Allcluster.markers_all_default0.25_q01

##step2 ：Construct CellDataSet object 10X的数据使用UMI count矩阵 
count_raw <- as(as.matrix(selected_cell@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = selected_cell@meta.data)
fData <- data.frame(gene_short_name = row.names(count_raw), row.names = row.names(count_raw))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(count_raw,phenoData = pd,featureData = fd,#lowerDetectionLimit = 0.5,
                              expressionFamily = uninormal())#expressionFamily = negbinomial.size())
#expressionFamily参数用于指定表达矩阵的数据类型，有几个选项可以选择：
#稀疏矩阵用negbinomial.size()，FPKM值用tobit()，logFPKM值用gaussianff()

##step3 Size factor and dispersion will be used to normalize data
HSMM<-monocle_cds
#hvg called by seurat：
length(VariableFeatures(selected_cell))#3701
hvg_genes<-VariableFeatures(selected_cell)
ordering_genes<-hvg_genes
hvg_HSMM <- setOrderingFilter(HSMM, ordering_genes)
## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
hvg_HSMM <- reduceDimension(hvg_HSMM,norm_method="none", 
                            reduction_method="DDRTree",
                            max_components=4,scaling=TRUE,verbose=TRUE,pseudo_expr=0)
## order cells change colors and theta to match your plot
hvg_HSMM <- orderCells(hvg_HSMM)
#Monocle does not know which the root state is for the tree, we need to define the root.
plot_cell_trajectory(hvg_HSMM, color_by = "prefinal_cluster_brief")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(hvg_HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)

#selected special HSMM
#DEGs_HSMM hvg_HSMM diff_HSMM disp_HSMM
HSMM<-hvg_HSMM
plot_cell_trajectory(HSMM, color_by = "monocle_names0",show_branch_points = FALSE,show_tree = TRUE,cell_size = 1) +
  scale_colour_manual(values=ppCor)+  theme(legend.position = "right")

#Step6 See what the trajectory looks like. for pre root identify
#Monocle does not know which the root state is for the tree, we need to define the root.
plot_cell_trajectory(HSMM, color_by = "prefinal_cluster_brief")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="State")+scale_colour_manual(values=rev(ppCor))
plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")

plot_cell_trajectory(HSMM, color_by = "monocle_names0",cell_size = 0.75)+scale_colour_manual(values=ppCor)+facet_wrap(~monocle_names0, ncol = 2) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief",cell_size = 0.75)+scale_colour_manual(values=ppCor)+facet_wrap(~final_major_subgroup_brief, ncol = 3) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)+facet_wrap(~Phase, ncol = 2)

plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)+facet_wrap(~sample_code, ncol = 3)
plot_cell_trajectory(HSMM, color_by="sample_code",cell_size=0.1)+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="Treat",cell_size=0.1)+scale_colour_manual(values=ppCor[c(3,5)])
plot_cell_trajectory(HSMM, color_by="Treat",cell_size=0.1)+scale_colour_manual(values=ppCor[c(3,5)])+facet_wrap(~Treat, ncol = 2)
plot_cell_trajectory(HSMM, color_by="monocle_names0",cell_size=0.1,show_branch_points = F)+scale_colour_manual(values=ppCor)+facet_wrap(~Treat, ncol = 2)
plot_cell_trajectory(HSMM, color_by="final_major_subgroup_brief",cell_size=0.1,show_branch_points = F)+scale_colour_manual(values=ppCor)+facet_wrap(~Treat, ncol = 2)

table(pData(HSMM)$final_major_subgroup_brief,pData(HSMM)$State)
# we clearly see that State1 should be the root # so now we can reorder the cells  again
plot_cell_trajectory(HSMM, color_by="State",cell_size=0.1)+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="State")+scale_colour_manual(values=rev(ppCor))+facet_wrap(~State, ncol = 3)


HSMM <- orderCells(HSMM, root_state=7)
# plot with coloring by pseudotime
plot_cell_trajectory(HSMM, color_by="Pseudotime")

#We can compare the inferred pseudotime to the known sampling timepoints.
meta_plot<-selected_cell@meta.data;head(meta_plot);dim(meta_plot)
pData_plot<-pData(HSMM);head(pData_plot);dim(pData_plot)
meta_plot2<-merge(meta_plot,pData_plot[,c("Pseudotime","State")],by=0)
dim(meta_plot2)
cell_distribution<-ggplot(meta_plot2,aes(x = Pseudotime,y = prefinal_cluster_brief, colour = prefinal_cluster_brief)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = ppCor) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("cell subtypes") +
  ggtitle("Cells ordered by monocle2 pseudotime")
cell_distribution+facet_wrap(~Treat, ncol = 3)

cell_distribution2<-ggplot(meta_plot2,aes(x = Pseudotime,y = prefinal_cluster_brief, colour = Treat)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = ppCor[c(3,5)]) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("cell subtypes") +
  ggtitle("Cells ordered by monocle2 pseudotime")
cell_distribution2

ggplot(meta_plot2,aes(x=Pseudotime,fill=Treat,alpha = 1/10))+geom_density(position="identity")+scale_fill_manual(values = ppCor[c(3,5)]) + theme_classic()

ggplot(meta_plot2,aes(x=Pseudotime,fill=Treat,alpha = 1/10))+geom_density()+scale_fill_manual(values = ppCor[c(3,5)]) + theme_classic()+facet_wrap(~prefinal_cluster_brief, ncol = 3)
ggplot(meta_plot2,aes(x=Pseudotime,colour=Treat))+geom_density()+scale_color_manual(values = ppCor[c(3,5)]) + theme_classic() +facet_wrap(~prefinal_cluster_brief, ncol = 3)


plotdf=pData(HSMM)
#monocle 山脊图绘制
ggplot(plotdf, aes(x=Pseudotime,y=prefinal_cluster_brief,fill=prefinal_cluster_brief))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
ggsave("tmp1.pdf",width = 13,height = 7,units = "cm")
#添加颜色
ggplot(plotdf, aes(x=Pseudotime,y=prefinal_cluster_brief,fill = stat(x))) +
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


## replace the exprs cds matrix value with orignal RNA counts after ordercells based on Seurat integrated matrix?
#Below I mae another cell_data_set object but this time using the RNA assay for the
#expression matrix. I also have to adjust the gene_metadata.
DefaultAssay(selected_cell) <- "RNA"
# Normalize RNA data for visualization purposes
selected_cell <- NormalizeData(selected_cell, verbose = FALSE)
#count_raw <- as(as.matrix(selected_cell@assays$SCT@data),'sparseMatrix')
count_raw <- as(as.matrix(selected_cell@assays$RNA@data),'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = selected_cell@meta.data)
fData <- data.frame(gene_short_name = row.names(count_raw), row.names = row.names(count_raw))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
SCT_monocle_cds <- newCellDataSet(count_raw,phenoData = pd,featureData = fd,#lowerDetectionLimit = 0.5,
                                  expressionFamily = uninormal())

HSMM@assayData <- SCT_monocle_cds@assayData
#cds_test now maintains the trajectory constructed using the integrated assay
#but contains expression data from the RNA assay.
proliferate_markers<-c("MKI67","TOP2A","TK1")
Myeloid_Cell_maker<-c("AIF1","CD14","HLA-DRA","S100A8","FCGR3A","CLEC9A","CD1C","CSF1R","CD163","CD68","CD209","LYVE1","CD86","FCN1","CCL2","MS4A3")

maker_gene<-Myeloid_Cell_maker
maker_plot<-plot_cell_trajectory(HSMM, markers =maker_gene, use_color_gradient = T,cell_size=0.1,show_branch_points = F)
maker_plot+scale_color_gradient2(low = "lightblue", mid = "darkred", high = "orange")+ theme(legend.position = "right")
maker_plot+geom_point(aes(color=value),size=0.5,na.rm = F) + theme(legend.position="top")+
  scale_color_gradient2(name = "RNA normalized Count",low = "lightblue",mid = "grey", high = "darkred")
#scale_color_plasma(name = "RNA normalized Count: value")

saveRDS(hvg_HSMM, file = "/mnt/data/chenwei/gongchen/6.monocle2/all_intergrated_monocle2_pseudotime_object_original_hvg.rds")
saveRDS(HSMM, file = "/mnt/data/chenwei/gongchen/6.monocle2/all_intergrated_monocle2_pseudotime_object_pre_gene_call.rds")
##以下未处理
#split analysis and plot by CTRL and Abortion
seurat_data_list<-readRDS(file = "/mnt/data/chenwei/gongchen/6.monocle2/All_Myeloid_Cell_seurat_sample_list.rds")
dtlist <- seurat_data_list[c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5")]
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 5000) #2976 
dtlist <- PrepSCTIntegration(object.list = dtlist,anchor.features = intfts)
anchors <- FindIntegrationAnchors(object.list = dtlist, normalization.method = "SCT", anchor.features = intfts)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated <- RunPCA(integrated)
ElbowPlot(object = integrated)
VizDimLoadings(object = integrated, dims = 1:2, reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:25)
integrated<- RunTSNE(integrated, dims = 1:25)
selected_cell <- FindNeighbors(object = integrated, dims = 1:25,verbose = FALSE)
selected_cell <- FindClusters(object = selected_cell,resolution = 0.6, verbose = FALSE)
DimPlot(object = selected_cell, reduction = "umap", group.by = "final_major_subgroup_brief",cols =ppCor)
DimPlot(object = selected_cell, reduction = "umap", group.by = "sample_code",cols =ppCor)
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
saveRDS(selected_cell, file = "/mnt/data/chenwei/gongchen/6.monocle2/CTRL_Myeloid_Cell_seurat_pbject.rds")

table(selected_cell$monocle_names0)
#CTBs EVTs STBs 
#8266 1671 1033
table(selected_cell$final_major_subgroup_brief)
#CTBs_0 CTBs_1 CTBs_2 CTBs_3 EVTs_1 EVTs_2 EVTs_3 STBs_1 STBs_2 
# 809   1405   5783    269    417    492    762    210    823 
##call DEGs inter group
#Idents(object = selected_cell) <- "final_major_subgroup_brief" #'merge_cluster'
#Allcluster.markers_all <- FindAllMarkers(object = selected_cell,min.pct = 0.25)
#Allcluster.markers_all_default0.25<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=0.25,]
#Allcluster.markers_all_default0.25_q01<-Allcluster.markers_all_default0.25[Allcluster.markers_all_default0.25$p_val_adj<0.1,]
#DEGs_seurat_call<-Allcluster.markers_all_default0.25_q01

##step2 ：Construct CellDataSet object 10X的数据使用UMI count矩阵 
count_raw <- as(as.matrix(selected_cell@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = selected_cell@meta.data)
fData <- data.frame(gene_short_name = row.names(count_raw), row.names = row.names(count_raw))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(count_raw,phenoData = pd,featureData = fd,#lowerDetectionLimit = 0.5,
                              expressionFamily = uninormal())#expressionFamily = negbinomial.size())
#expressionFamily参数用于指定表达矩阵的数据类型，有几个选项可以选择：
#稀疏矩阵用negbinomial.size()，FPKM值用tobit()，logFPKM值用gaussianff()

##step3 Size factor and dispersion will be used to normalize data
HSMM<-monocle_cds
#hvg called by seurat：
length(VariableFeatures(selected_cell))#3701
hvg_genes<-VariableFeatures(selected_cell)
ordering_genes<-hvg_genes
hvg_HSMM <- setOrderingFilter(HSMM, ordering_genes)
## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
hvg_HSMM <- reduceDimension(hvg_HSMM,norm_method="none", 
                            reduction_method="DDRTree",
                            max_components=4,scaling=TRUE,verbose=TRUE,pseudo_expr=0)
## order cells change colors and theta to match your plot
hvg_HSMM <- orderCells(hvg_HSMM)
#Monocle does not know which the root state is for the tree, we need to define the root.
plot_cell_trajectory(hvg_HSMM, color_by = "final_major_subgroup_brief")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(hvg_HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)

#selected special HSMM
#DEGs_HSMM hvg_HSMM diff_HSMM disp_HSMM
HSMM<-hvg_HSMM
plot_cell_trajectory(HSMM, color_by = "monocle_names0",show_branch_points = FALSE,show_tree = TRUE,cell_size = 1) +
  scale_colour_manual(values=ppCor)+  theme(legend.position = "right")

#Step6 See what the trajectory looks like. for pre root identify
#Monocle does not know which the root state is for the tree, we need to define the root.
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="State")+scale_colour_manual(values=rev(ppCor))
plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")

plot_cell_trajectory(HSMM, color_by = "monocle_names0",cell_size = 0.75)+scale_colour_manual(values=ppCor)+facet_wrap(~monocle_names0, ncol = 2) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief",cell_size = 0.75)+scale_colour_manual(values=ppCor)+facet_wrap(~final_major_subgroup_brief, ncol = 3) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)+facet_wrap(~Phase, ncol = 2)

plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)+facet_wrap(~sample_code, ncol = 3)
plot_cell_trajectory(HSMM, color_by="sample_code",cell_size=0.1)+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="monocle_names0",cell_size=0.1,show_branch_points = F)+scale_colour_manual(values=ppCor)+facet_wrap(~Treat, ncol = 2)
plot_cell_trajectory(HSMM, color_by="final_major_subgroup_brief",cell_size=0.1,show_branch_points = F)+scale_colour_manual(values=ppCor)+facet_wrap(~Treat, ncol = 2)

table(pData(HSMM)$final_major_subgroup_brief,pData(HSMM)$State)
# we clearly see that State1 should be the root # so now we can reorder the cells  again
plot_cell_trajectory(HSMM, color_by="State",cell_size=0.1)+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="State")+scale_colour_manual(values=rev(ppCor))+facet_wrap(~State, ncol = 3)


HSMM <- orderCells(HSMM, root_state=8)
# plot with coloring by pseudotime
plot_cell_trajectory(HSMM, color_by="Pseudotime")

#We can compare the inferred pseudotime to the known sampling timepoints.
meta_plot<-selected_cell@meta.data;head(meta_plot);dim(meta_plot)
pData_plot<-pData(HSMM);head(pData_plot);dim(pData_plot)
meta_plot2<-merge(meta_plot,pData_plot[,c("Pseudotime","State")],by=0)
dim(meta_plot2)
cell_distribution<-ggplot(meta_plot2,aes(x = Pseudotime,y = final_major_subgroup_brief, colour = final_major_subgroup_brief)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = ppCor) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("cell subtypes") +
  ggtitle("Cells ordered by monocle2 pseudotime")
cell_distribution

ggplot(meta_plot2,aes(x=Pseudotime,fill=Treat))+geom_density()+scale_fill_manual(values = ppCor[3])+ggtitle("CTRL Cells ordered by monocle2 pseudotime")
ggplot(meta_plot2,aes(x=Pseudotime,fill=final_major_subgroup_brief,alpha = 1/10))+geom_rug()+geom_density()+scale_fill_manual(values = ppCor)
ggplot(meta_plot2,aes(x=Pseudotime,fill=final_major_subgroup_brief,alpha = 1/10))+geom_rug()+geom_density()+scale_fill_manual(values = ppCor)+ theme_classic()+facet_wrap(~final_major_subgroup_brief, ncol = 3)

ggplot(meta_plot2,aes(x=Pseudotime,colour=final_major_subgroup_brief))+geom_density()+scale_color_manual(values = ppCor) + theme_classic() +facet_wrap(~final_major_subgroup_brief, ncol = 3)
ggplot(meta_plot2,aes(x=Pseudotime,colour=final_major_subgroup_brief,alpha = 1/10))+geom_rug()+geom_density()+scale_colour_manual(values = ppCor)

## replace the exprs cds matrix value with orignal RNA counts after ordercells based on Seurat integrated matrix?
#Below I mae another cell_data_set object but this time using the RNA assay for the
#expression matrix. I also have to adjust the gene_metadata.
DefaultAssay(selected_cell) <- "RNA"
# Normalize RNA data for visualization purposes
selected_cell <- NormalizeData(selected_cell, verbose = FALSE)
#count_raw <- as(as.matrix(selected_cell@assays$SCT@data),'sparseMatrix')
count_raw <- as(as.matrix(selected_cell@assays$RNA@data),'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = selected_cell@meta.data)
fData <- data.frame(gene_short_name = row.names(count_raw), row.names = row.names(count_raw))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
SCT_monocle_cds <- newCellDataSet(count_raw,phenoData = pd,featureData = fd,#lowerDetectionLimit = 0.5,
                                  expressionFamily = uninormal())

HSMM@assayData <- SCT_monocle_cds@assayData
#cds_test now maintains the trajectory constructed using the integrated assay
#but contains expression data from the RNA assay.
proliferate_markers<-c("MKI67","TOP2A","TK1")
Myeloid_Cell_maker<-c("HLA-B","VIM","KRT7","GATA3","PERP","XAGE2","CDH1","EGFR","PAGE4","HLA-G","PAPPA2","CGA","CSH2","ERVFRD-1")
maker_CBTs<- c("CDH1","EGFR","PAGE4","PEG10")
maker_EVTs<-c("HLA-G","MMP2","PAPPA2","DIO2","MMP11","FLT1","ITGA5","ADAM12","MCAM")
maker_STBs<-c("CSH1","CGA","CSH2","CSHL1","GH2","PSG2","HOPX","TFAP2A","CYP19A1","ERVFRD-1")

maker_gene<-maker_STBs
maker_plot<-plot_cell_trajectory(HSMM, markers =maker_gene, use_color_gradient = T,cell_size=0.1,show_branch_points = F)
maker_plot+scale_color_gradient2(low = "lightblue", mid = "darkred", high = "orange")+ theme(legend.position = "right")
maker_plot+geom_point(aes(color=value),size=0.5,na.rm = F) + theme(legend.position="top")+
  scale_color_gradient2(name = "RNA normalized Count",low = "lightblue",mid = "grey", high = "darkred")

saveRDS(hvg_HSMM, file = "/mnt/data/chenwei/gongchen/6.monocle2/CTRL_intergrated_monocle2_pseudotime_object_original_hvg.rds")
saveRDS(HSMM, file = "/mnt/data/chenwei/gongchen/6.monocle2/CTRL_intergrated_monocle2_pseudotime_object_pre_gene_call.rds")

names(seurat_data_list)
dtlist <- seurat_data_list[c("Arrest_1","Arrest_2","Arrest_3","Arrest_4")]
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 5000) #2976 
dtlist <- PrepSCTIntegration(object.list = dtlist,anchor.features = intfts)
anchors <- FindIntegrationAnchors(object.list = dtlist, normalization.method = "SCT", anchor.features = intfts)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated <- RunPCA(integrated)
ElbowPlot(object = integrated)
VizDimLoadings(object = integrated, dims = 1:2, reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:25)
integrated<- RunTSNE(integrated, dims = 1:25)
selected_cell <- FindNeighbors(object = integrated, dims = 1:25,verbose = FALSE)
selected_cell <- FindClusters(object = selected_cell,resolution = 0.6, verbose = FALSE)
DimPlot(object = selected_cell, reduction = "umap", group.by = "final_major_subgroup_brief",cols =ppCor)
DimPlot(object = selected_cell, reduction = "umap", group.by = "sample_code",cols =ppCor)
DimPlot(object = selected_cell, reduction = "umap", group.by = "sample_code",split.by ="sample_code",cols =ppCor,ncol=2)
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
saveRDS(selected_cell, file = "/mnt/data/chenwei/gongchen/6.monocle2/Abortion_Myeloid_Cell_seurat_pbject.rds")

table(selected_cell$monocle_names0)
#CTBs EVTs STBs 
#3017 4836  768 
table(selected_cell$final_major_subgroup_brief)
#CTBs_0 CTBs_1 CTBs_2 CTBs_3 EVTs_1 EVTs_2 EVTs_3 STBs_1 STBs_2 
#   298    444   1910    365   4094    353    389    388    380

##call DEGs inter group
DefaultAssay(selected_cell) <- "RNA"
selected_cell <- NormalizeData(selected_cell, verbose = FALSE)
Idents(object = selected_cell) <- "final_major_subgroup_brief" #'merge_cluster'
Allcluster.markers_all <- FindAllMarkers(object = selected_cell,min.pct = 0.25)
Allcluster.markers_all_default0.25<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=0.25,]
Allcluster.markers_all_default0.25_q01<-Allcluster.markers_all_default0.25[Allcluster.markers_all_default0.25$p_val_adj<0.1,]
DEGs_seurat_call<-Allcluster.markers_all_default0.25_q01

##step2 ：Construct CellDataSet object 10X的数据使用UMI count矩阵
DefaultAssay(selected_cell) <- "integrated"
count_raw <- as(as.matrix(selected_cell@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = selected_cell@meta.data)
fData <- data.frame(gene_short_name = row.names(count_raw), row.names = row.names(count_raw))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(count_raw,phenoData = pd,featureData = fd,#lowerDetectionLimit = 0.5,
                              expressionFamily = uninormal())#expressionFamily = negbinomial.size())
#expressionFamily参数用于指定表达矩阵的数据类型，有几个选项可以选择：
#稀疏矩阵用negbinomial.size()，FPKM值用tobit()，logFPKM值用gaussianff()

##step3 Size factor and dispersion will be used to normalize data
HSMM<-monocle_cds
#hvg called by seurat：
length(VariableFeatures(selected_cell))#4567
length(unique(DEGs_seurat_call$gene))
hvg_genes<-VariableFeatures(selected_cell)
ordering_genes<-hvg_genes[which(hvg_genes %in% unique(DEGs_seurat_call$gene))]
length(ordering_genes)#1712
#ordering_genes<-hvg_genes
hvg_HSMM <- setOrderingFilter(HSMM, ordering_genes)
## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
hvg_HSMM <- reduceDimension(hvg_HSMM,norm_method="none", 
                            reduction_method="DDRTree",
                            max_components=4,scaling=TRUE,verbose=TRUE,pseudo_expr=0)
## order cells change colors and theta to match your plot
hvg_HSMM <- orderCells(hvg_HSMM)
#Monocle does not know which the root state is for the tree, we need to define the root.
plot_cell_trajectory(hvg_HSMM, color_by = "final_major_subgroup_brief")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(hvg_HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)

#selected special HSMM
#DEGs_HSMM hvg_HSMM diff_HSMM disp_HSMM
HSMM<-hvg_HSMM
plot_cell_trajectory(HSMM, color_by = "monocle_names0",show_branch_points = FALSE,show_tree = TRUE,cell_size = 1) +
  scale_colour_manual(values=ppCor)+  theme(legend.position = "right")

#Step6 See what the trajectory looks like. for pre root identify
#Monocle does not know which the root state is for the tree, we need to define the root.
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="State")+scale_colour_manual(values=rev(ppCor))
plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")

plot_cell_trajectory(HSMM, color_by = "monocle_names0",cell_size = 0.75)+scale_colour_manual(values=ppCor)+facet_wrap(~monocle_names0, ncol = 3) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "final_major_subgroup_brief",cell_size = 0.75)+scale_colour_manual(values=ppCor)+facet_wrap(~final_major_subgroup_brief, ncol = 3) ##根据细胞聚类的进行着色
plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)+facet_wrap(~Phase, ncol = 3)

plot_cell_trajectory(HSMM, color_by = "Phase")+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by = "monocle_names0")+scale_colour_manual(values=ppCor)+facet_wrap(~sample_code, ncol =2)
plot_cell_trajectory(HSMM, color_by="sample_code",cell_size=0.1)+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="monocle_names0",cell_size=0.1,show_branch_points = F)+scale_colour_manual(values=ppCor)+facet_wrap(~Treat, ncol = 2)
plot_cell_trajectory(HSMM, color_by="final_major_subgroup_brief",cell_size=0.1,show_branch_points = F)+scale_colour_manual(values=ppCor)+facet_wrap(~final_major_subgroup_brief, ncol = 3)

table(pData(HSMM)$final_major_subgroup_brief,pData(HSMM)$State)
# we clearly see that State1 should be the root # so now we can reorder the cells  again
plot_cell_trajectory(HSMM, color_by="State",cell_size=0.1)+scale_colour_manual(values=ppCor)
plot_cell_trajectory(HSMM, color_by="State",show_branch_points = F)+scale_colour_manual(values=rev(ppCor))+facet_wrap(~State, ncol = 3)

HSMM <- orderCells(HSMM, root_state=1)
# plot with coloring by pseudotime
plot_cell_trajectory(HSMM, color_by="Pseudotime")

#We can compare the inferred pseudotime to the known sampling timepoints.
meta_plot<-selected_cell@meta.data;head(meta_plot);dim(meta_plot)
pData_plot<-pData(HSMM);head(pData_plot);dim(pData_plot)
meta_plot2<-merge(meta_plot,pData_plot[,c("Pseudotime","State")],by=0)
dim(meta_plot2)
cell_distribution<-ggplot(meta_plot2,aes(x = Pseudotime,y = final_major_subgroup_brief, colour = final_major_subgroup_brief)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = ppCor) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("cell subtypes") +
  ggtitle("Cells ordered by monocle2 pseudotime")
cell_distribution

ggplot(meta_plot2,aes(x=Pseudotime,fill=Treat))+geom_density()+scale_fill_manual(values = ppCor[5])+ggtitle("Abortion Cells ordered by monocle2 pseudotime")
ggplot(meta_plot2,aes(x=Pseudotime,fill=final_major_subgroup_brief,alpha = 1/10))+geom_rug()+geom_density()+scale_fill_manual(values = ppCor)
ggplot(meta_plot2,aes(x=Pseudotime,fill=final_major_subgroup_brief,alpha = 1/10))+geom_rug()+geom_density()+scale_fill_manual(values = ppCor)+ theme_classic()+facet_wrap(~final_major_subgroup_brief, ncol = 3)

ggplot(meta_plot2,aes(x=Pseudotime,colour=final_major_subgroup_brief))+geom_density()+scale_color_manual(values = ppCor) + theme_classic() +facet_wrap(~final_major_subgroup_brief, ncol = 3)
ggplot(meta_plot2,aes(x=Pseudotime,colour=final_major_subgroup_brief,alpha = 1/10))+geom_rug()+geom_density()+scale_colour_manual(values = ppCor)

## replace the exprs cds matrix value with orignal RNA counts after ordercells based on Seurat integrated matrix?
#Below I mae another cell_data_set object but this time using the RNA assay for the
#expression matrix. I also have to adjust the gene_metadata.
DefaultAssay(selected_cell) <- "RNA"
# Normalize RNA data for visualization purposes
selected_cell <- NormalizeData(selected_cell, verbose = FALSE)
#count_raw <- as(as.matrix(selected_cell@assays$SCT@data),'sparseMatrix')
count_raw <- as(as.matrix(selected_cell@assays$RNA@data),'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = selected_cell@meta.data)
fData <- data.frame(gene_short_name = row.names(count_raw), row.names = row.names(count_raw))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
SCT_monocle_cds <- newCellDataSet(count_raw,phenoData = pd,featureData = fd,#lowerDetectionLimit = 0.5,
                                  expressionFamily = uninormal())

HSMM@assayData <- SCT_monocle_cds@assayData
#cds_test now maintains the trajectory constructed using the integrated assay
#but contains expression data from the RNA assay.
proliferate_markers<-c("MKI67","TOP2A","TK1")
Myeloid_Cell_maker<-c("HLA-B","VIM","KRT7","GATA3","PERP","XAGE2","CDH1","EGFR","PAGE4","HLA-G","PAPPA2","CGA","CSH2","ERVFRD-1")
maker_CBTs<- c("CDH1","EGFR","PAGE4","PEG10")
maker_EVTs<-c("HLA-G","MMP2","PAPPA2","DIO2","MMP11","FLT1","ITGA5","ADAM12","MCAM")
maker_STBs<-c("CSH1","CGA","CSH2","CSHL1","GH2","PSG2","HOPX","TFAP2A","CYP19A1","ERVFRD-1")

maker_gene<-maker_STBs
maker_plot<-plot_cell_trajectory(HSMM, markers =maker_gene, use_color_gradient = T,cell_size=0.1,show_branch_points = F)
maker_plot+scale_color_gradient2(low = "lightblue", mid = "darkred", high = "orange")+ theme(legend.position = "right")
maker_plot+geom_point(aes(color=value),size=0.5,na.rm = F) + theme(legend.position="top")+
  scale_color_gradient2(name = "RNA normalized Count",low = "lightblue",mid = "grey", high = "darkred")

saveRDS(hvg_HSMM, file = "/mnt/data/chenwei/gongchen/6.monocle2/Abortion_intergrated_monocle2_pseudotime_object_original_hvg.rds")
saveRDS(HSMM, file = "/mnt/data/chenwei/gongchen/6.monocle2/Abortion_intergrated_monocle2_pseudotime_object_pre_gene_call.rds")


##Find genes that change along pseudotime
#Now we run a new differential gene test which test for changes along pseudotime.
#找随拟时序变化的差异基因，以及不同state之间的差异基因。
#这两个都是monocle里面的概念，与seurat里面找的差异基因不同。
#nocle中differentialGeneTest()函数可以按条件进行差异分析
#我们可以按一定的条件筛选基因后进行差异分析，全部基因都输入会耗费比较长的时间。建议使用cluster差异基因或高变基因输入函数计算。
my_pseudotime_de <- differentialGeneTest(HSMM, fullModelFormulaStr = "~sm.ns(Pseudotime)")
my_pseudotime_de <-my_pseudotime_de[order(my_pseudotime_de$qval,decreasing = F),]

##Plot some top genes along pseudotime
# take top 8 genes
plotgenes <- rownames(my_pseudotime_de[order(my_pseudotime_de$qval)[1:8],])
plot_genes_in_pseudotime(HSMM[plotgenes,], color_by="monocle_names0",ncol=3)+scale_colour_manual(values=ppCor)
#Clustering genes by pseudotemporal expression pattern
sig_pseudo_gene_names <- row.names(subset(my_pseudotime_de, qval < 0.01))# qval < 1e-5
length(sig_pseudo_gene_names)
plot_pseudotime_heatmap(HSMM[sig_pseudo_gene_names[1:50], ],num_clusters =3,show_rownames=TRUE)+ggtitle("Top 50 gene cluster by pseudotemporal expression pattern")

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