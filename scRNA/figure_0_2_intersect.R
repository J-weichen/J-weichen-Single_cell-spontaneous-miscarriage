rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

require("Matrix")
library(data.table)
library(Seurat)
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
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色) nejm，共8种
show_col(pal1)
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色) nejm，共8种
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

#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
library(future)
plan(strategy = "multicore", workers = 8)
options(future.globals.maxSize = 1500 * 1024^12)

########################
### Integrate
########################
integrated <- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713_integrated_seurat.rds")
######
plot0<-DimPlot(integrated, reduction = "pca")
plot1<-DimPlot(integrated, reduction = "umap")
plot2<-DimPlot(integrated, reduction = "tsne")
CombinePlots(plots = list(plot0,plot1,plot2),legend="top",ncol=3)

###build coldata for merge matrix
head(integrated@meta.data)
sample_info<-integrated@meta.data
sample_info$raw_id<-rownames(sample_info)
sample_info$gemgroup<-as.numeric(unlist(lapply(strsplit(rownames(sample_info),"_"), function(x) x[2])))
sample_info$BARCODE<-as.character(unlist(lapply(strsplit(rownames(sample_info),"_"), function(x) x[1])))
sample_info$BARCODE<-as.character(unlist(lapply(strsplit(sample_info$BARCODE,"-"), function(x) x[1])))
sample_info$BARCODE_new <-paste0(sample_info$BARCODE,"_",sample_info$gemgroup)
tail(sample_info)
table(sample_info$gemgroup)
#  1     2     3     4     5     6     7     8     9 
#11046  6524  4544  4354  3872  3092  8751  4215  3590 
colnames(sample_info)
sample_info_new<-sample_info[,-c(1,11:41)]
#mydata_list<-c("N1_out","N2_out","N3_out","N4_out","N5_out","A1_out","A3_out","A4_out","A8_out")

coldata<-data.frame(sample=c("N1","N2","N3","N4","N5","A1","A3","A4","A8"),
                    sample_code=c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","Arrest_1","Arrest_2","Arrest_3","Arrest_4"),
                    Treat=c(rep("CTRL",5),rep("Abortion",4)),gemgroup=c(1:9))

head(coldata);head(sample_info_new)
######merge meta information#####
meta_data_new<-merge(x =sample_info_new,y =coldata, by = "gemgroup")
rownames(meta_data_new)<-meta_data_new$raw_id
str(meta_data_new);tail(meta_data_new);dim(meta_data_new)##49988     16

########replace the old metadata
integrated2<-integrated#temp save
dim(meta_data_new);dim(integrated@meta.data)
integrated@meta.data<-meta_data_new
integrated$Treat<-factor(integrated$Treat,levels = c("Abortion","CTRL"))
integrated$sample<-factor(integrated$sample,levels = c("N1","N2","N3","N4","N5","A1","A3","A4","A8"))
integrated$sample_code<-factor(integrated$sample_code,levels = c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","Arrest_1","Arrest_2","Arrest_3","Arrest_4"))
table(integrated$Treat)
#Abortion     CTRL 
## 19648    30340 
P1<-DimPlot(integrated, group.by = "Treat",cols = ppCor)
P2<-DimPlot(integrated, group.by = "Treat",split.by ="Treat",cols = ppCor)
P3<-DimPlot(integrated, group.by = "sample",split.by ="sample",cols = ppCor,ncol = 3)
P4<-DimPlot(integrated, group.by = "sample_code",split.by ="sample_code",cols = ppCor,ncol = 3)

ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/Group_umap_Treat1.pdf"),P1,width=6, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/Group_umap_Treat2.pdf"),P2,width=12, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/Group_umap_sample.pdf"),P3,width=12, height=12)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/Group_umap_sample_code.pdf"),P4,width=12, height=12)

#DefaultAssay(integrated) <- "RNA"
# Normalize RNA data for visualization purposes
#integrated <- NormalizeData(integrated, verbose = FALSE)
DefaultAssay(integrated) <- "SCT"
#细胞周期评估
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
length(s.genes);length(g2m.genes)#43； 54
integrated <- CellCycleScoring(integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) #presumably running on SCT assay
##Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
##Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
integrated$CC.Difference <- integrated$S.Score-integrated$G2M.Score

Cell_score_plot1<-DimPlot(integrated,reduction = "umap",group.by= "Phase",split.by = "Phase",cols =ppCor)
Cell_score_plot2<-DimPlot(integrated,reduction = "pca",group.by= "Phase",split.by = "Phase",cols =ppCor)
Cell_score_plot1/Cell_score_plot2
DimPlot(integrated,reduction = "umap",group.by= "Phase",cols =ppCor)+ ggtitle('Cell cycle evaluation') 
DimPlot(integrated,reduction = "umap",group.by= "Phase",split.by= "Phase",cols =ppCor)+ ggtitle('Cell cycle evaluation') 
DimPlot(integrated,reduction = "umap",group.by= "Phase",split.by= "Treat",cols =ppCor)+ ggtitle('Cell cycle evaluation') 

# Visualize the distribution of cell cycle markers across
RidgePlot(integrated, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
FeaturePlot(object = integrated, features = c("PCNA", "TOP2A", "MCM6", "MKI67"),cols= c("grey", "red"),ncol=2)
table(integrated$Phase)
#  G1      G2M     S 
#30557  8085 11346 
table(integrated$Phase,integrated$Treat)

#高变基因中的cell cycle 基因
library(VennDiagram)
VRF_genes<-VariableFeatures(integrated)
length(VRF_genes)
cell_cycle_genes<-c(s.genes,g2m.genes)

grid.newpage(); #清空画板，开始画新图
grid.draw(venn.diagram(list(VRF_genes=VRF_genes,cell_cycle_genes=cell_cycle_genes),fill=c(ppCor[2],ppCor[1]), filename = NULL))
dev.off();grid.newpage(); #清空画板，开始画新图


#关注状态直观绘制
DimPlot(object = integrated, reduction = "umap",label = TRUE) + NoLegend() + ggtitle('Cell_sctransform') 
##sample information plot in generation
#plots <- DimPlot(integrated, group.by = c("Family","Tissue","cell.type","Age_group","gender","Age_gender_group"), combine = FALSE,cols = ppCor)
#plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(ncol =1, override.aes = list(size = 5))))
#CombinePlots(plots)

# Cluster the cells
#Seurat采用的是graph-based聚类方法，k-means方法在V3中已经不存在了
DefaultAssay(integrated) <- "integrated"
integrated <- FindNeighbors(object = integrated, dims = 1:50,verbose = FALSE)
integrated <- FindClusters(object = integrated, verbose = FALSE) #default=0.8
#Cell <- FindClusters(object = Cell, resolution = 1.5,verbose = FALSE

##different resolution
#integrated <- FindClusters( object = integrated,resolution = c(seq(.4,0.6,0.8,1.2,1.5,2,2.5,3,5)))

integrated <- FindClusters(object = integrated, resolution = 5,verbose = FALSE)
integrated <- FindClusters(object = integrated, resolution = 3,verbose = FALSE)
integrated <- FindClusters(object = integrated, resolution = 2.5,verbose = FALSE)
integrated <- FindClusters(object = integrated, resolution = 2,verbose = FALSE)
integrated <- FindClusters(object = integrated, resolution = 1.5,verbose = FALSE)
integrated <- FindClusters(object = integrated, resolution = 1.2,verbose = FALSE)
integrated <- FindClusters(object = integrated, resolution = 1,verbose = FALSE)
integrated <- FindClusters(object = integrated, verbose = FALSE)
integrated <- FindClusters(object = integrated,resolution = 0.6, verbose = FALSE)
integrated <- FindClusters(object = integrated, resolution = 0.4,verbose = FALSE)


plot1<-DimPlot(integrated, group.by ="integrated_snn_res.0.8",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.0.8') 
plot2<-DimPlot(integrated, group.by ="integrated_snn_res.1",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.1') 
plot3<-DimPlot(integrated, group.by ="integrated_snn_res.1.2",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.1.2') 
plot4<-DimPlot(integrated, group.by ="integrated_snn_res.1.5",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.1.5') 
plot5<-DimPlot(integrated, group.by ="integrated_snn_res.2",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.2') 
plot6<-DimPlot(integrated, group.by ="integrated_snn_res.2.5",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.2.5') 
plot7<-DimPlot(integrated, group.by ="integrated_snn_res.3",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.3') 
plot8<-DimPlot(integrated, group.by ="integrated_snn_res.5",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.5') 
plot_merge_eight<-CombinePlots(plots = list(plot1, plot2,plot3,plot4,plot5,plot6,plot7,plot8),ncol = 4,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/cluster_integrate_diff_res.pdf",plot_merge_eight,width=16, height=8)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/cluster_integrate_diff_res.png",plot_merge_eight,width=16, height=8)

plot9<-DimPlot(integrated, group.by ="integrated_snn_res.0.4",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.0.4') 
plot10<-DimPlot(integrated, group.by ="integrated_snn_res.0.6",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.0.6') 
plot_merge_two<-CombinePlots(plots = list(plot9,plot10),ncol = 2,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/cluster_integrate_diff_res2.pdf",plot_merge_two,width=10, height=5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/cluster_integrate_diff_res2.png",plot_merge_two,width=10, height=5)

length(unique(integrated@meta.data$integrated_snn_res.0.4));length(unique(integrated@meta.data$integrated_snn_res.0.6));#24 33
length(unique(integrated@meta.data$integrated_snn_res.0.8));length(unique(integrated@meta.data$integrated_snn_res.1));#40 43
length(unique(integrated@meta.data$integrated_snn_res.1.2));length(unique(integrated@meta.data$integrated_snn_res.1.5)); #48 52
length(unique(integrated@meta.data$integrated_snn_res.2));length(unique(integrated@meta.data$integrated_snn_res.2.5)); #62 70
length(unique(integrated@meta.data$integrated_snn_res.3));length(unique(integrated@meta.data$integrated_snn_res.5)) #76 107

#14=>20=>22=>30=>33=>36=>43=>52=>58=>82
#分群评估
#参考：https://www.jianshu.com/p/99cb6dc8de45  https://www.jianshu.com/p/a7a6b8b11e3c
library(clustree)
plot_clustree<-clustree(integrated@meta.data, prefix = "integrated_snn_res.")
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/clustree_dot_line.pdf",plot_clustree,width=10, height=20)

library(ggalluvial)
library(tidyverse)
##绘图时间长
head(integrated@meta.data)
plot_alluvial<-ggplot(data = integrated@meta.data,
       aes(axis1 = integrated_snn_res.0.4, axis2 = integrated_snn_res.0.6,axis3 = integrated_snn_res.0.8,axis4 = integrated_snn_res.1,axis5 = integrated_snn_res.1.2,axis6 = integrated_snn_res.1.5,axis7 = integrated_snn_res.2)) +
  scale_x_discrete(limits = c(paste0("integrated_snn_res.",seq(.4,1.5,.2))), expand = c(.01, .05)) +
  geom_alluvium(aes(fill = integrated_snn_res.1.5)) +
  geom_stratum() + geom_text(stat = "stratum", infer.label = TRUE) +
  #coord_polar()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("cell number in each cluster")
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/alluvial_stream.png",plot_alluvial,width=20, height=20)

saveRDS(integrated, file = "/mnt/data/chenwei/gongchen/3.seurat_result/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0723_integrated_seurat.rds")
DimPlot(object = integrated, reduction = "umap",label = TRUE) + NoLegend() + ggtitle('Cell_default0.8_sctransform') 

#合并注释大类群
target<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0723_integrated_seurat.rds")
#target<-integrated
DimPlot(object = target, reduction = "umap",label = TRUE) + NoLegend() + ggtitle('Cell_default0.6_sctransform') 

DefaultAssay(target) <- "RNA"
# Normalize RNA data for visualization purposes
target <- NormalizeData(target, verbose = FALSE)
selected_gene<-c("ERVFRD-1","RRM2","CCNB1","TAC3","PRG2","JAM2","SERPINE1","CD74","CTSS","LY2","MRC1")
#For selected gene
selected_gene<-c("BAX","BCL2","BCL2A1","BCL2L1", "BCL2L2","BCL3","BCL6","BID", "PUMA", "BIK","BMP4","BNIP1","BNIP2","BNIP3")  
selected_gene<-c("CASP1","CASP2","CASP3","CASP4", "CASP5","CASP6","CASP7","CASP8", "CASP9","CASP10","BAK","NOXA")  

#library(patchwork)
VlnPlot(target, features = selected_gene, split.by = "Treat",pt.size = 0,ncol=4)
plots <- VlnPlot(target, features = selected_gene, split.by = "Treat",pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

plot_select_gene<- FeaturePlot(object = target, features =selected_gene,cols= c("grey", "red"),ncol=4)
#variables were not found: BAK, NOXA ,PUMA
#ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/apoptosis_maker_for_intergrated_six_sample_after_filter_singlet_Cell.png", plot_major_cluster,width=20, height=20)
plot_select_gene2<- FeaturePlot(object = target, features =c("CASP8", "CASP4","BCL2A1","BNIP3"),cols= c("grey", "red"),split.by ="Treat" )
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/split_umap_apoptosis_maker_for_intergrated_nine_sample_after_filter_singlet_Cell.png", plot_select_gene2,width=10, height=20)

plot_select_gene3<- FeaturePlot(object = target, features =selected_gene,cols= c("grey", "red"),split.by ="Treat" )
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/BCL_gene_split_umap_apoptosis_maker_for_intergrated_nine_sample_after_filter_singlet_Cell.png", plot_select_gene3,width=10, height=5*length(selected_gene),limitsize = FALSE)
plot_select_gene4<- FeaturePlot(object = target, features =selected_gene,cols= c("grey", "red"),split.by ="Treat" )
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/CASP_gene_split_umap_apoptosis_maker_for_intergrated_nine_sample_after_filter_singlet_Cell.png", plot_select_gene4,width=10, height=5*length(selected_gene),limitsize = FALSE)

FeaturePlot(object = target, features =c("BNIP3"),cols= c("grey", "red"),split.by ="Treat" )

#For major cell population
plot_major_cluster<- FeaturePlot(object = target, features = c("VIM","HLA-B","KRT7","PERP","EPCAM","PAEP","DCN","DLK1","PECAM1","LYVE1","AIF1","CD3D","CD8A","IL7R","NKG7","CD79A","MS4A2","HBA1","PF4","PPBP"),cols= c("grey", "red"),ncol=4)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/major_maker_for_intergrated_nine_sample_after_filter_singlet_Cell.png", plot_major_cluster,width=20, height=20)
# Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe: LYVE1, MS4A2, HBA1, PPBP

#For no_immunine Cells 5*6
plot_no_immune_maker<-FeaturePlot(object = target, features = c("VIM","CDH1","EGFR","HLA-G","MMP11","CGA","CYP19A1","CSH2","ERVFRD-1",
                                          "DCN","DLK1","THY1","COL1A1","LAMA2","TIMP1",
                                          "MYH11","RGS5","PRL","IGFBP1","APOD","COL6A2",
                                          "PECAM1","LYVE1","EGFL7","PPBP","PF4","EPCAM","PAEP"),cols= c("grey", "red"),ncol=5)
FeaturePlot(object = target, features =c("LYVE1","AIF1"),cols= c("grey", "red"))

#For immunine Cells 6*7
plot_immune_maker<-FeaturePlot(object = target,features = c("MS4A1","CD79A","CD79B","JCHAIN","TCL1A","FCER2", 
                                         "CD3E", "CD3D", "CD8A","IL7R","SELL","CCR7","CCL5","S100A4","GZMK", 
                                         "NKG7", "GNLY","KLRB1","FGFBP2","XCL1",
                                         "AIF1","LYZ","CST3","CD14","FCGR3A","MS4A7","S100A8","FCN1","CLEC9A","CD1C",
                                         "CSF1R","CD163","CD209","CD69","LYVE1","APOE","MS4A3",
                                         "HBB","KIT"), pt.size = 0.2,cols= c("grey", "red"),ncol=6)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/no_immune_maker_for_intergrated_nine_sample_after_filter_singlet_Cell.png", plot_no_immune_maker,width=25, height=30)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/immune_maker_for_intergrated_nine_sample_after_filter_singlet_Cell.png", plot_immune_maker,width=30, height=35)
#troblast and no troblast
plot_maker0<-FeaturePlot(object = target, features = c("KRT7","VIM"),cols= c("grey", "red"))
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/troblast_maker_maker_for_intergrated_singlet_Cell_nine_sample.png", plot_maker0,width=10, height=5)

#immune maker 
plot_maker1<-FeaturePlot(object =  target, features = c("PTPRC","CD3E"),cols= c("grey", "red"))
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/whole_immune_maker_for_intergrated_singlet_Cell_nine_sample.png", plot_maker1,width=10, height=5)

#tissue resident markers 
plot_maker2<-FeaturePlot(object =  target, features = c("CD69","ITGA1","CD9"),cols= c("grey", "red"),ncol=3)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/tissue_resident_markers_for_intergrated_singlet_Cell_nine_sample.png", plot_maker2,width=15, height=5)

#Mitotic or  proliferating subpopulation: MKI67, TOP2A ,TK1;
plot_maker3<-FeaturePlot(object = target, features = c("MKI67","TOP2A","TK1"),cols= c("grey", "red"),ncol=3)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/proliferate_markers_for_intergrated_singlet_Cell_nine_sample.png", plot_maker3,width=15, height=5)

#plot for Trophoblast: 5 X 3
plot_maker_Trophoblast<-FeaturePlot(object = target, features = c("HLA-B","VIM","KRT7","GATA3","PERP","XAGE2","CDH1","EGFR","PAGE4","HLA-G","PAPPA2","CGA","CSH2","ERVFRD-1"),ncol=5,cols= c("grey", "red"))
#FeaturePlot(object = target, features = c("EPCAM","PAEP","CAPS"),cols= c("grey", "red"),ncol=2)
#Cytotrophoblast（CBTs）2X 2
plot_maker_CBTs<-FeaturePlot(object = target, features = c("CDH1","EGFR","PAGE4","PEG10"),ncol=2,cols= c("grey", "red"))
#Extravilous trophoblast（EVTs）3x3
plot_maker_EVTs<-FeaturePlot(object = target, features = c("HLA-G","MMP2","PAPPA2","DIO2","MMP11","FLT1","ITGA5","ADAM12","MCAM"),ncol=3,cols= c("grey", "red"))
#Syncytiotrophoblast (STBs) 4 x 3
plot_maker_STBs<-FeaturePlot(object = target, features = c("CSH1","CGA","CSH2","CSHL1","GH2","PSG2","HOPX","TFAP2A","CYP19A1","ERVFRD-1"),cols= c("grey", "red"))
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/whole_trophoblast_markers_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_Trophoblast,width=25, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/CBTs_markers_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_CBTs,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/EVTs_markers_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_EVTs,width=15, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/STBs_markers_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_STBs,width=20, height=15)


#for TNK cells
plot_maker_T_cell1<-FeaturePlot(object = target, features = c("CD3D","CD4","CD8A","KLRB1","NKG7","NCAM1","FGFBP2","FCGR3A","CXCR4","XCL1","PRF1","GNLY","GZMA","GZMB","GZMK"),cols= c("grey", "red"),ncol= 5)
plot_maker_T_cell2<-FeaturePlot(object = target, features = c("MS4A2","PTPRC","CD69","CD3D","CD4","CD8A","IL7R","CCR7","SELL","FOXP3", "KLRF1",
                                          "FGFBP2","NKG7","GNLY","GZMB","GZMK","GZMA","NCAM1","TRDC","TRGC2"),cols= c("grey", "red"),ncol= 5)
plot_maker_T_cell3<-FeaturePlot(object = target, features = c("PTPRC","HBA1","MS4A2","CD34","CD69","CD3D","CD8A","CD4","KLRB1","KLRF1","IL7R","CCR7","SELL","FOXP3",
                                          "NKG7","NCAM1", "FGFBP2","FCGR3A","GNLY","GZMB","GZMK","TRDV2","TRGV9","MKI67","TOP2A"),cols= c("grey", "red"),ncol= 5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/Tcell_markers_group1_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_T_cell1,width=25, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/Tcell_markers_group2_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_T_cell2,width=25, height=20)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/Tcell_markers_group3_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_T_cell3,width=25, height=25)

#T cell maker
plot_maker_T_cell<-FeaturePlot(object = target, features = c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B"),cols= c("grey", "red"),ncol=3)
plot_maker_TCR_cell<-FeaturePlot(object = target, features = c("TRAC","TRBC1","TRBC2","TRDC","TRGC1","TRGC2"),cols= c("grey", "red"),ncol=3)
plot_maker_TRV_cell<-FeaturePlot(object = target, features = c("TRAV1-2","TRDV1","TRDV2","TRDV3","TRGV3","TRGV4","TRGV9"),cols= c("grey", "red"),ncol=3)
#  The following requested variables were not found: TRDV1, TRDV3
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/Tcell_markers_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_T_cell,width=15, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/TCR_markers_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_TCR_cell,width=15, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/TRV_markers_for_intergrated_singlet_Cell_nine_sample.png", plot_maker_TRV_cell,width=15, height=10)

FeaturePlot(object = target, features = c("TK1","EGFL7","PAEP","MS4A2"),cols= c("grey", "red"))

#######identify cell subtypes############
DimPlot(object = target, group.by ="integrated_snn_res.0.6", label = TRUE) + NoLegend()# + ggtitle('Cell_default0.6_sctransform') 
Idents(object = target) <- "integrated_snn_res.0.6"
Idents(target)
target$seurat_cluster <- as.character(Idents(target))
target$merge_cluster <- as.character(Idents(target))

FeaturePlot(object = target, features = c("HBA1","HBB","KIT"),cols= c("grey", "red"))

#Erythrocyte
sub_cell<-subset(x = target,idents=c("15"))
target$merge_cluster[Cells(sub_cell)] <- "Erythrocyte_HBA1_pos_VIM_neg"
target$seurat_cluster[Cells(sub_cell)] <- "Ery"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 
#Mast cell: KIT
#sub_cell<-subset(x = target,idents=c("17"))
#target$merge_cluster[Cells(sub_cell)] <- "Mast_cells_KIT_pos" 
#DimPlot(object = target, group.by ="merge_cluster",label = TRUE) + NoLegend() 

#### Myeloid cells 

sub_cell<-subset(x = target,idents=c("5")) 
target$merge_cluster[Cells(sub_cell)] <- "MarcPhage_AIF1_pos_S100A8_pos"
target$seurat_cluster[Cells(sub_cell)] <- "dM1"
sub_cell<-subset(x = target,idents=c("8")) 
target$merge_cluster[Cells(sub_cell)] <- "MarcPhage_AIF1_pos_S100A8_neg"
target$seurat_cluster[Cells(sub_cell)] <- "dM2"

FeaturePlot(object = target, features = c("LYVE1","CCL4"),cols= c("grey", "red"))
sub_cell<-subset(x = target,idents=c("0","13")) 
target$merge_cluster[Cells(sub_cell)] <- "Hofbauer_Cell_AIF1_pos_LYVE1_pos"
target$seurat_cluster[Cells(sub_cell)] <- "HCs"

# TNKB_Cells:   CD3D+CD8A+
sub_cell<-subset(x = target,idents=c("14"))
target$merge_cluster[Cells(sub_cell)] <- "dNKsTB_Cells_CD3D_NKG7_CD79A"
target$seurat_cluster[Cells(sub_cell)] <- "TNKsBs"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 

# fibroblasts:   DLK1+COL1A1+
sub_cell<-subset(x = target,idents=c("6"))
target$merge_cluster[Cells(sub_cell)] <- "fibroblasts_DLK1_pos_COL1A1_pos"
target$seurat_cluster[Cells(sub_cell)] <- "FBs"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 

#Endothelial_Cell_PECAM1
sub_cell<-subset(x = target,idents=c("19"))
target$merge_cluster[Cells(sub_cell)] <- "Endothelial_Cell_PECAM1_pos"
target$seurat_cluster[Cells(sub_cell)] <- "Endo"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 

#Epithelial cells:PAEP 
#sub_cell<-subset(x = target,idents=c("24"))
#target$merge_cluster[Cells(sub_cell)] <- "Epithelial_cells_PAEP_pos"
#target$seurat_cluster[Cells(sub_cell)] <- "Epi"
#DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 

#Decidual_stromal_cells :DCN pos DLK1 neg
sub_cell<-subset(x = target,idents=c("18"))
target$merge_cluster[Cells(sub_cell)] <- "Decidual_stromal_cells_DCN_pos_DLK1_neg" 
target$seurat_cluster[Cells(sub_cell)] <- "dS"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 

##### Trophoblast_KRT7  
#Extravilous trophoblast（EVTs）
sub_cell<-subset(x = target,idents=c("2","3"))
target$merge_cluster[Cells(sub_cell)] <- "Extravilous_trophoblast_HLA-G_pos"
target$seurat_cluster[Cells(sub_cell)] <- "EVTs_1"
sub_cell<-subset(x = target,idents=c("11"))
target$merge_cluster[Cells(sub_cell)] <- "Extravilous_trophoblast_HLA-G_pos_PAPPA2_neg"
target$seurat_cluster[Cells(sub_cell)] <- "EVTs_2"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 

####Syncytiotrophoblast（STBs）
sub_cell<-subset(x = target,idents=c("17"))
target$merge_cluster[Cells(sub_cell)] <- "Syncytiotrophoblast_CGA_pos_CSH2_pos"
target$seurat_cluster[Cells(sub_cell)] <- "STBs_1"
sub_cell<-subset(x = target,idents=c("10"))
target$merge_cluster[Cells(sub_cell)] <- "Syncytiotrophoblast_CGA_pos_CSH2_neg"
target$seurat_cluster[Cells(sub_cell)] <- "STBs_2"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 
####Cytotrophoblast（CTBs）
sub_cell<-subset(x = target,idents=c("16"))
target$merge_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_neg_CDH1_pos"
target$seurat_cluster[Cells(sub_cell)] <- "CTBs_3"
sub_cell<-subset(x = target,idents=c("4","9"))
target$merge_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_neg"
target$seurat_cluster[Cells(sub_cell)] <- "CTBs_2"
sub_cell<-subset(x = target,idents=c("1"))
target$merge_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_pos_PCNA_pos"
target$seurat_cluster[Cells(sub_cell)] <- "CTBs_1"
sub_cell<-subset(x = target,idents=c("7"))
target$merge_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_pos_PCNA_neg"
target$seurat_cluster[Cells(sub_cell)] <- "CTBs_0"
###unclassical trophoblast（un_tros）
sub_cell<-subset(x = target,idents=c("12"))
target$merge_cluster[Cells(sub_cell)] <- "unclassical_trophoblasts_KRT7_pos"
target$seurat_cluster[Cells(sub_cell)] <- "un_tros"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 

#another round 
#plot merge group
DimPlot(object = target, group.by ="integrated_snn_res.0.8", label = TRUE) + NoLegend()# + ggtitle('Cell_default0.6_sctransform') 
DimPlot(object = target, group.by ="integrated_snn_res.1", label = TRUE) + NoLegend()# + ggtitle('Cell_default0.6_sctransform') 
DimPlot(object = target, group.by ="integrated_snn_res.3", label = TRUE) + NoLegend() 
DimPlot(object = target, group.by ="integrated_snn_res.5", label = TRUE) + NoLegend() 

Idents(object = target) <- "integrated_snn_res.1.2"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 
sub_cell<-subset(x = target,idents=c("24","27")) 
target$merge_cluster[Cells(sub_cell)] <- "Hofbauer_Cell_AIF1_pos_LYVE1_pos_MKI67_pos"
target$seurat_cluster[Cells(sub_cell)] <- "HCs_2"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 

#Extravilous trophoblast（EVTs）:PAPPA2 neg MKI67/TK1 pos
sub_cell<-subset(x = target,idents=c("21"))
target$merge_cluster[Cells(sub_cell)] <- "Extravilous_trophoblast_HLA-G_pos_PAPPA2_neg"
target$seurat_cluster[Cells(sub_cell)] <- "EVTs_2"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 

#Extravilous trophoblast（EVTs）
Idents(object = target) <- "integrated_snn_res.3"
sub_cell<-subset(x = target,idents=c("17")) 
target$merge_cluster[Cells(sub_cell)] <- "Hofbauer_Cell_AIF1_pos_LYVE1_pos"
target$seurat_cluster[Cells(sub_cell)] <- "HCs"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 
sub_cell<-subset(x = target,idents=c("36"))
target$merge_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_pos_PCNA_neg"
target$seurat_cluster[Cells(sub_cell)] <- "CTBs_0"
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE) + NoLegend() 


DimPlot(target, group.by = "merge_cluster",label = F,cols = ppCor_all2)
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE,cols = ppCor_all2) + NoLegend() 
DimPlot(object = target, reduction = "umap", split.by = "Treat",group.by = "seurat_cluster",label = TRUE,cols = ppCor_all2,ncol = 2)

#add nature refer meta
sample_info<-target@meta.data
head(sample_info)
table(sample_info$gemgroup)
sample_info2<-sample_info[,c("sample_code","gemgroup")]
colnames(sample_info2)<-c("sample_code","gemgroup_1")
head(sample_info2)
sample_info3<-distinct(sample_info2)

nature_ref_meta<-read.table(file="/mnt/data/chenwei/gongchen/3.seurat_result/metadata_gongchen_cell_integrated_after_nature_singlelet_all_nine_forced_cell_low_cell_double_both_remove_0723_seurat.txt",stringsAsFactors=F,row.names=1,header =T)
head(nature_ref_meta)
table(nature_ref_meta$gemgroup)

######merge meta information#####
meta_data_new<-merge(x =nature_ref_meta,y =sample_info3, by = "sample_code")
head(meta_data_new)

meta_data_new$raw_id2<-as.character(unlist(lapply(strsplit(meta_data_new$raw_id,"_"), function(x) x[1])))
meta_data_new$raw_id2 <-paste0(meta_data_new$raw_id2,"_",meta_data_new$gemgroup_1)
meta_data_new$BARCODE_new <-paste0(meta_data_new$BARCODE,"_",meta_data_new$gemgroup_1)
head(meta_data_new)
rownames(meta_data_new)<-meta_data_new$raw_id2
str(meta_data_new);tail(meta_data_new);dim(meta_data_new)##49988    44

######merge meta information#####
head(meta_data_new)
meta_data_new2<-meta_data_new[,c("seurat_cluster","merge_cluster")]
colnames(meta_data_new2)<-c("seurat_cluster_nature","merge_cluster_nature")
head(meta_data_new2);head(sample_info)
meta_data_new3<-merge(x =sample_info,y =meta_data_new2, by = 0)
head(meta_data_new3)
rownames(meta_data_new3)<-meta_data_new3$Row.names

meta_data_new3$major_group_raw<-as.character(meta_data_new3$seurat_cluster)
meta_data_new3$major_group_brief<-as.character(meta_data_new3$seurat_cluster)

table(meta_data_new3$major_group_raw)
table(meta_data_new3$major_group_brief)

meta_data_new3[which(meta_data_new3$seurat_cluster %in% c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","STBs_1","STBs_2","un_tros")),]$major_group_raw <- "Trophoblast_KRT7_pos_mix_Epithelial_Cell_EPCAM_PAEP"
meta_data_new3[which(meta_data_new3$seurat_cluster %in% c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","STBs_1","STBs_2","un_tros")),]$major_group_brief <- "Trophoblast"
meta_data_new3[which(meta_data_new3$seurat_cluster %in% c("dM1","dM2","HCs","HCs_2")),]$major_group_raw <- "Myeloid_Cell_AIF1_pos"
meta_data_new3[which(meta_data_new3$seurat_cluster %in% c("dM1","dM2","HCs","HCs_2")),]$major_group_brief <- "Myeloid_Cell"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("dS")),]$major_group_raw <- "Decidual_stromal_cells_DKK1_pos"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("dS")),]$major_group_brief <- "Decidual_stromal_cells"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("FBs")),]$major_group_raw <- "Fibroblasts_DLK1_pos_COL1A1_pos"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("FBs")),]$major_group_brief <- "Fibroblasts"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("Endo")),]$major_group_raw <- "Endothelial_Cell_PECAM1_pos"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("Endo")),]$major_group_brief <- "Endothelial_Cell"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("TNKsBs")),]$major_group_raw <- "T_NK_B_raw_Cell_CD4_NKG7_CD79A"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("TNKsBs")),]$major_group_brief <- "T_NK_B_Cell"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("Ery")),]$major_group_raw <- "Erythrocyte_HBA1"
meta_data_new3[which(meta_data_new3$seurat_cluster %in%  c("Ery")),]$major_group_brief <- "Erythrocyte"

head(meta_data_new3)

########replace the old metadata
target2<-target#temp save
dim(meta_data_new3);dim(target@meta.data)
target@meta.data<-meta_data_new3
table(target$seurat_cluster)
table(target$seurat_cluster_nature)
table(meta_data_new3$major_group_raw)
table(meta_data_new3$major_group_brief)

umap_plot1<-DimPlot(object = target, group.by ="seurat_cluster",label = TRUE,cols = ppCor_all2) +NoLegend()
umap_plot2<-DimPlot(object = target, group.by ="seurat_cluster_nature",label = TRUE,cols = ppCor_all2)  +NoLegend()
umap_plot1+umap_plot2

umap_plot3<-DimPlot(object = target, group.by ="major_group_raw",label = TRUE,cols = ppCor_all2)  +NoLegend()
umap_plot4<-DimPlot(object = target, group.by ="major_group_brief",label = TRUE,cols = ppCor_all2)  +NoLegend()
umap_plot3+umap_plot4

table(target$merge_cluster,target$Treat)
#           
#                                                                Abortion CTRL
#Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_neg              1450 5111
#Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_neg_CDH1_pos      460  638
#Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_pos_PCNA_neg      681 2789
#Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_pos_PCNA_pos      559 3099
#Decidual_stromal_cells_DCN_pos_DLK1_neg                              214  483
#dNKsTB_Cells_CD3D_NKG7_CD79A                                        1086  209
#Endothelial_Cell_PECAM1_pos                                           93   38
#Erythrocyte_HBA1_pos_VIM_neg                                         183  982
#Extravilous_trophoblast_HLA-G_pos                                   3620 2636
#Extravilous_trophoblast_HLA-G_pos_PAPPA2_neg                         529 2046
#fibroblasts_DLK1_pos_COL1A1_pos                                      689 2499
#Hofbauer_Cell_AIF1_pos_LYVE1_pos                                    4696 3044
#Hofbauer_Cell_AIF1_pos_LYVE1_pos_MKI67_pos                           295 1018
#MarcPhage_AIF1_pos_S100A8_neg                                        828 1606
#MarcPhage_AIF1_pos_S100A8_pos                                       3196  844
#Syncytiotrophoblast_CGA_pos_CSH2_neg                                 323 1492
#Syncytiotrophoblast_CGA_pos_CSH2_pos                                 395  585
#unclassical_trophoblasts_KRT7_pos                                    351 1221

saveRDS(target, file = "/mnt/data/chenwei/gongchen/3.seurat_result/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0812_integrated_seurat.rds")


target<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0812_integrated_seurat.rds")
Phase_data<-data.frame(table(target$Phase,target$Treat))

Phase_pieplot1<-ggplot(Phase_data[which(Phase_data$Var2=="CTRL"),],aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ 
  coord_polar("y",start=1) + 
  # geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ labs(title = "CTRL")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all)
Phase_pieplot2<-ggplot(Phase_data[which(Phase_data$Var2=="Abortion"),],aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ 
  coord_polar("y",start=1) + 
  # geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ labs(title = "Abortion")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all)
Phase_pieplot1+Phase_pieplot2 

#合并注释大类群
DimPlot(object = target, reduction = "umap",label = TRUE) 
DimPlot(object = target, group.by ="seurat_cluster",label = TRUE,cols = ppCor_all2) + NoLegend() + ggtitle('Cell_annotation') 
DimPlot(object = target, group.by ="seurat_cluster",cols = ppCor_all2)  + ggtitle('Cell_annotation') 

Metadate_plot<-target@meta.data
unique(Metadate_plot$seurat_cluster)
Metadate_plot$seurat_cluster<-factor(Metadate_plot$seurat_cluster,levels = c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","STBs_1","STBs_2","EVTs_1","EVTs_2","un_tros",
                                                               "FBs","dS","Endo","dM1","dM2","HCs","HCs_2","TNKsBs","Ery"))
unique(Metadate_plot$major_group_brief)
Metadate_plot$major_group_brief<-factor(Metadate_plot$major_group_brief,levels = c("Trophoblast","Fibroblasts","Decidual_stromal_cells","Endothelial_Cell","Myeloid_Cell","T_NK_B_Cell","Erythrocyte"))

#for CTRL without Erythrocyte_HBA1_pos_VIM_neg
Metadate_plot2<-subset(x = Metadate_plot,Treat == "CTRL" & seurat_cluster != "Ery")
Metadate_plot2$seurat_cluster<-factor(Metadate_plot2$seurat_cluster,levels = c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","STBs_1","STBs_2","EVTs_1","EVTs_2","un_tros",
                                                                             "FBs","dS","Endo","dM1","dM2","HCs","HCs_2","TNKsBs"))
dim(Metadate_plot2)
type_pieplot1<-ggplot(as.data.frame(table(Metadate_plot2$seurat_cluster)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ 
  coord_polar("y",start=1) + 
  # geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ labs(title = "CTRL")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all)
type_pieplot1 
#for Abortion   without Erythrocyte
Metadate_plot2<-subset(x = Metadate_plot,Treat == "Abortion" & seurat_cluster != "Ery")
Metadate_plot2$seurat_cluster<-factor(Metadate_plot2$seurat_cluster,levels = c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","STBs_1","STBs_2","EVTs_1","EVTs_2","un_tros",
                                                                               "FBs","dS","Endo","dM1","dM2","HCs","HCs_2","TNKsBs"))
dim(Metadate_plot2)
type_pieplot2<-ggplot(as.data.frame(table(Metadate_plot2$seurat_cluster)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ 
  coord_polar("y",start=1) + 
  #  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ labs(title = "Abortion")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all)
type_pieplot2
type_pieplot1+type_pieplot2

#propetion bar sperated without Ery
Metadate_plot3<-subset(x = Metadate_plot,seurat_cluster != "Ery")
Metadate_plot3$seurat_cluster<-factor(Metadate_plot3$seurat_cluster,levels = c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","STBs_1","STBs_2","EVTs_1","EVTs_2","un_tros",
                                                                               "FBs","dS","Endo","dM1","dM2","HCs","HCs_2","TNKsBs"))
cell.prop1<-as.data.frame(prop.table(table(Metadate_plot3$seurat_cluster, Metadate_plot3$Treat)))
colnames(cell.prop1)<-c("Cell_type","group","proportion")
ggplot(cell.prop1,aes(group,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=ppCor_all2)+coord_flip() 

#for CTRL without Erythrocyte_HBA1_pos_VIM_neg
Metadate_plot2<-subset(x = Metadate_plot,Treat == "CTRL" & major_group_brief != "Erythrocyte")
Metadate_plot2$major_group_brief<-factor(Metadate_plot2$major_group_brief,levels = c("Trophoblast","Fibroblasts","Decidual_stromal_cells","Endothelial_Cell","Myeloid_Cell","T_NK_B_Cell"))

dim(Metadate_plot2)
type_pieplot1<-ggplot(as.data.frame(table(Metadate_plot2$major_group_brief)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ 
  coord_polar("y",start=1) + 
  # geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ labs(title = "CTRL")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all)
type_pieplot1 
#for Abortion   without Erythrocyte
Metadate_plot2<-subset(x = Metadate_plot,Treat == "Abortion" & major_group_brief != "Erythrocyte")
Metadate_plot2$major_group_brief<-factor(Metadate_plot2$major_group_brief,levels = c("Trophoblast","Fibroblasts","Decidual_stromal_cells","Endothelial_Cell","Myeloid_Cell","T_NK_B_Cell"))

dim(Metadate_plot2)
type_pieplot2<-ggplot(as.data.frame(table(Metadate_plot2$major_group_brief)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ 
  coord_polar("y",start=1) + 
  #  geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ labs(title = "Abortion")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all)
type_pieplot2
type_pieplot1+type_pieplot2

#propetion bar sperated without Ery
Metadate_plot3<-subset(x = Metadate_plot,major_group_brief != "Erythrocyte")
Metadate_plot3$major_group_brief<-factor(Metadate_plot3$major_group_brief,levels = c("Trophoblast","Fibroblasts","Decidual_stromal_cells","Endothelial_Cell","Myeloid_Cell","T_NK_B_Cell"))

cell.prop1<-as.data.frame(prop.table(table(Metadate_plot3$major_group_brief, Metadate_plot3$Treat)))
colnames(cell.prop1)<-c("Cell_type","group","proportion")
ggplot(cell.prop1,aes(group,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=ppCor_all2)+coord_flip() 


#subset for each Cells groups
table(target$major_group_brief)
#Decidual_stromal_cells  Endothelial_Cell  Erythrocyte Fibroblasts  Myeloid_Cell  T_NK_B_Cell  Trophoblast
#      697                    131            1165         3188          15527           1295       27985

#For Trophoblast_mix_EECs_KRT7
sub_cell<-subset(x = target, subset = major_group_brief == "Trophoblast")
#sub_cell2<-subset(x = sub_cell, subset = cell.type  =="fetal")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Trophoblast <- sub_cell
table(Trophoblast$seurat_cluster)
#CTBs_0  CTBs_1  CTBs_2  CTBs_3  EVTs_1  EVTs_2  STBs_1  STBs_2 un_tros 
# 3470    3658    6561    1098    6256    2575     980    1815    1572

#For T_NK_B_Cell
sub_cell<-subset(x = target, subset = major_group_brief == "T_NK_B_Cell")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
T_NK_B_Cell <- sub_cell
table(T_NK_B_Cell$seurat_cluster)

#For Myeloid_Cell
sub_cell<-subset(x = target, subset = major_group_brief == "Myeloid_Cell")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Myeloid_Cell <- sub_cell
table(Myeloid_Cell$seurat_cluster)

#For Fibroblasts
sub_cell<-subset(x = target, subset = major_group_brief == "Fibroblasts")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Fibroblasts <- sub_cell
table(Fibroblasts$seurat_cluster)

#For Erythrocyte
sub_cell<-subset(x = target, subset = major_group_brief == "Erythrocyte")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Erythrocyte <- sub_cell
table(Erythrocyte$seurat_cluster)

#For Endothelial_Cell
sub_cell<-subset(x = target, subset = major_group_brief == "Endothelial_Cell")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Endothelial_Cell <- sub_cell
table(Endothelial_Cell$seurat_cluster)


#For Decidual_stromal_cells
sub_cell<-subset(x = target, subset = major_group_brief == "Decidual_stromal_cells")
DefaultAssay(sub_cell) <- "RNA"    ## very important
sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
ElbowPlot(object = sub_cell)
sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
Decidual_stromal_cells <- sub_cell
table(Decidual_stromal_cells$seurat_cluster)

Cell_submian_list<-c(list(Decidual_stromal_cells),list(Endothelial_Cell),list(Erythrocyte),list(Fibroblasts),list(Myeloid_Cell),list(T_NK_B_Cell),list(Trophoblast))
names(Cell_submian_list)<-c("Decidual_stromal_cells","Endothelial_Cell","Erythrocyte","Fibroblasts","Myeloid_Cell","T_NK_B_Cell","Trophoblast")
saveRDS(Cell_submian_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/Cell_submian_list_0813.rds")


## 查看每一类细胞的数目
sort(table(x = Idents(object = target)),decreasing = TRUE)
#查看某一cluster的具体细胞barcode
head(WhichCells(target,idents="2"))
# 提取某一cluster细胞查看。
##
head(subset(as.data.frame(target@active.ident),target@active.ident=="2"))
## or
head(subset(as.data.frame(target@meta.data),target@meta.data$merge_cluster=="Trophoblast_mix_EECs_KRT7"))
