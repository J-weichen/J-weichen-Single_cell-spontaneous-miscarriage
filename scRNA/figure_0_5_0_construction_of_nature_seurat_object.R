rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(Seurat)
library(pheatmap) 
library(dplyr)
#library(export)
require(gridExtra)
library(reshape2)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(scales)
library(stringr)
library(RColorBrewer)

#调颜色
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor_all)
ppCor<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
#ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)
length(ppCor)

#读入表达数据####除了第一次外其他次直接读入
data_trans<- as.data.frame(read.table("/mnt/data/chenwei/data4/chenwei/script_10X/other_placental_article/placenta_nature/10X/E-MTAB-6701.processed.1/raw_data_10x.txt",header =T,stringsAsFactors = FALSE,sep = "\t"))
coldata<- as.data.frame(read.table("/mnt/data/chenwei/qinmen_BR/00.ref_data/single_cell_placental_nature/E-MTAB-6701.processed.2/meta_10x.txt",header =T,stringsAsFactors = FALSE,sep = "\t"))
tail(data_trans[1:6,1:6])
tail(coldata[1:6,])

dim(data_trans)#31764 64735
rownames(data_trans)<- data_trans$Gene
data_trans<-data_trans[,-1]
#remove duplication
data_trans$median=apply(data_trans,1,median) 
data_trans$symbol <-  unlist(lapply(strsplit(rownames(data_trans),"_"), function(x) x[1]))
data_trans<-data_trans[order(data_trans$symbol,data_trans$median,decreasing = T),]
data_trans<-data_trans[!duplicated(data_trans$symbol),]
rownames(data_trans)<-data_trans$symbol
data_trans[1:4,c(ncol(data_trans),ncol(data_trans)-1)] 
data_trans<-data_trans[,-c(ncol(data_trans),ncol(data_trans)-1)]
data_trans[1:4,1:4] 
dim(data_trans)# 31743 64734

dim(data_trans);dim(coldata)

#构建SeuratObject对象
Cell <- CreateSeuratObject(counts = data_trans, meta.data = coldata,project ="nature_project")

Cell #31743 features across 64734 samples
head(Cell@meta.data)

#########################
####pre_clean data
#Pre-View the data quality for un-fielter Cells
MT_gene_4<-grep(pattern = "^MT-", x = rownames(Cell), value = TRUE)#13 MT genes
RPS_gene_4<-grep(pattern = "^RPS", x = rownames(Cell), value = TRUE)

#rRNA genes percent
Rps_features <-  grep(pattern =  "^RPS", x = rownames(x = Cell[["RNA"]]),value = TRUE)
Cell$nCount_Rps<-colSums(x = GetAssayData(object = Cell,assay = "RNA", slot = "counts")[Rps_features, , drop = FALSE])
meta.data_zero<-Cell@meta.data
meta.data_zero$percent.rps<-100*meta.data_zero$nCount_Rps/meta.data_zero$nCount_RNA
Cell@meta.data<-meta.data_zero
range(Cell$percent.rps)#  1.071429 28.320463

#perform normalization
Cell <- SCTransform(Cell, verbose = FALSE)
Cell <- RunPCA(object = Cell,verbose = FALSE)
ElbowPlot(object = Cell)
Cell <- RunUMAP(object = Cell, dims = 1:20,verbose = FALSE)
Cell<- RunTSNE(Cell, dims = 1:20)
plot0<-DimPlot(Cell, reduction = "pca")
plot1<-DimPlot(Cell, reduction = "umap")
plot2<-DimPlot(Cell, reduction = "tsne")
CombinePlots(plots = list(plot0,plot1, plot2),legend="top",ncol = 3)
saveRDS(Cell, file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_all2.rds")

#saveRDS(Cell, file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_all.rds")



#selected Villous cell
DefaultAssay(Cell) <- "RNA"
Vill_Cell<-subset(x = Cell,subset =location=="Placenta") 
Vill_Cell#18547 samples
#perform normalization
Vill_Cell <- SCTransform(Vill_Cell, verbose = FALSE)
Vill_Cell <- RunPCA(object = Vill_Cell,verbose = FALSE)
ElbowPlot(object = Vill_Cell)
Vill_Cell <- RunUMAP(object = Vill_Cell, dims = 1:20,verbose = FALSE)
Vill_Cell<- RunTSNE(Vill_Cell, dims = 1:20)
plot0<-DimPlot(Vill_Cell, reduction = "pca")
plot1<-DimPlot(Vill_Cell, reduction = "umap")
plot2<-DimPlot(Vill_Cell, reduction = "tsne")
CombinePlots(plots = list(plot0,plot1, plot2),legend="top",ncol = 3)
saveRDS(Vill_Cell, file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_Vill.rds")

#plot
#Cell<-readRDS(file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_all.rds")
Cell<-readRDS(file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_Vill.rds")

DimPlot(Cell, group.by = "annotation",reduction ="umap",cols = rev(ppCor_all))
names(table(Cell$annotation))
#Cell2<-subset(x = Cell,subset = annotation %in% c("dS1","dS2","dS3",),invert = TRUE)

metadata_new<-Cell@meta.data
metadata_new$reannotation<-metadata_new$annotation
metadata_new[which(metadata_new$annotation %in% c("MO","dM1","dM2","dM3","DC1","DC2","Granulocytes")),]$reannotation <- "MyCs"
metadata_new[which(metadata_new$annotation %in% c("dNK p","dNK1","dNK2","dNK3","Tcells", "ILC3","Plasma","NK CD16-","NK CD16+")),]$reannotation <- "TsBsNK"
#metadata_new[which(metadata_new$annotation %in% c("dP1","dP2")),]$reannotation <- "dPs"
metadata_new[which(metadata_new$annotation %in% c("dS1","dS2","dS3")),]$reannotation <- "dSs"
metadata_new[which(metadata_new$annotation %in% c("Endo (f)","Endo (m)","Endo L")),]$reannotation <- "Endos"
metadata_new[which(metadata_new$annotation %in% c("Epi1","Epi2")),]$reannotation <- "Epis"
metadata_new[which(metadata_new$annotation %in% c("fFB1","fFB2")),]$reannotation <- "fFBs"
metadata_new[which(metadata_new$annotation %in% c("VCT")),]$reannotation <- "CTBs"
metadata_new[which(metadata_new$annotation %in% c("EVT")),]$reannotation <- "EVTs"
metadata_new[which(metadata_new$annotation %in% c("SCT")),]$reannotation <- "STBs"
Cell@meta.data<-metadata_new
table(Cell$reannotation)
names(table(Cell$reannotation))

subgroup_order2<-c("CTBs","STBs","EVTs","fFBs","HB","MyCs","Endos","Epis","dSs","TsBsNK")
Cell$reannotation<- factor(Cell$reannotation, levels=subgroup_order2,ordered=TRUE)
Umap_anno_plot<-DimPlot(Cell,group.by = "reannotation",label = T,cols = ppCor_all)
Umap_anno_plot
ggsave(Umap_anno_plot,file="/mnt/data/chenwei/qinmen_BR/3.methy_result/AMA_project/single_cell_placenta/sc_Nature_villous_Umap_anno_plot.png",width = 8,height = 7)
ggsave(Umap_anno_plot,file="/mnt/data/chenwei/qinmen_BR/3.methy_result/AMA_project/single_cell_placenta/sc_Nature_villous_Umap_anno_plot.pdf",width = 8,height = 7)

Umap_anno_plot2<-DimPlot(Cell,group.by = "reannotation",label = F,cols = ppCor_all)+NoLegend()+theme(plot.title=element_blank(),axis.title=element_blank(), axis.text = element_blank(), legend.title = element_blank())
Umap_anno_plot2
ggsave(Umap_anno_plot2,file="/mnt/data/chenwei/qinmen_BR/3.methy_result/AMA_project/single_cell_placenta/sc_Nature_villous_Umap_anno_plot2.png",width = 7,height = 7)

saveRDS(Cell, file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_Vill_reanno.rds")

