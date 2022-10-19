rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(Seurat)
require("Matrix")
library(data.table)
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

single_cell_placent_all<-readRDS(file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_all.rds")
DefaultAssay(single_cell_placent_all) <- "SCT"
seurat_data_list<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713.rds")
seurat_data_list2<-c(single_cell_placent_all,seurat_data_list)
dtlist <- seurat_data_list2
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 5000) #2976 
dtlist <- PrepSCTIntegration(object.list = dtlist,anchor.features = intfts)
anchors <- FindIntegrationAnchors(object.list = dtlist, normalization.method = "SCT", anchor.features = intfts)

integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated <- RunPCA(integrated)
ElbowPlot(object = integrated)
VizDimLoadings(object = integrated, dims = 1:2, reduction = "pca")

integrated <- RunUMAP(integrated, dims = 1:25)
integrated<- RunTSNE(integrated, dims = 1:25)
saveRDS(integrated, file = "/mnt/data/chenwei/gongchen/3.seurat_result/sc_nature_individual_singlelet_nine_sample_integrated_seurat.rds")

integrated <- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/sc_nature_individual_singlelet_nine_sample_integrated_seurat.rds")
######
plot0<-DimPlot(integrated, reduction = "pca")
plot1<-DimPlot(integrated, reduction = "umap")
plot2<-DimPlot(integrated, reduction = "tsne")
CombinePlots(plots = list(plot0,plot1,plot2),legend="top",ncol=3)

head(integrated@meta.data);tail(integrated@meta.data)
table(integrated@meta.data$orig.ident)
DimPlot(object = integrated, group.by ="annotation",label = TRUE) + NoLegend() 
meta_data_new<-integrated@meta.data
meta_data_new$data_source<-meta_data_new$annotation
meta_data_new[which(!(is.na(meta_data_new$annotation))),]$data_source<-"nature_data"
meta_data_new[which(is.na(meta_data_new$annotation)),]$data_source<-"gong_chen_data"
meta_data_new[which(is.na(meta_data_new$annotation)),]$annotation<-"gong_chen_data"
table(meta_data_new$data_source)
table(meta_data_new$annotation)
meta_data_new$data_source<- factor(meta_data_new$data_source,levels = c("gong_chen_data","nature_data"))
integrated@meta.data<-meta_data_new
DimPlot(object = integrated, group.by ="annotation",split.by ="data_source",label = TRUE) + NoLegend() 
DimPlot(object = integrated, group.by ="annotation",label = TRUE) + NoLegend() 


###build coldata for merge matrix
head(integrated@meta.data)
sample_info<-integrated@meta.data
sample_info$raw_id<-rownames(sample_info)
sample_info$gemgroup<-as.numeric(unlist(lapply(strsplit(rownames(sample_info),"_"), function(x) x[length(x)])))

sample_info_gc_data<-sample_info[which(sample_info$orig.ident == "Gongchen_filter_by_hand"),]
sample_info_gc_data$BARCODE<-as.character(unlist(lapply(strsplit(sample_info_gc_data$raw_id,"_"), function(x) x[1])))
sample_info_gc_data$BARCODE<-as.character(unlist(lapply(strsplit(sample_info_gc_data$BARCODE,"-"), function(x) x[1])))
sample_info_gc_data$BARCODE_new <-paste0(sample_info_gc_data$BARCODE,"_",sample_info_gc_data$gemgroup)

sample_info_sc_nature<-sample_info[which(sample_info$orig.ident != "Gongchen_filter_by_hand"),]
sample_info_sc_nature$BARCODE<-sample_info_sc_nature$raw_id
sample_info_sc_nature$BARCODE_new<-sample_info_sc_nature$raw_id
sample_info2<-rbind(sample_info_sc_nature,sample_info_gc_data)

head(sample_info2)
table(sample_info2$gemgroup)
#   1     2     3     4     5     6     7     8     9    10 
#64734 11046  6524  4544  4354  3872  3092  8751  4215  3590 
colnames(sample_info2)
sample_info_new<-sample_info2[,-c(15:45)]
tail(sample_info_new)
coldata<-data.frame(sample=c("Nature_ref","N1","N2","N3","N4","N5","A1","A3","A4","A8"),
                    sample_code=c("Nature_ref","CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","Arrest_1","Arrest_2","Arrest_3","Arrest_4"),
                    Treat=c("Nature_ref",rep("CTRL",5),rep("Abortion",4)),
                    Treat2=c("Nature_ref",rep("GC_data",9)),
                    gemgroup=c(1:10))

head(coldata);head(sample_info_new)
######merge meta information#####
meta_data_new<-merge(x =sample_info_new,y =coldata, by = "gemgroup")
rownames(meta_data_new)<-meta_data_new$raw_id
str(meta_data_new);tail(meta_data_new);dim(meta_data_new)##114722     23
########replace the old metadata
integrated2<-integrated#temp save
dim(meta_data_new);dim(integrated@meta.data)
integrated@meta.data<-meta_data_new
integrated$Treat<-factor(integrated$Treat,levels = c("Nature_ref","Abortion","CTRL"))
integrated$sample<-factor(integrated$sample,levels = c("N1","N2","N3","N4","N5","A1","A3","A4","A8","Nature_ref"))
integrated$sample_code<-factor(integrated$sample_code,levels = c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","Arrest_1","Arrest_2","Arrest_3","Arrest_4","Nature_ref"))
table(integrated$Treat)
#Nature_ref   Abortion       CTRL 
#  64734      19648      30340 
DimPlot(object = integrated, group.by ="annotation",split.by = "Treat") + NoLegend() 


P1<-DimPlot(integrated, group.by = "Treat",cols = ppCor)
P2<-DimPlot(integrated, group.by = "Treat",split.by ="Treat",cols = ppCor)
P3<-DimPlot(integrated, group.by = "sample",split.by ="sample",cols = ppCor,ncol = 3)
P4<-DimPlot(integrated, group.by = "sample_code",split.by ="sample_code",cols = ppCor,ncol = 3)

ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/individual_sc_nature_intergrate_Group_umap_Treat1.pdf"),P1,width=6, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/individual_sc_nature_intergrate_Group_umap_Treat2.pdf"),P2,width=12, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/individual_sc_nature_intergrate_Group_umap_sample.pdf"),P3,width=12, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/individual_sc_nature_intergrate_Group_umap_sample_code.pdf"),P4,width=12, height=6)

head(meta_data_new)
integrated$Treat2<-as.character(integrated$Treat)
subset_cell<-subset(x = integrated,subset =Treat %in% c("Abortion","CTRL")) 
integrated$Treat2[Cells(subset_cell)] <- "GC_data"

integrated$annotation2<-integrated$annotation
subset_cell<-subset(x = integrated,subset =Treat %in% c("Abortion")) 
integrated$annotation2[Cells(subset_cell)] <- "GC_Abortion"
subset_cell<-subset(x = integrated,subset =Treat %in% c("CTRL")) 
integrated$annotation2[Cells(subset_cell)] <- "GC_CTRL"
plot_group_split1<-DimPlot(object = integrated, group.by ="annotation2",split.by = "Treat2",label = TRUE,cols=ppCor_all) + NoLegend() 
plot_group_split2<-DimPlot(object = integrated, group.by ="annotation2",split.by = "Treat",label = TRUE,cols=ppCor_all) + NoLegend() 
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/individual_sc_nature_intergrate_group_split1.png",plot_sample_split1,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/individual_sc_nature_intergrate_group_split2.png",plot_group_split2,width=15, height=10)

saveRDS(integrated, file = "/mnt/data/chenwei/gongchen/3.seurat_result/annotad_sc_nature_individual_singlelet_nine_sample_integrated_seurat.rds")

DefaultAssay(integrated) <- "RNA"
# Normalize RNA data for visualization purposes
integrated <- NormalizeData(integrated, verbose = FALSE)
#For major cell population
plot_major_cluster<- FeaturePlot(object = integrated, features = c("VIM","HLA-B","KRT7","PERP","EPCAM","PAEP","DCN","DLK1","PECAM1","LYVE1","AIF1","CD3D","CD8A","IL7R","NKG7","CD79A","MS4A2","HBA1","PF4","PPBP"),cols= c("grey", "red"),ncol=4)

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/major_maker_for_intergrated_six_sample_after_filter_singlet_Cell.png", plot_major_cluster,width=20, height=20)

#target_gene<-c(BAX,BCL2,BCL2A1,BCL2L1,BCL2L2,BCL3,BCL6,BID,PUMA,BIK,BMP4,BNIP1,BNIP2,BNIP3,CASP1,CASP2,CASP3,CASP4,CASP5,CASP6,CASP7,CASP8,CASP9,CASP10,BAK,NOXA)  
FeaturePlot(object = Cell, features = target_gene,cols= c("grey", "red"), ncol=6)
FeaturePlot(object = Cell, features = target_gene,cols= c("grey", "red"),max.cutoff =2, ncol=6)