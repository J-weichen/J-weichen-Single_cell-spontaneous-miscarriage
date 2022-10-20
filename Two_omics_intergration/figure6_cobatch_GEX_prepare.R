rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
dyn.load('/mnt/data/chenwei/software/hdf/hdf5/lib/libhdf5.so.200')
dyn.load('/mnt/data/chenwei/software/hdf/hdf5/lib/libhdf5_hl.so.200')

library(hdf5r)
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)

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

gtf_file<-"/mnt/data/chenwei/gongchen/chip_seq_data/scRNA_ATCA_data_test/Homo_sapiens.GRCh38.91.gtf.gz"
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","CTBs_3","CTBs_4","STBs_1","STBs_2","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

#plot  group
DimPlot(target_final, group.by = "re_annotation_TFac",cols =ppCor_all2)
DimPlot(target_final, group.by = "re_anno_TFac_major",label = F,cols = ppCor_all2)
DimPlot(object = target_final, group.by = "re_annotation_TFac", split.by = "Treat",reduction = "umap",cols = ppCor_all2,ncol = 2)

cobatch_list<-list()
sample_name<-"final_six_sample" ##test line 
print(as.character(sample_name))

#str(peaks)
cell_peak_file<-paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/",sample_name,"_k27ac_hg38_mat.csv")
peaks <- read.table(cell_peak_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
peaks[1:6,1:6]
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks,annotation.file = gtf_file, seq.levels = c(1:22, "X", "Y"), 
                                            include.body = TRUE,upstream = 10000,downstream = 10000,
                                            verbose = TRUE)

#设置Seurat对象，将原始的peak counts保存到assay中，命名为"Cobatch"
placenta.Cobatch <- CreateSeuratObject(counts = peaks, assay = "Cobatch", project = "GC_Cobatch")
#增加基因活跃度矩阵，命名为"ACTIVITY".
placenta.Cobatch[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
#质检
DefaultAssay(placenta.Cobatch) <- "Cobatch"
eval_plot<- VlnPlot(object = placenta.Cobatch,features = c("nFeature_Cobatch", "nCount_Cobatch", "nCount_ACTIVITY", "nFeature_ACTIVITY"),
                    ncol = 4, pt.size = 0)
pdf(paste("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/",sample_name,"_reads_evaluation.pdf"))
eval_plot
dev.off()

#过滤掉cobatch数据中总count数低于500,大于12000的细胞？，这个阈值需要根据具体实验设置
range(placenta.Cobatch$nCount_Cobatch)#  148 4828
#placenta.Cobatch <- subset(placenta.Cobatch, subset = nCount_Cobatch > 500 & nCount_Cobatch < 12000)
placenta.Cobatch2<-placenta.Cobatch
placenta.Cobatch$tech <- "Cobatch"

DefaultAssay(placenta.Cobatch) <- "Cobatch"
#VariableFeatures(placenta.Cobatch)
#placenta.Cobatch <- FindVariableFeatures(placenta.Cobatch)
VariableFeatures(placenta.Cobatch) <- names(which(Matrix::rowSums(placenta.Cobatch) > 100))
placenta.Cobatch <- RunLSI(placenta.Cobatch, n = 30, scale.max = NULL)
placenta.Cobatch <- RunTSNE(placenta.Cobatch, reduction = "lsi", dims =2:30)
placenta.Cobatch <- RunUMAP(placenta.Cobatch, reduction = "lsi", dims =2:30)
DimPlot(placenta.Cobatch, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
DimPlot(placenta.Cobatch, reduction = "tsne") + NoLegend() + ggtitle("scATAC-seq")


DefaultAssay(placenta.Cobatch) <- "ACTIVITY"
placenta.Cobatch <- FindVariableFeatures(placenta.Cobatch)
placenta.Cobatch <- NormalizeData(placenta.Cobatch)
placenta.Cobatch <- ScaleData(placenta.Cobatch)


placenta.rna <-target_final
DefaultAssay(placenta.rna) <- "RNA"
# Perform standard analysis of each modality independently RNA analysis
placenta.rna <- NormalizeData(placenta.rna)
placenta.rna <- FindVariableFeatures(placenta.rna)
placenta.rna <- ScaleData(placenta.rna)
placenta.rna <- RunPCA(placenta.rna)
placenta.rna<- RunTSNE(placenta.rna, dims = 1:30)
placenta.rna <- RunUMAP(placenta.rna, dims = 1:30)
placenta.rna$tech <- "rna"

p1 <- DimPlot(placenta.Cobatch, reduction = "tsne") + NoLegend() + ggtitle("Cobatch")
p2 <- DimPlot(placenta.rna, group.by = "re_annotation_TFac", reduction = "tsne", cols = ppCor_all2,label = TRUE, repel = TRUE) + NoLegend() +   ggtitle("scRNA-seq")
pcom1<-CombinePlots(plots = list(p1, p2))
p3 <- DimPlot(placenta.Cobatch, reduction = "umap") + NoLegend() + ggtitle("Cobatch")
p4 <- DimPlot(placenta.rna, group.by = "re_annotation_TFac", reduction = "umap",cols = ppCor_all2, label = TRUE, repel = TRUE) + NoLegend() +   ggtitle("scRNA-seq")
pcom2<-CombinePlots(plots = list(p3, p4))
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_tsne_annotasion.pdf"), pcom1,width=16, height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_umap_annotasion.pdf"), pcom2,width=16, height=8)

# 只选择高变动的基因作为参考
genes.use <- VariableFeatures(placenta.rna)
refdata <- GetAssayData(placenta.rna, assay = "RNA", slot = "data")[genes.use, ]

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = placenta.rna, query = placenta.Cobatch, features = genes.use,
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
##After identifying anchors,we can transfer annotations from the scRNA-seq dataset onto the scATAC-seq cells. The annotations are stored in the seurat_annotations field, and are provided as input to the refdata parameter. The output will contain a matrix with predictions and confidence scores for each ATAC-seq cell.
celltype.predictions <- TransferData(anchorset = transfer.anchors,refdata = as.character(placenta.rna$re_annotation_TFac),weight.reduction = placenta.Cobatch[["lsi"]],dims = 2:30)
placenta.Cobatch <- AddMetaData(placenta.Cobatch, metadata = celltype.predictions)
#placenta.Cobatch$annotation_correct <- placenta.Cobatch$predicted.id == placenta.Cobatch$seurat_annotations
p1 <- DimPlot(placenta.Cobatch, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_single_umap_annotasion.pdf"), p1,width=8, height=8)
#p2 <- DimPlot(placenta.Cobatch, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
#p1 | p2
saveRDS(placenta.Cobatch,file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_placenta_Cobatch_seurats.rds"))
write.csv(as.data.frame(placenta.Cobatch@meta.data)[,1:7], file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_cobatch_list_sample_predicted_metadata.csv"),quote=F, row.names=T) 
#placenta.Cobatch<-readRDS(file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_placenta_Cobatch_seurats.rds"))
p_placenta1<-DimPlot(placenta.Cobatch, group.by = "orig.ident",cols = ppCor_all2, label = TRUE)  + ggtitle("orig.ident")
p_placenta2<-DimPlot(placenta.Cobatch, group.by = "predicted.id",cols = ppCor_all2, label = TRUE)  + ggtitle("Predicted annotation")
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_orig_ident_placenta_Cobatch_umap.pdf"),p_placenta1, width=9, height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_Predicted_placenta_Cobatch_umap.pdf"),p_placenta2, width=9, height=8)


imputation <- TransferData(anchorset = transfer.anchors,refdata = refdata, weight.reduction = placenta.Cobatch[["lsi"]])
# this line adds the imputed data matrix to the placenta.Cobatch object
placenta.Cobatch[["RNA"]] <- imputation

#合并两个的结果，然后就是scRNA-seq的常规分析。
coembed <- merge(x = placenta.rna, y = placenta.Cobatch)
# 标准化分析
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunTSNE(coembed, dims = 1:30)# Remove duplicates before running TSNE.
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$orig.ident
coembed$re_annotation_TFac <- ifelse(!is.na(coembed$re_annotation_TFac), coembed$re_annotation_TFac, coembed$orig.ident)
saveRDS(coembed,file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_coembed_object.rds"))

coembed<-readRDS(file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_coembed_object.rds"))
#在t-SNE上绘制结果
p1 <- DimPlot(coembed, reduction="umap", group.by = "tech")
p2 <- DimPlot(coembed, reduction="umap", group.by = "re_annotation_TFac",cols = ppCor_all2, label = TRUE, repel = TRUE)
#CombinePlots(list(p1, p2))
p3 <- DimPlot(coembed, reduction="umap", group.by = "tech", split.by = "tech")
p4 <- DimPlot(coembed, reduction="umap", split.by = "tech", group.by = "re_annotation_TFac", label = TRUE, repel = TRUE) + NoLegend()

DimPlot(coembed, reduction="tsne", split.by = "tech", group.by = "re_annotation_TFac",cols = ppCor_all2, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(coembed, reduction="umap", split.by = "tech", group.by = "re_annotation_TFac",cols = ppCor_all2, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(coembed, reduction="umap", group.by = "tech", split.by = "tech")
DimPlot(coembed, reduction="tsne", group.by = "tech", split.by = "tech")
DimPlot(coembed, reduction="tsne", group.by = "tech")
DimPlot(coembed, reduction="umap", group.by = "tech")
DimPlot(coembed, reduction="umap", group.by = "sample")

ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap1.pdf"), p1,width=8, height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap2.pdf"), p2,width=8, height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap3.pdf"), p3,width=16,height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap4.pdf"), p4,width=16,height=8)
