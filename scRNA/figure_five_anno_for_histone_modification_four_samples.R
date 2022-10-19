#Seurat 3.X版本能够整合scRNA-seq和scATAC-seq, 主要体现在：
##1)基于scRNA-seq的聚类结果对scATAC-seq的细胞进行聚类
##2)scRNA-seq和scATAC-seq共嵌入(co-embed)分析
#整合步骤包括如下步骤:
#1)从ATAC-seq中估计RNA-seq表达水平，即从ATAC-seq reads定量基因表达活跃度
#2)使用LSI学习ATAC-seq数据的内部结构
#3)鉴定ATAC-seq和RNA-seq数据集的锚点
#4)数据集间进行转移，包括聚类的标签，在ATAC-seq数据中推测RNA水平用于共嵌入分析
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
#基因活跃度定量
##将peak矩阵转成基因活跃度矩阵。
##Seurat做了一个简单的假设，基因活跃度可以通过简单的将落在基因区和其上游2kb的 "count相加" 得到基因活跃度,
##并且这个结果Cicero等工具返回gene-by-cell矩阵是类似的。
#读取参考注释
gtf_file<-"/mnt/data/chenwei/gongchen/chip_seq_data/scRNA_ATCA_data_test/Homo_sapiens.GRCh38.91.gtf.gz"
#cell_peak_file<-"/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/A4_hg38_mat.csv"
#placenta.Cobatch<- readRDS("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/A4_H3K27ac_seurat.rds")
target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
target_final<-subset(x = target,subset = final_major_subgroup== "Troph_mix_group",invert = TRUE)
major_order<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN",
               "Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order<-c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2","MyCs","HCs","STCs","FBs","Endo","Epi", "NKs","Ts","Bs","Ery","Masts")
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief, levels=subgroup_order,ordered=TRUE)
target_final$final_major_group_brief<- factor(target_final$final_major_group_brief, levels=major_order,ordered=TRUE)

#plot  group
DimPlot(target_final, group.by = "final_major_subgroup_brief",cols =ppCor_all2)
DimPlot(target_final, group.by = "final_major_group_brief",label = F,cols = ppCor_all2)
DimPlot(object = target_final, group.by = "final_major_group_brief", split.by = "Treat",reduction = "umap",cols = ppCor_all2,ncol = 2)

cobatch_list<-list()
#sample_name_list<-"Four_used_sample"
sample_name<-"Four_used_sample"   ##test line 
print(as.character(sample_name))
  # 读取peak
  #peaks <- Read10X_h5("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/atac_v1_placenta_10k_filtered_peak_bc_matrix.h5")
  #str(peaks)
cell_peak_file<-paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/",sample_name,"_k27ac_hg38_mat.csv")
peaks <- read.table(cell_peak_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
peaks[1:6,1:6]
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks,annotation.file = gtf_file, seq.levels = c(1:22, "X", "Y"), 
                                              include.body = TRUE,upstream = 10000,downstream = 10000,
                                              verbose = TRUE)
#activity.matrix是一个dgCMatrix对象，其中行为基因，列为细胞。
#因此如果对于Cicero的输出结果，只要提供相应的矩阵数据结构即可。
activity.matrix[1:6,1:6]
  
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
  #增加细胞的meta，该信息来自于scATAC-seq的CellRanger处理的singlecell.csv
  #meta <- read.table(singlecell_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  #meta <- meta[colnames(placenta.Cobatch), ]
  #placenta.Cobatch <- AddMetaData(placenta.Cobatch, metadata = meta)
  
  #过滤掉cobatch数据中总count数低于500,大于12000的细胞？，这个阈值需要根据具体实验设置
  range(placenta.Cobatch$nCount_Cobatch)# 502 11984
  placenta.Cobatch <- subset(placenta.Cobatch, subset = nCount_Cobatch > 500 & nCount_Cobatch < 12000)
  placenta.Cobatch2<-placenta.Cobatch
  placenta.Cobatch$tech <- "Cobatch"
  
  #数据预处理
  #我们需要对peak矩阵进行处理。
  #这里我们用的隐语义(Latent semantic indexing, LSI)方法对scATAC-seq数据进行降维。
  #该步骤学习scRNA-seq数据的内部结构，并且在转换信息时对锚点恰当权重的决定很重要。
  #根据 Cusanovich et al, Science 2015提出的LSI方法，他们搞了一个RunLSI函数。LSI的计算方法为TF-IDF加SVD。
  #我们使用在所有细胞中至少有100个read的peak，然后降维到50。该参数的选择受之前的scATAC-seq研究启发，所以没有更改，当然你可以把它改了。
  ## We exclude the first dimension as this is typically correlated with sequencing depth
  
  DefaultAssay(placenta.Cobatch) <- "Cobatch"
  #VariableFeatures(placenta.Cobatch)
  #placenta.Cobatch <- FindVariableFeatures(placenta.Cobatch)
  VariableFeatures(placenta.Cobatch) <- names(which(Matrix::rowSums(placenta.Cobatch) > 100))
  placenta.Cobatch <- RunLSI(placenta.Cobatch, n = 30, scale.max = NULL)
  placenta.Cobatch <- RunTSNE(placenta.Cobatch, reduction = "lsi", dims =2:30)
  placenta.Cobatch <- RunUMAP(placenta.Cobatch, reduction = "lsi", dims =2:30)
  DimPlot(placenta.Cobatch, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
  DimPlot(placenta.Cobatch, reduction = "tsne") + NoLegend() + ggtitle("scATAC-seq")
  #注: 要将placenta.Cobatch的默认assay切换成"Cobatch",  非线性降维可以选择UMAP或者t-SNE。
  
  #为了找到scATAC-seq数据集和scRNA-seq数据集之间的锚定点，我们需要对基因活跃度矩阵进行预处理
  #设置placenta.Cobatch的默认Assay为"ACTIVITY"， 然后寻找高表达的基因，对基因活跃度矩阵进行标准化和Scale。
  DefaultAssay(placenta.Cobatch) <- "ACTIVITY"
  placenta.Cobatch <- FindVariableFeatures(placenta.Cobatch)
  placenta.Cobatch <- NormalizeData(placenta.Cobatch)
  placenta.Cobatch <- ScaleData(placenta.Cobatch)
  
  #我们之前使用过Seurat对scRNA-seq数据进行预处理和聚类，
  #placenta.rna <- readRDS("/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/target_final_A4_forced_cell_low_cell_double_both_remove_seurat.rds")
  #target_final_N4<-subset(x = target_final,subset = sample== "N4")
  sample_list<-c("A3","A4","N4")
  placenta.rna<-subset(x = target_final,subset = sample %in% sample_list)
  #placenta.rna <- readRDS("/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/target_final_N4_forced_cell_low_cell_double_both_remove_seurat.rds")
  DefaultAssay(placenta.rna) <- "RNA"
  # Perform standard analysis of each modality independently RNA analysis
  placenta.rna <- NormalizeData(placenta.rna)
  placenta.rna <- FindVariableFeatures(placenta.rna)
  placenta.rna <- ScaleData(placenta.rna)
  placenta.rna <- RunPCA(placenta.rna)
  placenta.rna<- RunTSNE(placenta.rna, dims = 1:30)
  placenta.rna <- RunUMAP(placenta.rna, dims = 1:30)
  placenta.rna$tech <- "rna"
  
  #将scRNA-seq和scATAC-seq共同展示，对一些骨髓(myeloid)和淋巴(lymphoid)placenta中比较常见的聚类，其实是能从图中看出来。
  p1 <- DimPlot(placenta.Cobatch, reduction = "tsne") + NoLegend() + ggtitle("Cobatch")
  p2 <- DimPlot(placenta.rna, group.by = "final_major_subgroup_brief", reduction = "tsne", cols = ppCor_all2,label = TRUE, repel = TRUE) + NoLegend() +   ggtitle("scRNA-seq")
  pcom1<-CombinePlots(plots = list(p1, p2))
  p3 <- DimPlot(placenta.Cobatch, reduction = "umap") + NoLegend() + ggtitle("Cobatch")
  p4 <- DimPlot(placenta.rna, group.by = "final_major_subgroup_brief", reduction = "umap",cols = ppCor_all2, label = TRUE, repel = TRUE) + NoLegend() +   ggtitle("scRNA-seq")
  pcom2<-CombinePlots(plots = list(p3, p4))
  ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_tsne_annotasion.pdf"), pcom1,width=16, height=8)
  ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_umap_annotasion.pdf"), pcom2,width=16, height=8)
  
  ##Annotate scATAC-seq cells via label transfer
  # 只选择高变动的基因作为参考
  genes.use <- VariableFeatures(placenta.rna)
  refdata <- GetAssayData(placenta.rna, assay = "RNA", slot = "data")[genes.use, ]
  
  # Identify anchors
  transfer.anchors <- FindTransferAnchors(reference = placenta.rna, query = placenta.Cobatch, features = genes.use,
                                          reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
  ##After identifying anchors,we can transfer annotations from the scRNA-seq dataset onto the scATAC-seq cells. The annotations are stored in the seurat_annotations field, and are provided as input to the refdata parameter. The output will contain a matrix with predictions and confidence scores for each ATAC-seq cell.
  celltype.predictions <- TransferData(anchorset = transfer.anchors,refdata = as.character(placenta.rna$final_major_subgroup_brief),weight.reduction = placenta.Cobatch[["lsi"]],dims = 2:30)
  placenta.Cobatch <- AddMetaData(placenta.Cobatch, metadata = celltype.predictions)
  #placenta.Cobatch$annotation_correct <- placenta.Cobatch$predicted.id == placenta.Cobatch$seurat_annotations
  p1 <- DimPlot(placenta.Cobatch, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
  ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_single_umap_annotasion.pdf"), p1,width=8, height=8)
  #p2 <- DimPlot(placenta.Cobatch, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
  #p1 | p2
  saveRDS(placenta.Cobatch,file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_placenta_Cobatch_seurats.rds"))
  write.csv(as.data.frame(placenta.Cobatch@meta.data)[,1:7], file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_cobatch_list_sample_predicted_metadata.csv"),quote=F, row.names=T) 
  
  #此外，该prediction.score.max字段量化了与我们预测的注释相关的不确定性。
  #正确注释的单元格通常与高预测分数（> 90％）相关，而错误注释的单元格与预测分数（<50％）相关联。
  #不正确的分配也往往反映密切相关的细胞类型（即中间与原始 B 细胞）
  #predictions <- table(placenta.Cobatch$seurat_annotations, placenta.Cobatch$predicted.id)
  #predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
  #predictions <- as.data.frame(predictions)
  #p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  # theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #correct <- length(which(placenta.Cobatch$seurat_annotations == placenta.Cobatch$predicted.id))
  #incorrect <- length(which(placenta.Cobatch$seurat_annotations != placenta.Cobatch$predicted.id))
  #data <- FetchData(placenta.Cobatch, vars = c("prediction.score.max", "annotation_correct"))
  #p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  #  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
  #labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
  #p1 + p2
  
  
  #预测后可进行的分析结果
  ##在转移细胞类型标记之后，我们就可以进行细胞特异水平上的下游分析。
  ##举个例子，我们可以去找一些某些细胞类型特异的增强子,寻找富集的motif。
  ##目前这些分析Seurat3还不直接支持，还在调试中。
  
  #共嵌入(co-embedding)
  #最后，如果你想将所有的细胞一同展示，可以将scRNA-seq和scATAC-seq数据嵌入到相同的低维空间。
  #我们使用之前的锚点从scATAC-seq细胞中推断RNA-seq的值，后续分析就相当于两个单细胞数据的分析流程。
  #注意: 这一步只是为了可视化，其实不做也行。
  #选择用于估计的基因，可以是高变动基因，也可以是所有基因。
  
  
  #之后利用TransferData推断scATAC-seq在这些基因中的可能值
  #输出结果就是ATAC细胞对应的scRNA-seq矩阵
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
  coembed$final_major_subgroup_brief <- ifelse(!is.na(coembed$final_major_subgroup_brief), coembed$final_major_subgroup_brief, coembed$orig.ident)
  saveRDS(coembed,file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_coembed_object.rds"))
  
  #在t-SNE上绘制结果
  p1 <- DimPlot(coembed, reduction="umap", group.by = "tech")
  p2 <- DimPlot(coembed, reduction="umap", group.by = "final_major_subgroup_brief",cols = ppCor_all2, label = TRUE, repel = TRUE)
  #CombinePlots(list(p1, p2))
  p3 <- DimPlot(coembed, reduction="umap", group.by = "tech", split.by = "tech")
  p4 <- DimPlot(coembed, reduction="umap", split.by = "tech", group.by = "final_major_subgroup_brief", label = TRUE, repel = TRUE) + NoLegend()
  
  #从上面的结果中，你可能会发现某些细胞只有在一类技术中存在。
  #从巨噬细胞(megakaryocytes)到成熟的血小板细胞(pletelet)由于没有细胞核，无法被scATAC-seq技术检测到。
  #我们可以单独展示每个技术，进行检查
  DimPlot(coembed, reduction="tsne", split.by = "tech", group.by = "final_major_subgroup_brief",cols = ppCor_all2, label = TRUE, repel = TRUE) + NoLegend()
  DimPlot(coembed, reduction="umap", split.by = "tech", group.by = "final_major_subgroup_brief",cols = ppCor_all2, label = TRUE, repel = TRUE) + NoLegend()
  DimPlot(coembed, reduction="umap", group.by = "tech", split.by = "tech")
  DimPlot(coembed, reduction="tsne", group.by = "tech", split.by = "tech")
  DimPlot(coembed, reduction="tsne", group.by = "tech")
  
  ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap1.pdf"), p1,width=8, height=8)
  ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap2.pdf"), p2,width=8, height=8)
  ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap3.pdf"), p3,width=16,height=8)
  ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap4.pdf"), p4,width=16,height=8)

#using all RNA cells
sample_name<-"All_RNA_cells"
placenta.rna<-target_final
#placenta.rna <- readRDS("/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/target_final_N4_forced_cell_low_cell_double_both_remove_seurat.rds")
DefaultAssay(placenta.rna) <- "RNA"
# Perform standard analysis of each modality independently RNA analysis
placenta.rna <- NormalizeData(placenta.rna)
placenta.rna <- FindVariableFeatures(placenta.rna)
placenta.rna <- ScaleData(placenta.rna)
placenta.rna <- RunPCA(placenta.rna)
placenta.rna<-  RunTSNE(placenta.rna, dims = 1:30)
placenta.rna <- RunUMAP(placenta.rna, dims = 1:30)
placenta.rna$tech <- "rna"

#将scRNA-seq和scATAC-seq共同展示，对一些骨髓(myeloid)和淋巴(lymphoid)placenta中比较常见的聚类，其实是能从图中看出来。
p1 <- DimPlot(placenta.Cobatch, reduction = "tsne") + NoLegend() + ggtitle("Cobatch")
p2 <- DimPlot(placenta.rna, group.by = "final_major_subgroup_brief", reduction = "tsne",cols = ppCor_all2, label = TRUE, repel = TRUE) + NoLegend() +   ggtitle("scRNA-seq")
pcom1<-CombinePlots(plots = list(p1, p2))
p3 <- DimPlot(placenta.Cobatch, reduction = "umap") + NoLegend() + ggtitle("Cobatch")
p4 <- DimPlot(placenta.rna, group.by = "final_major_subgroup_brief", reduction = "umap",cols = ppCor_all2, label = TRUE, repel = TRUE) + NoLegend() +   ggtitle("scRNA-seq")
pcom2<-CombinePlots(plots = list(p3, p4))
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_tsne_annotasion.pdf"), pcom1,width=16, height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_umap_annotasion.pdf"), pcom2,width=16, height=8)

##Annotate scATAC-seq cells via label transfer
# 只选择高变动的基因作为参考
genes.use <- VariableFeatures(placenta.rna)
refdata <- GetAssayData(placenta.rna, assay = "RNA", slot = "data")[genes.use, ]

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = placenta.rna, query = placenta.Cobatch, features = genes.use,
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
##After identifying anchors,we can transfer annotations from the scRNA-seq dataset onto the scATAC-seq cells. The annotations are stored in the seurat_annotations field, and are provided as input to the refdata parameter. The output will contain a matrix with predictions and confidence scores for each ATAC-seq cell.
celltype.predictions <- TransferData(anchorset = transfer.anchors,refdata = as.character(placenta.rna$final_major_subgroup_brief),weight.reduction = placenta.Cobatch[["lsi"]],dims = 2:30)
placenta.Cobatch <- AddMetaData(placenta.Cobatch, metadata = celltype.predictions)
#placenta.Cobatch$annotation_correct <- placenta.Cobatch$predicted.id == placenta.Cobatch$seurat_annotations
p1 <- DimPlot(placenta.Cobatch, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_single_umap_annotasion.pdf"), p1,width=8, height=8)
#p2 <- DimPlot(placenta.Cobatch, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
#p1 | p2
saveRDS(placenta.Cobatch,file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_placenta_Cobatch_seurats.rds"))
write.csv(as.data.frame(placenta.Cobatch@meta.data)[,1:7], file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_cobatch_list_sample_predicted_metadata.csv"),quote=F, row.names=T) 

#之后利用TransferData推断scATAC-seq在这些基因中的可能值
#输出结果就是ATAC细胞对应的scRNA-seq矩阵
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
coembed$final_major_subgroup_brief <- ifelse(!is.na(coembed$final_major_subgroup_brief), coembed$final_major_subgroup_brief, coembed$orig.ident)
saveRDS(coembed,file = paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_coembed_object.rds"))

#在t-SNE上绘制结果
p1 <- DimPlot(coembed, reduction="umap", group.by = "tech")
p2 <- DimPlot(coembed, reduction="umap", group.by = "final_major_subgroup_brief", label = TRUE, repel = TRUE)
#CombinePlots(list(p1, p2))
p3 <- DimPlot(coembed, reduction="umap", group.by = "tech", split.by = "tech")
p4 <- DimPlot(coembed, reduction="umap", split.by = "tech", group.by = "final_major_subgroup_brief", label = TRUE, repel = TRUE) + NoLegend()

#从上面的结果中，你可能会发现某些细胞只有在一类技术中存在。
#从巨噬细胞(megakaryocytes)到成熟的血小板细胞(pletelet)由于没有细胞核，无法被scATAC-seq技术检测到。
#我们可以单独展示每个技术，进行检查
DimPlot(coembed, reduction="tsne", split.by = "tech", group.by = "final_major_subgroup_brief", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(coembed, reduction="umap", split.by = "tech", group.by = "final_major_subgroup_brief",cols = ppCor_all2, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(coembed, reduction="umap", group.by = "tech")
DimPlot(coembed, reduction="tsne", group.by = "tech")

ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap1.pdf"), p1,width=8, height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap2.pdf"), p2,width=8, height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap3.pdf"), p3,width=16,height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",sample_name,"_intergrated_umap4.pdf"), p4,width=16,height=8)





##此外，你还可以发现B细胞前体类型附近有一类细胞只由scATAC-seq的细胞构成。
#通过展示这些细胞在CellRanger分析结果提供的黑名单位置所占read数，可以推测出这类细胞可能是死细胞，或者是其他scRNA-seq无法检测的人为因素。
coembed$blacklist_region_fragments[is.na(coembed$blacklist_region_fragments)] <- 0
FeaturePlot(coembed, features = "blacklist_region_fragments", max.cutoff = 500)

##Seurat这种基于基因活跃度得分进行细胞类型预测，是否靠谱，开发者提供了如下几个证据
##1）总体预测得分(placenta.Cobatch$prediction.score.max)高，意味着用scRNA-seq定义细胞类型比较可靠
# 2）我们可以在scATC-seq降维结果中利用相同锚点的共嵌入分析，发现两类形态能很好的混合
# 3）将ATAC-seq数据根据聚类结果构建pseduo bulk, 发现和真实的bulk数据近似
#参考资料：https://satijalab.org/seurat/v3.0/atacseq_integration_vignette.html