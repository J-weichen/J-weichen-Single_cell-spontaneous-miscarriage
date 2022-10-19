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
library(ComplexHeatmap)

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
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:46))]
show_col(ppCor_all2)

###实接figure_0_4
#reading expression data
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

metadata<-target_final@meta.data
head(metadata);dim(metadata)##45800    53

#reading activity matrix
data_TFs_activity<- as.data.frame(read.table("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/ALLcell_all_genes/s3_ALLcell_all_genes.AUCell.txt",header =T,stringsAsFactors = FALSE,sep="\t",as.is=T))
#metadata <-read.csv("G:/ACE2/Cell research/metadata_new.txt",header =T,stringsAsFactors = FALSE,sep = "\t", row.names = 1)
head(data_TFs_activity);dim(data_TFs_activity)#49988   545
rownames(data_TFs_activity)<-data_TFs_activity$X
data_TFs_activity<-data_TFs_activity[,-1]
data_TFs_activity2<-data_TFs_activity[rownames(metadata),]
dim(data_TFs_activity2)
colnames(data_TFs_activity2)<-gsub('.{3}$', '', colnames(data_TFs_activity2))
range(data_TFs_activity2[,"KLF9"])

data_TFs_activity2<-as.data.frame(t(data_TFs_activity2))
tail(data_TFs_activity2[1:6,1:10])
data_TFs_activity2
##save for uppublish
write.table(data_TFs_activity2, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RAS_matrix_regulon_activity_all_singlet_and_not_mix_cells.txt",sep="\t",row.names = T)

write.table(data_TFs_activity2, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/data_TFs_activity_All_Abortion_CTRL.txt",sep="\t",row.names = F)
write.table(metadata, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/metadata_TFs_activity_All_Abortion_CTRL.txt",sep=""),sep="\t",row.names = F)

#数据构建
##方法1
#TFs_seurat <- CreateSeuratObject(counts = data_TFs_activity2, meta.data = metadata,project ="GC_TF_project")
#TFs_seurat
#VariableFeatures(TFs_seurat) <- rownames(TFs_seurat)
#TFs_seurat <- RunLSI(TFs_seurat, n = 30, scale.max = NULL)
#TFs_seurat <- RunTSNE(TFs_seurat, reduction = "lsi", dims =1:30)
#TFs_seurat <- RunUMAP(TFs_seurat, reduction = "lsi", dims =1:30)
#DimPlot(TFs_seurat, reduction = "umap", group.by = "final_major_subgroup_brief",cols = ppCor_all2) + ggtitle("scTF activity")
#DimPlot(TFs_seurat, reduction = "tsne", group.by = "final_major_subgroup_brief",cols = ppCor_all2) + ggtitle("scTF activity")

#方法2
TFs_seurat <- CreateSeuratObject(counts = data_TFs_activity2, meta.data = metadata,project ="GC_TF_project")
VariableFeatures(TFs_seurat) <- rownames(TFs_seurat)
TFs_seurat <- ScaleData(TFs_seurat, features = rownames(data_TFs_activity2))
TFs_seurat <- RunPCA(object = TFs_seurat,verbose = FALSE)
TFs_seurat <- RunUMAP(object = TFs_seurat, dims = 1:15,verbose = FALSE)
TFs_seurat <- RunTSNE(TFs_seurat, dims = 1:15)
DimPlot(TFs_seurat, reduction = "umap", group.by = "final_major_subgroup_brief",cols = ppCor_all2) + ggtitle("scTF activity")
DimPlot(TFs_seurat, reduction = "tsne", group.by = "final_major_subgroup_brief",cols = ppCor_all2) + ggtitle("scTF activity")

TFs_seurat <- FindNeighbors(object = TFs_seurat, dims = 1:20,verbose = FALSE)
TFs_seurat <- FindClusters(object = TFs_seurat,resolution = 0.4,verbose = FALSE)
TFs_seurat <- FindClusters(object = TFs_seurat,resolution = 0.6,verbose = FALSE)
TFs_seurat <- FindClusters(object = TFs_seurat,resolution = 0.8,verbose = FALSE)
TFs_seurat <- FindClusters(object = TFs_seurat,resolution = 1,verbose = FALSE)
TFs_seurat <- FindClusters(object = TFs_seurat,resolution = 1.2,verbose = FALSE)
TFs_seurat <- FindClusters(object = TFs_seurat,resolution =1.5,verbose = FALSE)
TFs_seurat <- FindClusters(object = TFs_seurat,resolution = 2,verbose = FALSE)
TFs_seurat <- FindClusters(object = TFs_seurat,resolution = 2.5,verbose = FALSE)
TFs_seurat <- FindClusters(object = TFs_seurat,resolution =3,verbose = FALSE)

B1<-DimPlot(object = TFs_seurat, reduction = "umap",group.by = "RNA_snn_res.0.4",label = TRUE) + NoLegend() +ggtitle('Cell_default0.4_TFac')
B2<-DimPlot(object = TFs_seurat, reduction = "umap",group.by = "RNA_snn_res.0.6",label = TRUE) + NoLegend() +ggtitle("Cell_default0.6_TFac")
B3<-DimPlot(object = TFs_seurat, reduction = "umap",group.by = "RNA_snn_res.0.8",label = TRUE) + NoLegend() +ggtitle("Cell_default0.8_TFac")
B4<-DimPlot(object = TFs_seurat, reduction = "umap",group.by = "RNA_snn_res.1",label = TRUE) + NoLegend() +ggtitle("Cell_default1_TFac")
B5<-DimPlot(object = TFs_seurat, reduction = "umap",group.by = "RNA_snn_res.1.2",label = TRUE) + NoLegend() +ggtitle("Cell_default1.2_TFac")
B6<-DimPlot(object = TFs_seurat, reduction = "umap",group.by = "RNA_snn_res.1.5",label = TRUE) + NoLegend() +ggtitle("Cell_default1.5_TFac")
B7<-DimPlot(object = TFs_seurat, reduction = "umap",group.by = "RNA_snn_res.2",label = TRUE) + NoLegend() +ggtitle("Cell_default2_TFac")
B8<-DimPlot(object = TFs_seurat, reduction = "umap",group.by = "RNA_snn_res.2.5",label = TRUE) + NoLegend() +ggtitle("Cell_default2.5_TFac")
B9<-DimPlot(object = TFs_seurat, reduction = "umap",group.by = "RNA_snn_res.3",label = TRUE) + NoLegend() +ggtitle("Cell_default3_TFac")

plot_for_dif_res<-CombinePlots(plots = list(B1,B2,B3,B4,B5,B6,B7,B8,B9),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/dif_res_dim15_for_regulon_TFs.png", plot_for_dif_res,width=30, height=30)

DimPlot(TFs_seurat, reduction = "umap", group.by = "RNA_snn_res.0.8",cols = ppCor_all2,label = T) + ggtitle("scTF activity")
DimPlot(TFs_seurat, reduction = "umap", group.by = "final_major_subgroup_brief",cols = ppCor_all2,label = T) + ggtitle("scTF activity")
DimPlot(TFs_seurat, reduction = "umap", group.by = "sample",cols = ppCor_all2,split.by = "sample",ncol=3) + ggtitle("scTF activity")
saveRDS(TFs_seurat, file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/TFs_regulon_activity_seurat_object.rds")

#reading object from GRN activity
TFs_seurat<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/TFs_regulon_activity_seurat_object.rds")
coldata<-TFs_seurat@meta.data
dim(coldata)#45800    63
table(coldata$RNA_snn_res.0.8,coldata$final_major_subgroup_brief)
DimPlot(TFs_seurat, reduction = "umap", group.by = "RNA_snn_res.0.8",cols = ppCor_all2,label = T) + ggtitle("scTF activity")
DimPlot(TFs_seurat, reduction = "umap", group.by = "RNA_snn_res.1",label = T) + ggtitle("scTF activity")
DimPlot(TFs_seurat, reduction = "umap", group.by = "RNA_snn_res.0.6",label = T) + ggtitle("scTF activity")

#change mata
target_final@meta.data<-TFs_seurat@meta.data
#re-annotation for each cluster
target_final$re_annotation_TFac <- as.character(target_final$RNA_snn_res.0.8)
subgroup_order<-c("0","1","6","24","17","19","8","15","10","12","3","7","2","18",
                  "5","11","20","27","9","14","21","4","26","22","13","23","16","25")
length(unique(subgroup_order))
table(target_final$re_annotation_TFac)#STCs 3885 
Idents(object = target_final) <- "re_annotation_TFac"
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order,ordered=TRUE)
DimPlot(object = target_final, label = TRUE) + NoLegend() 

major_submaker<-c("KRT7","PAGE4","PCNA","MKI67","EGFR","HLA-G","PAPPA2","CYP19A1","CSH2","ERVFRD-1",
                  "VIM","AIF1","LYVE1","PECAM1","DCN","DLK1","PAEP","NCAM1","CD3D","CD79A","MS4A2","HBA1")
#手动提取数据绘制所有细胞热图
levels(target_final)
DefaultAssay(target_final)<-"SCT"
FeaturePlot(object = target_final, features = c("PAPPA2"),split.by = "sample",cols= c("grey", "purple"))
FeaturePlot(object = target_final, features = c("HLA-G"),split.by = "sample",cols= c("grey", "red"))

mat0<- GetAssayData(target_final, slot = "data")
#mat0 <- GetAssayData(target, slot = "counts")
#mat0 <- log2(mat0 + 1)
#其他SCT的数据说明：https://www.jianshu.com/p/e639cc257d51 
#get genes and cluster information
#target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief,levels = subgroup_order)
gene_features <- major_submaker
#cluster_info <- target_final$final_major_subgroup_brief
cluster_info <- target_final$re_annotation_TFac
#sort the expression matrix 
mat1 <- as.matrix(mat0[as.character(gene_features), names(cluster_info)])
mat1[1:6,1:6];dim(mat1)
name_mean1<-as.data.frame(t(unlist(apply(mat1,1,mean))))
name_mean1<-as.data.frame(rbind(c(1:length(rownames(mat1))),name_mean1))
name_mean1<-as.data.frame(rbind(colnames(name_mean1),name_mean1))
rownames(name_mean1)<-c("Gname","row","mean")
dim(name_mean1);str(name_mean1)
name_mean1<-as.data.frame(t(name_mean1))
name_mean1$row<-as.numeric(as.character(name_mean1$row))
name_mean1$mean<-as.numeric(as.character(name_mean1$mean))
name_mean1$Gname<-as.character(name_mean1$Gname)
dim(name_mean1)
head(name_mean1)
#去重
name_mean1<-name_mean1[order(name_mean1$Gname,name_mean1$mean,decreasing = T),]
name_mean1<-name_mean1[!duplicated(name_mean1$Gname),]
name_mean1<-name_mean1[order(name_mean1$row),]
length(unique(gene_features))#2105
#name_mean1$row
head(name_mean1)
gene_features_uniq<-name_mean1$Gname
gene_features_dup<-gene_features[which(duplicated(gene_features))]
length(gene_features);length(gene_features_uniq);length(gene_features_dup)
#重新去unique数据框
mat2 <- mat1[gene_features_uniq,]
#手动scale
mat3<-scale(t(mat2), center=T,scale=T)
range(mat3)#-1.715613 63.262483
mat4<-t(mat3)
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1.5 #小于负数时，加括号！
#set color for cell types
col <- ppCor_all[1:length(levels(cluster_info))]
#colourCount<-length(unique(target$sub_cluster))
names(col) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col),labels = levels(cluster_info), labels_gp = gpar(cex = 0.5, col = "black")))

#legend adjust
col_fun1  <- circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red"))
col_fun2  <- circlize::colorRamp2(c(-2,0,2), c("purple", "black", "yellow"))
col_fun<-col_fun1
#show_col(col_fun(seq(-3, 3)))
pdf("/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/major_submaker_heatmap1.pdf",width = 15,height = 8)

plot_HP<-Heatmap(mat4,cluster_rows = FALSE,cluster_columns = FALSE,
                 show_column_names = FALSE,
                 #show_row_names = FALSE,
                 show_row_names = TRUE,
                 column_split = cluster_info,
                 top_annotation = top_anno,
                 #  right_annotation = row_anno,
                 column_title = NULL,
                 col = col_fun,
                 heatmap_legend_param = list(legend_height = unit(4, "cm"),grid_width = unit(0.4, "cm"),title = "log2(count+1)",title_position = "leftcenter-rot",at = c(-2, 0, 2), labels = c("low", "median", "high") ))

plot_HP
dev.off()

#re-annotation for each cluster
target_final$re_annotation_TFac <- as.character(target_final$final_major_subgroup_brief)
Idents(object = target_final) <- "RNA_snn_res.0.8"
DimPlot(TFs_seurat, reduction = "umap", group.by = "RNA_snn_res.0.8",cols = ppCor_all2,label = T) + ggtitle("scTF activity")

# Trophoblast_PERP_pos
sub_cell<-subset(x = target_final,idents=c("0"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "CTBs_1"
sub_cell<-subset(x = target_final,idents=c("1","6","24"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "CTBs_2"
sub_cell<-subset(x = target_final,idents=c("17"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "CTBs_3"
sub_cell<-subset(x = target_final,idents=c("19"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "CTBs_4"
sub_cell<-subset(x = target_final,idents=c("8"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "STBs_1"
sub_cell<-subset(x = target_final,idents=c("15"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "STBs_2"

sub_cell<-subset(x = target_final,idents=c("10"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "EVTs_1"
sub_cell<-subset(x = target_final,idents=c("12","3","7"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "EVTs_2"
sub_cell<-subset(x = target_final,idents=c("18","2"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "EVTs_3"
sub_cell<-subset(x = target_final,idents=c("26"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "Epi"
sub_cell<-subset(x = target_final,idents=c("5","11","20","27"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "Mycs"

#The myeloid lineage-specifc gene MS4A3, a known signature of granulocytic-monocytic progenitors (GMPs) which give rise to mast cell progenitors (MCP) and basophil progenitors.
#Mast cells account for ~0.5% of total immune cells in the decidua,which has been reported that positive for tryptases, KIT and MS4A2.
sub_cell<-subset(x = target_final,idents=c("23"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "Masts" 
sub_cell<-subset(x = target_final,idents=c("9"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "HCs"
sub_cell<-subset(x = target_final,idents=c("22"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "Endo"
sub_cell<-subset(x = target_final,idents=c("14"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "STCs"
sub_cell<-subset(x = target_final,idents=c("4","21"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "FBs"
sub_cell<-subset(x = target_final,idents=c("16","25"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "Ery"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("CTBs_0","MyCs"))
target_final$re_annotation_TFac[Cells(sub_cell)] <- "Ts"
DimPlot(target_final, group.by = "re_annotation_TFac",label = F,cols = ppCor_all2)
table(target_final$re_annotation_TFac)
# Bs CTBs_1 CTBs_2 CTBs_3 CTBs_4   Endo    Epi    Ery EVTs_1 EVTs_2 EVTs_3    FBs    HCs  Masts   Mycs    NKs STBs_1 STBs_2   STCs     Ts 
# 37   7837   7893    853    372    133     61    944   1848   6705   4746   2859   1902    113   4574    625   1955    903   1020    420 
saveRDS(target_final, file="/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_add_TFs_cluster_annotation.rds")

#change umap coordinates in GEXs seurat using that in TFs activity 
target_final[["umap"]]<-TFs_seurat[["umap"]]
DimPlot(target_final, group.by = "re_annotation_TFac",label = T,cols = ppCor_all2)

subgroup_order<-c("CTBs_1","CTBs_2","CTBs_3","CTBs_4","STBs_1","STBs_2","EVTs_1","EVTs_2","EVTs_3",
                 "Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order,ordered=TRUE)
table(target_final$re_annotation_TFac)
#CTBs_1 CTBs_2 CTBs_3 CTBs_4 STBs_1 STBs_2 EVTs_1 EVTs_2 EVTs_3    Epi   Endo    HCs   Mycs   STCs    FBs    NKs     Ts     Bs  Masts    Ery 
#  7837   7893    853    372   1955    903   1848   6705   4746     61    133   1902   4574   1020   2859    625    420     37    113    944 
DimPlot(target_final, group.by = "re_annotation_TFac",label = T,cols = ppCor_all2)

target_final$re_anno_TFac_major<-as.character(target_final$re_annotation_TFac)
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("CTBs_1","CTBs_2","CTBs_3","CTBs_4","STBs_1","STBs_2","EVTs_1","EVTs_2","EVTs_3"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "Trophoblast_KRT7"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("Epi"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "Epithelial_Cell_EPCAM_PAEP"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("Endo"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "Endothelial_Cell_PECAM1"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("HCs","Mycs"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "Myeloid_Cell_AIF1"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("STCs","FBs"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "Stromal_cells_DCN"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("NKs"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "Natural_killer_Cell_NCAM1"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("Ts"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "T_Cell_CD3D"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("Bs"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "B_Cell_CD79A"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("Masts"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "Mast_Cell_MS4A2"
sub_cell<-subset(x = target_final,subset = re_annotation_TFac %in% c("Ery"))
target_final$re_anno_TFac_major[Cells(sub_cell)] <- "Erythrocyte_HBA1"

table(target_final$re_anno_TFac_major)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major, 
                                              levels=c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN",
                                                       "Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP",
                                                       "Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A",
                                                       "Erythrocyte_HBA1","Mast_Cell_MS4A2"),ordered=TRUE)
DimPlot(target_final, group.by = "re_anno_TFac_major",label = T,cols = ppCor_all2)

saveRDS(target_final, file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object.rds")
