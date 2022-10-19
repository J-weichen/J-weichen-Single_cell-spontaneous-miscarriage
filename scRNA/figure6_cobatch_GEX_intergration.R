rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ggplot2)
library(scales)
library(ggsci)
library(Seurat)

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

#reading GEX cell object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","CTBs_3","CTBs_4","STBs_1","STBs_2","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

DefaultAssay(target_final) <- "RNA"
target_final <- NormalizeData(target_final)
target_final <- FindVariableFeatures(target_final)

#reading cobatch data
placenta.Cobatch<-readRDS(file = "/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/final_six_sample_placenta_Cobatch_seurats.rds")

#########################
### Integrate for Cell2016_data
########################
dtlist <- list(placenta_RNA=target_final,placenta_Cobatch=placenta.Cobatch)
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 3000) #2976 
anchors <- FindIntegrationAnchors(object.list = dtlist, anchor.features = intfts)
integrated <- IntegrateData(anchorset = anchors)
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated)
ElbowPlot(object = integrated)
VizDimLoadings(object = integrated, dims = 1:2, reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:15)
integrated<- RunTSNE(integrated, dims = 1:15)
meta_new<-integrated@meta.data

meta_new[is.na(meta_new$tech),]$tech<-"scRNAseq"
meta_new[is.na(meta_new$tech),]$tech<-"scRNAseq"
meta_new$final_anno<-ifelse(is.na(meta_new$re_annotation_TFac),meta_new$predicted.id,meta_new$re_annotation_TFac)

table(meta_new$tech)
integrated@meta.data<-meta_new
integrated2<-integrated
saveRDS(integrated,file = "/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/placenta_Cobatch_and_RNA_intergrated_seurats.rds")
#integrated<-readRDS(file = "/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/placenta_Cobatch_and_RNA_intergrated_seurats.rds")
integrated@meta.data

metadata_for_target<-integrated@meta.data
target_coordinate<-data.frame(Embeddings(integrated[["umap"]]))
metadata_Cell_Umap_add <- merge(metadata_for_target,target_coordinate,by=0)
head(metadata_Cell_Umap_add)
metadata_Cell_Umap_add2<-metadata_Cell_Umap_add[,c("Row.names","tech","UMAP_1","UMAP_2")]
rownames(metadata_Cell_Umap_add2)<-metadata_Cell_Umap_add2$Row.names 
metadata_Cell_Umap_add2<-metadata_Cell_Umap_add2[,-1]
head(metadata_Cell_Umap_add2)
table(metadata_Cell_Umap_add2$tech)
write.table(as.data.frame(metadata_Cell_Umap_add2), "/mnt/data/chenwei/gongchen/manuscript/manu_table/UMAP_coordinate_and_annotation_information_for_coembedded_landscape.txt",quote = F, row.names = T,sep = "\t")

######
#plot
#####
p_placenta1<-DimPlot(integrated, reduction="umap", group.by = "tech", split.by = "tech",cols = ppCor_all2)
p_placenta2<-DimPlot(integrated, reduction="umap",group.by = "tech",cols = ppCor_all2)
p_placenta3<-DimPlot(integrated, reduction="umap",group.by = "final_anno",split.by = "tech",label = TRUE,cols = ppCor_all2)
p_placenta4<-DimPlot(integrated, reduction="umap",group.by = "re_annotation_TFac",label = TRUE,cols = ppCor_all2)
p_placenta5<-DimPlot(integrated, reduction="umap",group.by = "re_annotation_TFac", split.by = "tech",label = TRUE,cols = ppCor_all2)
p_placenta6<-DimPlot(integrated, reduction="umap", group.by = "orig.ident", split.by = "tech",label = TRUE,cols = ppCor_all2)
p_placenta7<-DimPlot(integrated, group.by = "orig.ident",cols = ppCor_all2, label = TRUE)  + ggtitle("orig.ident")

plot_name<-"three_vs_three_cobatch_GEXs_intergration"
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",plot_name,"_tech_umap_tech_split.pdf"),p_placenta1, width=13, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",plot_name,"_tech_umap.pdf"),p_placenta2, width=9, height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",plot_name,"_final_anno_umap_tech_split.pdf"),p_placenta3, width=13, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",plot_name,"_re_annotation_TFac_umap.pdf"),p_placenta4, width=9, height=8)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",plot_name,"_re_annotation_TFac_umap_tech_split.pdf"),p_placenta5, width=13, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",plot_name,"_orig.ident_umap_tech_split.pdf"),p_placenta6, width=13, height=6)
ggsave(paste0("/mnt/data/chenwei/gongchen/chip_seq_data/gongcheng_data/intergrated_result/",plot_name,"_orig.ident_TFac_umap.pdf"),p_placenta7, width=9, height=8)
