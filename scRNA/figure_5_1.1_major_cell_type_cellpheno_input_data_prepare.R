rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
#dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(Seurat)

#step1 ：read Seurat object and load initial coldata and matrix
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

target_final$re_annotation_TFac_Treat <- paste(target_final$re_annotation_TFac, target_final$Treat, sep = "_")
table(target_final$re_anno_TFac_major)
unique(target_final$re_annotation_TFac_Treat)

DefaultAssay(object = target_final)#"SCT"
DefaultAssay(target_final) <- "SCT"
#ALL_coldata<-target@meta.data
#exprMat0<-GetAssayData(object =target,slot = "counts")
#cellInfo <- ALL_coldata

#for cellphonedb_count
##CellPhoneDB requires data be be normalized but not log-transformed but Seurat LogNormalizes the data.
#The ‘corrected’ UMI counts are stored in pbmc[["SCT"]]@Counts. 
#We store log-normalized versions of these corrected counts in pbmc[["SCT"]]@DaTa, which are very helpful for visualization.
#I would note that in the CellPhoneDB paper, they recommend using normalized data. 
#In my analysis, I performed SCTransform on my data and used "scale.data" to run CellPhoneDB
#You could use either data slot, or scale.data slot after SCTransform.
#The scale.data slot contains the pearson residuals of variable genes, which is also corrected for the confounding effects you put.
#ref1::  https://github.com/satijalab/seurat/issues/2049
#ref2::  https://github.com/Teichlab/cellphonedb/issues/175
#ref3::  https://www.jianshu.com/p/38a9376f5286
table(target$Treat)
#Abortion     CTRL 
# 19648    30340 
mydata_list<-c("Abortion","CTRL")
for ( group_name in mydata_list){
  #group_name<-"CTRL"
  print(as.character(group_name))
  sub_cell<-subset(x = target_final, subset = Treat == group_name)
  #cellphonedb_count<-as.matrix(pbmc3k@assays$RNA@data)
  # cellphonedb_count <- as.matrix(sub_cell@assays$SCT@scale.data)
  count_raw <- target_final@assays$RNA@counts
  cellphonedb_count <- apply(count_raw, 2, function(x) (x/sum(x))*10000)#按照202004  nature protocol 文章进行
  cellphonedb_count <- cellphonedb_count[rowSums(cellphonedb_count)!=0,]
  dim(cellphonedb_count)
  meta_data <- sub_cell@meta.data[,c("BARCODE_new","re_anno_TFac_major")]
  meta_data$BARCODE_new<-rownames(meta_data)
  meta_data <- as.matrix(meta_data)
  meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
  head(meta_data)
  colnames(meta_data)<-c("Cell","cell_type")
  write.table(cellphonedb_count, paste0("//mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/major_celltype/",group_name,"/cellphonedb_count.txt"), sep='\t', quote=F)
  write.table(meta_data, paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/major_celltype/",group_name,"/cellphonedb_meta.txt"), sep='\t', quote=F, row.names=F)
}

##for all cell 
  count_raw <- target_final@assays$RNA@counts
  cellphonedb_count <- apply(count_raw, 2, function(x) (x/sum(x))*10000)#按照202004  nature protocol 文章进行
  cellphonedb_count <- cellphonedb_count[rowSums(cellphonedb_count)!=0,]

  meta_data <- target_final@meta.data[,c("BARCODE_new","re_anno_TFac_major")]
  meta_data$BARCODE_new<-rownames(meta_data)
  meta_data <- as.matrix(meta_data)
  meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
  #head(meta_data)
  colnames(meta_data)<-c("Cell","cell_type")
  write.table(cellphonedb_count, "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/major_celltype/all_cell/cellphonedb_count.txt", sep='\t', quote=F)
  write.table(meta_data,"/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/major_celltype/all_cell/cellphonedb_meta.txt", sep='\t', quote=F, row.names=F)
