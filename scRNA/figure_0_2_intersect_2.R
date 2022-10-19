rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib")
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

require("Matrix")
library(data.table)
library(Seurat)
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
saveRDS(seurat_data_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713.rds")

dtlist <- seurat_data_list
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 3000) #2976 
dtlist <- PrepSCTIntegration(object.list = dtlist,anchor.features = intfts)
anchors <- FindIntegrationAnchors(object.list = dtlist, normalization.method = "SCT", anchor.features = intfts)
saveRDS(anchors, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713_anchor.rds")
seurat_data_list<-anchors<-dtlist<-intfts<-0

########################
### Integrate
########################
anchors<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713_anchor.rds")
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated <- RunPCA(integrated)
ElbowPlot(object = integrated)
VizDimLoadings(object = integrated, dims = 1:2, reduction = "pca")
integrated <- RunUMAP(integrated, dims = 1:25)
integrated<- RunTSNE(integrated, dims = 1:25)
saveRDS(integrated, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0713_integrated_seurat.rds")

##for no fileter
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
