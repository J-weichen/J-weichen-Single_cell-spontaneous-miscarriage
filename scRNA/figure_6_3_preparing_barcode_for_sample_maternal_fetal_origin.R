rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(data.table)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(grid)
library(gridExtra)
library(viridis)
library(ggpubr)
library(purrr)
library(cowplot)
grid.newpage()
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:46))]


#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
##get metadata
meta_information<-target_final@meta.data
head(meta_information);tail(meta_information)
meta_tag<-meta_information[,c("sample","sample_code","raw_id","re_annotation_TFac")]
meta_tag$BARCODE<-as.character(unlist(lapply(strsplit(rownames(meta_tag),"_"), function(x) x[1])))
meta_tag$BARCODE_new<-paste0("CB:Z:",meta_tag$BARCODE)

meta_tag$re_annotation_TFac<-as.character(meta_tag$re_annotation_TFac)
names(table(meta_tag$re_annotation_TFac))
meta_tag$MF_type<-ifelse(meta_tag$re_annotation_TFac %in% c("Mycs","STCs","NKs","Ts","Bs"),"maternal","fetal")
head(meta_tag)
sample_name<-unique(as.character(meta_tag$sample))
for (population in sample_name) {
  #population<-"N1" ##test line
  print(as.character(population))
  sub_meta<-meta_tag[which(meta_tag$sample ==population),]
  
  sub_meta_maternal<-meta_tag[which(sub_meta$MF_type =="maternal"),]
  sub_meta_fetal<-meta_tag[which(sub_meta$MF_type =="fetal"),]
  write.table(as.data.frame(sub_meta_maternal[,"BARCODE_new"]), file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure6/GEX_Abortion_SNP_bam/barcode_split/",population,"_maternal_barcode.txt"),quote=F, row.names=F, col.names=F,sep="\t") 
  write.table(as.data.frame(sub_meta_fetal[,"BARCODE_new"]), file=paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure6/GEX_Abortion_SNP_bam/barcode_split/",population,"_fetal_barcode.txt"),quote=F, row.names=F, col.names=F,sep="\t") 
}  
table(meta_tag$sample,meta_tag$MF_type)
