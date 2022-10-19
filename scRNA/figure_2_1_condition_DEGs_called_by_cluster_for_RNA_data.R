rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(Seurat)
library(UpSetR)

library(future)
plan(strategy = "multicore", workers = 24)
options(future.globals.maxSize = 1500 * 1024^12)

#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)


DefaultAssay(target_final) <- "RNA"
# Normalize RNA data for visualization purposes
target_final <- NormalizeData(target_final, verbose = FALSE)

# set idents from a value in object metadata
colnames(x = target_final[[]])
plot<-as.data.frame(table(target_final$re_annotation_TFac))
cell_type<-as.character(plot[which(plot$Freq>=3),]$Var1)
cell_type

#for major group
target_final$re_anno_TFac_major_Treat <- paste(target_final$re_anno_TFac_major, target_final$Treat, sep = "_")
Idents(object = target_final) <- "re_anno_TFac_major_Treat"
levels(x = target_final)
length(unique(target_final$re_anno_TFac_major_Treat))

compare_name<-as.character(unique(target_final$re_anno_TFac_major))
length(compare_name)
compare_name2<-compare_name
compare_name2<-compare_name[which(!(compare_name %in% c("Mast_Cell_MS4A2","Erythrocyte_HBA1")))]
#Using default 参数
for (tag_name in compare_name2){
  #  tag_name <-"T_Cell_CD3D"
  print(paste0("RUN ",tag_name))
  Age_fault<- FindMarkers(object = target_final, ident.1 = paste(tag_name,"_Abortion",sep=""), ident.2 =paste(tag_name,"_CTRL",sep=""),verbose = FALSE,min.pct = 0.25, logfc.threshold = 0.25,min.cells.group = 3)
  Age_fault$cluster<-rep(tag_name,nrow(Age_fault))
  Age_fault$gene<-rownames(Age_fault)
  ##=============
  Age_fault %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  Allcluster.markers_pos <- Age_fault[Age_fault$avg_logFC>0,]
  Allcluster.markers_pos_01<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.1,]
  Allcluster.markers_pos_005<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.05,]
  Allcluster.markers_pos_001<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.01,]
  dim(Allcluster.markers_pos_001)
  dim(Allcluster.markers_pos_005)
  top10_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top30_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  top30_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  
  Allcluster.markers_neg <- Age_fault[Age_fault$avg_logFC<0,]
  Allcluster.markers_neg_01<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.1,]
  Allcluster.markers_neg_005<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.05,]
  Allcluster.markers_neg_001<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.01,]
  dim(Allcluster.markers_neg_001)
  dim(Allcluster.markers_neg_005)
  top10_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top30_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  top30_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  #==============
  
  write.table(Age_fault, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_All_Abortion_CTRL_RNA_data_defualt_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_pos_Abortion_CTRL_RNA_data_defualt_genes.xls",sep=""),sep="\t", row.names = F)
  write.table(Allcluster.markers_pos_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_pos_Abortion_CTRL_RNA_data_defualt_01_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_pos_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_pos_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_top10_pos_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_top10_pos_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_pos_001,file =paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_top30_pos_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_top30_pos_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  
  write.table(Allcluster.markers_neg, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_neg_Abortion_CTRL_RNA_data_defualt_genes.xls",sep=""), sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_neg_Abortion_CTRL_RNA_data_defualt_01_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_neg_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_neg_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_top10_neg_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_top10_neg_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_top30_neg_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/",tag_name,"_top30_neg_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
}


#对差基因的不同组关系进行比较
#For default 
#read DEGs between PPH list
files<- list.files("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/Pos_genes_01q")
#list.files命令将input文件夹下所有文件名输入
data_list<-list();cnane<-c()
for ( f in files[1:length(files)]){
  FF<-paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/Pos_genes_01q/",f,sep="")
  column <- read.delim(FF,header = T,sep = "\t")[,c("cluster","gene")]
  if (is.na(column[1,]$cluster)) {
    next
  }
  data_list<-c(data_list,list(as.character(column$gene)))
  print(unique(as.character(column$cluster)))
  cnane<-c(cnane,paste(unique(as.character(column$cluster)),"_up",sep=""))
}
length(data_list)
length(cnane)
names(data_list)<-cnane

files<- list.files("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/Neg_genes_01q")
data_list2<-list();cnane2<-c()
for ( f in files[1:length(files)]){
  FF<-paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/major_group/Neg_genes_01q/",f,sep="")
  column <- read.delim(FF,header = T,sep = "\t")[,c("cluster","gene")]
  if (is.na(column[1,]$cluster)) {
    next
  }
  data_list2<-c(data_list2,list(as.character(column$gene)))
  print(unique(as.character(column$cluster)))
  cnane2<-c(cnane2,paste(unique(as.character(column$cluster)),"_down",sep=""))
}
length(data_list2)
length(cnane2)
names(data_list2)<-cnane2

conflict_prefer("which", "Matrix")

listinput<-c(data_list,data_list2)
upset_dataframe<-as.data.frame(fromList(listinput))
dim(upset_dataframe)
upset_dataframe[1:10,1:10]
upset_dataframe_rowSums<-rowSums(upset_dataframe)
max(upset_dataframe_rowSums)
upset_dataframe_colSums<-colSums(upset_dataframe)
range(upset_dataframe_colSums)# 16 622
upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset_dataframe2<-upset_dataframe
upset_dataframe2$summ<-as.numeric(upset_dataframe_rowSums)
dim(upset_dataframe2)

write.table(upset_dataframe2, file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_major_group_DEG_01_default_relationship.txt",quote=F, row.names=F, col.names=T,sep="\t")

listinput_plot<-list()
cell_type<-names(upset_dataframe)
for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe) #2340   18
max(rowSums(upset_dataframe))#9
range(colSums(upset_dataframe))
#  16 622
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/RNA_data_Upset_plot_for_Abortion_related_DEGs_major_cluster.pdf",width = 10,height = 8)

#list2_q1<-list(query = intersects, params = list("B_Cell_neg", "Endothelial_Cell_neg","Epithelial_Cell_neg","Erythrocyte_neg","Myeloid_Cell_neg", "Stromal_Cell_neg","Trophoblast_neg"), color = "red", active = T,query.name = "common down_PPH_DEG")
#list2_q2<-list(query = intersects,params = list("B_Cell_pos","Endothelial_Cell_pos","Epithelial_Cell_pos","Erythrocyte_pos","Myeloid_Cell_pos","Stromal_Cell_pos","Trophoblast_pos"), color = "orange", active = T,query.name = "common up_PPH_DEG")
upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)),
      nintersects = 500,
      #      sets=c("dS1_neg","dS1_pos","dP1_neg","dP1_pos","dS2_neg","dS2_pos","SCT_neg","SCT_pos","dM3_pos","dP2_pos","dS3_pos","Epi2_pos","EVT_pos" ,
      #             "fFB1_pos","HB_pos","VCT_pos"),
      keep.order = TRUE,
      query.legend = "top",
      sets.bar.color = "purple",
      #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
      shade.color="pink",
      #      matrix.color="purple",
      order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
      point.size = 3,line.size = 1.3,
      mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
      mb.ratio = c(0.60, 0.40),
      # text.scale = c(1.5, 1.5, 1.2,1.5,1.5,1),
      show.numbers = 'yes'
      #  queries = list(list2_q1,list2_q2)
)

dev.off()

###########for subtypes
Idents(object = target_final) <- 're_annotation_TFac'
target_final$re_annotation_TFac_Treat <- paste(Idents(object = target_final), target_final$Treat, sep = "_")
levels(x = target_final)

#For genes change in different Age for cells of the same type.
Idents(object = target_final) <- "re_annotation_TFac_Treat"
levels(x = target_final)
length(unique(target_final$re_annotation_TFac_Treat))

compare_name<-as.character(unique(target_final$re_annotation_TFac))
#fine_test<-c(paste(compare_name,"_Abortion",sep=""),paste(compare_name,"_CTRL",sep=""))
#nocompare_name<-setdiff(fine_test,as.character(unique(target_final$re_annotation_TFac_Treat)))
#nocompare_name <- as.vector(unlist(sapply(strsplit(nocompare_name, split="_"), "[[", 1) ))
#nocompare_name
#compare_name<-compare_name[-match(nocompare_name,compare_name)]
length(compare_name)

compare_name2<-compare_name[which(!(compare_name %in% c("Masts","Ery")))]
#compare_name2<-compare_name
#
#Using default 参数
for (tag_name in compare_name2){
  #  tag_name <-"EVTs_2"
  print(paste0("RUN ",tag_name))
  Age_fault<- FindMarkers(object = target_final, ident.1 = paste(tag_name,"_Abortion",sep=""), ident.2 =paste(tag_name,"_CTRL",sep=""),verbose = FALSE,min.pct = 0.25, logfc.threshold = 0.25,min.cells.group = 3)
  Age_fault$cluster<-rep(tag_name,nrow(Age_fault))
  Age_fault$gene<-rownames(Age_fault)
  ##=============
  Age_fault %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  Allcluster.markers_pos <- Age_fault[Age_fault$avg_logFC>0,]
  Allcluster.markers_pos_01<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.1,]
  Allcluster.markers_pos_005<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.05,]
  Allcluster.markers_pos_001<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.01,]
  dim(Allcluster.markers_pos_001)
  dim(Allcluster.markers_pos_005)
  dim(Allcluster.markers_pos_01)
  
  top10_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top30_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  top30_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  
  Allcluster.markers_neg <- Age_fault[Age_fault$avg_logFC<0,]
  Allcluster.markers_neg_01<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.1,]
  Allcluster.markers_neg_005<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.05,]
  Allcluster.markers_neg_001<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.01,]
  dim(Allcluster.markers_neg_001)
  dim(Allcluster.markers_neg_005)
  top10_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top30_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  top30_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  #==============
  
  write.table(Age_fault, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_All_Abortion_CTRL_RNA_data_defualt_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_pos_Abortion_CTRL_RNA_data_defualt_genes.xls",sep=""),sep="\t", row.names = F)
  write.table(Allcluster.markers_pos_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_pos_Abortion_CTRL_RNA_data_defualt_01_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_pos_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_pos_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_top10_pos_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_top10_pos_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_pos_001,file =paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_top30_pos_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_top30_pos_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  
  write.table(Allcluster.markers_neg, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_neg_Abortion_CTRL_RNA_data_defualt_genes.xls",sep=""), sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_neg_Abortion_CTRL_RNA_data_defualt_01_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_neg_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_neg_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_top10_neg_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_top10_neg_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_top30_neg_Abortion_CTRL_RNA_data_defualt_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/",tag_name,"_top30_neg_Abortion_CTRL_RNA_data_defualt_005_genes.xls",sep=""),sep="\t",row.names = F)
}



#对差基因的不同组关系进行比较
#For default 
files<- list.files("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/Pos_genes_01q")
#list.files命令将input文件夹下所有文件名输入
data_list<-list();cnane<-c()
for ( f in files[1:length(files)]){
  FF<-paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/Pos_genes_01q/",f,sep="")
  column <- read.delim(FF,header = T,sep = "\t")[,c("cluster","gene")]
  if (is.na(column[1,]$cluster)) {
    next
  }
  data_list<-c(data_list,list(as.character(column$gene)))
  print(unique(as.character(column$cluster)))
  cnane<-c(cnane,paste(unique(as.character(column$cluster)),"_up",sep=""))
}
length(data_list)#15
length(cnane)#15
names(data_list)<-cnane

files<- list.files("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/Neg_genes_01q")
data_list2<-list();cnane2<-c()
for ( f in files[1:length(files)]){
  FF<-paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/RNA_data/Neg_genes_01q/",f,sep="")
  column <- read.delim(FF,header = T,sep = "\t")[,c("cluster","gene")]
  if (is.na(column[1,]$cluster)) {
    next
  }
  data_list2<-c(data_list2,list(as.character(column$gene)))
  print(unique(as.character(column$cluster)))
  cnane2<-c(cnane2,paste(unique(as.character(column$cluster)),"_down",sep=""))
}
length(data_list2)#17
length(cnane2)#17
names(data_list2)<-cnane2

listinput<-c(data_list,data_list2)
upset_dataframe<-as.data.frame(fromList(listinput))
dim(upset_dataframe)
upset_dataframe[1:10,1:10]
upset_dataframe_rowSums<-rowSums(upset_dataframe)
max(upset_dataframe_rowSums)#16
upset_dataframe_colSums<-colSums(upset_dataframe)
range(upset_dataframe_colSums)#16 660
upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset_dataframe2<-upset_dataframe
upset_dataframe2$summ<-as.numeric(upset_dataframe_rowSums)
dim(upset_dataframe2)# 2755   33

write.table(upset_dataframe2, file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_major_subgroup_DEG_01_default_relationship.txt",quote=F, row.names=F, col.names=T,sep="\t")

listinput_plot<-list()
cell_type<-names(upset_dataframe)
for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe) #2755   32
max(rowSums(upset_dataframe))#16
range(colSums(upset_dataframe))
# 16 660
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/RNA_data_Upset_plot_for_Abortion_related_DEGs_subgroup_cluster.pdf",width = 10,height = 8)

#list2_q1<-list(query = intersects, params = list("B_Cell_neg", "Endothelial_Cell_neg","Epithelial_Cell_neg","Erythrocyte_neg","Myeloid_Cell_neg", "Stromal_Cell_neg","Trophoblast_neg"), color = "red", active = T,query.name = "common down_PPH_DEG")
#list2_q2<-list(query = intersects,params = list("B_Cell_pos","Endothelial_Cell_pos","Epithelial_Cell_pos","Erythrocyte_pos","Myeloid_Cell_pos","Stromal_Cell_pos","Trophoblast_pos"), color = "orange", active = T,query.name = "common up_PPH_DEG")
upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)),
      nintersects = 500,
      #      sets=c("dS1_neg","dS1_pos","dP1_neg","dP1_pos","dS2_neg","dS2_pos","SCT_neg","SCT_pos","dM3_pos","dP2_pos","dS3_pos","Epi2_pos","EVT_pos" ,
      #             "fFB1_pos","HB_pos","VCT_pos"),
      keep.order = TRUE,
      query.legend = "top",
      sets.bar.color = "purple",
      #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
      shade.color="pink",
      #      matrix.color="purple",
      order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
      point.size = 3,line.size = 1.3,
      mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
      mb.ratio = c(0.60, 0.40),
      # text.scale = c(1.5, 1.5, 1.2,1.5,1.5,1),
      show.numbers = 'yes'
      #  queries = list(list2_q1,list2_q2)
)
dev.off()

##以下未跑
#plot for whole Cells
#target_final$re_annotation_TFac
#B_Cell_CD79A：5779      Decidual_stromal_cells_DKK1_ACTA2 ：12327
#Endometrial_Epithelial_Cell_EPCAM_PAEP：137   Endothelial_Cell_PECAM1 ：2635
#Fetal_Stromal_cells_Endothelial_cells_DLK1： 5519   Megakaryocytes_PPBP：211
#Myeloid_Cell_AIF1：40350    T_NK_raw_Erythrocyte_Mast_Cell：82775
#Trophoblast_mix_EECs_KRT7：6302

levels=c( "CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2","STBs_3","TroB_VIM",
          "F_STCs","F_FBs_1","F_FBs_2","F_FBs_3","F_FBs_VECs",
          "D_STCs","D_FBs_1","D_FBs_2","D_FBs_3","MFs","PVs",
          "D_LECs","D_VECs","F_VECs", "MKCs_1","MKCs_2","Epis",
          "GMPs","DCs","Mac_1","Mac_2","Mac_3","HBs","nc_Mon_1","nc_Mon_2","B_c_Mon","c_Mon_2","c_Mon_3",
          "B_Bs","Trm_Bs","MALT_Bs","Aty_Bs",
          "D_NKs","NB_NKs","D_Cyt_Ts","D_Mem_CD4_Ts","F_Nai_CD4_Ts",
          "B_NKs","B_Cyt_Ts","B_NKTs","B_Mem_CD4_Ts", "B_Nai_CD4_Ts",        
          "F_AIF1_NK","F_AIF1_T","Mast","Ery"),ordered=TRUE)

names(data_list)
names(data_list2)
listinput<-c(data_list,data_list2)

#For Trophoblast
listinput_plot<-list()
cell_type<-c("CTBs_1_up","CTBs_2_up","EVTs_1_up","EVTs_2_up","EVTs_3_up","STBs_1_up","STBs_2_up","STBs_3_up",
             "CTBs_1_down","CTBs_2_down","EVTs_1_down","EVTs_2_down","EVTs_3_down","STBs_1_down","STBs_2_down","STBs_3_down")
for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
Trophoblast_plot<-upset(upset_dataframe, main.bar.color = "black",nsets = 8, nintersects = 200,
                        sets=c(colnames(upset_dataframe)),
                        keep.order = TRUE,
                        query.legend = "top",
                        sets.bar.color = "purple",
                        #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                        shade.color="pink",
                        #      matrix.color="purple",
                        order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                        #      point.size = 3,line.size = 1,
                        mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                        mb.ratio = c(0.60, 0.40),
                        #      text.scale = c(2, 2, 1.2,2,2,1),
                        show.numbers = 'yes')
#      queries = list(list2_q1,list2_q2))

#For Fetal_Stromal_Cell
listinput_plot<-list()
cell_type<-c("F_STCs_up","F_FBs_1_up","F_FBs_2_up","F_FBs_3_up",
             "F_STCs_down","F_FBs_1_down","F_FBs_2_down","F_FBs_3_down")


for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
Fetal_Stromal_Cell_plot<-upset(upset_dataframe, main.bar.color = "black",
                               #nsets = 8,
                               nintersects = 400,
                               sets=c(colnames(upset_dataframe)),
                               keep.order = TRUE,
                               query.legend = "top",
                               sets.bar.color = "purple",
                               #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                               shade.color="pink",
                               #      matrix.color="purple",
                               order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                               #      point.size = 3,line.size = 1,
                               mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                               mb.ratio = c(0.60, 0.40),
                               #      text.scale = c(2, 2, 1.2,2,2,1),
                               show.numbers = 'yes')
#For Maternal_Stromal_Cell
listinput_plot<-list()
cell_type<-c("D_STCs_up","D_FBs_1_up","D_FBs_2_up","D_FBs_3_up","MFs_up","PVs_up",
             "D_STCs_down","D_FBs_1_down","D_FBs_2_down","D_FBs_3_down","MFs_down","PVs_down")

for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
Maternal_Stromal_Cell_plot<-upset(upset_dataframe, main.bar.color = "black",
                                  #nsets = 8,
                                  nintersects = 400,
                                  sets=c(colnames(upset_dataframe)),
                                  keep.order = TRUE,
                                  query.legend = "top",
                                  sets.bar.color = "purple",
                                  #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                                  shade.color="pink",
                                  #      matrix.color="purple",
                                  order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                                  #      point.size = 3,line.size = 1,
                                  mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                                  mb.ratio = c(0.60, 0.40),
                                  #      text.scale = c(2, 2, 1.2,2,2,1),
                                  show.numbers = 'yes')
#For Stromal_Cell
listinput_plot<-list()
cell_type<-c("F_STCs_up","F_FBs_1_up","F_FBs_2_up","F_FBs_3_up","D_STCs_up","D_FBs_1_up","D_FBs_2_up","D_FBs_3_up","MFs_up","PVs_up",
             "F_STCs_down","F_FBs_1_down","F_FBs_2_down","F_FBs_3_down","D_STCs_down","D_FBs_1_down","D_FBs_2_down","D_FBs_3_down","MFs_down","PVs_down")


for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
Fetal_maternal_Stromal_Cell_plot<-upset(upset_dataframe, main.bar.color = "black",
                                        #nsets = 8,
                                        nintersects = 400,
                                        sets=c(colnames(upset_dataframe)),
                                        keep.order = TRUE,
                                        query.legend = "top",
                                        sets.bar.color = "purple",
                                        #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                                        shade.color="pink",
                                        #      matrix.color="purple",
                                        order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                                        #      point.size = 3,line.size = 1,
                                        mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                                        mb.ratio = c(0.60, 0.40),
                                        #      text.scale = c(2, 2, 1.2,2,2,1),
                                        show.numbers = 'yes')
#      queries = list(list2_q1,list2_q2))

#For Endothelial_Cell
listinput_plot<-list()
cell_type<-c("D_LECs_up","D_VECs_up","F_VECs_up","D_LECs_down","D_VECs_down","F_VECs_down")

for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
Endothelial_Cell_plot<-upset(upset_dataframe, main.bar.color = "black",
                             #nsets = 8,
                             nintersects = 400,
                             sets=c(colnames(upset_dataframe)),
                             keep.order = TRUE,
                             query.legend = "top",
                             sets.bar.color = "purple",
                             #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                             shade.color="pink",
                             #      matrix.color="purple",
                             order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                             #      point.size = 3,line.size = 1,
                             mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                             mb.ratio = c(0.60, 0.40),
                             #      text.scale = c(2, 2, 1.2,2,2,1),
                             show.numbers = 'yes')

#For Myeloid_Cell
listinput_plot<-list()
cell_type<-c("GMPs_up","DCs_up","Mac_1_up","Mac_2_up","Mac_3_up","HBs_up","nc_Mon_1_up","nc_Mon_2_up","B_c_Mon_up","c_Mon_2_up","c_Mon_3_up",
             "GMPs_down","DCs_down","Mac_1_down","Mac_2_down","Mac_3_down","HBs_down","nc_Mon_1_down","nc_Mon_2_down","B_c_Mon_down","c_Mon_2_down","c_Mon_3_down")

for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
Myeloid_Cell_plot<-upset(upset_dataframe, main.bar.color = "black",
                         #nsets = 8,
                         nintersects = 500,
                         sets=c(colnames(upset_dataframe)),
                         keep.order = TRUE,
                         query.legend = "top",
                         sets.bar.color = "purple",
                         #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                         shade.color="pink",
                         #      matrix.color="purple",
                         order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                         #      point.size = 3,line.size = 1,
                         mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                         mb.ratio = c(0.60, 0.40),
                         #      text.scale = c(2, 2, 1.2,2,2,1),
                         show.numbers = 'yes')

#For T_NK_Cell
listinput_plot<-list()
cell_type<-c("D_NKs_up","NB_NKs_up","D_Cyt_Ts_up","D_Mem_CD4_Ts_up","F_Nai_CD4_Ts_up", "B_NKs_up","B_Cyt_Ts_up","B_NKTs_up","B_Mem_CD4_Ts_up", "B_Nai_CD4_Ts_up",
             "D_NKs_down","NB_NKs_down","D_Cyt_Ts_down","D_Mem_CD4_Ts_down","F_Nai_CD4_Ts_down","B_NKs_down","B_Cyt_Ts_down","B_NKTs_down","B_Mem_CD4_Ts_down", "B_Nai_CD4_Ts_down")

for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
T_NK_Cell_plot<-upset(upset_dataframe, main.bar.color = "black",
                      #nsets = 8,
                      nintersects = 500,
                      sets=c(colnames(upset_dataframe)),
                      keep.order = TRUE,
                      query.legend = "top",
                      sets.bar.color = "purple",
                      #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                      shade.color="pink",
                      #      matrix.color="purple",
                      order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                      #      point.size = 3,line.size = 1,
                      mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                      mb.ratio = c(0.60, 0.40),
                      #      text.scale = c(2, 2, 1.2,2,2,1),
                      show.numbers = 'yes')

#For T_Cell
listinput_plot<-list()
cell_type<-c("D_Cyt_Ts_up","D_Mem_CD4_Ts_up","F_Nai_CD4_Ts_up","B_Cyt_Ts_up","B_NKTs_up","B_Mem_CD4_Ts_up", "B_Nai_CD4_Ts_up",
             "D_Cyt_Ts_down","D_Mem_CD4_Ts_down","F_Nai_CD4_Ts_down","B_Cyt_Ts_down","B_NKTs_down","B_Mem_CD4_Ts_down","B_Nai_CD4_Ts_down")

for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
T_Cell_plot<-upset(upset_dataframe, main.bar.color = "black",
                   #nsets = 8,
                   nintersects = 500,
                   sets=c(colnames(upset_dataframe)),
                   keep.order = TRUE,
                   query.legend = "top",
                   sets.bar.color = "purple",
                   #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                   shade.color="pink",
                   #      matrix.color="purple",
                   order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                   #      point.size = 3,line.size = 1,
                   mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                   mb.ratio = c(0.60, 0.40),
                   #      text.scale = c(2, 2, 1.2,2,2,1),
                   show.numbers = 'yes')


#For NK_Cell
listinput_plot<-list()
cell_type<-c("D_NKs_up","NB_NKs_up","B_NKs_up","D_NKs_down","NB_NKs_down","B_NKs_down")

for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
NK_Cell_plot<-upset(upset_dataframe, main.bar.color = "black",
                    #nsets = 8,
                    nintersects = 500,
                    sets=c(colnames(upset_dataframe)),
                    keep.order = TRUE,
                    query.legend = "top",
                    sets.bar.color = "purple",
                    #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                    shade.color="pink",
                    #      matrix.color="purple",
                    order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                    #      point.size = 3,line.size = 1,
                    mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                    mb.ratio = c(0.60, 0.40),
                    #      text.scale = c(2, 2, 1.2,2,2,1),
                    show.numbers = 'yes')


#For B_Cell
listinput_plot<-list()
cell_type<-c("B_Bs_up","Trm_Bs_up","MALT_Bs_up","Aty_Bs_up"
             "B_Bs_down","Trm_Bs_down","MALT_Bs_down","Aty_Bs_down")

for (i in cell_type){
  listinput_plot[[i]]<- listinput[[i]]
}
names(listinput_plot)
upset_dataframe<-as.data.frame(fromList(listinput_plot))
which(colSums(upset_dataframe)== 0)
ifelse(length(which(colSums(upset_dataframe)==0))==0,upset_dataframe<-upset_dataframe,upset_dataframe<-upset_dataframe[,-which(colSums(upset_dataframe)== 0)])
head(upset_dataframe)
colnames(upset_dataframe)

dim(upset_dataframe)
max(rowSums(upset_dataframe))
range(colSums(upset_dataframe))
#upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]

upset(upset_dataframe,nsets = length(colnames(upset_dataframe)),order.by = "freq")
#list2_q1<-list(query = intersects, params = list("STBs0_neg","EVTs0_neg","VCTs0_neg","Trobt_Vim+_neg"), color = "red", active = T,query.name = "down_age_DEG in Trophoblast")
#list2_q2<-list(query = intersects,params = list("STBs0_pos","EVTs0_pos","VCTs0_pos","Trobt_Vim+_pos"), color = "orange", active = T,query.name = "up_age_DEG in Trophoblast")
NK_Cell_plot<-upset(upset_dataframe, main.bar.color = "black",
                    #nsets = 8,
                    nintersects = 500,
                    sets=c(colnames(upset_dataframe)),
                    keep.order = TRUE,
                    query.legend = "top",
                    sets.bar.color = "purple",
                    #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
                    shade.color="pink",
                    #      matrix.color="purple",
                    order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
                    #      point.size = 3,line.size = 1,
                    mainbar.y.label = "Gene number Intersections", sets.x.label = "Gene number per subset",
                    mb.ratio = c(0.60, 0.40),
                    #      text.scale = c(2, 2, 1.2,2,2,1),
                    show.numbers = 'yes')
