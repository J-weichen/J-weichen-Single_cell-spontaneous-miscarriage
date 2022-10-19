rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(Seurat)
#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
library(future)
plan(strategy = "multicore", workers = 50)

#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac,levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

DefaultAssay(target_final) <- "RNA"
# Normalize RNA data for visualization purposes
target_final <- NormalizeData(target_final, verbose = FALSE)

#1. sub cell type DEGs call 
Idents(object = target_final) <- "re_annotation_TFac" #'merge_cluster'
levels(x = target_final)
Allcluster.markers_all <- FindAllMarkers(object = target_final,min.pct = 0.25)
Allcluster.markers_all$cluster <- factor(x =Allcluster.markers_all$cluster,levels = as.character(levels(x = target_final)))
Allcluster.markers_all_default0.25<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=0.25,]
Allcluster.markers_all_FC1.5<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=log(1.5),]
dim(Allcluster.markers_all_default0.25);dim(Allcluster.markers_all_FC1.5)
DEGs<-c(list(Allcluster.markers_all_default0.25),list(Allcluster.markers_all_FC1.5))
DEGs_name<-c("default0.25","Foldchange1.5")
str(DEGs)
cell_name<-"All_final_subgroup"
for ( filter_number in c(1:2)){
  #filter_number=2#test
  Allcluster.markers_all<- data.frame(DEGs[filter_number])
  type<-DEGs_name[filter_number]
  print(type)
  #selected cell
  Allcluster.markers_all %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  Allcluster.markers_pos <- Allcluster.markers_all[Allcluster.markers_all$avg_logFC>0,]
  Allcluster.markers_pos_01<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.1,]
  Allcluster.markers_pos_005<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.05,]
  Allcluster.markers_pos_001<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.01,]
  dim(Allcluster.markers_pos_001)
  dim(Allcluster.markers_pos_005)
  top10_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top30_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  top30_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  
  Allcluster.markers_neg <- Allcluster.markers_all[Allcluster.markers_all$avg_logFC<0,]
  Allcluster.markers_neg_01<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.1,]
  Allcluster.markers_neg_005<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.05,]
  Allcluster.markers_neg_001<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.01,]
  dim(Allcluster.markers_neg_001)
  dim(Allcluster.markers_neg_005)
  top10_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top30_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  top30_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  
  write.table(Allcluster.markers_all, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_All_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_pos_genes.xls",sep=""),sep="\t", row.names = F)
  write.table(Allcluster.markers_pos_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_pos_01_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_pos_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_pos_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_top10_pos_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_top10_pos_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_pos_001,file =paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_top30_pos_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_top30_pos_005_genes.xls",sep=""),sep="\t",row.names = F)
  
  write.table(Allcluster.markers_neg, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_neg_genes.xls",sep=""), sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_neg_01_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_neg_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_neg_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_top10_neg_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_top10_neg_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_top30_neg_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/",cell_name,"_",type,"_top30_neg_005_genes.xls",sep=""),sep="\t",row.names = F)
  #mm<-read.table("/home/chenwei/10x_data/191125-merge_nonormalization_add/Stromal_Cell_pos_01_genes.xls",head = T,sep = ",")
  #diff_gene <- subset(res3, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1
}

#2. main cell type DEGs call 
Idents(object = target_final) <- "re_anno_TFac_major" #'merge_cluster'
levels(x = target_final)
Allcluster.markers_all <- FindAllMarkers(object = target_final,min.pct = 0.25)
Allcluster.markers_all$cluster <- factor(x =Allcluster.markers_all$cluster,levels = as.character(levels(x = target_final)))
Allcluster.markers_all_default0.25<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=0.25,]
Allcluster.markers_all_FC1.5<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=log(1.5),]
dim(Allcluster.markers_all_default0.25);dim(Allcluster.markers_all_FC1.5)
DEGs<-c(list(Allcluster.markers_all_default0.25),list(Allcluster.markers_all_FC1.5))
DEGs_name<-c("default0.25","Foldchange1.5")
str(DEGs)
cell_name<-"All_final_major_group"
for ( filter_number in c(1:2)){
  #filter_number=2#test
  Allcluster.markers_all<- data.frame(DEGs[filter_number])
  type<-DEGs_name[filter_number]
  print(type)
  #filter_number=2#test
  Allcluster.markers_all<- data.frame(DEGs[filter_number])
  type<-DEGs_name[filter_number]
  print(type)
  #selected cell
  Allcluster.markers_all %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  Allcluster.markers_pos <- Allcluster.markers_all[Allcluster.markers_all$avg_logFC>0,]
  Allcluster.markers_pos_01<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.1,]
  Allcluster.markers_pos_005<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.05,]
  Allcluster.markers_pos_001<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.01,]
  dim(Allcluster.markers_pos_001)
  dim(Allcluster.markers_pos_005)
  top10_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top30_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  top30_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  
  Allcluster.markers_neg <- Allcluster.markers_all[Allcluster.markers_all$avg_logFC<0,]
  Allcluster.markers_neg_01<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.1,]
  Allcluster.markers_neg_005<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.05,]
  Allcluster.markers_neg_001<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.01,]
  dim(Allcluster.markers_neg_001)
  dim(Allcluster.markers_neg_005)
  top10_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top30_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  top30_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  
  write.table(Allcluster.markers_all, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_All_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_pos_genes.xls",sep=""),sep="\t", row.names = F)
  write.table(Allcluster.markers_pos_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_pos_01_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_pos_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_pos_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_top10_pos_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_top10_pos_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_pos_001,file =paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_top30_pos_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_top30_pos_005_genes.xls",sep=""),sep="\t",row.names = F)
  
  write.table(Allcluster.markers_neg, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_neg_genes.xls",sep=""), sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_neg_01_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_neg_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(Allcluster.markers_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_neg_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_top10_neg_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top10_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_top10_neg_005_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_top30_neg_001_genes.xls",sep=""),sep="\t",row.names = F)
  write.table(top30_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/major_group/",cell_name,"_",type,"_top30_neg_005_genes.xls",sep=""),sep="\t",row.names = F)
  #mm<-read.table("/home/chenwei/10x_data/191125-merge_nonormalization_add/Stromal_Cell_pos_01_genes.xls",head = T,sep = ",")
  #diff_gene <- subset(res3, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1
}

##3. for each subtype cluster in  Trophoblast_KRT7 
cell_name<-"Trophoblast_KRT7"
Trophoblast <-subset(x = target_final, subset = re_anno_TFac_major == "Trophoblast_KRT7")
Idents(object = Trophoblast) <- "re_annotation_TFac" #'merge_cluster'
levels(x = Trophoblast)

Allcluster.markers_all <- FindAllMarkers(object = Trophoblast,min.pct = 0.25)
Allcluster.markers_all$cluster <- factor(x =Allcluster.markers_all$cluster,levels = as.character(levels(x = Trophoblast)))
Allcluster.markers_all_default0.25<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=0.25,]
Allcluster.markers_all_FC1.5<-Allcluster.markers_all[abs(Allcluster.markers_all$avg_logFC)>=log(1.5),]
dim(Allcluster.markers_all_default0.25);dim(Allcluster.markers_all_FC1.5)
DEGs<-c(list(Allcluster.markers_all_default0.25),list(Allcluster.markers_all_FC1.5))
DEGs_name<-c("default0.25","Foldchange1.5")
str(DEGs)

for ( filter_number in c(1:2)){
    #filter_number=2#test
    Allcluster.markers_all<- data.frame(DEGs[filter_number])
    type<-DEGs_name[filter_number]
    print(type)
    #selected cell
    Allcluster.markers_all %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
    Allcluster.markers_pos <- Allcluster.markers_all[Allcluster.markers_all$avg_logFC>0,]
    Allcluster.markers_pos_01<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.1,]
    Allcluster.markers_pos_005<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.05,]
    Allcluster.markers_pos_001<-Allcluster.markers_pos[Allcluster.markers_pos$p_val_adj<0.01,]
    dim(Allcluster.markers_pos_001)
    dim(Allcluster.markers_pos_005)
    top10_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    top10_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    top30_pos_001 <- Allcluster.markers_pos_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
    top30_pos_005 <- Allcluster.markers_pos_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
    
    Allcluster.markers_neg <- Allcluster.markers_all[Allcluster.markers_all$avg_logFC<0,]
    Allcluster.markers_neg_01<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.1,]
    Allcluster.markers_neg_005<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.05,]
    Allcluster.markers_neg_001<-Allcluster.markers_neg[Allcluster.markers_neg$p_val_adj<0.01,]
    dim(Allcluster.markers_neg_001)
    dim(Allcluster.markers_neg_005)
    top10_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    top10_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    top30_neg_001 <- Allcluster.markers_neg_001 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
    top30_neg_005 <- Allcluster.markers_neg_005 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
    
    write.table(Allcluster.markers_all, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_All_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(Allcluster.markers_pos, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_pos_genes.xls",sep=""),sep="\t", row.names = F)
    write.table(Allcluster.markers_pos_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_pos_01_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(Allcluster.markers_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_pos_001_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(Allcluster.markers_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_pos_005_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(top10_pos_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_top10_pos_001_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(top10_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_top10_pos_005_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(top30_pos_001,file =paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_top30_pos_001_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(top30_pos_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_top30_pos_005_genes.xls",sep=""),sep="\t",row.names = F)
    
    write.table(Allcluster.markers_neg, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_neg_genes.xls",sep=""), sep="\t",row.names = F)
    write.table(Allcluster.markers_neg_01,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_neg_01_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(Allcluster.markers_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_neg_001_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(Allcluster.markers_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_neg_005_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(top10_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_top10_neg_001_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(top10_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_top10_neg_005_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(top30_neg_001,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_top30_neg_001_genes.xls",sep=""),sep="\t",row.names = F)
    write.table(top30_neg_005,file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_1/inter_celltype_DEGs/Trophoblast/",cell_name,"_",type,"_top30_neg_005_genes.xls",sep=""),sep="\t",row.names = F)
    
  }  
