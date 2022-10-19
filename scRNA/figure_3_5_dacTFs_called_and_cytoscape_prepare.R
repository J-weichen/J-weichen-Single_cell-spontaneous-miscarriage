rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
options(stringsAsFactors = FALSE)
## Load packages
dyn.load('/usr/local/hdf5/lib/libhdf5_hl.so.100')
require("Matrix")
library(data.table)
library(Seurat)
library(dplyr)
library(data.table)
library(gridExtra)

library(ggplot2)
library(ggrepel)
library(grid)
library(ggsci)
library(cowplot)
library(scales)
library(ggpubr)
library(tidyverse)
library(pheatmap)
library(VennDiagram)
library(futile.logger)

#
library(philentropy)
#tidyverse_conflicts()

library(pbapply)
library(circlize)
library(plyr)
library(dendextend)
library(ComplexHeatmap)
#step0 :set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jAbortion("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(55)
show_col(Cells_col)
show_col(ppCor)


#load GenAge database genes

## DEGs called 
merge_data0<-read.table(file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/DEG_01q_default_merge.txt",header =T,sep="\t") 
head(merge_data0)
merge_data_rm0<-merge_data0[which(!(merge_data0$cluster %in% c("Ery"))),]
sc_Up_gene<-unique(as.character(merge_data_rm0[which(merge_data_rm0$avg_logFC > 0),]$gene))

sc_Down_gene<-unique(as.character(merge_data_rm0[which(merge_data_rm0$avg_logFC < 0),]$gene))
length(sc_Up_gene)# 130; length(sc_Down_gene)# 147


bulk_sc_DEGs_merge <- read.table(file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/bulk_DEGs_and_DEG_01q_default_merge.txt", sep = "\t", header = T)
head(bulk_sc_DEGs_merge);tail(bulk_sc_DEGs_merge);dim(bulk_sc_DEGs_merge)#2143    4
unique(bulk_sc_DEGs_merge$cluster)
bulk_sc_DEGs_up<-unique(bulk_sc_DEGs_merge[which(bulk_sc_DEGs_merge$avg_logFC >0),"gene"])
bulk_sc_DEGs_down<-unique(bulk_sc_DEGs_merge[which(bulk_sc_DEGs_merge$avg_logFC <0),"gene"])
length(bulk_sc_DEGs_up)# 449; 
length(bulk_sc_DEGs_down)# 264

##human TFs and TFs_cofactor 
Homo_sapiens_TF <- read.table("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/addiction_humanTFs_list/Homo_sapiens_TF.txt", sep = "\t", header = T)
head(Homo_sapiens_TF);tail(Homo_sapiens_TF);dim(Homo_sapiens_TF)
length(unique(Homo_sapiens_TF$Symbol))#1665
Homo_sapiens_coTF <- read.table("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/addiction_humanTFs_list/Homo_sapiens_TF_cofactors.txt", sep = "\t", header = T)
head(Homo_sapiens_coTF);tail(Homo_sapiens_coTF);dim(Homo_sapiens_coTF)
length(unique(Homo_sapiens_coTF$Symbol))#1025
Hm_TFs<-unique(Homo_sapiens_TF$Symbol)
Hm_coTFs<-unique(Homo_sapiens_coTF$Symbol)

#for homo genes
human_mouse_homo_gene0 <- read.table("/home/chenwei/10x_data/191125-merge_nonormalization_add/homo_gene_human_mouse/homo_human_mouse_YD_used_gene.txt", sep = "\t", header = T)
head(human_mouse_homo_gene0)
human_mouse_homo_gene_DEGs<-human_mouse_homo_gene0[which(human_mouse_homo_gene0$H_Gene %in% unique(DEGs_merge$gene)),]
dim(human_mouse_homo_gene_DEGs) #699   2
which(duplicated(human_mouse_homo_gene_DEGs$H_Gene))


#annalysis about acTFs
#read Seurat object and load initial cAbortionata and matrix
#target_final<-readRDS(file = "/home/chenwei/10x_data/191125-merge_nonormalization_add/Cell_0512.rds")
target_final<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
target_group="final_major_subgroup_brief"
unique(target_final[[target_group]][,target_group])
ALL_cAbortionata<-target_final@meta.data

#prepare for cell.info_target_final
cell_name="ALLcell_all_genes"
target_final_group="Treat"
emb.tsne <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".tsne.txt"), sep = "\t", row.names = 1, header = T)
emb.umap <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".umap.txt"), sep = "\t", row.names = 1, header = T)
colnames(emb.tsne) <- paste0("tSNE_", 1:2)
colnames(emb.umap) <- paste0("UMAP_", 1:2)
#for all cells
cell.info_target_final<-ALL_cAbortionata
cell.info_target_final <- cbind(cell.info_target_final, emb.tsne[rownames(cell.info_target_final), ])
cell.info_target_final <- cbind(cell.info_target_final, emb.umap[rownames(cell.info_target_final), ])
head(cell.info_target_final)
write.table(cell.info_target_final,file= paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s4_",cell_name,".info.xls"),sep="\t", quote=F, row.names=T,col.names=T)
saveRDS(cell.info_target_final,paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s4_",cell_name,".info.rds"))

##load pyscenic output
cell_name="ALLcell_all_genes"

## 读入RAS矩阵--具体见supplment script
##rasMat(rds可直接读取): rows=cellID, cols=regulon, values=Regulon Acticity Score
#rasMat0 <- fread(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".AUCell.txt"), sep = "\t",  header = T, data.table = F)
#rownames(rasMat0) <- rasMat0$V1
#colnames(rasMat0) <- sub("(+)", "", colnames(rasMat0), fixed = T)
#rasMat0 <- rasMat0[, -1]
#saveRDS(rasMat0, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".rasMat.rds"))

rasMat<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".rasMat.rds"))
rasMat_binary <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".binary_mtx.txt"),header=T)
colnames(rasMat_binary)<-gsub("[...]*$","",colnames(rasMat_binary))
colnames(rasMat_binary)<-gsub("[...]","-",colnames(rasMat_binary))
all(colnames(rasMat_binary) == colnames(rasMat_binary))

rasMat_active_thresholds <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".auc_thresholds.txt"), header=T)
rownames(rasMat_active_thresholds)<-gsub("[(+)]","",rownames(rasMat_active_thresholds))
rasMat_active_thresholds$V1<-rownames(rasMat_active_thresholds)
colnames(rasMat_active_thresholds)<-c("thresholds","lable")

cell.info_target_final<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s4_",cell_name,".info.rds"))
head2d <- function(x) x[1:5, 1:5]
head2d(rasMat);head2d(rasMat_binary);head(cell.info_target_final)

#tissue.use <- cell.info %>%  filter(CellType == "Dendritic cell") %>% count(Tissue) %>% filter(n > 20) %>% `[[`("Tissue")
#tissue.use
#cells.bm <- cell.info %>% filter(CellType == "Dendritic cell" & Tissue == "BoneMarrow") %>%  rownames()
#cells.pbmc <- cell.info %>%  filter(CellType == "Dendritic cell" & Tissue == "PeripheralBlood") %>% rownames()
#cells.use <- c(cells.bm, cells.pbmc)
#length(cells.use)

#(rsd already)::Identification of  regulon with differential activity between Abortion and CTRL
compare_group<-"final_major_subgroup_brief"
cell_types_name <- names(table(cell.info_target_final[compare_group]))
cell_types_name1<- cell_types_name[which(!(cell_types_name %in% c("Troph_mix","Ery")))]
#since none of them have active_percent more than 50 or just no care
D_ac_TFs_res_list<-list(); D_ac_TFs_lable_list<-list(); D_ac_TFs_active_high_list<-list(); Cell_group_name<-c()
for ( name_Cell in cell_types_name1){
  # name_Cell<-"MyCs"
  print(name_Cell)
  cell_identity<-name_Cell
  cells.CTRL <- cell.info_target_final %>% filter(final_major_subgroup_brief  == cell_identity) %>% filter(Treat  == "CTRL") %>%  rownames()
  cells.Abortion <- cell.info_target_final %>% filter(final_major_subgroup_brief  == cell_identity) %>% filter(Treat  == "Abortion") %>%  rownames()
  cells.use <- c(cells.CTRL, cells.Abortion)
  length(cells.use)#4413
  
  ##wilcox test on RAS
  dim(rasMat)#49988   545
  rasMat.CTRL<- rasMat[cells.CTRL, ];rasMat.Abortion <- rasMat[cells.Abortion, ]
  mean.CTRL <- colMeans(rasMat.CTRL);which(mean.CTRL==0)
  mean.Abortion <- colMeans(rasMat.Abortion);which(mean.Abortion==0)
  ##Abortion specific
  Abortion_spec_ac_TFs<-setdiff(names(which(mean.CTRL==0)), names(which(mean.Abortion==0)))
  CTRL_spec_ac_TFs<-setdiff(names(which(mean.Abortion==0)), names(which(mean.CTRL==0)))
  Abortion_spec_ac_TFs;CTRL_spec_ac_TFs
  
  #statistion 
  p.val <- sapply(seq_along(colnames(rasMat)), function(x) wilcox.test(rasMat.Abortion[,x],rasMat.CTRL[,x])$p.value)
  p.adj <- p.adjust(p.val, method = "BH")
  fc <-mean.Abortion/mean.CTRL;range(fc)##0.5206128 1.9991561
  res <- data.frame(mean.Abortion = mean.Abortion,mean.CTRL = mean.CTRL,fAbortionchange = fc,p.val = p.val,p.adj = p.adj)
  
  head(res);res[Abortion_spec_ac_TFs,]
  res<-res[which(!(is.na(res$fAbortionchange))),]
  
  #filter significant regulon
  cut_off_pvalue = 0.01 ;cut_off_FC = 1 
  res$change = ifelse(res$p.val < cut_off_pvalue,as.character(ifelse(res$fAbortionchange > cut_off_FC,'Up','Down')),'Stable')
  res$change <-factor(res$change,levels=c('Up','Down','Stable'))
  
  res.label <- res %>% filter(p.val < cut_off_pvalue);dim(res.label)
  res.label$label <- rownames(res.label)
  res.label<-res.label[order(res.label$p.val,decreasing = F),]
  head(res.label);dim(res.label)
  res.label_all<-res.label
  head(res.label_all)
  write.table(as.data.frame(res.label_all), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_identity,"_abortion_def_TFs_p001.txt"),quote = F, row.names = T,sep = "\t")
  #res.label_all2<-read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_identity,"_Age_def_TFs_p001.txt"),header=T,sep = "\t")
  write.table(as.data.frame(res.label_all[which(res.label_all$p.adj<0.05),]), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_identity,"_abortion_def_TFs_p001_q005.txt"),quote = F, row.names = T,sep = "\t")
  
  #Reduce(intersect,list(res.label_all$label,res.label$label))
  #setdiff(res.label_all$label,res.label$label)
  #setdiff(res.label$label,res.label_all$label)
  
  #add active pecentage
  rasMat_binary1<-rasMat_binary[,res.label$label]
  Abortion_active_percent<-colSums(rasMat_binary1[cells.Abortion,])/nrow(rasMat_binary1[cells.Abortion,])
  CTRL_active_percent<-colSums(rasMat_binary1[cells.CTRL,])/nrow(rasMat_binary1[cells.CTRL,])
  Active_decition<-data.frame(t(rbind(Abortion_active_percent,CTRL_active_percent)))
  Active_decition$ratio<-"low"
  Active_decition[which(Active_decition$Abortion_active_percent>=0.5 | Active_decition$CTRL_active_percent>=0.5),]$ratio<-"high"
  
  res.label2<-merge(res.label,Active_decition,by=0);dim(res.label2)
  rownames(res.label2)<-res.label2$label
  res.label2<-res.label2[,c(-1)]
  ##add active_thresholds
  res.label3<-merge(res.label2,rasMat_active_thresholds,by=0);dim(res.label3)
  rownames(res.label3)<-res.label3$label
  res.label3<-res.label3[,c(-1)]
  head(res.label3)
  TF_high<-rownames(Active_decition[which(Active_decition$ratio == "high"),])
  res.label2[TF_high,]
  
  TF_active_high<-res.label3[which(res.label3$mean.Abortion>=res.label3$thresholds|res.label3$mean.CTRL >=res.label3$thresholds),]
  #  D_TF_most_restrict<-rownames(TF_active_high[which(TF_active_high$ratio == "high"),])
  
  write.table(as.data.frame(TF_active_high), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_identity,"_Age_def_TFs_p001_active_mean_high.txt"),quote = F, row.names = T,sep = "\t")
  write.table(as.data.frame(TF_active_high[which(TF_active_high$ratio == "high"),]), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_identity,"_Age_def_TFs_most_strict.txt"),quote = F, row.names = T,sep = "\t")
  D_ac_TFs_res_list<-c(D_ac_TFs_res_list,list(res))
  D_ac_TFs_lable_list<-c(D_ac_TFs_lable_list,list(res.label_all))
  D_ac_TFs_active_high_list <-c(D_ac_TFs_active_high_list,list(TF_active_high))
  Cell_group_name<-c(Cell_group_name,cell_identity)
  
}

length(D_ac_TFs_res_list);length(D_ac_TFs_lable_list);length(D_ac_TFs_active_high_list);length(Cell_group_name)#76;76
names(D_ac_TFs_res_list)<-Cell_group_name
names(D_ac_TFs_lable_list)<-Cell_group_name
names(D_ac_TFs_active_high_list)<-Cell_group_name

saveRDS(D_ac_TFs_res_list, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_name,"_D_ac_TFs_res_list.rds"))
saveRDS(D_ac_TFs_lable_list, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_name,"_D_ac_TFs_lable_list.rds"))
saveRDS(D_ac_TFs_active_high_list, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_name,"_D_ac_TFs_active_high_list.rds"))

D_ac_TFs_res_list<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_name,"_D_ac_TFs_res_list.rds"))
D_ac_TFs_lable_list<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_name,"_D_ac_TFs_lable_list.rds"))
D_ac_TFs_active_high_list<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/",cell_name,"_D_ac_TFs_active_high_list.rds"))

D_ac_TF_matrix<-data.frame()
for ( f_name in names(D_ac_TFs_lable_list)){
  # f_name<-"D_Mem_CD4_Ts"
  target_final_TF_matrix<-D_ac_TFs_lable_list[[f_name]]
  target_final_TF_matrix$cell_identity<-rep(f_name,nrow(target_final_TF_matrix))
  D_ac_TF_matrix<-rbind(D_ac_TF_matrix,target_final_TF_matrix)
}
rownames(D_ac_TF_matrix)<-NULL
head(D_ac_TF_matrix);tail(D_ac_TF_matrix);dim(D_ac_TF_matrix)#5612    8
write.table(D_ac_TF_matrix, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/merge_",cell_name,"_D_ac_TF_matrix.txt"), sep = "\t")

D_ac_TFs_active_high_matrix<-data.frame()
for ( f_name in names(D_ac_TFs_active_high_list)){
  # f_name<-"D_Mem_CD4_Ts"
  target_final_TF_matrix<-D_ac_TFs_active_high_list[[f_name]]
  target_final_TF_matrix$cell_identity<-rep(f_name,nrow(target_final_TF_matrix))
  D_ac_TFs_active_high_matrix<-rbind(D_ac_TFs_active_high_matrix,target_final_TF_matrix)
}
rownames(D_ac_TFs_active_high_matrix)<-NULL
head(D_ac_TFs_active_high_matrix);tail(D_ac_TFs_active_high_matrix);dim(D_ac_TFs_active_high_matrix)#921  13
write.table(D_ac_TFs_active_high_matrix, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/merge_",cell_name,"_D_ac_TFs_active_high_matrix.txt"), sep = "\t")

#next step 
##files prepared for cytoscape to show the relationship between TFs-active change and Cell types
#file1 daTFs-Cells-change-importance
#file2 Cells-origin1(Fetal & Maternal & Mix)-origin2(Tissue:Base & Base-Peri & PBMCs & Others)-daTFs_number1(total)-daTFs_number(Up)-daTFs_number(Down)
#file3 daTFs-TFs_type(Up_only&Down_only&Contrary)-Cell_numbers1(Up)-Cell_number2(Down)--Cell_numbers3(Up-Down)-Cell_number4(Up+Down)-percent(UP/(Up+Down))
#D_ac_TF_matrix <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/merge_",cell_name,"_D_ac_TF_matrix.txt"), header=T)
D_ac_TFs_active_high_matrix <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/merge_",cell_name,"_D_ac_TFs_active_high_matrix.txt"), header=T)

TF_Cell_net<-D_ac_TFs_active_high_matrix[,c("lable","cell_identity","change","fAbortionchange")]
colnames(TF_Cell_net)<-c("daTFs_identity","cell_identity","change","fAbortionchange")

## De-TFs number for each Cells
num_cell_plot<-D_ac_TFs_active_high_matrix[,c("cell_identity","change")]
num_cell_plot$Num<-1
Cell_node <-dcast(num_cell_plot, cell_identity ~ change , value.var = "Num")
colnames(Cell_node)<-c("cell_identity","daTFs_num_Down","daTFs_num_Up")
Cell_node$daTFs_num_sum<-Cell_node$daTFs_num_Down+Cell_node$daTFs_num_Up
#Cell_node2<-data.frame(table(cell.info_target_final$sub_cluster_brief,cell.info_target_final$cell.type))
#Cell_node2<-dcast(Cell_node2, Var1 ~  Var2 , value.var = "Freq")
#Cell_node2$mix_ratio<-Cell_node2$fetal/(Cell_node2$fetal+Cell_node2$maternal)*100
#table(Cell_node2$mix_ratio)
#Cell_node2$origin1<- ifelse(Cell_node2$mix_ratio>1  & Cell_node2$mix_ratio<100, "Mix",ifelse(Cell_node2$mix_ratio == 0, "maternal","Fetal"))
#Cell_node2<-Cell_node2[,c(1,5,4)]
#colnames(Cell_node2)<-c("cell_identity","origin1","mix_ratio")
#Cell_node<-merge(Cell_node,Cell_node2)

#Cell_node3<-data.frame(table(cell.info_target_final$sub_cluster_brief,cell.info_target_final$Tissue))
#Cell_node3<-dcast(Cell_node3, Var1 ~  Var2 , value.var = "Freq")
#Cell_node3<-merge(Cell_node2,Cell_node3)
com_trend_freq_active_high<-data.frame(table(D_ac_TFs_active_high_matrix$label,D_ac_TFs_active_high_matrix$change))

daTFs_node <-dcast(com_trend_freq_active_high, Var1 ~ Var2 , value.var = "Freq")
colnames(daTFs_node)<-c("daTFs_identity","Cell_num_Down","Cell_num_Up")
daTFs_node$Cell_num_sum<-daTFs_node$Cell_num_Down+daTFs_node$Cell_num_Up
daTFs_node$Cell_Down_percent<-daTFs_node$Cell_num_Down/daTFs_node$Cell_num_sum*100
daTFs_node$TFs_type<- ifelse(daTFs_node$Cell_num_Down> 0  & daTFs_node$Cell_num_Up > 0, "Contrary",ifelse(daTFs_node$Cell_num_Down == 0, "Up_only","Down_only"))
table(daTFs_node$TFs_type)
#Contrary Down_only   Up_only 
#  87        79       148 
head(TF_Cell_net);head(Cell_node);head(daTFs_node)
write.table(as.data.frame(TF_Cell_net), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/cytoscape_file/",cell_name,"_TF_Cell_Network.txt"),quote = F, row.names = F,sep = "\t")
write.table(as.data.frame(Cell_node), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/cytoscape_file/",cell_name,"_Cell_node_stat.txt"),quote = F, row.names = F,sep = "\t")
write.table(as.data.frame(daTFs_node), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/cytoscape_file/",cell_name,"_daTFs_node_stat.txt"),quote = F, row.names = F,sep = "\t")

TF_Cell_net2<-merge(TF_Cell_net,Cell_node)
TF_Cell_net2<-merge(TF_Cell_net2,daTFs_node)
#TF_Cell_net2$origin_change<-paste(TF_Cell_net2$origin1,TF_Cell_net2$change,sep="-" )
head(TF_Cell_net2);dim(TF_Cell_net2)# 921  12

#for 79  Down_TFs
Down_TFs_matrix<-data.frame(table(TF_Cell_net2[which(TF_Cell_net2$TFs_type == "Down_only"),]$daTFs_identity))
#Down_TFs_matrix2<-data.frame(dcast(Down_TFs_matrix, Var1 ~ Var2 , value.var = "Freq"))
#Down_TFs_matrix2[which(Down_TFs_matrix2$Mix.Down > 0),]
#Down_TFs_matrix2$TFs_type2<- ifelse(Down_TFs_matrix2$Mix.Down > 0, "Mix_related_Down",ifelse(Down_TFs_matrix2$Fetal.Down > 0 & Down_TFs_matrix2$maternal.Down > 0, "Both_MF_Down",ifelse(Down_TFs_matrix2$Fetal.Down > 0, "Fetal_Down","Maternal_Down")))
#table(Down_TFs_matrix2$TFs_type2)
##Both_MF_Down    Fetal_Down Maternal_Down   Mix_related_Down
##  53               31               49               30 
dim(Down_TFs_matrix)
#for 148   Up_TFs
Up_TFs_matrix<-data.frame(table(TF_Cell_net2[which(TF_Cell_net2$TFs_type == "Up_only"),]$daTFs_identity))
dim(Up_TFs_matrix)

#Up_TFs_matrix2<-data.frame(dcast(Up_TFs_matrix, Var1 ~ Var2 , value.var = "Freq"))
#Up_TFs_matrix2[which(Up_TFs_matrix2$Mix.Up > 0),]
#Up_TFs_matrix2$TFs_type2<- ifelse(Up_TFs_matrix2$Mix.Up > 0, "Mix_related_Up",ifelse(Up_TFs_matrix2$Fetal.Up > 0 & Up_TFs_matrix2$maternal.Up > 0, "Both_MF_Up",ifelse(Up_TFs_matrix2$Fetal.Up > 0, "Fetal_Up","Maternal_Up")))
#table(Up_TFs_matrix2$TFs_type2)
##Both_MF_Up    Fetal_Up Maternal_Up Mix_related_Up
#    4             14        5              8 

#for 87  Contrary TFs
Contrary_TFs_matrix<-data.frame(table(TF_Cell_net2[which(TF_Cell_net2$TFs_type == "Contrary"),]$daTFs_identity))
dim(Contrary_TFs_matrix)
#Contrary_TFs_matrix2<-data.frame(dcast(Contrary_TFs_matrix, Var1 ~ Var2 , value.var = "Freq"))
#Contrary_TFs_matrix2$TFs_type2<-"Change"
#Contrary_TFs_matrix2$TFs_type2<-ifelse(Contrary_TFs_matrix2$Mix.Down > 0 & Contrary_TFs_matrix2$Mix.Up > 0, "Mix_inne", ifelse(Contrary_TFs_matrix2$Mix.Down > 0 | Contrary_TFs_matrix2$Mix.Up > 0, "Mix_out",ifelse(Contrary_TFs_matrix2$Fetal.Up > 0 & Contrary_TFs_matrix2$Fetal.Down > 0 | Contrary_TFs_matrix2$maternal.Up > 0 & Contrary_TFs_matrix2$maternal.Down > 0,"inne_origin",ifelse(Contrary_TFs_matrix2$Fetal.Up > 0, "FM_Up_down","MF_Up_down"))))
#Contrary_TFs_matrix2$TFs_type2<-ifelse(Contrary_TFs_matrix2$Mix.Down > 0 & Contrary_TFs_matrix2$Mix.Up > 0, "Mix_inne", ifelse(Contrary_TFs_matrix2$Mix.Down > 0 | Contrary_TFs_matrix2$Mix.Up > 0, "Mix_out",ifelse(Contrary_TFs_matrix2$Fetal.Up > 0 & Contrary_TFs_matrix2$Fetal.Down > 0 & Contrary_TFs_matrix2$maternal.Up == 0 & Contrary_TFs_matrix2$maternal.Down == 0, "Fetal_inne_origin",ifelse(Contrary_TFs_matrix2$Fetal.Up == 0 & Contrary_TFs_matrix2$Fetal.Down == 0 & Contrary_TFs_matrix2$maternal.Up > 0 & Contrary_TFs_matrix2$maternal.Down > 0,"Maternal_inne_origin",ifelse(Contrary_TFs_matrix2$Fetal.Up > 0 & Contrary_TFs_matrix2$Fetal.Down > 0 | Contrary_TFs_matrix2$maternal.Up > 0 & Contrary_TFs_matrix2$maternal.Down > 0,"other_inne",ifelse(Contrary_TFs_matrix2$Fetal.Up > 0, "FM_Up_down","MF_Up_down"))))))
#table(Contrary_TFs_matrix2$TFs_type2)
#FM_Up_down inne_origin  MF_Up_down    Mix_inne     Mix_out 
#     4          95           1          13          56 

#daTFs_node3<-rbind(Down_TFs_matrix2[,c(1,ncol(Down_TFs_matrix2))],Up_TFs_matrix2[,c(1,ncol(Up_TFs_matrix2))],Contrary_TFs_matrix2[,c(1,ncol(Contrary_TFs_matrix2))])
#colnames(daTFs_node3)<-c("daTFs_identity","TFs_type2")
#daTFs_node3<-merge(daTFs_node,daTFs_node3)
#write.table(as.data.frame(daTFs_node3), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/cytoscape_file/",cell_name,"_daTFs_node_stat2.txt"),quote = F, row.names = F,sep = "\t")
##then detail result for Cell-TFs_network showed in cytoscape
#daTFs_node3<-read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/cytoscape_file/",cell_name,"_daTFs_node_stat2.txt"),header = T)

#save information file
#dim(Down_TFs_matrix2);dim(Up_TFs_matrix2);dim(Contrary_TFs_matrix2)#163   5  #31  5  #169   8
#daTFs_node4<-merge(TF_Cell_net2,daTFs_node3)
#write.csv(as.data.frame(daTFs_node4), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/cytoscape_file/",cell_name,"_daTFs_Cell_relationship.csv"),quote = F, row.names = F,sep = ",")

#next step show relationship for def_acTFs and Cell subtypes
##evalute the intercell type relationship for changes in TFs
dim(TF_Cell_net)# 921   4
TF_Cell_evaluate <-TF_Cell_net[,1:3]
TF_Cell_evaluate$Num<-ifelse(TF_Cell_evaluate$change == "Up",1,-1)
TF_Cell_evaluate2 <-dcast(TF_Cell_evaluate, cell_identity ~ daTFs_identity , value.var = "Num")
TF_Cell_evaluate2[is.na(TF_Cell_evaluate2)]<- 0
rownames(TF_Cell_evaluate2)<-TF_Cell_evaluate2$cell_identity
TF_Cell_evaluate2<-TF_Cell_evaluate2[,-1]
head(TF_Cell_evaluate2)
##http://www.bio-info-trainee.com/903.html
TF_Cell_dend = as.dendrogram(hclust(dist(TF_Cell_evaluate2,method="euclidean"), method = "ward.D"))
clusters <- cutree(TF_Cell_dend, h =9) # dendextend::cutree()
TF_Cell_dend = color_branches(TF_Cell_dend, h = 9, col = pal_d3("category20")(20))
plot(TF_Cell_dend)
TF_Cell_dend2 = as.dendrogram(hclust(dist(TF_Cell_evaluate2,method="euclidean"), method = "complete"))
clusters <- cutree(TF_Cell_dend2, h =8) # dendextend::cutree()
TF_Cell_dend2 = color_branches(TF_Cell_dend2, h = 8, col = pal_d3("category20")(20))
plot(TF_Cell_dend2)

### plot for cell type sepecify of de-acTFs
##group de-acTFs to evaluate the cor relationship 
TF_type<-"de_ac_TFs"
Up_acTF<-as.character(com_trend_freq_active_high[which(com_trend_freq_active_high$Var2 == "Up" & com_trend_freq_active_high$Freq>0),"Var1"])
Down_acTF<-as.character(com_trend_freq_active_high[which(com_trend_freq_active_high$Var2 == "Down" & com_trend_freq_active_high$Freq>0 ),"Var1"])
length(unique(Up_acTF));length(unique(Down_acTF))#235 166
length(unique(com_trend_freq_active_high$Var1))#314

regulon.names<-unique(c(Up_acTF,Down_acTF))

rasMat_de_high_acTF <- rasMat[, regulon.names]
dim(rasMat_de_high_acTF)#49988   314
pccMat <- cor(rasMat_de_high_acTF)

CSI <- function(r1, r2) {
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat), c(r1,r2))
  N <- sum(pccMat[r1, r.others] < delta) + sum(pccMat[r2, r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}
csiMat <- pblapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)))
csiMat <- do.call(rbind, csiMat)
rownames(csiMat) <- rownames(pccMat)

round(pccMat[1:10,1:10], 2)
round(csiMat[1:10,1:10], 2)
csiMat.binary <- matrix(as.numeric(csiMat >= 0.7), nrow = nrow(csiMat))
colnames(csiMat.binary) <- colnames(csiMat)
rownames(csiMat.binary) <- rownames(csiMat)
csiMat.binary[1:10,1:10]
#接下来即可对CSI matrix进行层次聚类；而CSI>0.7对CSI matrix进行二值化,所得csiMat.binary主要用于cytoscape中构建相关性网络。
saveRDS(csiMat, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s6_",cell_name,".",TF_type,"_csiMat.rds"))

mat = readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s6_",cell_name,".",TF_type,"_csiMat.rds"))

### Dendrogram
h = 5.5##这个值比较重要,用于确定分组，此处使用h法，还可以用k法固定分类数目。需要根据热图以及层次聚类图的结果来迭代确定此值
tree_height<-"h5.5"

row_dend = as.dendrogram(hclust(dist(mat), method = "complete"))
clusters <- cutree(row_dend, h = h) # dendextend::cutree()
row_dend = color_branches(row_dend, h = h, col = pal_d3("category20")(20))
plot(row_dend)

### Heatmap
col_range = c(0.7, 1)
col_fun <- colorRamp2(col_range, c("#FCF8DE", "#253177"))

ht <- Heatmap(matrix = mat,col = col_fun, name = "ht1",
              cluster_rows = TRUE,cluster_columns = TRUE,
              show_column_names = FALSE,show_row_names = FALSE, show_heatmap_legend = FALSE)
lgd <- Legend(col_fun = col_fun,title = "",  at = col_range, labels = c("low", "high"), 
              direction = "horizontal", legend_width = unit(1, "in"), border = FALSE)
draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))

decorate_heatmap_body("ht1", {
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(l) which(l)[1]
  last_index = function(l) { x = which(l); x[length(x)] }
  clusters <- names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  x1 = x1/length(ind)
  x2 = x2/length(ind)
  grid.rect(x = x1, width = (x2 - x1), y = 1-x1, height = (x1 - x2), 
            hjust = 0, vjust = 0, default.units = "npc", 
            gp = gpar(fill=NA, col="#FCB800", lwd=3))
  grid.text(label = paste0("M",clusters),x = x2-length(clusters)/length(ind), y = 1-x1-(x2-x1)/2,
            default.units = "npc",hjust = 1, vjust = 0.5, gp = gpar(fontsize=12, fontface="bold"))
  ##此处为添加字符方便查看，最后定下图则可以注释掉
})
decorate_column_dend("ht1", {
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(l) which(l)[1]
  last_index = function(l) { x = which(l); x[length(x)] }
  clusters <- names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
            default.units = "npc", gp = gpar(fill = pal_d3("category20")(20), alpha=.5, col = NA))
})

#用于与热图上的区块进行对应，方便查找motif
tree = column_dend(ht)
ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
clusters <- names(table(ind))
regulon.clusters <- data.frame(regulon=names(ind), cluster=paste0("M",ind))

write.table(regulon.clusters, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s6_",cell_name,".",TF_type,"_regulon_clusters_",tree_height,".txt"), sep = "\t", quote = F, row.names = F)
table(regulon.clusters$cluster)
##至此，我们确认了相关性较高的regulion群，接下来还需要确定这些regulon群是否与某一群或者某一状态的cell亚群有直接相关

regulon.clusters<-read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s6_",cell_name,".",TF_type,"_regulon_clusters_",tree_height,".txt"),header = T)
#regulon.clusters<-regulon.clusters[which(regulon.clusters$cluster %in% names(which(table(regulon.clusters$cluster)>1))),]
cell.info_target_final<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s4_",cell_name,".info.rds"))

moduleRasMat <- lapply(names(which(table(regulon.clusters$cluster)>1)), function(x) {
  regulon.use <- subset(regulon.clusters, cluster == x)$regulon
  rowMeans(rasMat[, regulon.use])
})

names(moduleRasMat) <- names(which(table(regulon.clusters$cluster)>1))
moduleRasMat2 <- do.call(cbind, moduleRasMat)
cell.info_target_final <- cbind(cell.info_target_final, moduleRasMat2[rownames(cell.info_target_final), ])

p.list <- lapply(names(moduleRasMat), function(module){
  data.use <- cell.info_target_final
  expression.color <- c("darkblue", "lightblue", "green", "yellow", "red")
  max.val <- quantile(data.use[, module], 0.99)
  low.val <- quantile(data.use[, module], 0.1)
  data.use[, module] <- ifelse(data.use[, module] > max.val, max.val, data.use[, module])
  ggplot(data.use, aes(tSNE_1, tSNE_2, color=get(module))) + 
    geom_point(size=0.05) + 
    theme_bw(base_size = 15) + 
    ggtitle(module) + 
    scale_color_gradientn(name = NULL, colors = expression.color) + 
    theme(legend.position = "right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = .5, face = "bold", size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black")
    )
})

cowplot::plot_grid(plotlist = p.list, ncol = ceiling(length(names(moduleRasMat))/5))

## 具体查看各个regulon对应的最有可能的Cell Type--在每个细胞的平均regulon活性的基础上面计算每个regulon群在各个细胞类群中的均值
#cellType.info <- read.table("../data/CellType.Info.txt", sep = "\t", header = T)
## group module score by clusterID
clusters <- names(table(cell.info_target_final$final_major_subgroup_brief))

moduleScoreInCellType <- lapply(clusters, function(x) {
  cells.use <- rownames(subset(cell.info_target_final, final_major_subgroup_brief == x))
  colMeans(cell.info_target_final[cells.use, names(moduleRasMat)])
})
names(moduleScoreInCellType) <- clusters
moduleScoreInCellType <- do.call(rbind, moduleScoreInCellType)
moduleScoreInCellType <- as.data.frame(moduleScoreInCellType)
#moduleScoreInCellType$CellType <- mapvalues(x = rownames(moduleScoreInCellType),from = cellType.info$Cluster,to = cellType.info$Cell.Type)
moduleScoreInCellType$CellType <- rownames(moduleScoreInCellType)

plotCellTypeRank <- function(data, module, topn=5){
  data.use <- data
  data.use <- data.use[order(data.use[, module], decreasing = TRUE), ]
  data.use$Rank <- 1:nrow(data.use)
  data.use$pt.col <- ifelse(data.use$Rank <= topn, "#007D9B", "#BECEE3")
  data.label <- head(data.use, n = topn)
  data.label$delta <- c(Inf, abs(diff(data.label[, module])))
  
  ggplot(data.use, aes(Rank, get(module))) + 
    geom_point(size=3, color=data.use$pt.col) + 
    geom_text_repel(inherit.aes = FALSE, data = data.label, aes(Rank, get(module), label=CellType), size=4, max.iter = 2e4) + 
    ggtitle(module) + 
    #  title("Module for Cell types")+
    ylab("Regulon activity score") + xlab("Cell type") + 
    theme_bw(base_size = 12) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5, face = "bold")
    )
}

p.list <- lapply(names(moduleRasMat), function(x) plotCellTypeRank(moduleScoreInCellType, module = x, topn = 10))
cowplot::plot_grid(plotlist = p.list, ncol = ceiling(length(names(moduleRasMat))/5))
cowplot::plot_grid(plotlist = p.list[1:20], ncol =5)
cowplot::plot_grid(plotlist = p.list[21:22], ncol =2)

##select presentive de_TFs for corresponding module and Cell subgroups
head(moduleScoreInCellType)
cell_type_list<-list()
module_name<-c()
for ( modules in colnames(moduleScoreInCellType[,-ncol(moduleScoreInCellType)])){
  # modules<-"M1"
  print(modules)
  cell_type<-moduleScoreInCellType[order(moduleScoreInCellType[, modules], decreasing = TRUE),]$CellType
  cell_type_list<-c(cell_type_list,list(cell_type))
  module_name<-c(module_name,modules)
}
length(cell_type_list);length(module_name)
names(cell_type_list)<-module_name

#TFs for selective Cells in each module
regulon.clusters<-read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s6_",cell_name,".",TF_type,"_regulon_clusters_",tree_height,".txt"),header = T)
TF_Cell_evaluate0 <-TF_Cell_net[,1:3]
#rssMat <- readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s5_",cell_name,"_cell_type.rssMat.rds"))
rasMat<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".rasMat.rds"))
binMat <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".binary_mtx.txt"), sep = "\t", header = T, row.names = 1, check.names = FALSE)
colnames(binMat) <- sub("(+)", "", colnames(binMat), fixed = T)
cell.info_target_final2 <- cbind(cell.info_target_final, binMat[rownames(cell.info_target_final), ])

cell_target_final<-unique(as.character(target_final$final_major_subgroup_brief))
#cell_target_final<-cell_types_name1
## Regulon Rank Plot
source("/mnt/data/chenwei/jianghai/0.script/SCENIC/plotRegulonRank.R")
source("/mnt/data/chenwei/jianghai/0.script/SCENIC/DimPlot_module.R")
p_module_list<-list();Module_name<-c();TFs_list<-list()
##number of TFs in M7 only 4 
#for (modules in module_name[which(!(module_name %in% c("M12", "M4", "M6","M7")))]){

for (modules in module_name){
  
  # modules<-"M7"
  print(modules)
  plot_group<-cell_target_final 
  Cell_select<-cell_type_list[[modules]][1:3]
  TF_Cell_select<-unique(TF_Cell_evaluate0[TF_Cell_evaluate0$cell_identity %in% Cell_select,]$daTFs_identity)
  regulon_select<-regulon.clusters[which(regulon.clusters$cluster == modules),]$regulon
  TFs_select<-Reduce(intersect,list(TF_Cell_select,regulon_select))
  #Regulon Specificity Score calculation 
  compare_group<-"final_major_subgroup_brief"
  ctMat_select<-cell.info_target_final[compare_group]
  ctMat_select$num<-ifelse(ctMat_select$final_major_subgroup_brief %in% Cell_select,1,0)
  colnames(ctMat_select) <- c("cell_type",paste0(modules,"_top4_cell"))
  rasMat_select<-rasMat[,TFs_select]
  
  rssMat_select <- pblapply(colnames(rasMat_select), function(i) {  1 - JSD(rbind(rasMat_select[, i], ctMat_select[, 2]), unit = 'log2', est.prob = "empirical")})
  rssMat_select2 <- do.call(rbind, rssMat_select)
  rownames(rssMat_select2) <- colnames(rasMat_select)
  colnames(rssMat_select2) <-paste0(modules,"_top3_cell")
  TF_top5<-head(rownames(rssMat_select2)[order(rssMat_select2, decreasing = TRUE)],n=5)
  if (length(TF_top5)<5) {
    print(paste0(modules,":TFs number smaller than 5"))
    next
  }
  
  TFs_list<-c(TFs_list,list(TF_top5))
  fig2Plot_module <- function(Modules,group,cell.type, Regulon) {
    p.list <- list(
      PlotRegulonRank(rssMat_select2, paste0(Modules,"_top3_cell")),
      DimPlot_module(cell.info_target_final2,module_NAME=Modules,group=group, cell.type = cell.type),
      DimPlot_module(cell.info_target_final2,regulon = Regulon) )
    cowplot::plot_grid(plotlist = p.list, ncol = 3,rel_widths = c(3,5,5))
  }
  pp1<-fig2Plot_module(modules,compare_group,Cell_select,TF_top5[1])
  pp2<-(DimPlot_module(cell.info_target_final2,regulon = TF_top5[2])+DimPlot_module(cell.info_target_final2,regulon = TF_top5[3]))/(DimPlot_module(cell.info_target_final2,regulon = TF_top5[4])+DimPlot_module(cell.info_target_final2,regulon = TF_top5[5]))
  p_module_combin<-cowplot::plot_grid(pp1,pp2, ncol = 2,rel_widths = c(8,4))
  p_module_list<-c(p_module_list,list(p_module_combin))
  Module_name<-c(Module_name,modules)
}

length(p_module_list);length(TFs_list);length(Module_name)#14
names(p_module_list)<-Module_name
names(TFs_list)<-Module_name

saveRDS(p_module_list, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s5_",cell_name,"_cell_TFs_module_plot.rds"))
saveRDS(TFs_list, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s5_",cell_name,"_cell_TFs_specific_TOP5.rds"))

p_module_list <- readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s5_",cell_name,"_cell_TFs_module_plot.rds"))
TFs_list <- readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s5_",cell_name,"_cell_TFs_specific_TOP5.rds"))

#Module_name
length(TFs_list);length(p_module_list)
TFs_list[["M9"]];p_module_list[["M9"]]
p_module_list[["M1"]]
##get special TFs matrix
#manual web:http://jaspar.genereg.net/
library(seqLogo)
#step0 在JASPAR数据库中下载给定转录因子的motif矩阵,存放至特定文件夹中
#1）在JASPAR数据库中检索Arx基因 2）找到对应的Logo，点击Logo 3）找到Frequency matrix，下载PFM矩阵
#step1 绘制motif logo##or using ggseqlogo绘制seqlogo图 ref: https://www.plob.org/article/12630.html
pfm.files <- list.files("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/seqLogo", "*.pfm", full.names = T)
for(file in pfm.files){
  print(file)
  pdf.file <- paste0(sub("(/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/seqLogo/.+)pfm", "\\1", file), "pdf")
  m <- read.table(file, skip = 1)
  m <- m/rep(colSums(m), each=4)
  motif<-makePWM(m)
  pdf(pdf.file, width = 5, height = 2)
  seqLogo(motif, xaxis=F, yaxis=F)
  dev.off()
}

#https://bioconductor.org/packages/release/bioc/vignettes/TFBSTools/inst/doc/TFBSTools.html
#library(TFBSTools)
#library(JASPAR2020)
#opts <- list()
#opts[["species"]] <- 9606
#opts[["type"]] <- "SELEX"
#opts[["all_versions"]] <- TRUE
#PFMatrixList <- getMatrixSet(JASPAR2020, opts)
#data("TBX19")
#####==待写出

##for age group
cell.info_target_final$sub_cluster_brief_Age<-paste(cell.info_target_final$sub_cluster_brief,cell.info_target_final$Age_group,sep="_")
clusters <- names(table(cell.info_target_final$sub_cluster_brief_Age))
#clusters <- names(table(cell.info_target_final$sub_cluster_brief))
moduleScoreInCellType <- lapply(clusters, function(x) {
  cells.use <- rownames(subset(cell.info_target_final, sub_cluster_brief_Age == x))
  colMeans(cell.info_target_final[cells.use, names(moduleRasMat)])
})
names(moduleScoreInCellType) <- clusters
moduleScoreInCellType <- do.call(rbind, moduleScoreInCellType)
moduleScoreInCellType <- as.data.frame(moduleScoreInCellType)
#moduleScoreInCellType$CellType <- mapvalues(x = rownames(moduleScoreInCellType),from = cellType.info$Cluster,to = cellType.info$Cell.Type)
moduleScoreInCellType$CellType <- rownames(moduleScoreInCellType)

p.list <- lapply(names(moduleRasMat), function(x) plotCellTypeRank(moduleScoreInCellType, module = x, topn = 10))
cowplot::plot_grid(plotlist = p.list, ncol = ceiling(length(names(moduleRasMat))/5))

##intergaration of DEGs and De-acTFs
#next discusion about DEGs with d-acTFs

Reduce(intersect,list(sc_Up_gene,Up_acTF))#57
Reduce(intersect,list(sc_Down_gene,Down_acTF))#25
TFs_DEGs_overlap<-Reduce(intersect,list(c(sc_Up_gene,sc_Down_gene),c(Up_acTF,Down_acTF)))
#78

grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(sc_Up_gene=sc_Up_gene,sc_Down_gene=sc_Down_gene,Up_acTF=Up_acTF,Down_acTF=Down_acTF),
                                   alpha=c(0.7,0.7,0.7,0.7),lwd=1,lty=1,col="black" ,  fill=ppCor[c(5:6,1:2)], 
                                   cex = 1.5, cat.col=ppCor[c(5:6,1:2)],cat.fontface=4,  cat.cex = 1.5,    
                                   main = paste0(cell_name,":",":p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)
grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(bulk_sc_DEGs_up=bulk_sc_DEGs_up,bulk_sc_DEGs_down=bulk_sc_DEGs_down,Up_acTF=Up_acTF,Down_acTF=Down_acTF),
                                   alpha=c(0.7,0.7,0.7,0.7),lwd=1,lty=1,col="black" ,  fill=ppCor[c(5:6,1:2)], 
                                   cex = 1.5, cat.col=ppCor[c(5:6,1:2)],cat.fontface=4,  cat.cex = 1.5,    
                                   main = paste0(cell_name,":",":p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)

##prepare for TFs_target_final  and corresponding cytoscape
#first preparing TFs_genes pairs
regulon_file0 <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".regulons.txt"), header=F)
rownames(regulon_file0)<-unlist(lapply(strsplit(regulon_file0$V1,"[(]"), function(x) x[1]))
regulon_file0$V1<-rownames(regulon_file0)
regulon_file1<-regulon_file0[,c("V1","V3")]

TFs_name <- regulon_file1$V1
TF_gene<-data.frame()
for ( f_name in TFs_name){
  # f_name<-"ALX4"
  target_final_gene<-unlist(strsplit(regulon_file1[f_name,]$V3,","))
  TFs<-rep(f_name,length(target_final_gene))
  TF_gene_raw<-data.frame(TFs=TFs,target_final_gene=target_final_gene)
  TF_gene<-rbind(TF_gene,TF_gene_raw)
}
head(TF_gene);tail(TF_gene);dim(TF_gene)
write.table(TF_gene, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".Cyto_Network.txt"), sep = "\t")

TF_gene_ALL <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".Cyto_Network.txt"), header=T)
head(TF_gene_ALL)

##filter TFs for downstream
head(TF_Cell_net);head(daTFs_node3)
special_TFs_DEGs_merge<-DEGs_merge[which(DEGs_merge$gene %in% unique(c(Reduce(intersect,list(unique(DEGs_merge$gene),Up_acTF)),Reduce(intersect,list(genAge_database,unique(DEGs_merge$gene),Down_acTF))))),]
special_TFs_DEGs_overlap<-DEGs_merge[which(DEGs_merge$gene %in% TFs_DEGs_overlap),]
table(special_TFs_DEGs_merge$cluster)
table(special_TFs_DEGs_overlap$cluster)

special_TFs_Cell0<-TF_Cell_net[which(TF_Cell_net$daTFs_identity %in% TFs_DEGs_overlap),]
special_daTFs_node0<-daTFs_node3[which(daTFs_node3$daTFs_identity %in%TFs_DEGs_overlap),]

special_TFs_Cell<-TF_Cell_net[which(TF_Cell_net$daTFs_identity %in% unique(c(Reduce(intersect,list(genAge_database,unique(DEGs_merge$gene),Up_acTF)),Reduce(intersect,list(genAge_database,unique(DEGs_merge$gene),Down_acTF))))),]
special_daTFs_node<-daTFs_node3[which(daTFs_node3$daTFs_identity %in% unique(c(Reduce(intersect,list(genAge_database,unique(DEGs_merge$gene),Up_acTF)),Reduce(intersect,list(genAge_database,unique(DEGs_merge$gene),Down_acTF))))),]
table(special_TFs_Cell$daTFs_identity,special_TFs_Cell$cell_identity)
table(special_TFs_Cell$daTFs_identity); table(special_TFs_Cell$cell_identity)

##以下待做 20210919
##select target_final regulon and corresponding weight
TF_gene_weight <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s2_",cell_name,".adj.txt"), header=T)  
head(TF_gene_weight);dim(TF_gene_weight)
TF_gene_ALL <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".Cyto_Network.txt"), header=T)
head(TF_gene_ALL);dim(TF_gene_ALL)

TF_gene_weight1<-TF_gene_weight[which(TF_gene_weight$TF %in% TFs_DEGs_overlap),]
head(TF_gene_weight1);dim(TF_gene_weight1)
TF_gene_weight1$pair<-paste0(TF_gene_weight1$TF,"_",TF_gene_weight1$target)
nrow(TF_gene_weight1)#216822
head(TF_gene_weight1)

TF_gene_final<-TF_gene_ALL[which(TF_gene_ALL$TFs %in% TFs_DEGs_overlap),]
head(TF_gene_final);dim(TF_gene_final)
TF_gene_final$pair<-paste0(TF_gene_final$TFs,"_",TF_gene_final$target_final_gene)
head(TF_gene_final)

TF_gene_final1<-TF_gene_weight1[which(TF_gene_weight1$pair %in% TF_gene_final$pair),]
nrow(TF_gene_final1);nrow(TF_gene_final)#24994 25008
range(TF_gene_final1$importance) ##8.500642e-03 1.882429e+02
table(TF_gene_final1$TF);length(table(TF_gene_final1$target))

setdiff(TF_gene_final$pair,TF_gene_final1$pair)
#[1] "ATF3_ATF3"     "CEBPB_CEBPB"   "CREM_CREM"     "EGR1_EGR1"     "ETS2_ETS2"     "FOS_FOS"       "FOSB_FOSB"     "HIF1A_HIF1A"  
#[9] "IRF1_IRF1"     "JUN_JUN"       "JUNB_JUNB"     "KLF4_KLF4"     "KLF5_KLF5"     "REL_REL"       "TFAP2A_TFAP2A" "YBX1_YBX1" 

TFs_target_DEGs_overlap<-Reduce(intersect,list(unique(DEGs_merge$gene),unique(TF_gene_final1$target)))
length(TFs_target_DEGs_overlap)#572

TF_gene_final_stat<-data.frame(table(c(TF_gene_final1$TF,TF_gene_final1$target)))
colnames(TF_gene_final_stat)<-c("gene","count")
TF_gene_final_stat<-TF_gene_final_stat[order(TF_gene_final_stat$count,decreasing = T),]
TF_gene_final_stat$type<-ifelse(TF_gene_final_stat$gene %in% TFs_DEGs_overlap,"TFs","target_genes")
head(TF_gene_final_stat)
table(TF_gene_final_stat$type)
#target_genes             TFs 
#           8028                 19 
write.table(as.data.frame(TF_gene_final1[,1:3]), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/cytoscape_file/",cell_name,"_DEGs_overlap_TFs_Network.txt"),quote = F, row.names = F,sep = "\t")
write.table(as.data.frame(TF_gene_final_stat), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/cytoscape_file/",cell_name,"_DEGs_overlap_TFs_stat.txt"),quote = F, row.names = F,sep = "\t")
write.table(as.data.frame(special_daTFs_node0[,c(1,6,7)]), paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/cytoscape_file/",cell_name,"_DEGs_overlap_TFs_stat2.txt"),quote = F, row.names = F,sep = "\t")

grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(genAge_database=genAge_database,All_DEGs=unique(DEGs_merge$gene),target_final_acTF_DEGs_overlap=unique(as.character(TF_gene_final_stat[TF_gene_final_stat$type == "target_genes",]$gene))),
                                   alpha=c(0.7,0.7,0.7),lwd=1,lty=1,col="black" ,  fill=ppCor[c(5:6,3)], 
                                   cex = 1.5, cat.col=ppCor[c(5:6,3)],cat.fontface=4,  cat.cex = 1.5,    
                                   main = paste0(cell_name,":",":p<0.01(active_high TFs):target_final gene"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)
grid.newpage()