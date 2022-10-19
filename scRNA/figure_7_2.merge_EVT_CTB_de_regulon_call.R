rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
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
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)

ALL_coldata<-target_final@meta.data
head(ALL_coldata)
ALL_coldata_CTBs_EVTs<-ALL_coldata[which(ALL_coldata$re_annotation_TFac %in% c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3")),]
dim(ALL_coldata_CTBs_EVTs)#29882    65
ALL_coldata_CTBs_EVTs$cluster_type<- as.character(ALL_coldata_CTBs_EVTs$re_annotation_TFac)
ALL_coldata_CTBs_EVTs[which(ALL_coldata_CTBs_EVTs$re_annotation_TFac  %in% c("CTBs_1","CTBs_2")),]$cluster_type<-"CTBs"
ALL_coldata_CTBs_EVTs[which(ALL_coldata_CTBs_EVTs$re_annotation_TFac  %in% c("EVTs_1","EVTs_2","EVTs_3")),]$cluster_type<-"EVTs"

##load pyscenic output
cell_name="ALLcell_all_genes"

## 读入RAS矩阵--supplment script
rasMat<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".rasMat.rds"))
rasMat_binary <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".binary_mtx.txt"),header=T)
colnames(rasMat_binary)<-gsub("[...]*$","",colnames(rasMat_binary))
colnames(rasMat_binary)<-gsub("[...]","-",colnames(rasMat_binary))
all(colnames(rasMat_binary) == colnames(rasMat_binary))

rasMat_active_thresholds <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".auc_thresholds.txt"), header=T)
rownames(rasMat_active_thresholds)<-gsub("[(+)]","",rownames(rasMat_active_thresholds))
rasMat_active_thresholds$V1<-rownames(rasMat_active_thresholds)
colnames(rasMat_active_thresholds)<-c("thresholds","lable")

cell.info_target_final<-ALL_coldata_CTBs_EVTs
head2d <- function(x) x[1:5, 1:5]
head2d(rasMat);head2d(rasMat_binary);head(cell.info_target_final)

#(rsd already)::Identification of  regulon with differential activity between Abortion and CTRL
compare_group<-"cluster_type"
cell_types_name <- names(table(cell.info_target_final[compare_group]))
# "CTBs" "EVTs"
#since none of them have active_percent more than 50 or just no care
D_ac_TFs_res_list<-list(); D_ac_TFs_lable_list<-list(); D_ac_TFs_active_high_list<-list(); Cell_group_name<-c()
for ( name_Cell in cell_types_name1){
  # name_Cell<-"EVTs"
  print(name_Cell)
  cell_identity<-name_Cell
  cells.CTRL <- cell.info_target_final %>% filter(cluster_type  == cell_identity) %>% filter(Treat  == "CTRL") %>%  rownames()
  cells.Abortion <- cell.info_target_final %>% filter(cluster_type  == cell_identity) %>% filter(Treat  == "Abortion") %>%  rownames()
  cells.use <- c(cells.CTRL, cells.Abortion)
  length(cells.use)#7837
  
  ##wilcox test on RAS
  dim(rasMat)#49988   545
  rasMat.CTRL<- rasMat[cells.CTRL, ];rasMat.Abortion <- rasMat[cells.Abortion, ]
  mean.CTRL <- colMeans(rasMat.CTRL);which(mean.CTRL==0)
  mean.Abortion <- colMeans(rasMat.Abortion);which(mean.Abortion==0)
  ##Abortion specific
  Abortion_spec_ac_TFs<-setdiff(names(which(mean.CTRL==0)), names(which(mean.Abortion==0)))
  CTRL_spec_ac_TFs<-setdiff(names(which(mean.Abortion==0)), names(which(mean.CTRL==0)))
  Abortion_spec_ac_TFs;CTRL_spec_ac_TFs
  ##For CTBs
  #[1] "TBX15"  "ZNF214"
  #[1] "AR"     "FOXD4"  "LMX1B"  "SOX11"  "SOX12"  "TAL2"   "ZNF239"
  ####For EVTs
  ##"AR"     "EGR3"   "ESR2"   "FOXD4"  "MEIS3"  "RARB"   "SOX11"  "SOX14"  "TAL2"   "TCF23"  "TGIF2"  "VAX2"   "ZNF214"
  
  #statistion 
  p.val <- sapply(seq_along(colnames(rasMat)), function(x) wilcox.test(rasMat.Abortion[,x],rasMat.CTRL[,x])$p.value)
  p.adj <- p.adjust(p.val, method = "BH")
  fc <-mean.Abortion/mean.CTRL;range(fc)
  res <- data.frame(mean.Abortion = mean.Abortion,mean.CTRL = mean.CTRL,foldchange = fc,p.val = p.val,p.adj = p.adj)
  
  head(res);res[Abortion_spec_ac_TFs,]
  res<-res[which(!(is.na(res$foldchange))),]
  
  #filter significant regulon
  cut_off_pvalue = 0.05 ;cut_off_FC = 1 
  res$change = ifelse(res$p.val < cut_off_pvalue,as.character(ifelse(res$foldchange > cut_off_FC,'Up','Down')),'Stable')
  res$change <-factor(res$change,levels=c('Up','Down','Stable'))
  
  res.label <- res %>% filter(p.val < cut_off_pvalue);dim(res.label)
  res.label$label <- rownames(res.label)
  res.label<-res.label[order(res.label$p.val,decreasing = F),]
  head(res.label);dim(res.label)
  res.label_all<-res.label
  head(res.label_all)
  
  write.table(as.data.frame(res.label_all), paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/CTBs_EVTs_diff_regulon/",cell_identity,"_USM_def_TFs_p005.txt"),quote = F, row.names = T,sep = "\t")
  #res.label_all2<-read.table(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Different_TFs/",cell_identity,"_Age_def_TFs_p001.txt"),header=T,sep = "\t")
  write.table(as.data.frame(res.label_all[which(res.label_all$p.adj<0.05),]), paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/CTBs_EVTs_diff_regulon/",cell_identity,"_USM_def_TFs_p005_q005.txt"),quote = F, row.names = T,sep = "\t")
  
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
  
  write.table(as.data.frame(TF_active_high), paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/CTBs_EVTs_diff_regulon/",cell_identity,"_USM_def_TFs_p005_active_mean_high.txt"),quote = F, row.names = T,sep = "\t")
  #write.table(as.data.frame(TF_active_high[which(TF_active_high$ratio == "high"),]), paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Different_TFs/",cell_identity,"_USM_def_TFs_most_strict.txt"),quote = F, row.names = T,sep = "\t")
  D_ac_TFs_res_list<-c(D_ac_TFs_res_list,list(res))
  D_ac_TFs_lable_list<-c(D_ac_TFs_lable_list,list(res.label_all))
  D_ac_TFs_active_high_list <-c(D_ac_TFs_active_high_list,list(TF_active_high))
  Cell_group_name<-c(Cell_group_name,cell_identity)
  
}

length(D_ac_TFs_res_list);length(D_ac_TFs_lable_list);length(D_ac_TFs_active_high_list);length(Cell_group_name)#76;76
names(D_ac_TFs_res_list)<-Cell_group_name
names(D_ac_TFs_lable_list)<-Cell_group_name
names(D_ac_TFs_active_high_list)<-Cell_group_name

De_regulon_active_high_matrix<-data.frame()
for ( f_name in names(D_ac_TFs_active_high_list)){
  # f_name<-"D_Mem_CD4_Ts"
  target_final_TF_matrix<-D_ac_TFs_active_high_list[[f_name]]
  target_final_TF_matrix$cell_identity<-rep(f_name,nrow(target_final_TF_matrix))
  De_regulon_active_high_matrix<-rbind(De_regulon_active_high_matrix,target_final_TF_matrix)
}
rownames(De_regulon_active_high_matrix)<-NULL
head(De_regulon_active_high_matrix);tail(De_regulon_active_high_matrix);dim(De_regulon_active_high_matrix)#127   13
write.table(De_regulon_active_high_matrix, "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/CTBs_EVTs_diff_regulon/merge_CTB_EVT_De_Regulons_active_high_matrix.txt", sep = "\t")


###perform comparing the TFs of H3K27ac and regulon
##read TFs selected with activity and expression trend between USM and CTRL
EVT_all_TF <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/evt_all_TF.csv", header = T,row.names= 1)
EVT_up_TF  <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/evt_up_TF.csv", header = T,row.names= 1)
EVT_down_TF  <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/evt_down_TF.csv", header = T,row.names= 1)
length(EVT_all_TF$x);length(EVT_up_TF$x);length(EVT_down_TF$x)##101  27  30
CTB_all_TF <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/ctb_all_TF.csv", header = T,row.names= 1)
CTB_up_TF  <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/ctb_up_TF.csv", header = T,row.names= 1)
CTB_down_TF  <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/ctb_down_TF.csv", header = T,row.names= 1)
length(CTB_all_TF$x);length(CTB_up_TF$x);length(CTB_down_TF$x) #101 26 20

#read regulon differential in TFs
De_regulon_EVT_CTB <- read.table( "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/CTBs_EVTs_diff_regulon/merge_CTB_EVT_De_Regulons_active_high_matrix.txt", header=T)
dim(De_regulon_EVT_CTB)# 127  13
CTB_De_regulon<-De_regulon_EVT_CTB[which(De_regulon_EVT_CTB$cell_identity =="CTBs"),]
dim(CTB_De_regulon)#39 13
CTB_up_regulon<-CTB_De_regulon[which(CTB_De_regulon$change =="Up"),]
CTB_down_regulon<-CTB_De_regulon[which(CTB_De_regulon$change =="Down"),]

EVT_De_regulon<-De_regulon_EVT_CTB[which(De_regulon_EVT_CTB$cell_identity  =="EVTs"),]
dim(EVT_De_regulon)#88 13
EVT_up_regulon<-EVT_De_regulon[which(EVT_De_regulon$change =="Up"),]
EVT_down_regulon<-EVT_De_regulon[which(EVT_De_regulon$change =="Down"),]

###perform venn diagram
#EVT
venn <-venn.diagram(list(regulon_Up_EVT=unique(EVT_up_regulon$lable),regulon_Down_EVT=unique(EVT_down_regulon$lable),
                         Up_EVT=EVT_up_TF$x,Down_EVT=EVT_down_TF$x),
                    alpha=c(0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=ppCor[c(1,2,3,6)], 
                    cex = 1.5,cat.col=ppCor[c(1,2,3,6)], cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)

grid.newpage(); grid.draw(venn)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/merge_EVT_common_de_regulon_TF_H3K37ac_venn.pdf",width = 8,height = 6)
grid.draw(venn)
dev.off()

#CTB
venn <-venn.diagram(list(regulon_Up_CTB=unique(CTB_up_regulon$lable),regulon_Down_CTB=unique(CTB_down_regulon$lable),
                         Up_CTB=CTB_up_TF$x,Down_CTB=CTB_down_TF$x),
                    alpha=c(0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=ppCor[c(1,2,3,6)], 
                    cex = 1.5,cat.col=ppCor[c(1,2,3,6)], cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)

grid.newpage(); grid.draw(venn)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/merge_CTBs_common_de_regulon_TF_H3K37ac_venn.pdf",width = 8,height = 6)
grid.draw(venn)
dev.off()

#selected unique trend TFs：CTBs
CTB_up_target<-Reduce(intersect,list(unique(CTB_up_regulon$lable),CTB_up_TF$x))##"TEF"    "TFAP2A"
CTB_down_target<-Reduce(intersect,list(unique(CTB_down_regulon$lable),CTB_down_TF$x))#
CTB_target_data<-CTB_De_regulon[which(CTB_De_regulon$label %in% c(CTB_up_target,CTB_down_target)),]
CTB_target_data[,c("lable", "cell_identity","change","foldchange","p.val","p.adj")]
#    lable cell_identity change foldchange         p.val         p.adj
#29    TEF          CTBs     Up   1.427734 2.728799e-122 1.052813e-121
#30 TFAP2A          CTBs     Up   1.309563 7.011616e-203 4.540856e-202

#selected unique trend DEGs
##EVTs
EVT_up_target<-Reduce(intersect,list(unique(EVT_up_regulon$lable),EVT_up_TF$x))##"ARNT2" "ATF3"  "SMAD3" "TEAD1" 
EVT_down_target<-Reduce(intersect,list(unique(EVT_down_regulon$lable),EVT_down_TF$x))## "ELF3"   "GATA2"  "HSF1"   "SNAI1"  "TFAP2A"

EVT_target_data<-EVT_De_regulon[which(EVT_De_regulon$label %in% c(EVT_up_target,EVT_down_target)),]
EVT_target_data[,c("lable", "cell_identity","change","foldchange","p.val","p.adj")]
#  lable cell_identity change foldchange         p.val         p.adj
#  ARNT2          EVTs     Up 19.1559032  0.000000e+00  0.000000e+00
#   ATF3          EVTs     Up  1.0874811  0.000000e+00  0.000000e+00
#   ELF3          EVTs   Down  0.8792147  0.000000e+00  0.000000e+00
#  GATA2          EVTs   Down  0.8924066 2.247207e-277 4.909560e-277
#   HSF1          EVTs   Down  0.3757229  0.000000e+00  0.000000e+00
#  SMAD3          EVTs     Up  2.9715904  0.000000e+00  0.000000e+00
#  SNAI1          EVTs   Down  0.8166780  0.000000e+00  0.000000e+00
#  TEAD1          EVTs     Up  1.2349047 1.331109e-176 2.514317e-176
# TFAP2A          EVTs   Down  0.8800474 3.371902e-171 6.303487e-171

####save activity of target regulon for yaxi
rasMat[1:4,1:4]
gene_signature<-c("ARNT2","ATF3","SMAD3","TEAD1","ELF3","GATA2","HSF1","SNAI1","TFAP2A","TEF")
head(cell.info_target_final)
Target_rasMat <- rasMat[rownames(cell.info_target_final), gene_signature]
Target_rasMat[1:4,];dim(Target_rasMat)
ALL_coldata_CTBs_EVTs[1:4,];dim(ALL_coldata_CTBs_EVTs)
Target_rasMat2<-merge(Target_rasMat,ALL_coldata_CTBs_EVTs[,c("Treat","cluster_type")],by=0)
Target_rasMat2[1:4,]
rownames(Target_rasMat2)<-Target_rasMat2$Row.names
Target_rasMat2$group<-paste(Target_rasMat2$cluster_type,Target_rasMat2$Treat,sep="_")
Target_rasMat3<-Target_rasMat2[,-c(1,12,13)]
Target_rasMat3[1:4,];dim(Target_rasMat3)#29882    11
write.table(Target_rasMat3,file ="/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/activity_level_regulon_CTBs_EVTs_TFs_regulon_activity.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

###plot violin plot for activity of regulon
head(Target_rasMat2)
gene_list<-list()
for ( genename in gene_signature){
  # genename <-"ATF3"
  print(genename)
  rasMat_data_target<- Target_rasMat2[,c("cluster_type","Treat",genename)]
  colnames(rasMat_data_target)<-c("cell_type","group","Activity")
  plot_Cell_exp<- ggplot(rasMat_data_target,aes(x = cell_type,y=Activity,fill = group)) + 
    geom_violin(scale = "width") + ylim(0,max(rasMat_data_target$Activity)*1.05)+
    theme_bw()+scale_fill_manual(values=ppCor[c(10,6)])+
    labs(y = "Activity level", x = "Cell_type", title = genename)+theme_classic()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(plot.title = element_text(hjust = 0.5),axis.title = element_text(size=10),axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),legend.position = "right")+
    theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))+
    stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
    stat_compare_means(method = "wilcox.test",data=rasMat_data_target,aes(x=cell_type,y = Activity,group = group,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((rasMat_data_target$Activity)))
    gene_list<-c(gene_list,list(plot_Cell_exp))
}

gene_list_merge<-grid.arrange(gene_list[[1]], gene_list[[2]],gene_list[[3]], gene_list[[4]],gene_list[[5]],
                              gene_list[[6]], gene_list[[7]],gene_list[[8]], gene_list[[9]],gene_list[[10]],ncol=5)
Group_tag<-"merge_USM_TF_select_0713"
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/regulon_activity_",Group_tag,"_violin_target_regulon_plot_stat.pdf"),gene_list_merge,width=15, height=6)


###plot the network between  TFs and their target USM-related DEGs in seleted regulon

options(stringsAsFactors = FALSE)
cell_name="ALLcell_all_genes"

TF_gene_ALL <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".Cyto_Network.txt"), header=T)
head(TF_gene_ALL)
TF_gene_weight <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s2_",cell_name,".adj.tsv"), header=T)  
head(TF_gene_weight);dim(TF_gene_weight)

##read single cell DEGs
#get overlapped gene between developped genes and abortion DEGs 
merge_data0 <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge_new.txt", header=T)
head(merge_data0)
cell_subtype<-c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3")
USM_CTBs_EVTs_DEGs_data1<-merge_data0[which(merge_data0$cluster %in% cell_subtype),];dim(USM_CTBs_EVTs_DEGs_data1)#1976    6

##filter TF_regulon with only DEGs in CTBs and EVTs
USM_CTBs_EVTs_DEGs<-unique(as.character(USM_CTBs_EVTs_DEGs_data1$gene))
TF_gene_filter<-TF_gene_ALL[which(TF_gene_ALL$target_final_gene %in% USM_CTBs_EVTs_DEGs),];dim(TF_gene_filter)#16510     2
head(TF_gene_filter)
head(USM_CTBs_EVTs_DEGs_data1)

###plot network for all target TFs and DEGs
##for CTBs
TFs_name<-c("TEF","TFAP2A")
USM_CTBs_DEGs_data1<-USM_CTBs_EVTs_DEGs_data1[which(USM_CTBs_EVTs_DEGs_data1$cluster %in% c("CTBs_1","CTBs_2")),]
USM_CTBs_DEGs<-unique(as.character(USM_CTBs_DEGs_data1$gene))
TF_gene_filter_CTBs<-TF_gene_ALL[which(TF_gene_ALL$target_final_gene %in% USM_CTBs_DEGs),];dim(TF_gene_filter_CTBs)#5119    2

dim(USM_CTBs_DEGs_data1)# 598   6
TFs_name<-c("TEF","TFAP2A")
target_TF_filter_CTBs<-TF_gene_filter_CTBs[which(TF_gene_filter_CTBs$TFs %in% TFs_name),]
targets_trend<-USM_CTBs_DEGs_data1[which(USM_CTBs_DEGs_data1$gene %in% unique(c(unique(target_TF_filter_CTBs$target_final_gene),TFs_name))),]

####select target_final regulon and corresponding weight
TF_gene_weight1<-TF_gene_weight[which(TF_gene_weight$TF %in% TFs_name),]
head(TF_gene_weight1);dim(TF_gene_weight1)
TF_gene_weight1$pair<-paste0(TF_gene_weight1$TF,"_",TF_gene_weight1$target)
nrow(TF_gene_weight1)#216822
head(TF_gene_weight1)

target_TF_gene_final<-target_TF_filter_CTBs
target_TF_gene_final$pair<-paste0(target_TF_gene_final$TFs,"_",target_TF_gene_final$target_final_gene)
head(target_TF_gene_final)
target_TF_gene_final_1<-TF_gene_weight1[which(TF_gene_weight1$pair %in% target_TF_gene_final$pair),]
target_TF_gene_final_1
nrow(target_TF_gene_final_1);nrow(target_TF_gene_final_1)#7
range(target_TF_gene_final_1$importance) ##1.369561 61.109176
table(target_TF_gene_final_1$TF);length(table(target_TF_gene_final_1$target))
#TEF TFAP2A 
# 1      6 
##get gene pair link number
TF_gene_final_stat<-data.frame(table(c(target_TF_gene_final_1$TF,target_TF_gene_final_1$target)))
colnames(TF_gene_final_stat)<-c("gene","count")
TF_gene_final_stat<-TF_gene_final_stat[order(TF_gene_final_stat$count,decreasing = T),]
TF_gene_final_stat$type<-ifelse(TF_gene_final_stat$gene %in% unique(target_TF_gene_final_1$TF),"TFs","targets")
head(TF_gene_final_stat)
table(TF_gene_final_stat$type)
#targets     TFs 
#    7       2 
TF_gene_final_stat$node<-as.character(TF_gene_final_stat$gene)

##get Trend and cell number of DEGs
head(targets_trend)
targets_trend2<-targets_trend[which(!(targets_trend$gene %in% TFs_name)),]
targets_split <-dcast(targets_trend2, gene ~ trend)
#targets_split$Cell_sum<-targets_split$Down+targets_split$Up
#targets_split$Up_percent<-targets_split$Up/targets_split$Cell_sum*100
#targets_split$DEG_type<- ifelse(targets_split$Down> 0  & targets_split$Up > 0, "Contrary",ifelse(targets_split$Down == 0, "Up_only","Down_only"))
targets_split$Cell_sum<-targets_split$Up
targets_split$DEG_type<- "Up_only"
table(targets_split$DEG_type)
#Up_only 
# 7 

head(targets_split)

##merge for node information
node_info_merge<-merge(TF_gene_final_stat,targets_split,by="gene",all=T)

algorithm_plot<-"nicely";regulon_DEG<-"merge_CTBs_TEF_TFAP2A"
node_info_merge[which(is.na(node_info_merge$DEG_type)),]$DEG_type<-"TFs"
node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$Cell_sum<-node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$count
node_info_merge[which(node_info_merge$DEG_type=="TFs"),]
node_info_merge$color<-node_info_merge$DEG_type
node_info_merge$color<-ifelse(node_info_merge$color =="Up_only","indianred3",ifelse(node_info_merge$color =="Down_only","steelblue",
                                                                                    ifelse(node_info_merge$color =="Contrary","lightgreen","purple")))

nodes<-node_info_merge
edges<-target_TF_gene_final_1[,c("TF","target","importance")]

net<-graph_from_data_frame(d=edges,vertices = nodes,directed = T)
node_color<-node_info_merge$color

set.seed(19921010)
layout <- create_layout(net, layout ="nicely")
head(layout)
CTBs_all_network<-ggraph(layout)+
  geom_edge_fan(aes(edge_width=importance),color="lightblue",alpha=0.7) + #,show.legend=FALSE
  geom_node_point(aes(size=Cell_sum,color=as.factor(DEG_type)))+
  geom_node_text(aes(label=node),size=3)+
  #geom_node_text(aes(filter= DEG_type =="TFs",label=node),size=3)+
  #geom_node_text(aes(filter= (DEG_type =="TFs" | Cell_sum>1),label=node),size=3)+
  scale_color_manual(limits = as.factor(layout$DEG_type), values =node_color)+ 
  #  scale_edge_color_continuous(low = "cyan",high = "red")+
  scale_size_area(max_size = 20)+guides(size="none")+
  scale_edge_width(range = c(0.5, 3))+ 
  guides(color=guide_legend(order=3))+  theme_graph()+ 
  ggtitle(paste(regulon_DEG,":",algorithm_plot)) + 
  theme(plot.title = element_text(hjust = 0.5))
CTBs_all_network
#pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/merge_CTBs_TEF_TFAP2A_all_network_plot.pdf",width = 6,height=6)
#ggsave(CTBs_all_network,"/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/merge_CTBs_TEF_TFAP2A_all_network_plot.pdf", width=6, height =6)

##for EVTs
TFs_name<-c("ARNT2","ATF3","SMAD3","TEAD1","ELF3","GATA2","HSF1","SNAI1","TFAP2A")

USM_EVTs_DEGs_data1<-USM_CTBs_EVTs_DEGs_data1[which(USM_CTBs_EVTs_DEGs_data1$cluster %in% c("EVTs_1","EVTs_2","EVTs_3")),]
USM_EVTs_DEGs<-unique(as.character(USM_EVTs_DEGs_data1$gene))
TF_gene_filter_EVTs<-TF_gene_ALL[which(TF_gene_ALL$target_final_gene %in% USM_EVTs_DEGs),];dim(TF_gene_filter)#16510     2
target_TF_filter_EVTs<-TF_gene_filter_EVTs[which(TF_gene_filter_EVTs$TFs %in% TFs_name),]
targets_trend<-USM_EVTs_DEGs_data1[which(USM_EVTs_DEGs_data1$gene %in% unique(c(unique(target_TF_filter_EVTs$target_final_gene),TFs_name))),]

####select target_final regulon and corresponding weight
TF_gene_weight1<-TF_gene_weight[which(TF_gene_weight$TF %in% TFs_name),]
head(TF_gene_weight1);dim(TF_gene_weight1)
TF_gene_weight1$pair<-paste0(TF_gene_weight1$TF,"_",TF_gene_weight1$target)
nrow(TF_gene_weight1)#74783     3
head(TF_gene_weight1)

target_TF_gene_final<-target_TF_filter_EVTs
target_TF_gene_final$pair<-paste0(target_TF_gene_final$TFs,"_",target_TF_gene_final$target_final_gene)
head(target_TF_gene_final)
target_TF_gene_final_1<-TF_gene_weight1[which(TF_gene_weight1$pair %in% target_TF_gene_final$pair),]
target_TF_gene_final_1
nrow(target_TF_gene_final_1);nrow(target_TF_gene_final_1)#647
range(target_TF_gene_final_1$importance) ## 0.3865055 147.3493267
table(target_TF_gene_final_1$TF);length(table(target_TF_gene_final_1$target))

##get gene pair link number
TF_gene_final_stat<-data.frame(table(c(target_TF_gene_final_1$TF,target_TF_gene_final_1$target)))
colnames(TF_gene_final_stat)<-c("gene","count")
TF_gene_final_stat<-TF_gene_final_stat[order(TF_gene_final_stat$count,decreasing = T),]
TF_gene_final_stat$type<-ifelse(TF_gene_final_stat$gene %in% unique(target_TF_gene_final_1$TF),"TFs","targets")
head(TF_gene_final_stat)
table(TF_gene_final_stat$type)
TF_gene_final_stat$node<-as.character(TF_gene_final_stat$gene)

##get Trend and cell number of DEGs
head(targets_trend)
targets_trend2<-targets_trend[which(!(targets_trend$gene %in% TFs_name)),]
targets_split <-dcast(targets_trend2, gene ~ trend)
targets_split$Cell_sum<-targets_split$Down+targets_split$Up
targets_split$Up_percent<-targets_split$Up/targets_split$Cell_sum*100
targets_split$DEG_type<- ifelse(targets_split$Down> 0  & targets_split$Up > 0, "Contrary",ifelse(targets_split$Down == 0, "Up_only","Down_only"))
table(targets_split$DEG_type)
#Contrary Down_only   Up_only 
#  10       274       201 

head(targets_split)

##merge for node information
node_info_merge<-merge(TF_gene_final_stat,targets_split,by="gene",all=T)
algorithm_plot<-"nicely";regulon_DEG<-"merge_EVTs_ATF3_TEAD1_GATA2_SNAI1_TFAP2A_HSF1_ELF3_SMAD3_ARNT2"
node_info_merge[which(is.na(node_info_merge$DEG_type)),]$DEG_type<-"TFs"
node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$Cell_sum<-node_info_merge[which(is.na(node_info_merge$Cell_sum)),]$count
node_info_merge[which(node_info_merge$DEG_type=="TFs"),]
node_info_merge$color<-node_info_merge$DEG_type
node_info_merge$color<-ifelse(node_info_merge$color =="Up_only","indianred3",ifelse(node_info_merge$color =="Down_only","steelblue",
                                                                                    ifelse(node_info_merge$color =="Contrary","lightgreen","purple")))
node_info_merge$lable<-node_info_merge$gene
node_info_merge[which(!(node_info_merge$DEG_type =="TFs" | node_info_merge$count>1)),]$lable<-NA

nodes<-node_info_merge
edges<-target_TF_gene_final_1[,c("TF","target","importance")]

net<-graph_from_data_frame(d=edges,vertices = nodes,directed = T)
node_color<-node_info_merge$color

set.seed(19921010)
layout <- create_layout(net, layout ="nicely");algorithm_plot<-"nicely"

#igraph_layouts <- c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl', 'lgl')
head(layout)
EVTs_all_network<-ggraph(layout)+
  geom_edge_fan(aes(edge_width=importance),color="lightblue",alpha=0.7) + #,show.legend=FALSE
  geom_node_point(aes(size=Cell_sum,color=as.factor(DEG_type)))+
  #geom_node_text(aes(filter= DEG_type =="TFs",label=node),size=3)+
  #geom_node_text(aes(filter= (DEG_type =="TFs" | count>1),label=node),size=3)+
  geom_text_repel(aes(x, y, label=lable),force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),segment.color="grey20",segment.size=0.2,segment.alpha=0.8,nudge_y=1,max.overlaps = Inf)+
  scale_color_manual(limits = as.factor(layout$DEG_type), values =node_color)+ 
  #  scale_edge_color_continuous(low = "cyan",high = "red")+
  scale_size_area(max_size = 20)+guides(size="none")+
  scale_edge_width(range = c(0.5, 3))+ 
  guides(color=guide_legend(order=3))+  theme_graph()+ 
  ggtitle(paste(regulon_DEG,":",algorithm_plot)) + 
  theme(plot.title = element_text(hjust = 0.5))
EVTs_all_network
#pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/merge_EVTs_ATF3_TEAD1_GATA2_SNAI1_TFAP2A_HSF1_ELF3_SMAD3_ARNT2_all_network_plot.pdf",width = 6,height=6)

####get the gene expression matrix
#reading expression data
Troph_group<-c("CTBs_1","CTBs_2","EVTs_1","EVTs_2","EVTs_3")#"STBs_1","STBs_2","STBs_3",
Troph_target<-subset(x = target_final, subset = re_annotation_TFac %in% Troph_group)
Troph_target$re_annotation_TFac<- factor(Troph_target$re_annotation_TFac, levels=Troph_group,ordered=TRUE)

##reannotation for CTBs and EVTs
Troph_target$merge_cluster <- as.character(Troph_target$re_annotation_TFac)
sub_cell<-subset(x = Troph_target,subset = re_annotation_TFac %in% c("CTBs_1","CTBs_2"))
Troph_target$merge_cluster[Cells(sub_cell)] <- "CTBs"
sub_cell<-subset(x = Troph_target,subset = re_annotation_TFac %in% c("EVTs_1","EVTs_2","EVTs_3"))
Troph_target$merge_cluster[Cells(sub_cell)] <- "EVTs"
DimPlot(Troph_target, group.by = "merge_cluster",label = F,cols = ppCor_all2)
Troph_target$merge_cluster<- factor(Troph_target$merge_cluster, levels=c("CTBs","EVTs"),ordered=TRUE)

##calculation 
DefaultAssay(Troph_target) <- "RNA"
# Normalize RNA data for visualization purposes
Troph_target <- NormalizeData(Troph_target, verbose = FALSE)
express_data <- as.matrix(GetAssayData(Troph_target, slot = "data"))

##for merge CTBs and EVTs
Group_tag<-"merge_USM_TF_select_0713"
gene_signature<-c("ARNT2","ATF3","SMAD3","TEAD1","ELF3","GATA2","HSF1","SNAI1","TFAP2A","TEF")

##expression matrix for yaxi for figure four
target_exprs <- data.frame(FetchData(object = Troph_target,slot = "data",vars = c("Treat","merge_cluster",gene_signature)))
target_exprs$Treat<-as.character(target_exprs$Treat)
target_exprs[which(target_exprs$Treat =="Abortion"),]$Treat<-"USM"
target_exprs$group<-paste(target_exprs$merge_cluster,target_exprs$Treat,sep="_")
target_exprs<-target_exprs[,-c(1:2)]
target_exprs[1:4,];dim(target_exprs)#29882    11
write.table(target_exprs,file ="/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/normalized_gene_expression_CTBs_EVTs_TFs_regulon_activity.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

###plot the gene expression
Idents(object = Troph_target) <- "merge_cluster"
VlnPlot(Troph_target, features =gene_signature, pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)],ncol = 5)

plot1<-VlnPlot(Troph_target, features = gene_signature[1], y.max =2.6,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot1_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =2.5,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[2], y.max =4.5,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot2_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =4,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[3], y.max =2.6,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot3_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =2.5,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[4], y.max =4,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot4_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =3.5,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[5], y.max =4,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot5_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =3.5,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[6], y.max =3.5,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot6_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =3,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[7], y.max =2.8,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot7_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =2.6,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[8], y.max =2.6,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot8_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =2.5,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[9], y.max =3.1,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot9_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =3,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
plot1<-VlnPlot(Troph_target, features = gene_signature[10], y.max =2,  group.by = "merge_cluster",pt.size = 0, split.by = "Treat",cols = ppCor[c(10,6)], combine =F) 
plot10_1<-plot1[[1]]+ stat_compare_means(method = "wilcox.test",data=plot1$data,label.y =1.7,label = "p.signif",hide.ns = F, label.x = 2)+stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))

gene_list_merge<-grid.arrange(plot1_1,plot2_1,plot3_1,plot4_1,plot5_1,plot6_1,plot7_1,plot8_1,plot9_1,plot10_1,ncol=5)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/merge_USM_TF_select_0713_expression_violin_target_gene_plot_stat_add.pdf"),gene_list_merge,width=15, height=6)

##or manuplot for the expression of target gene 
target_exprs <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/normalized_gene_expression_CTBs_EVTs_TFs_regulon_activity.txt", header=T)
head(target_exprs)
target_exprs$group_2 <- unlist(lapply(strsplit(target_exprs$group,"_"), function(x) x[2]))
target_exprs$cell <- unlist(lapply(strsplit(target_exprs$group,"_"), function(x) x[1]))
head(target_exprs)
##SMAD3和TEAD1
###plot violin plot for activity of regulon
genename<-c("SMAD3","TEAD1")
target_data_plot<-target_exprs[,c(genename,"group_2","cell")]
target_data_plot[target_data_plot==0]<-NA
head(target_data_plot)

target_data_plot2 <- melt(target_data_plot,id=c("group_2","cell"),variable.name="gene",value.name="expression",na.rm = TRUE)
colnames(target_data_plot2)<-c("group","cell","gene","expression")
head(target_data_plot2)
range(target_data_plot2$expression)#0.01782394 3.35799612
table(target_data_plot2$cell,target_data_plot2$group)
plot_Cell_exp<- ggplot(target_data_plot2,aes(x = cell,y=expression,fill = group)) + 
  geom_violin(scale = "width") +#ylim(min(target_data_plot2$expression)*1.05,max(target_data_plot2$expression)*1.05)+
  theme_bw()+scale_fill_manual(values=ppCor[c(6,10)])+
  labs(y = "expression level", x = "Cell_type", title = genename)+theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5),axis.title = element_text(size=10),axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),legend.position = "right")+
  theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))+
  stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
# stat_compare_means(method = "wilcox.test",data=TFs_active_high,aes(x=cell,y = Activity,group = group,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((TFs_active_high$Activity)))
  facet_wrap(~gene,scales = "free_y", ncol = 2) 
plot_Cell_exp
plot_Cell_exp2<-plot_Cell_exp+  stat_compare_means(method = "wilcox.test",data=target_data_plot2,aes(x=cell,y = expression,group = group,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =2.5)
plot_Cell_exp3<-plot_Cell_exp+  stat_compare_means(method = "t.test",data=target_data_plot2,aes(x=cell,y = expression,group = group,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =2.5)

ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/figure4c_final_used_SMAD3_TEAD1_violin_target_regulon_plot_stat_wilcox.pdf",plot_Cell_exp2,width=8, height=4)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/figure4c_final_used_SMAD3_TEAD1_violin_target_regulon_plot_stat_ttest.pdf",plot_Cell_exp3,width=8, height=4)


##calculated H3K27ac activity
TFs_active_high <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/STAT3_k27_activity.txt", header=T)
table(TFs_active_high$ident)
#CTB_A CTB_C EVT_A EVT_C 
# 172   151   774   191
TFs_active_high$group <- unlist(lapply(strsplit(TFs_active_high$ident,"_"), function(x) x[2]))
TFs_active_high$cell <- unlist(lapply(strsplit(TFs_active_high$ident,"_"), function(x) x[1]))
head(TFs_active_high)
###plot violin plot for activity of regulon
genename<-"STAT3"
colnames(TFs_active_high)<-c("Activity","cell_type","group","cell")
range(TFs_active_high$Activity)

plot_Cell_exp<- ggplot(TFs_active_high,aes(x = cell,y=Activity,fill = group)) + 
  geom_violin(scale = "width") + ylim(min(TFs_active_high$Activity)*1.05,max(TFs_active_high$Activity)*1.05)+
  theme_bw()+scale_fill_manual(values=ppCor[c(10,6)])+
  labs(y = "Activity level", x = "Cell_type", title = genename)+theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5),axis.title = element_text(size=10),axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),legend.position = "right")+
  theme(strip.text=element_text(colour = "black", face="bold", size=rel(1)), strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1))+
  stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))
  
# stat_compare_means(method = "wilcox.test",data=TFs_active_high,aes(x=cell,y = Activity,group = group,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((TFs_active_high$Activity)))
plot_Cell_exp2<-plot_Cell_exp+  stat_compare_means(method = "wilcox.test",data=TFs_active_high,aes(x=cell,y = Activity,group = group,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((TFs_active_high$Activity)))
plot_Cell_exp3<-plot_Cell_exp+  stat_compare_means(method = "t.test",data=TFs_active_high,aes(x=cell,y = Activity,group = group,label= ..p.signif..),hide.ns = TRUE, label.x = 2,label.y =max((TFs_active_high$Activity)))


ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/H3K27ac_activity_STAT3_violin_target_regulon_plot_stat_wilcox.pdf",plot_Cell_exp2,width=8, height=4)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/H3K27ac_activity_STAT3_violin_target_regulon_plot_stat_ttest.pdf",plot_Cell_exp3,width=8, height=4)
