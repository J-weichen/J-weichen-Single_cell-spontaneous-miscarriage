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
#reading expression data
#########################
###KS_panc_fl_data
########################
#reading final seurat object
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)
Umap_coordinate<-data.frame(Embeddings(target_final[["umap"]]))
target_final <- AddMetaData(target_final, Umap_coordinate)
DimPlot(target_final, group.by = "re_annotation_TFac",label = F,cols = ppCor_all2)

##substract Trophoblast
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

## Figure 6b plot 
##calculation 
DefaultAssay(Troph_target) <- "RNA"#counts 
# Normalize RNA data for visualization purposes
Troph_target <- NormalizeData(Troph_target, verbose = FALSE)
express_data <- as.matrix(GetAssayData(Troph_target, slot = "data"))

##for pre identified genes
gene_signature<-c("ZEB1", "ZEB2","ZBTB6","TBX18","SOX17","PLAGL1","HIC1")


Idents(object = Troph_target) <- "merge_cluster"
VlnPlot(Troph_target, features =gene_signature, pt.size = 0.01,cols = ppCor)
VlnPlot(Troph_target, features ="ZEB2", pt.size = 0.01,cols = ppCor)
VlnPlot(Troph_target, features ="ZEB1", pt.size = 0.01,cols = ppCor)

Idents(object = Troph_target) <- "re_annotation_TFac"
VlnPlot(Troph_target, features =gene_signature, pt.size = 0.01,cols = ppCor)

VlnPlot(Troph_target, features ="ZEB2", pt.size = 0.01,cols = ppCor)
VlnPlot(Troph_target, features ="ZEB1", pt.size = 0.01,cols = ppCor)

###for final genes
DefaultAssay(Troph_target) <- "RNA"#counts 
Group_tag<-"Twelve_USM_predict_TFs"
gene_signature<-c("ZEB1", "ZEB2","TBX3","TBX1","TBX19","SOX4","SOX13","ZBTB4","ZBTB6","PLAGL1","PLAG1","PLAGL2")
exprs <- data.frame(FetchData(object = Troph_target,slot = "data",vars = c(gene_signature,"Treat","re_annotation_TFac","merge_cluster")))
data_plot <- melt(exprs)
data_plot$re_annotation_TFac<- factor(data_plot$re_annotation_TFac,levels = Troph_group)
head(data_plot);range(data_plot$value)

###
#cell_group<-"merge_cluster";split_class<-"merge_cluster"
cell_group<-"re_annotation_TFac";split_class<-"Treat"

maker_plot<-ggplot(data_plot,aes(x = re_annotation_TFac,y=value,fill = Treat)) +  
  geom_point(size = 1,position = 'jitter',alpha = 0.1)+ 
  geom_violin(aes(fill = Treat),scale = "width",alpha = 0.7) + #geom_violin(scale = "width") +  count area
  facet_wrap(~variable,scales = "free_y", ncol = 3) +
   theme_bw()+ #scale_colour_gradient(low = 'lightblue', high = 'darkblue')
   scale_fill_manual(values=ppCor[c(3,9)])+  
  labs(x = "Group", y = "LogNormalized count", title ="SNPs correlated TFs")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x.bottom = element_blank(),legend.position = "right")
maker_plot
#ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_contain_cell_",Group_tag,"_violin_target_gene_plot.pdf"),maker_plot,width=9, height=12)
plot1_13<-plot1[[13]]+ stat_compare_means(method = "wilcox.test",data=plot1[[13]]$data,label = "p.signif",hide.ns = F, label.x = 2,label.y =2)+geom_boxplot(width=0.15,outlier.shape = NA,notch=T,position=position_dodge(0.9),color="grey", alpha=0.4) +stat_summary(fun=mean, geom="point", shape=18, size=2, col="orange", position = position_dodge(0.9))

ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_contain_cell_sub_and_split_",Group_tag,"_violin_target_gene_plot.pdf"),maker_plot,width=12, height=12)

##remove zero expression
data_plot2<-data_plot[which(data_plot$value>0),]
#maker_plot<-ggplot(data_plot2,aes(x = merge_cluster,y=value,fill = merge_cluster)) +
maker_plot<-ggplot(data_plot2,aes(x = re_annotation_TFac,y=value,fill = Treat)) +  
 # geom_point(size = 1,position = 'jitter',alpha = 0.1)+  
  geom_violin(aes(fill = Treat),scale = "width",alpha = 0.7) +  
  facet_wrap(~variable,scales = "free_y", ncol = 3) +
  scale_fill_manual(values=ppCor[c(3,9)])+  theme_bw()+ 
  labs(x = "Group", y = "LogNormalized count", title ="SNPs correlated TFs")+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x.bottom = element_blank(),legend.position = "right")+
  stat_summary(fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))

maker_plot
#ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_",Group_tag,"_violin_target_gene_plot.pdf"),maker_plot,width=9, height=12)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_sub_and_split_",Group_tag,"_violin_target_gene_plot.pdf"),maker_plot,width=12, height=12)

##add p value
gene_signature<-c("ZEB1", "ZEB2","TBX3","TBX1","TBX19","SOX4","SOX13","ZBTB4","ZBTB6","PLAGL1","PLAG1","PLAGL2")
##ZEB1，ZEB2，PLAGL1，TBX3，ZBTB6，SOX4
##For SOX4
genename <-"SOX4"
print(genename)
merge_data_target<- data_plot2[which(data_plot2$variable==genename),]
head(merge_data_target);dim(merge_data_target)#17723     5
merge_data_target_original<-merge_data_target
range(merge_data_target_original$value)# 0.1479054 4.6395716
max_value<-4.6
merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot3<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  scale_y_continuous(limits = c(0,4.7),breaks = c(0,1,2,3,4,5))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot3
#ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_sub_and_split_max_value_adjust_",genename,"_violin_target_gene_plot.pdf"),maker_plot3,width=5, height=4)

##For ZBTB6
genename <-"ZBTB6"
print(genename)
merge_data_target<- data_plot2[which(data_plot2$variable==genename),]
head(merge_data_target);dim(merge_data_target)#263   5
merge_data_target_original<-merge_data_target
range(merge_data_target_original$value)# 0.03473989 1.90326545
max_value<-2
#merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot3<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  scale_y_continuous(limits = c(0,2),breaks = c(0,0.5,1,1.5,2))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot3
#ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_sub_and_split_max_value_adjust_",genename,"_violin_target_gene_plot.pdf"),maker_plot3,width=5, height=4)


##For TBX3
genename <-"TBX3"
print(genename)
merge_data_target<- data_plot2[which(data_plot2$variable==genename),]
head(merge_data_target);dim(merge_data_target)#263   5
merge_data_target_original<-merge_data_target
range(merge_data_target_original$value)# 0.01782394 2.90499725
max_value<-2.9
merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot3<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  scale_y_continuous(limits = c(0,3),breaks = c(0,0.5,1,1.5,2,2.5,3))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot3
#ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_sub_and_split_max_value_adjust_",genename,"_violin_target_gene_plot.pdf"),maker_plot3,width=5, height=4)


##For PLAGL1
genename <-"PLAGL1"
print(genename)
merge_data_target<- data_plot2[which(data_plot2$variable==genename),]
head(merge_data_target);dim(merge_data_target)#263   5
merge_data_target_original<-merge_data_target
range(merge_data_target_original$value)#0.04635717 2.28346065
max_value<-2.2
merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot3<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  scale_y_continuous(limits = c(0,2.5),breaks = c(0,0.5,1,1.5,2,2.5))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot3
#ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_sub_and_split_max_value_adjust_",genename,"_violin_target_gene_plot.pdf"),maker_plot3,width=5, height=4)



##For ZEB2
genename <-"ZEB2"
print(genename)
merge_data_target<- data_plot2[which(data_plot2$variable==genename),]
head(merge_data_target);dim(merge_data_target)#263   5
merge_data_target_original<-merge_data_target
range(merge_data_target_original$value)#0.01752079 2.88355061
max_value<-2.8
merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot3<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  scale_y_continuous(limits = c(0,3),breaks = c(0,0.5,1,1.5,2,2.5,3))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot3
#ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_sub_and_split_max_value_adjust_",genename,"_violin_target_gene_plot.pdf"),maker_plot3,width=5, height=4)


##For ZEB1
genename <-"ZEB1"
print(genename)
merge_data_target<- data_plot2[which(data_plot2$variable==genename),]
head(merge_data_target);dim(merge_data_target)#263   5
merge_data_target_original<-merge_data_target
max_value<-1.6
merge_data_target[which(merge_data_target$value> max_value),]$value<-max_value
maker_plot3<-ggplot(merge_data_target,aes(x = re_annotation_TFac,y=value,fill = Treat)) + geom_violin(scale = "width",trim = TRUE) + 
  stat_summary(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
  stat_compare_means(data=merge_data_target_original,mapping =aes(x = re_annotation_TFac,y=value), method = "wilcox.test",label.y =max_value,label = "p.signif",hide.ns = F, label.x = 1.5)+
  scale_y_continuous(limits = c(0,1.8),breaks = c(0,0.5,1,1.5))+
  scale_fill_manual(values=ppCor[c(10,6)])+ 
  #   stat_compare_means(method = "wilcox.test",label.y =max(merge_data_target$value),label = "p.signif",hide.ns = F, label.x = 1.5)+
  labs(x = "Group", y = "LogNormalized count", title =genename)+
  theme_bw()+ theme( panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),legend.position = "right")
maker_plot3
#ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Three_USM_SNP_overlapped_max_value_adjust_",genename,"_violin_target_gene_plot_stat.pdf"),maker_plot,width=5, height=4)
ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/zero_remove_cell_sub_and_split_max_value_adjust_",genename,"_violin_target_gene_plot.pdf"),maker_plot3,width=5, height=4)



rownames(Troph_target)[grep("ZEB",rownames(Troph_target))]# "ZEB2"     "ZEB1-AS1" "ZEB1"     "ZEB2-AS1"

##for TBX
gene_signature<-rownames(Troph_target)[grep("TBX",rownames(Troph_target))]# "TBX19"    "TBX20"    "TBXAS1"   "TBX5"     "TBX3"     "TBX21"    "TBX2-AS1" "TBX2"     "TBXA2R"   "TBX1"     "TBX15"    "TBX18"  "TBX5-AS1" "TBX6"     "TBX4"     "TBX10"    "TBXT"    
DotPlot(Troph_target, features = gene_signature,group.by ="re_annotation_TFac",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()
DotPlot(Troph_target, features = gene_signature,group.by ="merge_cluster",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()

##for PLAG
gene_signature<-rownames(Troph_target)[grep("PLAG",rownames(Troph_target))]#"PLAGL1" "PLAG1"  "PLAGL2"
DotPlot(Troph_target, features = gene_signature,group.by ="re_annotation_TFac",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()
DotPlot(Troph_target, features = gene_signature,group.by ="merge_cluster",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()

#SOX17
gene_signature<-rownames(Troph_target)[grep("SOX",rownames(Troph_target))]#"QSOX1"     "SOX13"     "SOX14"     "SOX30"     "SOX4"      "SOX7"      "QSOX2"     "SOX6"      "SOX15"     "SOX12"     "SOX2" "SOX17"     "SOX5"      "SOX8"      "SOX18"     "SOX11"     "SOX9"      "SOX21"     "SOX21-AS1" "SOX1"      
DotPlot(Troph_target, features = gene_signature,group.by ="re_annotation_TFac",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()
DotPlot(Troph_target, features = gene_signature,group.by ="merge_cluster",cols =  c("blue","red"),assay = "RNA")# + RotatedAxis()

#ZBTB
gene_signature<-rownames(Troph_target)[grep("ZBTB",rownames(Troph_target))]# 48 gene  
DotPlot(Troph_target, features = gene_signature,group.by ="re_annotation_TFac",cols =  c("blue","red"),assay = "RNA")+ RotatedAxis()
DotPlot(Troph_target, features = gene_signature,group.by ="merge_cluster",cols =  c("blue","red"),assay = "RNA") + RotatedAxis()

##plot dot
Group_tag<-"Twelve_USM_predict_TFs"
gene_signature<-c("ZEB1", "ZEB2","TBX3","TBX1","TBX19","SOX4","SOX13","ZBTB4","ZBTB6","PLAGL1","PLAG1","PLAGL2")

DefaultAssay(Troph_target) <- "SCT" #counts ##https://www.jianshu.com/p/21b2f08652ed corrected count
exprs <- data.frame(FetchData(object = Troph_target,slot = "counts",vars = c(gene_signature,"Treat","re_annotation_TFac","merge_cluster")))
data_plot <- melt(exprs)
head(data_plot);range(data_plot$value)# 0 157

head(data_plot);range(data_plot$value)
data_plot$value2<-data_plot$value
data_plot[which(data_plot$value !=0),]$value2<-1
##calculated expression level
data_preplot<-aggregate(value ~  merge_cluster + variable, data = data_plot, mean)#median
head(data_preplot)
data_preplot$expression<-log2(data_preplot$value+1)
data_preplot2<-dcast(data_plot,merge_cluster + variable~ value2)
data_preplot2$percent<-round((data_preplot2$`1`/(data_preplot2$`0`+data_preplot2$`1`))*100,2)

data_preplot$cell_class<-paste(data_preplot$merge_cluster,data_preplot$variable,sep="-")
data_preplot2$cell_class<-paste(data_preplot2$merge_cluster,data_preplot2$variable,sep="-")
data_preplot3<-merge(data_preplot,data_preplot2[,c("cell_class","percent")])

data_preplot3$merge_cluster<-factor(data_preplot3$merge_cluster,levels = c("CTBs","EVTs"))
data_preplot3$variable<-factor(data_preplot3$variable,levels =gene_signature)

head(data_preplot3)
range(data_preplot3$expression)# 0.0009566655 3.9667915089
Gene_expression_plot<-ggplot(data_preplot3,aes(x=merge_cluster,y=variable,size=percent,colour=expression))+
  geom_point(alpha =0.8,na.rm = TRUE)+
  scale_size(breaks = c(1,10,20,40,60,80,100),range = c(1,6),name='Postive Percentage')+
  scale_color_gradient("Log2(corrected count+1)", low = "darkred", high = "yellow", limits=c(0.0001,4))+#, breaks = c(0.00001,0.001,0.01,0.1,0.5,1,1.5,2,2.1,2.2))
  #scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,0.2,0.5,1,1.2,1.5,2,2.5),name='LogNormalized count')+  
  theme_classic()+labs(x="",y="",title="Target average gene expression")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),
        legend.position ="right",legend.direction = "vertical")
Gene_expression_plot

ggsave(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/Manu_",Group_tag,"_genes_dotplot_all_celltype_split.pdf"),Gene_expression_plot,width=6, height=6)

