rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(Seurat)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ComplexHeatmap)
library(pheatmap) #https://www.cnblogs.com/Mao1518202/p/11269595.html
library(scales)
library(circlize)
library(ggsci)
library(psych)
library(pheatmap)
library(UpSetR)
library(DESeq2)
library(VennDiagram)
library(grid)
library(futile.logger)

pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
show_col(pal)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
show_col(nejm)
ppCor <-c(nejm,pal[c(2,5,9)])

#step1 ：read Seurat object and load initial coldata and matrix
target_final<-readRDS(file = "/mnt/data/chenwei/gongchen/4.SCENIC/reannotation_file/Target_final_TFs_regulon_activity_seurat_object_new.rds")
subgroup_order0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
major_order0<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN","Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
target_final$re_annotation_TFac<- factor(target_final$re_annotation_TFac, levels=subgroup_order0,ordered=TRUE)
target_final$re_anno_TFac_major<- factor(target_final$re_anno_TFac_major,levels=major_order0,ordered=TRUE)

cell_prop<-as.data.frame(table(target_final$re_annotation_TFac))
colnames(cell_prop)<-c("celltype","Cell_number")
cell_prop<-cell_prop[which(cell_prop$Cell_number >0),]
cell_prop$proportion<-cell_prop$Cell_number *100/sum(cell_prop$Cell_number)
colnames(cell_prop)<-c("celltype","Cell_number","proportion")
head(cell_prop);dim(cell_prop);range(cell_prop$Cell_number);range(cell_prop$proportion)
#19 3 # 37 8746 # 0.08078603 19.09606987

#important step selection
#cell_exclude<-c("Troph_mix_group","Erythrocyte_HBA1","Mast_Cell_MS4A2","B_Cell_CD79A")
#cell_prop_select<-cell_prop
Used_cell<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Mycs","HCs","STCs","FBs","Endo","NKs","Ts")
cell_prop_select_2<-cell_prop[which(cell_prop$celltype %in% Used_cell),]
cell_prop_select_2[order(cell_prop_select_2$proportion,decreasing = T),]
dim(cell_prop_select_2)#15  3

#
CTRL<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/CTRL/significant_means.txt",head = T,sep = "\t")
Abortion<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/Abortion/significant_means.txt",head = T,sep = "\t")
str(CTRL);str(Abortion);dim(CTRL);dim(Abortion)#[1]1199  373#[1]1199  373

Cell_pair_CTRL1<-unlist(lapply(strsplit(as.character(colnames(CTRL)),"[.]"), function(x) x[1]))
Cell_pair_CTRL2<-unlist(lapply(strsplit(as.character(colnames(CTRL)),"[.]"), function(x) x[2]))
Cell_pair_Abortion1<-unlist(lapply(strsplit(as.character(colnames(Abortion)),"[.]"), function(x) x[1]))
Cell_pair_Abortion2<-unlist(lapply(strsplit(as.character(colnames(Abortion)),"[.]"), function(x) x[2]))

unique(c(Cell_pair_CTRL1,Cell_pair_CTRL2)); unique(c(Cell_pair_Abortion1,Cell_pair_Abortion2))


CTRL_rm<-CTRL[,c(1:12,which(Cell_pair_CTRL1 %in% Used_cell & Cell_pair_CTRL2 %in% Used_cell))]
Abortion_rm<-Abortion[,c(1:12,which(Cell_pair_Abortion1 %in% Used_cell & Cell_pair_Abortion2 %in% Used_cell))]
CTRL[1:6,1:6];Abortion_rm[1:6,1:6]
dim(CTRL_rm);dim(Abortion_rm)#1199  237# 1199  237

#merge CTRL and Abortion
CTRL_rm[1:6,1:13]
Abortion1<-Abortion_rm[,c(1,13:ncol(Abortion_rm))]
CTRL1<-CTRL_rm[,c(1,13:ncol(CTRL_rm))]
Abortion1[1:6,1:6];CTRL1[1:6,1:6]
rownames(Abortion1)<-Abortion1$id_cp_interaction;Abortion1<-Abortion1[,-1]
rownames(CTRL1)<-CTRL1$id_cp_interaction;CTRL1<-CTRL1[,-1]
Abortion1[!(is.na(Abortion1))] <- 1;Abortion1[is.na(Abortion1)] <- 0
CTRL1[!(is.na(CTRL1))] <- 1;CTRL1[is.na(CTRL1)] <- 0

#for number stastic
dim(Abortion1);dim(CTRL1)#1199    225 # 1199    225
Abortion_st<-Abortion1;CTRL_st<-CTRL1
upset_rowSums<-rowSums(Abortion_st);range(upset_rowSums)#  0 91
upset_colSums<-colSums(Abortion_st);range(upset_colSums)#  1 163

ifelse(length(which(rowSums(Abortion_st)==0))==0,Abortion_st<-Abortion_st,Abortion_st<-Abortion_st[-which(rowSums(Abortion_st)== 0),])
ifelse(length(which(colSums(Abortion_st)==0))==0,Abortion_st<-Abortion_st,Abortion_st<-Abortion_st[,-which(colSums(Abortion_st)== 0)])
ifelse(length(which(rowSums(CTRL_st)==0))==0,CTRL_st<-CTRL_st,CTRL_st<-CTRL_st[-which(rowSums(CTRL_st)== 0),])
ifelse(length(which(colSums(CTRL_st)==0))==0,CTRL_st<-CTRL_st,CTRL_st<-CTRL_st[,-which(colSums(CTRL_st)== 0)])
dim(Abortion_st);dim(CTRL_st)#633 225#615 225
Abortion_st[1:6,1:6];CTRL_st[1:6,1:6]

grid.newpage();
grid.draw(venn.diagram(list(CTRL=rownames(CTRL_st),Abortion=rownames(Abortion_st)),fill=c("blue","red"),main="Cell_communication_molecular_pair", filename = NULL))
grid.newpage(); 
grid.draw(venn.diagram(list(CTRL=colnames(CTRL_st),Abortion=colnames(Abortion_st)),fill=c(ppCor[4],ppCor[3]),main="communication_Cell_pair",  filename = NULL))

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/common_molecular_or_cell_pair.pdf")
grid.draw(venn.diagram(list(CTRL=rownames(CTRL_st),Abortion=rownames(Abortion_st)),fill=c("blue","red"),main="Cell_communication_molecular_pair", filename = NULL))
grid.draw(venn.diagram(list(CTRL=colnames(CTRL_st),Abortion=colnames(Abortion_st)),fill=c(ppCor[4],ppCor[3]),main="communication_Cell_pair",  filename = NULL))
dev.off()

length(intersect(rownames(CTRL_st),rownames(Abortion_st)))#   582
length(union(rownames(CTRL_st),rownames(Abortion_st)))# 666
setdiff(rownames(CTRL_st),rownames(Abortion_st))

CTRL_special_MP<-setdiff(rownames(CTRL_st),rownames(Abortion_st))
Abortion_special_MP<-setdiff(rownames(Abortion_st),rownames(CTRL_st))

write.table(as.data.frame(CTRL_special_MP),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/CTRL_special_MP.csv',row.names = F,sep = ",")
write.table(as.data.frame(Abortion_special_MP),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_special_MP.csv',row.names = F,sep = ",")

#for indentify molecular pair name
pre_enrich_name<-merge(Abortion_rm[,1:2], CTRL_rm[,1:2], all = T, by = "id_cp_interaction")
pre_enrich_name$interacting_pair.x<-as.character(pre_enrich_name$interacting_pair.x)
pre_enrich_name$interacting_pair.y<-as.character(pre_enrich_name$interacting_pair.y)
pre_enrich_name$interacting_pair<-ifelse(is.na(pre_enrich_name$interacting_pair.x),pre_enrich_name$interacting_pair.y, pre_enrich_name$interacting_pair.x)
head(pre_enrich_name);dim(pre_enrich_name);which(is.na(pre_enrich_name$interacting_pair))
# 1199    4

CTRL_Mpair<-pre_enrich_name[which(pre_enrich_name$id_cp_interaction %in% CTRL_special_MP),]$interacting_pair
CTRL_pair_pre_enrich_genes<-unique(c(unlist(lapply(strsplit(as.character(CTRL_Mpair),"_"), function(x) x[1])),unlist(lapply(strsplit(as.character(CTRL_Mpair),"_"), function(x) x[2]))))

Abortion_Mpair<-pre_enrich_name[which(pre_enrich_name$id_cp_interaction %in% Abortion_special_MP),]$interacting_pair
Abortion_pair_pre_enrich_genes<-unique(c(unlist(lapply(strsplit(as.character(Abortion_Mpair),"_"), function(x) x[1])),unlist(lapply(strsplit(as.character(Abortion_Mpair),"_"), function(x) x[2]))))

write.table(as.data.frame(CTRL_pair_pre_enrich_genes),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/CTRL_pair_pre_enrich_genes_list.csv',row.names = F,sep = ",")
write.table(as.data.frame(Abortion_pair_pre_enrich_genes),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_pair_pre_enrich_genes_list.csv',row.names = F,sep = ",")

write.table(as.data.frame(pre_enrich_name[which(pre_enrich_name$id_cp_interaction %in% CTRL_special_MP),]),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/CTRL_pair_special_genes_pair.csv',row.names = F,sep = ",")
write.table(as.data.frame(pre_enrich_name[which(pre_enrich_name$id_cp_interaction %in% Abortion_special_MP),]),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_pair_special_genes_pair.csv',row.names = F,sep = ",")

##statistic two
merge_data<-merge(Abortion_st, CTRL_st, all = T, by=0)
dim(merge_data);dim(Abortion_st);dim(CTRL_st)#666 451 #633 225 #615 225
head(merge_data[,1:6]);tail(merge_data[,1:6])
rownames(merge_data)<-merge_data$Row.names;merge_data<-merge_data[,-1]
plot_num1<-ggplot(data =as.data.frame(table(rowSums(merge_data,na.rm=TRUE)))) + 
  geom_bar(aes(x =  Var1, y = Freq), stat = "identity",fill="#69b3a2",color="#e9ecef", alpha=0.9)+
  theme_bw()+labs(x="Cell_pairs_number",y="molecule_type_number")
plot_num2<-ggplot(data =as.data.frame(table(colSums(merge_data,na.rm=TRUE)))) + 
  geom_bar(aes(x =  Var1, y = Freq), stat = "identity",fill="orange",color="#e9ecef", alpha=0.9)+
  theme_bw()+labs(x="molecule_type_number",y="Cell_pairs_number")
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/relationship_molecular_or_cell_pair.pdf")
plot_num1+plot_num2
dev.off()

#for calculate next
Abortion1[1:6,1:6];CTRL1[1:6,1:6]
npair_Abortion<-ncol(Abortion1);npair_CTRL<-ncol(CTRL1)
colnames(Abortion1)<-paste0(colnames(Abortion1),".Abortion");colnames(CTRL1)<-paste0(colnames(CTRL1),".CTRL")

merge_data<-merge(Abortion1, CTRL1, all = T, by=0)
dim(merge_data);dim(Abortion1);dim(CTRL1)#1199   451 #1199   225 #1199    225
head(merge_data[,1:6]);tail(merge_data[,1:6])
rownames(merge_data)<-merge_data$Row.names;merge_data<-merge_data[,-1]

merge_data1<-merge_data
merge_data2<-as.data.frame(t(merge_data1))
merge_data2[1:6,1:6];dim(merge_data2)# 450 1199
merge_data2[is.na(merge_data2)] <- 0;range(merge_data2)# 0 1

upset_rowSums<-rowSums(merge_data2);range(upset_rowSums)# 1 189
upset_colSums<-colSums(merge_data2);range(upset_colSums)#  0 200
upset_colSums[which(upset_colSums== 200)]#CPI-SS034D36D2F 

merge_data_write<-merge_data2
merge_data_write$summ<-as.numeric(upset_rowSums)

which(colSums(merge_data2)== 0)
ifelse(length(which(colSums(merge_data2)==0))==0,merge_data2<-merge_data2,merge_data2<-merge_data2[,-which(colSums(merge_data2)== 0)])
head(merge_data2);dim(merge_data2)# 450 666
colnames(merge_data2)
range(colSums(merge_data2));range(rowSums(merge_data2))#1 200#1 189
length(which(rowSums(merge_data2)== 0));length(which(rowSums(merge_data2)>10)) #0; 410

dim(merge_data2)#450 666
cell_type_pair<-rownames(merge_data2);length(cell_type_pair) # 450
#unlist(lapply(strsplit(rownames(merge_data2), "[.]"), function(x) paste(x[1],x[2],sep="-")))
annotation_cell_type = data.frame(origin =c(rep("Abortion",npair_Abortion),rep("CTRL",npair_CTRL)))
rownames(annotation_cell_type) = cell_type_pair

merge_data3<-t(merge_data2)
dim(merge_data3) ## 666 450
range(merge_data3)
merge_data3[1:6,1:6]

anno_colors = list(origin=c(Abortion =ppCor[1],CTRL =ppCor[2]))
labels_col = c("")
p1<-pheatmap(merge_data3,cluster_rows=T, cluster_cols =F, show_rownames=F,show_colnames=F,
         annotation_col =annotation_cell_type,labels_row  = labels_col,
         annotation_colors = anno_colors,
         gaps_col = npair_Abortion,cutree_col = 2,
         main = "inter_pair",
         legend_breaks = c(0.2,0.8),legend_labels = c("unfunc","func"),color = c(ppCor[6],ppCor[3]))

p2<-pheatmap(merge_data3,cluster_rows=T, cluster_cols =T, show_rownames=F,show_colnames=F,
         annotation_col =annotation_cell_type,labels_row  = labels_col,
         annotation_colors = anno_colors,main = "inter_pair",
         legend_breaks = c(0.2,0.8),legend_labels = c("unfunc","func"),color = c(ppCor[6],ppCor[3]))

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/heatmap_interaction_between_Abortion_and_CTRL_1.pdf")
print(p1)
dev.off()
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/heatmap_interaction_between_Abortion_and_CTRL_2.pdf")
print(p2)
dev.off()

##1.prepare file for Cytocyte
merge_data2[1:6,1:6];dim(merge_data2)#450 666
merge_data_cyto<-merge_data2

merge_data_cyto$RowS<-rowSums(merge_data_cyto)
merge_data_cyto$RowS<-as.numeric(merge_data_cyto$RowS)

merge_data_cyto$Treat_group<-c(rep("Abortion",npair_Abortion),rep("CTRL",npair_CTRL))
merge_data_cyto$Cell_pair<-unlist(lapply(strsplit(as.character(rownames(merge_data_cyto)),"[.]"), function(x) paste(x[1],x[2],sep="-")))
dim(merge_data_cyto)# 450 669

merge_data_cyto1<-dcast(merge_data_cyto[,c("Treat_group","Cell_pair","RowS")], Cell_pair~Treat_group)
merge_data_cyto1[is.na(merge_data_cyto1)]<-0
head(merge_data_cyto1);dim(merge_data_cyto1)#[1] 256   3
merge_data_cyto2<-merge_data_cyto1[which(!(merge_data_cyto1$Abortion==0 & merge_data_cyto1$CTRL==0)),]
dim(merge_data_cyto2)#[1] 225   3
merge_data_cyto2$Count<- merge_data_cyto2$Abortion-merge_data_cyto2$CTRL 
range(merge_data_cyto2$Count)#-71  24

merge_data_cyto2$Ligand_celltype= unlist(strsplit(as.character(merge_data_cyto2$Cell_pair) ,split='-') )[seq(1,by=2,length.out = nrow(merge_data_cyto2))] 
merge_data_cyto2$Receptor_celltype=unlist(strsplit(as.character(merge_data_cyto2$Cell_pair) ,split='-'))[seq(2,by=2,length.out = nrow(merge_data_cyto2))]

merge_data_cyto2$Regulated<-'none'
merge_data_cyto2[merge_data_cyto2$Count>0,]$Regulated<-'up'
merge_data_cyto2[merge_data_cyto2$Count<0,]$Regulated<-'down'
merge_data_cyto2$Count<-abs(merge_data_cyto2$Count)
head(merge_data_cyto2)

#for direction of interaction
#Cell_type0<-as.character(merge_data_cyto2$Cell_pair)
#Cell_type0<-Cell_type0[order(Cell_type0)]
Cell_type<-unique(c(unlist(lapply(strsplit(as.character(merge_data_cyto2$Cell_pair),"-"), function(x) x[1])),unlist(lapply(strsplit(as.character(merge_data_cyto2$Cell_pair),"-"), function(x) x[2]))))
Cell_type<-Cell_type[order(Cell_type)]
Cell_type2<-as.data.frame(t(combn(Cell_type,2)))
merge_data_cyto2$Dir<-'self'
merge_data_cyto2[which(merge_data_cyto2$Cell_pair %in% paste(Cell_type2$V1,Cell_type2$V2,sep="-")),]$Dir<-'dir_a'
merge_data_cyto2[which(merge_data_cyto2$Cell_pair %in% paste(Cell_type2$V2,Cell_type2$V1,sep="-")),]$Dir<-'dir_b'

merge_data_cyto2<-merge_data_cyto2[,c('Ligand_celltype','Receptor_celltype','Count','Regulated','Dir')]
head(merge_data_cyto2)

tapply.df1=data.frame(celltype=names(tapply(abs(merge_data_cyto2$Count),merge_data_cyto2$Ligand_celltype, sum)),L_Total_Count=unname(tapply(abs(merge_data_cyto2$Count) ,merge_data_cyto2$Ligand_celltype, sum)))
tapply.df2=data.frame(celltype=names(tapply(abs( merge_data_cyto2$Count),merge_data_cyto2$Receptor_celltype, sum)),R_Total_Count=unname(tapply(abs(merge_data_cyto2$Count) ,merge_data_cyto2$Receptor_celltype, sum)))
stat.df=merge(tapply.df1,tapply.df2,by='celltype',all=T)
stat.df$Total_Count=rowSums(stat.df[,2:3],na.rm = T)
head(stat.df)
stat.df2=merge(stat.df,cell_prop_select_2,by='celltype',all=T)
head(stat.df2);dim(stat.df2)
#install.packages("xlsx")
#library(xlsx)
saveroad<-"/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/"
sample_tag<-"Cytocyte_select_Cell"
write.table(as.data.frame(merge_data_cyto2),paste(saveroad,sample_tag,'_Abortion_CTRL_network.txt',sep = ''),row.names = F,sep = "\t")
write.table(as.data.frame(stat.df2),paste(saveroad,sample_tag,'_Abortion_CTRL_stat.txt',sep = ''),row.names = F,sep = "\t")

head(merge_data_cyto2)
head(stat.df2)

#2. dotplot for inter_pair
Abortion_mean<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/Abortion/means.txt",head = T,sep = "\t")
CTRL_mean<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/CTRL/means.txt",head = T,sep = "\t")
Abortion_pvalue<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/Abortion/pvalues.txt",head = T,sep = "\t")
CTRL_pvalue<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/CTRL/pvalues.txt",head = T,sep = "\t")
dim(CTRL_pvalue);dim(Abortion_pvalue)##1199  372 ##1199  372
dim(CTRL_mean);dim(Abortion_mean)##1199  372 ##1199  372
Abortion_mean[1:6,10:16];Abortion_pvalue[1:6,10:16]

Cell_pair_CTRL1<-unlist(lapply(strsplit(as.character(colnames(CTRL_mean)),"[.]"), function(x) x[1]))
Cell_pair_CTRL2<-unlist(lapply(strsplit(as.character(colnames(CTRL_mean)),"[.]"), function(x) x[2]))
CTRL_mean_rm<-CTRL_mean[,c(1:11,which(Cell_pair_CTRL1 %in% Used_cell & Cell_pair_CTRL2 %in% Used_cell))]
CTRL_mean[10:16,10:16];CTRL_mean_rm[10:16,10:16]

Cell_pair_Abortion1<-unlist(lapply(strsplit(as.character(colnames(Abortion_mean)),"[.]"), function(x) x[1]))
Cell_pair_Abortion2<-unlist(lapply(strsplit(as.character(colnames(Abortion_mean)),"[.]"), function(x) x[2]))
Abortion_mean_rm<-Abortion_mean[,c(1:11,which(Cell_pair_Abortion1 %in% Used_cell & Cell_pair_Abortion2 %in% Used_cell))]
Abortion_rm[1:6,1:6];Abortion_mean_rm[1:16,1:16]
dim(CTRL_mean_rm);dim(Abortion_mean_rm)#1199    236 #1199    236

Cell_pair_CTRL1<-unlist(lapply(strsplit(as.character(colnames(CTRL_pvalue)),"[.]"), function(x) x[1]))
Cell_pair_CTRL2<-unlist(lapply(strsplit(as.character(colnames(CTRL_pvalue)),"[.]"), function(x) x[2]))
CTRL_pvalue_rm<-CTRL_pvalue[,c(1:11,which(Cell_pair_CTRL1 %in% Used_cell & Cell_pair_CTRL2 %in% Used_cell))]
CTRL_pvalue[10:16,10:16];CTRL_pvalue_rm[10:16,10:16]

Cell_pair_Abortion1<-unlist(lapply(strsplit(as.character(colnames(Abortion_pvalue)),"[.]"), function(x) x[1]))
Cell_pair_Abortion2<-unlist(lapply(strsplit(as.character(colnames(Abortion_pvalue)),"[.]"), function(x) x[2]))
Abortion_pvalue_rm<-Abortion_pvalue[,c(1:11,which(Cell_pair_Abortion1 %in% Used_cell & Cell_pair_Abortion2 %in% Used_cell))]
Abortion_rm[1:6,1:6];Abortion_pvalue_rm[1:16,1:16]
dim(CTRL_pvalue_rm);dim(Abortion_pvalue_rm)#1199  236 #1199  236

CTRL_mean_new<-CTRL_mean_rm[which(CTRL_mean_rm$id_cp_interaction %in% colnames(merge_data2)),c(1,12:ncol(CTRL_mean_rm))]
Abortion_mean_new<-Abortion_mean_rm[which(Abortion_mean_rm$id_cp_interaction %in% colnames(merge_data2)),c(1,12:ncol(Abortion_mean_rm))]
CTRL_pvalue_new<-CTRL_pvalue_rm[which(CTRL_pvalue_rm$id_cp_interaction %in% colnames(merge_data2)),c(1,12:ncol(CTRL_pvalue_rm))]
Abortion_pvalue_new<-Abortion_pvalue_rm[which(Abortion_pvalue_rm$id_cp_interaction %in% colnames(merge_data2)),c(1,12:ncol(Abortion_pvalue_rm))]
Abortion_mean_new[1:6,1:6];Abortion_pvalue_new[1:6,1:6]

merge_mean_data<-merge(Abortion_mean_new, CTRL_mean_new, all = T, by = "id_cp_interaction")
dim(merge_mean_data);dim(Abortion_mean_new);dim(CTRL_mean_new)#666 451 #666 226#666 226
head(merge_mean_data[,1:6]);tail(merge_mean_data[,1:6])
merge_pvalue_data<-merge(Abortion_pvalue_new, CTRL_pvalue_new, all = T, by = "id_cp_interaction")
dim(merge_pvalue_data);dim(Abortion_pvalue_new);dim(CTRL_pvalue_new)#666 451 #666 226#666 226
head(merge_pvalue_data[,1:6]);tail(merge_pvalue_data[,1:6])

head(melt(merge_mean_data))
head(melt(merge_pvalue_data))
merge_mean_melt<-melt(merge_mean_data)
merge_pvalue_melt<-melt(merge_pvalue_data)
colnames(merge_mean_melt)<-c("id_cp_interaction","variable","mean")
colnames(merge_pvalue_melt)<-c("id_cp_interaction","variable","pvalue")

merge_mean_pvalue<-merge(merge_mean_melt, merge_pvalue_melt, by = c("id_cp_interaction","variable"),all = T)
head(merge_mean_pvalue);dim(merge_mean_pvalue)

merge_mean_pvalue$Treat_group<-unlist(lapply(strsplit(as.character(merge_mean_pvalue$variable),"[.]"), function(x) x[3]))
merge_mean_pvalue$Treat_group<-as.factor(ifelse(grepl("x",merge_mean_pvalue$Treat_group), "Abortion", "CTRL"))
merge_mean_pvalue$Cell_pair<-unlist(lapply(strsplit(as.character(merge_mean_pvalue$variable),"[.]"), function(x) paste(x[1],x[2],sep="-")))
merge_mean_pvalue$plogtran<-c(-log((merge_mean_pvalue$pvalue+0.00001),10))
head(merge_mean_pvalue)

summary(merge_mean_pvalue$plogtran);summary(merge_mean_pvalue$mean)
#   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.000004 -0.000004 -0.000004  0.373954 -0.000004  5.000000 
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.0000   0.0060   0.0630   0.9646   0.3510 433.8970 

#mtcars2 <- transform(merge_mean_pvalue, plogtran = ifelse(pvalue >= 0.05, NA, plogtran))

merge_mean_pvalue2 <- transform(merge_mean_pvalue, plogtran = ifelse(pvalue >= 0.05, NA, plogtran), mean = ifelse(mean ==0, NA, mean))
head(merge_mean_pvalue2)

## for all cell with exist inter_pairs
whole_Cells_dcast<-dcast(merge_mean_pvalue2[,c(1,5,6,3)], id_cp_interaction+Cell_pair~Treat_group)
head(whole_Cells_dcast);dim(whole_Cells_dcast)#149850      4
whole_Cells_dcast_rm1<-whole_Cells_dcast[which(!(is.na(whole_Cells_dcast$Abortion) & is.na(whole_Cells_dcast$CTRL))),]
head(whole_Cells_dcast_rm1);dim(whole_Cells_dcast_rm1)#136225      4

#decat cells
whole_Cells_dcast<-dcast(merge_mean_pvalue2[,c(1,5,6,7)], id_cp_interaction+Cell_pair~Treat_group)
head(whole_Cells_dcast);dim(whole_Cells_dcast)
whole_Cells_dcast_rm2<-whole_Cells_dcast[which(!(is.na(whole_Cells_dcast$Abortion) & is.na(whole_Cells_dcast$CTRL))),]
head(whole_Cells_dcast_rm2);dim(whole_Cells_dcast_rm2)#15238     4

whole_Cells_dcast_rm_pair<-intersect(whole_Cells_dcast_rm1$id_cp_interaction,whole_Cells_dcast_rm2$id_cp_interaction)
length(whole_Cells_dcast_rm_pair);dim(whole_Cells_dcast_rm1);dim(whole_Cells_dcast_rm2) #666  136225  15238 
head(merge_mean_pvalue2)
merge_mean_pvalue3<-merge_mean_pvalue2[which(as.character(merge_mean_pvalue2$id_cp_interaction) %in% whole_Cells_dcast_rm_pair),]
head(merge_mean_pvalue3);dim(merge_mean_pvalue3)#299700      7

#plot
merge_name<-merge(Abortion_mean_rm[,1:2], CTRL_mean_rm[,1:2], all = T, by = "id_cp_interaction")
str(merge_name);head(merge_name)
merge_name$interacting_pair.x<-as.character(merge_name$interacting_pair.x)
merge_name$interacting_pair.y<-as.character(merge_name$interacting_pair.y)

merge_name$interacting_pair<-ifelse(is.na(merge_name$interacting_pair.x),merge_name$interacting_pair.y, merge_name$interacting_pair.x)
head(merge_name);dim(merge_name);which(is.na(merge_name$interacting_pair))# 1199    4

merge_name2<-merge_name[which(merge_name$id_cp_interaction %in% merge_mean_pvalue3$id_cp_interaction),][,c(1,4)]
dim(merge_name2);head(merge_name2)#666   2

merge_mean_pvalue4<-merge(merge_name2, merge_mean_pvalue3, all = T, by = "id_cp_interaction")
merge_mean_pvalue4[1:6,]
Cell_type<-unique(c(unlist(lapply(strsplit(as.character(merge_mean_pvalue4$Cell_pair),"-"), function(x) x[1])),unlist(lapply(strsplit(as.character(merge_mean_pvalue4$Cell_pair),"-"), function(x) x[2]))))
Cell_type2<-as.data.frame(t(combn(Cell_type,2) ))

Cell_type_Dir_a<-c(paste(Cell_type2$V1,Cell_type2$V2,sep="-"))
Cell_type_Dir_b<-c(paste(Cell_type2$V2,Cell_type2$V1,sep="-"))
Cell_type_Dir_c<-c(paste(Cell_type,Cell_type,sep="-"))

Cell_type_Dir_a[order(Cell_type_Dir_a)]

#pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Three_types_all_cell_interaction_between_Abortion_and_CTRL.pdf")
#dev.off()
#dir_a
merge_mean_pvalue_dir_a<-merge_mean_pvalue4[which(merge_mean_pvalue4$Cell_pair %in% Cell_type_Dir_a),]
heatmap1<-ggplot(merge_mean_pvalue_dir_a[order(merge_mean_pvalue_dir_a$Cell_pair),],aes(Cell_pair,interacting_pair))+
  geom_point(aes(size=plogtran,color=mean),alpha =0.9,na.rm = TRUE)+
  facet_wrap( ~ Treat_group )+
  theme(axis.text.x = element_text(vjust=1,hjust=1,angle = 90,size=4),panel.border = element_rect(colour="black",fill=NA))+
  #scale_color_manual(values=ppCor[4:3])+
  labs(color="Mean_expression",
       size="-log10(p_value)",
       x="Cell Pairs",y="Ligand_Receptor",title="Cell Communiction::Cell_type_Dir_a")+
  scale_color_gradient(low="cyan",high ="red")+
  scale_size_continuous(range=c(0.1,3))

#dir_b
merge_mean_pvalue_dir_b<-merge_mean_pvalue4[which(merge_mean_pvalue4$Cell_pair %in% Cell_type_Dir_b),]
heatmap2<-ggplot(merge_mean_pvalue_dir_b[order(merge_mean_pvalue_dir_b$Cell_pair),],aes(Cell_pair,interacting_pair))+
  geom_point(aes(size=plogtran,color=mean),alpha =0.9,na.rm = TRUE)+
  facet_wrap( ~ Treat_group )+
  theme(axis.text.x = element_text(vjust=1,hjust=1,angle = 90,size=4),panel.border = element_rect(colour="black",fill=NA))+
  #scale_color_manual(values=ppCor[4:3])+
  labs(color="Mean_expression",
       size="-log10(p_value)",
       x="Cell Pairs",y="Ligand_Receptor",title="Cell Communiction::Cell_type_Dir_b")+
  scale_color_gradient(low="cyan",high ="red")+
  scale_size_continuous(range=c(0.1,3))

#dir_c
merge_mean_pvalue_dir_c<-merge_mean_pvalue4[which(merge_mean_pvalue4$Cell_pair %in% Cell_type_Dir_c),]
heatmap3<-ggplot(merge_mean_pvalue_dir_c[order(merge_mean_pvalue_dir_c$Cell_pair),],aes(Cell_pair,interacting_pair))+
  geom_point(aes(size=plogtran,color=mean),alpha =0.9,na.rm = TRUE)+
  facet_wrap( ~ Treat_group )+
  theme(axis.text.x = element_text(vjust=1,hjust=1,angle = 90,size=4),panel.border = element_rect(colour="black",fill=NA))+
  #scale_color_manual(values=ppCor[4:3])+
  labs(color="Mean_expression",
       size="-log10(p_value)",
       x="Cell Pairs",y="Ligand_Receptor",title="Cell Communiction :: self_interaction")+
  scale_color_gradient(low="cyan",high ="red")+
  scale_size_continuous(range=c(0.1,3))

write.table(as.data.frame(merge_mean_pvalue_dir_c),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/self_interaction_mean_pvalue_merge.csv',row.names = F,sep = ",")
write.table(as.data.frame(merge_mean_pvalue_dir_a),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/dir_a_interaction_mean_pvalue_merge.csv',row.names = F,sep = ",")
write.table(as.data.frame(merge_mean_pvalue_dir_b),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/dir_b_interaction_mean_pvalue_merge.csv',row.names = F,sep = ",")

ggsave(heatmap1,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Interaction_heatmap_dir_a.pdf",width = 24,height = 40)
ggsave(heatmap2,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Interaction_heatmap_dir_b.pdf",width = 24,height = 40)
ggsave(heatmap3,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Interaction_heatmap_dir_c.pdf",width = 10,height = 40)
