rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)
library(ggsci)
library(grid)
library(futile.logger)

pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
show_col(pal)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])

#2. dotplot for inter_pair
Used_cell<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Mycs","HCs","STCs","FBs","Endo","NKs","Ts")

Abortion_mean<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/Abortion/means.txt",head = T,sep = "\t")
CTRL_mean<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/CTRL/means.txt",head = T,sep = "\t")
Abortion_pvalue<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/Abortion/pvalues.txt",head = T,sep = "\t")
CTRL_pvalue<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/CTRL/pvalues.txt",head = T,sep = "\t")
dim(CTRL_pvalue);dim(Abortion_pvalue)##1199  372 ##1199  372
dim(CTRL_mean);dim(Abortion_mean)##1199  372 ##1199  372
Abortion_mean[1:6,10:16];Abortion_pvalue[1:6,10:16]
##mean
Cell_pair_CTRL1<-unlist(lapply(strsplit(as.character(colnames(CTRL_mean)),"[.]"), function(x) x[1]))
Cell_pair_CTRL2<-unlist(lapply(strsplit(as.character(colnames(CTRL_mean)),"[.]"), function(x) x[2]))
CTRL_mean_rm<-CTRL_mean[,c(1:11,which(Cell_pair_CTRL1 %in% Used_cell & Cell_pair_CTRL2 %in% Used_cell))]
CTRL_mean[10:16,10:16];CTRL_mean_rm[10:16,10:16]
Cell_pair_Abortion1<-unlist(lapply(strsplit(as.character(colnames(Abortion_mean)),"[.]"), function(x) x[1]))
Cell_pair_Abortion2<-unlist(lapply(strsplit(as.character(colnames(Abortion_mean)),"[.]"), function(x) x[2]))
Abortion_mean_rm<-Abortion_mean[,c(1:11,which(Cell_pair_Abortion1 %in% Used_cell & Cell_pair_Abortion2 %in% Used_cell))]
Abortion_mean[1:6,1:6];Abortion_mean_rm[1:6,1:6]
dim(CTRL_mean_rm);dim(Abortion_mean_rm)#1199    236 #1199    236
##p value 
Cell_pair_CTRL1<-unlist(lapply(strsplit(as.character(colnames(CTRL_pvalue)),"[.]"), function(x) x[1]))
Cell_pair_CTRL2<-unlist(lapply(strsplit(as.character(colnames(CTRL_pvalue)),"[.]"), function(x) x[2]))
CTRL_pvalue_rm<-CTRL_pvalue[,c(1:11,which(Cell_pair_CTRL1 %in% Used_cell & Cell_pair_CTRL2 %in% Used_cell))]
CTRL_pvalue[1:6,1:6];CTRL_pvalue_rm[1:6,1:6]
Cell_pair_Abortion1<-unlist(lapply(strsplit(as.character(colnames(Abortion_pvalue)),"[.]"), function(x) x[1]))
Cell_pair_Abortion2<-unlist(lapply(strsplit(as.character(colnames(Abortion_pvalue)),"[.]"), function(x) x[2]))
Abortion_pvalue_rm<-Abortion_pvalue[,c(1:11,which(Cell_pair_Abortion1 %in% Used_cell & Cell_pair_Abortion2 %in% Used_cell))]
Abortion_pvalue[1:6,1:6];Abortion_pvalue_rm[1:6,1:6]
dim(CTRL_pvalue_rm);dim(Abortion_pvalue_rm)#1199  236

##selected for target intermolecular pairs
merge_mean_data<-merge(Abortion_mean_rm[,c(1,12:ncol(Abortion_mean_rm))], CTRL_mean_rm[,c(1,12:ncol(CTRL_mean_rm))], all = T, by = "id_cp_interaction")
merge_pvalue_data<-merge(Abortion_pvalue_rm[,c(1,12:ncol(Abortion_pvalue_rm))], CTRL_pvalue_rm[,c(1,12:ncol(CTRL_pvalue_rm))], all = T, by = "id_cp_interaction")
merge_pvalue_data[1:6,1:6];merge_mean_data[1:6,1:6]

merge_mean_melt<-melt(merge_mean_data);merge_pvalue_melt<-melt(merge_pvalue_data)
colnames(merge_mean_melt)<-c("id_cp_interaction","variable","mean")
colnames(merge_pvalue_melt)<-c("id_cp_interaction","variable","pvalue")
merge_mean_melt[1:6,]

merge_mean_pvalue<-merge(merge_mean_melt, merge_pvalue_melt, by = c("id_cp_interaction","variable"),all = T)
head(merge_mean_pvalue);dim(merge_mean_pvalue)

merge_mean_pvalue$Treat_group<-unlist(lapply(strsplit(as.character(merge_mean_pvalue$variable),"[.]"), function(x) x[3]))
merge_mean_pvalue$Treat_group<-as.factor(ifelse(grepl("x",merge_mean_pvalue$Treat_group), "Abortion", "CTRL"))
merge_mean_pvalue$Cell_pair<-unlist(lapply(strsplit(as.character(merge_mean_pvalue$variable),"[.]"), function(x) paste(x[1],x[2],sep="-")))
merge_mean_pvalue$plogtran<-c(-log((merge_mean_pvalue$pvalue+0.00001),10))
head(merge_mean_pvalue)
merge_mean_pvalue2 <- transform(merge_mean_pvalue, plogtran = ifelse(pvalue >= 0.05, NA, plogtran), mean = ifelse(mean ==0, NA, mean))
head(merge_mean_pvalue2)

## for all cell with exist inter_pairs
whole_Cells_dcast<-dcast(merge_mean_pvalue2[,c(1,5,6,3)], id_cp_interaction+Cell_pair~Treat_group)
head(whole_Cells_dcast);dim(whole_Cells_dcast)# 269775      4
whole_Cells_dcast_rm1<-whole_Cells_dcast[which(!(is.na(whole_Cells_dcast$Abortion) & is.na(whole_Cells_dcast$CTRL))),]
head(whole_Cells_dcast_rm1);dim(whole_Cells_dcast_rm1)#194645      4
whole_Cells_dcast<-dcast(merge_mean_pvalue2[,c(1,5,6,7)], id_cp_interaction+Cell_pair~Treat_group)
whole_Cells_dcast_rm2<-whole_Cells_dcast[which(!(is.na(whole_Cells_dcast$Abortion) & is.na(whole_Cells_dcast$CTRL))),]
head(whole_Cells_dcast_rm2);dim(whole_Cells_dcast_rm2)#15238     4

##get final mean pvalue matrix
whole_Cells_dcast_rm_pair<-intersect(whole_Cells_dcast_rm1$id_cp_interaction,whole_Cells_dcast_rm2$id_cp_interaction)
length(whole_Cells_dcast_rm_pair);dim(whole_Cells_dcast_rm1);dim(whole_Cells_dcast_rm2) #666  194645    15238 
head(merge_mean_pvalue2)

merge_mean_pvalue3<-merge_mean_pvalue2[which(as.character(merge_mean_pvalue2$id_cp_interaction) %in% whole_Cells_dcast_rm_pair),]
head(merge_mean_pvalue3);dim(merge_mean_pvalue3)#299700      7

merge_name<-merge(Abortion_mean_rm[,1:2], CTRL_mean_rm[,1:2], all = T, by = "id_cp_interaction")
merge_name$interacting_pair.x<-as.character(merge_name$interacting_pair.x);merge_name$interacting_pair.y<-as.character(merge_name$interacting_pair.y)
merge_name$interacting_pair<-ifelse(is.na(merge_name$interacting_pair.x),merge_name$interacting_pair.y, merge_name$interacting_pair.x)
merge_name2<-merge_name[which(merge_name$id_cp_interaction %in% merge_mean_pvalue3$id_cp_interaction),][,c(1,4)]
dim(merge_name2);head(merge_name2)#666   2
merge_mean_pvalue4<-merge(merge_name2, merge_mean_pvalue3, all = T, by = "id_cp_interaction")
merge_mean_pvalue4[1:6,]
write.table(as.data.frame(merge_mean_pvalue4), file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/cellpheno_merge_mean_pvalue.xls",sep = "\t", row.names=T, col.names=T) 

Cell_type<-unique(c(unlist(lapply(strsplit(as.character(merge_mean_pvalue4$Cell_pair),"-"), function(x) x[1])),unlist(lapply(strsplit(as.character(merge_mean_pvalue4$Cell_pair),"-"), function(x) x[2]))))
Cell_type2<-as.data.frame(t(combn(Cell_type,2) ))
Cell_type_Dir_a<-c(paste(Cell_type2$V1,Cell_type2$V2,sep="-"))
Cell_type_Dir_b<-c(paste(Cell_type2$V2,Cell_type2$V1,sep="-"))
Cell_type_Dir_c<-c(paste(Cell_type,Cell_type,sep="-"))
Cell_type_Dir_a[order(Cell_type_Dir_a)]

####
merge_mean_pvalue_CTRL<-na.omit(merge_mean_pvalue4[which(merge_mean_pvalue4$Treat_group =="CTRL" ),])
merge_mean_pvalue_Abortion<-na.omit(merge_mean_pvalue4[which(merge_mean_pvalue4$Treat_group =="Abortion" ),])
CPI_CTRL<-unique(merge_mean_pvalue_CTRL$id_cp_interaction);CPI_Abortion<-unique(merge_mean_pvalue_Abortion$id_cp_interaction)
CPI_CTRL_sp<-setdiff(CPI_CTRL,CPI_Abortion)#33
CPI_Abortion_sp<-setdiff(CPI_Abortion,CPI_CTRL)##51


#pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Three_types_all_cell_interaction_between_Abortion_and_CTRL.pdf")
#dev.off()
merge_mean_pvalue_plot<-merge_mean_pvalue4[which(merge_mean_pvalue4$id_cp_interaction  %in% c(CPI_Abortion_sp,CPI_CTRL_sp)),]
merge_mean_pvalue_plot2<-na.omit(merge_mean_pvalue_plot)
merge_mean_pvalue_plot2$id_cp_interaction<- factor(merge_mean_pvalue_plot2$id_cp_interaction,levels = c(CPI_Abortion_sp,CPI_CTRL_sp))
merge_mean_pvalue_plot2<-merge_mean_pvalue_plot2[order(merge_mean_pvalue_plot2$id_cp_interaction,decreasing = T),]
merge_mean_pvalue_plot2$interacting_pair<- factor(merge_mean_pvalue_plot2$interacting_pair,levels = as.character(unique(merge_mean_pvalue_plot2$interacting_pair)))
merge_mean_pvalue_plot2<-merge_mean_pvalue_plot2[order(merge_mean_pvalue_plot2$interacting_pair,decreasing = T),]

head(merge_mean_pvalue_plot2)
write.table(merge_mean_pvalue_plot2, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/group_unique_CPI.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

heatmap1<-ggplot(merge_mean_pvalue_plot2[order(merge_mean_pvalue_plot2$Cell_pair),],aes(Cell_pair,interacting_pair))+
  geom_point(aes(size=plogtran,color=mean),alpha =0.9,na.rm = TRUE)+
  theme(axis.text.x = element_text(vjust=1,hjust=1,angle = 90,size=4),panel.border = element_rect(colour="black",fill=NA))+
  labs(color="Mean_expression",size="-log10(p_value)",
       x="Cell Pairs",y="Ligand_Receptor",title="Cell Communiction::Cell_type_all")+
  scale_color_gradient(low="cyan",high ="red")+ scale_size_continuous(range=c(0.1,3)) + facet_wrap( ~ Treat_group )
heatmap1

plot_dira<-merge_mean_pvalue_plot2[which(merge_mean_pvalue_plot2$Cell_pair %in% Cell_type_Dir_a),]
plot_dira$Cell_pair<-factor(plot_dira$Cell_pair,levels = unique(plot_dira$Cell_pair))
heatmap_dira<-ggplot(plot_dira,aes(Cell_pair,interacting_pair))+
  geom_point(aes(size=plogtran,color=mean),alpha =0.9,na.rm = TRUE)+
  theme(axis.text.x = element_text(vjust=0.5,hjust=1,angle = 90,size=8),panel.border = element_rect(colour="black",fill=NA))+
  labs(color="Mean_expression",size="-log10(p_value)",
       x="Cell Pairs",y="Ligand_Receptor",title="Cell Communiction::Cell_type_Dir_a")+
  scale_color_gradient(low="cyan",high ="red")+ scale_size_continuous(range=c(0.1,3)) + facet_wrap( ~ Treat_group )
heatmap_dira

plot_dirb<-merge_mean_pvalue_plot2[which(merge_mean_pvalue_plot2$Cell_pair %in% Cell_type_Dir_b),]
plot_dirb$Cell_pair<-factor(plot_dirb$Cell_pair,levels = unique(plot_dirb$Cell_pair))
heatmap_dirb<-ggplot(plot_dirb,aes(Cell_pair,interacting_pair))+
  geom_point(aes(size=plogtran,color=mean),alpha =0.9,na.rm = TRUE)+
  theme(axis.text.x = element_text(vjust=0.5,hjust=1,angle = 90,size=8),panel.border = element_rect(colour="black",fill=NA))+
  labs(color="Mean_expression",size="-log10(p_value)",
       x="Cell Pairs",y="Ligand_Receptor",title="Cell Communiction::Cell_type_Dir_b")+
  scale_color_gradient(low="cyan",high ="red")+ scale_size_continuous(range=c(0.1,3)) + facet_wrap( ~ Treat_group )
heatmap_dirb

plot_dirc<-merge_mean_pvalue_plot2[which(merge_mean_pvalue_plot2$Cell_pair %in% Cell_type_Dir_c),]
plot_dirc$Cell_pair<-factor(plot_dirc$Cell_pair,levels = unique(plot_dirc$Cell_pair))
heatmap_dirc<-ggplot(plot_dirc,aes(Cell_pair,interacting_pair))+
  geom_point(aes(size=plogtran,color=mean),alpha =0.9,na.rm = TRUE)+
  theme(axis.text.x = element_text(vjust=0.5,hjust=1,angle = 90,size=8),panel.border = element_rect(colour="black",fill=NA))+
  labs(color="Mean_expression",size="-log10(p_value)",
       x="Cell Pairs",y="Ligand_Receptor",title="Cell Communiction::Cell_type_Dir_c")+
  scale_color_gradient(low="cyan",high ="red")+ scale_size_continuous(range=c(0.1,3)) + facet_wrap( ~ Treat_group )
heatmap_dirc

write.table(as.data.frame(merge_mean_pvalue_plot2),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_CTRL_special_whole_cell_pair_interaction_mean_pvalue_merge.csv',row.names = F,sep = ",")
write.table(as.data.frame(plot_dirc),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_CTRL_special_self_cell_pair_interaction_mean_pvalue_merge.csv',row.names = F,sep = ",")
write.table(as.data.frame(plot_dira),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_CTRL_special_dir_a_cell_pair_interaction_mean_pvalue_merge.csv',row.names = F,sep = ",")
write.table(as.data.frame(plot_dirb),'/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_CTRL_special_dir_b_cell_pair_interaction_mean_pvalue_merge.csv',row.names = F,sep = ",")

ggsave(heatmap1,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_CTRL_special_whole_cell_pair_Interaction_heatmap.pdf",width = 30,height = 12)
ggsave(heatmap_dira,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_CTRL_special_dir_a_cell_pair_Interaction_heatmap.pdf",width = 20,height = 12)
ggsave(heatmap_dirb,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_CTRL_special_dir_b_cell_pair_Interaction_heatmap.pdf",width = 20,height = 12)
ggsave(heatmap_dirc,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_CTRL_special_dir_c_cell_pair_Interaction_heatmap.pdf",width = 10,height = 12)
