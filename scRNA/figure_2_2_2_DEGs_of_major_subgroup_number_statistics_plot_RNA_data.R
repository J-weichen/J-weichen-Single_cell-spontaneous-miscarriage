rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ggsci)
library(ggplot2)
library(scales)
library(dplyr)
library(data.table)
library(tidyverse)
library(ggpubr)
#set colors

##sample information plot in generation
#set colors
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)

pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:46))]
show_col(ppCor_all2)
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(80)
show_col(Cells_col)
cell_type_all0<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","STCs","FBs","Mycs","HCs","NKs","Ts","Bs")
conflict_prefer("which", "Matrix")
DEGs_merge_data <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge_new.txt", header=T)
head(DEGs_merge_data)
range(table(DEGs_merge_data$cluster))#16 962
range(table(DEGs_merge_data$cluster,DEGs_merge_data$trend))#0  660

DEGs_merge_data$Num<-1
DEGs_merge_data2<-aggregate(DEGs_merge_data$Num, by=list(DEGs_merge_data$cluster,DEGs_merge_data$trend), FUN=sum)
colnames(DEGs_merge_data2)<-c("cell_identity","tendency","DEGs_num")
DEGs_merge_data2$cell_identity<-as.character(DEGs_merge_data2$cell_identity)
DEGs_merge_data2$tendency<-factor(DEGs_merge_data2$tendency,levels=c("Up","Down"))

DEGs_merge_data02<-aggregate(DEGs_merge_data$Num, by=list(DEGs_merge_data$cluster), FUN=sum)
DEGs_merge_data02$Group.1<-as.character(DEGs_merge_data02$Group.1)
colnames(DEGs_merge_data02)<-c("cell_identity","DEGs_num2")
DEGs_merge_data02<-DEGs_merge_data02[order(DEGs_merge_data02$DEGs_num2,decreasing = T),]
cell_type_order<-DEGs_merge_data02$cell_identity

DEGs_merge_data3<-merge(DEGs_merge_data2,DEGs_merge_data02)
DEGs_merge_data3$cell_identity<-factor(DEGs_merge_data3$cell_identity,levels=cell_type_order)
str(DEGs_merge_data3)
head(DEGs_merge_data3)
DEGs_merge_data4 <-dcast(DEGs_merge_data3, cell_identity+DEGs_num2 ~ tendency , value.var = "DEGs_num")
head(DEGs_merge_data4)
DEGs_merge_data5<-merge(DEGs_merge_data3[,1:3],DEGs_merge_data4,by="cell_identity")
DEGs_merge_data5[is.na(DEGs_merge_data5)]<-0
range(DEGs_merge_data5$DEGs_num2)
#rose maps2
plot_mn<-ggplot(DEGs_merge_data5, aes(x = cell_identity, y = DEGs_num, fill = tendency)) +  
  geom_bar(stat = "identity", alpha = 0.6) + coord_polar() +
  theme_grey() + 
  labs(x = "", y = "", title = "The number of DEGs for each cell subtype") + 
  theme(legend.position = "right") +  
  theme(title = element_text(vjust = -56, face = "bold"))+
  scale_fill_manual(name="Active trend",values=ppCor_all[1:2])+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),#panel.grid.major.y  = element_blank(),
        axis.text.x = element_text(size = 6,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())


plot_mn2<-plot_mn+geom_text(aes(x= cell_identity, y= 1000, label =  DEGs_num2),data =DEGs_merge_data3,hjust =0.5,vjust =0.5,color="black", size=2 )
plot_mn2
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/roseplot_DEGs_number_decrease_in_all_celltypes.pdf",plot_mn2,width=10, height=10)


DEGs_merge_data5$cell_identity<-factor(DEGs_merge_data5$cell_identity,levels=cell_type_all0)
str(DEGs_merge_data5)
head(DEGs_merge_data5)

#rose maps1 
plot_mn3<-ggplot(DEGs_merge_data5, aes(x = cell_identity, y = DEGs_num, fill = tendency)) + 
  geom_bar(stat = "identity", alpha = 0.6) + coord_polar() +  scale_fill_manual(name="Active trend",values=ppCor_all[1:2])+
  geom_text(aes(x= cell_identity, y= 1000, label =  DEGs_num2),data =DEGs_merge_data3,hjust =0.5,vjust =0.5,color="black", size=2 )+
  theme_grey() + labs(x = "", y = "", title = "The number of DEGs for each cell subtype") + 
  theme(legend.position = "right") +  
  theme(title = element_text(vjust = -56, face = "bold"))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 6,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
plot_mn3
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/MF1i_roseplot_DEGs_number_order_in_all_celltypes.pdf",plot_mn3,width=10, height=10)

plot_mn4<-ggplot(DEGs_merge_data5, aes(x = cell_identity, y = DEGs_num, fill = tendency)) + 
  geom_bar(stat = "identity", alpha = 0.6) + coord_polar() +  scale_fill_manual(name="Active trend",values=ppCor_all[1:2])+
  geom_text(aes(x=cell_identity, y= 1000, label = Up),data =DEGs_merge_data5,hjust =0.5,vjust =0.5,color=ppCor[1], size=2 )+
  geom_text(aes(x=cell_identity, y= 900, label = Down),data =DEGs_merge_data5,hjust =0.5,vjust =0.5,color=ppCor[2], size=2 )+
  theme_grey() + labs(x = "", y = "", title = "The number of DEGs for each cell subtype") + 
  theme(legend.position = "right") +  
  theme(title = element_text(vjust = -56, face = "bold"))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 6,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
plot_mn4
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/MF1i_roseplot_DEGs_number_order_in_all_celltypes_anno_split.pdf",plot_mn4,width=10, height=10)

##=barplot ===========
num_cell_plot3 <-dcast(DEGs_merge_data2, cell_identity ~ tendency , value.var = "DEGs_num")
num_cell_plot3[is.na(num_cell_plot3)]<-0
num_cell_plot3$DEGs_total_num<-num_cell_plot3$Up+num_cell_plot3$Down
num_cell_plot3<-num_cell_plot3[order(num_cell_plot3$DEGs_total_num,decreasing = T),]
conflict_prefer("melt", "reshape2")
#barplot1
num_cell_plot3$cell_identity<-factor(num_cell_plot3$cell_identity,levels = as.character(num_cell_plot3$cell_identity))
num_cell_plot4<-melt(num_cell_plot3[,c(1:3)])
str(num_cell_plot4);head(num_cell_plot4)
bar_DEG_plot<-ggplot(data=num_cell_plot4, mapping=aes(x= cell_identity,y= value,fill=variable))+
  geom_bar(stat="identity",width=0.8,position= 'dodge', alpha = 0.6)+
  scale_fill_manual(name="Active trend",values=ppCor[1:2])+
  geom_text(aes(label=value,y=value+2),size=3,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="DEGs Number",x="Cell types",title="The number of Abortion related DEGs in each Cells")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 10),  axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")
bar_DEG_plot
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/SF4c_barplot_Up_Down_DEGs_number_decrease_in_all_celltypes.pdf",bar_DEG_plot,width=10, height=8)

#barplot2
num_cell_plot3$cell_identity<-factor(num_cell_plot3$cell_identity,levels = cell_type_all0)
num_cell_plot4<-melt(num_cell_plot3[,c(1:3)])
str(num_cell_plot4);head(num_cell_plot4)
bar_DEG_plot2<-ggplot(data=num_cell_plot4, mapping=aes(x= cell_identity,y= value,fill=variable))+
  geom_bar(stat="identity",width=0.8,position= 'dodge', alpha = 0.6)+
  scale_fill_manual(name="Active trend",values=ppCor[1:2])+
  geom_text(aes(label=value,y=value+2),size=3,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="DEGs Number",x="Cell types",title="The number of AMA related DEGs in each Cells")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 10),axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")
bar_DEG_plot2
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/SF4a_left_barplot_Up_Down_DEGs_number_order_in_all_celltypes.pdf",bar_DEG_plot2,width=10, height=8)

#top common genes
common_freq<-data.frame(table(DEGs_merge_data$gene,DEGs_merge_data$trend))
common_freq2 <-dcast(common_freq, Var1  ~ Var2  , value.var = "Freq")
common_freq2$add_freq<-common_freq2$Down + common_freq2$Up
head(common_freq2)
write.table(as.data.frame(common_freq2), file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Cell_type_number_for_Abortion_DEGs.xls",sep = "\t", row.names=T, col.names=T) 
#common_freq2<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Cell_type_number_for_Abortion_DEGs.xls",sep = "\t", header=T) 
head(common_freq2)
#distinct(DEGs_merge_data[,c("gene","DEGs_inter_type")])
common_freq2_down<-common_freq2[which(common_freq2$Down>0),]
common_freq2_Up<-common_freq2[which(common_freq2$Up>0),]
common_freq2_down<-common_freq2_down[order(common_freq2_down$Down,decreasing = T),]
common_freq2_Up<-common_freq2_Up[order(common_freq2_Up$Up,decreasing = T),]
head(common_freq2_down);head(common_freq2_Up)
data_down<-head(common_freq2_down,n=10)
data_down$Var1<-factor(data_down$Var1,levels = as.character(data_down$Var1))
data_Up<-head(common_freq2_Up,n=10)
data_Up$Var1<-factor(data_Up$Var1,levels = as.character(data_Up$Var1))
plot_Up<-ggplot(data_Up, aes(x = Var1, y =  Up)) +  geom_bar(stat = "identity", alpha = 0.9,fill=ppCor_all[1]) + coord_polar() +
  geom_text(aes(x= Var1, y= Up+2, label =  Up),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of cell subtype for top 10 common Up_DEGs") + 
  theme(title = element_text(vjust = -56, face = "bold"))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
plot_down<-ggplot(data_down, aes(x = Var1, y =  Down)) +  geom_bar(stat = "identity", alpha = 0.9,fill=ppCor_all[2]) + coord_polar() +
  geom_text(aes(x= Var1, y= Down+1, label =  Down),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of cell subtype for top 10 common Down_DEGs") + 
  theme(title = element_text(vjust = -56, face = "bold"))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
plot_Up
plot_down
top10_rose<-grid.arrange(plot_Up, plot_down,ncol=2)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Top10_common_number_for_Abortion_related_DEGs.pdf",top10_rose,width=12, height=6)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Top10_common_number_for_Abortion_related_DEGs.pdf",width =6,height = 6)
grid.arrange(plot_Up, plot_down,ncol=2)
dev.off()
distinct(DEGs_merge_data[which(DEGs_merge_data$gene %in% as.character(data_down$Var1)),-c(2:5)])

#plot for selected DEGs
gene_target<-unique(c(as.character(data_down$Var1),as.character(data_Up$Var1)))
data_target<-common_freq[which(common_freq$Var1 %in% gene_target),]
data_target$Var1<-factor(data_target$Var1,levels=rev(gene_target))
data_target$Var2<-factor(data_target$Var2,levels=c("Up","Down"))

data_target

Top10_DEGs_plot<-ggplot(data_target, aes(x = Var1, y = Freq, fill = Var2)) +  
  geom_bar(stat = "identity", alpha = 0.7) + coord_polar() +
  geom_text(aes(x= Var1, y= max(Freq+1), label =  Freq),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of Cells type for DEGs(Top10)") + 
  theme(legend.position = "right") +  
  theme(title = element_text(vjust = -56, face = "bold"))+
  scale_fill_manual(name="Active trend",values=ppCor_all[1:2])+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
Top10_DEGs_plot
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/S4d_rose_plot_top10_number_cell_type_Up_Down_common_DEGs.pdf",Top10_DEGs_plot,width=10, height=10)


Distribution1<-ggbarplot(data.frame(table(common_freq2$add_freq)), x="Var1", y="Freq",fill="navy",color = "white")+#,rotate=TRUE
  geom_text(aes(label=Freq),size=2,position = position_dodge(0.5),vjust=-0.5)+
  labs(y="DEGs  Number",x="Cell type Number",title="Distribution of Abortion DEGs among Cells ")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),axis.ticks = element_blank(),
        legend.position ="top",legend.direction = "horizontal")
DEGs_freq0<-common_freq[which(common_freq$Freq>0),]
DEGs_freq<-data.frame(table(DEGs_freq0$Freq,DEGs_freq0$Var2))
head(DEGs_freq)
colnames(DEGs_freq)<-c("cell_types_number","trend","TFs_type_freq")
DEGs_freq$trend<-factor(DEGs_freq$trend,levels=c("Up","Down"))

Distribution2<-ggplot(data=DEGs_freq, mapping=aes(x= cell_types_number,y=TFs_type_freq,fill=trend))+
  geom_bar(stat="identity",width=0.8,position= 'dodge')+
  scale_fill_manual(name="Change trend",values=ppCor[1:2])+
  geom_text(aes(y=TFs_type_freq+6,label=TFs_type_freq),size=2,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="DEGs  Number",x="Cell type Number",title="Distribution of Abortion  DEGs among Cells ")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")# + coord_flip()

Distribution1
Distribution2
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Distribution_of_common_number_for_Abortion_related_DEGs.pdf",width =12,height = 6)
grid.arrange(Distribution1,Distribution2,ncol=2)
dev.off()

######以下未跑：特定基因的细胞频率绘制
#plot for all DEGs
head(DEGs_merge_data)
common_freq<-data.frame(table(DEGs_merge_data$gene,DEGs_merge_data$trend))
length(unique(as.character(common_freq$Var1)))#702
common_freq2 <-dcast(common_freq, Var1  ~ Var2  , value.var = "Freq")
common_freq2$add_freq<-common_freq2$Down + common_freq2$Up
common_freq2<-common_freq2[order(common_freq2$add_freq,decreasing = T),]
head(common_freq2);nrow(common_freq2)#363
common_freq$Var1<-factor(common_freq$Var1,levels = common_freq2$Var1)
common_freq$Var2<-factor(common_freq$Var2,levels =c("Up","Down"))

ggplot(data=common_freq, mapping=aes(x= Var1,y= Freq,fill=Var2))+
  geom_bar(stat="identity",width=0.8,position= 'dodge')+
  scale_fill_manual(name="Active trend",values=ppCor[1:2])+
  geom_text(aes(label=Freq,y=Freq+0.5),size=3,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="Cell type Number",x="DEGs",title="The number of Cells type for each DEGs")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 10),  axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")

#plot for selected DEGs
gene_target<-c("TRAF4","ASCL2","YBX1","JUN","TFAP2A","REL","HIF1A","IRF1","ETS2","CEBPD","CREM","FOSB","JUNB","CEBPB","KLF5","KLF4","EGR1","ATF3","FOS")
data_target<-common_freq[which(common_freq$Var1 %in% gene_target),]
data_target$Var1<-factor(data_target$Var1,levels=rev(gene_target))

ggplot(data=data_target, mapping=aes(x= Var1,y= Freq,fill=Var2))+
  geom_bar(stat="identity",width=0.8,position= 'dodge')+
  scale_fill_manual(name="Active trend",values=ppCor[1:2])+
  geom_text(aes(label=Freq,y=Freq+0.2),size=3,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="Cell type Number",x="DEGs",title="The number of Cells type for DEGs overlaped with daTFs")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),  axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")

ggplot(data_target, aes(x = Var1, y = Freq, fill = Var2)) +  
  geom_bar(stat = "identity", alpha = 0.7) + coord_polar() +
  geom_text(aes(x= Var1, y= 10, label =  Freq),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of Cells type for DEGs overlaped with daTFs") + 
  theme(legend.position = "right") +  
  theme(title = element_text(vjust = -56, face = "bold"))+
  scale_fill_manual(name="Active trend",values=ppCor_all[1:2])+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())

