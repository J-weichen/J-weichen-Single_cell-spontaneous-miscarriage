rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
options(stringsAsFactors = FALSE)
## Load packages
library(reshape2)
library(dendextend)
library(ggsci)
library(ggpubr)
library(gridExtra)
library(scales)
library(VennDiagram)

##sample information plot in generation
#set colors
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

#for plot evaluation

#for plot evaluation
D_ac_TFs_active_high_matrix <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Different_TFs/merge_ALLcell_all_genes_D_ac_TFs_active_high_matrix.txt", header=T)
range(table(D_ac_TFs_active_high_matrix$label))#1 13
head(D_ac_TFs_active_high_matrix)

conflict_prefer("rbind", "spam")

cell_order<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","STCs","FBs","Mycs","HCs","NKs","Ts","Bs")
up_matrix<-D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$change =="Up"),]
down_matrix<-D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$change =="Down"),]
up_matrix$cell_identity<-factor(up_matrix$cell_identity,levels=cell_order)
up_matrix<-up_matrix[order(up_matrix$cell_identity),]
down_matrix$cell_identity<-factor(down_matrix$cell_identity,levels=cell_order)
down_matrix<-down_matrix[order(down_matrix$cell_identity),]

merge_data<-as.data.frame(rbind(up_matrix,down_matrix))
gene_order<-unique(merge_data$label)

##formal plot for major figure three 
conflict_prefer("dcast", "reshape2")
merge_data1<-merge_data[,c("cell_identity","label","foldchange")]
merge_data1[which(merge_data1$foldchange ==1),]
merge_data2<- dcast(merge_data1,label~cell_identity)
merge_data2$label<- factor(merge_data2$label,levels=gene_order,ordered=TRUE)
merge_data2<-merge_data2[order(merge_data2$label),]

head(merge_data2);dim(merge_data2)# 364  18
rownames(merge_data2)<-merge_data2$label
merge_data2<-merge_data2[,-(1)]
merge_data2[is.na(merge_data2)] <- 0
merge_data2[merge_data2>1] <- 1
merge_data2[merge_data2<1 & merge_data2>0] <- c(-1)
head(merge_data2);tail(merge_data2)

merge_data3<-t(merge_data2)
dim(merge_data3) ##17 364
range(merge_data3)# -1  1
merge_data3[1:6,1:6]

#Building anno 
cell_type_all2<-c(rep("Cytotrophoblasts_PAGE4",2),rep("Syncytiotrophoblast_CGA",3),rep("Extravilous_trophoblast_HLA_G",3),
                  rep("Epithelial_Cell_EPCAM",1),rep("Endothelial_Cell_PECAM1",1),rep("Stromal_cells_DCN",2),
                  rep("Myeloid_Cell_AIF1",2),rep("Leukomonocyte",3))
cell_type_all3<-c(rep("Trophoblasts_KRT7",8),rep("Epithelial_Cell_EPCAM",1),rep("Endothelial_Cell_PECAM1",1),rep("Stromal_cells_DCN",2),
                  rep("Myeloid_Cell_AIF1",2),rep("Leukomonocyte",3))

annotation_cell_type <-data.frame(cluster=cell_order,main_group=cell_type_all2,main_group2=cell_type_all3)
rownames(annotation_cell_type) = annotation_cell_type$cluster
head(annotation_cell_type)

anno_colors = list(
  main_group =c(Cytotrophoblasts_PAGE4=Cells_col[1],Extravilous_trophoblast_HLA_G=Cells_col[18],
                Syncytiotrophoblast_CGA=Cells_col[9],Epithelial_Cell_EPCAM=Cells_col[30],Endothelial_Cell_PECAM1=Cells_col[27],Stromal_cells_DCN=Cells_col[14],
                Myeloid_Cell_AIF1=Cells_col[22], Leukomonocyte=Cells_col[54]),
  main_group2 =c(Trophoblasts_KRT7=Cells_col[1], Epithelial_Cell_EPCAM=Cells_col[30],Endothelial_Cell_PECAM1=Cells_col[27],Stromal_cells_DCN=Cells_col[14],
                 Myeloid_Cell_AIF1=Cells_col[22], Leukomonocyte=Cells_col[54]))
labels_col = c("")

p1<-pheatmap(merge_data3,cluster_rows=F, cluster_cols =F, 
             annotation_row=annotation_cell_type[,c("main_group","main_group2")],
             labels_col = labels_col,
             annotation_colors = anno_colors,
             #gaps_row = c(18),
             main = "Disease all DEG_defualt(q=0.1)",
             legend_breaks = c(-0.6,0,0.6),legend_labels = c("Down","Unchange","Up"),
             color = colorRampPalette(colors = c(ppCor_all2[22],"lightgray",ppCor_all2[21]))(3))

p2<-pheatmap(merge_data3,cluster_rows=F, cluster_cols =T, 
             annotation_row=annotation_cell_type[,c("main_group","main_group2")],labels_col = labels_col,
             annotation_colors = anno_colors,
             main = "Disease all DEG_defualt(q=0.1)",
             legend_breaks = c(-0.6,0,0.6),legend_labels = c("Down","Unchange","Up"),
             color = colorRampPalette(colors = c(ppCor_all2[22],"lightgray",ppCor_all2[21]))(3))

p3<-pheatmap(merge_data3,cluster_rows=T, cluster_cols =T, 
             annotation_row=annotation_cell_type[,c("main_group","main_group2")],labels_col = labels_col,
             annotation_colors = anno_colors,
             main = "Disease all DEG_defualt(q=0.1)",
             legend_breaks = c(-0.6,0,0.6),legend_labels = c("Down","Unchange","Up"),
             color = colorRampPalette(colors = c(ppCor_all2[22],"lightgray",ppCor_all2[21]))(3))

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/SF4b_heat_plot_for_Abortion_related_DacTFs_all_subcluster1.pdf",width =10,height = 6)
print(p1)
dev.off()
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/SF4b_heat_plot_for_Abortion_related_DacTFs_all_subcluster2.pdf",width =10,height = 6)
print(p2)
dev.off()
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/SF4b_heat_plot_for_Abortion_related_DacTFs_all_subcluster3.pdf",width =10,height = 6)
print(p3)
dev.off()
#存储相关关系文件
head(merge_data2)
write.table(merge_data2, file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/relationship_for_Cellsubtypes_DacTFs_active_high_matrix.txt",quote=F, row.names=F, col.names=T,sep="\t") 

#selected candiidated group for venn plot
merge_data0 <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Different_TFs/merge_ALLcell_all_genes_D_ac_TFs_active_high_matrix.txt", header=T)
head(merge_data0)
merge_data_rm0<-merge_data0[which(!(merge_data0$cell_identity %in% c("Ery"))),]

Vagina_up<-unique(as.character(merge_data_rm0[which(merge_data_rm0$change =="Up"),]$lable))
length(Vagina_up)# 1528
Vagina_down<-unique(as.character(merge_data_rm0[which(merge_data_rm0$change =="Down"),]$lable))
length(Vagina_down)# 1627
Up_gene<-unique(c(Vagina_up));Down_gene<-unique(c(Vagina_down))
#FOR UP_gene VERSUS Down_gene
venn <-venn.diagram(list(Up_TF=Up_gene,Down_TF= Down_gene),
                    alpha=c(0.8,0.8),lwd=1,lty=1, col="white",fill=c(ppCor_all2[21],ppCor_all2[22]), 
                    cex = 1.5,cat.col=c(ppCor_all2[21],ppCor_all2[22]), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=90,filename = NULL)
grid.newpage(); grid.draw(venn)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/SF4b_venn_plot_for_Abortion_related_DacTFs_all_subcluster.pdf",width = 6,height = 6)
grid.draw(venn)
dev.off()


cell_name="ALLcell_all_genes"
#target_final_group="Treat"
D_ac_TF_matrix <- read.table(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Different_TFs/merge_",cell_name,"_D_ac_TF_matrix.txt"), header=T)
D_ac_TFs_active_high_matrix <- read.table(paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Different_TFs/merge_",cell_name,"_D_ac_TFs_active_high_matrix.txt"), header=T)

range(table(D_ac_TF_matrix$label))#1 17
range(table(D_ac_TFs_active_high_matrix$label))#1 13
head(D_ac_TF_matrix)
head(D_ac_TFs_active_high_matrix)
##formal plot for major figure three 

## De-TFs number for each Cells
num_cell_plot<-D_ac_TFs_active_high_matrix[,c("cell_identity","change")]
num_cell_plot$Num<-1
num_cell_plot2<-aggregate(num_cell_plot$Num, by=list(num_cell_plot$cell_identity,num_cell_plot$change), FUN=sum)
#table(D_ac_TFs_active_high_matrix$cell_identity,D_ac_TFs_active_high_matrix$change)
colnames(num_cell_plot2)<-c("cell_identity","tendency","TFs_num")
cell_order<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","STCs","FBs","Mycs","HCs","NKs","Ts","Bs")

num_cell_plot2$cell_identity<-factor(num_cell_plot2$cell_identity,levels=cell_order)
num_cell_plot2$tendency<-factor(num_cell_plot2$tendency,levels=c("Up","Down"))
str(num_cell_plot2)
head(num_cell_plot2)
#rose maps 
TF_roseplot<-ggplot(num_cell_plot2, aes(x = cell_identity, y = TFs_num, fill = tendency)) +  
  geom_bar(stat = "identity", alpha = 0.5) + coord_polar() +
  theme_grey() + labs(x = "", y = "", title = "The number of TFs for each cell subtype") + 
  geom_text(aes(x= cell_identity, y= 100, label =  TFs_num),data =num_cell_plot2,hjust =0.5,vjust =0.5,color="black", size=2 )+
  # theme(axis.text.y = element_blank()) +    theme(axis.ticks = element_blank()) +      ## 去掉左上角的刻度标签## 去掉左上角的刻度线
  # theme(panel.border = element_blank()) +   ## 去掉外层边框
  theme(legend.position = "right") +   ## 去掉图例
  theme(title = element_text(vjust = -56, face = "bold"))+
  scale_fill_manual(name="Active trend",values=c(ppCor_all2[21],ppCor_all2[22]))
TF_roseplot
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/MF1j_TF_roseplot_Up_Down_dacTFs_number_order_in_all_celltypes.pdf",TF_roseplot,width=10, height=8)

#conflict_prefer("melt", "reshape2")
num_cell_plot3 <-dcast(num_cell_plot2, cell_identity ~ tendency , value.var = "TFs_num")
num_cell_plot3[is.na(num_cell_plot3)]<-0
num_cell_plot3$TFs_total_num<-num_cell_plot3$Up+num_cell_plot3$Down
num_cell_plot3<-num_cell_plot3[order(num_cell_plot3$TFs_total_num,decreasing = T),]
num_cell_plot3$cell_identity<-factor(num_cell_plot3$cell_identity,levels = as.character(num_cell_plot3$cell_identity))
head(num_cell_plot2)
TFs_merge_data<-merge(num_cell_plot2,num_cell_plot3,by="cell_identity")
TFs_merge_data[is.na(TFs_merge_data)]<-0
range(TFs_merge_data$TFs_total_num)#8 109

plot_mn4<-ggplot(TFs_merge_data, aes(x = cell_identity, y = TFs_num, fill = tendency)) + 
  geom_bar(stat = "identity", alpha = 0.6) + coord_polar() +  scale_fill_manual(name="Active trend",values=ppCor_all2[21:22])+
  geom_text(aes(x=cell_identity, y= 115, label = Up),data =TFs_merge_data,hjust =0.5,vjust =0.5,color=ppCor_all2[21], size=2 )+
  geom_text(aes(x=cell_identity, y= 105, label = Down),data =TFs_merge_data,hjust =0.5,vjust =0.5,color=ppCor_all2[22], size=2 )+
  theme_grey() + labs(x = "", y = "", title = "The number of DacTFs for each cell subtype") + 
  theme(legend.position = "right") +  
  theme(title = element_text(vjust = -56, face = "bold"))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 6,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
plot_mn4
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/MF1j_TF_roseplot_Up_Down_dacTFs_number_order_in_all_celltypes_number_split.pdf",plot_mn4,width=10, height=8)




num_cell_plot4<-melt(num_cell_plot3[,c(1:3)])
num_cell_plot4$cell_identity<-factor(num_cell_plot4$cell_identity,levels=cell_order)
num_cell_plot4$variable<-factor(num_cell_plot4$variable,levels=c("Up","Down"))
str(num_cell_plot4);head(num_cell_plot4)
TF_number_barplot<-ggplot(data=num_cell_plot4, mapping=aes(x= cell_identity,y=value,fill=variable))+
  geom_bar(stat="identity",width=0.8,position= 'dodge',alpha=0.5)+
  scale_fill_manual(name="Active trend",values=c(ppCor_all2[21],ppCor_all2[22]))+
  geom_text(aes(label=value),size=3,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="def_acTFs Number",x="Cell types",title="The number of differential AMA active_high TFs in each Cells ")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0,hjust=1,angle = 90),
        legend.title = element_text(size = 10),  axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")
TF_number_barplot
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/SF4e_barplot_Up_Down_dac_TF_number_order_in_all_celltypes.pdf",TF_number_barplot,width=10, height=8)

#进一步统计组间重复
common_freq<-data.frame(table(D_ac_TF_matrix$label))
common_freq<-common_freq[order(common_freq$Freq,decreasing = T),]

common_freq_active_high<-data.frame(table(D_ac_TFs_active_high_matrix$label))
common_freq_active_high<-common_freq_active_high[order(common_freq_active_high$Freq,decreasing = T),]

length(unique(common_freq$Var1));length(unique(common_freq_active_high$Var1))##539 364   
#for common TFs
merge_data0<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge.txt",header =T,sep="\t") 

Reduce(intersect,list(unique(merge_data0$gene), common_freq$Var1))#88
Reduce(intersect,list(unique(merge_data0$gene), common_freq_active_high$Var1))#70

venn_Abortion_DEGs_daTFs1<-venn.diagram(list(All_DEGs=unique(merge_data0$gene),All_sig_DATF= common_freq$Var1),
                                        alpha=c(0.7,0.7),
                                        lwd=1,lty=1, #lwd用于设定圆弧的宽度，lty用于设定圆弧的线型
                                        col="black" , #指定图形的圆周边缘颜色  transparent 透明   
                                        fill=ppCor[c(5:6)], #参数fill表示各个集合对应的圆的填充颜色,
                                        cex = 1.5,    #每个区域label名称的大小
                                        cat.col=ppCor[c(5:6)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                        # cat.fontface=4,  #字体格式
                                        cat.cex = 1.5,      #每个分类名称大小
                                        main = paste0(cell_name,":p<0.05(unselect)"),
                                        #  main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                        filename = NULL)

venn_Abortion_DEGs_daTFs2<-venn.diagram(list(All_DEGs=unique(merge_data0$gene),All_sig_DATF= common_freq_active_high$Var1),
                                        alpha=c(0.7,0.7),
                                        lwd=1,lty=1, #lwd用于设定圆弧的宽度，lty用于设定圆弧的线型
                                        col="black" , #指定图形的圆周边缘颜色  transparent 透明   
                                        fill=ppCor[c(5:6)], #参数fill表示各个集合对应的圆的填充颜色,
                                        cex = 1.5,    #每个区域label名称的大小
                                        cat.col=ppCor[c(5:6)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                        #cat.fontface=4,  #字体格式
                                        cat.cex = 1.5,      #每个分类名称大小
                                        main = paste0(cell_name,":p<0.05(active_high)"),
                                        # main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                        filename = NULL)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Venn_between_dacTFs_scDEGs_number_for_Abortion.pdf",width =10,height = 8)
grid.newpage(); grid.draw(venn_Abortion_DEGs_daTFs1)
grid.newpage();grid.draw(venn_Abortion_DEGs_daTFs2)
dev.off()

head(D_ac_TFs_active_high_matrix)
com_trend_freq_active_high<-data.frame(table(D_ac_TFs_active_high_matrix$label,D_ac_TFs_active_high_matrix$change))
com_trend_freq_active_high<-com_trend_freq_active_high[order(com_trend_freq_active_high$Var2,com_trend_freq_active_high$Freq,decreasing = T),]
dim(com_trend_freq_active_high)## 728     3
com_trend_freq_active_high<-com_trend_freq_active_high[which(com_trend_freq_active_high$Freq>0),]
dim(com_trend_freq_active_high)##500    3
range(com_trend_freq_active_high$Freq)# 1 10

Up_acTF<-as.character(com_trend_freq_active_high[which(com_trend_freq_active_high$Var2 == "Up"),"Var1"])
Down_acTF<-as.character(com_trend_freq_active_high[which(com_trend_freq_active_high$Var2 == "Down" ),"Var1"])

#for common TFs
length(Up_acTF);length(Down_acTF)#290 210
length(Reduce(intersect,list(Up_acTF,Down_acTF)))#136

#plot different types of TFs with different number among Cells
dim(common_freq);dim(common_freq_active_high)#539   2 #364    2

Distribution1<-ggbarplot(data.frame(table(common_freq$Freq)), x="Var1", y="Freq",fill="Pink",color = "white")+#,rotate=TRUE
  geom_text(aes(label=Freq),size=2,position = position_dodge(0.5),vjust=-0.5)+
  labs(y="TFs Number",x="Cell type Number",title="Distribution of Abortion TFs among Cells ")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 10),axis.ticks = element_blank(),
        legend.position ="top",legend.direction = "horizontal")
Distribution2<-ggbarplot(data.frame(table(common_freq_active_high$Freq)), x="Var1", y="Freq",fill="navy",color = "white")+#,rotate=TRUE
  geom_text(aes(label=Freq),size=2,position = position_dodge(0.5),vjust=-0.5)+
  labs(y="TFs Number",x="Cell type Number",title="Distribution of Abortion active_high TFs among Cells ")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 10),axis.ticks = element_blank(),
        legend.position ="top",legend.direction = "horizontal")
Distribution1
Distribution2

head(com_trend_freq_active_high)
TFs_freq<-data.frame(table(com_trend_freq_active_high$Freq,com_trend_freq_active_high$Var2))
colnames(TFs_freq)<-c("cell_types_number","trend","TFs_type_freq")
TFs_freq$trend<-factor(TFs_freq$trend,levels=c("Up","Down"))

TF_number_disribution<-ggplot(data=TFs_freq, mapping=aes(x= cell_types_number,y=TFs_type_freq,fill=trend))+
  geom_bar(stat="identity",width=0.8,position= 'dodge')+
  scale_fill_manual(name="Active trend",values=c(ppCor_all2[21],ppCor_all2[22]))+
  geom_text(aes(label=TFs_type_freq),size=3,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="TFs type Number",x="Cell type Number",title="Distribution of Abortion active_high TFs among Cells ")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")# + coord_flip()
TF_number_disribution
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/barplot_Up_Down_common_TFs_number_order_in_common_celltype_number.pdf",TF_number_disribution,width=10, height=8)

#top common genes
head(com_trend_freq_active_high)
#common_freq<-data.frame(table(com_trend_freq_active_high$Var1,TFs_freq$trend))
common_freq2 <-dcast(com_trend_freq_active_high, Var1  ~ Var2 , value.var = "Freq")
common_freq2[is.na(common_freq2)]<-0
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
plot_Up<-ggplot(data_Up, aes(x = Var1, y =  Up)) +  geom_bar(stat = "identity", alpha = 0.8,fill=ppCor_all2[21]) + coord_polar() +
  geom_text(aes(x= Var1, y= Up+2, label =  Up),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of cell subtype for top 10 common Up_dacTFs") + 
  theme(title = element_text(vjust = -56, face = "bold"))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
plot_down<-ggplot(data_down, aes(x = Var1, y =  Down)) +  geom_bar(stat = "identity", alpha = 0.8,fill=ppCor_all2[22]) + coord_polar() +
  geom_text(aes(x= Var1, y= Down+1, label =  Down),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of cell subtype for top 10 common Down_dacTFs") + 
  theme(title = element_text(vjust = -56, face = "bold"))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
plot_Up
plot_down
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Top10_common_number_for_Abortion_related_dacTFs.pdf",width =12,height = 6)
grid.arrange(plot_Up, plot_down,ncol=2)
dev.off()

#plot for selected dacTFs
gene_target<-unique(c(as.character(data_down$Var1),as.character(data_Up$Var1)))
data_target<-com_trend_freq_active_high[which(com_trend_freq_active_high$Var1 %in% gene_target),]
data_target$Var1<-factor(data_target$Var1,levels=rev(gene_target))
data_target$Var2<-factor(data_target$Var2,levels=c("Up","Down"))

data_target

Top10_dacTFs<-ggplot(data_target, aes(x = Var1, y = Freq, fill = Var2)) +  
  geom_bar(stat = "identity", alpha = 0.7) + coord_polar() +
  geom_text(aes(x= Var1, y= 10, label =  Freq),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of Cells type for daTFs(Top10)") + 
  theme(legend.position = "right") +  
  theme(title = element_text(vjust = -56, face = "bold"))+
  scale_fill_manual(name="Active trend",values=c(ppCor_all2[21],ppCor_all2[22]))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
Top10_dacTFs
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/SF4f_rose_plot_top10_number_cell_type_Up_Down_common_TFs.pdf",Top10_dacTFs,width=10, height=10)
