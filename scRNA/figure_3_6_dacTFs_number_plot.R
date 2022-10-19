#for plot evaluation
cell_name="ALLcell_all_genes"
#target_final_group="Treat"
D_ac_TF_matrix <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/merge_",cell_name,"_D_ac_TF_matrix.txt"), header=T)
D_ac_TFs_active_high_matrix <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/Def_ac_TFs/merge_",cell_name,"_D_ac_TFs_active_high_matrix.txt"), header=T)

range(table(D_ac_TF_matrix$label))#1 17
range(table(D_ac_TFs_active_high_matrix$label))#1 14
head(D_ac_TF_matrix)
head(D_ac_TFs_active_high_matrix)
##formal plot for major figure three 

## De-TFs number for each Cells
num_cell_plot<-D_ac_TFs_active_high_matrix[,c("cell_identity","change")]
num_cell_plot$Num<-1
num_cell_plot2<-aggregate(num_cell_plot$Num, by=list(num_cell_plot$cell_identity,num_cell_plot$change), FUN=sum)
#table(D_ac_TFs_active_high_matrix$cell_identity,D_ac_TFs_active_high_matrix$change)
colnames(num_cell_plot2)<-c("cell_identity","tendency","TFs_num")
cell_order<-c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2",
              "STCs","FBs","Endo","Epi","HCs","MyCs","NKs","Ts","Bs","Masts")
num_cell_plot2$cell_identity<-factor(num_cell_plot2$cell_identity,levels=cell_order)
num_cell_plot2$tendency<-factor(num_cell_plot2$tendency,levels=c("Up","Down"))
str(num_cell_plot2)
#rose maps 
ggplot(num_cell_plot2, aes(x = cell_identity, y = TFs_num, fill = tendency)) +  geom_bar(stat = "identity", alpha = 0.7) + coord_polar() +
  theme_grey() + labs(x = "", y = "", title = "The number of TFs for each cell subtype") + 
  # geom_text(aes(y = TFs_num/2 + max(TFs_num)/4, label = TFs_num, color = tendency), size = 5) +    ## 加上数字
  # theme(axis.text.y = element_blank()) +    theme(axis.ticks = element_blank()) +      ## 去掉左上角的刻度标签## 去掉左上角的刻度线
  # theme(panel.border = element_blank()) +   ## 去掉外层边框
  theme(legend.position = "right") +   ## 去掉图例
  theme(title = element_text(vjust = -56, face = "bold"))+
  scale_fill_manual(name="Active trend",values=ppCor_all[1:2])

num_cell_plot3 <-dcast(num_cell_plot2, cell_identity ~ tendency , value.var = "TFs_num")
num_cell_plot3[is.na(num_cell_plot3)]<-0
num_cell_plot3$TFs_total_num<-num_cell_plot3$Up+num_cell_plot3$Down
num_cell_plot3<-num_cell_plot3[order(num_cell_plot3$TFs_total_num,decreasing = T),]
num_cell_plot3$cell_identity<-factor(num_cell_plot3$cell_identity,levels = as.character(num_cell_plot3$cell_identity))
num_cell_plot4<-melt(num_cell_plot3[,c(1:3)],)
cell_order<-c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2",
              "STCs","FBs","Endo","Epi","HCs","MyCs","NKs","Ts","Bs","Masts")
num_cell_plot4$cell_identity<-factor(num_cell_plot4$cell_identity,levels=cell_order)
num_cell_plot4$variable<-factor(num_cell_plot4$variable,levels=c("Up","Down"))
str(num_cell_plot4);head(num_cell_plot4)
ggplot(data=num_cell_plot4, mapping=aes(x= cell_identity,y=value,fill=variable))+
  geom_bar(stat="identity",width=0.8,position= 'dodge')+
  scale_fill_manual(name="Active trend",values=ppCor[1:2])+
  geom_text(aes(label=value),size=3,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="def_acTFs Number",x="Cell types",title="The number of differential AMA active_high TFs in each Cells ")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 10),  axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")


#进一步统计组间重复
common_freq<-data.frame(table(D_ac_TF_matrix$label))
common_freq<-common_freq[order(common_freq$Freq,decreasing = T),]

common_freq_active_high<-data.frame(table(D_ac_TFs_active_high_matrix$label))
common_freq_active_high<-common_freq_active_high[order(common_freq_active_high$Freq,decreasing = T),]

length(unique(common_freq$Var1));length(unique(common_freq_active_high$Var1))##538  314   
#for common TFs
Reduce(intersect,list(unique(merge_data_rm0$gene), common_freq$Var1))#106
Reduce(intersect,list(unique(merge_data_rm0$gene), common_freq_active_high$Var1))#78

grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(All_DEGs=unique(merge_data_rm0$gene),All_sig_DATF= common_freq$Var1),
                                   alpha=c(0.7,0.7),
                                   lwd=1,lty=1, #lwd用于设定圆弧的宽度，lty用于设定圆弧的线型
                                   col="black" , #指定图形的圆周边缘颜色  transparent 透明   
                                   fill=ppCor[c(5:6)], #参数fill表示各个集合对应的圆的填充颜色,
                                   cex = 1.5,    #每个区域label名称的大小
                                   cat.col=ppCor[c(5:6)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,  #字体格式
                                   cat.cex = 1.5,      #每个分类名称大小
                                   main = paste0(cell_name,":p<0.01(unselect)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)
grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(All_DEGs=unique(merge_data_rm0$gene),All_sig_DATF= common_freq_active_high$Var1),
                                   alpha=c(0.7,0.7),
                                   lwd=1,lty=1, #lwd用于设定圆弧的宽度，lty用于设定圆弧的线型
                                   col="black" , #指定图形的圆周边缘颜色  transparent 透明   
                                   fill=ppCor[c(5:6)], #参数fill表示各个集合对应的圆的填充颜色,
                                   cex = 1.5,    #每个区域label名称的大小
                                   cat.col=ppCor[c(5:6)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,  #字体格式
                                   cat.cex = 1.5,      #每个分类名称大小
                                   main = paste0(cell_name,":p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)
grid.newpage(); #清空画板，开始画新图

com_trend_freq_active_high<-data.frame(table(D_ac_TFs_active_high_matrix$label,D_ac_TFs_active_high_matrix$change))
com_trend_freq_active_high<-com_trend_freq_active_high[order(com_trend_freq_active_high$Var2,com_trend_freq_active_high$Freq,decreasing = T),]
dim(com_trend_freq_active_high)## 628   3
com_trend_freq_active_high<-com_trend_freq_active_high[which(com_trend_freq_active_high$Freq>0),]
dim(com_trend_freq_active_high)##401   3
range(com_trend_freq_active_high$Freq)# 1 13

Up_acTF<-as.character(com_trend_freq_active_high[which(com_trend_freq_active_high$Var2 == "Up"),"Var1"])
Down_acTF<-as.character(com_trend_freq_active_high[which(com_trend_freq_active_high$Var2 == "Down" ),"Var1"])

#for common TFs
length(Up_acTF);length(Down_acTF)#235 166
length(Reduce(intersect,list(Up_acTF,Down_acTF)))#87

grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(Up_acTF=Up_acTF,Down_acTF= Down_acTF),
                                   alpha=c(0.7,0.7),lwd=1,lty=1,col="black" ,  fill=ppCor[c(1:2)], 
                                   cex = 1.5, cat.col=ppCor[c(1:2)],cat.fontface=4,  cat.cex = 1.5,    
                                   main = paste0(cell_name,":",":p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)

Reduce(intersect,list(sc_Up_gene,sc_Down_gene,Up_acTF,Down_acTF))
#"FOS"    "JUND"   "HMGB3"  "FOSB"   "SOX4"   "TFDP2"  "GATA2"  "TFAP2A"
Reduce(intersect,list(sc_Up_gene,Up_acTF,Down_acTF))
#[1] "FOS"    "JUND"   "HMGB3"  "FOSB"   "CEBPB"  "ATF3"   "KLF6"   "SOX4"   "ASCL2"  "GATA3"  "STAT2" 
#[12] "TFDP2"  "FOXP1"  "REL"    "GATA2"  "SNAI1"  "TFAP2A" "ELF3"   "ESRRA"  "NR3C1"  "HIF1A"  "ELF1" 
Reduce(intersect,list(sc_Up_gene,Up_acTF))
Reduce(intersect,list(sc_Down_gene,Up_acTF,Down_acTF))
# [1] "MYCN"    "TFDP2"   "BHLHE40" "TFAP2A"  "UGP2"    "FOS"     "HMGB3"   "GATA2"   "NR2F6"  
#[10] "SOX4"    "JUND"    "BCLAF1"  "FOSB"   
Reduce(intersect,list(sc_Down_gene,Down_acTF))

setdiff(Reduce(intersect,list(sc_Up_gene,Up_acTF)),c(sc_Down_gene,Down_acTF))
setdiff(Reduce(intersect,list(sc_Down_gene,Down_acTF)),c(sc_Up_gene,Up_acTF))

grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(sc_Up_gene=sc_Up_gene,sc_Down_gene=sc_Down_gene,Up_DATF= Up_acTF,Down_DATF=Down_acTF ),
                                   alpha=c(0.9,0.9,0.9,0.9),
                                   lwd=1,lty=1, #lwd用于设定圆弧的宽度，lty用于设定圆弧的线型
                                   col="black" , #指定图形的圆周边缘颜色  transparent 透明   
                                   fill=ppCor[c(1,2,3,5)], #参数fill表示各个集合对应的圆的填充颜色,
                                   cex = 1.5,    #每个区域label名称的大小
                                   cat.col=ppCor[c(1,2,3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,  #字体格式
                                   cat.cex = 1.5,      #每个分类名称大小
                                   main = paste0(cell_name,"::DEG_005_default_and_DE_acTFq005_p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)
grid.newpage(); #清空画板，开始画新图
#reading DEGs from bulk
bulk_merge_data<-read.table(file="/mnt/data/chenwei/gongchen/2.map_result/count_file/4.count/Abortion_vs_CTRL.DEG_information_pvalue005_FC1.5.txt",header =T,sep="\t") 
head(bulk_merge_data);dim(bulk_merge_data)#5316    4
bulk_Up_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange > 0),]$ID))
length(bulk_Up_gene)# 2748
bulk_Down_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange < 0),]$ID))
length(bulk_Down_gene)# 2568

bulk_sc_Up_DEGs<- setdiff(Reduce(intersect,list(bulk_Up_gene,sc_Up_gene)),c(sc_Down_gene,bulk_Down_gene))
bulk_sc_Down_DEGs<-setdiff(Reduce(intersect,list(sc_Down_gene,bulk_Down_gene)),c(bulk_Up_gene,sc_Up_gene))
setdiff(Reduce(intersect,list(bulk_sc_Up_DEGs,Up_acTF)),c(bulk_sc_Down_DEGs,Down_acTF))
setdiff(Reduce(intersect,list(bulk_sc_Down_DEGs,Down_acTF)),c(bulk_sc_Up_DEGs,Up_acTF))

grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_genAge<-venn.diagram(list(bulk_sc_Up_DEGs=bulk_sc_Up_DEGs,bulk_sc_Down_DEGs=bulk_sc_Down_DEGs,Up_DATF= Up_acTF,Down_DATF=Down_acTF ),
                                   alpha=c(0.9,0.9,0.9,0.9),
                                   lwd=1,lty=1, #lwd用于设定圆弧的宽度，lty用于设定圆弧的线型
                                   col="black" , #指定图形的圆周边缘颜色  transparent 透明   
                                   fill=ppCor[c(1,2,3,5)], #参数fill表示各个集合对应的圆的填充颜色,
                                   cex = 1.5,    #每个区域label名称的大小
                                   cat.col=ppCor[c(1,2,3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,  #字体格式
                                   cat.cex = 1.5,      #每个分类名称大小
                                   main = paste0(cell_name,"::DEG_005_default_and_DE_acTFq005_p<0.01(active_high)"),
                                   main.cex = 1.5, main.fontface = 1.5, main.fontfamily =1.5, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)
                        
#plot different types of TFs with different number among Cells
dim(common_freq);dim(common_freq_active_high)#538   2 #314    2

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
Distribution1+Distribution2

TFs_freq<-data.frame(table(com_trend_freq_active_high$Freq,com_trend_freq_active_high$Var2))
colnames(TFs_freq)<-c("cell_types_number","trend","TFs_type_freq")
TFs_freq$trend<-factor(TFs_freq$trend,levels=c("Up","Down"))

ggplot(data=TFs_freq, mapping=aes(x= cell_types_number,y=TFs_type_freq,fill=trend))+
  geom_bar(stat="identity",width=0.8,position= 'dodge')+
  scale_fill_manual(name="Active trend",values=ppCor[1:2])+
  geom_text(aes(label=TFs_type_freq),size=3,position = position_dodge(1),vjust=0.5)+
  theme_classic()+labs(y="TFs type Number",x="Cell type Number",title="Distribution of Abortion active_high TFs among Cells ")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")# + coord_flip()