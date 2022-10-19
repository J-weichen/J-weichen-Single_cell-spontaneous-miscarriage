rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(future)
plan(strategy = "multicore", workers = 24)
#富集分析
library(stringr)
library(ggplot2)
library(ggpubr) 
library(ggsci)
library(scales)
library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)#org.Hs.eg.db 的数据类型以及使用简介：可用于数据类型转换
keytypes(org.Hs.eg.db)#查看org.Hs.eg.db数据对象里面包含着各大主流数据库的数据

pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
show_col(pal)
pal
#[1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF"
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
show_col(nejm)
nejm
#[1] "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"
#decide color
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)
#for selected  
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/enrichment_result/sc_bulk_abortion_related_overlap_DEG_GO_BP_selected.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
Enrich_all_plot_select$GeneList<-factor(Enrich_all_plot_select$GeneList,levels = c("down","up"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$GeneList,decreasing = T ),]
Enrich_all_plot_select$GO<-factor(Enrich_all_plot_select$GO,levels = rev(c(unique(as.character(Enrich_all_plot_select$GO)))))
range(Enrich_all_plot_select$LogP)#4 22
range(Enrich_all_plot_select$GeneInHitList)# 11 41
Enrich_all_plot_select$LogP2<-Enrich_all_plot_select$LogP
Enrich_all_plot_select[which(Enrich_all_plot_select$LogP>10),]$LogP2<-10

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=GeneList,y=GO,size=GeneInHitList,colour=LogP2))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,10,20,30,40),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,2,4,6,8,10),name='-log(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="Abortion related bulk single cell overlapped DEGs: metascape BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/default_DEGs/enrichment_result/bulk_single_cell_overlapped_DEGs_Enrich_all_BP_select.pdf",plot_BP,width = 8, height =8)

#未跑
p <- ggplot(Enrich_all_plot_select,aes(LogP,reorder(GO,LogP)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=X.GeneInGOAndHitList,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log(p.adjust)",y="Description",title="KEGG enrichment")+
  scale_size_continuous(range=c(8,28))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=55))+
  geom_text(aes(label=qlogtran), size=4)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")

p

