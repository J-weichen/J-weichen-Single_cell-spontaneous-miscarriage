library(ggplot2)
library(ggpubr)
library(grid)
library(ggsci)
library(data.table)
library(gridExtra)
library(cowplot)
library(scales)
library(reshape2)
library(stringr)
#SET colors
pal <- pal_npg("nrc", alpha=0.6)(9)
nejm<-pal_nejm("default",alpha = 0.6)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)


row_anno2 <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/DEGs_EVT2_EVT3_overlapped_development_genes_along_sleceted_EVT_branchment_pseudotime.txt", header=T)
head(row_anno2)
table(row_anno2$Cluster)
#cluster1 cluster2 cluster3 
#  153       58      366

plot_data0<-data.frame(table(row_anno2$type))
plot_data0$Cluster<-unlist(lapply(strsplit(as.character(plot_data0$Var1),"-"), function(x) x[1]))
plot_data0$class<-unlist(lapply(strsplit(as.character(plot_data0$Var1),"-"), function(x) x[2]))
head(plot_data0)

p_ring_list<-list();cnane_region<-c()

for ( region_name in unique(plot_data0$Cluster)){
  #  region_name<-"LTR"  ##test line 
  print(as.character(region_name))
  cluster_number<-region_name
 # cluster_number<-"cluster2"
  plot_data<-plot_data0[which(plot_data0$Cluster ==cluster_number),]
  plot_data$fraction<-plot_data$Freq/sum(plot_data$Freq)
  plot_data$ymax<-cumsum(plot_data$fraction)
  plot_data$ymin<-c(0,head(plot_data$ymax,n=-1))
  plot_data$labelPosition<-(plot_data$ymax + plot_data$ymin)/2
  plot_data$label<-paste0(plot_data$class,"\n number: ",plot_data$Freq)
  p_ring <- ggplot(plot_data,aes(ymax=ymax,ymin=ymin,xmax=4,xmin=3))+
  geom_rect(aes(fill=class))+geom_label(x=3.5,aes(y=labelPosition,label=label),size=4)+
  #scale_fill_brewer(palette = 4)+
  scale_fill_manual(values=my_morandi_colors)+
  coord_polar(theta = "y")+  xlim(2,4)+theme_void()#+theme(legend.position = "none")
  p_ring_list<-c(p_ring_list,list(p_ring))
  cnane_region<-c(cnane_region,region_name)
}

length(p_ring_list);length(cnane_region)
names(p_ring_list)<-cnane_region

p_ring_group <- ggpubr::ggarrange(p_ring_list[[1]],p_ring_list[[2]],p_ring_list[[3]],nrow = 3, ncol = 1)
                                  #labels = c('A', 'B', 'C', 'D'),  font.label = list(color = 'black')
p_ring_group2<-annotate_figure(p_ring_group, top=text_grob("number of DEGs in diff classes", color = "black",face = "bold", size=12))
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/ring_plot_for_DEGs_in_diff_class.pdf",p_ring_group2,width = 6, height =19)
