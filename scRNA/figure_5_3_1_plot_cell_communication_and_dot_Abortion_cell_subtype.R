#https://www.jianshu.com/p/f196c98e0954
rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(Seurat)
library(reshape2)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
show_col(allcolour)
conflict_prefer("dcast", "reshape2")

Used_cell2<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Mycs","HCs","STCs","FBs","Endo","NKs","Ts")
#cellphone 结果的可视化:在这一步中，我一般只用到上述的means.txt和pvalues.txt文件
#非对称热图：将A-->B 和B-->A之间的互作数目分别讨论，不合并
#reading network from cellpheno
pheno_file='/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/Abortion/' #outs下的文件放在这里了
pvalues <- read.delim(paste0(pheno_file,"pvalues.txt"), check.names = FALSE,header = T,sep = "\t",stringsAsFactors = F) #读取数据
head(pvalues);dim(pvalues)#1199  372
pvalues=pvalues[,12:dim(pvalues)[2]] #此时不关注前11列
statdf=as.data.frame(colSums(pvalues < 0.05)) #统计在某一种细胞pair的情况之下，显著的受配体pair的数目；阈值可以自己选
colnames(statdf)=c("number")
head(statdf)
#排在前面的分子定义为indexa；排在后面的分子定义为indexb
statdf$indexb<-unlist(lapply(strsplit(as.character(rownames(statdf)),"[|]"), function(x) x[2]))
statdf$indexa<-unlist(lapply(strsplit(as.character(rownames(statdf)),"[|]"), function(x) x[1]))
#selected cells in target cell
statdf<-statdf[which(statdf$indexa %in% Used_cell2 & statdf$indexb %in% Used_cell2),]
head(statdf);dim(statdf)#225   3

#设置合适的细胞类型的顺序
rankname=Used_cell2
#转成因子类型，画图时，图形将按照预先设置的顺序排列
statdf$indexa=factor(statdf$indexa,levels = rev(rankname))
statdf$indexb=factor(statdf$indexb,levels = rev(rankname))
range(statdf$number)# 1 163
Abortion_heatmap<-ggplot(statdf,aes(x=indexa,y=indexb,fill=number))+geom_tile(color="white")+
  geom_text(aes(label = round(number, 1))) +
  #scale_fill_gradientn(colours = c("blue","lightblue","orange","red","darkred"),limits=c(0,180))+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,180))+
  scale_x_discrete("cluster 1 produces molecule 1")+
  scale_y_discrete("cluster 2 produces molecule 2")+
  theme_minimal()+labs(title = "Abortion interaction number")+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
    panel.grid = element_blank())
Abortion_heatmap
ggsave(Abortion_heatmap,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_interaction_number_no_merge.pdf",width = 15,height = 14)
#对称热图：将两细胞之间的互作分子对数目合并
statdf$total_number=0
for (i in 1:dim(statdf)[1]) {
  tmp_indexb=statdf[i,"indexb"]
  tmp_indexa=statdf[i,"indexa"]
  if (tmp_indexa == tmp_indexb) {
    statdf[i,"total_number"] = statdf[i,"number"]
  } else {
    statdf[i,"total_number"] = statdf[statdf$indexb==tmp_indexb & statdf$indexa==tmp_indexa,"number"]+
      statdf[statdf$indexa==tmp_indexb & statdf$indexb==tmp_indexa,"number"]
  }
}
head(statdf)
rankname=Used_cell2
statdf$indexa=factor(statdf$indexa,levels =rev(rankname))
statdf$indexb=factor(statdf$indexb,levels =rev(rankname))
range(statdf$total_number)
Abortion_heatmap2<-ggplot(statdf,aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  geom_text(aes(label = round(total_number, 1))) +
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,250))+
  scale_x_discrete("cluster 1")+scale_y_discrete("cluster 2")+
  theme_minimal()+labs(title = "Abortion interaction number")+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
    panel.grid = element_blank()
  )
Abortion_heatmap2
ggsave(Abortion_heatmap2,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_interaction_number_merge.pdf",width = 15,height = 14)

#获取上三角热图
statdf_dc<-dcast(statdf[,c("indexb","indexa","total_number")],indexa~ indexb)
rownames(statdf_dc)<-statdf_dc$indexa;statdf_dc<-statdf_dc[,-1]
head(statdf_dc)
#获取相关性矩阵的上三角矩阵
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(statdf_dc)
dim(upper_tri)
upper_tri$indexa<-rownames(upper_tri)
statdf_dc_tri <- melt(upper_tri,"indexa",na.rm = TRUE)
rankname=rev(Used_cell2)
statdf_dc_tri$indexa=factor(statdf_dc_tri$indexa,levels = rankname)
statdf_dc_tri$variable=factor(statdf_dc_tri$variable,levels = rankname)
head(statdf_dc_tri)
range(statdf_dc_tri$value)# 1 247
Abortion_heatmap3<-ggplot(statdf_dc_tri, aes(indexa, variable)) + 
  geom_tile(aes(fill = value),colour = "white") +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,250))+
  scale_x_discrete("cluster 1")+scale_y_discrete("cluster 2")+
  theme_bw()+labs(title = "Abortion interaction number")+
  theme(
    axis.text.x.bottom  = element_text(hjust = 1, vjust = NULL, angle = 45),
    panel.grid = element_blank()
  )
Abortion_heatmap3
ggsave(Abortion_heatmap3,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_interaction_number_merge_up_tri.pdf",width = 15,height = 14)


#网络图表示合并后细胞互作关系的数量
mynet=statdf_dc_tri[statdf_dc_tri$value > 0,] #过滤掉值为0的观测
head(mynet);dim(mynet)
rankname=Used_cell2
mynet$indexa=factor(mynet$indexa,levels = rankname)
mynet<-mynet[order(mynet$indexa),]
head(mynet) 
length(unique(as.character(mynet$indexa)))
net<- graph_from_data_frame(mynet) #构建net对象
#plot network by R
length(allcolour)
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
karate_groups <- cluster_optimal(net) #统计每个端点的和
coords <- layout_in_circle(net, order = order(membership(karate_groups)))  # 设置网络布局

#设置节点和连线的颜色
color1=allcolour[1:length(rankname)]
names(color1)=rankname
color2=colorRampPalette(brewer.pal(9, "Reds")[3:7])(10)#将颜色分成多少份，取决于互作关系数目的最大值
names(color2)=1:10 #每一份颜色用对应的数字命名
range(E(net)$value/20)
#设置网络性质
E(net)$width  <- E(net)$value/20  #根据count值设置边的宽 
E(net)$label = E(net)$value #连线的标注
E(net)$color <- color2[as.character(ifelse(E(net)$width > 10,10,ceiling(E(net)$width)))] #用前面设置好的颜色赋给连线，颜色深浅对应数值大小
E(net)$label = E(net)$value #连线的标注
E(net)$label.color <- "black" #连线标注的颜色
V(net)$label.color <- "black" #节点标注的颜色
V(net)$color <- color1[names(V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息

#调整节点位置的线条角度
##如果没有这两行代码，节点位置的圆圈是向右的
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_subcell_net_merge_plot.pdf",width=12,height=12)
plot(net,
     edge.arrow.size = 0, #连线不带箭头
     edge.curved = 0.2, #连线不弯曲
     vertex.frame.color = "black", #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = 30) #节点大小
dev.off()


###网络图和贝壳图表示非合并细胞互作关系的数量
#selected cells in target cell
mynet=statdf[statdf$number > 0,c("indexa","indexb","number")] #过滤掉值为0的观测
colnames(mynet)<-c("SOURCE","Target","count")
head(mynet);dim(mynet)#225   3
table(mynet$count)
mynet$SOURCE<-factor(mynet$SOURCE,levels = Used_cell2)
mynet<-mynet[order(mynet$SOURCE),]
head(mynet) 
length(unique(as.character(mynet$SOURCE)))

net<- graph_from_data_frame(mynet) #构建net对象
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
karate_groups <- cluster_optimal(net) #统计每个端点的和
coords <- layout_in_circle(net, order = order(membership(karate_groups)))  # 设置网络布局
E(net)$width  <- E(net)$count/10  #根据count值设置边的宽 
range(E(net)$width)
#调整节点位置的线条角度 ##如果没有这两行代码，节点位置的圆圈是向右的
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

net2 <- net; mynet2<-mynet# 复制一份备用

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_subcell_not_merge_network_R_plot.pdf",width=12,height=12)
mynet$SOURCE<-as.character(mynet$SOURCE)
for (i in 1: length(unique(mynet$SOURCE)) ){ #配置发出端的颜色
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x)) })%>% unlist()]$color <- allcolour[i] }  # 这波操作谁有更好的解决方案？ 

plot(net, edge.arrow.size=.1,  #箭头大小设置为0.1
     edge.curved=0.2, # 只是调了曲率这个参数
     vertex.color=allcolour,
     vertex.frame.color="#555555", #圆圈颜色
     vertex.label.color="black", #标签颜色
     layout = coords, #网络布局位点
     vertex.label.cex=.7 )#标签大小) 
dev.off()

#贝壳图
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Abortion_subcell_net_nomerge_plot_sperated.pdf",width=30,height=30)
length(unique(mynet2$SOURCE))#15
par(mfrow=c(6,3), mar=c(.3,.3,.3,.3))
for (i in 1: length(unique(as.character(mynet$SOURCE))) ){
  net1<-net2
  E(net1)$count <- ""
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  # 故技重施
  
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.4,
       edge.label = E(net1)$count, # 绘制边的权重
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1
  ) 
  
}
dev.off()


##气泡图——具体的互作关系--未细化
##以上几种图，只是用来展示数量，具体的两种细胞之间的互作关系可以用如下的代码展示：
source("CCC.bubble.R")
CCC(
  pfile="./test/pvalues.txt",
  mfile="./test/means.txt",
  #neg_log10_th= -log10(0.05),
  #means_exp_log2_th=1,
  #neg_log10_th2=3,
  #means_exp_log2_th2=c(-4,6),
  #notused.cell=c("Bcell","Gcell"),
  #used.cell=c("Mcell"),
  #cell.pair=c("Mcell.Scell","Mcell.NKcell","Mcell.Tcell","Scell.Mcell","NKcell.Mcell","Tcell.Mcell"),#这里是自定义的顺序，若是可选细胞对的子集，则只展示子集，若有交集则只展示交集；空值情况下，会根据可选细胞对自动排序
  #gene.pair=c("MIF_TNFRSF14","FN1_aVb1 complex","EGFR_MIF")#作用同上
)
ggsave(filename = "interaction.detail.1.pdf",device = "pdf",width =20,height = 12,units = "cm")

#参数解释：
#neg_log10_th和means_exp_log2_th两个参数用来筛选显著的互作关系；
#neg_log10_th2和means_exp_log2_th2两个参数用来限定最终气泡图的数值范围；
#notused.cell不包含的细胞类型；
#used.cell必须包含的细胞类型；
#cell.pair必须包含的细胞pair，以及它们的顺序；
#gene.pair必须包含的基因pair，以及它们的顺序。
#后面四个参数在细化气泡图的时候，很有用。
CCC(
  pfile="./test/pvalues.txt",
  mfile="./test/means.txt",
  cell.pair=c("Mcell.Scell","Mcell.NKcell","Mcell.Tcell","Scell.Mcell","NKcell.Mcell","Tcell.Mcell"),#这里是自定义的顺序，若是可选细胞对的子集，则只展示子集，若有交集则只展示交集；空值情况下，会根据可选细胞对自动排序
  gene.pair=c("MIF_TNFRSF14","FN1_aVb1 complex","EGFR_MIF")#作用同上
)
ggsave(filename = "interaction.detail.2.pdf",device = "pdf",width =16,height = 10,units = "cm")
