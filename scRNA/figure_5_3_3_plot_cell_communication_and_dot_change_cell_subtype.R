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
#conflict_prefer("dcast", "reshape2")

Used_cell2<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Mycs","HCs","STCs","FBs","Endo","NKs","Ts")
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)}

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)}
#reading network from cellpheno
pheno_file='/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/Abortion/' #outs下的文件放在这里了
pvalues <- read.delim(paste0(pheno_file,"pvalues.txt"), check.names = FALSE,header = T,sep = "\t",stringsAsFactors = F) #读取数据
head(pvalues);dim(pvalues)# 1199  411
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05)) 
colnames(statdf)=c("number")
head(statdf)

#排在前面的分子定义为indexa；排在后面的分子定义为indexb
statdf$indexb<-unlist(lapply(strsplit(as.character(rownames(statdf)),"[|]"), function(x) x[2]))
statdf$indexa<-unlist(lapply(strsplit(as.character(rownames(statdf)),"[|]"), function(x) x[1]))
#selected cells in target cell
statdf<-statdf[which(statdf$indexa %in% Used_cell2 & statdf$indexb %in% Used_cell2),]
head(statdf);dim(statdf)#225   3
statdf_Abortion<-statdf
colnames(statdf_Abortion)=c("number_Abortion","indexb_Abortion","indexa_Abortion")

#reading network from cellpheno of CTRL
pheno_file='/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/2.statistical_result/CTRL/' #outs下的文件放在这里了
pvalues <- read.delim(paste0(pheno_file,"pvalues.txt"), check.names = FALSE,header = T,sep = "\t",stringsAsFactors = F) #读取数据
head(pvalues);dim(pvalues)# 1199  411
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05)) 
colnames(statdf)=c("number")
head(statdf)

#排在前面的分子定义为indexa；排在后面的分子定义为indexb
statdf$indexb<-unlist(lapply(strsplit(as.character(rownames(statdf)),"[|]"), function(x) x[2]))
statdf$indexa<-unlist(lapply(strsplit(as.character(rownames(statdf)),"[|]"), function(x) x[1]))
#selected cells in target cell
statdf<-statdf[which(statdf$indexa %in% Used_cell2 & statdf$indexb %in% Used_cell2),]
head(statdf);dim(statdf)#225   3
statdf_CTRL<-statdf
colnames(statdf_CTRL)=c("number_CTRL","indexb_CTRL","indexa_CTRL")

#合并数值
head(statdf_CTRL);head(statdf_Abortion)
statdf_nomerge<-merge(statdf_CTRL,statdf_Abortion,by=0, all = T)
head(statdf_nomerge)
statdf_nomerge$diff_count<-statdf_nomerge$number_Abortion-statdf_nomerge$number_CTRL
rownames(statdf_nomerge)<-statdf_nomerge$Row.names;statdf_nomerge<-statdf_nomerge[,-1]

statdf_nomerge<-statdf_nomerge[,c("indexa_CTRL","indexb_CTRL","diff_count")]
colnames(statdf_nomerge)<-c("indexa","indexb","diff_count")
head(statdf_nomerge);dim(statdf_nomerge)

statdf_nomerge$total_number=0
for (i in 1:dim(statdf_nomerge)[1]) {
  #i=181
  tmp_indexb=as.character(statdf_nomerge[i,"indexb"])
  tmp_indexa=as.character(statdf_nomerge[i,"indexa"])
  if (tmp_indexa == tmp_indexb) {
    statdf_nomerge[i,"total_number"] = statdf_nomerge[i,"diff_count"]
  } else {
    statdf_nomerge[i,"total_number"] = statdf_nomerge[statdf_nomerge$indexb==tmp_indexb & statdf_nomerge$indexa==tmp_indexa,"diff_count"]+
      statdf_nomerge[statdf_nomerge$indexa == tmp_indexb & statdf_nomerge$indexb == tmp_indexa,"diff_count"]
  }
}

rankname=rev(Used_cell2)
statdf_nomerge$indexa=factor(statdf_nomerge$indexa,levels = rankname)
statdf_nomerge$indexb=factor(statdf_nomerge$indexb,levels = rankname)
head(statdf_nomerge)

#获取相关性矩阵的上三角矩阵
statdf_dc<-dcast(statdf_nomerge[,c("indexb","indexa","total_number")],indexa~ indexb)
rownames(statdf_dc)<-statdf_dc$indexa;statdf_dc<-statdf_dc[,-1]
head(statdf_dc)
upper_tri <- get_upper_tri(statdf_dc)
upper_tri$indexa<-rownames(upper_tri)
statdf_merge <- melt(upper_tri,"indexa",na.rm = TRUE)
rownames(statdf_merge)<-paste0(statdf_merge$indexa,"-",statdf_merge$variable)
colnames(statdf_merge)<-c("indexa","indexb","diff_count")
head(statdf_merge);dim(statdf_merge)#136   3

#转成因子类型，画图时，图形将按照预先设置的顺序排列
statdf_merge$indexa=factor(statdf_merge$indexa,levels = rankname)
statdf_merge$indexb=factor(statdf_merge$indexb,levels = rankname)

range(statdf_merge$diff_count)#-114   22
ggplot(statdf_merge, aes(indexa, indexb)) + 
  geom_tile(aes(fill = diff_count),colour = "white") +
  geom_text(aes(label = round(diff_count, 1))) +
 # scale_fill_gradient2(low="navy",mid = "white",high="red",limits=c(-120,40),midpoint=0)+
  scale_fill_gradient2(low="cornflowerblue",mid = "gray100",high="red",limits=c(-120,40),midpoint=0)+
  #  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(-125,50))+
  scale_x_discrete("cluster 1")+scale_y_discrete("cluster 2")+
  theme_bw()+labs(title = "Abortion changed interaction number")+
  theme(
    axis.text.x.bottom  = element_text(hjust = 1, vjust = NULL, angle = 45),
    panel.grid = element_blank()
  )
ggsave(filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/changed interaction_number_merge_up_tri.pdf",device = "pdf",width = 22,height = 20,units = c("cm"))

#获取相关性矩阵的上三角矩阵
statdf_dc<-dcast(statdf_nomerge[,c("indexb","indexa","total_number")],indexa~ indexb)
rownames(statdf_dc)<-statdf_dc$indexa;statdf_dc<-statdf_dc[,-1]
head(statdf_dc)
lower_tri <- get_lower_tri(statdf_dc)
lower_tri$indexa<-rownames(lower_tri)
statdf_merge <- melt(lower_tri,"indexa",na.rm = TRUE)
rownames(statdf_merge)<-paste0(statdf_merge$indexa,"-",statdf_merge$variable)
colnames(statdf_merge)<-c("indexa","indexb","diff_count")
head(statdf_merge);dim(statdf_merge)#136   3

#转成因子类型，画图时，图形将按照预先设置的顺序排列
statdf_merge$indexa=factor(statdf_merge$indexa,levels = rankname)
statdf_merge$indexb=factor(statdf_merge$indexb,levels = rankname)
#不对称热图
range(statdf_merge$diff_count)#-114   22
ggplot(statdf_merge, aes(indexa, indexb)) + 
  geom_tile(aes(fill = diff_count),colour = "white") +
  geom_text(aes(label = round(diff_count, 1))) +
  # scale_fill_gradient2(low="navy",mid = "white",high="red",limits=c(-120,40),midpoint=0)+
  scale_fill_gradient2(low="steelblue3",mid = "gray100",high="red",limits=c(-120,40),midpoint=0)+
  #  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(-125,50))+
  scale_x_discrete("cluster 1")+scale_y_discrete("cluster 2")+
  theme_bw()+labs(title = "Abortion changed interaction number")+
  theme(
    axis.text.x.bottom  = element_text(hjust = 1, vjust = NULL, angle = 45),
    panel.grid = element_blank()
  )
ggsave(filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/changed interaction_number_merge_down_tri.pdf",device = "pdf",width = 22,height = 20,units = c("cm"))

#不对称热图
range(statdf_nomerge$diff_count)#-74    24
change_heatmap<-ggplot(statdf_nomerge,aes(x=indexa,y=indexb,fill=diff_count))+geom_tile(color="white")+
  geom_text(aes(label = round(diff_count, 1))) +
  scale_fill_gradient2(low="navy",mid = "white",high="red",limits=c(-75,25),midpoint=0)+
  scale_x_discrete("cluster 1 produces molecule 1")+
  scale_y_discrete("cluster 2 produces molecule 2")+
  theme_bw()+labs(title = "USM changed interaction number")+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
        panel.grid = element_blank())
change_heatmap
ggsave(change_heatmap,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/changed_interaction_number_no_merge.pdf",width = 16,height = 15)


#网络图表示合并后细胞互作关系的数量
mynet=statdf_merge[statdf_merge$diff_count != 0,] #过滤掉值为0的观测
colnames(mynet)<-c("indexa","indexb","value")
head(mynet);dim(mynet)#119   3
unique(as.character((mynet$indexa)))
rankname=rev(Used_cell2)[Used_cell2 %in%unique(as.character((mynet$indexa)))]
mynet$indexa=factor(mynet$indexa,levels = rankname)
rankname2=Used_cell2[Used_cell2 %in%unique(as.character((mynet$indexb)))]
mynet$indexb=factor(mynet$indexb,levels = rankname2)
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
rankname<-rankname2
color1=allcolour[1:length(rankname)]
names(color1)=rankname
color2=rev(colorRampPalette(brewer.pal(10,"RdBu"))(110))#将颜色分成多少份，取决于互作关系数目的最大值
show_col(color2)
names(color2)=c((-50:-1),1:60) #每一份颜色用对应的数字命名
range(E(net)$value[(E(net)$value < 0)])#-114   -1
range(E(net)$value[(E(net)$value > 0)])# 3    22

#设置网络性质
E(net)$width  <- abs(E(net)$value)/10  #根据count值设置边的宽 
E(net)$label = E(net)$value #连线的标注

E(net)$color <- color2[as.character(ifelse(E(net)$value <c(-50),-50,E(net)$value))]#用前面设置好的颜色赋给连线，颜色深浅对应数值大小
E(net)$label.color <- "black" #连线标注的颜色
V(net)$label.color <- "black" #节点标注的颜色
V(net)$color <- color1[names(V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息

#调整节点位置的线条角度
##如果没有这两行代码，节点位置的圆圈是向右的
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/subcell_changed_net_merge_plot.pdf",width=12,height=12)
plot(net,
     edge.arrow.size = 0, #连线不带箭头
     edge.curved = 0.2, #连线不弯曲
     vertex.frame.color = "black", #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = 30) #节点大小
title(main = " merged change count", sub = "Abortion versus CTRL")
dev.off()

#网络图表示非合并细胞互作关系的数量
mynet=statdf_nomerge[statdf_nomerge$diff_count != 0,c(1:3)] #过滤掉值为0的观测
colnames(mynet)<-c("indexa","indexb","value")
head(mynet);dim(mynet)#221   3
unique(as.character((mynet$indexa)))
rankname=rev(Used_cell2)[Used_cell2 %in%unique(as.character((mynet$indexa)))]
mynet$indexa=factor(mynet$indexa,levels = rankname)
rankname2=Used_cell2[Used_cell2 %in%unique(as.character((mynet$indexb)))]
mynet$indexb=factor(mynet$indexb,levels = rankname2)

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
rankname<-rankname2
color1=allcolour[1:length(rankname)]
names(color1)=rankname
range(E(net)$value[(E(net)$value < 0)])# -71  -1
range(E(net)$value[(E(net)$value > 0)])#  1 24
color2=rev(colorRampPalette(brewer.pal(10,"RdBu"))(110))#将颜色分成多少份，取决于互作关系数目的最大值
color2<-color2[c(1:50,70:100)]
show_col(color2)
names(color2)=c((-50:-1),1:30) #每一份颜色用对应的数字命名

#设置网络性质
E(net)$width  <- abs(E(net)$value)/5  #根据count值设置边的宽 
E(net)$label = E(net)$value #连线的标注
E(net)$color <- color2[as.character(ifelse(E(net)$value <c(-50),-50,E(net)$value))]#用前面设置好的颜色赋给连线，颜色深浅对应数值大小
E(net)$label.color <- "black" #连线标注的颜色
V(net)$label.color <- "black" #节点标注的颜色
V(net)$color <- color1[names(V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息

#调整节点位置的线条角度
##如果没有这两行代码，节点位置的圆圈是向右的
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/subcell_changed_net_not_merge_plot.pdf",width=12,height=12)
plot(net,
     edge.arrow.size = 0.1, #连线不带箭头
     edge.curved = 0.2, #连线不弯曲
     vertex.frame.color = "black", #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = 30) #节点大小
title(main = "no merged change count", sub = "Abortion versus CTRL")
dev.off()
