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
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)

merge_mean_pvalue4<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/cellpheno_merge_mean_pvalue.xls",sep = "\t", row.names=1, header =T) 
head(merge_mean_pvalue4);dim(merge_mean_pvalue4)
Cell_type<-unique(c(unlist(lapply(strsplit(as.character(merge_mean_pvalue4$Cell_pair),"-"), function(x) x[1])),unlist(lapply(strsplit(as.character(merge_mean_pvalue4$Cell_pair),"-"), function(x) x[2]))))
Cell_type2<-as.data.frame(t(combn(Cell_type,2) ))
Cell_type_Dir_a<-c(paste(Cell_type2$V1,Cell_type2$V2,sep="-"))
Cell_type_Dir_b<-c(paste(Cell_type2$V2,Cell_type2$V1,sep="-"))
Cell_type_Dir_c<-c(paste(Cell_type,Cell_type,sep="-"))
Cell_type_Dir_a[order(Cell_type_Dir_a)]


merge_mean_pvalue4[1:6,]
merge_mean_pvalue5<-na.omit(merge_mean_pvalue4)
merge_mean_pvalue_CTRL<-na.omit(merge_mean_pvalue4[which(merge_mean_pvalue4$Treat_group =="CTRL" ),])
merge_mean_pvalue_Abortion<-na.omit(merge_mean_pvalue4[which(merge_mean_pvalue4$Treat_group =="Abortion" ),])

Cell_pair<-c(Cell_type_Dir_a,Cell_type_Dir_b,Cell_type_Dir_c)
sp_CPI_dataframe<-data.frame()
for (sp_Cell_pair in Cell_pair){
  #sp_Cell_pair<-"CTBs_1-CTBs_2"
  print(sp_Cell_pair)
  CPI_CTRL<-unique(merge_mean_pvalue_CTRL[which(merge_mean_pvalue_CTRL$Cell_pair ==sp_Cell_pair),]$interacting_pair)
  CPI_Abortion<-unique(merge_mean_pvalue_Abortion[which(merge_mean_pvalue_Abortion$Cell_pair ==sp_Cell_pair),]$interacting_pair)
  CPI_CTRL_sp<-setdiff(CPI_CTRL,CPI_Abortion)
  CPI_Abortion_sp<-setdiff(CPI_Abortion,CPI_CTRL)
  merge_mean_pvalue6<-merge_mean_pvalue5[which((merge_mean_pvalue5$Cell_pair == sp_Cell_pair) & (merge_mean_pvalue5$interacting_pair %in% c(CPI_CTRL_sp,CPI_Abortion_sp))),]
  sp_CPI_dataframe<-rbind(sp_CPI_dataframe,merge_mean_pvalue6)
}
head(sp_CPI_dataframe)
write.table(as.data.frame(sp_CPI_dataframe), file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Changed_CPI_details_in_each_selected_celltype.xls",sep = "\t", row.names=T, col.names=T) 

data_sum<-as.data.frame(table(sp_CPI_dataframe$interacting_pair,sp_CPI_dataframe$Treat_group))
merge_mean_pvalue5[1:6,]

#top common CPI
common_freq<-data.frame(table(sp_CPI_dataframe$interacting_pair,sp_CPI_dataframe$Treat_group))
common_freq2 <-dcast(common_freq, Var1  ~ Var2  , value.var = "Freq")
common_freq2$add_freq<-common_freq2$Abortion + common_freq2$CTRL
head(common_freq2);dim(common_freq2)#613   4
write.table(as.data.frame(common_freq2), file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Changed_CPI_number_summary.xls",sep = "\t", row.names=T, col.names=T) 
#common_freq2<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Changed_CPI_number_summary.xls",sep = "\t", header=T) 
head(common_freq2)
#distinct(DEGs_merge_data[,c("gene","DEGs_inter_type")])
common_freq2_CTRL<-common_freq2[which(common_freq2$CTRL>0),]
common_freq2_Abortion<-common_freq2[which(common_freq2$Abortion>0),]
common_freq2_CTRL<-common_freq2_CTRL[order(common_freq2_CTRL$CTRL,decreasing = T),]
common_freq2_Abortion<-common_freq2_Abortion[order(common_freq2_Abortion$Abortion,decreasing = T),]
head(common_freq2_CTRL);head(common_freq2_Abortion)
data_CTRL<-head(common_freq2_CTRL,n=10)
data_CTRL$Var1<-factor(data_CTRL$Var1,levels = as.character(data_CTRL$Var1))
data_Abortion<-head(common_freq2_Abortion,n=10)
data_Abortion$Var1<-factor(data_Abortion$Var1,levels = as.character(data_Abortion$Var1))
plot_Abortion<-ggplot(data_Abortion, aes(x = Var1, y =  Abortion)) +  geom_bar(stat = "identity", alpha = 0.9,fill=ppCor[1]) + coord_polar() +
  geom_text(aes(x= Var1, y= Abortion+2, label =  Abortion),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of cell_pair for top 10 common Abortion_special_CPI") + 
  theme(title = element_text(vjust = -56, face = "bold"))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
plot_CTRL<-ggplot(data_CTRL, aes(x = Var1, y =  CTRL)) +  geom_bar(stat = "identity", alpha = 0.9,fill=ppCor[2]) + coord_polar() +
  geom_text(aes(x= Var1, y= CTRL+1, label =  CTRL),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of cell_pair for top 10 common CTRL_special_CPI") + 
  theme(title = element_text(vjust = -56, face = "bold"))+
  theme(axis.text.y = element_blank(),panel.grid.minor.y  = element_blank(),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1),
        axis.ticks = element_blank(), axis.title = element_blank())
plot_Abortion
plot_CTRL
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Changed_CPI_Top10_common_number_for_cell_pair.pdf",width =12,height = 6)
grid.arrange(plot_Abortion, plot_CTRL,ncol=2)
dev.off()

##summary for cell_paires with most changed CPI
## Unique CPI number for each Cells
length(unique(merge_mean_pvalue5$Cell_pair))#225
num_cell_plot<-merge_mean_pvalue5[,c("Cell_pair","Treat_group")]
num_cell_plot$Num<-1
num_cell_plot2<-aggregate(num_cell_plot$Num, by=list(num_cell_plot$Cell_pair,num_cell_plot$Treat_group), FUN=sum)
#table(D_ac_TFs_active_high_matrix$cell_identity,D_ac_TFs_active_high_matrix$change)
colnames(num_cell_plot2)<-c("Cell_pair","tendency","TFs_num")
num_cell_plot2$cell_one<- as.character(unlist(lapply(strsplit(as.character(num_cell_plot2$Cell_pair),"-"), function(x) x[1])))
num_cell_plot2$cell_two<- as.character(unlist(lapply(strsplit(as.character(num_cell_plot2$Cell_pair),"-"), function(x) x[2])))
head(num_cell_plot2)
length(unique(c(num_cell_plot2$cell_one,num_cell_plot2$cell_two)))#15

num_cell_plot3 <-dcast(num_cell_plot2, Cell_pair ~ tendency , value.var = "TFs_num")
head(num_cell_plot3)
num_cell_plot3[is.na(num_cell_plot3)]<-0
num_cell_plot3$cell_one<- as.character(unlist(lapply(strsplit(as.character(num_cell_plot3$Cell_pair),"-"), function(x) x[1])))
num_cell_plot3$cell_two<- as.character(unlist(lapply(strsplit(as.character(num_cell_plot3$Cell_pair),"-"), function(x) x[2])))
head(num_cell_plot3)

head(num_cell_plot2)
num_cell_plot4<-merge(num_cell_plot2,num_cell_plot3[, c("Cell_pair","Abortion","CTRL")],by="Cell_pair")
num_cell_plot4$total_num<-num_cell_plot4$Abortion+num_cell_plot4$CTRL
head(num_cell_plot4)
Used_cell<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Mycs","HCs","STCs","FBs","Endo","NKs","Ts")
num_cell_plot4$cell_one <-factor(num_cell_plot4$cell_one,levels=Used_cell)
num_cell_plot4$tendency<-factor(num_cell_plot4$tendency,levels=c("Abortion","CTRL"))

TF_roseplot3<-ggplot(num_cell_plot4, aes(x = Cell_pair, y = TFs_num, fill = tendency)) +  
  geom_bar(stat = "identity", alpha = 0.5) + coord_polar() +
  theme_grey() + labs(x = "", y = "", title = "The number of changed CPI for each cellpair") + 
  geom_text(aes(x= Cell_pair, y= 400, label = Abortion),data =num_cell_plot4,hjust =0.5,vjust =0.5,color=ppCor[1], size=2 )+
  geom_text(aes(x= Cell_pair, y= 380, label = CTRL),data =num_cell_plot4,hjust =0.5,vjust =0.5,color=ppCor[2], size=2 )+
  # theme(panel.border = element_blank()) +   ## 去掉外层边框
  theme(legend.position = "right") +   ## 去掉图例
  theme(title = element_text(vjust = -56, face = "bold"))+
  scale_fill_manual(name="Change class",values=c(ppCor[1],ppCor[2]))
TF_roseplot4 <- TF_roseplot3+ theme(axis.text.x = element_blank()) 
TF_roseplot3
TF_roseplot4
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Changed_CPI_number_details_in_each_cellpair.pdf",TF_roseplot3,width=12, height=12)
ggsave("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Changed_CPI_number_details_in_each_cellpair2.pdf",TF_roseplot4,width=12, height=12)
##plot for top change cell pairs
head(num_cell_plot4)
unique(num_cell_plot4[which(num_cell_plot4$total_num==max(num_cell_plot4$total_num)),]$Cell_pair) 
unique(num_cell_plot4[which(num_cell_plot4$Abortion==max(num_cell_plot4$Abortion)),]$Cell_pair) 
unique(num_cell_plot4[which(num_cell_plot4$CTRL ==max(num_cell_plot4$CTRL )),]$Cell_pair) 

num_cell_plot4[order(num_cell_plot4$total_num,decreasing = T),]

########绘制目标细胞对中具体变化的分子对

##绘制互作频率变化最大的细胞对的具体变化分子对
CPI_down_most<-"FBs-EVTs_3"##-71
table(merge_mean_pvalue5[which(merge_mean_pvalue5$Cell_pair ==CPI_down_most),]$Treat_group)#Abortion     CTRL  66      137 
table(sp_CPI_dataframe[which(sp_CPI_dataframe$Cell_pair ==CPI_down_most),]$Treat_group)    #Abortion     CTRL  7       78 
CPI_up_most<-"FBs-EVTs_1"##24
table(merge_mean_pvalue5[which(merge_mean_pvalue5$Cell_pair ==CPI_up_most),]$Treat_group) #Abortion     CTRL  117       93 
table(sp_CPI_dataframe[which(sp_CPI_dataframe$Cell_pair ==CPI_up_most),]$Treat_group)     #Abortion     CTRL  49       25 

##for plot CPI
##down
data_down_most<-sp_CPI_dataframe[which(sp_CPI_dataframe$Cell_pair ==CPI_down_most),]
data_up_most<-sp_CPI_dataframe[which(sp_CPI_dataframe$Cell_pair ==CPI_up_most),]
head(data_down_most)

library(tidyverse)
library(tidygraph)
library(ggraph)
library(ggdark)
library(ggpubr)
##netplot
CPI_pair<-as.character(data_down_most$interacting_pair)[grep("complex",as.character(data_down_most$interacting_pair),invert =T)]
data_down_most2<-data_down_most[grep("complex",as.character(data_down_most$interacting_pair),invert =T),]
CPI_pair<-gsub("TGFbeta receptor", "TGFBR", CPI_pair) 
CPI_pair<-gsub("ACVR_", "ACVR-", CPI_pair) 
down_most_CPI_pair<-gsub(" receptor", "-receptor", CPI_pair) 
CPI_pair<-as.character(data_up_most$interacting_pair)[grep("complex",as.character(data_up_most$interacting_pair),invert =T)]
data_up_most2<-data_up_most[grep("complex",as.character(data_up_most$interacting_pair),invert =T),]
CPI_pair<-gsub("TGFbeta receptor", "TGFBR", CPI_pair) 
CPI_pair<-gsub("ACVR_", "ACVR-", CPI_pair) 
CPI_pair<-gsub(" receptor", "-receptor", CPI_pair) 
Up_most_CPI_pair<-gsub("_TGFR_AVR2A", "_TGFR-AVR2A", CPI_pair) 

down_send<-unlist(lapply(strsplit(down_most_CPI_pair,"_"), function(x) x[1]))
down_receptor<-unlist(lapply(strsplit(down_most_CPI_pair,"_"), function(x) x[2]))
UP_send<-unlist(lapply(strsplit(Up_most_CPI_pair,"_"), function(x) x[1]))
UP_receptor<-unlist(lapply(strsplit(Up_most_CPI_pair,"_"), function(x) x[2]))

##for down 
edges<-data.frame(from=down_send,to=down_receptor,edge.colour =data_down_most2$Treat_group,edge.width = data_down_most2$mean)
node <- unique(c(edges$from, edges$to)) %>% sort()
nodes <- data.frame(node, node.size = as.numeric(table(c(down_send,down_receptor))),stringsAsFactors = FALSE)
graph_data <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
graph_data

##我们实际使用的是这个数据框来绘制的，因此，我们可以在 aes() 中使用布局算法返回的数据框中的变量
plot<-ggraph(graph = graph_data, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width),edge_alpha=0.7,arrow = arrow(length = unit(4, 'mm')), start_cap = circle(3, 'mm'),end_cap = circle(3, 'mm')) +
  scale_colour_manual(values = c("Abortion" = ppCor[3], "CTRL" = ppCor[2])) +
  #geom_edge_fan(aes(edge_colour = edge.colour,edge_width =edge.width)) +
  geom_node_point(aes(size = node.size),colour="grey") +
  scale_size_continuous(range = c(5,10)) +
  guides(colour = guide_legend(order = 1),size = guide_legend(order = 2),colour = guide_edge_colourbar(order = 3))+
 # guides(edge_alpha = guide_edge_direction())+
  #scale_edge_colour_gradient2(low = "#155F83FF", mid = "white", high = "#800000FF") +
  theme_void()
  #ggdark::dark_theme_void()
node_name<-node
angle <- 360 * (c(1:length(node_name)) - 0.5)/length(node_name)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)
plot2<-plot + geom_node_text(aes(x = x * 1.05,y = y * 1.15,label = node_name),angle = angle, hjust = hjust,colour = "black",size = 3)+ 
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
plot_down<-plot2+labs(title=paste0("Cell Communiction CPI pairs in most down  changed: ",CPI_down_most))
ggsave(plot_down,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/circle_plot_for_most_down_change_CPI_pairs.pdf",width = 11,height = 10)

##for UP 
edges<-data.frame(from=UP_send,to=UP_receptor,edge.colour =data_up_most2$Treat_group,edge.width = data_up_most2$mean)
node <- unique(c(edges$from, edges$to)) %>% sort()
nodes <- data.frame(node, node.size = as.numeric(table(c(UP_send,UP_receptor))),stringsAsFactors = FALSE)
graph_data <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
graph_data

##我们实际使用的是这个数据框来绘制的，因此，我们可以在 aes() 中使用布局算法返回的数据框中的变量
plot<-ggraph(graph = graph_data, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width),edge_alpha=0.7,arrow = arrow(length = unit(4, 'mm')), start_cap = circle(3, 'mm'),end_cap = circle(3, 'mm')) +
  scale_colour_manual(values = c("Abortion" = ppCor[3], "CTRL" = ppCor[2])) +
  #geom_edge_fan(aes(edge_colour = edge.colour,edge_width =edge.width)) +
  geom_node_point(aes(size = node.size),colour="grey") +
  scale_size_continuous(range = c(5,10)) +
  guides(colour = guide_legend(order = 1),size = guide_legend(order = 2),colour = guide_edge_colourbar(order = 3))+
  # guides(edge_alpha = guide_edge_direction())+
  #scale_edge_colour_gradient2(low = "#155F83FF", mid = "white", high = "#800000FF") +
  theme_void()
#ggdark::dark_theme_void()
node_name<-node
angle <- 360 * (c(1:length(node_name)) - 0.5)/length(node_name)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)
plot2<-plot + geom_node_text(aes(x = x * 1.05,y = y * 1.15,label = node_name),angle = angle, hjust = hjust,colour = "black",size = 3)+ 
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
plot_UP<-plot2+labs(title=paste0("Cell Communiction CPI pairs in most UP  changed: ",CPI_up_most))
ggsave(plot_UP,filename = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/circle_plot_for_most_up_change_CPI_pairs.pdf",width = 11,height = 10)

#ggarrange(plot_UP,plot_down,widths = c(2,1), labels = c('a', 'b'))


##其他变形待尝试https://min.news/zh-cn/technique/5812887eae2456997cba2b02f3abe9ae.html
##https://www.jianshu.com/p/059b6a9402b8
##其他美化http://www.17bigdata.com/r%E8%AF%AD%E8%A8%80%E5%88%A9%E7%94%A8igraph%E5%92%8Cnetworkd3%E5%8C%85%E5%BF%AB%E9%80%9F%E5%85%A5%E9%97%A8%E5%81%9A%E5%87%BA%E7%82%AB%E9%85%B7%E7%9A%84%E7%A4%BE%E4%BA%A4%E7%BD%91%E7%BB%9C%E5%9B%BE/
networkD3