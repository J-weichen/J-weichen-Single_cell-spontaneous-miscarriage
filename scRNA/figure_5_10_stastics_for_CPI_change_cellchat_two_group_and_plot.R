rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

##Load the required libraries
library(reticulate)
use_python('/mnt/data/chenwei/software/miniconda2_new/envs/scenic/bin/',required = T)
library(CellChat)
library(patchwork)
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
library(ggsci)
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

#Create a directory to save figures
data.dir <- '/mnt/data/chenwei/gongchen/7.cellchat/Merge'
#dir.create(data.dir)
setwd(data.dir)

#Load CellChat object of each dataset and then merge together
cellchat.CTRL <- readRDS("/mnt/data/chenwei/gongchen/7.cellchat/CTRL/ctrl_cellchat.rds")
cellchat.Abortion <- readRDS( "/mnt/data/chenwei/gongchen/7.cellchat/Abortion/Abortion_cellchat.rds")
object.list <- list(CTRL = cellchat.CTRL, Abortion = cellchat.Abortion)
##提取配体和受体互作的通讯表格
##note:结果存储于@net下,增加一列概率值和对应的pvalue
##cellchat 中的computecommunprob函数 根据表达值推测细胞互作的概率；这一步也是cellchat比cellphenoDB多出来的一步
##cellphenoDB 仅仅用互作分子的平均表达值作为互作强度
##信号通路水平的通信概率： 进一步计算与每一个信号通路相关的所有配体-受体相互作用的通信概率

#step one:extract interaction information in each group for each interpairs
CTRL_net <- subsetCommunication(cellchat.CTRL)
USM_net <- subsetCommunication(cellchat.Abortion)
write.table(CTRL_net,file="/mnt/data/chenwei/gongchen/7.cellchat/Merge/CTRL_interaction_communprob.txt",sep="\t", quote=F, row.names=T,col.names=T)
write.table(USM_net,file="/mnt/data/chenwei/gongchen/7.cellchat/Merge/USM_interaction_communprob.txt",sep="\t", quote=F, row.names=T,col.names=T)
CTRL_net$class<-"CTRL";USM_net$class<-"USM"
all_nat<-rbind(CTRL_net,USM_net)
write.table(all_nat,file="/mnt/data/chenwei/gongchen/7.cellchat/Merge/All_interaction_communprob.txt",sep="\t", quote=F, row.names=T,col.names=T)

#step two:extract interaction information in each group for each pathway
CTRL_netpath <- subsetCommunication(cellchat.CTRL,slot.name = "netP")
USM_netpath  <- subsetCommunication(cellchat.Abortion,slot.name = "netP")
write.table(CTRL_netpath,file="/mnt/data/chenwei/gongchen/7.cellchat/Merge/CTRL_inter_pathway_communprob.txt",sep="\t", quote=F, row.names=T,col.names=T)
write.table(USM_netpath,file="/mnt/data/chenwei/gongchen/7.cellchat/Merge/USM_inter_pathway_communprob.txt",sep="\t", quote=F, row.names=T,col.names=T)
CTRL_netpath$class<-"CTRL";USM_netpath$class<-"USM"
all_netpath<-rbind(CTRL_netpath,USM_netpath)
write.table(all_netpath,file="/mnt/data/chenwei/gongchen/7.cellchat/Merge/All_inter_pathway_communprob.txt",sep="\t", quote=F, row.names=T,col.names=T)


#step three:calculate unique interaction information in each cell pairs
all_nat<-read.table(file="/mnt/data/chenwei/gongchen/7.cellchat/Merge/All_interaction_communprob.txt",sep="\t", row.names=1,header =T)
head(all_nat)

all_nat$cell_pairs<-paste(all_nat$source,all_nat$target,sep="-")
sp_CPI_dataframe<-data.frame()
for (sp_Cell_pair in unique(all_nat$cell_pairs)){
  #sp_Cell_pair<- "STBs_1-STBs_1"
  print(sp_Cell_pair)
  CPI_CTRL<-unique(all_nat[which(all_nat$cell_pairs ==sp_Cell_pair & all_nat$class == "CTRL"),]$interaction_name)
  CPI_Abortion<-unique(all_nat[which(all_nat$cell_pairs ==sp_Cell_pair & all_nat$class == "USM"),]$interaction_name)
  CPI_CTRL_sp<-setdiff(CPI_CTRL,CPI_Abortion)
  CPI_Abortion_sp<-setdiff(CPI_Abortion,CPI_CTRL)
  Selected_net<-all_nat[which(all_nat$cell_pairs == sp_Cell_pair & all_nat$interaction_name %in% c(CPI_CTRL_sp,CPI_Abortion_sp)),]
  sp_CPI_dataframe<-rbind(sp_CPI_dataframe,Selected_net)
}
head(sp_CPI_dataframe)
write.table(as.data.frame(sp_CPI_dataframe), file="/mnt/data/chenwei/gongchen/7.cellchat/Merge/Changed_CPI_details_in_all_cellpairs.txt",sep = "\t", row.names=T, col.names=T) 

EVT3_interaction<-sp_CPI_dataframe[which(sp_CPI_dataframe$source =="EVTs_3"| sp_CPI_dataframe$target =="EVTs_3"),]
table(EVT3_interaction$target)

data_sum<-as.data.frame(table(sp_CPI_dataframe$interaction_name,sp_CPI_dataframe$class))

#top common CPI
common_freq<-data.frame(table(sp_CPI_dataframe$interaction_name,sp_CPI_dataframe$class))
common_freq2 <-dcast(common_freq, Var1  ~ Var2  , value.var = "Freq")
common_freq2$add_freq<-common_freq2$USM  + common_freq2$CTRL
head(common_freq2);dim(common_freq2)#613   4
write.table(as.data.frame(common_freq2), file="/mnt/data/chenwei/gongchen/7.cellchat/Merge/Changed_CPI_number_summary_cellchat.txt",sep = "\t", row.names=T, col.names=T) 
#common_freq2<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_4/subcelltype/Changed_CPI_number_summary.xls",sep = "\t", header=T) 
head(common_freq2)
#distinct(DEGs_merge_data[,c("gene","DEGs_inter_type")])
common_freq2_CTRL<-common_freq2[which(common_freq2$CTRL>0),]
common_freq2_USM<-common_freq2[which(common_freq2$USM >0),]
common_freq2_CTRL<-common_freq2_CTRL[order(common_freq2_CTRL$CTRL,decreasing = T),]
common_freq2_USM<-common_freq2_USM[order(common_freq2_USM$USM,decreasing = T),]
head(common_freq2_CTRL);head(common_freq2_USM)
data_CTRL<-head(common_freq2_CTRL,n=10)
data_CTRL$Var1<-factor(data_CTRL$Var1,levels = as.character(data_CTRL$Var1))
data_USM<-head(common_freq2_USM,n=10)
data_USM$Var1<-factor(data_USM$Var1,levels = as.character(data_USM$Var1))
plot_USM<-ggplot(data_USM, aes(x = Var1, y =  USM)) +  geom_bar(stat = "identity", alpha = 0.9,fill=ppCor[1]) + coord_polar() +
  geom_text(aes(x= Var1, y= USM+2, label =  USM),hjust =0.5,vjust =0.5,color="black", size=3 )+
  theme_grey() + labs(x = "", y = "", title = "The number of cell_pair for top 10 common USM_special_CPI") + 
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
plot_USM
plot_CTRL
pdf("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Changed_CPI_Top10_common_number_for_cell_pair.pdf",width =12,height = 6)
grid.arrange(plot_USM, plot_CTRL,ncol=2)
dev.off()

##summary for cell_paires with most changed CPI
## Unique CPI number for each Cells
length(unique(sp_CPI_dataframe$cell_pairs))#225
num_cell_plot<-sp_CPI_dataframe[,c("cell_pairs","class")]
num_cell_plot$Num<-1
num_cell_plot2<-aggregate(num_cell_plot$Num, by=list(num_cell_plot$cell_pairs,num_cell_plot$class), FUN=sum)
#table(D_ac_TFs_active_high_matrix$cell_identity,D_ac_TFs_active_high_matrix$change)
colnames(num_cell_plot2)<-c("Cell_pair","tendency","CPI_num")
num_cell_plot2$cell_one<- as.character(unlist(lapply(strsplit(as.character(num_cell_plot2$Cell_pair),"-"), function(x) x[1])))
num_cell_plot2$cell_two<- as.character(unlist(lapply(strsplit(as.character(num_cell_plot2$Cell_pair),"-"), function(x) x[2])))
head(num_cell_plot2)
length(unique(c(num_cell_plot2$cell_one,num_cell_plot2$cell_two)))#15

num_cell_plot3 <-dcast(num_cell_plot2, Cell_pair ~ tendency , value.var = "CPI_num")
head(num_cell_plot3)
num_cell_plot3[is.na(num_cell_plot3)]<-0
num_cell_plot3$cell_one<- as.character(unlist(lapply(strsplit(as.character(num_cell_plot3$Cell_pair),"-"), function(x) x[1])))
num_cell_plot3$cell_two<- as.character(unlist(lapply(strsplit(as.character(num_cell_plot3$Cell_pair),"-"), function(x) x[2])))
head(num_cell_plot3)

head(num_cell_plot2)
num_cell_plot4<-merge(num_cell_plot2,num_cell_plot3[, c("Cell_pair","USM","CTRL")],by="Cell_pair")
num_cell_plot4$total_num<-num_cell_plot4$USM+num_cell_plot4$CTRL
head(num_cell_plot4)
Used_cell<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Mycs","HCs","STCs","FBs","Endo","NKs","Ts")
num_cell_plot4$cell_one <-factor(num_cell_plot4$cell_one,levels=Used_cell)
num_cell_plot4$tendency<-factor(num_cell_plot4$tendency,levels=c("USM","CTRL"))

range(num_cell_plot4$total_num)
TF_roseplot3<-ggplot(num_cell_plot4, aes(x = Cell_pair, y = CPI_num, fill = tendency)) +  
  geom_bar(stat = "identity", alpha = 0.5) + coord_polar() +
  theme_grey() + labs(x = "", y = "", title = "The number of changed CPI for each cellpair") + 
  geom_text(aes(x= Cell_pair, y= 100, label = USM),data =num_cell_plot4,hjust =0.5,vjust =0.5,color=ppCor[1], size=2 )+
  geom_text(aes(x= Cell_pair, y= 95, label = CTRL),data =num_cell_plot4,hjust =0.5,vjust =0.5,color=ppCor[2], size=2 )+
  # theme(panel.border = element_blank()) +   ## 去掉外层边框
  theme(legend.position = "right") +   ## 去掉图例
  theme(title = element_text(vjust = -56, face = "bold"))+
  scale_fill_manual(name="Change class",values=c(ppCor[1],ppCor[2]))
TF_roseplot4 <- TF_roseplot3+ theme(axis.text.x = element_blank()) 
TF_roseplot4

ggsave("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Changed_CPI_number_details_in_each_cellpair_cellchat.pdf",TF_roseplot3,width=12, height=12)
ggsave("/mnt/data/chenwei/gongchen/7.cellchat/Merge/Changed_CPI_number_details_in_each_cellpair_cellchat2.pdf",TF_roseplot4,width=12, height=12)
##plot for top change cell pairs
head(num_cell_plot4)
unique(num_cell_plot4[which(num_cell_plot4$total_num==max(num_cell_plot4$total_num)),]$Cell_pair) 
unique(num_cell_plot4[which(num_cell_plot4$Abortion==max(num_cell_plot4$Abortion)),]$Cell_pair) 
unique(num_cell_plot4[which(num_cell_plot4$CTRL ==max(num_cell_plot4$CTRL )),]$Cell_pair) 

num_cell_plot4[order(num_cell_plot4$total_num,decreasing = T),]
