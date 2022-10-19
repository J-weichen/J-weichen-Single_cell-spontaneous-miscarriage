#For default 
files<- list.files("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup/BP_GO")
Enrich_all_plot<-data_frame()
for ( f in files[1:length(files)]){
  # f=files[1]
  samplenames <- gsub("_q005_unselect_BP_GO.txt","",f)
  print(samplenames)
  FF<-paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup/BP_GO/",f,sep="")
  column<-read.table(FF,head = T)
  head(column)
  column2<-column[which(column$Count>=2),]
  if (is.na(column2$Description)) {
    next
  }
  column2$group<-rep(samplenames,nrow(column2))
  Enrich_all_plot<- rbind(Enrich_all_plot,column2)
}
str(Enrich_all_plot);
Enrich_all_plot2<-Enrich_all_plot[which(!(Enrich_all_plot$group %in% c("Ery_down","Ery_up"))),]
dim(Enrich_all_plot);dim(Enrich_all_plot2)#64190   58732
table(Enrich_all_plot2$group)
length(unique(Enrich_all_plot2$Description))#5882
length(unique(Enrich_all_plot2[which(Enrich_all_plot2$p.adjust<0.05),]$Description))#3127
length(unique(Enrich_all_plot2[which(Enrich_all_plot2$p.adjust<0.01),]$Description))#1710

Enrich_all_plot2$LogPadjust<- c(-log10(Enrich_all_plot2$p.adjust))
range(Enrich_all_plot2$LogPadjust)#2.017767e-10 5.480113e+01
range(Enrich_all_plot2$Count)#2 81

Enrich_all_number <- data.frame(table(as.character(Enrich_all_plot2$Description)))
Enrich_all_number2 <- Enrich_all_number[which(Enrich_all_number$Freq>1),]
dim(Enrich_all_number2);dim(Enrich_all_number)#4895 5882 

colnames(Enrich_all_number)<-c("Description","freq")
colnames(Enrich_all_number2)<-c("Description","freq")

Enrich_all_plot3<-merge(Enrich_all_plot2,Enrich_all_number)
tail(Enrich_all_plot3);dim(Enrich_all_plot3)#58732    12
write.table(Enrich_all_plot3, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup_all_GO_BP_merge.txt",row.names=T, col.names=T) 
Enrich_all_plot3<-merge(Enrich_all_plot2,Enrich_all_number2)
tail(Enrich_all_plot3);dim(Enrich_all_plot3)#57745    12
write.table(Enrich_all_plot3, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup_common_GO_BP_merge.txt",row.names=T, col.names=T) 

#for no selected GO terms of change protein
Enrich_all_plot_noselect<-read.table("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup_all_GO_BP_merge.txt",header =T) 
head(Enrich_all_plot_noselect)
length(unique(Enrich_all_plot_noselect$Description))#5882 
length(unique(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>2),]$Description))#4063 
Enrich_all_plot_noselect$trend <- unlist(lapply(strsplit(as.character(Enrich_all_plot_noselect$group), "_"), function(x) x[length(x)]))

names(table(Enrich_all_plot_noselect$group))
group_order<-c("CTBs_0_down","CTBs_1_down","CTBs_2_down","CTBs_3_down","STBs_1_down","STBs_2_down",
               "EVTs_1_down","EVTs_2_down","EVTs_3_down","STCs_down","FBs_down",
               "Endo_down","Epi_down", "HCs_down","MyCs_down","NKs_down","Ts_down","Bs_down",
               "CTBs_0_up","CTBs_1_up","CTBs_2_up","CTBs_3_up","STBs_1_up","STBs_2_up",
               "EVTs_1_up","EVTs_2_up","EVTs_3_up","STCs_up","FBs_up","Endo_up","HCs_up","MyCs_up","NKs_up","Ts_up")

head(Enrich_all_plot_noselect)
Enrich_all_plot_noselect$group<-factor(Enrich_all_plot_noselect$group,levels = group_order)
Enrich_all_plot_noselect <- Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$group),]

Enrich_all_plot_noselect<-Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$freq,decreasing = F),]

Enrich_all_plot_noselect$Description<-factor(Enrich_all_plot_noselect$Description,levels = c(unique(as.character(Enrich_all_plot_noselect$Description))))
head(Enrich_all_plot_noselect)
range(Enrich_all_plot_noselect$LogPadjust) #2.017767e-10 5.480113e+01
range(Enrich_all_plot_noselect$Count) #2 81
library(RColorBrewer)
library(stringr)

#ggplot(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>40),], 
#ggplot(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$freq>20),], 
ggplot(Enrich_all_plot_noselect,
  aes(x=group,y=Description,size=Count,colour=LogPadjust))+
  geom_point(alpha =0.8,na.rm = TRUE)+
  scale_size(breaks = c(5,10,15,20,90),range = c(1,3),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,2,4,6,8),name='-log10(p_value_adjust)')+  
 # scale_y_discrete(labels=function(x) str_wrap(x, width=150))+
  theme_classic()+labs(x="",y="GO terms",title="AMA related change protein BP enrichment")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 2,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1,angle = 90),
        legend.title = element_text(size = 10),
        legend.position ="right",legend.direction = "vertical")
write.table(Enrich_all_plot_noselect, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/wait_selected_major_subgroup_all_GO_BP_merge.txt",row.names=T, col.names=T) 

#For default 
UP_GO_BP_files<- "/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/merge_DEG_enrichment/Up_uniq_q005_unselect_BP_GO.txt"
up_column<-read.table(UP_GO_BP_files,head = T)
up_column2<-up_column[which(up_column$Count>=2),]
up_column2$group<-rep("Up_uniq",nrow(up_column2))
head(up_column2)
Down_GO_BP_files<- "/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/merge_DEG_enrichment/Down_uniq_q005_unselect_BP_GO.txt"
Down_column<-read.table(Down_GO_BP_files,head = T)
Down_column2<-Down_column[which(Down_column$Count>=2),]
Down_column2$group<-rep("Down_uniq",nrow(Down_column2))
head(Down_column2)
Enrich_all_plot<- rbind(up_column2,Down_column2)

str(Enrich_all_plot);dim(Enrich_all_plot)
table(Enrich_all_plot$group)
#Down_uniq   Up_uniq 
# 4374      5166 
length(unique(Enrich_all_plot$Description))#6181
length(unique(Enrich_all_plot[which(Enrich_all_plot$p.adjust<0.05),]$Description))# 1827
length(unique(Enrich_all_plot[which(Enrich_all_plot$p.adjust<0.01),]$Description))#1050

Enrich_all_plot$LogPadjust<- c(-log10(Enrich_all_plot$p.adjust))
range(Enrich_all_plot$LogPadjust)#8.987532e-14 3.984749e+01
range(Enrich_all_plot$Count)#2 113

Enrich_all_number <- data.frame(table(as.character(Enrich_all_plot$Description)))
Enrich_all_number2 <- Enrich_all_number[which(Enrich_all_number$Freq>1),]
dim(Enrich_all_number2);dim(Enrich_all_number)#3359   6181 

colnames(Enrich_all_number)<-c("Description","freq")
colnames(Enrich_all_number2)<-c("Description","freq")

Enrich_all_plot2<-merge(Enrich_all_plot,Enrich_all_number)
tail(Enrich_all_plot2);dim(Enrich_all_plot2)#9540   12
write.table(Enrich_all_plot2, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/UP_Down_unique_DEGs_all_GO_BP_merge.txt",row.names=T, col.names=T) 
Enrich_all_plot3<-merge(Enrich_all_plot,Enrich_all_number2)
tail(Enrich_all_plot3);dim(Enrich_all_plot3)# 6718   12
write.table(Enrich_all_plot3, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/UP_Down_unique_DEGs_common_GO_BP_merge.txt",row.names=T, col.names=T) 

#for no selected GO terms of change protein
Enrich_all_plot_noselect<-read.table("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/UP_Down_unique_DEGs_all_GO_BP_merge.txt",header =T) 
head(Enrich_all_plot_noselect)
length(unique(Enrich_all_plot_noselect$Description))#6181 
length(unique(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>2),]$Description))# 4659 

head(Enrich_all_plot_noselect)
table(Enrich_all_plot_noselect$group)
Enrich_all_plot_noselect$group<-factor(Enrich_all_plot_noselect$group,levels = c("Down_uniq","Up_uniq"))
Enrich_all_plot_noselect <- Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$group),]
Enrich_all_plot_noselect<-Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$freq,decreasing = F),]
Enrich_all_plot_noselect$Description<-factor(Enrich_all_plot_noselect$Description,levels = c(unique(as.character(Enrich_all_plot_noselect$Description))))
head(Enrich_all_plot_noselect)
range(Enrich_all_plot_noselect$LogPadjust) #8.987532e-14 3.984749e+01
range(Enrich_all_plot_noselect$Count) #2 81
library(RColorBrewer)
library(stringr)

#ggplot(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>40),], 
ggplot(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$p.adjust<0.05),], 
       #ggplot(Enrich_all_plot_noselect,
       aes(x=group,y=Description,size=Count,colour=LogPadjust))+
  geom_point(alpha =0.8,na.rm = TRUE)+
  scale_size(breaks = c(5,10,15,20,90),range = c(1,3),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,2,4,6,8),name='-log10(p_value_adjust)')+  
  # scale_y_discrete(labels=function(x) str_wrap(x, width=150))+
  theme_classic()+labs(x="",y="GO terms",title="Abortion related unique change DEGs BP enrichment")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 2,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1,angle = 90),
        legend.title = element_text(size = 10),
        legend.position ="right",legend.direction = "vertical")
write.table(Enrich_all_plot_noselect, file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/wait_selected_UP_Down_unique_DEGs_all_GO_BP_merge0.05q.txt",row.names=T, col.names=T) 
