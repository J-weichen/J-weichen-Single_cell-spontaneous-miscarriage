rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)#org.Hs.eg.db 的数据类型以及使用简介：可用于数据类型转换
keytypes(org.Hs.eg.db)#查看org.Hs.eg.db数据对象里面包含着各大主流数据库的数据
library(RColorBrewer)

Abortion_celltalk_special_gene<-read.csv("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_CTRL_special_interacting_pair_genelist.csv",head = T)
head(Abortion_celltalk_special_gene)
Abortion_gene_list<-na.omit(unique(Abortion_celltalk_special_gene$Abortion))
CTRL_gene_list<-na.omit(unique(Abortion_celltalk_special_gene$Control))
gene_list<-c(list(Abortion_gene_list),list(CTRL_gene_list))
names(gene_list)<-c("Abortion","CTRL")
#list.files命令将input文件夹下所有文件名输入
for ( fname in c("Abortion","CTRL")){
  # fname<-"Abortion"     
  gene_list0<-gene_list[[fname]]
  print(fname)
  gene_list1<- bitr(gene_list0, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) 
  BP <- enrichGO(gene_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2, readable=T)  
  MF <- enrichGO(gene_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, minGSSize = 2, readable=T) 
  CC <- enrichGO(gene_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2,  readable=T) 
  head(summary(BP))
  CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
  BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
  MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
  head(as.data.frame(BP_simp@result))
  kk <- clusterProfiler::enrichKEGG(gene = gene_list1$ENTREZID,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
  kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
  head(as.data.frame(kk@result))
  ## Reactome pathway enrichment analysis
  rpea <- enrichPathway(gene= gene_list1$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
  head(as.data.frame(rpea@result))
  
  write.table(as.data.frame(CC@result), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/",fname,"_special_CPI_gene_unselect_CC_GO.txt",sep = ""), row.names=T, col.names=T) 
  write.table(as.data.frame(BP@result), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/",fname,"_special_CPI_gene_unselect_BP_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(MF@result), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/",fname,"_special_CPI_gene_unselect_MF_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(CC_simp@result), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/",fname,"_special_CPI_gene_unselect_CC_GO_simp.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(BP_simp@result), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/",fname,"_special_CPI_gene_unselect_BP_GO_simp.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(MF_simp@result), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/",fname,"_special_CPI_gene_unselect_MF_GO_simp.txt",sep = ""), row.names=T, col.names=T) 
  write.table((as.data.frame(kk@result)), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/",fname,"_special_CPI_gene_unselect_KEGG.txt",sep = ""),row.names=T, col.names=T) 
  write.table((as.data.frame(rpea@result)), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/",fname,"_special_CPI_gene_unselect_RPEA.txt",sep = ""),row.names=T, col.names=T) 
  gene_GO_RPEA_erichment_results<-list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
  names(gene_GO_RPEA_erichment_results)<-c("CC","BP","MF","CC_simp","BP_simp","MF_simp","kkeg","rpea")
  saveRDS(gene_GO_RPEA_erichment_results, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/",fname,"_q005_unselect_GO_RPEA_result.rds",sep = ""))
}

##plot enrichment
#construct dataframe for all GO terms 
Abortion_BP <-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/Abortion_special_CPI_gene_unselect_BP_GO_simp.txt",head = T)
CTRL_BP <-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/CTRL_special_CPI_gene_unselect_BP_GO_simp.txt",head = T)

dim(Abortion_BP[which(Abortion_BP$Count>=2),]);dim(CTRL_BP[which(CTRL_BP$Count>=3),])#145  9  126  9
Enrich_all_plot<-as.data.frame(rbind(Abortion_BP,CTRL_BP))
Enrich_all_plot$group<-c(rep("Abortion",nrow(Abortion_BP)),rep("CTRL",nrow(CTRL_BP)))
dim(Enrich_all_plot)#446  10
Enrich_all_plot$LogPadjust<- c(-log10(Enrich_all_plot$p.adjust))
range(Enrich_all_plot$LogPadjust)# 1.30158 10.14360
range(Enrich_all_plot$Count)# 1 16

Enrich_all_number<-data.frame(table(as.character(Enrich_all_plot$Description)))
#Enrich_all_number2<-Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot2<-merge(Enrich_all_plot,Enrich_all_number)
write.table(Enrich_all_plot2, file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/ALL_GO_BP_merge.txt",row.names=T, col.names=T) 

#no selected 
Enrich_all_plot_noselect<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/ALL_GO_BP_merge.txt",header =T) 
head(Enrich_all_plot_noselect)
length(unique(as.character(Enrich_all_plot_noselect$Description)))#415
length(unique(as.character(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>=3),]$Description)))#199
range(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>=3),]$Count)#3 16

#Enrich_all_plot_noselect$Description<-factor(Enrich_all_plot_noselect$Description,levels = rev(c(unique(as.character(Enrich_all_plot_noselect$Description)))))
Enrich_all_plot_noselect2<-Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count >= 3),]
Enrich_all_number <-data.frame(table(as.character(Enrich_all_plot_noselect2$Description)))
colnames(Enrich_all_number)<-c("Description","freq2")
Enrich_all_plot2<-merge(Enrich_all_plot_noselect2,Enrich_all_number)
head(Enrich_all_plot2)
length(unique(as.character(Enrich_all_plot2$Description)))#199
write.table(Enrich_all_plot2, file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/GO_BP_merge_count_more_than_one.txt",row.names=T, col.names=T) 


Enrich_all_plot_noselect$group<-factor(Enrich_all_plot_noselect$group,levels = c("Abortion","CTRL"))
Enrich_all_plot_noselect2<-Enrich_all_plot_noselect2[order(Enrich_all_plot_noselect2$freq,Enrich_all_plot_noselect2$group,decreasing = T ),]
Enrich_all_plot_noselect2$Description<-factor(Enrich_all_plot_noselect2$Description,levels = rev(c(unique(as.character(Enrich_all_plot_noselect2$Description)))))
length(unique(as.character(Enrich_all_plot_noselect2[which(Enrich_all_plot_noselect2$freq>1),]$Description)))#29

plot_BP<-ggplot(Enrich_all_plot_noselect2, aes(x=group,y=Description,size=Count,colour=LogPadjust))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(3,6,9,12,15,18),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,4,8,12),name='-log10(p_adjust)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="Abortion related ICP gene  BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/Enrich_all_BP_noselect.pdf",plot_BP,width = 12, height =30)
write.table(Enrich_all_plot_noselect2, file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/GO_BP_merge_count_more_than_one_wait_selected.txt",row.names=T, col.names=T) 

#for selected  
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/GO_BP_merge_count_more_than_one_selected.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
length(unique(Enrich_all_plot_select$Description))#26

Enrich_all_plot_select$group<-factor(Enrich_all_plot_select$group,levels = c("Abortion","CTRL"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$freq,Enrich_all_plot_select$group,decreasing =F),]
Enrich_all_plot_select$Description<-factor(Enrich_all_plot_select$Description,levels = unique(as.character(Enrich_all_plot_select$Description)))

length(unique(as.character(Enrich_all_plot_select[which(Enrich_all_plot_select$freq>1),]$Description)))#11

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=group,y=Description,size=Count,colour=LogPadjust))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(3,6,9,12,15,18),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,4,8,12),name='-log10(p_adjust)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="Abortion related ICP gene  BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/enrichment_result/Enrich_all_BP_select.pdf",plot_BP,width = 8, height =8)

