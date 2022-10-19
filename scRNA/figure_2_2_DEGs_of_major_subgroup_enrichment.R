rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)#org.Hs.eg.db 的数据类型以及使用简介：可用于数据类型转换
keytypes(org.Hs.eg.db)#查看org.Hs.eg.db数据对象里面包含着各大主流数据库的数据
library(KEGG.db)

#analysis for each genelist
#For default 
files<- list.files("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/default_DEGs/Pos_genes_01q")
#list.files命令将input文件夹下所有文件名输入
data_list<-list();cnane<-c()
for ( f in files[1:length(files)]){
  FF<-paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/default_DEGs/Pos_genes_01q/",f,sep="")
  column <- read.delim(FF,header = T,sep = "\t")[,c("cluster","gene")]
  if (is.na(column[1,]$cluster)) {
    next
  }
  data_list<-c(data_list,list(as.character(column$gene)))
  print(unique(as.character(column$cluster)))
  cnane<-c(cnane,paste(unique(as.character(column$cluster)),"_up",sep=""))
}
length(data_list);length(cnane)
names(data_list)<-cnane

files<- list.files("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/default_DEGs/Neg_genes_01q")
data_list2<-list();cnane2<-c()
for ( f in files[1:length(files)]){
  FF<-paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/default_DEGs/Neg_genes_01q/",f,sep="")
  column <- read.delim(FF,header = T,sep = "\t")[,c("cluster","gene")]
  if (is.na(column[1,]$cluster)) {
    next
  }
  data_list2<-c(data_list2,list(as.character(column$gene)))
  print(unique(as.character(column$cluster)))
  cnane2<-c(cnane2,paste(unique(as.character(column$cluster)),"_down",sep=""))
}
length(data_list2);length(cnane2)
names(data_list2)<-cnane2
data_list_two<-c(data_list,data_list2)
#enrichment
for ( DEGs_name in names(data_list_two)){
  #DEGs_name <- "HCs_up"   ##test line
  DEGs_list<-data_list_two[[DEGs_name]]
  print(DEGs_name)
  DEGs_list1<- bitr(DEGs_list, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) #translate into other types ID
  #  gene_GO_RPEA_erichment_results<-list()
  BP <- enrichGO(DEGs_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2, readable=T)  
  MF <- enrichGO(DEGs_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, minGSSize = 2, readable=T) 
  CC <- enrichGO(DEGs_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2,  readable=T) 
  head(summary(BP))
  kk <- clusterProfiler::enrichKEGG(gene = DEGs_list1$ENTREZID,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
  kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
  head(as.data.frame(kk@result))
  ## Reactome pathway enrichment analysis
  rpea <- enrichPathway(gene= DEGs_list1$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
  head(as.data.frame(rpea@result))
  write.table(as.data.frame(CC@result), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup/",DEGs_name,"_q005_unselect_CC_GO.txt",sep = ""), row.names=T, col.names=T) 
  write.table(as.data.frame(BP@result), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup/",DEGs_name,"_q005_unselect_BP_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(MF@result), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup/",DEGs_name,"_q005_unselect_MF_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table((as.data.frame(kk@result)), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup/",DEGs_name,"_q005_unselect_KEGG.txt",sep = ""),row.names=T, col.names=T) 
  write.table((as.data.frame(rpea@result)), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup/",DEGs_name,"_q005_unselect_RPEA.txt",sep = ""),row.names=T, col.names=T) 
  gene_GO_RPEA_erichment_results<-list(CC,BP,MF,kk,rpea)
  names(gene_GO_RPEA_erichment_results)<-c("CC","BP","MF","kkeg","rpea")
  saveRDS(gene_GO_RPEA_erichment_results, file = paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/major_subgroup/",DEGs_name,"_q005_unselect_GO_RPEA_result.rds",sep = ""))
}

##analysis for whole genelist
DEGs_merge_data <- read.table(file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/DEG_01q_default_merge.txt", header=T)
DEGs_merge_data$cluster<-as.character(DEGs_merge_data$cluster)
table(DEGs_merge_data$cluster)
# Bs CTBs_0 CTBs_1 CTBs_2 CTBs_3   Endo    Epi EVTs_1 EVTs_2 EVTs_3    FBs    HCs   MyCs    NKs STBs_1 STBs_2   STCs     Ts 
# 16    659    587    508    280    164     33    932   1035    653    448   1171     96     82    810    631    297     41 
DEGs_merge_data$trend<-"Up"
DEGs_merge_data[which(DEGs_merge_data$avg_logFC <0),]$trend<-"Down"
DEGs_merge_data$cluster_trend <- paste0(DEGs_merge_data$cluster,"_",DEGs_merge_data$trend)

#seperated genes
#merge_data_rm0<-DEGs_merge_data[which(!(DEGs_merge_data$cluster %in% c("Ery"))),]
Vagina_up<-unique(as.character(DEGs_merge_data[which(DEGs_merge_data$avg_logFC > 0),]$gene))
length(Vagina_up)# 130
Vagina_down<-unique(as.character(DEGs_merge_data[which(DEGs_merge_data$avg_logFC < 0),]$gene))
length(Vagina_down)# 147

Up_gene<-unique(c(Vagina_up));Down_gene<-unique(c(Vagina_down))
length(Up_gene);length(Down_gene)
# 1551 #1601
Up_uniq<-setdiff(Up_gene,Down_gene)
Down_uniq<-setdiff(Down_gene,Up_gene)
length(Up_uniq);length(Down_uniq)
#1258 #1308
DEGs_merge_data$trend2<-DEGs_merge_data$trend
DEGs_merge_data$trend2<-ifelse(DEGs_merge_data$gene %in% Up_uniq,"Up_uniq",ifelse(DEGs_merge_data$gene %in% Down_uniq,"Down_uniq","contract"))
#save result
write.table(as.data.frame(DEGs_merge_data), file="/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/DEG_01q_default_merge_new.txt",row.names=T, col.names=T) 
Venn_gene_list<-c(list(Up_uniq),list(Down_uniq))
names(Venn_gene_list)<-c("Up_uniq","Down_uniq")
saveRDS(Venn_gene_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/figure_two_Venn_gene_list.rds")
figure_two_Venn_gene_list<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/figure_two_Venn_gene_list.rds")
#enrichment
for ( DEGs_name in names(figure_two_Venn_gene_list)){
  #DEGs_name <- "Up_uniq"  ##test line
  DEGs_list<-figure_two_Venn_gene_list[[DEGs_name]]
  print(DEGs_name)
  DEGs_list1<- bitr(DEGs_list, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) #translate into other types ID
  #  gene_GO_RPEA_erichment_results<-list()
  BP <- enrichGO(DEGs_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2, readable=T)  
  MF <- enrichGO(DEGs_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, minGSSize = 2, readable=T) 
  CC <- enrichGO(DEGs_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2,  readable=T) 
  head(summary(BP))
  kk <- clusterProfiler::enrichKEGG(gene = DEGs_list1$ENTREZID,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
  kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
  head(as.data.frame(kk@result))
  ## Reactome pathway enrichment analysis
  rpea <- enrichPathway(gene= DEGs_list1$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
  head(as.data.frame(rpea@result))
  write.table(as.data.frame(CC@result), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/merge_DEG_enrichment/",DEGs_name,"_q005_unselect_CC_GO.txt",sep = ""), row.names=T, col.names=T) 
  write.table(as.data.frame(BP@result), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/merge_DEG_enrichment/",DEGs_name,"_q005_unselect_BP_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(MF@result), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/merge_DEG_enrichment/",DEGs_name,"_q005_unselect_MF_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table((as.data.frame(kk@result)), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/merge_DEG_enrichment/",DEGs_name,"_q005_unselect_KEGG.txt",sep = ""),row.names=T, col.names=T) 
  write.table((as.data.frame(rpea@result)), file=paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/merge_DEG_enrichment/",DEGs_name,"_q005_unselect_RPEA.txt",sep = ""),row.names=T, col.names=T) 
  gene_GO_RPEA_erichment_results<-list(CC,BP,MF,kk,rpea)
  names(gene_GO_RPEA_erichment_results)<-c("CC","BP","MF","kkeg","rpea")
  saveRDS(gene_GO_RPEA_erichment_results, file = paste("/mnt/data/chenwei/gongchen/3.seurat_result/figure_two_result/enrichment_result/merge_DEG_enrichment/",DEGs_name,"_q005_unselect_GO_RPEA_result.rds",sep = ""))
}

