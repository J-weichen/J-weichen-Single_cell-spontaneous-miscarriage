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

CTRL_pair_pre_enrich_genes<-as.character(as.data.frame(read.csv('/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/CTRL_pair_pre_enrich_genes_list.csv',sep = ","))[,1])
Abortion_pair_pre_enrich_genes<-as.character(as.data.frame(read.csv('/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_pair_pre_enrich_genes_list.csv',sep = ","))[,1])
gene_list<-c(list(CTRL_pair_pre_enrich_genes),list(Abortion_pair_pre_enrich_genes))
names(gene_list)<-c("CTRL_special","Abortion_special")

#enrichment
for ( gene_name in names(gene_list)){
  #gene_name <- "CTRL_special"   ##test line
  genes_list<-gene_list[[gene_name]]
  print(gene_name)
  genes_list1<- bitr(genes_list, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) #translate into other types ID
  #  gene_GO_RPEA_erichment_results<-list()
  BP <- enrichGO(genes_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2, readable=T)  
  MF <- enrichGO(genes_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, minGSSize = 2, readable=T) 
  CC <- enrichGO(genes_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2,  readable=T) 
  head(summary(BP))
  kk <- clusterProfiler::enrichKEGG(gene = genes_list1$ENTREZID,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
  kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
  head(as.data.frame(kk@result))
  ## Reactome pathway enrichment analysis
  rpea <- enrichPathway(gene= genes_list1$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
  head(as.data.frame(rpea@result))
  write.table(as.data.frame(CC@result), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/",gene_name,"_q005_unselect_CC_GO.txt",sep = ""), row.names=T, col.names=T) 
  write.table(as.data.frame(BP@result), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/",gene_name,"_q005_unselect_BP_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(MF@result), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/",gene_name,"_q005_unselect_MF_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table((as.data.frame(kk@result)), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/",gene_name,"_q005_unselect_KEGG.txt",sep = ""),row.names=T, col.names=T) 
  write.table((as.data.frame(rpea@result)), file=paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/",gene_name,"_q005_unselect_RPEA.txt",sep = ""),row.names=T, col.names=T) 
  gene_GO_RPEA_erichment_results<-list(CC,BP,MF,kk,rpea)
  names(gene_GO_RPEA_erichment_results)<-c("CC","BP","MF","kkeg","rpea")
  saveRDS(gene_GO_RPEA_erichment_results, file = paste("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/",gene_name,"_q005_unselect_GO_RPEA_result.rds",sep = ""))
}

