rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
## reading gene in diff clusters
k27_USM_gene_module<-read.csv(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/k27_CTB2EVT_usm_gene_4cluster.csv", sep = ',', row.names = 1, header = F)
k27_USM_gene_module$V2<-paste0("cluster",k27_USM_gene_module$V2)
colnames(k27_USM_gene_module)<-c("k27_USM_cluster")

k27_CTRL_gene_module<-read.csv(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/k27_CTB2EVT_ct_gene_4cluster.csv", sep = ',', row.names = 1, header = F)
k27_CTRL_gene_module$V2<-paste0("cluster",k27_CTRL_gene_module$V2)
colnames(k27_CTRL_gene_module)<-c("k27_CTRL_cluster")

RNA_USM_gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/gene_Cluster_USM_CRB2EVT_RNA_K27ac.txt", sep = '\t', row.names = 1, header = T)
rownames(RNA_USM_gene_module)<-RNA_USM_gene_module$row_genes
colnames(RNA_USM_gene_module)<-c("gene_name","RNA_USM_cluster")

RNA_CTRL_gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/gene_Cluster_CTRL_CRB2EVT_RNA_K27ac.txt", sep = '\t', row.names = 1, header = T)
rownames(RNA_CTRL_gene_module)<-RNA_CTRL_gene_module$row_genes
colnames(RNA_CTRL_gene_module)<-c("gene_name","RNA_CTRL_cluster")

dim(k27_USM_gene_module);dim(k27_CTRL_gene_module);dim(RNA_USM_gene_module);dim(RNA_CTRL_gene_module)
k27_gene_module<-merge(k27_USM_gene_module,k27_CTRL_gene_module,by=0)
RNA_gene_module<-merge(RNA_USM_gene_module,RNA_CTRL_gene_module,by=0,all=T)
dim(RNA_gene_module)
All_gene_module<-merge(k27_gene_module,RNA_gene_module[,c(1,3,5)],by="Row.names")
rownames(All_gene_module)<-All_gene_module$Row.names
All_gene_module<-All_gene_module[,-1]
head(All_gene_module)
gene_signature<-c("TIMP3","SERPINE1","CDH1","STMN1","THBD")
All_gene_module[gene_signature,]

All_gene_module$k27_USM_cluster<-ifelse(All_gene_module$k27_USM_cluster =="cluster1","Up_Down_Up",
                                        ifelse(All_gene_module$k27_USM_cluster =="cluster2","Up_along",
                                               ifelse(All_gene_module$k27_USM_cluster =="cluster3","Down_Up_Down","Down_along")))
All_gene_module$k27_CTRL_cluster<-ifelse(All_gene_module$k27_CTRL_cluster =="cluster1","Up_along",
                                         ifelse(All_gene_module$k27_CTRL_cluster =="cluster2","Up_pre_along",
                                                ifelse(All_gene_module$k27_CTRL_cluster =="cluster3","Down_Up_Down","Down_along")))
All_gene_module$RNA_USM_cluster <-ifelse(All_gene_module$RNA_USM_cluster  =="cluster1","Up_along",
                                         ifelse(All_gene_module$RNA_USM_cluster  =="cluster2","Down_along",
                                                ifelse(All_gene_module$RNA_USM_cluster  =="cluster3","Down_Up","Up_Down_Up")))
All_gene_module$RNA_CTRL_cluster <-ifelse(All_gene_module$RNA_CTRL_cluster  =="cluster1","Down_along",
                                          ifelse(All_gene_module$RNA_CTRL_cluster  =="cluster2","Up_along",
                                                 ifelse(All_gene_module$RNA_CTRL_cluster  =="cluster3","Down_mini","Up_mini")))

head(All_gene_module)
gene_signature<-c("MKI67","SERPINE1","TIMP3","STMN1","TIMP2","TIMP1","HLA-G","THBD","MALAT1","TGFB1","LAMA4")
gene_signature<-c("PAPPA2","CDH1","KRT19","KRT7","CYP19A1","EGFR","ERVFRD-1","RPL39","TIMP4","MMP11","MMP15","MMP17")
gene_signature<-c("KRT19","MMP15","SERPINE1","MALAT1","TGFB1")
All_gene_module[gene_signature,]

###using USM CTRL toghether cluster
d_num<-16
select_RNA_gene_module<-read.table(file =  paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/gene_Cluster_",d_num,"_ALL_CRB2EVT_RNA_K27ac.txt"), sep = '\t', row.names = 1, header = T)
head(select_RNA_gene_module)
select_RNA_gene_module$module<-select_RNA_gene_module$Cluster
table(select_RNA_gene_module$Cluster)
select_RNA_gene_module[which(select_RNA_gene_module$module %in% c("cluster10", "cluster11")),]$module<-"cluster10"
select_RNA_gene_module[which(select_RNA_gene_module$module %in% c("cluster12", "cluster13","cluster14", "cluster15", "cluster16")),]$module<-"cluster11"


select_RNA_gene_module$module<-factor(select_RNA_gene_module$module,levels = paste0("cluster",1:11))
select_RNA_gene_module<-select_RNA_gene_module[order(select_RNA_gene_module$module,decreasing = F),]
write.table(select_RNA_gene_module, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_rearrange_gene_Cluster_ALL_CRB2EVT_RNA_K27ac.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
