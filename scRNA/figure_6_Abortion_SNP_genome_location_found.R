rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(vcfR)
vcf <- read.vcfR("/mnt/data/chenwei/gongchen/manuscript/manu_table/dbsnp_144.hg38_nohead.vcf", verbose = FALSE )
rs_ID<- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/SNP_number_355.txt",header = T)
SNP_matrix<-as.data.frame(vcf@fix)
conflict_prefer("which", "Matrix")
target_vcf<-SNP_matrix[which(SNP_matrix$ID %in% rs_ID$SNP_names),]
rs_ID_unknown<-setdiff(rs_ID$SNP_names,target_vcf$ID)
dim(target_vcf)# 302   8
length(rs_ID_unknown)#28
write.table(target_vcf, file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure6/Abortion_related_SNP_302_genome_location_hg38_dbsnp_144.txt",quote=F, row.names=F, col.names=T,sep="\t")
write.table(rs_ID_unknown, file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure6/no_found_28_Abortion_related_SNP_genome_location_in_hg38_dbsnp_144.txt",quote=F, row.names=F, col.names=T,sep="\t")
