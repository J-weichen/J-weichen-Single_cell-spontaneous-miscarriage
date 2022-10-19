rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(scales)
library(ggsci)
library(pheatmap)
library(VennDiagram)
library(UpSetR)
library(grid)
library(futile.logger)
library(reshape2)
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:46))]
show_col(ppCor_all2)
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(80)
show_col(Cells_col)



#selected candiidated group for venn plot
#gene file0 reading peak annotation files
CTRL_EVT1<-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/cobatch_abortion_related_peak/CTRL-EVT1_specific_peak.csv",header =T,sep=",",check.names = F)
CTRL_CTB1<-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/cobatch_abortion_related_peak/CTRL-CTB1_specific_peak.csv",header =T,sep=",",check.names = F)
CTRL_other<-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/cobatch_abortion_related_peak/CTRL-other_specific_peak.csv",header =T,sep=",",check.names = F)
Abortion_EVT2<-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/cobatch_abortion_related_peak/Abortion-EVT2_specific_peak.csv",header =T,sep=",",check.names = F)
Abortion_CTB2<-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/cobatch_abortion_related_peak/Abortion-CTB2_specific_peak.csv",header =T,sep=",",check.names = F)
Abortion_other<-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/cobatch_abortion_related_peak/Abortion-other_specific_peak.csv",header =T,sep=",",check.names = F)

#merge alldata
CTRL_EVT1$cell_type<-"CTRL_EVT1";CTRL_CTB1$cell_type<-"CTRL_CTB1";CTRL_other$cell_type<-"CTRL_other"
Abortion_EVT2$cell_type<-"Abortion_EVT2";Abortion_CTB2$cell_type<-"Abortion_CTB2";Abortion_other$cell_type<-"Abortion_other"
All_peak_file<-rbind(CTRL_EVT1,CTRL_CTB1,CTRL_other,Abortion_EVT2,Abortion_CTB2,Abortion_other)

All_CTRL_peak_file<-rbind(CTRL_EVT1,CTRL_CTB1,CTRL_other)
All_Abortion_peak_file<-rbind(Abortion_EVT2,Abortion_CTB2,Abortion_other)
All_peak_file$location<- unlist(lapply(strsplit(as.character(All_peak_file$Annotation)," "), function(x) x[1]))
All_CTRL_peak_file$location<- unlist(lapply(strsplit(as.character(All_CTRL_peak_file$Annotation)," "), function(x) x[1]))
All_Abortion_peak_file$location<- unlist(lapply(strsplit(as.character(All_Abortion_peak_file$Annotation)," "), function(x) x[1]))

dim(All_peak_file);dim(All_CTRL_peak_file);dim(All_Abortion_peak_file)# 5608 3020  2588 
length(unique(All_peak_file$`Gene Name`))# 3190
length(unique(All_CTRL_peak_file$`Gene Name`))# 1822
length(unique(All_Abortion_peak_file$`Gene Name`))#1752

#selected candiidated group for venn plot
#gene file 1: single cell DEGs
merge_data0<-read.table(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge.txt",header =T,sep="\t",check.names = F)
head(merge_data0)
sc_Up_gene<-unique(as.character(merge_data0[which(merge_data0$avg_logFC > 0),]$gene))
length(sc_Up_gene)# 1528
sc_Down_gene<-unique(as.character(merge_data0[which(merge_data0$avg_logFC < 0),]$gene))
length(sc_Down_gene)# 1627

#gene file 2: DEGs from bulk
bulk_merge_data<-read.table(file="/mnt/data/chenwei/gongchen/2.map_result/count_file/4.count/Abortion_vs_CTRL.DEG_information_pvalue005_FC1.5.txt",header =T,sep="\t") 
head(bulk_merge_data);dim(bulk_merge_data)#5316    4
bulk_Up_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange > 0),]$ID))
length(bulk_Up_gene)# 2748
bulk_Down_gene<-unique(as.character(bulk_merge_data[which(bulk_merge_data$log2FoldChange < 0),]$ID))
length(bulk_Down_gene)# 2568

#gene file 3: single cell dacTFs
cell_name="ALLcell_all_genes"
D_ac_TFs_active_high_matrix <- read.table( paste0("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Different_TFs/merge_",cell_name,"_D_ac_TFs_active_high_matrix.txt"), header=T)
head(D_ac_TFs_active_high_matrix);dim(D_ac_TFs_active_high_matrix)# 1026   13
Up_acTF<-as.character(unique(D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$change == "Up"),]$label))
Down_acTF<-as.character(unique(D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$change == "Down"),]$label))

#gene file 4：interaction special moleculer
Abortion_inter<-as.character(read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/Abortion_pair_pre_enrich_genes_list.csv",header =T,sep=",",check.names = F)[,1])
CTRL_inter<-as.character(read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_4/subcelltype/CTRL_pair_pre_enrich_genes_list.csv",header =T,sep=",",check.names = F)[,1])

listinput<-c(list(unique(All_CTRL_peak_file$`Gene Name`)),list(unique(All_Abortion_peak_file$`Gene Name`)),
             list(sc_Up_gene),list(sc_Down_gene),list(bulk_Up_gene),list(bulk_Down_gene),list(Up_acTF),list(Down_acTF),list(Abortion_inter),list(CTRL_inter))
names(listinput)<-c("CTRL_peak","Abortion_peak","sc_Up","sc_Down","bulk_Up","bulk_Down","Up_acTF","Down_acTF","Abortion_inter","CTRL_inter" )


#FOR UP_gene VERSUS Down_gene
venn <-venn.diagram(listinput[c("Abortion_peak","CTRL_peak","bulk_Up","bulk_Down")],
                    alpha=c(0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","navy","orange","pink"), 
                    cex = 1.5,cat.col=c("red","navy","orange","pink"), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage(); grid.draw(venn)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Venn_peak_and_bulk_DEGs_for_Abortion.pdf",width = 8,height = 6)
grid.draw(venn)
dev.off()

overlap1<-setdiff(Reduce(intersect,listinput[c("Abortion_peak","bulk_Up")]),c(as.character(unlist(listinput["CTRL_peak"])),as.character(unlist(listinput["bulk_Down"]))))
overlap2<-setdiff(Reduce(intersect,listinput[c("CTRL_peak","bulk_Down")]),c(as.character(unlist(listinput["Abortion_peak"])),as.character(unlist(listinput["bulk_Up"]))))
table(All_peak_file[which(All_peak_file$`Gene Name` %in% overlap1),]$location)
#3'         exon   Intergenic       intron   non-coding promoter-TSS          TTS 
#2            3           89          162            3           15            7 
table(All_peak_file[which(All_peak_file$`Gene Name` %in% overlap2),]$location)
# 3'         exon   Intergenic       intron   non-coding promoter-TSS          TTS 
# 8           25           77          181            2           29           23 

venn <-venn.diagram(listinput[c("Abortion_peak","CTRL_peak","sc_Up","sc_Down")],
                    alpha=c(0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","navy","orange","pink"), 
                    cex = 1.5,cat.col=c("red","navy","orange","pink"), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage(); grid.draw(venn)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Venn_peak_and_sc_DEGs_for_Abortion.pdf",width = 8,height = 6)
grid.draw(venn)
dev.off()
overlap3<-setdiff(Reduce(intersect,listinput[c("Abortion_peak","sc_Up")]),c(as.character(unlist(listinput["CTRL_peak"])),as.character(unlist(listinput["sc_Down"]))))
overlap4<-setdiff(Reduce(intersect,listinput[c("CTRL_peak","sc_Down")]),c(as.character(unlist(listinput["Abortion_peak"])),as.character(unlist(listinput["sc_Up"]))))
table(All_peak_file[which(All_peak_file$`Gene Name` %in% overlap3),]$location)
#  3'           5'         exon   Intergenic       intron promoter-TSS          TTS 
#  1            1            2           65          130            9            1 

table(All_peak_file[which(All_peak_file$`Gene Name` %in% overlap4),]$location)
#  3'           5'         exon   Intergenic       intron promoter-TSS          TTS 
#  6            2           20           49           82           14           14 


venn <-venn.diagram(listinput[c("Abortion_inter","CTRL_inter","sc_Up","sc_Down")],
                    alpha=c(0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","navy","orange","pink"), 
                    cex = 1.5,cat.col=c("red","navy","orange","pink"), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage(); grid.draw(venn)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Venn_inter_molecule_and_sc_DEGs_for_Abortion.pdf",width = 8,height = 6)
grid.draw(venn)
dev.off()

overlap5<-setdiff(Reduce(intersect,listinput[c("Abortion_inter","sc_Up")]),c(as.character(unlist(listinput["CTRL_inter"])),as.character(unlist(listinput["sc_Down"]))))
overlap6<-setdiff(Reduce(intersect,listinput[c("CTRL_inter","sc_Down")]),c(as.character(unlist(listinput["Abortion_inter"])),as.character(unlist(listinput["sc_Up"]))))
overlap5;overlap6
#[1] "FCGR2A" "PDGFRB" "EGFR"   "CCL13" 
#[1] "SELP" "CD24" "TNF" 


venn <-venn.diagram(listinput[c("Abortion_inter","CTRL_inter","bulk_Up","bulk_Down")],
                    alpha=c(0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","navy","orange","pink"), 
                    cex = 1.5,cat.col=c("red","navy","orange","pink"), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage(); grid.draw(venn)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_2/Venn_inter_molecule_and_bulk_DEGs_for_Abortion.pdf",width = 8,height = 6)
grid.draw(venn)
dev.off()

overlap7<-setdiff(Reduce(intersect,listinput[c("Abortion_inter","bulk_Up")]),c(as.character(unlist(listinput["CTRL_inter"])),as.character(unlist(listinput["bulk_Down"]))))
overlap8<-setdiff(Reduce(intersect,listinput[c("CTRL_inter","bulk_Down")]),c(as.character(unlist(listinput["Abortion_inter"])),as.character(unlist(listinput["bulk_Up"]))))
overlap7;overlap8
#[1] "CCL8"     "ADORA1"   "CCL21"    "FCGR2A"   "CD274"    "ENTPD1"   "CCL13"    "HRH1"     "TGFA"     "PDCD1LG2"
#[1] "COL2A1" "EFNB2"  "WNT7B"  "FGFR3"  "CD34"   "AGTR1"  "CD24"   "TNF"  


#upset plot for all six lists of DMRs
#for upset plot
upset_dataframe<-as.data.frame(fromList(listinput))
dim(upset_dataframe)# 7338   10
upset_dataframe[1:10,1:5]
upset_dataframe_rowSums<-rowSums(upset_dataframe)
range(upset_dataframe_rowSums)#1:6
upset_dataframe_colSums<-colSums(upset_dataframe)
range(upset_dataframe_colSums)# 45 2748
upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]
upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400)

##以下未跑

list2_q1<-list(query = intersects, params = list("Kids_hyper", "Mother_hyper","Father_hyper"), color = "red", active = T,query.name = "common hyper DMRs")
list2_q2<-list(query = intersects, params = list("Kids_hypo","Mother_hypo","Father_hypo"), color = "navy", active = T,query.name = "common hypo DMRs")
list2_q3<-list(query = intersects, params = list("Kids_hyper", "Mother_hyper"), color = "orange", active = T,query.name = "mother_SF_DMRs")
list2_q4<-list(query = intersects, params = list("Kids_hypo","Mother_hypo"), color = "orange", active = T,query.name = "mother_SF_DMRs")
list2_q5<-list(query = intersects, params = list("Kids_hyper", "Father_hyper"), color = "purple", active = T,query.name = "father_SF_DMRs")
list2_q6<-list(query = intersects, params = list("Kids_hypo","Father_hypo"), color = "purple", active = T,query.name = "father_SF_DMRs")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/upsetplot_six_lists_DMRs_AMA_vs_Young_200bin_q005_dif15.pdf",width = 10,height=8)

upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400,
      sets=c("Kids_hyper","Mother_hyper","Father_hyper","Kids_hypo","Mother_hypo","Father_hypo"),
      keep.order = TRUE,
      query.legend = "top",
      sets.bar.color = "brown",
      #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
      shade.color="pink",
      #      matrix.color="purple",
      order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
      point.size = 3,line.size = 1.3,
      mainbar.y.label = "DMRs number Intersections", sets.x.label = "DMRs number per subset",
      mb.ratio = c(0.60, 0.40),
      text.scale = c(1.5, 1.5,1.2,1.5,1.5,1),
      show.numbers = 'yes',
      queries = list(list2_q1,list2_q2,list2_q3,list2_q4,list2_q5,list2_q6)
)
dev.off()
