rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
options(stringsAsFactors = FALSE)
## Load packages
library(reshape2)
library(dendextend)
library(ggsci)
library(ggpubr)
library(gridExtra)
library(scales)
library(VennDiagram)

##sample information plot in generation
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

#read regulon differential in TFs
D_ac_TFs_active_high_matrix <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/Different_TFs/merge_ALLcell_all_genes_D_ac_TFs_active_high_matrix.txt", header=T)
range(table(D_ac_TFs_active_high_matrix$label))#1 13
head(D_ac_TFs_active_high_matrix)

CTB_D_ac_TFs<-D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$cell_identity %in% c("CTBs_1","CTBs_2")),]
dim(CTB_D_ac_TFs)#68 13
CTB_up_matrix<-CTB_D_ac_TFs[which(CTB_D_ac_TFs$change =="Up"),]
CTB_down_matrix<-CTB_D_ac_TFs[which(CTB_D_ac_TFs$change =="Down"),]


EVT_D_ac_TFs<-D_ac_TFs_active_high_matrix[which(D_ac_TFs_active_high_matrix$cell_identity %in% c("EVTs_1","EVTs_2","EVTs_3")),]
dim(EVT_D_ac_TFs)#206  13
EVT_up_matrix<-EVT_D_ac_TFs[which(EVT_D_ac_TFs$change =="Up"),]
EVT_down_matrix<-EVT_D_ac_TFs[which(EVT_D_ac_TFs$change =="Down"),]

##read TFs selected with activity and expression trend between USM and CTRL
EVT_all_TF <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/evt_all_TF.csv", header = T,row.names= 1)
EVT_up_TF  <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/evt_up_TF.csv", header = T,row.names= 1)
EVT_down_TF  <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/evt_down_TF.csv", header = T,row.names= 1)
length(EVT_all_TF$x);length(EVT_up_TF$x);length(EVT_down_TF$x)##101  27  30
CTB_all_TF <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/ctb_all_TF.csv", header = T,row.names= 1)
CTB_up_TF  <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/ctb_up_TF.csv", header = T,row.names= 1)
CTB_down_TF  <-read.csv(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_5/USM_differential_TFs/ctb_down_TF.csv", header = T,row.names= 1)
length(CTB_all_TF$x);length(CTB_up_TF$x);length(CTB_down_TF$x) #101 26 20

##perform intersection

#EVT
venn <-venn.diagram(list(regulon_Up_EVT=unique(EVT_up_matrix$lable),regulon_Down_EVT=unique(EVT_down_matrix$lable),
                         Up_EVT=EVT_up_TF$x,Down_EVT=EVT_down_TF$x),
                    alpha=c(0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=ppCor[c(1,2,3,6)], 
                    cex = 1.5,cat.col=ppCor[c(1,2,3,6)], cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)

grid.newpage(); 
grid.draw(venn)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/no_merge_EVTs_common_de_regulon_TF_H3K37ac_venn.pdf",width = 8,height = 6)
grid.draw(venn)
dev.off()

#CTB
venn <-venn.diagram(list(regulon_Up_CTB=unique(CTB_up_matrix$lable),regulon_Down_CTB=unique(CTB_down_matrix$lable),
                         Up_CTB=CTB_up_TF$x,Down_CTB=CTB_down_TF$x),
                    alpha=c(0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=ppCor[c(1,2,3,6)], 
                    cex = 1.5,cat.col=ppCor[c(1,2,3,6)], cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)

grid.newpage(); grid.draw(venn)
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/selected_genes/no_merge_CTBs_common_de_regulon_TF_H3K37ac_venn.pdf",width = 8,height = 6)
grid.draw(venn)
dev.off()

#selected unique trend DEGs：CTBs
CTB_up_target<-Reduce(intersect,list(unique(CTB_up_matrix$lable),CTB_up_TF$x))##"TEF"    "TFAP2A"
CTB_down_target<-Reduce(intersect,list(unique(CTB_down_matrix$lable),CTB_down_TF$x))## "TFAP2C"
CTB_target_data<-CTB_D_ac_TFs[which(CTB_D_ac_TFs$label %in% c(CTB_up_target,CTB_down_target)),]
CTB_target_data[,c("lable", "cell_identity","change","foldchange","p.val","p.adj")]


#selected unique trend DEGs
##EVTs
EVT_up_target<-Reduce(intersect,list(unique(EVT_up_matrix$lable),EVT_up_TF$x))##"TFAP2C" "ARNT2"  "SMAD3"  
EVT_down_target<-Reduce(intersect,list(unique(EVT_down_matrix$lable),EVT_down_TF$x))##  "KLF13"  "ELF3"   "GATA2"  "HSF1"   "SNAI1"  "TFAP2A"

Up_special_TF<-setdiff(EVT_up_target,unique(EVT_down_matrix$lable))#"ARNT2" "SMAD3"
Down_special_TF<-setdiff(EVT_down_target,unique(EVT_up_matrix$lable))#"ELF3" "HSF1"

EVT_specia_data<-EVT_D_ac_TFs[which(EVT_D_ac_TFs$label %in% c(Up_special_TF,Down_special_TF)),]
EVT_specia_data[,c("lable", "cell_identity","change","foldchange","p.val","p.adj")]

Up_conflict_TF<-setdiff(EVT_up_target,Up_special_TF)# "TFAP2C"
Down_conflict_TF<-setdiff(EVT_down_target,Down_special_TF)#  "KLF13"  "GATA2"  "SNAI1"  "TFAP2A"
EVT_conflict_data<-EVT_D_ac_TFs[which(EVT_D_ac_TFs$label %in% c(Up_conflict_TF,Down_conflict_TF)),]
EVT_conflict_data[,c("lable", "cell_identity","change","foldchange","p.val","p.adj")]

EVT_target_data<-EVT_D_ac_TFs[which(EVT_D_ac_TFs$label %in% c(EVT_up_target,EVT_down_target)),]
EVT_target_data[,c("lable", "cell_identity","change","foldchange","p.val","p.adj")]
