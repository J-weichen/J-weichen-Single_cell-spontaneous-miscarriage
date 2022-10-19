rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(Seurat)
library(future)

library(data.table)
library(scales)
library(dplyr)
library(psych)

library(ggplot2)
library(cowplot)
library(ggsci)
library(grid)
library(gridExtra)

library(pheatmap)
library(ComplexHeatmap)

grid.newpage()

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
#n=5;barplot(rep(1,n), col=colorRampPalette(colors = c('red', 'white'))( n ))
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
show_col(ppCor_all2)

##plot cells
Cells_col<-colorRampPalette(colors = rev(Cells_col_raw))(80)
length(unique(Cells_col))
barplot(rep(1,80), col=Cells_col)
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(80)
barplot(rep(1,80), col=Cells_col)

#多核并行运算。加速差异计算https://github.com/satijalab/seurat/issues/1278
plan(strategy = "multicore", workers = 42)
options(future.globals.maxSize = 1500 * 1024^12)

target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
target_final<-subset(x = target,subset = final_major_subgroup== "Troph_mix_group",invert = TRUE)
major_order<-c("Trophoblast_KRT7","Myeloid_Cell_AIF1","Stromal_cells_DCN",
               "Endothelial_Cell_PECAM1","Epithelial_Cell_EPCAM_PAEP","Natural_killer_Cell_NCAM1","T_Cell_CD3D","B_Cell_CD79A","Erythrocyte_HBA1","Mast_Cell_MS4A2")
subgroup_order<-c("CTBs_0","CTBs_1","CTBs_2","CTBs_3","EVTs_1","EVTs_2","EVTs_3","STBs_1","STBs_2","MyCs","HCs","STCs","FBs","Endo","Epi", "NKs","Ts","Bs","Ery","Masts")
target_final$final_major_subgroup_brief<- factor(target_final$final_major_subgroup_brief, levels=subgroup_order,ordered=TRUE)
target_final$final_major_group_brief<- factor(target_final$final_major_group_brief, levels=major_order,ordered=TRUE)

#plot  group
DimPlot(target_final, group.by = "final_major_subgroup_brief",cols =Cells_col)
DimPlot(target_final, group.by = "final_major_group_brief",label = F,cols = ppCor_all2)
DimPlot(object = target_final, group.by = "final_major_group_brief", split.by = "Treat",reduction = "umap",cols = ppCor_all2,ncol = 2)

#target_final_N4<-subset(x = target_final,subset = sample== "N4")
target_final_A4<-subset(x = target_final,subset = sample== "A4")

target_final_selected<-target_final_A4
DefaultAssay(target_final_selected) <- "RNA"    ## very important
#normalization
target_final_selected <- SCTransform(target_final_selected, verbose = FALSE)
target_final_selected <- RunPCA(target_final_selected)
ElbowPlot(object = target_final_selected)
target_final_selected <- RunUMAP(target_final_selected, dims = 1:25)
target_final_selected<- RunTSNE(target_final_selected, dims = 1:25)
saveRDS(target_final_selected, file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_0713/target_final_A4_forced_cell_low_cell_double_both_remove_seurat.rds")

table(target_final_selected$final_major_group_brief)
DimPlot(target_final_selected, group.by = "final_major_subgroup_brief",cols =Cells_col)
DimPlot(target_final_selected, group.by = "final_major_group_brief",label = F,cols = ppCor_all2)