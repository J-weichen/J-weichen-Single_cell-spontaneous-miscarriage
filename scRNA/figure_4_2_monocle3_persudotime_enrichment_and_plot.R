rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so')
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
dyn.load('/mnt/data/chenwei/software/share_road/libudunits2/libudunits2.so.0')
dyn.load('/mnt/data/chenwei/software/share_road/libproj/libproj.so.0')
dyn.load('/mnt/data/chenwei/software/share_road/libgdal.so.20')
dyn.load('/mnt/data/chenwei/software/share_road/libgeos-3.4.2.so')
dyn.load('/mnt/data/chenwei/software/share_road/libgeos_c.so.1')

library(Seurat)
library(monocle3)
library(scales)
library(ggsci)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dendextend)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)

pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

##plot maker genes for cells
seurat_scRNA_final<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_seurat_object_rearrange_monocle3_pseudotime_add.rds")
seurat_scRNA_final
#metadata_for_all_cells<-seurat_scRNA_final@meta.data
#Umap_coordinate<-data.frame(Embeddings(seurat_scRNA_final[["umap"]]))
#target_coordinate2<-merge(metadata_for_all_cells,Umap_coordinate,by=0)
#target_coordinate3<-target_coordinate2[,c("Row.names","sample_code","Treat","re_annotation_TFac","re_anno_TFac_major","pseudotime","UMAP_1","UMAP_2")]
#write.table(target_coordinate3, file = "/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_seurat_object_rearrange_monocle3_UMAP_distribution.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

DefaultAssay(seurat_scRNA_final) <- "RNA"
seurat_scRNA_final <- NormalizeData(seurat_scRNA_final, verbose = FALSE)

final_anno_maker_plot <- FeaturePlot(object = seurat_scRNA_final, features = "MKI67",cols= c("lightgrey", "red"))
final_anno_plot01<-final_anno_maker_plot+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                                plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                                axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                                legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
final_anno_plot02<-final_anno_maker_plot+theme_bw()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                     plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                     axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                     legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_seurat_MKI67_Trophoblast.pdf",final_anno_maker_plot,width =20, height =20)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_seurat_Trophoblast_monocle3_MKI67_no_text.png",final_anno_plot01,width = 8, height =8)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_seurat_Trophoblast_monocle3_MKI67_no_text2.png",final_anno_plot02,width = 8, height =8)


#major_submaker<-c("MTHFR")
#"VEGFA","MTHFR"
major_submaker<-c("KRT7","PCNA","PAGE4","PEG10","CDH1","ERVFRD-1","CYP19A1","HLA-G","PAPPA2")
#Trophoblast_maker<-c("KRT7","GATA3","PERP","XAGE2","CDH1","EGFR","PAGE4","HLA-G","PAPPA2","CGA","CSH2","ERVFRD-1")
final_anno_maker_plot <- FeaturePlot(object = seurat_scRNA_final, features = major_submaker,cols= c("grey", "red"),ncol=3)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_makers_Trophoblast_monocle3_seurat.pdf",final_anno_maker_plot,width =20, height =20)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_makers_Trophoblast_monocle3_seurat.png",final_anno_maker_plot,width =20, height =20)

##selected final cell rds
cds_used<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/All_Trophoblast_harmony_monocle3_pseudotime_object_final_used.rds")
cds_final <- choose_cells(cds_used)
saveRDS(cds_final,file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_fielt_final_monocle3_object.rds")

##formal building
cds_final<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_fielt_final_monocle3_object.rds")
##24099 genes in 32226 cells
##plot maker genes for cells
major_submaker<-c("KRT7","PCNA","PAGE4","PEG10","CDH1","ERVFRD-1","CYP19A1","HLA-G","PAPPA2")
#Trophoblast_maker<-c("KRT7","GATA3","PERP","XAGE2","CDH1","EGFR","PAGE4","HLA-G","PAPPA2","CGA","CSH2","ERVFRD-1")
final_anno_maker_plot<-plot_cells(cds_final,genes= major_submaker,label_cell_groups=FALSE,show_trajectory_graph=FALSE)+ scale_color_gradient2(low = "gray", mid = "white", high =  "red")
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_makers_Trophoblast_monocle3.pdf",final_anno_maker_plot,width =20, height =20)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_makers_Trophoblast_monocle3.png",final_anno_maker_plot,width =20, height =20)

                  
final_anno_plot0<-plot_cells(cds_final, color_cells_by = "pseudotime",show_trajectory_graph =FALSE,  label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
final_anno_plot1<-plot_cells(cds_final, color_cells_by = "re_annotation_TFac",show_trajectory_graph =FALSE,   label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=my_morandi_colors)
final_anno_plot2<-plot_cells(cds_final, color_cells_by = "Treat", show_trajectory_graph =FALSE,  label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=my_morandi_colors[c(1,21)])
final_anno_plot1_split<-final_anno_plot1+facet_wrap(~ Treat)
final_anno_plot2_split<-final_anno_plot2+facet_wrap(~ Treat)

ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_pseudotime.pdf",final_anno_plot0,width = 8, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_annotation.pdf",final_anno_plot1,width = 8, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_Treat.pdf",final_anno_plot2,width = 8, height =6)

ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_pseudotime.png",final_anno_plot0,width = 8, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_annotation.png",final_anno_plot1,width = 8, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_Treat.png",final_anno_plot2,width = 8, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_Treat_split.pdf",final_anno_plot2_split,width = 13, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_annotation_split.pdf",final_anno_plot1_split,width = 13, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_Treat_split.png",final_anno_plot2_split,width = 13, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_annotation_split.png",final_anno_plot1_split,width = 13, height =6)

final_anno_plot3<-plot_cells(cds_final,genes= "MKI67",label_cell_groups=FALSE,show_trajectory_graph=FALSE)+ scale_color_gradient2(low = "gray", mid = "white", high =  "red")
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_MKI67.pdf",final_anno_plot3,width = 8, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_monocle3_MKI67.png",final_anno_plot3,width = 8, height =6)

#######cell annotation without lables
final_anno_plot01<-final_anno_plot0+theme_bw()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                     plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                     axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                     legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
final_anno_plot02<-final_anno_plot01+NoLegend()

final_anno_plot11<-final_anno_plot1+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                                plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                                axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                                legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
final_anno_plot21<-final_anno_plot2+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                                plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                                axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                                legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
final_anno_plot1_split1<-final_anno_plot1_split+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                                            plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                                            axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                                            legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
final_anno_plot2_split1<-final_anno_plot2_split+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                                            plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                                            axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                                            legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
final_anno_plot31<-final_anno_plot3+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                                plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                                axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                                legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())
final_anno_plot32<-final_anno_plot3+theme_bw()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                                plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                                axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                                legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())

ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_monocle3_pseudotime_no_text.png",final_anno_plot01,width = 8.5, height =8)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_monocle3_pseudotime_no_text2.png",final_anno_plot02,width = 8, height =8)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_monocle3_annotation_no_text.png",final_anno_plot11,width = 8, height =8)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_monocle3_Treat_no_text.png",final_anno_plot21,width = 8, height =8)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_monocle3_Treat_split_no_text.png",final_anno_plot2_split1,width = 12, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_monocle3_annotation_split_no_text.png",final_anno_plot1_split1,width = 12, height =6)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_monocle3_MKI67_no_text.png",final_anno_plot31,width = 8, height =8)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_monocle3_MKI67_no_text2.png",final_anno_plot32,width = 8, height =8)

##for revison
final_anno_plot2<-plot_cells(cds_final, color_cells_by = "Treat", show_trajectory_graph =FALSE,  label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=c("#FFCBCD","#7FC3E2"))

final_anno_plot2_split<-final_anno_plot2+facet_wrap(~ Treat)

final_anno_plot2_split1<-final_anno_plot2_split+theme_bw()+NoLegend()+theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
                                                                            plot.title=element_blank(),axis.title =element_blank(),axis.text.y=element_blank(),
                                                                            axis.text.x=element_blank(), axis.text = element_blank(),legend.text=element_blank(),
                                                                            legend.title = element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),panel.border= element_blank())

ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_monocle3_Treat_split_no_text.png",final_anno_plot2_split1,width = 12, height =6)

#plot for density along persudotime
metadata<-as.data.frame(colData(cds_final))
pseudotime <- pseudotime(cds_final, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(metadata)]
metadata$pseudotime <- pseudotime
head(metadata)
cell_distribution<-ggplot(metadata,aes(x=pseudotime,colour=Treat))+geom_density()+scale_color_manual(values = ppCor[c(3,5)]) + theme_classic()   +facet_wrap(~re_annotation_TFac, ncol = 2)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_subcyte_cell_distribution.pdf",cell_distribution,width = 6, height =10)
cell_distribution2<-ggplot(metadata,aes(x=pseudotime,colour=Treat))+geom_density()+scale_color_manual(values = ppCor[c(3,5)]) + theme_classic()   +facet_wrap(~re_annotation_TFac, ncol = 4)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_subcyte_cell_distribution2.pdf",cell_distribution2,width = 12, height =7)

##plot for module genes
gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_gene_module_along_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
head(gene_module)
all_module_expression<-plot_cells(cds_final, genes=gene_module ,label_cell_groups=FALSE,show_trajectory_graph=FALSE)
all_module_expression2<-all_module_expression + scale_color_gradient2(low = "blue", mid = "gray", high =  "red")
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/final_Trophoblast_all_module_expression2.png",all_module_expression2,width = 20, height =20,limitsize = FALSE)
plot_cells(cds_final, genes=gene_module %>% filter(module %in% c(27, 10, 7, 30)),label_cell_groups=FALSE,show_trajectory_graph=FALSE)

#######################################注意线
###analysis for EVTs branch
##heatmap plot for persudotime related genes
cds_sub_EVT <- choose_graph_segments(cds_final)
cds_branch<-cds_final[,rownames(colData(cds_final)) %in% rownames(colData(cds_sub_EVT))]
EVT_branch_plot<-plot_cells(cds_branch,  reduction_method="UMAP",color_cells_by = "pseudotime", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_EVT_branch_umap.pdf",EVT_branch_plot,width = 8, height =7)
ggsave(file="/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_EVT_branch_umap.png",EVT_branch_plot,width = 8, height =7)

plot_cells(cds_branch, color_cells_by = "re_annotation_TFac", show_trajectory_graph =FALSE, label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)+scale_colour_manual(values=ppCor)
saveRDS(cds_branch,file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_EVT_selected_branch_monocle3_object.rds")

###analysis for EVTs branch
cds_branch<-readRDS(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_EVT_selected_branch_monocle3_object.rds")
pseu_time<-data.frame(pseudotime(cds_branch))
colnames(pseu_time)<-"pseudotime"
head(pseu_time)
#pseu_time$number<-ifelse(rownames(pseu_time) %in% rownames(colData(cds_branch)),1,-1)
pseu_time$cell_name<-rownames(pseu_time)
pseu_time<-pseu_time[order(pseu_time$pseudotime,decreasing = F),]
head(pseu_time);tail(pseu_time)

gene_module<-read.table(file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/final_Trophoblast_gene_module_along_pseudotime.txt", sep = '\t', row.names = 1, header = TRUE)
head(gene_module)
gene_module$id<-as.character(gene_module$id)
gene_module<-gene_module[order(gene_module$module,decreasing = F),]
gap_number0<-as.numeric(table(gene_module$module))
gap_number<-gap_number0[-length(gap_number0)]
head(gene_module)

pt.matrix <- exprs(cds_branch)[match(gene_module$id,rownames(rowData(cds_branch))),]
pt.matrix[1:4,1:4];dim(pt.matrix)#1296 8232
#reorder cell by persudotime
pt.matrix<-pt.matrix[,pseu_time$cell_name]
pt.matrix[1:4,1:4];dim(pt.matrix)#1296 8232
##recalculation for smooth expression along psudotime
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
#rownames(pt.matrix) <- pr_DEGs_branch
pt.matrix[1:4,1:4];dim(pt.matrix)#1296 8232
#plot gene trend in pheatmap
plot_matrix<-pt.matrix
colnames(plot_matrix) = 1:ncol(plot_matrix)

htkm <- Heatmap(plot_matrix,name= "z-score",#border = TRUE,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(paste0("Module_",1:length(gap_number0)),gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_selected_direction_EVTs_all_module.pdf",height=15,width=8)
print(htkm)
dev.off()

gene_module$type<-"Early_time"
gene_module[which(gene_module$module %in% c(10,12)),]$type<-"Later_time"
gene_module[which(gene_module$module %in% c(11,4,7)),]$type<-"Mid_low"
#gene_module[which(gene_module$module %in% c(4,7)),]$type<-"Mid_high"
#gene_module$type<-factor(gene_module$type,levels = c("Early_time","Mid_low","Later_time","Mid_high"))
gene_module$type<-factor(gene_module$type,levels = c("Early_time","Mid_low","Later_time"))
gene_module<-gene_module[order(gene_module$type,decreasing = F),]
gap_number0<-as.numeric(table(gene_module$type))
gap_number<-gap_number0[-length(gap_number0)]
head(gene_module)

plot_matrix2<-plot_matrix[as.character(gene_module$id),]
htkm <- Heatmap(plot_matrix2,name= "z-score",border = TRUE,
                # row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(c("Early_time","Mid_low","Later_time","Mid_high"),gap_number0),
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(c("Early_time","Mid_low","Later_time"),gap_number0),
                cluster_column_slices = FALSE,column_title = "Developmental Gene Expression in direction of EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_rows= FALSE,cluster_row_slices= FALSE,cluster_columns= FALSE)
print(htkm)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_selected_direction_EVTs_rearrange_module.pdf",height=15,width=8)
print(htkm)
dev.off()

#GO BP enrichment
allcluster_BP_GO=data.frame();allcluster_BP_GO_whole=data.frame()

for (cluster_name in levels(gene_module$type)) {
  #cluster_name<-"Mid_low"
  print(paste0("cluster number: ",cluster_name))
  small_gene_group=unique(as.character(gene_module[which(gene_module$type == cluster_name),]$id))
  gene_change=bitr(small_gene_group, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  BP_GO <- enrichGO(gene= unique(gene_change$ENTREZID),OrgDb= org.Hs.eg.db,keyType= 'ENTREZID',ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff= 0.1,readable= TRUE)
  BP_GO2 <- simplify(BP_GO, cutoff=0.7, by="p.adjust", select_fun=min)
  go_res=BP_GO2@result
  go_res_whole=BP_GO@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=cluster_name;allcluster_BP_GO=rbind(allcluster_BP_GO,go_res)}
  if (dim(go_res_whole)[1] != 0) {
    go_res_whole$cluster=cluster_name;allcluster_BP_GO_whole=rbind(allcluster_BP_GO_whole,go_res_whole)}
}
head(allcluster_BP_GO[,c("ID","Description","qvalue","cluster")])
table(allcluster_BP_GO$cluster)
#Early_time Later_time    Mid_low 
#   202         55         72 
write.table(allcluster_BP_GO, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_for_gene_module_along_slected_EVT_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(allcluster_BP_GO_whole, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_whole_for_module_along_slected_EVT_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

#relationship between DEGs and development genes
#selected candiidated group for venn plot
#gene file 1: single cell DEGs
merge_data0 <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_2/RNA_data_DEG_01q_default_merge_new.txt", header=T)
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


#FOR BULK UP/Down_gene
venn <-venn.diagram(list(bulk_Up=bulk_Up_gene,bulk_Down=bulk_Down_gene,develop_gene=gene_module$id),
                    alpha=c(0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","navy","grey"), 
                    cex = 1.5,cat.col=c("red","navy","grey"), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage()
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Venn_development_gene_and_bulk_DEGs_for_Abortion.pdf",width = 6,height = 6)
grid.draw(venn)
dev.off()

table(gene_module$type,gene_module$id %in% c(bulk_Up_gene,bulk_Down_gene))
table(gene_module$type,gene_module$id %in% c(bulk_Up_gene))

#           FALSE TRUE  Up_TRUE
#Early_time   630  305  97
#Mid_low       11   19  19
#Later_time    45   38  37
#Mid_high     161   87  59

subgroup_order2<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3")
merge_data1<-distinct(merge_data0[which(merge_data0$cluster %in% subgroup_order2),c("gene","cluster_trend")])
merge_data1$cluster_trend<-factor(merge_data1$cluster_trend,levels =c(paste0(subgroup_order2,"_Up"),paste0(subgroup_order2,"_Down")))
table(merge_data1$gene %in% gene_module$id,merge_data1$cluster_trend)
##         CTBs_1_Up CTBs_2_Up STBs_1_Up STBs_2_Up STBs_3_Up EVTs_1_Up EVTs_2_Up EVTs_3_Up CTBs_1_Down CTBs_2_Down STBs_1_Down STBs_2_Down STBs_3_Down EVTs_1_Down EVTs_2_Down EVTs_3_Down
##FALSE        80        80       128       164       136        52       209        21          56          65         170          85         136          27         132          67
##TRUE         87       111       107       100       256        70       221        76          52          67         153         105         195          77         292         134

#FOR sc EVT2 and EVT3 UP/Down_gene

venn <-venn.diagram(list(EVTs_2_Up=unique(as.character(merge_data1[which(merge_data1$cluster_trend == "EVTs_2_Up"),]$gene)),
                         EVTs_3_Up=unique(as.character(merge_data1[which(merge_data1$cluster_trend == "EVTs_3_Up"),]$gene)),
                         EVTs_2_Down =unique(as.character(merge_data1[which(merge_data1$cluster_trend == "EVTs_2_Down"),]$gene)),
                         EVTs_3_Down =unique(as.character(merge_data1[which(merge_data1$cluster_trend == "EVTs_3_Down"),]$gene)),
                         develop_genes =unique(gene_module$id)),
                    alpha=c(0.8,0.8,0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c(ppCor[1:5]), 
                    cex = 1.5,cat.col=ppCor[1:5], cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage()
pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/Venn_development_gene_and_sc_EVTs_2_3_DEGs_for_Abortion.pdf",width = 8,height = 8)
grid.draw(venn)
dev.off()

#complexheatmap for Arrest changed_development genes 
merge_data2<-merge_data1[which(merge_data1$cluster_trend %in% c("EVTs_2_Up","EVTs_3_Up","EVTs_2_Down","EVTs_3_Down")),]

#building expression matrix
Abortion_EVTs_2_3_DEGs<-unique(as.character(merge_data1[which(merge_data1$cluster_trend %in% c("EVTs_2_Up","EVTs_3_Up","EVTs_2_Down","EVTs_3_Down")),]$gene))
overlap_genes<-intersect(Abortion_EVTs_2_3_DEGs,unique(gene_module$id))

table(gene_module[which(gene_module$id %in% overlap_genes),]$module)
##1   2   3   4   5   6   7   8   9  10  11  12 
#129  53  84  86  15  85  15   6  16  60   7  21 

row_anno <-data.frame(row_genes=overlap_genes)
row_anno$EVTs_2_Up<-"No";row_anno[which(row_anno$row_genes %in% as.character(merge_data2[which(merge_data2$cluster_trend == "EVTs_2_Up"),]$gene)),]$EVTs_2_Up<-"EVTs_2_Up"
row_anno$EVTs_3_Up<-"No";row_anno[which(row_anno$row_genes %in% as.character(merge_data2[which(merge_data2$cluster_trend == "EVTs_3_Up"),]$gene)),]$EVTs_3_Up<-"EVTs_3_Up"
row_anno$EVTs_2_Down<-"No";row_anno[which(row_anno$row_genes %in% as.character(merge_data2[which(merge_data2$cluster_trend == "EVTs_2_Down"),]$gene)),]$EVTs_2_Down<-"EVTs_2_Down"
row_anno$EVTs_3_Down<-"No";row_anno[which(row_anno$row_genes %in% as.character(merge_data2[which(merge_data2$cluster_trend == "EVTs_3_Down"),]$gene)),]$EVTs_3_Down<-"EVTs_3_Down"
head(row_anno);dim(row_anno)# 577   5

plot_matrix2<-plot_matrix[overlap_genes,]
dim(plot_matrix2)#577 8232
plot_matrix2[1:4,1:4]

#set annotation
EVTs_2_Up<-row_anno$EVTs_2_Up;EVTs_3_Up<-row_anno$EVTs_3_Up;EVTs_2_Down<-row_anno$EVTs_2_Down;EVTs_3_Down<-row_anno$EVTs_3_Down
ann_colors = list(EVTs_2_Up=c(EVTs_2_Up=ppCor[1],No="darkgray"),EVTs_3_Up=c(EVTs_3_Up=ppCor[2],No="darkgray"),
                  EVTs_2_Down=c(EVTs_2_Down=ppCor[3],No="darkgray"),EVTs_3_Down=c(EVTs_3_Down=ppCor[5],No="darkgray"))
row_ha = rowAnnotation(EVTs_2_Up =EVTs_2_Up,EVTs_3_Up =EVTs_3_Up,EVTs_2_Down =EVTs_2_Down,EVTs_3_Down =EVTs_3_Down,col = ann_colors)

#Heatmap(plot_matrix2, name = "mat", left_annotation = row_ha)


dend = as.dendrogram(hclust(dist(plot_matrix2),method="ward.D2"))
#d_num<-6;dend = color_branches(dend, k =d_num)
d_num<-3;dend = color_branches(dend, k =d_num)

htkm <- Heatmap(plot_matrix2,name= "z-score", border = TRUE,
                cluster_column_slices = FALSE,column_title = "Abortion DEG & Developmental Gene Expression in direction of EVTs", 
                col= colorRamp2(seq(from=-5,to=5,length=11),rev(brewer.pal(11, "Spectral"))),
                show_row_names= FALSE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,#cluster_rows= FALSE,
                cluster_rows = dend,
                row_split = d_num,row_gap = unit(rep(2, d_num-1), "mm"),row_dend_width = unit(4, "cm"),
                right_annotation = row_ha)
print(htkm)

pdf("/mnt/data/chenwei/gongchen/manuscript/manu_figure/major_3/pheatmap_overlap_DEGs_selected_direction_EVTs_rearrange_module.pdf",height=9,width=8)
print(htkm)
dev.off()

#exact genes names in each cluster
clusterlist = row_order(htkm)
names(clusterlist)<-1:d_num
htkm_module <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(plot_matrix2[clusterlist[[i]],]),Cluster = paste0("cluster", i),stringsAsFactors = FALSE)
  return(out)}) %>%   do.call(rbind, .)

head(htkm_module);colnames(htkm_module)<-c("row_genes","Cluster")
table(htkm_module$Cluster)
#cluster1 cluster2 cluster3 cluster4 cluster5 cluster6 
#  45      108       58      161      137       68 
#cluster1 cluster2 cluster3 
#  153       58      366 

##add cluster information for genes annotation
row_anno2<-merge(row_anno,htkm_module,by="row_genes")
head(row_anno2)
row_anno2$trend<-ifelse(row_anno2$EVTs_2_Up != "No" & row_anno2$EVTs_3_Up != "No","Both_Up",
                        ifelse(row_anno2$EVTs_2_Down != "No" & row_anno2$EVTs_3_Down != "No","Both_Down",
                               ifelse(row_anno2$EVTs_2_Up != "No" & row_anno2$EVTs_3_Down != "No","Up_Down",
                                      ifelse(row_anno2$EVTs_2_Down != "No" & row_anno2$EVTs_3_Up != "No","Down_Up",
                                             ifelse(row_anno2$EVTs_2_Down != "No","EVTs_2_Down", ifelse(row_anno2$EVTs_2_Up != "No","EVTs_2_Up", ifelse(row_anno2$EVTs_3_Down != "No","EVTs_3_Down","EVTs_3_Up")))))))

row_anno2$type<-paste0(row_anno2$Cluster,"-",row_anno2$trend)
table(row_anno2$type)
head(row_anno2)

##rearange genes by different classes
row_anno2$row_genes<-as.character(row_anno2$row_genes)
row_anno2$type<-factor(row_anno2$type,levels = names(table(row_anno2$type)))
row_anno2<-row_anno2[order(row_anno2$type,decreasing = F),]
head(row_anno2)
write.table(row_anno2, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/DEGs_EVT2_EVT3_overlapped_development_genes_along_sleceted_EVT_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
table(row_anno2$Cluster)
#cluster1 cluster2 cluster3 
# 153       58      366 
row_anno2 <- read.table("/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/DEGs_EVT2_EVT3_overlapped_development_genes_along_sleceted_EVT_branchment_pseudotime.txt", header=T)
#GO BP enrichment1
allcluster_BP_GO=data.frame();allcluster_BP_GO_whole=data.frame()

for (cluster_name in unique(row_anno2$Cluster)) {
  #cluster_name<-"cluster1"
  print(paste0("cluster number: ",cluster_name))
  small_gene_group=unique(as.character(row_anno2[which(row_anno2$Cluster == cluster_name),]$row_genes))
  gene_change=bitr(small_gene_group, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  BP_GO <- enrichGO(gene= unique(gene_change$ENTREZID),OrgDb= org.Hs.eg.db,keyType= 'ENTREZID',ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff= 0.1,readable= TRUE)
  BP_GO2 <- simplify(BP_GO, cutoff=0.7, by="p.adjust", select_fun=min)
  go_res=BP_GO2@result
  go_res_whole=BP_GO@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=cluster_name;allcluster_BP_GO=rbind(allcluster_BP_GO,go_res)}
  if (dim(go_res_whole)[1] != 0) {
    go_res_whole$cluster=cluster_name;allcluster_BP_GO_whole=rbind(allcluster_BP_GO_whole,go_res_whole)}
}
head(allcluster_BP_GO[,c("ID","Description","qvalue","cluster")])
table(allcluster_BP_GO$cluster)
write.table(allcluster_BP_GO, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_for_DEGs_EVT2_EVT3_overlapped_development_genes_for_different_cluster_along_slected_EVT_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(allcluster_BP_GO_whole, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_whole_for_DEGs_EVT2_EVT3_overlapped_development_genes_for_different_cluster_along_slected_EVT_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)

#GO BP enrichment2
allcluster_BP_GO=data.frame();allcluster_BP_GO_whole=data.frame()
cluster_type<-c("cluster1-Both_Up","cluster1-EVTs_2_Up","cluster1-EVTs_3_Up",
                "cluster2-Both_Up","cluster2-EVTs_2_Up","cluster2-EVTs_3_Up","cluster3-Both_Up",   
                "cluster3-EVTs_2_Down","cluster3-EVTs_2_Up","cluster3-EVTs_3_Down","cluster3-EVTs_3_Up", 
                "cluster3-Up_Down","cluster4-Both_Down","cluster4-Down_Up","cluster4-EVTs_2_Down",
                "cluster4-EVTs_3_Down","cluster5-Both_Down","cluster5-EVTs_2_Down","cluster5-EVTs_3_Down",
                "cluster5-Up_Down","cluster6-Both_Down","cluster6-EVTs_2_Down","cluster6-EVTs_2_Up",  
                "cluster6-EVTs_3_Down","cluster6-EVTs_3_Up","cluster6-Up_Down")
#cluster1-Up_Down
for (cluster_name in cluster_type) {
  #cluster_name<-"cluster1-Both_Up"
  print(paste0("cluster number: ",cluster_name))
  small_gene_group=unique(as.character(row_anno2[which(row_anno2$type == cluster_name),]$row_genes))
  gene_change=bitr(small_gene_group, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  BP_GO <- enrichGO(gene= unique(gene_change$ENTREZID),OrgDb= org.Hs.eg.db,keyType= 'ENTREZID',ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff= 0.1,readable= TRUE)
  BP_GO2 <- simplify(BP_GO, cutoff=0.7, by="p.adjust", select_fun=min)
  go_res=BP_GO2@result
  go_res_whole=BP_GO@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=cluster_name;allcluster_BP_GO=rbind(allcluster_BP_GO,go_res)}
  if (dim(go_res_whole)[1] != 0) {
    go_res_whole$cluster=cluster_name;allcluster_BP_GO_whole=rbind(allcluster_BP_GO_whole,go_res_whole)}
}
head(allcluster_BP_GO[,c("ID","Description","qvalue","cluster")])
table(allcluster_BP_GO$cluster)
write.table(allcluster_BP_GO, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_for_DEGs_EVT2_EVT3_overlapped_development_genes_along_slected_EVT_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
write.table(allcluster_BP_GO_whole, file = "/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure_3/Enrichment_GO_BP_whole_for_DEGs_EVT2_EVT3_overlapped_development_genes_along_slected_EVT_branchment_pseudotime.txt", quote = FALSE, sep = '\t', row.names = T, col.names = TRUE)
