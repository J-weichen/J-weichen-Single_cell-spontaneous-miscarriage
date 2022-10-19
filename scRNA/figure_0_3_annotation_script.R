#合并注释大类群
###
target<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0723_integrated_seurat.rds")
#target<-integrated
DimPlot(object = target, reduction = "umap",label = TRUE) + NoLegend() + ggtitle('Cell_default0.6_sctransform') 

DefaultAssay(target) <- "RNA"
# Normalize RNA data for visualization purposes
target <- NormalizeData(target, verbose = FALSE)

##major population 
# "Trophoblast"  "Stromal_cells" "Epithelial_Cell" "Endothelial_Cell" "Erythrocyte" "Myeloid_Cell" "Bs_Cell"  "Ts_Cell"  "NKs_Cell"   
#maker list
#maker list
major_maker<-c("VIM","HLA-B","KRT7","PERP","EPCAM","PAEP","DCN","DLK1","PECAM1","LYVE1","AIF1","CD3D","CD8A","IL7R","NKG7","CD79A","MS4A2","HBA1","PF4","PPBP")
immune_maker<-c("AIF1","CD14","S100A8","FCN1","FCGR3A","CLEC9A","CD1C",
                "CSF1R","CD163","CD209","CD69","LYVE1","LYZ","APOE","MS4A3",
                "MS4A1","CD79A","CD79B","JCHAIN",
                "CD3D","CD8A","IL7R","NKG7","KLRB1","FGFBP2","S100A4","HBB","KIT")

no_immune_maker<-c("VIM","CDH1","EGFR","HLA-G","MMP11","CGA","CYP19A1","CSH2","ERVFRD-1","DCN","DLK1","THY1","COL1A1","LAMA2","TIMP1",
                   "MYH11","RGS5","PRL","IGFBP1","APOD","COL6A2","PECAM1","LYVE1","EGFL7","PPBP","PF4","EPCAM","PAEP")

HBs_maker<-c("LYVE1","AIF1")
troblast_maker<- c("KRT7","VIM")
global_immune_maker <-c("PTPRC","CD3E")
tissue_resident_markers<-c("CD69","ITGA1","CD9")
proliferate_markers<-c("MKI67","TOP2A","TK1","PCNA")
Trophoblast_maker<-c("HLA-B","VIM","KRT7","GATA3","PERP","XAGE2","CDH1","EGFR","PAGE4","HLA-G","PAPPA2","CGA","CSH2","ERVFRD-1")
maker_CBTs<- c("CDH1","EGFR","PAGE4","PEG10")
maker_EVTs<-c("HLA-G","MMP2","PAPPA2","DIO2","MMP11","FLT1","ITGA5","ADAM12","MCAM")
maker_STBs<-c("CSH1","CGA","CSH2","CSHL1","GH2","PSG2","HOPX","TFAP2A","CYP19A1","ERVFRD-1")

global_TNK_maker<-c("PTPRC","HBA1","MS4A2","CD34","CD69","CD3D","CD8A","CD4","KLRB1","KLRF1","IL7R","CCR7","SELL","FOXP3",
                    "NKG7","NCAM1", "FGFBP2","FCGR3A","GNLY","GZMB","GZMK","TRDV2","TRGV9","MKI67","TOP2A")
effect_TNK_maker<-c("CD3D","CD4","CD8A","KLRB1","NKG7","NCAM1","FGFBP2","FCGR3A","CXCR4","XCL1","PRF1","GNLY","GZMA","GZMB","GZMK")
naive_TNK_maker<-c("MS4A2","PTPRC","CD69","CD3D","CD4","CD8A","IL7R","CCR7","SELL","FOXP3", "KLRF1","FGFBP2","NKG7","GNLY","GZMB","GZMK","GZMA","NCAM1","TRDC","TRGC2")
maker_Ts<-c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B")
maker_TCR<-c("TRAC","TRBC1","TRBC2","TRDC","TRGC1","TRGC2")
maker_TRV<-c("TRAV1-2","TRDV1","TRDV2","TRDV3","TRGV3","TRGV4","TRGV9")
Treg_maker<-c("FOXP3","TNFRSF18","IL2RA","CTLA4","TNFRSF9","TIGIT","IL10","IL23R","IL17A","RORA","CCR6")
Bs_maker<-c("MS4A1","CD79A","CD79B","CD19","IGHM", "IGHD","IGHA1","JCHAIN","CD69")
maker_Erythrocyte<-c("HBA1","HBA2","HBG1","HBG2","HBA1","HBB","KIT","ALAS2","SNCA","GYPA","GYPB","GATA1","TK1","TOP2A", "MKI67")
maker_immature_erythrocytic<-c("CD34","GATA2", "GATA1", "TFRC","CD47","HBA1")
maker_Mast<-c("KIT","MS4A2","CPA3","IL7R")
maker_Endothelial<- c("CD34","PECAM1","CLDN5","PLVAP","CDH5","ICAM1","LYVE1","CD9","CCL21","EGFL7")
maker_Stromal_cells<-c("DCN","COL6A2","COL3A1","THY1","DLK1")
maker_Fibroblasts<-c("COL1A1","COL3A1","COL1A2","LAMA2","LAMC3","TIMP1","APOD","PCOLCE")
maker_Perivascular_cells<-c("ACTA2","MYH11","RGS5","NDUFA4L2")
maker_Fibroblasts2<- c("PRL","IGFBP1","IGFBP2","IGFBP6","CYP11A1","APOA1","CHI3L2","SERPINA3","IL1B","PROK1")

#plot makers  For major cell population
plot_raw_cluster<- FeaturePlot(object = target, features = major_maker,cols= c("grey", "red"),ncol=4)
#For no_immunine Cells 5*6
plot_no_immune_maker<-FeaturePlot(object = target, features = no_immune_maker,cols= c("grey", "red"),ncol=5)
plot_HBs_maker<-FeaturePlot(object = target, features =HBs_maker,cols= c("grey", "red"))
#For immunine Cells 6*7
plot_immune_maker<-FeaturePlot(object = target,features = immune_maker, pt.size = 0.2,cols= c("grey", "red"),ncol=6)

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/major_maker.png", plot_raw_cluster,width=20, height=20)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/no_immune_maker.png", plot_no_immune_maker,width=25, height=30)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/immune_maker.png", plot_immune_maker,width=30, height=35)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/HBs_maker.png", plot_HBs_maker,width=10, height=5)

#troblast and no troblast
plot_maker0<-FeaturePlot(object = target, features = troblast_maker,cols= c("grey", "red"))
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/troblast_maker_maker.png", plot_maker0,width=10, height=5)
#global_immune maker 
plot_global_immune_maker<-FeaturePlot(object =  target, features =global_immune_maker,cols= c("grey", "red"))
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/global_immune_maker.png", plot_global_immune_maker,width=10, height=5)
#tissue resident markers 
plot_maker2<-FeaturePlot(object =  target, features = tissue_resident_markers,cols= c("grey", "red"),ncol=3)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/tissue_resident_markers.png", plot_maker2,width=15, height=5)
#Mitotic or  proliferating subpopulation: MKI67, TOP2A ,TK1;
plot_maker3<-FeaturePlot(object = target, features = proliferate_markers,cols= c("grey", "red"),ncol=3)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/proliferate_markers.png", plot_maker3,width=15, height=5)

#plot for Trophoblast: 5 X 3
plot_maker_Trophoblast<-FeaturePlot(object = target, features = Trophoblast_maker,ncol=5,cols= c("grey", "red"))
#FeaturePlot(object = target, features = c("EPCAM","PAEP","CAPS"),cols= c("grey", "red"),ncol=2)
#Cytotrophoblast（CBTs）2X 2
plot_maker_CBTs<-FeaturePlot(object = target, features = maker_CBTs,ncol=2,cols= c("grey", "red"))
#Extravilous trophoblast（EVTs）3x3
plot_maker_EVTs<-FeaturePlot(object = target, features =maker_EVTs,ncol=3,cols= c("grey", "red"))
#Syncytiotrophoblast (STBs) 4 x 3
plot_maker_STBs<-FeaturePlot(object = target, features = maker_STBs,cols= c("grey", "red"))
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/whole_trophoblast_markers.png", plot_maker_Trophoblast,width=25, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/CBTs_markers.png", plot_maker_CBTs,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/EVTs_markers.png", plot_maker_EVTs,width=15, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/STBs_markers.png", plot_maker_STBs,width=20, height=15)


#for TNK cells
plot_maker_T_cell1<-FeaturePlot(object = target, features =global_TNK_maker,cols= c("grey", "red"),ncol= 5)
plot_maker_T_cell2<-FeaturePlot(object = target, features =effect_TNK_maker,cols= c("grey", "red"),ncol= 5)
plot_maker_T_cell3<-FeaturePlot(object = target, features = naive_TNK_maker,cols= c("grey", "red"),ncol= 5)

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/global_TNK_maker.png", plot_maker_T_cell1,width=25, height=25)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/effect_TNK_maker.png", plot_maker_T_cell2,width=25, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/naive_TNK_maker.png", plot_maker_T_cell3,width=25, height=20)

#T cell maker
plot_maker_T_cell<-FeaturePlot(object = target, features =maker_Ts,cols= c("grey", "red"),ncol=3)
plot_maker_TCR_cell<-FeaturePlot(object = target, features =maker_TCR,cols= c("grey", "red"),ncol=3)
plot_maker_TRV_cell<-FeaturePlot(object = target, features = maker_TRV,cols= c("grey", "red"),ncol=3)
#  The following requested variables were not found: TRDV1, TRDV3
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/Tcell_markers.png", plot_maker_T_cell,width=15, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/TCR_markers.png", plot_maker_TCR_cell,width=15, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/TRV_markers.png", plot_maker_TRV_cell,width=15, height=10)

#######roughly identify for major population############
##major population 
# "Trophoblast"  "Stromal_cells" "Epithelial_Cell" "Endothelial_Cell" "Erythrocyte" "Myeloid_Cell" "Bs_Cell"  "Ts_Cell"  "NKs_Cell" 

DimPlot(object = target, group.by ="integrated_snn_res.0.6", label = TRUE) + NoLegend()# + ggtitle('Cell_default0.6_sctransform') 
Idents(object = target) <- "integrated_snn_res.0.6"
target$raw_cluster <- as.character(Idents(target))
target$raw_cluster_brif <- as.character(Idents(target))

#Erythrocyte
sub_cell<-subset(x = target,idents=c("15"))
target$raw_cluster[Cells(sub_cell)] <- "Erythrocyte_HBA1_pos"
target$raw_cluster_brif[Cells(sub_cell)] <- "Ery"
#Myeloid_Cell
sub_cell<-subset(x = target,idents=c("0","5","8","13")) 
target$raw_cluster[Cells(sub_cell)] <- "Myeloid_Cell_AIF1_pos"
target$raw_cluster_brif[Cells(sub_cell)] <- "MyCs"

# TNKB_Cells:  CD79A+ CD3D+CD8A+
sub_cell<-subset(x = target,idents=c("14"))
target$raw_cluster[Cells(sub_cell)] <- "NKsTsBs_Cells_NCAM1_CD3D_CD79A"
target$raw_cluster_brif[Cells(sub_cell)] <- "NKsTsBs"

# Stromal_cells_DCN_pos
sub_cell<-subset(x = target,idents=c("6","18"))
target$raw_cluster[Cells(sub_cell)] <- "Stromal_cells_DCN_pos"
target$raw_cluster_brif[Cells(sub_cell)] <- "STCs"

#Endothelial_Cell_PECAM1
sub_cell<-subset(x = target,idents=c("19"))
target$raw_cluster[Cells(sub_cell)] <- "Endothelial_Cell_PECAM1_pos"
target$raw_cluster_brif[Cells(sub_cell)] <- "Endo"

##### Trophoblast_KRT7_pos
sub_cell<-subset(x = target,idents=c("1","2","3","4","7","9","10","11","12","16","17"))
target$raw_cluster[Cells(sub_cell)] <- "Trophoblast_KRT7_pos"
target$raw_cluster_brif[Cells(sub_cell)] <- "Troph"
DimPlot(object = target, group.by ="raw_cluster_brif",label = TRUE) + NoLegend() 

plot_for_raw_cluster_brif<-DimPlot(object = target, group.by ="raw_cluster_brif",label = TRUE,cols = ppCor_all2)
plot_for_raw_cluster_brif_split<-DimPlot(object = target, reduction = "umap", split.by = "Treat",group.by = "raw_cluster_brif",label = TRUE,cols = ppCor_all2,ncol = 2)

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_cluster_brif.png", plot_for_raw_cluster_brif,width=12, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_cluster_brif_split.png", plot_for_raw_cluster_brif_split,width=21, height=10)

#recluster for each Cells groups
table(target$raw_cluster)
#Endothelial_Cell_PECAM1_pos   Erythrocyte_HBA1_pos         Myeloid_Cell_AIF1_pos NKsTsBs_Cells_NCAM1_CD3D_CD79A          Stromal_cells_DCN_pos 
#         131                           1166                          13129                           1296                           3885 
#Trophoblast_KRT7_pos 
#       30381
raw_population<-names(table(target$raw_cluster)) 
Cell_submian_list<-list()
raw_population
for (population in raw_population) {
    #population<-"Erythrocyte_HBA1_pos"  ##test line
    print(as.character(population))
    sub_cell<-subset(x = target, subset = raw_cluster == population)
    DefaultAssay(sub_cell) <- "RNA"    ## very important
    sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
    sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
    ElbowPlot(object = sub_cell)
    sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
    sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
    
    #PCs determination
    ElbowPlot(object = sub_cell)
    sub_cell <- RunUMAP(object = sub_cell, dims = 1:10,verbose = FALSE)
    sub_cell <- FindNeighbors(object = sub_cell, dims = 1:10,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
    A1<-DimPlot(sub_cell,label = TRUE)
    
    sub_cell <- RunUMAP(object = sub_cell, dims = 1:16,verbose = FALSE)
    sub_cell <- FindNeighbors(object = sub_cell, dims = 1:16,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
    A2<-DimPlot(sub_cell,label = TRUE)
    
    sub_cell <- RunUMAP(object = sub_cell, dims = 1:20,verbose = FALSE)
    sub_cell <- FindNeighbors(object = sub_cell, dims = 1:20,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
    A3<-DimPlot(sub_cell,label = TRUE)
    
    sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
    sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
    A4<-DimPlot(sub_cell,label = TRUE)
    plot_for_dif_dims<-CombinePlots(plots = list(A1,A2,A3,A4),ncol = 2,legend = NULL)
    ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/dif_dim_for_",population,".png"), plot_for_dif_dims,width=24, height=20)
           
    #determination of resolution
    sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
    sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
    
    sub_cell <- FindClusters(object = sub_cell,resolution = 0.4,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,resolution = 0.6,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,resolution = 0.8,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,resolution = 1,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,resolution = 1.2,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,resolution =1.5,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,resolution = 2,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,resolution = 2.5,verbose = FALSE)
    sub_cell <- FindClusters(object = sub_cell,resolution =3,verbose = FALSE)
    
    B1<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.4",label = TRUE) + NoLegend() +ggtitle('Cell_default0.4_sctransform')
    B2<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() +ggtitle("Cell_default0.6_sctransform")
    B3<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.8",label = TRUE) + NoLegend() +ggtitle("Cell_default0.8_sctransform")
    B4<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1",label = TRUE) + NoLegend() +ggtitle("Cell_default1_sctransform")
    B5<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1.2",label = TRUE) + NoLegend() +ggtitle("Cell_default1.2_sctransform")
    B6<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1.5",label = TRUE) + NoLegend() +ggtitle("Cell_default1.5_sctransform")
    B7<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.2",label = TRUE) + NoLegend() +ggtitle("Cell_default2_sctransform")
    B8<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.2.5",label = TRUE) + NoLegend() +ggtitle("Cell_default1.5_sctransform")
    B9<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.3",label = TRUE) + NoLegend() +ggtitle("Cell_default3_sctransform")
    
    plot_for_dif_res<-CombinePlots(plots = list(B1,B2,B3,B4,B5,B6,B7,B8,B9),ncol = 3,legend = NULL)
    ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/dif_res_dim25_for",population,".png"), plot_for_dif_res,width=30, height=30)
           
    DefaultAssay(sub_cell) <- "RNA"
    sub_cell <- NormalizeData(sub_cell, verbose = FALSE)
    plot_major_maker_subcluster<- FeaturePlot(object = sub_cell, features = major_maker,cols= c("grey", "red"),ncol=4)
    ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/major_maker_for_",population,".png"), plot_major_maker_subcluster,width=20, height=20)
    
    P1<-DimPlot(object = sub_cell, reduction = "umap", group.by = "Phase",cols = ppCor)
    P2<-DimPlot(object = sub_cell, reduction = "umap", group.by = "Treat",cols =ppCor_all)
    P3<-DimPlot(object = sub_cell, reduction = "umap",group.by = "sample",cols = ppCor)
    P4<-DimPlot(object = sub_cell, reduction = "umap",group.by = "sample_code",cols = ppCor)
    plot_group<-CombinePlots(plots = list(P1,P2,P3,P4),ncol = 2,legend = NULL)
    ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/anno_sample_group_for_",population,".png"), plot_group,width=20, height=20)
    plot_group_split<-DimPlot(object = sub_cell, reduction = "umap", split.by="Treat",group.by = "Treat",cols = ppCor,ncol=2)
    ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/treat_split_for_",population,".png"), plot_group_split,width=20, height=10)
                         
    Cell_submian_list<-c(Cell_submian_list,list(sub_cell))
  }
names(Cell_submian_list)<-raw_population
saveRDS(Cell_submian_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_population_list.rds")


#Cell type determination
Cell_submian_list<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_population_list.rds")
names(Cell_submian_list)
#Erythrocyte_HBA1_pos
target_new<-Cell_submian_list[["Erythrocyte_HBA1_pos_VIM_neg"]]
#target_new<-Cell_submian_list[["Erythrocyte_HBA1_pos"]]
table(target_new$raw_cluster_brif)#Ery 1166
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 
target_new$sub_cluster_brief <- as.character(target_new$raw_cluster_brif)
table(target_new$sub_cluster_brief )#Ery 1166 

Idents(object = target_new) <- "SCT_snn_res.2"
DimPlot(object = target_new, label = TRUE) + NoLegend() 

# Trophoblast_PERP_pos
sub_cell<-subset(x = target_new,idents=c("5","16","9"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "Troph(Er)"
#Hofbauer_Cell_AIF1_pos_LYVE1_pos
sub_cell<-subset(x = target_new,idents=c("13"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "HCs(Er)"
DimPlot(object = target_new, group.by = "sub_cluster_brief",cols = ppCor_all2,label = TRUE)  

table(target_new$sub_cluster_brief)
# Ery   HCs(Er) Troph(Er) 
#968        37       161 
DimPlot(object = target_new, reduction = "umap",group.by = "sub_cluster_brief",cols = ppCor)
Erythrocyte <-target_new

#plot
P1<-DimPlot(object = Erythrocyte, reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = Erythrocyte, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P3<-DimPlot(object = Erythrocyte, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P4<-DimPlot(object = Erythrocyte, reduction = "umap",group.by = "sample",cols = ppCor)
P5<-DimPlot(object = Erythrocyte, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_Ery<-CombinePlots(plots = list(P1,P2,P3,P4,P5),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/subset_type_for_Erythrocyte.png", type_plot_Ery,width=16, height=10)
split_plot_Ery<-DimPlot(object = Erythrocyte, reduction = "umap", split.by="Treat",group.by = "sub_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/plit_subset_type_for_Erythrocyte.png", split_plot_Ery,width=13, height=6)


#2.Endothelial_Cell_PECAM1_pos
names(Cell_submian_list)
target_new<-Cell_submian_list[["Endothelial_Cell_PECAM1_pos"]]
table(target_new$raw_cluster)#Endo 131 
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 
#Cell type determination
target_new$sub_cluster_brief <- as.character(target_new$raw_cluster_brif)
table(target_new$sub_cluster_brief)
#Endo 
#131 
DimPlot(object = target_new, reduction = "umap",group.by = "sub_cluster_brief",cols = ppCor)
Endothelial_cells <-target_new

#plot for Endothelial_cells
P1<-DimPlot(object = Endothelial_cells, reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = Endothelial_cells, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P3<-DimPlot(object = Endothelial_cells, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P4<-DimPlot(object = Endothelial_cells, reduction = "umap",group.by = "sample",cols = ppCor)
P5<-DimPlot(object = Endothelial_cells, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_Endo<-CombinePlots(plots = list(P1,P2,P3,P4,P5),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/subset_type_for_Endothelial_cells.png", type_plot_Endo,width=15, height=10)
split_plot_Endo<-DimPlot(object = Endothelial_cells, reduction = "umap", split.by="Treat",group.by = "sub_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/split_subset_type_for_Endothelial_cells.png", split_plot_Endo,width=13, height=6)

#3.Stromal_cells_DCN_pos   
names(Cell_submian_list)
#[1] "Endothelial_Cell_PECAM1_pos"    "Erythrocyte_HBA1_pos_VIM_neg"   "Myeloid_Cell_AIF1_pos"         
#[4] "NKsTsBs_Cells_NCAM1_CD3D_CD79A" "Stromal_cells_DCN_pos"          "Trophoblast_KRT7_pos"  
target_new<-Cell_submian_list[["Stromal_cells_DCN_pos"]]
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 
DimPlot(object = target_new, reduction = "umap",label = TRUE) + NoLegend() 

#Cell type determination
target_new$sub_cluster_brief <- as.character(target_new$raw_cluster_brif)
table(target_new$sub_cluster_brief)#STCs 3885 
Idents(object = target_new) <- "SCT_snn_res.1"
DimPlot(object = target_new, label = TRUE) + NoLegend() 
# Trophoblast_PERP_pos
sub_cell<-subset(x = target_new,idents=c("2","15","13"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "Troph(SCTs)"
# Myeloid_Cell_AIF1_pos
sub_cell<-subset(x = target_new,idents=c("16"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "MyCs(SCTs)"
table(target_new$sub_cluster_brief)
#  MyCs(SCTs)        STCs Troph(SCTs) 
#         20        3388         477 
DimPlot(object = target_new, reduction = "umap",group.by = "sub_cluster_brief",cols = ppCor)
Stromal_cells <-target_new

P1<-DimPlot(object = Stromal_cells, reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = Stromal_cells, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P3<-DimPlot(object = Stromal_cells, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P4<-DimPlot(object = Stromal_cells, reduction = "umap",group.by = "sample",cols = ppCor)
P5<-DimPlot(object = Stromal_cells, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_SCTs<-CombinePlots(plots = list(P1,P2,P3,P4,P5),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/subset_type_for_Stromal_cells.png", type_plot_FBs,width=15, height=10)
split_plot_SCTs<-DimPlot(object = Stromal_cells, reduction = "umap", split.by="Treat",group.by = "sub_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/plit_subset_type_for_Stromal_cells.png", split_plot_SCTs,width=13, height=6)

#4.NKsTsBs_Cells_NCAM1_CD3D_CD79A
names(Cell_submian_list)
target_new<-Cell_submian_list[["NKsTsBs_Cells_NCAM1_CD3D_CD79A"]]
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 
target_new$sub_cluster_brief <- as.character(target_new$raw_cluster_brif)
table(target_new$sub_cluster_brief)#NKsTsBs 1296 
Idents(object = target_new) <- "SCT_snn_res.2.5"
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.2.5",label = TRUE) + NoLegend() 
#immature erythrocytic cells:CD34,GATA2
FeaturePlot(object = target_new,features = c("CD34","GATA2", "GATA1", "TFRC","CD47","HBA1"),cols= c("grey", "red"))
#Mast_Cells :MS4A2
FeaturePlot(object = target_new, features = c("KIT","MS4A2","CPA3","IL7R"),cols= c("grey", "red"),ncol=2)
#B_Cell_CD79A
sub_cell<-subset(x = target_new,idents=c("19"))
target_new$raw_cluster[Cells(sub_cell)] <- "B_Cell_CD79A"
target_new$raw_cluster_brif[Cells(sub_cell)] <- "Bs"
target_new$sub_cluster_brief[Cells(sub_cell)] <- "Bs"
# Trophoblast_PERP_pos
sub_cell<-subset(x = target_new,idents=c("6","7","21","17"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "Troph(TsNKs)"
# 3.Stromal_cells_DCN_pos_DLK1_pos
sub_cell<-subset(x = target_new,idents=c("22"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "STCs(TsNKs)"
#Myeloid cells 
sub_cell<-subset(x = target_new,idents=c("5"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "MyCs(TsNKs)"
#Immature_Erythrocytic_cells_AIF1_pos_CD34_pos_GATA2_pos
#sub_cell<-subset(x = target_new,idents=c("5"))
#target_new$sub_cluster_brief[Cells(sub_cell)] <- "Ery_im"
table(target_new$sub_cluster_brief)
#   Bs  MyCs(TsNKs)      NKsTsBs  STCs(TsNKs) Troph(TsNKs) 
#   32           76          971           21          196 
DimPlot(object = target_new, group.by = "sub_cluster_brief",cols = ppCor_all2,label = TRUE)
NKsTsBs_Cells<-target_new

P1<-DimPlot(object = NKsTsBs_Cells, reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = NKsTsBs_Cells, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P3<-DimPlot(object = NKsTsBs_Cells, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P4<-DimPlot(object = NKsTsBs_Cells, reduction = "umap",group.by = "sample",cols = ppCor)
P5<-DimPlot(object = NKsTsBs_Cells, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_NKsTsBs<-CombinePlots(plots = list(P1,P2,P3,P4,P5),ncol = 3,legend = NULL)

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/subset_type_for_NKsTsBs_Cells.png", type_plot_NKsTsBs,width=15, height=10)
split_plot_NKsTsBs<-DimPlot(object = NKsTsBs_Cells, reduction = "umap", split.by="Treat",group.by = "sub_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/plit_subset_type_for_NKsTsBs_Cells.png", split_plot_NKsTsBs,width=13, height=6)

#5.Trophoblast_KRT7
names(Cell_submian_list)
target_new<-Cell_submian_list[["Trophoblast_KRT7_pos"]]
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 
#Cell type determination
target_new$sub_cluster_brief <- as.character(target_new$raw_cluster_brif)
table(target_new$sub_cluster_brief)#Troph 30381 
Idents(object = target_new) <- "SCT_snn_res.1.5"
DimPlot(object = target_new, label = TRUE) + NoLegend()

#Hofbauer_Cell_AIF1_pos_LYVE1_pos
sub_cell<-subset(x = target_new,idents=c("20","30"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "MyCs(Troph)"
#dNKsTB_Cells_CD3D_NKG7_CD79A
sub_cell<-subset(x = target_new,idents=c("32"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "TsNKsBs(Troph)"
#Epithelial cells:PAEP 
sub_cell<-subset(x = target_new,idents=c("39"))  
target_new$sub_cluster_brief[Cells(sub_cell)] <- "Epi(Troph)"
#Stromal_cells_DCN_pos_DLK1_pos
sub_cell<-subset(x = target_new,idents=c("25","27","40"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "STCs(Troph)"
#Erythrocyte_HBA1_pos
sub_cell<-subset(x = target_new,idents=c("38"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "Ery(Troph)"

table(target_new$sub_cluster_brief)
# Epi(Troph)     Ery(Troph)    MyCs(Troph)    STCs(Troph)          Troph TsNKsBs(Troph) 
#   54             79            835            810          28344            259 
DimPlot(object = target_new, group.by = "sub_cluster_brief",cols = ppCor_all2,label = TRUE) 
Trophoblast<-target_new

P1<-DimPlot(object = Trophoblast, reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = Trophoblast, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P3<-DimPlot(object = Trophoblast, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P4<-DimPlot(object = Trophoblast, reduction = "umap",group.by = "sample",cols = ppCor)
P5<-DimPlot(object = Trophoblast, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_Troph<-CombinePlots(plots = list(P1,P2,P3,P4,P5),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/subset_type_for_Trophoblast.png", type_plot_Troph,width=15, height=10)
split_plot_Troph<-DimPlot(object = Trophoblast, reduction = "umap", split.by="Treat",group.by = "sub_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/plit_subset_type_for_Trophoblast.png", split_plot_Troph,width=13, height=6)

#6.Myeloid_Cell_AIF1
names(Cell_submian_list)
target_new<-Cell_submian_list[["Myeloid_Cell_AIF1_pos"]]
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 
#Cell type determination
target_new$sub_cluster_brief <- as.character(target_new$raw_cluster_brif)
table(target_new$sub_cluster_brief)#MyCs 13129 

Idents(object = target_new) <- "SCT_snn_res.2"
DimPlot(object = target_new, label = TRUE) + NoLegend() 

# Trophoblast_PERP_pos
sub_cell<-subset(x = target_new,idents=c("1","8","6","9","12","14","19","24","30"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "Troph(MyCs)"
#dNKsTB_Cells_CD3D_NKG7_CD79A
sub_cell<-subset(x = target_new,idents=c("27"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "TNKsBs(MyCs)"
#Epithelial cells:PAEP 
#sub_cell<-subset(x = target_new,idents=c("25"))  
#target_new$sub_cluster_brief[Cells(sub_cell)] <- "wait_Epi(Mye)"
#Stromal_cells_DCN_pos_DLK1_pos
sub_cell<-subset(x = target_new,idents=c("21","22","23"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "SCTs(MyCs)"
#Trophoblast_PERP_pos mixed Myeloid cells
sub_cell<-subset(x = target_new,idents=c("2","7","15","16","17","20","26","28","32"))
target_new$sub_cluster_brief[Cells(sub_cell)] <- "Troph_MyCs_mix"

table(target_new$sub_cluster_brief)
#  MyCs     SCTs(MyCs)   TNKsBs(MyCs) Troph_MyCs_mix    Troph(MyCs) 
#  5345            609            104           2968           4103 
DimPlot(object = target_new, group.by = "sub_cluster_brief",cols = ppCor_all2,label = TRUE) 
Myeloid_cells<-target_new

P1<-DimPlot(object = Myeloid_cells, reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = Myeloid_cells, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P3<-DimPlot(object = Myeloid_cells, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P4<-DimPlot(object = Myeloid_cells, reduction = "umap",group.by = "sample",cols = ppCor)
P5<-DimPlot(object = Myeloid_cells, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_MyCs<-CombinePlots(plots = list(P1,P2,P3,P4,P5),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/subset_type_for_Myeloid_cells.png", type_plot_MyCs,width=15, height=10)
split_plot_MyCs<-DimPlot(object = Myeloid_cells, reduction = "umap", split.by="Treat",group.by = "sub_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/raw_sub_population/split_subset_type_for_Myeloid_cells.png", split_plot_MyCs,width=13, height=6)

####
Cell_submian_list<-c(list(Erythrocyte),list(Endothelial_cells),list(Stromal_cells),list(Trophoblast),
                     list(Myeloid_cells),list(NKsTsBs_Cells))
names(Cell_submian_list)<-c("Erythrocyte","Endothelial_cells","Stromal_cells","Trophoblast","Myeloid_cells","NKsTsBs_Cells")
saveRDS(Cell_submian_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/selected_subpopulation_list.rds")


P1<-DimPlot(object = Cell_submian_list[["Erythrocyte"]], reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Erythrocyte')
P2<-DimPlot(object =Cell_submian_list[["Endothelial_cells"]], reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Endothelial_cells')
P3<-DimPlot(object =Cell_submian_list[["Stromal_cells"]], reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Fibroblasts')
P4<-DimPlot(object =Cell_submian_list[["Trophoblast"]], reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Trophoblast_KRT7')
P5<-DimPlot(object =Cell_submian_list[["Myeloid_cells"]], reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('Myeloid_cells_AIF1')
P6<-DimPlot(object =Cell_submian_list[["NKsTsBs_Cells"]], reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() +ggtitle('BsTsNKs_cells')

type_plot_each_cells<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/selected_type_for_each_major_cell.png", type_plot_each_cells,width=26, height=16)

#next round for fine adjustment of cell population
Cell_submian_list<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/selected_subpopulation_list.rds")

#merge metadata
dim(Cell_submian_list[["Erythrocyte"]]@meta.data);dim(Cell_submian_list[["Endothelial_cells"]]@meta.data);#1166   131
dim(Cell_submian_list[["Stromal_cells"]]@meta.data);dim(Cell_submian_list[["Trophoblast"]]@meta.data);#3885 30381 
dim(Cell_submian_list[["Myeloid_cells"]]@meta.data);dim(Cell_submian_list[["NKsTsBs_Cells"]]@meta.data);#13129  1296
metadata_new<-rbind(Cell_submian_list[["Erythrocyte"]]@meta.data,Cell_submian_list[["Endothelial_cells"]]@meta.data,
                    Cell_submian_list[["Stromal_cells"]]@meta.data,Cell_submian_list[["Trophoblast"]]@meta.data,
                    Cell_submian_list[["Myeloid_cells"]]@meta.data,Cell_submian_list[["NKsTsBs_Cells"]]@meta.data)
dim(metadata_new)#49988    46
head(metadata_new)

table(metadata_new$sub_cluster_brief)
names(table(metadata_new$sub_cluster_brief))
  
metadata_new$major_group_brief<-metadata_new$sub_cluster_brief
metadata_new[which(metadata_new$sub_cluster_brief %in% c("Epi(Troph)","Troph","Troph(Er)","Troph(MyCs)","Troph(SCTs)","Troph(TsNKs)")),]$major_group_brief <- "Trophoblast_mix_Epi"
metadata_new[which(metadata_new$sub_cluster_brief %in% c("HCs(Er)","MyCs","Troph_MyCs_mix","MyCs(SCTs)","MyCs(Troph)","MyCs(TsNKs)")),]$major_group_brief <- "Myeloid_Cell"
metadata_new[which(metadata_new$sub_cluster_brief %in%  c("FBs(Mye)","STCs","SCTs(MyCs)","STCs(Troph)","STCs(TsNKs)")),]$major_group_brief <- "Stromal_cells"
metadata_new[which(metadata_new$sub_cluster_brief %in%  c("Endo")),]$major_group_brief <- "Endothelial_Cell"
metadata_new[which(metadata_new$sub_cluster_brief %in%  c("Bs","NKsTsBs","TNKsBs(MyCs)","TsNKsBs(Troph)")),]$major_group_brief <- "NKsTsBs_Cells"
metadata_new[which(metadata_new$sub_cluster_brief %in%  c("Ery","Ery(Troph)")),]$major_group_brief <- "Erythrocyte"

table(metadata_new$major_group_brief)
#Endothelial_Cell         Erythrocyte        Myeloid_Cell       NKsTsBs_Cells       Stromal_cells Trophoblast_mix_Epi 
#      131                1047                9281                1366                4828               33335 

########replace the old metadata
Metadate_old<-target@meta.data
target2<-target#temp save
dim(metadata_new);dim(Metadate_old)#47 47
target@meta.data<-metadata_new
#saveRDS(target, file = "/mnt/data/chenwei/gongchen/3.seurat_result/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0818_integrated_seurat.rds")
umap_p1<-DimPlot(object = target, reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all,label =F) +ggtitle('sub_cluster_brief') 
umap_p2<-DimPlot(object = target, reduction = "umap", group.by = "major_group_brief",cols =ppCor_all,label =T) +ggtitle('major_group_brief') + NoLegend() 
umap_p3<-DimPlot(object = target, reduction = "umap", group.by = "raw_cluster_brif",cols =ppCor_all,label =T) +ggtitle('raw_cluster_brif') + NoLegend() 
umap_p1+umap_p2+umap_p3


###recalculation for pure subgroups
raw_population<-names(table(target$major_group_brief)) 
Cell_submian_list<-list()
raw_population
for (population in raw_population) {
  #population<-"Erythrocyte"  ##test line
  print(as.character(population))
  sub_cell<-subset(x = target, subset = major_group_brief == population)
  DefaultAssay(sub_cell) <- "RNA"    ## very important
  sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
  sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
  ElbowPlot(object = sub_cell)
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
  
  #PCs determination
  ElbowPlot(object = sub_cell)
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:10,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:10,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A1<-DimPlot(sub_cell,label = TRUE)
  
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:16,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:16,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A2<-DimPlot(sub_cell,label = TRUE)
  
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:20,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:20,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A3<-DimPlot(sub_cell,label = TRUE)
  
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A4<-DimPlot(sub_cell,label = TRUE)
  plot_for_dif_dims<-CombinePlots(plots = list(A1,A2,A3,A4),ncol = 2,legend = NULL)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/dif_dim_for_",population,".png"), plot_for_dif_dims,width=24, height=20)
  
  #determination of resolution
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
  
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.4,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.6,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.8,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 1,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 1.2,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution =1.5,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 2,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 2.5,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution =3,verbose = FALSE)
  
  B1<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.4",label = TRUE) + NoLegend() +ggtitle('Cell_default0.4_sctransform')
  B2<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() +ggtitle("Cell_default0.6_sctransform")
  B3<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.8",label = TRUE) + NoLegend() +ggtitle("Cell_default0.8_sctransform")
  B4<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1",label = TRUE) + NoLegend() +ggtitle("Cell_default1_sctransform")
  B5<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1.2",label = TRUE) + NoLegend() +ggtitle("Cell_default1.2_sctransform")
  B6<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1.5",label = TRUE) + NoLegend() +ggtitle("Cell_default1.5_sctransform")
  B7<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.2",label = TRUE) + NoLegend() +ggtitle("Cell_default2_sctransform")
  B8<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.2.5",label = TRUE) + NoLegend() +ggtitle("Cell_default2.5_sctransform")
  B9<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.3",label = TRUE) + NoLegend() +ggtitle("Cell_default3_sctransform")
  
  plot_for_dif_res<-CombinePlots(plots = list(B1,B2,B3,B4,B5,B6,B7,B8,B9),ncol = 3,legend = NULL)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/dif_res_dim25_for",population,".png"), plot_for_dif_res,width=30, height=30)
  
  DefaultAssay(sub_cell) <- "RNA"
  sub_cell <- NormalizeData(sub_cell, verbose = FALSE)
  plot_major_maker_subcluster<- FeaturePlot(object = sub_cell, features = major_maker,cols= c("grey", "red"),ncol=4)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/major_maker_for_",population,".png"), plot_major_maker_subcluster,width=20, height=20)
  
  #plot
  P1<-DimPlot(object = sub_cell, reduction = "umap", group.by = "sub_cluster_brief",cols =ppCor_all2)
  P2<-DimPlot(object = sub_cell, reduction = "umap", group.by = "major_group_brief",cols =ppCor_all2)
  P3<-DimPlot(object = sub_cell, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
  P4<-DimPlot(object = sub_cell, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
  P5<-DimPlot(object = sub_cell, reduction = "umap",group.by = "sample",cols = ppCor)
  P6<-DimPlot(object = sub_cell, reduction = "umap",group.by = "sample_code",cols = ppCor)
  type_plot<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/anno_sample_group_for_",population,".png"), type_plot,width=26, height=16)
  plot_group_split<-DimPlot(object = sub_cell, reduction = "umap", split.by="Treat",group.by = "Treat",cols = ppCor,ncol=2)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/treat_split_for_",population,".png"), plot_group_split,width=20, height=10)
  
  Cell_submian_list<-c(Cell_submian_list,list(sub_cell))
}
names(Cell_submian_list)<-raw_population
saveRDS(Cell_submian_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/recal_selected_population_list.rds")

##final identification for cells 
#Cell type determination
names(Cell_submian_list)
#"Endothelial_Cell"    "Erythrocyte"         "Myeloid_Cell"        "NKsTsBs_Cells"       "Stromal_cells"       "Trophoblast_mix_Epi"
##1.Erythrocyte
target_new<-Cell_submian_list[["Erythrocyte"]]
table(target_new$major_group_brief)#Erythrocyte 1047 
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 
DimPlot(object = target_new, reduction = "umap",label = TRUE) + NoLegend() 
#plot for Erythrocyte
maker_plot_Ery<-FeaturePlot(object = target_new, features = maker_Erythrocyte,cols= c("grey", "red"),ncol=5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_markers_for_Erythrocyte.png", maker_plot_Ery,width=25, height=15)
#immature erythrocytic cells:CD34,GATA2
maker_plot_imma_ery<-FeaturePlot(object = target_new,features =maker_immature_erythrocytic ,cols= c("grey", "red"))
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_markers_for_maker_Erythrocyte_imma_ery.png", maker_plot_imma_ery,width=10, height=15)

#Mast_Cells :MS4A2
maker_plot_ery_Mast_Cells<-FeaturePlot(object = target_new, features = maker_Mast,cols= c("grey", "red"),ncol=2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_markers_for_maker_Erythrocyte_Mast_Cells.png", maker_plot_ery_Mast_Cells,width=10, height=10)

target_new$prefinal_cluster <- as.character(target_new$major_group_brief)
target_new$prefinal_cluster_brief <- as.character(target_new$major_group_brief)

table(target_new$major_group_brief)#5316 

Idents(object = target_new) <- "SCT_snn_res.0.6"
DimPlot(object = target_new, label = TRUE) + NoLegend() 

sub_cell<-subset(x = target_new,idents=c("5","9","2","8","1","7"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Erythrocyte_HBA1_pos_MKI67_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Ery_prol"

sub_cell<-subset(x = target_new,subset = prefinal_cluster_brief %in% c("Ery_prol"),invert = TRUE)
target_new$prefinal_cluster[Cells(sub_cell)] <- "Erythrocyte_HBA1_pos_MKI67_neg"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Ery_no_prol"
DimPlot(object = target_new, group.by = "prefinal_cluster_brief",cols = ppCor_all2,label = TRUE)  

table(target_new$prefinal_cluster_brief)
#Ery_no_prol    Ery_prol 
#  536         511
DimPlot(object = target_new, reduction = "umap",group.by = "prefinal_cluster_brief",cols = ppCor)
Erythrocyte <-target_new

#plot
P1<-DimPlot(object = Erythrocyte, reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = Erythrocyte, reduction = "umap", group.by = "prefinal_cluster",cols = ppCor[c(2:10)])+NoLegend()
P3<-DimPlot(object = Erythrocyte, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P4<-DimPlot(object = Erythrocyte, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P5<-DimPlot(object = Erythrocyte, reduction = "umap",group.by = "sample",cols = ppCor)
P6<-DimPlot(object = Erythrocyte, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_Ery<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/prefinal_subset_type_for_Erythrocyte.png", type_plot_Ery,width=16, height=10)
split_plot_Ery<-DimPlot(object = Erythrocyte, reduction = "umap", split.by="Treat",group.by = "sub_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/prefinal_plit_subset_type_for_Erythrocyte.png", split_plot_Ery,width=13, height=6)


#2.Endothelial_Cell_PECAM1_pos
names(Cell_submian_list)
target_new<-Cell_submian_list[["Endothelial_Cell"]]

table(target_new$raw_cluster)#Endo 131 
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 
DimPlot(object = target_new, reduction = "umap",label = TRUE) + NoLegend() 
maker_plot_Endo<-FeaturePlot(object = target_new, features =maker_Endothelial,cols= c("grey", "red"),ncol=4)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_markers_for_Endothelial_cells.png", maker_plot_Endo,width=20, height=15)

#Cell type determination
target_new$prefinal_cluster <- as.character(target_new$major_group_brief)
target_new$prefinal_cluster_brief <- as.character(target_new$major_group_brief)
table(target_new$major_group_brief)#Endothelial_Cell 131 

Idents(object = target_new) <- "SCT_snn_res.0.4"
DimPlot(object = target_new, label = TRUE) + NoLegend() 

#LECs
sub_cell<-subset(x = target_new,idents=c("1"))  
target_new$prefinal_cluster[Cells(sub_cell)] <- "Lymphatic_Endothelial_cells_CCL21_pos_cd34_neg"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Endo_1"

sub_cell<-subset(x = target_new,subset = prefinal_cluster_brief %in% c("Endo_1"),invert = TRUE)
target_new$prefinal_cluster[Cells(sub_cell)] <- "Endothelial_cells_CCL21_neg_cd34_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Endo_2"
DimPlot(object = target_new, group.by = "prefinal_cluster_brief",cols = ppCor_all2,label = TRUE) + NoLegend() 

table(target_new$prefinal_cluster_brief)
#Endo_1 Endo_2 
#32     99
DimPlot(object = target_new, reduction = "umap",group.by = "prefinal_cluster_brief",cols = ppCor)
Endothelial_cells <-target_new

#plot for Endothelial_cells
P1<-DimPlot(object = Endothelial_cells, reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = Endothelial_cells, reduction = "umap", group.by = "prefinal_cluster",cols = ppCor[c(2:10)])+NoLegend()
P3<-DimPlot(object = Endothelial_cells, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P4<-DimPlot(object = Endothelial_cells, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P5<-DimPlot(object = Endothelial_cells, reduction = "umap",group.by = "sample",cols = ppCor)
P6<-DimPlot(object = Endothelial_cells, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_Endo<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/prefinal_subset_type_for_Endothelial_cells.png", type_plot_Endo,width=15, height=10)
split_plot_Endo<-DimPlot(object = Endothelial_cells, reduction = "umap", split.by="Treat",group.by = "sub_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/prefinal_split_subset_type_for_Endothelial_cells.png", split_plot_Endo,width=13, height=6)


#3.Stromal_cells_DCN_pos   
names(Cell_submian_list)
# "Endothelial_Cell"    "Erythrocyte"      "Stromal_cells"      "Myeloid_Cell"        "NKsTsBs_Cells" 
#"Trophoblast_mix_Epi"
target_new<-Cell_submian_list[["Stromal_cells"]]
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 
DimPlot(object = target_new, reduction = "umap",label = TRUE) + NoLegend() 


#plot for Trophoblast_KRT7
maker_plot_SCT1<-FeaturePlot(object = target_new, features = Trophoblast_maker,ncol=5,cols= c("grey", "red"))
maker_plot_SCT2<-FeaturePlot(object = target_new, features = c("EPCAM","PAEP","CAPS"),cols= c("grey", "red"),ncol=2)
#Cytotrophoblast（VCTs）
maker_plot_SCT3<-FeaturePlot(object = target_new, features = maker_CBTs,ncol=2,cols= c("grey", "red"))
#Extravilous trophoblast（EVTs）
maker_plot_SCT4<-FeaturePlot(object = target_new, features = maker_EVTs,ncol=3,cols= c("grey", "red"))
#Syncytiotrophoblast (STBs)
maker_plot_SCT5<-FeaturePlot(object = target_new, features = maker_STBs,cols= c("grey", "red"))

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_SCTs_trophoblast1.png", maker_plot_SCT1,width=25, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_SCTs_trophoblast2.png", maker_plot_SCT2,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_SCTs_CBTs.png", maker_plot_SCT3,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_SCTs_EVTs.png", maker_plot_SCT4,width=15, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_SCTs_STBs.png", maker_plot_SCT5,width=20, height=15)

#plot for Fibroblasts
maker_plot_FBs_1<-FeaturePlot(object = target_new, features =maker_Stromal_cells,cols= c("grey", "red"))
maker_plot_FBs_2<-FeaturePlot(object = target_new, features = maker_Fibroblasts,cols= c("grey", "red"))
maker_plot_FBs_3<-FeaturePlot(object = target_new, features = maker_Perivascular_cells,cols= c("grey", "red"),ncol=2)
maker_plot_FBs_4<-FeaturePlot(object = target_new, features = c("CD34","PECAM1","CLDN5","PLVAP","CDH5","ICAM1","LYVE1","CD9","CCL21","EGFL7"),cols= c("grey", "red"),ncol=4)
maker_plot_FBs_5<-FeaturePlot(object = target_new, features = maker_Fibroblasts2,cols= c("grey", "red"),ncol=3)

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_markers_for_Stromal_cells.png", maker_plot_FBs_1,width=10, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_markers_for_Fibroblasts1.png", maker_plot_FBs_2,width=15, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_markers_for_Myofibroblast_Perivascular_cells.png", maker_plot_FBs_3,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_markers_for_Fibroblasts2.png", maker_plot_FBs_4,width=20, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_markers_for_Fibroblasts3.png", maker_plot_FBs_5,width=15, height=20)

#Cell type determination
target_new$prefinal_cluster <- as.character(target_new$major_group_brief)
target_new$prefinal_cluster_brief <- as.character(target_new$major_group_brief)
table(target_new$major_group_brief)#Stromal_cells 4828  
Idents(object = target_new) <- "SCT_snn_res.2.5"
DimPlot(object = target_new, label = TRUE) + NoLegend() 
FeaturePlot(object = target_new, features = proliferate_markers,cols= c("grey", "red"),ncol=3)
# Erythrocyte_HBA1_pos
sub_cell<-subset(x = target_new,idents=c("21"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Erythrocyte_HBA1_pos_MKI67_neg"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Ery_no_prol"

##Fibroblasts_DCN_pos_DLK1_pos_mix_trophoblast_KRT7 :DCN pos DLK1 neg
sub_cell<-subset(x = target_new,idents=c("2","12","13","17","20","23","24","29"))  
target_new$prefinal_cluster[Cells(sub_cell)] <- "Fibroblasts_mix_trophoblast_DCN_pos_DLK1_pos__KRT7"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "FBs_Troph_mix"
#Stromal_cells_DCN_pos_DLK1_neg :DCN pos DLK1 neg
sub_cell<-subset(x = target_new,idents=c("5","8","10","19","22","26","28","30"))  
target_new$prefinal_cluster[Cells(sub_cell)] <- "Stromal_cells_DCN_pos_DLK1_neg"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "STCs"
#Perivascular_cells_RGS5_pos_NDUFA4L2_pos
sub_cell<-subset(x = target_new,idents=c("14","16","27")) 
target_new$prefinal_cluster[Cells(sub_cell)] <- "Perivascular_cells_RGS5_pos_NDUFA4L2_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "PVs"
#Myofibroblast_ACTA2_pos
sub_cell<-subset(x = target_new,idents=c("0","6","15","18"))  
target_new$prefinal_cluster[Cells(sub_cell)] <- "Myofibroblast_ACTA2_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "MFs"
#Fibroblasts_PCOLCE_pos_ACTA2_neg
sub_cell<-subset(x = target_new,idents=c("1","3","4","7","9","11","25"))  
target_new$prefinal_cluster[Cells(sub_cell)] <- "Fibroblasts_PCOLCE_pos_ACTA2_neg_NDUFA4L2_neg"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "FBs"

DimPlot(object = target_new, group.by = "prefinal_cluster_brief",cols = ppCor_all2,label = TRUE) + NoLegend() 
table(target_new$prefinal_cluster_brief)
#Ery_no_prol           FBs FBs_Troph_mix           MFs           PVs          STCs 
#         89          1579           994           818           315          1033 

DimPlot(object = target_new, reduction = "umap",group.by = "prefinal_cluster_brief",cols = ppCor)
Stromal_cells <-target_new

P1<-DimPlot(object = Stromal_cells, reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = Stromal_cells, reduction = "umap", group.by = "prefinal_cluster",cols = ppCor[c(2:10)])+NoLegend()
P3<-DimPlot(object = Stromal_cells, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P4<-DimPlot(object = Stromal_cells, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P5<-DimPlot(object = Stromal_cells, reduction = "umap",group.by = "sample",cols = ppCor)
P6<-DimPlot(object = Stromal_cells, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_SCTs<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/subset_type_for_Stromal_cells.png", type_plot_FBs,width=15, height=10)
split_plot_SCTs<-DimPlot(object = Stromal_cells, reduction = "umap", split.by="Treat",group.by = "prefinal_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/split_subset_type_for_split_plot_Stromal_cells.png", split_plot_SCTs,width=13, height=6)


#4.NKsTsBs_Cells
names(Cell_submian_list)
target_new<-Cell_submian_list[["NKsTsBs_Cells"]]
table(target_new$sub_cluster_brief)
#Bs        NKsTsBs   TNKsBs(MyCs) TsNKsBs(Troph) 
#32            971            104            259 
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 

#plot for Trophoblast_mix_EECs_KRT7
maker_plot_NKsTsBs1<-FeaturePlot(object = target_new, features = Trophoblast_maker,ncol=5,cols= c("grey", "red"))
maker_plot_NKsTsBs2<-FeaturePlot(object = target_new, features = c("EPCAM","PAEP","CAPS"),cols= c("grey", "red"),ncol=2)
#Cytotrophoblast（VCTs）
maker_plot_NKsTsBs3<-FeaturePlot(object = target_new, features = maker_CBTs,ncol=2,cols= c("grey", "red"))
#Extravilous trophoblast（EVTs）
maker_plot_NKsTsBs4<-FeaturePlot(object = target_new, features = maker_EVTs,ncol=3,cols= c("grey", "red"))
#Syncytiotrophoblast (STBs)
maker_plot_NKsTsBs5<-FeaturePlot(object = target_new, features = maker_STBs,cols= c("grey", "red"))

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_NKsTsBs_trophoblast_1.png", maker_plot_NKsTsBs1,width=25, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_NKsTsBs_trophoblast_2.png", maker_plot_NKsTsBs2,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_NKsTsBs_CBTs.png", maker_plot_NKsTsBs3,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_NKsTsBs_EVTs.png", maker_plot_NKsTsBs4,width=15, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_NKsTsBs_STBs.png", maker_plot_NKsTsBs5,width=20, height=15)


#for maker of TNK cells
plot_maker_T_cell1<-FeaturePlot(object = target_new, features = effect_TNK_maker,cols= c("grey", "red"),ncol= 5)
plot_maker_T_cell2<-FeaturePlot(object = target_new, features = naive_TNK_maker,cols= c("grey", "red"),ncol= 5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/Tcell_markers_group1_for_effect_TNK_maker.png", plot_maker_T_cell1,width=25, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/Tcell_markers_group2_for_naive_TNK_maker.png", plot_maker_T_cell2,width=25, height=20)

#plot 
#mucosa-associated lymphoid tissue-derived B(MALT-B): IGHM, IgA, JCHAIN
maker_plot_Bs<-FeaturePlot(object = target_new, features =Bs_maker,cols= c("grey", "red"),ncol=3)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Bs.png", maker_plot_Bs,width=15, height=15)
##plot for global T cells
plot_maker_global_T_cell<-FeaturePlot(object = target_new, features = global_TNK_maker,cols= c("grey", "red"),ncol= 5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_global_TsNKs.png", plot_maker_global_T_cell,width=25, height=25)

#immune maker 
FeaturePlot(object =  target_new, features = c("PTPRC"),cols= c("grey", "red"))
#tissue resident markers 
maker_plot_tissue_resident<-FeaturePlot(object =  target_new, features = tissue_resident_markers,cols= c("grey", "red"),ncol=3)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_tissue_resident_cell.png", maker_plot_tissue_resident,width=15, height=6)

#Mitotic or  proliferating subpopulation: MKI67, TOP2A ,TK1
maker_plot_TNKB_prolif<-FeaturePlot(object = target_new, features = c("CD3D","CD4","CD8A","KLRB1","KLRF1",proliferate_markers),cols= c("grey", "red"),ncol=4)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_TNKB_prolif.png", maker_plot_TNKB_prolif,width=20, height=10)

#T cell maker
plot_maker_T_cell<-FeaturePlot(object = target_new, features = maker_Ts,cols= c("grey", "red"),ncol=3)
plot_maker_TCR_cell<-FeaturePlot(object = target_new, features =maker_TCR,cols= c("grey", "red"),ncol=3)
plot_maker_TRV_cell<-FeaturePlot(object = target_new, features =maker_TRV,cols= c("grey", "red"),ncol=3)
#  The following requested variables were not found: TRDV1, TRDV3
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/Tcell_markers.png", plot_maker_T_cell,width=15, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/TCR_markers.png", plot_maker_TCR_cell,width=15, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/TRV_markers.png", plot_maker_TRV_cell,width=15, height=10)

#immature erythrocytic cells:CD34,GATA2
maker_plot_BsTsNKs_imma_ery<-FeaturePlot(object = target_new,features = c("CD34","GATA2", "GATA1", "TFRC","CD47","HBA1"),cols= c("grey", "red"))
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_BsTsNKs_imma_ery.png", maker_plot_BsTsNKs_imma_ery,width=10, height=15)

#Mast_Cells :MS4A2
maker_plot_BsTsNKs_Mast_Cells<-FeaturePlot(object = target_new, features = c("KIT","MS4A2","CPA3","IL7R"),cols= c("grey", "red"),ncol=2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_BsTsNKs_Mast_Cells.png", maker_plot_BsTsNKs_Mast_Cells,width=10, height=10)

#Natural killer Cells:  KLRB1 也存在一定的极化现象
maker_plot_BsTsNKs_NK_Cells<-FeaturePlot(object = target_new, features = c("KLRD1","KLRB1","KLRF1"),cols= c("grey", "red"),ncol= 3)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_NK_Cells.png", maker_plot_BsTsNKs_NK_Cells,width=15, height=6)

maker_plot_BsTsNKs_TNK_Cells<-FeaturePlot(object = target_new, features = c(effect_TNK_maker,"ITGAE","TBX21","CX3CR1","IL7R","SELL"),cols= c("grey", "red"),ncol= 5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_BsTsNKs_TNK_Cells.png", maker_plot_BsTsNKs_TNK_Cells,width=25, height=20)

#circulating NK (cNK) cell markers FCGR3A, CX3CR1, and the transcription factor T-bet (TBX21)
maker_plot_BsTsNKs_cNK_Cells<-FeaturePlot(object = target_new, features = c("CCR1","CCR7","XCL2","ENTPD1","CYP26A1","B4GALNT1","ANXA1", "ITGB2","ITGAE","FCGR3A", "CX3CR1","TBX21","IFNG"),cols= c("grey", "red"),ncol=5)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_BsTsNKs_cNK_Cells.png", maker_plot_BsTsNKs_cNK_Cells,width=25, height=15)

#ILC1-like CELL or CD56bright NK cells : IL7R+TCF7+SELL+KLRD1+NCAM1+ :  innate lymphocyte cell marker CD127 (also known as IL7R)
#unclassical circulating NK (ucNK)   lacked FCGR3A and FGFBP2 
maker_plot_BsTsNKs_ucNK_Cells<-FeaturePlot(object = target_new, features = c("IL7R","TCF7","SELL","KLRD1","NCAM1","CD160","ITGAE","TBX21"),cols= c("grey", "red"),ncol=4)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_ucNK_Cells.png", maker_plot_BsTsNKs_ucNK_Cells,width=20, height=10)

#blood γδT Cells :CD3D nagetive; KLRB1+NCAM1, CX3CR1+ TRDV2 and TRGV9 pos
maker_plot_BsTsNKs_γδT_Cells<-FeaturePlot(object =  target_new, features = c("KLRG1","KLRD1","KLRB1","KLRF1","TRDV2","TRGV9"),cols= c("grey", "red"),ncol=2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_BsTsNKs_γδT_Cells.png", maker_plot_BsTsNKs_γδT_Cells,width=10, height=15)

##Naïve CD4 CD3D+CD4+CCR7 SELL  naive marker genes: SELL, TCF7,CCR7 and ID3  #Memory CD4+ T: IL7R, S100A4, CCR10
maker_plot_BsTsNKs_CD4T_Cells<-FeaturePlot(object = target_new, features = c("CD3D","CD4","CD8A","IL7R","CCR7","ID3","SELL","TCF7","S100A4", "CCR10"),cols= c("grey", "red"),ncol= 4)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_BsTsNKs_CD4T_Cells.png", maker_plot_BsTsNKs_CD4T_Cells,width=20, height=15)

#Active CD4T: CD3D+CD4+LTBhiGATA3+; cTh1: CD3D+CD4+CXCR3+TBX21+STAT4+IFNG+   ex-cTh17: CD3D+CD4+RORA+IL23R+CCR6+STAT4+EOMES+IFNG+
maker_plot_BsTsNKs_other_CD4T_Cells<-FeaturePlot(object = target_new, features = c("LTB","GATA3","CXCR3","TBX21","STAT1","STAT4","IFNG","RORA","RORC","CCR6","IL23R","EOMES"),cols= c("grey", "red"),ncol= 3)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_other_CD4T_Cells.png", maker_plot_BsTsNKs_other_CD4T_Cells,width=15, height=20)

### Regulatory T cells(Treg): FOXP3, CTLA-4, TNFRSF18; IL2RA (CD25 +)
maker_plot_BsTsNKs_Treg_Cells<-FeaturePlot(object = target_new, features = Treg_maker,cols= c("grey", "red"),ncol=3)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_BsTsNKs_Treg_Cells.png", maker_plot_BsTsNKs_Treg_Cells,width=15, height=20)

#Th17 like cells
maker_plot_BsTsNKs_Th17_like_Cells<-FeaturePlot(object = target_new, features = c("RORC","KIT","IL23R","RORA"),cols= c("grey", "red"),ncol=2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_BsTsNKs_Th17_like_Cells.png", maker_plot_BsTsNKs_Th17_like_Cells,width=20, height=20)

FeaturePlot(object = target_new, features = c("TK1","EGFL7","PAEP","MS4A2"),cols= c("grey", "red"))

target_new$prefinal_cluster <- as.character(target_new$major_group_brief)
target_new$prefinal_cluster_brief <- as.character(target_new$major_group_brief)
table(target_new$major_group_brief)#NKsTsBs_Cells 1366 
#plot for NKsTsBs_Cells
Idents(object = target_new) <- "SCT_snn_res.3"
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.3",label = TRUE) + NoLegend() 

#B_Cell_CD79A
sub_cell<-subset(x = target_new,idents=c("20"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "B_Cell_CD79A"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Bs"
#Mast cell: KIT: "Mast_cells_KIT_pos" 
sub_cell<-subset(x = target_new,idents=c("22"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Mast_cells_MS4A2_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Masts"
#Th17_like_Cells_RORC_pos
sub_cell<-subset(x = target_new,idents=c("13"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Th17_like_Cells_RORC_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Th17_L"

#naive_CD4_T_Cells_CD4_pos_SELL_pos
sub_cell<-subset(x = target_new,idents=c("5"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "naive_CD4_T_Cells_CD4_pos_SELL_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "CD4_Ts_nai"
#memory_CD4_T_Cells_CD4_pos_S100A4_pos
sub_cell<-subset(x = target_new,idents=c("12"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "memory_CD4_T_Cells_CD4_pos_S100A4_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "CD4_Ts_mem"

### Regulatory T cells(Treg): FOXP3, CTLA-4, TNFRSF18; IL2RA (CD25 +)
sub_cell<-subset(x = target_new,idents=c("10"))
#sub_cell2<-subset(x = sub_cell,FOXP3>0) 
target_new$prefinal_cluster[Cells(sub_cell)] <- "Regulatory_T_Cells_CD4_pos_FOXP3_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <-"Treg"

#naive_Cytotoxic_T_Cells_CD8A_pos_IL7R_pos
#sub_cell<-subset(x = target_new,idents=c("12"))
#target_new$prefinal_cluster[Cells(sub_cell)] <- "naive_Cytotoxic_T_Cells_CD8A_pos_IL7R_pos"
#target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Cyt_Ts_nai"

#Effective_memory_Cytotoxic_T_Cells_CD8A_pos_GZMK_pos_GNLY_neg
sub_cell<-subset(x = target_new,idents=c("15","6","19"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Effective_memory_Cytotoxic_T_Cells_CD8A_pos_GZMK_pos_GNLY_neg"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Cyt_Ts_em"
#Effective_Cytotoxic_T_Cells_CD8A_pos_GZMB_pos_KLRB1_low_GNLY_pos
sub_cell<-subset(x = target_new,idents=c("14"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Effective_Cytotoxic_T_Cells_CD8A_pos_GZMB_pos_GNLY_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Cyt_Ts_eff"

#Natural_killer_Cells_CD3D_neg_NCAM1_mix_trophoblast
sub_cell<-subset(x = target_new,idents=c("2","7","8"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Natural_killer_Cells_mix_Trophoblast_CD3D_neg_NCAM1_KRT7_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "NKs_Troph_mix"

#Natural_killer_Cells_CD3D_neg_NCAM1_low_FCGR3A_pos【CD56 dim CD16 pos】
sub_cell<-subset(x = target_new,idents=c("16","17","21"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Natural_killer_Cells_CD3D_neg_NCAM1_low_FCGR3A_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "NKs_1"
#Natural_killer_Cells_CD3D_neg_NCAM1_low_FCGR3A_neg【CD56 dim CD16 neg】
sub_cell<-subset(x = target_new,idents=c("11"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Natural_killer_Cells_CD3D_neg_NCAM1_low_FCGR3A_neg"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "NKs_2"
#Natural_killer_Cells_CD3D_neg_NCAM1_high_FCGR3A_neg【CD56 bright CD16 neg】
sub_cell<-subset(x = target_new,idents=c("0","1","3","4"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Natural_killer_Cells_CD3D_neg_NCAM1_high_FCGR3A_neg"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "NKs_3"
#Natural_killer_Cells_CD3D_neg_NCAM1_high_MKI67_pos
sub_cell<-subset(x = target_new,idents=c("9","18"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Natural_killer_Cells_CD3D_neg_NCAM1_high_MKI67_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "NKs_p"
table(target_new$prefinal_cluster_brief)
# Bs    CD4_Ts_mem    CD4_Ts_nai    Cyt_Ts_eff     Cyt_Ts_em         Masts         NKs_1 
# 37            52            69            48           146            19           124 
#NKs_2         NKs_3         NKs_p NKs_Troph_mix        Th17_L          Treg 
#56           387            97           227            48            56 
DimPlot(object = target_new, group.by = "prefinal_cluster_brief",cols = ppCor_all2,label = TRUE)
NKsTsBs_Cells<-target_new

P1<-DimPlot(object = NKsTsBs_Cells, reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = NKsTsBs_Cells, reduction = "umap", group.by = "prefinal_cluster",cols = ppCor_all)+NoLegend()
P3<-DimPlot(object = NKsTsBs_Cells, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P4<-DimPlot(object = NKsTsBs_Cells, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P5<-DimPlot(object = NKsTsBs_Cells, reduction = "umap",group.by = "sample",cols = ppCor)
P6<-DimPlot(object = NKsTsBs_Cells, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_BsTsNKs<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/prefinal_subset_type_for_BsTsNKs_cells.png", type_plot_BsTsNKs,width=15, height=10)
split_plot_BsTsNKs<-DimPlot(object = NKsTsBs_Cells, reduction = "umap", split.by="Treat",group.by = "prefinal_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/prefinal_split_subset_type_for_BsTsNKs_cells.png", split_plot_BsTsNKs,width=13, height=6)

#5.Myeloid_Cell_AIF1
names(Cell_submian_list)
target_new<-Cell_submian_list[["Myeloid_Cell"]]
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 

#plot for Trophoblast_mix_EECs_KRT7
maker_plot_MyC1<-FeaturePlot(object = target_new, features = c("HLA-B","VIM","KRT7","GATA3","PERP","XAGE2","CDH1","EGFR","PAGE4","HLA-G","PAPPA2","CGA","CSH2","ERVFRD-1"),ncol=5,cols= c("grey", "red"))
maker_plot_MyC2<-FeaturePlot(object = target_new, features = c("EPCAM","PAEP","CAPS"),cols= c("grey", "red"),ncol=2)
#Cytotrophoblast（VCTs）
maker_plot_MyC3<-FeaturePlot(object = target_new, features = c("CDH1","EGFR","PAGE4","PEG10"),ncol=2,cols= c("grey", "red"))
#Extravilous trophoblast（EVTs）
maker_plot_MyC4<-FeaturePlot(object = target_new, features = c("HLA-G","MMP2","PAPPA2","DIO2","MMP11","FLT1","ITGA5","ADAM12","MCAM"),ncol=3,cols= c("grey", "red"))
#Syncytiotrophoblast (STBs)
maker_plot_MyC5<-FeaturePlot(object = target_new, features = c("CSH1","CGA","CSH2","CSHL1","GH2","PSG2","HOPX","TFAP2A","CYP19A1","ERVFRD-1"),cols= c("grey", "red"))

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_trophoblast_1.png", maker_plot_MyC1,width=25, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_trophoblast_2.png", maker_plot_MyC2,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_trophoblast_3.png", maker_plot_MyC3,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_trophoblast_4.png", maker_plot_MyC4,width=15, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_trophoblast_5.png", maker_plot_MyC5,width=20, height=15)


#plot for Myeloid_Cell
#	myeloid population: monocytes, DCs (cDC1 and cDC2, pDCs), macrophages (M1, M2), and mast Cells
maker_plot_Myc_1<-FeaturePlot(object = target_new, features = c("AIF1","CD14","HLA-DRA","S100A8","FCGR3A","CLEC9A","CD1C","CSF1R","CD163","CD68","CD209","LYVE1","CD86","FCN1","CCL2","MS4A3"),cols= c("grey", "red"),ncol= 4)
#granulocytic-monocytic progenitors (GMPs)
maker_plot_Myc_GMPs<-FeaturePlot(object = target_new, features = c("MS4A3"),cols= c("grey", "red"))
#	Dendritic Cells(DC): FCER1A, CST3, CLEC9A
maker_plot_Myc_DCs<-FeaturePlot(object = target_new, features = c("CLEC9A","CD1C","TSPAN13","XCR1","THBD","CADM1","CLEC10A","FCER1A"),cols= c("grey", "red"))
###Hofbauer Cells(fetal macrophages): 
#CD14+CSF1R+CD68+CD163+CD209+（resident marker）; #LYVE1+CD86+HLA-DR-HLA-DP-HLA-DQ-
maker_plot_Myc_HCs1<-FeaturePlot(object = target_new, features = c("CD14","CSF1R","CD163","CD68","CD209","LYVE1","CD86"),cols= c("grey", "red")) 
maker_plot_Myc_HCs2<-FeaturePlot(object = target_new, features = c("CCL2","CCL3","CCL4","CCL13","CCL4L2","CCL3L3","CXCL8"),cols= c("grey", "red")) 
##macrophage:
maker_plot_Mac_1<-FeaturePlot(object = target_new, features = c("CD14","CSF1R","CD163","CD68"),cols= c("grey", "red")) 
maker_plot_Mac_2<-FeaturePlot(object = target_new, features = c("F13A1","SLC40A1","GPNMB","S100A8","S100A9","CXCL8","C1QA","GPNMB","APOE","GPX3","VCAN","MARCO"),cols= c("grey", "red")) 
#Monocytes: 
#  Positive: ITGAX (CD11c) and HLA-DRA, FCN1
#	 CD14+ Monocytes (classical monocytes): CD14, LYZ, FTL, S100A9
#  FCGR3A+ monocytes (nonclassical monocytes): FCGR3A, MS4A7, VMO1，CCL2
#  granulocytic-monocytic progenitors (GMPs)：MS4A3
maker_plot_Myc_mono<-FeaturePlot(object = target_new, features = c("CD14","FCGR3A","LYZ","FCN1","CCL2","S100A8","S100A9","MS4A7","MS4A3"),cols= c("grey", "red"),ncol= 3)

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_whole.png", maker_plot_Myc_1,width=20, height=20)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_GMPs_cells.png", maker_plot_Myc_GMPs,width=6, height=6)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_DCs_cells.png", maker_plot_Myc_DCs,width=15, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_HC1_cells.png", maker_plot_Myc_HCs1,width=15, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_HC2_cells.png", maker_plot_Myc_HCs2,width=10, height=15)

ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_Mac1_cells.png", maker_plot_Mac_1,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_Mac2_cells.png", maker_plot_Mac_2,width=20, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/pre_final_subset_markers_for_Myeloid_Cell_Mono_cells.png", maker_plot_Myc_mono,width=15, height=15)

#Cell type determination
target_new$prefinal_cluster <- as.character(target_new$major_group_brief)
target_new$prefinal_cluster_brief <- as.character(target_new$major_group_brief)

Idents(object = target_new) <- "SCT_snn_res.2.5"
DimPlot(object = target_new, label = TRUE) + NoLegend() 
#	Dendritic Cells(DC):CLEC9A,CD1C
sub_cell<-subset(x = target_new,idents=c("31"))
target_new$prefinal_cluster[Cells(sub_cell)] <-  "Dendritic_Cells_AIF1_pos_CLEC9A_CD1C_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "DCs"
#classical monocytes: FCN1 pos CD14, LYZ, FTL, S100A9 
sub_cell<-subset(x = target_new,idents=c("25"))
target_new$prefinal_cluster[Cells(sub_cell)] <-  "Classical_monocytes_FCN1_CD14_S100A9_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "cMon"
#FCGR3A+ monocytes (nonclassical monocytes): FCGR3A, MS4A7, VMO1，CCL2
#MarcPhage_AIF1_pos_LYVE1_neg_CD68_pos
sub_cell<-subset(x = target_new,idents=c("3","4","5","8","9","10","12","15","18","20","21","22","28","30","32","35"))
#target_new$prefinal_cluster[Cells(sub_cell)] <-  "Nonclassical_monocytes_AIF1_pos_FCGR3A_pos"
#target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "ncMon"
target_new$prefinal_cluster[Cells(sub_cell)] <-  "MarcPhage_AIF1_pos_LYVE1_neg_CD68_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Mac"

#Hofbauer_Cell_AIF1_pos_LYVE1_pos
sub_cell<-subset(x = target_new,idents=c("1","2","11","17","23","26","29"))
target_new$prefinal_cluster[Cells(sub_cell)] <-  "Hofbauer_Cell_AIF1_pos_LYVE1_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "HCs"
sub_cell<-subset(x = target_new,subset = prefinal_cluster_brief %in% c("DCs","cMon","Mac","HCs"),invert = TRUE)
target_new$prefinal_cluster[Cells(sub_cell)] <-  "Myeloid_Trophoblast_mix_AIF1_KRT7_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "MyCs_Troph_mix"

DimPlot(object = target_new, group.by = "prefinal_cluster_brief",cols = ppCor_all2,label = TRUE)
Myeloid_Cell<-target_new

P1<-DimPlot(object = Myeloid_Cell, reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all)
P2<-DimPlot(object = Myeloid_Cell, reduction = "umap", group.by = "prefinal_cluster",cols = ppCor[c(2:10)])+NoLegend()
P3<-DimPlot(object = Myeloid_Cell, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
P4<-DimPlot(object = Myeloid_Cell, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
P5<-DimPlot(object = Myeloid_Cell, reduction = "umap",group.by = "sample",cols = ppCor)
P6<-DimPlot(object = Myeloid_Cell, reduction = "umap",group.by = "sample_code",cols = ppCor)
type_plot_MyCs<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/prefinal_subset_type_for_Myeloid_Cell.png", type_plot_Endo,width=15, height=10)
split_plot_MyCs<-DimPlot(object = Myeloid_Cell, reduction = "umap", split.by="Treat",group.by = "prefinal_cluster_brief",cols = ppCor_all2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/prefinal_split_subset_type_for_Myeloid_Cell.png", split_plot_Endo,width=13, height=6)

#6.For Trophoblast
names(Cell_submian_list)
target_new<-Cell_submian_list[["Trophoblast_mix_Epi"]]
target_new$prefinal_cluster <- as.character(target_new$major_group_brief)
target_new$prefinal_cluster_brief <- as.character(target_new$major_group_brief)
DimPlot(object = target_new, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() 

Idents(object = target_new) <- "SCT_snn_res.0.4"
DimPlot(object = target_new, label = TRUE) + NoLegend() 

#Decidual_stromal_cells :DCN pos DLK1 neg
sub_cell<-subset(x = target_new,idents=c("15"))  
target_new$prefinal_cluster[Cells(sub_cell)] <- "Stromal_cells_DCN_pos_DLK1_neg_KRT7_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "STCs_Troph_mix"
#Decidual_stromal_cells :DCN pos DLK1 neg
sub_cell<-subset(x = target_new,idents=c("16"))
target_new$prefinal_cluster[Cells(sub_cell)] <- "Epithelial_cells_PAEP_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Epi"
sub_cell<-subset(x = target_new,subset = prefinal_cluster_brief %in% c("STCs_Troph_mix","Epi"),invert = TRUE)
target_new$prefinal_cluster[Cells(sub_cell)] <- "Trophoblast_KRT7_pos"
target_new$prefinal_cluster_brief[Cells(sub_cell)] <- "Troph"
Trophoblast_mix_Epi<-target_new

##replace old metadata
dim(Erythrocyte@meta.data);dim(Endothelial_cells@meta.data)
dim(Stromal_cells@meta.data);dim(Trophoblast_mix_Epi@meta.data)
dim(Myeloid_Cell@meta.data);dim(NKsTsBs_Cells@meta.data)

metadata_new<-rbind(Erythrocyte@meta.data,Endothelial_cells@meta.data,
                    Stromal_cells@meta.data,Trophoblast_mix_Epi@meta.data,
                    Myeloid_Cell@meta.data,NKsTsBs_Cells@meta.data)
dim(metadata_new)#49988    49

table(metadata_new$prefinal_cluster_brief)
names(table(metadata_new$prefinal_cluster_brief))

metadata_new$final_major_group_brief<-metadata_new$prefinal_cluster_brief
metadata_new[which(metadata_new$prefinal_cluster_brief %in% c("Troph")),]$final_major_group_brief <- "Trophoblast_KRT7"
metadata_new[which(metadata_new$prefinal_cluster_brief %in% c("HCs","Mac","cMon","DCs","MyCs_Troph_mix")),]$final_major_group_brief <- "Myeloid_Cell_AIF1"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("FBs","MFs","PVs","STCs","STCs_Troph_mix")),]$final_major_group_brief <- "Stromal_cells_DCN"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("NKs_1","NKs_2","NKs_3","NKs_p")),]$final_major_group_brief <- "Natural_killer_Cell_NCAM1"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Th17_L","CD4_Ts_mem","CD4_Ts_nai","Treg","Cyt_Ts_eff","Cyt_Ts_em")),]$final_major_group_brief <- "T_Cell_CD3D"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Ery_no_prol","Ery_prol")),]$final_major_group_brief <- "Erythrocyte_HBA1"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Bs")),]$final_major_group_brief <- "B_Cell_CD79A"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Epi")),]$final_major_group_brief <- "Epithelial_Cell_EPCAM_PAEP"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Endo_1","Endo_2")),]$final_major_group_brief <- "Endothelial_Cell_PECAM1"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Masts")),]$final_major_group_brief <- "Mast_Cell_MS4A2"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("FBs_Troph_mix","MyCs_Troph_mix","NKs_Troph_mix","STCs_Troph_mix")),]$final_major_group_brief <- "Troph_mix_group"

table(metadata_new$final_major_group_brief)
#B_Cell_CD79A    Endothelial_Cell_PECAM1 Epithelial_Cell_EPCAM_PAEP           Erythrocyte_HBA1            Mast_Cell_MS4A2 
#      37                        131                         71                       1136                         19 
#Myeloid_Cell_AIF1  Natural_killer_Cell_NCAM1          Stromal_cells_DCN          T_Cell_CD3D            Troph_mix_group 
#     6387                        664                       3745                        419                       4188 
#Trophoblast_KRT7 
#     33191 

########replace the old metadata
Metadate_old<-target@meta.data
target2<-target#temp save
dim(metadata_new);dim(Metadate_old)
target@meta.data<-metadata_new
table(target$final_major_group_brief)

#saveRDS(target, file = "/mnt/data/chenwei/gongchen/3.seurat_result/singlelet_all_nine_forced_cell_low_cell_double_both_remove_sample_list_0818_integrated_seurat.rds")
umap_p1<-DimPlot(object = target, reduction = "umap", group.by = "final_major_group_brief",cols =ppCor_all2,label =F) +ggtitle('final_major_group_brief') 
umap_p2<-DimPlot(object = target, reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all2,label =T) +ggtitle('prefinal_cluster_brief') + NoLegend() 
umap_p3<-DimPlot(object = target, reduction = "umap", group.by = "raw_cluster_brif",cols =ppCor,label =T) +ggtitle('raw_cluster_brif') + NoLegend() 
umap_p1+umap_p2+umap_p3

#6.For Trophoblast
Trophoblast_raw<-subset(x = target, subset = final_major_group_brief == "Trophoblast_KRT7")
table(Trophoblast_raw$final_major_group_brief)
table(Trophoblast_raw$prefinal_cluster_brief)

DimPlot(object = Trophoblast_raw, reduction = "umap", group.by = "integrated_snn_res.0.6",cols =ppCor_all,label = TRUE)  +ggtitle('Trophoblast_KRT7')
DimPlot(object = Trophoblast_raw, reduction = "umap", group.by = "integrated_snn_res.0.6",split.by = "major_group_brief", cols =ppCor_all,label = TRUE)  +ggtitle('Trophoblast_KRT7')
DimPlot(object = Trophoblast_raw, reduction = "umap", group.by = "sample",split.by="sample", cols =ppCor_all,label = TRUE,ncol = 3)  +ggtitle('Trophoblast_KRT7')

DefaultAssay(Trophoblast_raw) <- "RNA"
Trophoblast_raw <- NormalizeData(Trophoblast_raw, verbose = FALSE)
#plot for Trophoblast: 5 X 3
plot_maker_Trophoblast<-FeaturePlot(object = Trophoblast_raw, features = c("HLA-B","VIM","KRT7","GATA3","PERP","XAGE2","CDH1","EGFR","PAGE4","HLA-G","PAPPA2","CGA","CSH2","ERVFRD-1"),ncol=5,cols= c("grey", "red"))
#FeaturePlot(object = target, features = c("EPCAM","PAEP","CAPS"),cols= c("grey", "red"),ncol=2)
#Cytotrophoblast（CBTs）2X 2
plot_maker_CBTs<-FeaturePlot(object = Trophoblast_raw, features = c("CDH1","EGFR","PAGE4","PEG10"),ncol=2,cols= c("grey", "red"))
#Extravilous trophoblast（EVTs）3x3
plot_maker_EVTs<-FeaturePlot(object = Trophoblast_raw, features = c("HLA-G","MMP2","PAPPA2","DIO2","MMP11","FLT1","ITGA5","ADAM12","MCAM"),ncol=3,cols= c("grey", "red"))
#Syncytiotrophoblast (STBs) 4 x 3
plot_maker_STBs<-FeaturePlot(object = Trophoblast_raw, features = c("CSH1","CGA","CSH2","CSHL1","GH2","PSG2","HOPX","TFAP2A","CYP19A1","ERVFRD-1"),cols= c("grey", "red"))
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/whole_trophoblast_markers_for_Trophoblast_raw.png", plot_maker_Trophoblast,width=25, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/CBTs_markers_for_Trophoblast_raw.png", plot_maker_CBTs,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/EVTs_markers_for_Trophoblast_raw.png", plot_maker_EVTs,width=15, height=15)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/STBs_markers_for_Trophoblast_raw.png", plot_maker_STBs,width=20, height=15)
plot_maker_profile<-FeaturePlot(object = Trophoblast_raw, features = c("PCNA", "TOP2A", "MCM6", "MKI67"),cols= c("grey", "red"),ncol=2)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/profile_markers_for_Trophoblast_raw.png", plot_maker_profile,width=10, height=10)

plot1<-DimPlot(Trophoblast_raw, group.by ="integrated_snn_res.0.8",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.0.8') 
plot2<-DimPlot(Trophoblast_raw, group.by ="integrated_snn_res.1",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.1') 
plot3<-DimPlot(Trophoblast_raw, group.by ="integrated_snn_res.1.2",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.1.2') 
plot4<-DimPlot(Trophoblast_raw, group.by ="integrated_snn_res.1.5",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.1.5') 
plot5<-DimPlot(Trophoblast_raw, group.by ="integrated_snn_res.2",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.2') 
plot6<-DimPlot(Trophoblast_raw, group.by ="integrated_snn_res.2.5",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.2.5') 
plot7<-DimPlot(Trophoblast_raw, group.by ="integrated_snn_res.3",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.3') 
plot8<-DimPlot(Trophoblast_raw, group.by ="integrated_snn_res.5",label = TRUE)+ NoLegend()+ ggtitle('integrated_snn_res.5') 
plot_merge_eight<-CombinePlots(plots = list(plot1, plot2,plot3,plot4,plot5,plot6,plot7,plot8),ncol = 4,legend = NULL)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/filter_sub_population/cluster_for_Trophoblast_raw.png",plot_merge_eight,width=16, height=8)

#Idents(object = Trophoblast_raw) <- "integrated_snn_res.0.6"
plot_06<-DimPlot(object = Trophoblast_raw, reduction = "umap", group.by = "integrated_snn_res.0.6",cols =ppCor_all,label = TRUE)  +ggtitle('Trophoblast_KRT7')
DimPlot(object = Trophoblast_raw, reduction = "umap", group.by = "integrated_snn_res.0.8",cols =ppCor_all,label = TRUE)  +ggtitle('Trophoblast_KRT7')
DimPlot(object = Trophoblast_raw, reduction = "umap", group.by = "integrated_snn_res.1",cols =ppCor_all,label = TRUE)  +ggtitle('Trophoblast_KRT7')
plot_12<-DimPlot(object = Trophoblast_raw, reduction = "umap", group.by = "integrated_snn_res.1.2",cols =ppCor_all,label = TRUE)  +ggtitle('Trophoblast_KRT7')
plot_06+plot_12
Idents(object = Trophoblast_raw) <- "integrated_snn_res.1.2"
Trophoblast_raw$prefinal_cluster <- as.character(Idents(Trophoblast_raw))
Trophoblast_raw$prefinal_cluster_brief <- as.character(Idents(Trophoblast_raw))

##### Trophoblast_KRT7  
DimPlot(object = Trophoblast_raw, group.by ="prefinal_cluster_brief",label = TRUE) + NoLegend() 

#Extravilous trophoblast（EVTs）
sub_cell<-subset(x = Trophoblast_raw,idents=c("6","10","14","16","29","31","21","18"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Extravilous_trophoblast_HLA-G_pos_PAPPA2_pos"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "EVTs_1"
sub_cell<-subset(x = Trophoblast_raw,idents=c("24","27"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Extravilous_trophoblast_HLA-G_pos_PAPPA2_neg"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "EVTs_2"
DimPlot(object = Trophoblast_raw, group.by ="prefinal_cluster_brief",cols =ppCor_all,label = TRUE) 
sub_cell<-subset(x = Trophoblast_raw,idents=c("2","8","15","17","13","20"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Extravilous_trophoblast_HLA-G_pos_PAPPA2_neg"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "EVTs_3"
DimPlot(object = Trophoblast_raw, group.by ="prefinal_cluster_brief",cols =ppCor_all,label = TRUE) 
####Syncytiotrophoblast（STBs）
sub_cell<-subset(x = Trophoblast_raw,idents=c("25","28"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Syncytiotrophoblast_CGA_pos_CSH2_pos"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "STBs_1"
sub_cell<-subset(x = Trophoblast_raw,idents=c("9"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Syncytiotrophoblast_CGA_pos_CSH2_neg"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "STBs_2"
DimPlot(object = Trophoblast_raw, group.by ="prefinal_cluster_brief",label = TRUE) + NoLegend() 
####Cytotrophoblast（CTBs）
sub_cell<-subset(x = Trophoblast_raw,idents=c("23"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_neg_CDH1_pos"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "CTBs_3"
sub_cell<-subset(x = Trophoblast_raw,idents=c("3","4","5","7","12","19","30","1","26","32"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_neg"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "CTBs_2"
sub_cell<-subset(x = Trophoblast_raw,idents=c("0","22"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_pos_PCNA_pos"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "CTBs_1"
sub_cell<-subset(x = Trophoblast_raw,idents=c("11"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_pos_PCNA_neg"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "CTBs_0"

DimPlot(object = Trophoblast_raw, group.by ="prefinal_cluster_brief",label = TRUE) + NoLegend() 

#another fine
#Extravilous trophoblast（EVTs）
Idents(object = Trophoblast_raw) <- "integrated_snn_res.3"
DimPlot(object = Trophoblast_raw, group.by ="integrated_snn_res.3", label = TRUE) + NoLegend() 

sub_cell<-subset(x = Trophoblast_raw,idents=c("32"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Extravilous_trophoblast_HLA-G_pos_PAPPA2_neg"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "EVTs_2"
sub_cell<-subset(x = Trophoblast_raw,idents=c("36"))
Trophoblast_raw$prefinal_cluster[Cells(sub_cell)] <- "Cytotrophoblasts_PAGE4_pos_HLA-G_neg_CGA_neg_MKI67_pos_PCNA_neg"
Trophoblast_raw$prefinal_cluster_brief[Cells(sub_cell)] <- "CTBs_0"
DimPlot(object = Trophoblast_raw, group.by ="prefinal_cluster_brief",label = TRUE) + NoLegend() 
#plot merge group
DimPlot(object = Trophoblast_raw, group.by ="integrated_snn_res.3", label = TRUE) + NoLegend() 

DimPlot(object = Trophoblast_raw,group.by ="prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() 
DimPlot(object = Trophoblast_raw,group.by ="prefinal_cluster_brief",cols =ppCor_all,label = TRUE)

##build final metadata
sub_cell<-subset(x = target,subset = final_major_group_brief %in% c("Trophoblast_KRT7"),invert = TRUE)
dim(Trophoblast_raw@meta.data);dim(sub_cell@meta.data)
metadata_new<-rbind(sub_cell@meta.data,Trophoblast_raw@meta.data)
########replace the old metadata
Metadate_old<-target@meta.data
target2<-target#temp save
dim(metadata_new);dim(Metadate_old)
target@meta.data<-metadata_new
DimPlot(object = target,group.by ="prefinal_cluster_brief",cols =ppCor_all,label = TRUE) + NoLegend() 
DimPlot(object = target,group.by ="prefinal_cluster",cols =ppCor_all,label = TRUE)
DimPlot(object = target,group.by ="final_major_group_brief",cols =ppCor_all2,label = TRUE) 

###recalculation for pure subgroups
raw_population<-names(table(target$final_major_group_brief)) 
Cell_submian_list<-list()
raw_population2 <-c("Endothelial_Cell_PECAM1","Erythrocyte_HBA1","Myeloid_Cell_AIF1","Natural_killer_Cell_NCAM1",
                    "Stromal_cells_DCN","T_Cell_CD3D","Troph_mix_group","Trophoblast_KRT7")
for (population in raw_population2) {
  #population<-"B_Cell_CD79A"  ##test line
  print(as.character(population))
  sub_cell<-subset(x = target, subset = final_major_group_brief == population)
  DefaultAssay(sub_cell) <- "RNA"    ## very important
  sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
  sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
  ElbowPlot(object = sub_cell)
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
  
  #PCs determination
  ElbowPlot(object = sub_cell)
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:10,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:10,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A1<-DimPlot(sub_cell,label = TRUE)
  
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:16,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:16,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A2<-DimPlot(sub_cell,label = TRUE)
  
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:20,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:20,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A3<-DimPlot(sub_cell,label = TRUE)
  
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A4<-DimPlot(sub_cell,label = TRUE)
  plot_for_dif_dims<-CombinePlots(plots = list(A1,A2,A3,A4),ncol = 2,legend = NULL)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population/dif_dim_for_",population,".png"), plot_for_dif_dims,width=24, height=20)
  
  #determination of resolution
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
  
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.4,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.6,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.8,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 1,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 1.2,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution =1.5,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 2,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 2.5,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution =3,verbose = FALSE)
  
  B1<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.4",label = TRUE) + NoLegend() +ggtitle('Cell_default0.4_sctransform')
  B2<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() +ggtitle("Cell_default0.6_sctransform")
  B3<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.8",label = TRUE) + NoLegend() +ggtitle("Cell_default0.8_sctransform")
  B4<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1",label = TRUE) + NoLegend() +ggtitle("Cell_default1_sctransform")
  B5<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1.2",label = TRUE) + NoLegend() +ggtitle("Cell_default1.2_sctransform")
  B6<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1.5",label = TRUE) + NoLegend() +ggtitle("Cell_default1.5_sctransform")
  B7<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.2",label = TRUE) + NoLegend() +ggtitle("Cell_default2_sctransform")
  B8<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.2.5",label = TRUE) + NoLegend() +ggtitle("Cell_default2.5_sctransform")
  B9<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.3",label = TRUE) + NoLegend() +ggtitle("Cell_default3_sctransform")
  
  plot_for_dif_res<-CombinePlots(plots = list(B1,B2,B3,B4,B5,B6,B7,B8,B9),ncol = 3,legend = NULL)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population/dif_res_dim25_for",population,".png"), plot_for_dif_res,width=30, height=30)
  
  DefaultAssay(sub_cell) <- "RNA"
  sub_cell <- NormalizeData(sub_cell, verbose = FALSE)
  plot_major_maker_subcluster<- FeaturePlot(object = sub_cell, features = major_maker,cols= c("grey", "red"),ncol=4)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population/major_maker_for_",population,".png"), plot_major_maker_subcluster,width=20, height=20)
  
  #plot
  P1<-DimPlot(object = sub_cell, reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all2)
  P2<-DimPlot(object = sub_cell, reduction = "umap", group.by = "final_major_group_brief",cols =ppCor_all2)
  P3<-DimPlot(object = sub_cell, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
  P4<-DimPlot(object = sub_cell, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
  P5<-DimPlot(object = sub_cell, reduction = "umap",group.by = "sample",cols = ppCor)
  P6<-DimPlot(object = sub_cell, reduction = "umap",group.by = "sample_code",cols = ppCor)
  type_plot<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population/anno_sample_group_for_",population,".png"), type_plot,width=26, height=16)
  plot_group_split<-DimPlot(object = sub_cell, reduction = "umap", split.by="Treat",group.by = "Treat",cols = ppCor,ncol=2)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population/treat_split_for_",population,".png"), plot_group_split,width=20, height=10)
  
  Cell_submian_list<-c(Cell_submian_list,list(sub_cell))
}
names(Cell_submian_list)<-raw_population2
saveRDS(Cell_submian_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_major_recal_selected_population_list.rds")
raw_population1<-c("B_Cell_CD79A","Epithelial_Cell_EPCAM_PAEP","Mast_Cell_MS4A2")


#including mixpopulation
metadata_new$final_major_group_brief_mix<-metadata_new$final_major_group_brief
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("FBs_Troph_mix","STCs_Troph_mix")),]$final_major_group_brief_mix <- "Stromal_cells_DCN"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("MyCs_Troph_mix")),]$final_major_group_brief_mix <- "Myeloid_Cell_AIF1"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("NKs_Troph_mix")),]$final_major_group_brief_mix <- "Natural_killer_Cell_NCAM1"
table(metadata_new$final_major_group_brief_mix)
table(metadata_new$final_major_group_brief)

########replace the old metadata
Metadate_old<-target@meta.data
target2<-target#temp save
dim(metadata_new);dim(Metadate_old)
target@meta.data<-metadata_new
DimPlot(object = target,group.by ="final_major_group_brief",cols =ppCor_all2,label = TRUE) 
DimPlot(object = target,group.by ="final_major_group_brief_mix",cols =ppCor_all,label = TRUE) + NoLegend() 
DimPlot(object = target,group.by ="prefinal_cluster_brief",cols =ppCor_all,label = TRUE)

###recalculation for pure subgroups
raw_population<-names(table(target$final_major_group_brief_mix)) 
Cell_submian_list<-list()
raw_population3 <-c("Endothelial_Cell_PECAM1","Erythrocyte_HBA1","Myeloid_Cell_AIF1","Natural_killer_Cell_NCAM1",
                    "Stromal_cells_DCN","T_Cell_CD3D","Trophoblast_KRT7")
for (population in raw_population3) {
  #population<-"B_Cell_CD79A"  ##test line
  print(as.character(population))
  sub_cell<-subset(x = target, subset = final_major_group_brief_mix == population)
  DefaultAssay(sub_cell) <- "RNA"    ## very important
  sub_cell <- SCTransform(object = sub_cell,verbose = FALSE) 
  sub_cell <- RunPCA(object = sub_cell,verbose = FALSE)
  ElbowPlot(object = sub_cell)
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.6, verbose = FALSE)
  
  #PCs determination
  ElbowPlot(object = sub_cell)
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:10,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:10,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A1<-DimPlot(sub_cell,label = TRUE)
  
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:16,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:16,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A2<-DimPlot(sub_cell,label = TRUE)
  
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:20,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:20,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A3<-DimPlot(sub_cell,label = TRUE)
  
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,verbose = FALSE)
  A4<-DimPlot(sub_cell,label = TRUE)
  plot_for_dif_dims<-CombinePlots(plots = list(A1,A2,A3,A4),ncol = 2,legend = NULL)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population_mix/dif_dim_for_",population,".png"), plot_for_dif_dims,width=24, height=20)
  
  #determination of resolution
  sub_cell <- RunUMAP(object = sub_cell, dims = 1:25,verbose = FALSE)
  sub_cell <- FindNeighbors(object = sub_cell, dims = 1:25,verbose = FALSE)
  
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.4,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.6,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 0.8,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 1,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 1.2,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution =1.5,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 2,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution = 2.5,verbose = FALSE)
  sub_cell <- FindClusters(object = sub_cell,resolution =3,verbose = FALSE)
  
  B1<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.4",label = TRUE) + NoLegend() +ggtitle('Cell_default0.4_sctransform')
  B2<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.6",label = TRUE) + NoLegend() +ggtitle("Cell_default0.6_sctransform")
  B3<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.0.8",label = TRUE) + NoLegend() +ggtitle("Cell_default0.8_sctransform")
  B4<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1",label = TRUE) + NoLegend() +ggtitle("Cell_default1_sctransform")
  B5<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1.2",label = TRUE) + NoLegend() +ggtitle("Cell_default1.2_sctransform")
  B6<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.1.5",label = TRUE) + NoLegend() +ggtitle("Cell_default1.5_sctransform")
  B7<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.2",label = TRUE) + NoLegend() +ggtitle("Cell_default2_sctransform")
  B8<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.2.5",label = TRUE) + NoLegend() +ggtitle("Cell_default2.5_sctransform")
  B9<-DimPlot(object = sub_cell, reduction = "umap",group.by = "SCT_snn_res.3",label = TRUE) + NoLegend() +ggtitle("Cell_default3_sctransform")
  
  plot_for_dif_res<-CombinePlots(plots = list(B1,B2,B3,B4,B5,B6,B7,B8,B9),ncol = 3,legend = NULL)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population_mix/dif_res_dim25_for",population,".png"), plot_for_dif_res,width=30, height=30)
  
  DefaultAssay(sub_cell) <- "RNA"
  sub_cell <- NormalizeData(sub_cell, verbose = FALSE)
  plot_major_maker_subcluster<- FeaturePlot(object = sub_cell, features = major_maker,cols= c("grey", "red"),ncol=4)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population_mix/major_maker_for_",population,".png"), plot_major_maker_subcluster,width=20, height=20)
  
  #plot
  P1<-DimPlot(object = sub_cell, reduction = "umap", group.by = "prefinal_cluster_brief",cols =ppCor_all2)
  P2<-DimPlot(object = sub_cell, reduction = "umap", group.by = "final_major_group_brief_mix",cols =ppCor_all2)
  P3<-DimPlot(object = sub_cell, reduction = "umap",group.by = "Phase",cols = ppCor_all2)
  P4<-DimPlot(object = sub_cell, reduction = "umap",  group.by = "Treat",cols =ppCor_all)
  P5<-DimPlot(object = sub_cell, reduction = "umap",group.by = "sample",cols = ppCor)
  P6<-DimPlot(object = sub_cell, reduction = "umap",group.by = "sample_code",cols = ppCor)
  type_plot<-CombinePlots(plots = list(P1,P2,P3,P4,P5,P6),ncol = 3,legend = NULL)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population_mix/anno_sample_group_for_",population,".png"), type_plot,width=26, height=16)
  plot_group_split<-DimPlot(object = sub_cell, reduction = "umap", split.by="Treat",group.by = "Treat",cols = ppCor,ncol=2)
  ggsave(paste0("/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_sub_population_mix/treat_split_for_",population,".png"), plot_group_split,width=20, height=10)
  
  Cell_submian_list<-c(Cell_submian_list,list(sub_cell))
}
names(Cell_submian_list)<-raw_population3
saveRDS(Cell_submian_list, file = "/mnt/data/chenwei/gongchen/3.seurat_result/1.global_cell_maker/final_mix_major_recal_selected_population_list.rds")

##final major subgroup
metadata_new<-target@meta.data
dim(metadata_new)
table(metadata_new$prefinal_cluster_brief)
table(metadata_new$prefinal_cluster)

table(metadata_new$prefinal_cluster_brief,metadata_new$sample_code)

names(table(metadata_new$prefinal_cluster_brief))
names(table(metadata_new$prefinal_cluster))
metadata_new$final_major_subgroup<-metadata_new$prefinal_cluster
metadata_new[which(metadata_new$prefinal_cluster_brief %in% c("HCs")),]$final_major_subgroup <- "Hofbauer_Cell_AIF1_LYVE1_pos"
metadata_new[which(metadata_new$prefinal_cluster_brief %in% c("Mac","cMon","DCs")),]$final_major_subgroup <- "Myeloid_Cell_AIF1_pos_LYVE1_neg"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("FBs","MFs","PVs")),]$final_major_subgroup <- "Fibroblasts_DCN_DLK1_pos"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("STCs")),]$final_major_subgroup <- "Stromal_cells_DCN_pos_DLK1_neg"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("NKs_1","NKs_2","NKs_3","NKs_p")),]$final_major_subgroup <- "Natural_killer_Cell_NCAM1_pos"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Th17_L","CD4_Ts_mem","CD4_Ts_nai","Treg","Cyt_Ts_eff","Cyt_Ts_em")),]$final_major_subgroup <- "T_Cell_CD3D_pos"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Ery_no_prol","Ery_prol")),]$final_major_subgroup <- "Erythrocyte_HBA1_pos"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Bs")),]$final_major_subgroup <- "B_Cell_CD79A_pos"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Epi")),]$final_major_subgroup <- "Epithelial_Cell_EPCAM_PAEP_pos"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Endo_1","Endo_2")),]$final_major_subgroup <- "Endothelial_Cell_PECAM1_pos"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Masts")),]$final_major_subgroup <- "Mast_Cell_MS4A2_pos"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("FBs_Troph_mix","MyCs_Troph_mix","NKs_Troph_mix","STCs_Troph_mix")),]$final_major_subgroup <- "Troph_mix_group"

table(metadata_new$final_major_subgroup)
names(table(metadata_new$prefinal_cluster_brief))
metadata_new$final_major_subgroup_brief<-metadata_new$prefinal_cluster_brief
metadata_new[which(metadata_new$prefinal_cluster_brief %in% c("HCs")),]$final_major_subgroup_brief <- "HCs"
metadata_new[which(metadata_new$prefinal_cluster_brief %in% c("Mac","cMon","DCs")),]$final_major_subgroup_brief <- "MyCs"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("FBs","MFs","PVs")),]$final_major_subgroup_brief <- "FBs"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("NKs_1","NKs_2","NKs_3","NKs_p")),]$final_major_subgroup_brief <- "NKs"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Th17_L","CD4_Ts_mem","CD4_Ts_nai","Treg","Cyt_Ts_eff","Cyt_Ts_em")),]$final_major_subgroup_brief <- "Ts"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Ery_no_prol","Ery_prol")),]$final_major_subgroup_brief <- "Ery"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("Endo_1","Endo_2")),]$final_major_subgroup_brief <- "Endo"
metadata_new[which(metadata_new$prefinal_cluster_brief %in%  c("FBs_Troph_mix","MyCs_Troph_mix","NKs_Troph_mix","STCs_Troph_mix")),]$final_major_subgroup_brief <- "Troph_mix"

table(metadata_new$final_major_subgroup_brief)
# Bs    CTBs_0    CTBs_1    CTBs_2    CTBs_3      Endo       Epi       Ery    EVTs_1    EVTs_2    EVTs_3       FBs       HCs     Masts      MyCs       NKs    STBs_1    STBs_2 
# 37      2144      3283     11650       913       131        71      1136      7058      1858      3422      2712      1974        19      4413       664       963      1900 
#STCs  Troph_mix        Ts 
# 1033      4188       419 

#add nature refer meta
sample_info<-metadata_new
table(sample_info$gemgroup)
sample_info2<-sample_info[,c("sample_code","gemgroup")]
colnames(sample_info2)<-c("sample_code","gemgroup_1")
sample_info3<-distinct(sample_info2)
nature_ref_meta<-read.table(file="/mnt/data/chenwei/gongchen/3.seurat_result/metadata_gongchen_cell_integrated_after_nature_singlelet_all_nine_forced_cell_low_cell_double_both_remove_0723_seurat.txt",stringsAsFactors=F,row.names=1,header =T)
head(nature_ref_meta)
table(nature_ref_meta$gemgroup)
table(nature_ref_meta$sample_code)

######merge meta information#####
meta_data_new<-merge(x =nature_ref_meta,y =sample_info3, by = "sample_code")
head(meta_data_new)
meta_data_new$raw_id2<-as.character(unlist(lapply(strsplit(meta_data_new$raw_id,"_"), function(x) x[1])))
meta_data_new$raw_id2 <-paste0(meta_data_new$raw_id2,"_",meta_data_new$gemgroup_1)
meta_data_new$BARCODE_new <-paste0(meta_data_new$BARCODE,"_",meta_data_new$gemgroup_1)
head(meta_data_new)
rownames(meta_data_new)<-meta_data_new$raw_id2
str(meta_data_new);tail(meta_data_new);dim(meta_data_new)##49988    44
meta_data_new2<-meta_data_new[,c("seurat_cluster","merge_cluster")]
colnames(meta_data_new2)<-c("seurat_cluster_nature","merge_cluster_nature")
head(meta_data_new2);head(sample_info)

##nature merge by our metadata
meta_data_new3<-merge(x =sample_info,y =meta_data_new2, by = 0)
head(meta_data_new3)
rownames(meta_data_new3)<-meta_data_new3$Row.names

########replace the old metadata
target2<-target
target@meta.data<-meta_data_new3
table(target$seurat_cluster_nature)
table(target$final_major_subgroup_brief)
saveRDS(target, file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
write.table(as.data.frame(target@meta.data), file="/mnt/data/chenwei/gongchen/3.seurat_result/metadata_for_final_GC_seurat_0906.txt",quote=F, row.names=T, col.names=T,sep="\t") 
umap_tx = integrated@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tx = so@meta.data$Treat)

ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, color=tx)) + geom_point() + 
  scale_color_manual(values=c("group1_untreated" = "darkblue", "group1_treated" = "darkred"))