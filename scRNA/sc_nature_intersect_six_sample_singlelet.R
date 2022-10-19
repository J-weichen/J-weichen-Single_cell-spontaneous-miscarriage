single_cell_placent_all<-readRDS(file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_all.rds")
gc_single_cell<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/intergrated_six_sample_after_filter_singlet_Cell_0525.rds")
head(single_cell_placent_all@meta.data)
gc_single_cell$annotation<-gc_single_cell$integrated_snn_res.0.6

##perform normalization
DefaultAssay(gc_single_cell) <- "RNA"
gc_single_cell2 <- SCTransform(gc_single_cell, verbose = FALSE)
DefaultAssay(single_cell_placent_all) <- "SCT"

dtlist <- list(SC_nature=single_cell_placent_all,GC_cell = gc_single_cell2)
intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 4000);length(intfts)#3401
dtlist <- PrepSCTIntegration(object.list = dtlist,anchor.features = intfts)
anchors <- FindIntegrationAnchors(object.list = dtlist, normalization.method = "SCT", anchor.features = intfts)
#saveRDS(anchors, file = "/mnt/data/chenwei/gongchen/3.seurat_result/gc_nature_single_cell_anchors_210525.rds")
#anchors<-readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/gc_nature_single_cell_anchors_210525.rds")
integrated2 <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated2 <- RunPCA(integrated2)
ElbowPlot(object = integrated2)
VizDimLoadings(object = integrated2, dims = 1:2, reduction = "pca")
integrated2 <- RunUMAP(integrated2, dims = 1:25)
integrated2<- RunTSNE(integrated2, dims = 1:25)      
######
plot0<-DimPlot(integrated2, reduction = "pca")
plot1<-DimPlot(integrated2, reduction = "umap")
plot2<-DimPlot(integrated2, reduction = "tsne")
CombinePlots(plots = list(plot0,plot1,plot2),legend="top",ncol=3)

head(integrated2@meta.data);tail(integrated2@meta.data)
table(integrated2@meta.data$orig.ident)
DimPlot(object = integrated2, group.by ="annotation",label = TRUE) + NoLegend() 

subset<-subset(x = integrated2,cells = names(which(is.na(integrated2$Treat)))) 
integrated2$Treat[Cells(subset)] <- "sc_nature"
table(integrated2@meta.data$Treat)
DimPlot(object = integrated2, group.by ="annotation",split.by = "Treat",label = TRUE) + NoLegend() 

integrated2$Treat2<-integrated2$Treat
subset_cell<-subset(x = integrated2,subset =Treat2 %in% c("Abortion","CTRL")) 
integrated2$Treat2[Cells(subset_cell)] <- "GC_data"
DimPlot(object = integrated2, group.by ="annotation",split.by = "Treat2",label = TRUE) + NoLegend() 

integrated2$annotation2<-integrated2$annotation
integrated2$annotation2[Cells(subset_cell)] <- "GC_data"
plot_group_split_anno3<-DimPlot(object = integrated2, group.by ="annotation2",split.by = "Treat2",label = TRUE,cols=ppCor_all) + NoLegend() 
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/sc_nature_intergrate_group_anno.png",plot_group_split_anno3,width=10, height=10)

#split by three group
integrated2$annotation3<-integrated2$annotation2
subset_cell<-subset(x = integrated2,subset =Treat %in% c("Abortion")) 
integrated2$annotation3[Cells(subset_cell)] <- "GC_Abortion"
subset_cell<-subset(x = integrated2,subset =Treat %in% c("CTRL")) 
integrated2$annotation3[Cells(subset_cell)] <- "GC_CTRL"
plot_group_split1<-DimPlot(object = integrated2, group.by ="annotation3",split.by = "Treat2",label = TRUE,cols=ppCor_all) + NoLegend() 
plot_group_split2<-DimPlot(object = integrated2, group.by ="annotation3",split.by = "Treat",label = TRUE,cols=ppCor_all) + NoLegend() 
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/sc_nature_intergrate_group_split1.png",plot_sample_split,width=10, height=10)
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/sc_nature_intergrate_group_split2.png",plot_group_split2,width=15, height=10)

#split by sample
head(integrated2@meta.data);tail(integrated2@meta.data)
subset<-subset(x = integrated2,cells = names(which(is.na(integrated2$sample_code)))) 
integrated2$sample_code[Cells(subset)] <- "sc_nature"
table(integrated2@meta.data$sample_code)
plot_sample_split<-DimPlot(object = integrated2, group.by ="annotation3",split.by = "sample_code",label = TRUE,cols=ppCor_all,ncol=4) + NoLegend() 
ggsave("/mnt/data/chenwei/gongchen/3.seurat_result/sc_nature_intergrate_sample_split.png"),plot_sample_split,width=20, height=10)
saveRDS(integrated2, file = "/mnt/data/chenwei/gongchen/3.seurat_result/sc_nature_intergrated_whole_six_sample_singlet_Cell.rds")