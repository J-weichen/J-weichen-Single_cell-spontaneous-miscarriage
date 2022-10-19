setwd("/media/helab/data1/xiaoge/tropho/filter_2000/tropho.k27ac")
library(cisTopic)
library(ggExtra)   
#start from bam files and peka file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
pathToBams <- './bam_hg38/all/'                                                                            
bamFiles <- paste(pathToBams, list.files(pathToBams), sep='')
regions <- './3v3/final_tri_peak_For_seurat.bed'
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, paired = T,regions, project.name='tropho.k27ac')
dat <- cbind(cisTopicObject@cell.data[which(cisTopicObject@cell.data$Total_reads<20000),]$Total_reads,cisTopicObject@cell.data[which(cisTopicObject@cell.data$Total_reads<20000),]$nCounts)
colnames(dat) <- c("fragments","fragments_in_peaks")
dat <- as.data.frame(dat)
pdf("scatter_qc.pdf")
p <- ggplot(dat,aes(x=fragments,y=fragments_in_peaks)) + geom_point(shape=19) +
  xlab("fragments") + ylab("fragments_in_peaks") 
p1 <- ggMarginal(p, type = "histogram")
p1
dev.off()

rna_anno <- read.csv("all_cobatch_list_sample_predicted_metadata.csv",header = T,row.names = 1)
rna_anno1 <- read.csv("A5_cobatch_list_sample_predicted_metadata.csv",header = T,row.names = 1)
rna_anno <- rbind(rna_anno,rna_anno1)
library(stropho.k27acngr)
new_name <- str_replace_all(as.character(colnames(tropho.k27ac)), "-", ".")

colnames(count_mat) <- new_name
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(harmony)
set.seed(1234)

count_mat <- as.matrix(cisTopicObject@count.matrix[,which(cisTopicObject@cell.data$Total_reads<12000&cisTopicObject@cell.data$nCounts<5000)]) 
counts <- counts[-c(93613, 93614, 93615, 93616, 93617, 93618, 93619, 93620, 93621, 93622, 93623, 93624, 93625, 93626, 93627, 93628),]
rm(tropho.k27ac)
# create object
k27ac.assay <- CreateChromatinAssay(
  counts = count_mat,
  sep = c(":", "-"),
  min.cells = 10,
  max.cells = "q90",
  min.features = 10,
  validate.fragments = F
)


tropho.k27ac <- CreateSeuratObject(counts = k27ac.assay, assay = "RNA")
tropho.k27ac <- subset(tropho.k27ac, nCount_peaks > 500 & nCount_peaks < 12000)
write.csv(tropho.k27ac@assays$RNA@counts,"final_3v3_cobatch_mat.csv")
# compute LSI
tropho.k27ac <- FindTopFeatures(tropho.k27ac, min.cutoff = "q10")
tropho.k27ac <- RunTFIDF(tropho.k27ac)
tropho.k27ac <- RunSVD(tropho.k27ac)

tropho.k27ac <- RunUMAP(object = tropho.k27ac, reduction = 'lsi', dims = 2:30)
tropho.k27ac <- FindNeighbors(object = tropho.k27ac, reduction = 'lsi', dims = 2:30)
tropho.k27ac <- FindClusters(object = tropho.k27ac, verbose = FALSE, resolution = 0.8,algorithm = 3)

label <- rownames(as.data.frame(tropho.k27ac@active.ident))
label <- as.data.frame(label)

label[grep("A3",label[,1]),2] <- "A3"
label[grep("A4",label[,1]),2] <- "A4"
label[grep("A5",label[,1]),2] <- "A5"
label[grep("N1",label[,1]),2] <- "N1"
label[grep("N2",label[,1]),2] <- "N2"
label[grep("N4",label[,1]),2] <- "N4"

label[grep("A3",label[,1]),3] <- "Abortion"
label[grep("A4",label[,1]),3] <- "Abortion"
label[grep("A5",label[,1]),3] <- "Abortion"
label[grep("N1",label[,1]),3] <- "Normal"
label[grep("N2",label[,1]),3] <- "Normal"
label[grep("N4",label[,1]),3] <- "Normal"

tropho.k27ac@meta.data$orig.ident <- as.factor(label[,2])
tropho.k27ac@meta.data$group <- as.factor(label[,3])
tropho.k27ac@meta.data$predict <- rna_anno[match(colnames(tropho.k27ac),rownames(rna_anno),nomatch = NA_integer_ ),7]
new_name <- str_replace_all(as.character(rownames(cisTopicObject@cell.data)), "-", ".")
rownames(cisTopicObject@cell.data) <- new_name
tropho.k27ac@meta.data$depth <- cisTopicObject@cell.data[which(rownames(cisTopicObject@cell.data) %in% colnames(count_mat)),]$Total_reads
tropho.k27ac <- ScaleData(tropho.k27ac)
tropho.k27ac <- RunPCA(
  object = tropho.k27ac, npcs = 20, verbose = FALSE
)

tropho.k27ac <- tropho.k27ac %>% RunHarmony("orig.ident", 
                                            plot_convergence = F, theta=0, lambda=0.01, asssy.use = "RNA")

tropho.k27ac <- RunSVD(
  object = tropho.k27ac,
  assay = 'RNA',
  reduction.key = 'harmony_',
  reduction.name = 'harmony'
)

pdf("tropho.k27ac_depth.pdf")
DepthCor(tropho.k27ac, reduction = "harmony",n=30)
dev.off()
tropho.k27ac <- RunUMAP(
  object = tropho.k27ac,
  reduction = 'harmony',
  dims = c(2:5,7:30)
)
tropho.k27ac <- FindNeighbors(
  object = tropho.k27ac,
  reduction = 'harmony',
  dims = c(2:5,7:30)
)
tropho.k27ac <- FindClusters(
  object = tropho.k27ac,
  algorithm = 3,
  resolution = 0.8,
  verbose = FALSE
)
identity(tropho.k27ac)

tropho.k27ac@meta.data$predict[which(tropho.k27ac$seurat_clusters=="9")] <- 'STB'
tropho.k27ac@active.ident <- tropho.k27ac@meta.data$predict
tropho.k27ac$seurat_clusters <- tropho.k27ac@meta.data$predict
tropho.k27ac@meta.data$seurat_clusters[which(tropho.k27ac$seurat_clusters=="13")] <- 'Ery'
double_id <- cbind(tropho.k27ac$group,tropho.k27ac$predict)
write.csv(double_id,"double_id.csv")

pdf("tropho.k27ac_umap.pdf",width = 20, height = 5)
P1 <- DimPlot(object = tropho.k27ac, group.by = "predict" , label.size = 5,label = TRUE) 
P2 <- DimPlot(object = tropho.k27ac ,group.by = "group" ,label.size = 5,label = TRUE) 
P3 <- FeaturePlot(
  object = tropho.k27ac,
  features = "depth",
  pt.size = 0.1,
  max.cutoff = 'q99',
  min.cutoff = 'q1'
)
P4 <- DimPlot(object = tropho.k27ac, group.by = "orig.ident" , label.size = 5,label = TRUE) 
P1 |P2 |P3 |P4
dev.off()

pdf("tropho.k27ac_raw_cluster_umap.pdf",width = 5, height = 5)
P4 <- DimPlot(object = tropho.k27ac ,label.size = 5,label = T, pt.size = 0.5) 
P4
dev.off()

write.csv(rownames(tropho.k27ac@assays$RNA),"3v3_peaks_seurat.csv")
## cell proportion
cell_number <- read.csv("cell_number.csv",header = T)
library(reshape2)
mydata<-reshape2::melt(cell_number,id.vars="celltype",variable.name="group",value.name="number")
ab_cell <- mydata[1:12,]
ct_cell <- mydata[13:24,]
ab_cell$proportion <- ab_cell$number/sum(ab_cell$number)
ct_cell$proportion <- ct_cell$number/sum(ct_cell$number)
cell_proportion <- rbind(ab_cell,ct_cell)
class(cell_proportion)
cell_proportion$proportion2 <- ifelse(cell_proportion$group == "control", 
                                      cell_proportion$proportion *1,
                                      cell_proportion$proportion * -1)
pdf("cell_proportion.pdf")
ggplot(data = cell_proportion) + geom_col(aes(x = factor(celltype), 
                                              y = proportion2, 
                                              fill = group)) + 
  scale_y_continuous(breaks = seq(from = -1, to = 1,
                                  by = 0.1),
                     labels = c(seq(1, 0, -0.1), 
                                seq(0.1, 1, 0.1)),
                     limits = c(NA,0.5)) + 
  coord_flip() + theme_bw()
dev.off()


##signature gene
CCG <- read.table("anno/3v3_CCG_signature_uniq.txt")
CCG_count <-tropho.k27ac@assays$RNA@counts[rownames(tropho.k27ac@assays$RNA@counts) %in% CCG$V1, ]
CCG_mean <- colMeans(CCG_count)
CCG_mean <- as.matrix(CCG_mean)

is.na(CCG_mean) <- 0
CCG_nor <- CCG_mean/tropho.k27ac@meta.data$nCount_RNA
tropho.k27ac@meta.data$CCG <- CCG_nor

plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'CCG',
  pt.size = 0.5,
  min.cutoff = 'q30',
  max.cutoff = 'q99',
  cols = c("lightgrey", "magenta")
)
plot1 <- VlnPlot(
  object = tropho.k27ac,
  group.by = "group" ,
  features = "CCG",
  pt.size = 0.1,
)
library(ggsignif)
library(ggpubr)
pdf("CCG_vln.pdf")
plot1
dev.off()
my_comparisons <- list(c("Abortion", "Normal"))
iris <- as.data.frame(plot1$data)
ggboxplot(iris, x="ident", y="CCG",add = "jitter")+ 
  stat_compare_means(label.y = 9.5)+
  ### 添加每个两两比较的显著性标记位置信息
  stat_compare_means(comparisons=my_comparisons, label ="p.signif")



# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
tropho.k27ac <- AddMotifs(
  object = tropho.k27ac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

CTB2_peaks <- FindMarkers(
  object = tropho.k27ac,
  group.by = "group",
  ident.1 = c('0'),
  ident.2 = c('3'),
  only.pos = TRUE,
  test.use = "LR",
  min.pct = 0.05
)
write.table(da_peaks,"EVT2_specific_peak.txt")

EVT_both_peaks <- FindMarkers(
  object = tropho.k27ac,
  ident.1 = c('4','5'),
  ident.2 = c('12','1','2','13','0','9','3','11','7','8','6','10'),
  only.pos = TRUE,
  test.use = "LR",
  min.pct = 0.05
)
write.table(rownames(CTB_both_peaks),"cluster_specific_peak/CTB_both_specific_peak.txt")
# test enrichment
EVT1_enriched.motifs <- FindMotifs(
  object = tropho.k27ac,
  features = rownames(EVT1_peaks),
  background = rownames(tropho.k27ac@assays$RNA@counts)
)
write.csv(CTB1_enriched.motifs,"CTB1_enriched.motifs.csv")
pdf("EVT1_motif.pdf",width = 8,height = 5)
MotifPlot(
  object = tropho.k27ac,
  motifs = head(rownames(EVT1_enriched.motifs))
)
dev.off()

DefaultAssay(tropho.k27ac) <- 'RNA'
library(magrittr)
library(dplyr)
tropho.k27ac.markers <- FindAllMarkers(object = tropho.k27ac, only.pos = TRUE, min.pct = 0.25, 
                                       thresh.use = 0.25)
top5 <- tropho.k27ac.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top10 <- tropho.k27ac.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

data <- tropho.k27ac@assays$RNA@data[which(rownames(tropho.k27ac@assays$RNA@data) %in% top10$gene),]
data <- as.data.frame(data)
library(pheatmap)
pdf("top10_heatmap.pdf",width = 12, height = 8)
p <- pheatmap(data, cluster_rows = F,show_rownames = F, show_colnames = F,
              cluster_cols = F,color = colorRampPalette(c("navy", "white", "firebrick3"))(10))
DoHeatmap(object = tropho.k27ac, features = top10$gene, slot = "data", disp.min = 0, disp.max = 2, 
          raster = T)
p
dev.off()

c3_markers <- tropho.k27ac.markers[c(which(tropho.k27ac.markers$cluster=='3')),]
write.table(c3_markers$gene,"c3_diff_peaks.txt")
normal_peaks <- FindMarkers(
  object = tropho.k27ac,
  group.by = 'group',
  ident.1 = c("Abortion"),
  ident.2 = c("Normal"),
  subset.ident = c("1","2","6","11","12","13"),
  only.pos = T,
  test.use = "LR",
  min.pct = 0.05
)

c0_u <- da_peaks[da_peaks$avg_logFC>0,]
write.table(normal_peaks[which(normal_peaks$p_val<0.05),],"./3v3/abortion_other_specific_peak.txt")
#c1
pdf("c1_Alu.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr8:142124913-142195063',
  pt.size = 0.5,
  max.cutoff = 'q70'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr14:104162971-104210271',
  pt.size = 0.5,
  max.cutoff = 'q70'
)
plot1 | plot2
dev.off()

pdf("c1_GALNT9_CACNA1I.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr12:132300773-132358192',
  pt.size = 0.5,
  max.cutoff = 'q70'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr22:39618130-39725030',
  pt.size = 0.5,
  max.cutoff = 'q80'
)
plot1 | plot2
dev.off()
#c2
pdf("c2_CD148_CD93.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr11:47979720-48057245',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q50'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr20:23073946-23147945',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q50'
)
plot1 | plot2
dev.off()
#c3
pdf("c3_ISM2_GPR78.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr14:77460079-77520288',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q20'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr4:8535701-8563444',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q20'
)
plot1 | plot2
dev.off()
#c4
pdf("c4_GFI1_ZBTB7A.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr1:92480148-92486883',
  pt.size = 0.5,
  max.cutoff = 'q70',
  min.cutoff = 'q10'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr19:4027086-4128644',
  pt.size = 0.5,
  max.cutoff = 'q90',
  min.cutoff = 'q50'
)
plot1 | plot2
dev.off()

#c5
pdf("c5_RXRA_LRG1.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr9:134320673-134393570',
  pt.size = 0.5,
  max.cutoff = 'q99',
  min.cutoff = 'q80'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr19:4532247-4552674',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q50'
)
plot1 | plot2
dev.off()
#c6
pdf("c6_CDH4_FGFR2.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr20:61414210-61472251',
  pt.size = 0.5,
  max.cutoff = 'q90',
  min.cutoff = 'q20'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr10:121260456-121304733',
  pt.size = 0.5,
  max.cutoff = 'q50',
  min.cutoff = 'q5'
)
plot1 | plot2
dev.off()
#c7
pdf("c7_TCF7_FGFR3.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr9:134568662-134626402',
  pt.size = 0.1,
  max.cutoff = 'q80',
  min.cutoff = 'q30'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr4:1735008-1818531',
  pt.size = 0.1,
  max.cutoff = 'q80',
  min.cutoff = 'q50'
)
plot1 | plot2
dev.off()
#c8
pdf("c8_PTPRCAP_TCF7.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr11:67427190-67448558',
  pt.size = 0.5,
  max.cutoff = 'q98',
  min.cutoff = 'q40'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr5:134119426-134126819',
  pt.size = 0.5,
  max.cutoff = 'q50',
  min.cutoff = 'q10'
)
plot1 | plot2
dev.off()
#c9
pdf("c9_CBFA2T3_GFI1B.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr16:88983283-88995490',
  pt.size = 0.1,
  max.cutoff = 'q70',
  min.cutoff = 'q10'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr9:132994673-133002895',
  pt.size = 0.1,
  max.cutoff = 'q80',
  min.cutoff = 'q10'
)
plot1 | plot2
dev.off()

#c10
pdf("c10_TEAD3_SYDE1.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr6:35483107-35498969',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q50'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr19:15106498-15120267',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q50'
)
plot1 | plot2
dev.off()

#c0
pdf("c0_CSF1R_VASH1.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr5:149425025-149468390',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q50'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr14:77218991-77251517',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q50'
)
plot1 | plot2
dev.off()

#c10
pdf("c10_PLAC4_RHOV.pdf",width = 10, height = 5)
plot1 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr21:42585332-42596977',
  pt.size = 0.5,
  max.cutoff = 'q80',
  min.cutoff = 'q10'
)
plot2 <- FeaturePlot(
  object = tropho.k27ac,
  features = 'chr15:41163851-41172018',
  pt.size = 0.5,
  max.cutoff = 'q60',
  min.cutoff = 'q10'
)
plot1 | plot2
dev.off()
# Get a list of motif position frequency matropho.k27acces from the JASPAR database
# Get a list of motif position frequency matropho.k27acces from the JASPAR database
pfm <- getMatropho.k27acxSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

chrom_counts <- GetAssayData(tropho.k27ac[["RNA"]])
chromatinassay <- CreateChromatinAssay(counts = chrom_counts, genome = "hg38")
#object <- CreateSeuratObject(counts = chromatinassay)
tropho.k27ac[["RNA"]] <- CreateChromatinAssay(counts = chrom_counts, genome = "hg38")

# add motif information
tropho.k27ac <- AddMotifs(
  object = tropho.k27ac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
# Create a new Motif object to store the results
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(tropho.k27ac), sep = c("-", "-")),
  pwm = pfm,
  genome = 'hg38',
  sep = c("-", "-"),
  use.counts = FALSE
)
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

tropho.k27ac[['RNA']] <- AddMotifs(
  object = tropho.k27ac[['RNA']],
  motif.object = motif,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
tropho.k27ac <- RegionStats(
  object = tropho.k27ac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c("-", "-")
)
class(CTB1_enriched.motifs$motif.name)
TFname = read.table("top_motif_uniq.txt",header=F,sep="\t")
TFname <- as.matrix(TFname)
EVT2_motif <- EVT2_enriched.motifs[match(TFname,EVT2_enriched.motifs$motif.name),c(6:7)]
all_TF_dat <- cbind(TFname,CTB1_motif,CTB2_motif,EVT1_motif,EVT2_motif)
write.csv(all_TF_dat,"all_TF_bubble.csv")

bubble <- read.csv("all_TF_bubble.csv",header = T)
bubble <- as.data.frame(bubble)
pp <- ggplot(bubble,aes(group,motif.name))
pbubble <-  pp + geom_point(aes(size=fold.enrichment,color=-1*log10(pvalue)))
plot <- pbubble + scale_colour_gradient(low="blue",high="yellow")
pdf("TF_top200.pdf",height =30)
plot
dev.off()
library(chromVAR)
#Computing motif activities
bgpeaks <- getBackgroundPeaks(tropho.k27ac)
tropho.k27ac <- RunChromVAR(
  object = tropho.k27ac,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

## pseudotime analysis
library(cicero)
library(monocle3)
library(SeuratWrappers)
library(reshape2)

rm(label_all)


tropho.k27ac <- readRDS("tropho.k27ac.rds")
DefaultAssay(tropho.k27ac) <- "RNA"
normal <- subset(tropho.k27ac,cells = colnames(tropho.k27ac[,which(tropho.k27ac@meta.data$group=="Normal")]))
Abortion <- subset(tropho.k27ac,cells = colnames(tropho.k27ac[,which(tropho.k27ac@meta.data$group=="Abortion")]))
normal_tropho <- normal[, normal$predict %in% c("CTB2", "STB", "EVT2")]
abortion_tropho <- Abortion[, Abortion$predict %in% c("CTB1", "STB", "EVT1")]

saveRDS(tropho.k27ac, file="tropho.k27ac.rds")
saveRDS(normal_tropho, file="normal_tropho.rds")
saveRDS(abortion_tropho, file="abortion_tropho.rds")

normal.cds <- as.cell_data_set(normal_tropho)
normal.cds <- cluster_cells(cds = normal.cds, reduction_method = "UMAP")
normal.cds <- learn_graph(normal.cds, use_partition = F)

abortion.cds <- as.cell_data_set(abortion_tropho)
abortion.cds <- cluster_cells(cds = abortion.cds, reduction_method = "UMAP")
abortion.cds <- learn_graph(abortion.cds, use_partition = F)

normal_umap <- normal_tropho@reductions$umap@cell.embeddings
abortion_umap <- abortion_tropho@reductions$umap@cell.embeddings
# order cells
normal.cds <- order_cells(normal.cds, reduction_method = "UMAP", root_cells = "N1_Trop.501_702_FKDL210008195.1a.D702.AK1680_11_6.rmdup_hg38.bam.sorted.bam" )
abortion.cds <- order_cells(abortion.cds, reduction_method = "UMAP", root_cells ="A3_tro.1.9_FKDL210056153.1a.D709.AK1680_2_7.rmdup_hg38.bam.sorted.bam")


# plot trajectories colored by pseudotime
pdf("normal_evt_pseudo.pdf")
plot_cells(
  cds = normal.evt.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
dev.off()

pdf("abortion_evt_pseudo.pdf")
plot_cells(
  cds = abortion.evt.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
dev.off()

tropho.k27ac <- AddMetaData(
  object = tropho.k27ac,
  metadata = normal.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "normal"
)

tropho.k27ac <- AddMetaData(
  object = tropho.k27ac,
  metadata = abortion.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "abortion"
)

pdf("abortion_normal_pseudo.pdf",width = 10, height = 5)
FeaturePlot(tropho.k27ac, c("normal", "abortion"), pt.size = 0.1) & scale_color_viridis_c()
dev.off()
###plot dynamics
normal_EVT <- normal[, normal$predict %in% c("CTB2",  "EVT2")]
abortion_EVT <- Abortion[, Abortion$predict %in% c("CTB1",  "EVT1")]

normal.evt.cds <- as.cell_data_set(normal_EVT)
normal.evt.cds <- cluster_cells(cds = normal.evt.cds, reduction_method = "UMAP")
normal.evt.cds <- learn_graph(normal.evt.cds, use_partition = F)

abortion.evt.cds <- as.cell_data_set(abortion_EVT)
abortion.evt.cds <- cluster_cells(cds = abortion.evt.cds, reduction_method = "UMAP")
abortion.evt.cds <- learn_graph(abortion.evt.cds, use_partition = F)

normal.evt.umap <- normal_EVT@reductions$umap@cell.embeddings
abortion.evt.umap <- abortion_EVT@reductions$umap@cell.embeddings
# order cells
normal.evt.cds <- order_cells(normal.evt.cds, reduction_method = "UMAP", 
                              root_cells = "N1_Trop.502_702_FKDL210008195.1a.D702.AK1681_10_2.rmdup_hg38.bam.sorted.bam" )
abortion.evt.cds <- order_cells(abortion.evt.cds, reduction_method = "UMAP", 
                                root_cells ="A3_tro.1.9_FKDL210056153.1a.D709.AK1680_2_7.rmdup_hg38.bam.sorted.bam")


tropho.k27ac <- AddMetaData(
  object = tropho.k27ac,
  metadata = normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "normal.evt"
)

tropho.k27ac <- AddMetaData(
  object = tropho.k27ac,
  metadata = abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "abortion.evt"
)

pdf("abortion_normal_evt_pseudo.pdf",width = 10, height = 5)
FeaturePlot(tropho.k27ac, c("normal.evt", "abortion.evt"), pt.size = 0.1) & scale_color_viridis_c()
dev.off()
dim(order_normal_count)
order_normal_count <- normal.evt.cds@assays@data@listData$logcounts[,order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
order_abortion_count <- abortion.evt.cds@assays@data@listData$logcounts[,order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
order_normal_count <- as.data.frame(order_normal_count)
order_abortion_count <- as.data.frame(order_abortion_count)
grep('chr10-128',rownames(order_abortion_count))

krt19_n <- order_normal_count[c(18002),]
krt19_a <- order_abortion_count[c(18002),]
krt19_n[2,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
krt19_a[2,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(krt19_n) <- c("krt19","pseudotime")
rownames(krt19_a) <- c("krt19","pseudotime")
krt19_n <- t(krt19_n)
krt19_n <- as.data.frame(krt19_n)
krt19_a <- t(krt19_a)
krt19_a <- as.data.frame(krt19_a)
pdf("krt19_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=krt19_n, aes(x=pseudotime, y=krt19, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=krt19_a, aes(x=pseudotime, y=krt19, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

krt7_n <- order_normal_count[c(6281 ,6282,  6283,  6284, 28816, 28819, 28820, 28821),]
krt7_a <- order_abortion_count[c(6281 ,6282,  6283,  6284, 28816, 28819, 28820, 28821),]
krt7_n[9,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
krt7_a[9,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(krt7_n) <- c(rownames(krt7_n[1:8,]),"pseudotime")
rownames(krt7_a) <- c(rownames(krt7_a[1:8,]),"pseudotime")
krt7_n <- t(krt7_n)
krt7_n <- as.data.frame(krt7_n)
krt7_n$krt7_avg <- rowMeans(krt7_n[,1:8])
krt7_a <- t(krt7_a)
krt7_a <- as.data.frame(krt7_a)
krt7_a$krt7_avg <- rowMeans(krt7_a[,1:8])
pdf("krt7_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=krt7_n, aes(x=pseudotime, y=krt7_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=krt7_a, aes(x=pseudotime, y=krt7_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()


cyp19a1_n <- order_normal_count[c(15286),]
cyp19a1_a <- order_abortion_count[c(15286),]
cyp19a1_n[2,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
cyp19a1_a[2,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(cyp19a1_n) <- c("cyp19a1","pseudotime")
rownames(cyp19a1_a) <- c("cyp19a1","pseudotime")
cyp19a1_n <- t(cyp19a1_n)
cyp19a1_n <- as.data.frame(cyp19a1_n)
cyp19a1_a <- t(cyp19a1_a)
cyp19a1_a <- as.data.frame(cyp19a1_a)
pdf("cyp19a1_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=cyp19a1_n, aes(x=pseudotime, y=cyp19a1, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=cyp19a1_a, aes(x=pseudotime, y=cyp19a1, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

egfr_n <- order_normal_count[c(16989,16992,16993,2088,17003,17004,17006,17007,17032,30572,17034,17037
                               ,17039,17040,17041,17063,2232,17065),]
egfr_a <- order_abortion_count[c(16989,16992,16993,2088,17003,17004,17006,17007,17032,30572,17034,17037
                                 ,17039,17040,17041,17063,2232,17065),]
egfr_n[19,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
egfr_a[19,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(egfr_n) <- c(rownames(egfr_n[1:18,]),"pseudotime")
rownames(egfr_a) <- c(rownames(egfr_a[1:18,]),"pseudotime")
egfr_n <- t(egfr_n)
egfr_n <- as.data.frame(egfr_n)
egfr_n$egfr_avg <- rowMeans(egfr_n[,1:18])
egfr_a <- t(egfr_a)
egfr_a <- as.data.frame(egfr_a)
egfr_a$egfr_avg <- rowMeans(egfr_a[,1:18])
pdf("egfr_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=egfr_n, aes(x=pseudotime, y=egfr_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=egfr_a, aes(x=pseudotime, y=egfr_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

SERPINE1_n <- order_normal_count[c(13046,13045,29947),]
SERPINE1_a <- order_abortion_count[c(13046,13045,29947),]
SERPINE1_n[4,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
SERPINE1_a[4,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(SERPINE1_n) <- c("SERPINE1_1","SERPINE1_2","SERPINE1_3","pseudotime")
rownames(SERPINE1_a) <- c("SERPINE1_1","SERPINE1_2","SERPINE1_3","pseudotime")
SERPINE1_n <- t(SERPINE1_n)
SERPINE1_n <- as.data.frame(SERPINE1_n)
SERPINE1_n$SERPINE1_avg <- rowMeans(SERPINE1_n[,1:3])

SERPINE1_a <- t(SERPINE1_a)
SERPINE1_a <- as.data.frame(SERPINE1_a)
SERPINE1_a$SERPINE1_avg <- rowMeans(SERPINE1_a[,1:3])
pdf("SERPINE1_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=SERPINE1_n, aes(x=pseudotime, y=SERPINE1_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=SERPINE1_a, aes(x=pseudotime, y=SERPINE1_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

MALAT1_n <- order_normal_count[c(2749,2750,2844,18149,18150,18159,18160,2845,2846,18161,30747),]
MALAT1_a <- order_abortion_count[c(2749,2750,2844,18149,18150,18159,18160,2845,2846,18161,30747),]
MALAT1_n[12,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
MALAT1_a[12,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(MALAT1_n) <- c(rownames(MALAT1_n[1:11,]),"pseudotime")
rownames(MALAT1_a) <- c(rownames(MALAT1_a[1:11,]),"pseudotime")
MALAT1_n <- t(MALAT1_n)
MALAT1_n <- as.data.frame(MALAT1_n)
MALAT1_n$MALAT1_avg <- rowMeans(MALAT1_n[,1:11])
MALAT1_a <- t(MALAT1_a)
MALAT1_a <- as.data.frame(MALAT1_a)
MALAT1_a$MALAT1_avg <- rowMeans(MALAT1_a[,1:11])
pdf("MALAT1_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=MALAT1_n, aes(x=pseudotime, y=MALAT1_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=MALAT1_a, aes(x=pseudotime, y=MALAT1_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

HLAG_n <- order_normal_count[c(28424),]
HLAG_a <- order_abortion_count[c(28424),]
HLAG_n[2,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
HLAG_a[2,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(HLAG_n) <- c("HLAG","pseudotime")
rownames(HLAG_a) <- c("HLAG","pseudotime")
HLAG_n <- t(HLAG_n)
HLAG_n <- as.data.frame(HLAG_n)
HLAG_a <- t(HLAG_a)
HLAG_a <- as.data.frame(HLAG_a)
pdf("HLAG_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=HLAG_n, aes(x=pseudotime, y=HLAG, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=HLAG_a, aes(x=pseudotime, y=HLAG, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

TIMP4_n <- order_normal_count[c(16663,16665,16681),]
TIMP4_a <- order_abortion_count[c(16663,16665,16681),]
TIMP4_n[4,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
TIMP4_a[4,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(TIMP4_n) <- c("TIMP4_1","TIMP4_2","TIMP4_3","pseudotime")
rownames(TIMP4_a) <- c("TIMP4_1","TIMP4_2","TIMP4_3","pseudotime")
TIMP4_n <- t(TIMP4_n)
TIMP4_n <- as.data.frame(TIMP4_n)
TIMP4_n$TIMP4_avg <- rowMeans(TIMP4_n[,1:3])
TIMP4_a <- t(TIMP4_a)
TIMP4_a <- as.data.frame(TIMP4_a)
TIMP4_a$TIMP4_avg <- rowMeans(TIMP4_a[,1:3])
pdf("TIMP4_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=TIMP4_n, aes(x=pseudotime, y=TIMP4_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=TIMP4_a, aes(x=pseudotime, y=TIMP4_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

TIMP3_n <- order_normal_count[c(12007,12012,12046),]
TIMP3_a <- order_abortion_count[c(12007,12012,12046),]
TIMP3_n[4,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
TIMP3_a[4,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(TIMP3_n) <- c("TIMP3_1","TIMP3_2","TIMP3_3","pseudotime")
rownames(TIMP3_a) <- c("TIMP3_1","TIMP3_2","TIMP3_3","pseudotime")
TIMP3_n <- t(TIMP3_n)
TIMP3_n <- as.data.frame(TIMP3_n)
TIMP3_n$TIMP3_avg <- rowMeans(TIMP3_n[,1:3])
TIMP3_a <- t(TIMP3_a)
TIMP3_a <- as.data.frame(TIMP3_a)
TIMP3_a$TIMP3_avg <- rowMeans(TIMP3_a[,1:3])
pdf("TIMP3_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=TIMP3_n, aes(x=pseudotime, y=TIMP3_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=TIMP3_a, aes(x=pseudotime, y=TIMP3_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

TIMP2_n <- order_normal_count[c(640 ,6523,  28851,  1008, 8481, 28852, 28859, 28860),]
TIMP2_a <- order_abortion_count[c(640 ,6523,  28851,  1008, 8481, 28852, 28859, 28860),]
TIMP2_n[9,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
TIMP2_a[9,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(TIMP2_n) <- c(rownames(TIMP2_n[1:8,]),"pseudotime")
rownames(TIMP2_a) <- c(rownames(TIMP2_a[1:8,]),"pseudotime")
TIMP2_n <- t(TIMP2_n)
TIMP2_n <- as.data.frame(TIMP2_n)
TIMP2_n$TIMP2_avg <- rowMeans(TIMP2_n[,1:8])
TIMP2_a <- t(TIMP2_a)
TIMP2_a <- as.data.frame(TIMP2_a)
TIMP2_a$TIMP2_avg <- rowMeans(TIMP2_a[,1:8])
pdf("TIMP2_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=TIMP2_n, aes(x=pseudotime, y=TIMP2_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=TIMP2_a, aes(x=pseudotime, y=TIMP2_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

TIMP1_n <- order_normal_count[c(9330,9352),]
TIMP1_a <- order_abortion_count[c(9330,9352),]
TIMP1_n[3,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
TIMP1_a[3,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(TIMP1_n) <- c("TIMP1_1","TIMP1_2","pseudotime")
rownames(TIMP1_a) <- c("TIMP1_1","TIMP1_2","pseudotime")
TIMP1_n <- t(TIMP1_n)
TIMP1_n <- as.data.frame(TIMP1_n)
TIMP1_n$TIMP1_avg <- rowMeans(TIMP1_n[,1:2])
TIMP1_a <- t(TIMP1_a)
TIMP1_a <- as.data.frame(TIMP1_a)
TIMP1_a$TIMP1_avg <- rowMeans(TIMP1_a[,1:2])
pdf("TIMP1_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=TIMP1_n, aes(x=pseudotime, y=TIMP1_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=TIMP1_a, aes(x=pseudotime, y=TIMP1_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

MMP11_n <- order_normal_count[c(31491),]
MMP11_a <- order_abortion_count[c(31491),]
MMP11_n[2,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
MMP11_a[2,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(MMP11_n) <- c("MMP11","pseudotime")
rownames(MMP11_a) <- c("MMP11","pseudotime")
MMP11_n <- t(MMP11_n)
MMP11_n <- as.data.frame(MMP11_n)
MMP11_a <- t(MMP11_a)
MMP11_a <- as.data.frame(MMP11_a)
pdf("MMP11_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=MMP11_n, aes(x=pseudotime, y=MMP11, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=MMP11_a, aes(x=pseudotime, y=MMP11, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

MMP17_n <- order_normal_count[c(11981,11982,11984,11986,11987),]
MMP17_a <- order_abortion_count[c(11981,11982,11984,11986,11987),]
MMP17_n[6,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
MMP17_a[6,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(MMP17_n) <- c("MMP17_1","MMP17_2","MMP17_3","MMP17_4","MMP17_5","pseudotime")
rownames(MMP17_a) <- c("MMP17_1","MMP17_2","MMP17_3","MMP17_4","MMP17_5","pseudotime")
MMP17_n <- t(MMP17_n)
MMP17_n <- as.data.frame(MMP17_n)
MMP17_n$MMP17_avg <- rowMeans(MMP17_n[,1:5])
MMP17_a <- t(MMP17_a)
MMP17_a <- as.data.frame(MMP17_a)
MMP17_a$MMP17_avg <- rowMeans(MMP17_a[,1:5])
pdf("MMP17_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=MMP17_n, aes(x=pseudotime, y=MMP17_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=MMP17_a, aes(x=pseudotime, y=MMP17_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()



RPL39_n <- order_normal_count[c(5576,6684),]
RPL39_a <- order_abortion_count[c(5576,6684),]
RPL39_n[3,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
RPL39_a[3,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(RPL39_n) <- c("RPL39_1","RPL39_2","pseudotime")
rownames(RPL39_a) <- c("RPL39_1","RPL39_2","pseudotime")
RPL39_n <- t(RPL39_n)
RPL39_n <- as.data.frame(RPL39_n)
RPL39_n$RPL39_avg <- rowMeans(RPL39_n[,1:2])
RPL39_a <- t(RPL39_a)
RPL39_a <- as.data.frame(RPL39_a)
RPL39_a$RPL39_avg <- rowMeans(RPL39_a[,1:2])
pdf("RPL39_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=RPL39_n, aes(x=pseudotime, y=RPL39_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=RPL39_a, aes(x=pseudotime, y=RPL39_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

TGFB1_n <- order_normal_count[c(11155,29209),]
TGFB1_a <- order_abortion_count[c(11155,29209),]
TGFB1_n[3,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
TGFB1_a[3,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(TGFB1_n) <- c("TGFB1_1","TGFB1_2","pseudotime")
rownames(TGFB1_a) <- c("TGFB1_1","TGFB1_2","pseudotime")
TGFB1_n <- t(TGFB1_n)
TGFB1_n <- as.data.frame(TGFB1_n)
TGFB1_n$TGFB1_avg <- rowMeans(TGFB1_n[,1:2])
TGFB1_a <- t(TGFB1_a)
TGFB1_a <- as.data.frame(TGFB1_a)
TGFB1_a$TGFB1_avg <- rowMeans(TGFB1_a[,1:2])
pdf("TGFB1_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=TGFB1_n, aes(x=pseudotime, y=TGFB1_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=TGFB1_a, aes(x=pseudotime, y=TGFB1_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

THBD_n <- order_normal_count[c(25415,25417,25455),]
THBD_a <- order_abortion_count[c(25415,25417,25455),]
THBD_n[4,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
THBD_a[4,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(THBD_n) <- c("THBD_1","THBD_2","THBD_3","pseudotime")
rownames(THBD_a) <- c("THBD_1","THBD_2","THBD_3","pseudotime")
THBD_n <- t(THBD_n)
THBD_n <- as.data.frame(THBD_n)
THBD_n$THBD_avg <- rowMeans(THBD_n[,1:3])
THBD_a <- t(THBD_a)
THBD_a <- as.data.frame(THBD_a)
THBD_a$THBD_avg <- rowMeans(THBD_a[,1:3])
pdf("THBD_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=THBD_n, aes(x=pseudotime, y=THBD_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=THBD_a, aes(x=pseudotime, y=THBD_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

stmn1_n <- order_normal_count[c(29724),]
stmn1_a <- order_abortion_count[c(29724),]
stmn1_n[2,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
stmn1_a[2,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(stmn1_n) <- c("stmn1","pseudotime")
rownames(stmn1_a) <- c("stmn1","pseudotime")
stmn1_n <- t(stmn1_n)
stmn1_n <- as.data.frame(stmn1_n)
stmn1_a <- t(stmn1_a)
stmn1_a <- as.data.frame(stmn1_a)
pdf("stmn1_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=stmn1_n, aes(x=pseudotime, y=stmn1, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=stmn1_a, aes(x=pseudotime, y=stmn1, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

LAMA4_n <- order_normal_count[c(8910,16427,30470,16459,17457, 8912, 8913,17461),]
LAMA4_a <- order_abortion_count[c(8910,16427,30470,16459,17457, 8912, 8913,17461),]
LAMA4_n[9,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
LAMA4_a[9,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(LAMA4_n) <- c(rownames(LAMA4_n[1:8,]),"pseudotime")
rownames(LAMA4_a) <- c(rownames(LAMA4_a[1:8,]),"pseudotime")
LAMA4_n <- t(LAMA4_n)
LAMA4_n <- as.data.frame(LAMA4_n)
LAMA4_n$LAMA4_avg <- rowMeans(LAMA4_n[,1:8])
LAMA4_a <- t(LAMA4_a)
LAMA4_a <- as.data.frame(LAMA4_a)
LAMA4_a$LAMA4_avg <- rowMeans(LAMA4_a[,1:8])
pdf("LAMA4_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=LAMA4_n, aes(x=pseudotime, y=LAMA4_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=LAMA4_a, aes(x=pseudotime, y=LAMA4_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

cdh1_n <- order_normal_count[c(25057,25080,25082,25083,25084,25085,25087,25119,25120,25142),]
cdh1_a <- order_abortion_count[c(25057,25080,25082,25083,25084,25085,25087,25119,25120,25142),]
cdh1_n[11,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
cdh1_a[11,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(cdh1_n) <- c("cdh1_1","cdh1_2","cdh1_3","cdh1_4","cdh1_5","cdh1_6","cdh1_7","cdh1_8","cdh1_9","cdh1_10","pseudotime")
rownames(cdh1_a) <- c("cdh1_1","cdh1_2","cdh1_3","cdh1_4","cdh1_5","cdh1_6","cdh1_7","cdh1_8","cdh1_9","cdh1_10","pseudotime")
cdh1_n <- t(cdh1_n)
cdh1_n <- as.data.frame(cdh1_n)
cdh1_n$cdh1_avg <- rowMeans(cdh1_n[,1:10])
cdh1_a <- t(cdh1_a)
cdh1_a <- as.data.frame(cdh1_a)
cdh1_a$cdh1_avg <- rowMeans(cdh1_a[,1:10])
pdf("cdh1_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=cdh1_n, aes(x=pseudotime, y=cdh1_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=cdh1_a, aes(x=pseudotime, y=cdh1_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

pappa2_n <- order_normal_count[c(2178,5379,5380,5381,28680),]
pappa2_a <- order_abortion_count[c(2178,5379,5380,5381,28680),]
pappa2_n[6,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
pappa2_a[6,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(pappa2_n) <- c("pappa2_1","pappa2_2","pappa2_3","pappa2_4","pappa2_5","pseudotime")
rownames(pappa2_a) <- c("pappa2_1","pappa2_2","pappa2_3","pappa2_4","pappa2_5","pseudotime")
pappa2_n <- t(pappa2_n)
pappa2_n <- as.data.frame(pappa2_n)
pappa2_n$pappa2_avg <- rowMeans(pappa2_n[,1:5])
pappa2_a <- t(pappa2_a)
pappa2_a <- as.data.frame(pappa2_a)
pappa2_a$pappa2_avg <- rowMeans(pappa2_a[,1:5])
pdf("pappa2_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=pappa2_n, aes(x=pseudotime, y=pappa2_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=pappa2_a, aes(x=pseudotime, y=pappa2_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

mki67_n <- order_normal_count[c(24687,24688,24694,24697,31804),]
mki67_a <- order_abortion_count[c(24687,24688,24694,24697,31804),]
mki67_n[6,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
mki67_a[6,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(mki67_n) <- c("mki67_1","mki67_2","mki67_3","mki67_4","mki67_5","pseudotime")
rownames(mki67_a) <- c("mki67_1","mki67_2","mki67_3","mki67_4","mki67_5","pseudotime")
mki67_n <- t(mki67_n)
mki67_n <- as.data.frame(mki67_n)
mki67_n$mki_avg <- rowMeans(mki67_n[,1:5])
mki67_a <- t(mki67_a)
mki67_a <- as.data.frame(mki67_a)
mki67_a$mki_avg <- rowMeans(mki67_a[,1:5])
pdf("mki67_avg_dynamic.pdf", width = 10, height = 5)
pn <- ggplot(data=mki67_n, aes(x=pseudotime, y=mki_avg, color=pseudotime))+
  geom_smooth()
pa <- ggplot(data=mki67_a, aes(x=pseudotime, y=mki_avg, color=pseudotime))+
  geom_smooth()
pn|pa
dev.off()

rna_3gene_dat <- read.table("./rna_pseudo/USM_EVT1EVT3_smooth_spline_scale_selected_genes_expression_for_Yaxi3.txt")
rna_pseudo <- read.table("./rna_pseudo/final_pseudotime_info_for_cells_in_USM_EVT1EVT3_branch.txt")
usm_pseudo <- t(as.matrix(rna_pseudo$pseudotime[order(rna_pseudo$pseudotime)]))
#colnames(usm_pseudo) <- colnames(rna_3gene_dat)
mki67_rna <- rbind(rna_3gene_dat[1,],usm_pseudo)
serpine1_rna <- rbind(rna_3gene_dat[2,],usm_pseudo)
timp3_rna <- rbind(rna_3gene_dat[3,],usm_pseudo)
rownames(mki67_rna) <- c("mki67_rna","pseudotime")
rownames(serpine1_rna) <- c("serpine1_rna","pseudotime")
rownames(timp3_rna) <- c("timp3_rna","pseudotime")
mki67_rna <- as.data.frame(t(mki67_rna))
serpine1_rna <- as.data.frame(t(serpine1_rna))
timp3_rna <- as.data.frame(t(timp3_rna))

pdf("3gene_rna_USM_0607.pdf", width = 15, height = 5)
p_mki67_rna <- ggplot(data=mki67_rna, aes(x=pseudotime, y=mki67_rna))+
  geom_smooth(method = "loess")
p_serpine1_rna <- ggplot(data=serpine1_rna, aes(x=pseudotime, y=serpine1_rna))+
  geom_smooth(method = "loess")
p_timp3_rna <- ggplot(data=timp3_rna, aes(x=pseudotime, y=timp3_rna, color=pseudotime))+
  geom_smooth(method = "loess")
p_mki67_rna | p_serpine1_rna | p_timp3_rna
dev.off()

## cell distribution
cell_number <- read.csv("cell_number.csv",header = T)
library(reshape2)
mydata<-melt(cell_number,id.vars="celltype",variable.name="group",value.name="number")
pdf("cell_distribution.pdf")

ggplot(mydata,aes(celltype,number,fill=group))+
  geom_bar(stat="identity",position="fill")
dev.off()
## cell density
pseudo_abortion <- cbind(tropho.k27ac$abortion,as.matrix(rep('abortion',7611)))
pseudo_normal <- cbind(tropho.k27ac$normal,as.matrix(rep('normal',7611)))
pseudo_density <- rbind(pseudo_normal,pseudo_abortion)
pseudo_density <- na.omit(pseudo_density)
colnames(pseudo_density) <- c("pseudotime","group")
pseudo_density <- as.data.frame(pseudo_density)
pseudo_density$pseudotime <- round(as.numeric(pseudo_density$pseudotime),3)
pdf("tropho_density.pdf",width = 10, height = 5)
ggplot(pseudo_density, aes(x = pseudotime)) +
  geom_density(aes(fill = group), trim = T,alpha=0.4)
dev.off()
# co-accessible
tropho.abortion.cds <- as.CellDataSet(x = abortion_tropho)
tropho.abortion.cicero <- make_cicero_cds(tropho.abortion.cds, reduced_coordinates = abortion_tropho@reductions$umap@cell.embeddings)
tropho.normal.cds <- as.CellDataSet(x = normal_tropho)
tropho.normal.cicero <- make_cicero_cds(tropho.normal.cds, reduced_coordinates = normal_tropho@reductions$umap@cell.embeddings)


# run cicero
normal_conns <- run_cicero(tropho.normal.cicero, genomic_coords = 'hg38.chrom.sizes', sample_num = 100)
head(normal_conns)
normal_ccans <- generate_ccans(normal_conns)
write.csv(normal_ccans,"control_ccans.csv")
write.table(rownames(normal_ccans),"control_ccans.txt")

normal_links <- ConnectionsToLinks(conns = normal_conns, ccans = normal_ccans)
Links(normal_tropho) <- normal_links
head(normal_links)

abortion_conns <- run_cicero(tropho.abortion.cicero, genomic_coords = 'hg38.chrom.sizes', sample_num = 100)
head(abortion_conns)
abortion_ccans <- generate_ccans(abortion_conns)
write.csv(abortion_ccans,"abortion_ccans.csv")
write.table(rownames(abortion_ccans),"abortion_ccans.txt")

abortion_links <- ConnectionsToLinks(conns = abortion_conns, ccans = abortion_ccans)
Links(abortion_tropho) <- abortion_links
head(abortion_links)

##top co-accessible links 
abortion_high_chr <- abortion_links@seqnames[which(abortion_links@elementMetadata@listData[["score"]] >= 0.8),]
normal_high_chr <- normal_links@seqnames[which(normal_links@elementMetadata@listData[["score"]] >= 0.6),]
abortion_high <- abortion_links@ranges[which(abortion_links@elementMetadata@listData[["score"]] >= 0.8),]
normal_high <- normal_links@ranges[which(normal_links@elementMetadata@listData[["score"]] >= 0.6),]
head(abortion_high_chr@lengths)
tail(normal_high)
normal_high_chr
head(abortion_high)
write.csv(normal_high,"normal_high_144_links.csv")

normal_chr14_link <- normal_links@ranges[which(normal_links@seqnames=='chr14'),]
View(normal_chr3_link[grep('138',normal_chr3_link),])

pdf("link_SNAI1_normal.pdf",height = 3)
link_plot <- LinkPlot(
  object = normal_tropho,
  region = c("chr20-49936808-50034714"),min.cutoff = 0.1
)
link_plot
dev.off()

pdf("link_SNAI1_abortion.pdf",height = 3)
link_plot <- LinkPlot(
  object = abortion_tropho,
  region = c("chr20-49936808-50034714"),min.cutoff = 0.1
)
link_plot
dev.off()

pdf("link_ATF3_normal.pdf",height = 3)
link_plot <- LinkPlot(
  object = normal_tropho,
  region = c("chr1-212520000-212660000"),min.cutoff = 0.001
)
link_plot
dev.off()

pdf("link_ATF3_abortion.pdf",height = 3)
link_plot <- LinkPlot(
  object = abortion_tropho,
  region = c("chr1-212520000-212660000"),min.cutoff = 0.001
)
link_plot
dev.off()


pdf("link_STAT3_abortion.pdf",height = 3)
link_plot <- LinkPlot(
object = abortion_tropho,
region = c("chr17-42194573-42447782"),min.cutoff = 0.05
)
link_plot
dev.off()

pdf("link_TFAP2A_normal.pdf",height = 3)
link_plot <- LinkPlot(
  object = normal_tropho,
  region = c("chr6-10371110-10435983"),min.cutoff = 0.05
)
link_plot
dev.off()
pdf("link_TFAP2A_abortion.pdf",height = 3)
link_plot <- LinkPlot(
  object = abortion_tropho,
  region = c("chr6-10371110-10435983"),min.cutoff = 0.05
)
link_plot
dev.off()
#ETV5(0) EGR1(0) TEAD4 NFE2L3 TFAP2C(0)/A NFYB(0) TBPL1(0) HEY1(0) BHLHE40(0)
#SOX4 MYCN(0) ELF3 ELF1(0) NFIL3

##co-embedding
inte <- readRDS("../placenta_Cobatch_and_RNA_intergrated_seurats.rds")
inte@meta.data[["Treat"]][grep("A5",colnames(inte))] <- 'abortion_Cobatch'
pdf("3v3_inte_raw_cluster.pdf",width = 20,height = 5)
DimPlot(inte, group.by = c("seurat_clusters", "Treat","final_anno"))
dev.off()
pdf("3v3_feature4.pdf",width = 15,height = 15)
DefaultAssay(inte) <- "ACTIVITY"
FeaturePlot(
  object = inte,
  features = c('DMPK', 'CCDC106', 'ASCL2', 'DNAJB6', 'MIR1229', 'RPS6KA4',
               'MAP7D1', 'RSPH6A', 'DLL4', 'MRPL58', 'CARMIL2', 'C1QL4',
               'GYG2', 'LRFN4', 'CATSPERZ', 'CSPG5'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4
)
dev.off()
FeaturePlot(
  object = inte,
  features = c('MGME1', 'PTGES3', 'ID3', 'HARBI1', 'NRF1', 'HNRNPA2B1',
               'TIA1','PAN2','H3F3B','CLIC4','NOG','HIST1H4H','DAP',
               'CEP85','FLJ37453','HIST1H2BE'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4
)

FeaturePlot(
  object = inte,
  features = c('C2CD4C', 'FIBCD1', 'NKX1-1', 'PTGER4', 'NKX1-1', 'MIR153-2',
               'CFAP46','SENCR','MIR1257','LDLRAD4','HMX1','PALD1','JPH3',
               'NUP210','MSX1','OLFM1'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4
)


RidgePlot(inte, group.by =  "final_anno",  assay = "ACTIVITY",
          features = c("C2CD4C"), ncol = 1)
DefaultAssay(inte) <- 'RNA'
ctrl_gene_myc <- FindMarkers(
  object = inte,
  slot = "data",
  subset.ident = c("Mycs"),
  group.by = 'Treat',
  logfc.threshold = 0.25,
  ident.2 = c("Abortion"),
  ident.1 = c("CTRL"),
  only.pos = T,
  test.use = "wilcox",
  min.pct = 0.1
)

control_c12_peak <- FindMarkers(
  object = tropho.k27ac,
  slot = "data",
  group.by = 'active.ident',
  ident.1 = c("12"),
  ident.2 = c("0","1","2","3","4","5","6","7","8","9","10","11","13"),
  subset.ident = c("Normal"),
  only.pos = T,
  test.use = "LR",
  min.pct = 0.01
)
DefaultAssay(inte) <- "ACTIVITY"
activity_myc_ab <- as.matrix(inte@assays$ACTIVITY@data[which(rownames(inte@assays$ACTIVITY@data) %in% rownames(abortion_gene_myc)),1]) 

pdf("3v3_myc_ab_activity.pdf",width = 15,height = 10)
FeaturePlot(
  object = inte,
  features = c("C3","ALDH2","ID2","TMEM176A","CHI3L1","LTC4S"),
  pt.size = 0.1,
  max.cutoff = 'q55',
  ncol = 3
)
dev.off()

activity_myc_ct <- as.matrix(inte@assays$ACTIVITY@data[which(rownames(inte@assays$ACTIVITY@data) %in% rownames(ctrl_gene_myc)),1]) 
pdf("3v3_myc_ct_activity.pdf",width = 20,height = 15)
FeaturePlot(
  object = inte,
  features = c("CLIP1","SMCHD1","TNFRSF1B","HSPA1B"
               ,"IL4I1","CD300E"),
  pt.size = 0.1,
  max.cutoff = 'q55',
  ncol = 3
)
dev.off()
## myc
myc_data <- rbind(ctrl_gene_myc,abortion_gene_myc)
myc_data[,6] <- c(rep("ctrl",158),rep("abortion",95))
colnames(myc_data) <- c("p_val","log2fc","pct.1","pct.2","p_val_adj","sig")
myc_data <- as.data.frame(myc_data)
myc_data[159:253,2] <- -myc_data[159:253,2]
library(ggplot2)
library(RColorBrewer)

pdf(file = "myc_volcano.pdf",width=8,height=8)
theme_set(theme_bw())
p <- ggplot(myc_data,aes(log2fc,-1*log10(p_val_adj),
                         color = sig))+geom_point()+
  xlim(-2.5,2.5) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-0.5,0.5),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))+ylim(0,15)
p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p
dev.off()
print(p)




label_rna[grep("Epi",label_rna[,1]),2] <- "non"
label_rna[grep("Endo",label_rna[,1]),2] <- "non"
label_rna[grep("Ery",label_rna[,1]),2] <- "non"
label_rna[grep("FBs",label_rna[,1]),2] <- "non"
label_rna[grep("HCs",label_rna[,1]),2] <- "non"
label_rna[grep("Masts",label_rna[,1]),2] <- "non"
label_rna[grep("MyCs",label_rna[,1]),2] <- "non"
label_rna[grep("NKs",label_rna[,1]),2] <- "non"
label_rna[grep("STCs",label_rna[,1]),2] <- "non"
label_rna[grep("Ts",label_rna[,1]),2] <- "non"
label_rna[grep("Bs",label_rna[,1]),2] <- "non"
label_rna[grep("CTBs",label_rna[,1]),2] <- "tropho"
label_rna[is.na(label_rna)]<-"tropho"

##figure5
label_rna <- as.data.frame(inte_rna@meta.data$final_anno)
label_rna[grep("CTBs_3",label_rna[,1]),2] <- "CTB_A"
label_rna[grep("CTBs_4",label_rna[,1]),2] <- "CTB_A"
label_rna[grep("CTBs_1",label_rna[,1]),2] <- "CTB_C"
label_rna[grep("EVTs_3",label_rna[,1]),2] <- "EVT_A"
label_rna[grep("EVTs_1",label_rna[,1]),2] <- "EVT_C"
label_rna[is.na(label_rna)]<-"other"
inte_rna@meta.data$label1 <- as.factor(label_rna[,2])

label_k27 <- as.data.frame(inte_k27@meta.data$final_anno)
label_k27[grep("CTBs_1",label_k27[,1]),2] <- "CTB_A"
label_k27[grep("CTBs_0",label_k27[,1]),2] <- "CTB_C"
label_k27[grep("EVTs_3",label_k27[,1]),2] <- "EVT_C"
label_k27[grep("EVTs_1",label_k27[,1]),2] <- "EVT_A"
label_k27[is.na(label_k27)]<-"other"
inte_k27@meta.data$label1 <- as.factor(label_k27[,2])


ctb_A_rna <- rowMeans(inte_rna@assays$RNA@data[c("TFDP1","ZBTB14","TFAP2A","EGR3","TFAP2C",
                                                 "TFAP2E","TFAP2A","TFAP2B","TFAP2B","TFAP2B",
                                                 "TFAP2C","EGR2","EGR4","TCFL5", "ZFP57",
                                                 "YY2","VEZF1","ZBTB7A","ZBTB7B","ZBTB7C",
                                                 "MYCN","MYC","CREB3L2","NPAS2","MAX","MLXIP",
                                                 "HEY2","HEY1","CREB3L4","TFEC","FOSL2","ELK4",
                                                 "ELF4","ETV6","ETV5","ETV2","ETS1","ETS2",
                                                 "PROX1","ELF2"),
                                               which(inte_rna@meta.data$label1 == c("CTBs_1A"))])

tf_rna <- cbind(as.matrix(ctb_c_rna),as.matrix(evt_c_rna),as.matrix(ctb_a_rna),as.matrix(evt_a_rna))
pdf("tf_rna_1.pdf",width = 5,height = 10)
p <- pheatmap(tf_rna, cluster_rows = F,cluster_cols = F)
p
dev.off()
library(pheatmap)
TF_rna <- inte_rna@assays$ACTIVITY@data[which(rownames(inte_rna@assays$ACTIVITY@data) %in% TF.id2[,]),which(inte_rna@meta.data$label1!="other")]
ctb_a_rna <- rowMeans(inte_rna@assays$RNA@data[match(rownames(k27_merge),rownames(inte_rna@assays$RNA@data)),which(inte_rna@meta.data$label1 == c("CTB_A"))])
ctb_c_rna <- rowMeans(inte_rna@assays$RNA@data[match(rownames(k27_merge),rownames(inte_rna@assays$RNA@data)),which(inte_rna@meta.data$label1 == c("CTB_C"))])
evt_a_rna <- rowMeans(inte_rna@assays$RNA@data[match(rownames(k27_merge),rownames(inte_rna@assays$RNA@data)),which(inte_rna@meta.data$label1 == c("EVT_A"))])
evt_c_rna <- rowMeans(inte_rna@assays$RNA@data[match(rownames(k27_merge),rownames(inte_rna@assays$RNA@data)),which(inte_rna@meta.data$label1 == c("EVT_C"))])
tf_combine_rna <- cbind(as.matrix(ctb_c_rna),as.matrix(ctb_a_rna),as.matrix(evt_c_rna),as.matrix(evt_a_rna))
colnames(tf_combine_rna) <- c("CTB_C_rna","CTB_A_rna","EVT_C_rna","EVT_A_rna") 
tf_combine_rna_nor <- cbind(tf_combine_rna[,1]/rowMeans(tf_combine_rna),
                            tf_combine_rna[,2]/rowMeans(tf_combine_rna),
                            tf_combine_rna[,3]/rowMeans(tf_combine_rna),
                            tf_combine_rna[,4]/rowMeans(tf_combine_rna))
tf_combine_rna_evt <- cbind(tf_combine_rna_nor[,3],tf_combine_rna_nor[,4])
tf_combine_rna_evt_nor <- tf_combine_rna_evt[apply(tf_combine_rna_evt, 1, function(x) sd(x)!=0),]
pdf("combine_TF_evt_rna_heatmap.pdf",width = 5, height = 15)
p_ctb <- pheatmap(tf_combine_rna_ctb_nor, cluster_rows = T,cluster_cols = F,scale = "row",
              color = colorRampPalette(colors = c("grey","grey","white","cyan"))(100),)
p_evt
dev.off()

k27_ctb_101 <- k27_ctb_nor[match(rownames(tf_combine_rna_ctb_nor),rownames(k27_ctb_nor)),]
k27_evt_101 <- k27_evt_nor[match(rownames(tf_combine_rna_evt_nor),rownames(k27_evt_nor)),]
k27_ctb_exp = k27_ctb_101[p_ctb$tree_row$order,]
k27_evt_exp = k27_evt_101[p_evt$tree_row$order,]
colnames(k27_ctb_exp) <- c("CTB_c_k27","CTB_u_k27")
colnames(k27_evt_exp) <- c("EVT_c_k27","EVT_u_k27")
k27_ctb_exp_rna = tf_combine_rna_ctb_nor[p_ctb$tree_row$order,]
k27_evt_exp_rna = tf_combine_rna_evt_nor[p_evt$tree_row$order,]
colnames(k27_ctb_exp_rna) <- c("CTB_c_rna","CTB_u_rna")
colnames(k27_evt_exp_rna) <- c("EVT_c_rna","EVT_u_rna")
ctb_all <- cbind(k27_ctb_exp_rna,k27_ctb_exp)
evt_all <- cbind(k27_evt_exp_rna,k27_evt_exp)
ctb_down <- ctb_all[which((ctb_all[,1]-ctb_all[,2])>0 & (ctb_all[,3]-ctb_all[,4])>0),]
ctb_up <- ctb_all[which((ctb_all[,1]-ctb_all[,2])<0 & (ctb_all[,3]-ctb_all[,4])<0),]
evt_down <- evt_all[which((evt_all[,1]-evt_all[,2])>0 & (evt_all[,3]-evt_all[,4])>0),]
evt_up <- evt_all[which((evt_all[,1]-evt_all[,2])<0 & (evt_all[,3]-evt_all[,4])<0),]
ctb_select <- rbind(ctb_down,ctb_up)
evt_select <- rbind(evt_down,evt_up)
pdf("select_TF_ctb_heatmap1.pdf",width = 5, height = 10)
p <- pheatmap(ctb_select, cluster_rows = T,cluster_cols = F,scale = "row",
              clustering_method="centroid",clustering_distance_rows	="euclidean",
              color = colorRampPalette(colors = c("grey","lightgrey","white","cyan","cyan"))(50),
              show_rownames = T)
p
dev.off()
ctb_select.sort <- ctb_select[match(rownames(ctb_select.km.cluster.sort),rownames(ctb_select)),]
ctb_select.km<-kmeans(ctb_select, centers=2,iter.max = 20)
ctb_select.km.cluster <- as.matrix(ctb_select.km$cluster)
ctb_select.km.cluster.sort <- as.matrix(ctb_select.km.cluster[order(ctb_select.km.cluster[,1]),])
pdf("select_TF_evt_heatmap.pdf",width = 5, height = 10)
p <- pheatmap(evt_select, cluster_rows = T,cluster_cols = F,scale = "row",
              clustering_method="centroid",clustering_distance_rows	="canberra",
              color = colorRampPalette(colors = c("grey","white","magenta"))(50),
              show_rownames = T)
p
dev.off()

ctb_a_k27 <- rowMeans(inte_k27@assays$ACTIVITY@data[match(rownames(tf_rna),rownames(inte_k27@assays$ACTIVITY@data)),which(inte_k27@meta.data$label1 == c("CTB_A"))])
ctb_c_k27 <- rowMeans(inte_k27@assays$ACTIVITY@data[match(rownames(tf_rna),rownames(inte_k27@assays$ACTIVITY@data)),which(inte_k27@meta.data$label1 == c("CTB_C"))])
evt_a_k27 <- rowMeans(inte_k27@assays$ACTIVITY@data[match(rownames(tf_rna),rownames(inte_k27@assays$ACTIVITY@data)),which(inte_k27@meta.data$label1 == c("EVT_A"))])
evt_c_k27 <- rowMeans(inte_k27@assays$ACTIVITY@data[match(rownames(tf_rna),rownames(inte_k27@assays$ACTIVITY@data)),which(inte_k27@meta.data$label1 == c("EVT_C"))])
tf_k27 <- cbind(as.matrix(ctb_c_k27),as.matrix(evt_c_k27),as.matrix(ctb_a_k27),as.matrix(evt_a_k27))
tf_k27_diff <- 
colnames(tf_k27) <- c("CTB_C_k27","EVT_C_k27","CTB_A_k27","EVT_A_k27") 
pdf("inte_TF_k27_heatmap.pdf",width = 5, height = 10)
p <- pheatmap(tf_k27, cluster_rows = F,cluster_cols = F)
p
dev.off()

tf_merge <- cbind(tf_rna,tf_k27)
pdf("inte_TF_merge_heatmap.pdf",width = 5, height = 30)
p <- pheatmap(tf_merge, cluster_rows = F,cluster_cols = T)
p
dev.off()


###20220507
## overlap of coaccessible genes and TFs
ccans_gene <- read.csv("ccan_gene.csv",header = F)
ccans_gene <- as.matrix(ccans_gene)
overlap_TF_gene <- ccans_gene[which(ccans_gene[,] %in%  TF.id2[,] ),]
table(inte_rna@meta.data$Treat)
subgroup_order<-c("CTBs_1","CTBs_2","STBs_1","STBs_2","STBs_3","EVTs_1","EVTs_2","EVTs_3","Epi","Endo","HCs","Mycs","STCs","FBs", "NKs","Ts","Bs","Masts","Ery")
inte_rna@meta.data$re_annotation_TFac1<- factor(inte_rna@meta.data$re_annotation_TFac, levels=subgroup_order,ordered=TRUE)
ctb_1_C_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                               which(inte_rna@meta.data$re_annotation_TFac1 == c("CTBs_1") &inte_rna@meta.data$Treat == c("CTRL"))])
ctb_1_A_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                                 which(inte_rna@meta.data$re_annotation_TFac1 == c("CTBs_1")&inte_rna@meta.data$Treat == c("Abortion"))])
ctb_2_C_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                                 which(inte_rna@meta.data$re_annotation_TFac1 == c("CTBs_2") &inte_rna@meta.data$Treat == c("CTRL"))])
ctb_2_A_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                                 which(inte_rna@meta.data$re_annotation_TFac1 == c("CTBs_2")&inte_rna@meta.data$Treat == c("Abortion"))])
evt_1_C_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                                 which(inte_rna@meta.data$re_annotation_TFac1 == c("EVTs_1") &inte_rna@meta.data$Treat == c("CTRL"))])
evt_1_A_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                                 which(inte_rna@meta.data$re_annotation_TFac1 == c("EVTs_1")&inte_rna@meta.data$Treat == c("Abortion"))])
evt_2_C_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                                 which(inte_rna@meta.data$re_annotation_TFac1 == c("EVTs_2") &inte_rna@meta.data$Treat == c("CTRL"))])
evt_2_A_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                                 which(inte_rna@meta.data$re_annotation_TFac1 == c("EVTs_2")&inte_rna@meta.data$Treat == c("Abortion"))])
evt_3_C_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                                 which(inte_rna@meta.data$re_annotation_TFac1 == c("EVTs_3") &inte_rna@meta.data$Treat == c("CTRL"))])
evt_3_A_rna <- rowMeans(inte_rna@assays$RNA@data[match( TF_select[,],rownames(inte_rna@assays$RNA@data)),
                                                 which(inte_rna@meta.data$re_annotation_TFac1 == c("EVTs_3")&inte_rna@meta.data$Treat == c("Abortion"))])

tf_rna <- cbind(as.matrix(ctb_1_C_rna),as.matrix(ctb_1_A_rna),as.matrix(ctb_2_C_rna),as.matrix(ctb_2_A_rna),
                as.matrix(evt_1_C_rna),as.matrix(evt_1_A_rna),as.matrix(evt_2_C_rna),
                as.matrix(evt_2_A_rna),as.matrix(evt_3_C_rna),as.matrix(evt_3_A_rna))
tf_rna_sort <- tf_rna[order(row.names(tf_rna)),]
colnames(tf_rna) <- c("CTBs_1_C","CTBs_1_A","CTBs_2_C","CTBs_2_A","EVTs_1_C",
                           "EVTs_1_A","EVTs_2_C","EVTs_2_A","EVTs_3_C","EVTs_3_A")
pdf("tf_rna_select.pdf",width = 6,height = 20)
library(pheatmap)
p <- pheatmap(tf_rna, scale = "row",cluster_rows = F,cluster_cols = F,show_colnames = T,show_rownames = T)
p
dev.off()
tf_ctb <- cbind(as.matrix(ctb_c_k27),as.matrix(ctb_a_k27),as.matrix(ctb_1_C_rna),
                as.matrix(ctb_2_C_rna),as.matrix(ctb_1_A_rna),as.matrix(ctb_2_A_rna))
colnames(tf_ctb) <- c("CTB_C_k27","CTB_A_k27","CTBs_1_C","CTBs_2_C","CTBs_1_A","CTBs_2_A") 
tf_evt <- cbind(as.matrix(evt_c_k27),as.matrix(evt_a_k27),as.matrix(evt_1_C_rna),
                as.matrix(evt_2_C_rna),as.matrix(evt_3_C_rna),as.matrix(evt_1_A_rna),
                as.matrix(evt_2_A_rna),as.matrix(evt_3_A_rna))
colnames(tf_evt) <- c("EVT_C_k27","EVT_A_k27","EVTs_1_C","EVTs_2_C","EVTs_3_C",
                      "EVTs_1_A","EVTs_2_A","EVTs_3_A") 
ctb_a_k27 <- rowMeans(inte_k27@assays$ACTIVITY@data[match(TF_select[,],rownames(inte_k27@assays$ACTIVITY@data)),which(inte_k27@meta.data$label1 == c("CTB_A"))])
ctb_c_k27 <- rowMeans(inte_k27@assays$ACTIVITY@data[match(TF_select[,],rownames(inte_k27@assays$ACTIVITY@data)),which(inte_k27@meta.data$label1 == c("CTB_C"))])
evt_a_k27 <- rowMeans(inte_k27@assays$ACTIVITY@data[match(TF_select[,],rownames(inte_k27@assays$ACTIVITY@data)),which(inte_k27@meta.data$label1 == c("EVT_A"))])
evt_c_k27 <- rowMeans(inte_k27@assays$ACTIVITY@data[match(TF_select[,],rownames(inte_k27@assays$ACTIVITY@data)),which(inte_k27@meta.data$label1 == c("EVT_C"))])
tf_k27 <- cbind(as.matrix(ctb_c_k27),as.matrix(ctb_a_k27),as.matrix(evt_c_k27),as.matrix(evt_a_k27))
colnames(tf_k27) <- c("CTB_C_k27","CTB_A_k27","EVT_C_k27","EVT_A_k27") 

pdf("TF_k27_select.pdf",width = 6, height = 20)
p <- pheatmap(tf_k27, scale = "row", cluster_rows = F,cluster_cols = F,show_colnames = T,show_rownames = T)
p
dev.off()
write.csv(tf_rna_sort,"TF_rna.csv")

tf_merge <- cbind(tf_ctb,tf_evt)
pdf("TF_k27_rna_merge_heatmap.pdf",width = 5, height = 20)
p <- pheatmap(tf_merge, cluster_rows = T,cluster_cols = F,scale = "row")
p
dev.off()

###heatmap_0526

k27_ctb <- cbind(as.matrix(ctb_c_k27),as.matrix(ctb_a_k27))
k27_ctb_nor <- cbind(k27_ctb[,1]/rowMeans(k27_ctb),k27_ctb[,2]/rowMeans(k27_ctb))
is.na(k27_ctb_nor) <- 0

k27_evt <- cbind(as.matrix(evt_c_k27),as.matrix(evt_a_k27))
k27_evt_nor <- cbind(k27_evt[,1]/rowMeans(k27_evt),k27_evt[,2]/rowMeans(k27_evt))
is.na(k27_evt_nor) <- 0
k27_evt <- k27_evt[apply(k27_evt, 1, function(x) sd(x)!=0),]
k27_merge <- cbind(k27_ctb_nor,k27_evt_nor)

colnames(k27_ctb) <- c("CTB_C_k27","CTB_A_k27")
colnames(k27_evt) <- c("EVT_C_k27","EVT_A_k27")

colnames(k27_merge) <- c("CTB_C_k27","CTB_A_k27","EVT_C_k27","EVT_A_k27") 
pdf("TF_k27_combine_nor_heatmap.pdf",width = 5, height = 20)
p <- pheatmap(k27_merge_exp, cluster_rows = F,cluster_cols = T,scale = "row",
              color = colorRampPalette(colors = c("grey","grey","magenta"))(100),
              show_rownames = T)
p
dev.off()
k27_merge <- k27_merge[-100,]
k27_merge_exp = k27_merge[p$tree_row$order,]
rownames(TF_order) <- rownames(k27_merge)
TF_order[,2] <- rownames(k27_merge)
TF_order_new <- as.data.frame(TF_order[order(TF_order[,1]),])

rna_ctb <- cbind(as.matrix(ctb_1_C_rna),
                as.matrix(ctb_2_C_rna),as.matrix(ctb_1_A_rna),as.matrix(ctb_2_A_rna))
colnames(rna_ctb) <- c("CTBs_1_C","CTBs_2_C","CTBs_1_A","CTBs_2_A") 
rna_evt <- cbind(as.matrix(evt_1_C_rna),
                 as.matrix(evt_2_C_rna),as.matrix(evt_3_C_rna),as.matrix(evt_1_A_rna),
                 as.matrix(evt_2_A_rna),as.matrix(evt_3_A_rna))
colnames(rna_evt) <- c("EVTs_1_C","EVTs_2_C","EVTs_3_C",
                       "EVTs_1_A","EVTs_2_A","EVTs_3_A")
rna_ctb_sort <- rna_ctb[pmatch(rownames(k27_merge_exp),rownames(rna_ctb)),]
rna_evt_sort <- rna_evt[pmatch(rownames(k27_merge_exp),rownames(rna_evt)),]
pdf("TF_rna_evt_heatmap.pdf",width = 5, height = 20)
p <- pheatmap(rna_evt_sort, cluster_rows = F,cluster_cols = T,scale = "row",
              color = colorRampPalette(colors = c("grey","grey","white","cyan","cyan"))(100),
              show_rownames = T)
p
dev.off()

##TF violin plot 20220606
inte_k27@active.ident <- inte_k27@meta.data$label1
pdf("TF_1_k27_violin.pdf",width = 30,height = 20)
VlnPlot(
  object = inte_k27,
  features = c("SNAI1","IRF8","HES5","IRF4","TFAP2A",
               "HNF1B","GLI2","CREB3L1","SOX9","NFIC"),
  group.by = "label1" ,
  split.by = "label1",
  pt.size = 0.1,
  ncol = 4,
  flip = T,
  log = T,
  slot = "scale.data",
  same.y.lims = T
)
dev.off()
my_comparisons=list(c("CTB_C","CTB_A"),
                    c("EVT_C","EVT_A"))
my_comparisons
ident <- c("CTB_C","CTB_A","EVT_C","EVT_A","other") 
library(ggpubr)
pdf("TEAD1_inte_k27_vln.pdf")
plot_TEAD1 <- VlnPlot(
  object = inte_k27,
  group.by = "label1" ,
  features = "TEAD1",
  pt.size = 0.1,
  ncol = 4,
  slot = "scale.data"
)
plot_TEAD1
dev.off()
pdf("TEAD1_inte_k27_vln.pdf")
ggviolin(plot_TEAD1$data, x="ident", y="integrated_TEAD1", fill = "ident",add = "mean", 
         add.params = list(color="white") ,  order= ident )+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")
dev.off()
pdf("TEAD1_inte_rna_vln.pdf",width = 20,height = 10)
plot_TEAD1 <- VlnPlot(
  object = inte_rna,
  group.by = "label1" ,
  features = "TEAD1",
  pt.size = 0.1,
  ncol = 4,
  slot = "scale.data"
)
plot_TEAD1
dev.off()
pdf("TEAD1_inte_rna_vln.pdf",width = 20,height = 20)
ggviolin(plot_TEAD1$data, x="ident", y="integrated_TEAD1", fill = "ident",add = "mean", 
         add.params = list(color="white") ,  order= ident )+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")
dev.off()
ScaleData(inte_rna@assays$RNA)
ScaleData(inte_k27@assays$ACTIVITY)
pdf("TF_1_rna_violin.pdf",width = 30,height = 20)
VlnPlot(
  object = inte_rna,
  features = c("SNAI1","IRF8","HES5","IRF4","TFAP2A",
               "HNF1B","GLI2","CREB3L1","SOX9","NFIC"),
  group.by = "label1" ,
  split.by = "label1",
  pt.size = 0.1,
  ncol = 4,
  flip = T,
  log = T,
  slot = "scale.data",
  same.y.lims = T
)
dev.off()
DefaultAssay(inte_k27) <-"ACTIVITY"
DefaultAssay(inte_rna) <-"RNA"
pdf("TF_2_k27_violin.pdf",width = 20,height = 10)
VlnPlot(
  object = inte_k27,
  features = c("PAX5","KLF2","TEAD1","LMX1B","ZBTB32",
               "GATA3","EBF3","PKNOX2","TBX21","CUX2"),
  group.by = "label1" ,
  pt.size = 0.1,
  ncol = 4
)
dev.off()
pdf("TF_2_rna_violin.pdf",width = 20,height = 10)
VlnPlot(
  object = inte_rna,
  features = c("PAX5","KLF2","TEAD1","LMX1B","ZBTB32",
               "GATA3","EBF3","PKNOX2","TBX21","CUX2"),
  group.by = "label1" ,
  pt.size = 0.1,
  ncol = 4
)
dev.off()

pdf("TF_3_k27_violin.pdf",width = 20,height = 10)
VlnPlot(
  object = inte_k27,
  features = c("GATA2",
               "ELF3","SPDEF"),
  group.by = "label1" ,
  pt.size = 0.1,
  ncol = 4
)
dev.off()
pdf("TF_3_rna_violin.pdf",width = 20,height = 10)
VlnPlot(
  object = inte_rna,
  features = c("GATA2",
               "ELF3","SPDEF"),
  group.by = "label1" ,
  pt.size = 0.1,
  ncol = 4
)
dev.off()
####
TF_select <- read.table("TF_select.txt")
library(dplyr)
TF_select <- distinct(TF_select)

pdf("TF_rna_dot.pdf",width = 80, height = 30)
dot_rna <- DotPlot(object = inte_rna, scale.by = "size",col.min = -2,col.max = 2 
                   ,assay = 'RNA',cols = c("blue","yellow"),dot.scale = 12,
        features = rownames(tf_rna[which(rownames(tf_rna)!="RELA"),]), group.by = "final_anno"
          )
dot_rna
dev.off()
pdf("TF_k27_dot.pdf",width = 80, height = 20)
dot_k27 <- DotPlot(object = inte_k27,scale.by = "size",col.min = -1.5,col.max = 1.5 ,
                   dot.scale = 12,assay = 'ACTIVITY',cols = c("blue","yellow"),
        features = rownames(tf_rna[which(rownames(tf_rna)!="RELA"),]), group.by = "label1"
)
dot_k27
dev.off()

pdf("inte_TF_k27_heatmap.pdf",width = 10, height = 10)
DoHeatmap(object = inte_k27, features = rownames(TF_rna), cells = colnames(inte_k27[,which(inte_k27@meta.data$label1!="other")]) ,group.by = "label1",
          slot = "scale.data",   draw.lines = T,raster = F)
dev.off()

###figure6
##

inte_rna <- subset(inte, tech == "scRNAseq")
pdf("tropho.rna_umap.pdf",width = 10, height = 5)
P4 <- DimPlot(object = inte_rna, split.by = "Treat" ,group.by = "final_anno", label.size = 5,label = TRUE) 
P4
dev.off()

inte_k27 <- subset(inte, tech != "scRNAseq")
pdf("tropho.k27_umap.pdf",width = 10, height = 5)
P4 <- DimPlot(object = inte_k27, split.by = "Treat" , group.by = "final_anno",label.size = 5,label = TRUE) 
P4
dev.off()


pdf("tropho.k27_label.pdf",width = 10, height = 5)
P4 <- DimPlot(object = inte_k27, group.by = "label1",label.size = 5,label = TRUE) 
P4
dev.off()
pdf("tropho.rna_label.pdf",width = 10, height = 5)
P4 <- DimPlot(object = inte_rna, group.by = "label",label.size = 5,label = TRUE) 
P4
dev.off()
#non-trophoblasts
DefaultAssay(inte_rna) <- 'integrated'
control_non_rna <- FindMarkers(
  object = inte_rna,
  slot = "data",
  group.by = 'Treat',
  ident.2 = c("Abortion"),
  ident.1 = c("CTRL"),
  subset.ident = c("Bs","Endo","Epi","Ery","Mycs","FBs","HCs","Masts",
                   "NKs","Ts","STCs"),
  only.pos = T,
  test.use = "LR",
  min.pct = 0.01
)

DefaultAssay(inte_k27) <- 'integrated'
k27_active_ident <- inte_k27@meta.data[["final_anno"]]
rownames(k27_active_ident) <- colnames(inte_k27@assays$integrated@data)
names(k27_active_ident) <- names(inte_k27@active.ident)
inte_k27@active.ident <- as.factor(k27_active_ident)
control_non_k27 <- FindMarkers(
  object = inte_k27,
  slot = "data",
  group.by = 'Treat',
  ident.2 = c("abortion_Cobatch"),
  ident.1 = c("normal_Cobatch"),
  subset.ident = c("Ery","MyCs","FBs","HCs","NKs","Ts","STCs"),
  only.pos = T,
  test.use = "LR",
  min.pct = 0.01
)

heatmap_DAG <- c(rownames(control_uniq_gene),rownames(abortion_non_k27))

DAG_rna <- inte_rna@assays$integrated@data[heatmap_DAG,]
pdf("inte_rna_heatmap.pdf",width = 20, height = 10)
DoHeatmap(object = inte_rna, features = heatmap_DAG, cells = colnames(inte_rna[,which(inte_rna@meta.data$label=="non")]) ,group.by = "final_anno_1",
          slot = "data",  disp.min = -2,disp.max = 2, draw.lines =F,raster = T)
dev.off()


DAG_k27 <- inte_k27@assays$integrated@data[k27_heatmap_DAG,]
pdf("inte_k27_heatmap.pdf",width = 20, height = 10)
DoHeatmap(object = inte_k27, features = heatmap_DAG, cells = colnames(inte_k27[,which(inte_k27@meta.data$label1 == "non")]) ,group.by = "final_anno_1",
          slot = "data",  disp.min = 0,disp.max = 2, draw.lines =F,raster = T)
dev.off()


library(tidyr)
rna_anno_mat <- as.data.frame(cbind(inte_rna@meta.data$final_anno,inte_rna@meta.data$Treat))
colnames(rna_anno_mat) <- c("anno","treat")
rna_anno_mat <- unite(rna_anno_mat, "final_anno_treat", "anno", "treat",
                      sep = "_", remove = FALSE)
rna_anno_mat[rna_anno_mat$final_anno_treat == "NKs_Abortion",]$final_anno_treat <- "Lyms_Abortion"
inte_rna@meta.data$final_anno_1 <- rna_anno_mat$final_anno_treat

k27_anno_mat[k27_anno_mat$final_anno_treat == "Masts_abortion_Cobatch",]$final_anno_treat <- "MyCs_abortion_Cobatch"
inte_k27@meta.data$final_anno_1 <- k27_anno_mat$final_anno_treat

abortion_uniq_gene <- rbind(abortion_non_k27,abortion_non_rna)
abortion_uniq_gene$name <- rownames(abortion_uniq_gene)
abortion_uniq_gene <- abortion_uniq_gene[!duplicated(abortion_uniq_gene$name),]
write.csv(abortion_uniq_gene,"abortion_non_uniq_gene.csv")
control_uniq_gene <- rbind(control_non_k27,control_non_rna)
control_uniq_gene$name <- rownames(control_uniq_gene)


control_uniq_gene <- control_uniq_gene[!duplicated(control_uniq_gene$name),]
write.csv(control_uniq_gene,"control_non_uniq_gene.csv")

##TF enhancer H3K27ac activity
TF.id <- read.csv("EVT1_enriched.motifs.csv",header = T)
TF.id <- as.matrix(TF.id)
TF.id2 <- as.matrix(TF.id1[-grep("::",TF.id1),])
write.table(toupper(TF.id2),"TF_name.txt")
hg38_ref <- read.table("hg38.refGene.uniq.-1k+100.bed",header = F)
TF_hg38 <- hg38_ref[which(hg38_ref[,4] %in% TF.id2[,] ),]
TF_hg38 <- TF_hg38[-grep("_",TF_hg38[,1]),]
write.csv(TF_hg38,"TF_hg38_raw.csv")

##rose plot
data0510 <- read.table("diff_peak_number.txt")
data0510 <- as.data.frame(data0510)
colnames(data0510) <- c("number","celltype","group")
pdf("peak_rose.pdf")
rose=ggplot(data0510)+
  geom_bar(stat="identity",aes(y=number,x=celltype,fill=group))+ 
  coord_polar()
rose
dev.off()

##figure 3 heatmap
k27_CTB_EVT_norm_mat <- inte_k27@assays$ACTIVITY@data[,match(colnames(order_normal_count),colnames(inte_k27@assays$ACTIVITY@data),)]
k27_CTB_EVT_usm_mat <- inte_k27@assays$ACTIVITY@data[,match(colnames(order_abortion_count),colnames(inte_k27@assays$ACTIVITY@data),)]
rna_gene <- read.table("developped_genes_rearranged_along_CTBs2EVTs_branchment_pseudotime.txt")
rownames(rna_gene) <- rna_gene$id
k27_CTB_EVT_ct_1174_mat <- k27_CTB_EVT_norm_mat[match(rownames(rna_gene),rownames(k27_CTB_EVT_norm_mat),nomatch=0),]
k27_CTB_EVT_usm_1174_mat <- k27_CTB_EVT_usm_mat[match(rownames(rna_gene),rownames(k27_CTB_EVT_usm_mat),nomatch=0),]
write.table(rownames(k27_CTB_EVT_usm_1174_mat),"k27_CTB2EVT_gene.txt")
write.csv(k27_CTB_EVT_ct_1174_mat,"k27_CTB_EVT_ct_501_mat.csv")
write.csv(k27_CTB_EVT_usm_1174_mat,"k27_CTB_EVT_usm_501_mat.csv")

# 构建列注释信息

CTRL.k27 <- read.csv("k27_CTB_EVT_ct_501_mat.csv",header = T,row.names = 1)
CTRL.k27.norm <- scale(CTRL.k27)
pdf("k27_CTB2EVT_ct_501_heatmp.pdf",width = 15,height = 20)
p <- pheatmap(CTRL.k27.norm, scale = "row",cluster_rows = T,
              cluster_cols = F,show_colnames = F,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),show_rownames = T)
p
dev.off()

USM.k27 <- read.csv("k27_CTB_EVT_usm_501_mat.csv",header = T,row.names = 1)
USM.k27.norm <- scale(USM.k27)
pdf("k27_CTB2EVT_usm_501_heatmp.pdf",width = 15,height = 20)
p <- pheatmap(USM.k27.norm, scale = "row",cluster_rows = T,
              cluster_cols = F,show_colnames = F,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),show_rownames = T)
p
dev.off()
write.csv(p$kmeans$cluster,"k27_CTB2EVT_usm_gene_4cluster.csv")
USM.k27.norm.order <- USM.k27.norm[match(rownames(CTRL.k27.norm),rownames(USM.k27.norm)),]
combine.k27.norm <- cbind(CTRL.k27.norm,USM.k27.norm)
combine.k27 <- cbind(k27_CTB_EVT_ct_501_mat,k27_CTB_EVT_usm_501_mat)
write.csv(combine.k27,"combine.k27.raw.csv")
write.table(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)],"USM_CTB2EVT_pseudotime.txt")
dim(combine.k27)
dim(CTRL.k27.norm)
dim(USM.k27.norm)
k27_CTB_EVT_ct_501_mat <- as.matrix(k27_CTB_EVT_ct_1174_mat)
signature1_n <- CTRL.k27.norm[c("ACTN4","CEBPB","COL4A2","EFNA1","MYH9","ARHGDIB","ASCL2",
                                "KRT19","ORC6","SNAI1","SOCS3","SPINT2","STK4","TPD52L1",
                                "ACAN","ADAM12","ADAM19","ADAM9","ANGPT4"),]
signature1_n <- as.data.frame(signature1_n)
signature1_n[20,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
signature1_n <- t(signature1_n)
colnames(signature1_n) <- c("ACTN4","CEBPB","COL4A2","EFNA1","MYH9","ARHGDIB","ASCL2",
                            "KRT19","ORC6","SNAI1","SOCS3","SPINT2","STK4","TPD52L1",
                            "ACAN","ADAM12","ADAM19","ADAM9","ANGPT4","pseudotime")

k27_CTB_EVT_usm_501_mat <- as.matrix(k27_CTB_EVT_usm_1174_mat)
signature1_a <- USM.k27.norm[c("ACTN4","CEBPB","COL4A2","EFNA1","MYH9","ARHGDIB","ASCL2",
                            "KRT19","ORC6","SNAI1","SOCS3","SPINT2","STK4","TPD52L1",
                            "ACAN","ADAM12","ADAM19","ADAM9","ANGPT4"),]
signature1_a <- as.data.frame(signature1_a)
signature1_a[20,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(signature1_a) <- c("ACTN4","CEBPB","COL4A2","EFNA1","MYH9","ARHGDIB","ASCL2",
                         "KRT19","ORC6","SNAI1","SOCS3","SPINT2","STK4","TPD52L1",
                         "ACAN","ADAM12","ADAM19","ADAM9","ANGPT4","pseudotime")
signature1_a <- t(signature1_a)
signature1_a <- as.data.frame(signature1_a)
signature1_n[signature1_n == 0] <- NA
signature1_a[signature1_a == 0] <- NA
pdf("signature1_k27.pdf", width = 10, height = 5)

pn1 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ACTN4, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa1 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ACTN4, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn1|pa1

pn2 <- ggplot(data=signature1_n, aes(x=pseudotime, y=CEBPB, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa2 <- ggplot(data=signature1_a, aes(x=pseudotime, y=CEBPB, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn2|pa2

pn3 <- ggplot(data=signature1_n, aes(x=pseudotime, y=COL4A2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa3 <- ggplot(data=signature1_a, aes(x=pseudotime, y=COL4A2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn3|pa3

pn4 <- ggplot(data=signature1_n, aes(x=pseudotime, y=EFNA1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa4 <- ggplot(data=signature1_a, aes(x=pseudotime, y=EFNA1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn4|pa4

pn5 <- ggplot(data=signature1_n, aes(x=pseudotime, y=MYH9, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa5 <- ggplot(data=signature1_a, aes(x=pseudotime, y=MYH9, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn5|pa5

pn6 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ARHGDIB, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa6 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ARHGDIB, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn6|pa6

pn7 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ASCL2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa7 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ASCL2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn7|pa7

pn8 <- ggplot(data=signature1_n, aes(x=pseudotime, y=KRT19, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa8 <- ggplot(data=signature1_a, aes(x=pseudotime, y=KRT19, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn8|pa8

pn9 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ORC6, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa9 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ORC6, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn9|pa9

pn10 <- ggplot(data=signature1_n, aes(x=pseudotime, y=SNAI1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa10 <- ggplot(data=signature1_a, aes(x=pseudotime, y=SNAI1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn10|pa10

pn11 <- ggplot(data=signature1_n, aes(x=pseudotime, y=SOCS3, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa11 <- ggplot(data=signature1_a, aes(x=pseudotime, y=SOCS3, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn11|pa11

pn12 <- ggplot(data=signature1_n, aes(x=pseudotime, y=SPINT2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa12 <- ggplot(data=signature1_a, aes(x=pseudotime, y=SPINT2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn12|pa12

pn13 <- ggplot(data=signature1_n, aes(x=pseudotime, y=STK4, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa13 <- ggplot(data=signature1_a, aes(x=pseudotime, y=STK4, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn13|pa13

pn14 <- ggplot(data=signature1_n, aes(x=pseudotime, y=TPD52L1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa14 <- ggplot(data=signature1_a, aes(x=pseudotime, y=TPD52L1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn14|pa14

pn15 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ACAN, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa15 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ACAN, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn15|pa15

pn16 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ADAM12, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa16 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ADAM12, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn16|pa16

pn17 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ADAM19, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa17 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ADAM19, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn17|pa17

pn18 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ADAM9, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa18 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ADAM9, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn18|pa18

pn19 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ANGPT4, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa19 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ANGPT4, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn19|pa19

dev.off()
signature1_n <- CTRL.k27.norm[c("AKAP12","AURKA","AURKB","CD74","DNMT1","DUSP1","EDF1",
                                "EEF1E1","FLT1","GADD45B","GMNN","GTSE1","HMGB1",
                                "HSP90AB1","IQGAP1","ITGA5","KLF6","NPM1","NRN1","NUSAP1","PTMS",
                                "PTPRF","QSOX1","RPL29","RRM1","SERF2","SPINT1","TCF7L2",
                                "TOP2A","WWC3"),]
signature1_n <- as.data.frame(signature1_n)
signature1_n[31,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
signature1_n <- t(signature1_n)
colnames(signature1_n) <- c("AKAP12","AURKA","AURKB","CD74","DNMT1","DUSP1","EDF1",
                            "EEF1E1","FLT1","GADD45B","GMNN","GTSE1","HMGB1",
                            "HSP90AB1","IQGAP1","ITGA5","KLF6","NPM1","NRN1","NUSAP1","PTMS",
                            "PTPRF","QSOX1","RPL29","RRM1","SERF2","SPINT1","TCF7L2",
                            "TOP2A","WWC3","pseudotime")

k27_CTB_EVT_usm_501_mat <- as.matrix(k27_CTB_EVT_usm_1174_mat)
signature1_a <- USM.k27.norm[c("AKAP12","AURKA","AURKB","CD74","DNMT1","DUSP1","EDF1",
                               "EEF1E1","FLT1","GADD45B","GMNN","GTSE1","HMGB1",
                               "HSP90AB1","IQGAP1","ITGA5","KLF6","NPM1","NRN1","NUSAP1","PTMS",
                               "PTPRF","QSOX1","RPL29","RRM1","SERF2","SPINT1","TCF7L2",
                               "TOP2A","WWC3"),]
signature1_a <- as.data.frame(signature1_a)
signature1_a[31,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(signature1_a) <- c("AKAP12","AURKA","AURKB","CD74","DNMT1","DUSP1","EDF1",
                            "EEF1E1","FLT1","GADD45B","GMNN","GTSE1","HMGB1",
                            "HSP90AB1","IQGAP1","ITGA5","KLF6","NPM1","NRN1","NUSAP1","PTMS",
                            "PTPRF","QSOX1","RPL29","RRM1","SERF2","SPINT1","TCF7L2",
                            "TOP2A","WWC3","pseudotime")
signature1_a <- t(signature1_a)
signature1_a <- as.data.frame(signature1_a)
pdf("signature1_k27.pdf", width = 10, height = 5)

pn1 <- ggplot(data=signature1_n, aes(x=pseudotime, y=AKAP12, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa1 <- ggplot(data=signature1_a, aes(x=pseudotime, y=AKAP12, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn1|pa1

pn2 <- ggplot(data=signature1_n, aes(x=pseudotime, y=AURKA, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa2 <- ggplot(data=signature1_a, aes(x=pseudotime, y=AURKA, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn2|pa2

pn3 <- ggplot(data=signature1_n, aes(x=pseudotime, y=AURKB, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa3 <- ggplot(data=signature1_a, aes(x=pseudotime, y=AURKB, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn3|pa3

pn4 <- ggplot(data=signature1_n, aes(x=pseudotime, y=CD74, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa4 <- ggplot(data=signature1_a, aes(x=pseudotime, y=CD74, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn4|pa4

pn5 <- ggplot(data=signature1_n, aes(x=pseudotime, y=DNMT1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa5 <- ggplot(data=signature1_a, aes(x=pseudotime, y=DNMT1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn5|pa5

pn6 <- ggplot(data=signature1_n, aes(x=pseudotime, y=DUSP1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa6 <- ggplot(data=signature1_a, aes(x=pseudotime, y=DUSP1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn6|pa6

pn7 <- ggplot(data=signature1_n, aes(x=pseudotime, y=EDF1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa7 <- ggplot(data=signature1_a, aes(x=pseudotime, y=EDF1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn7|pa7

pn8 <- ggplot(data=signature1_n, aes(x=pseudotime, y=EEF1E1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa8 <- ggplot(data=signature1_a, aes(x=pseudotime, y=EEF1E1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn8|pa8

pn9 <- ggplot(data=signature1_n, aes(x=pseudotime, y=FLT1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa9 <- ggplot(data=signature1_a, aes(x=pseudotime, y=FLT1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn9|pa9

pn10 <- ggplot(data=signature1_n, aes(x=pseudotime, y=GADD45B, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa10 <- ggplot(data=signature1_a, aes(x=pseudotime, y=GADD45B, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn10|pa10

pn11 <- ggplot(data=signature1_n, aes(x=pseudotime, y=GMNN, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa11 <- ggplot(data=signature1_a, aes(x=pseudotime, y=GMNN, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn11|pa11

pn12 <- ggplot(data=signature1_n, aes(x=pseudotime, y=GTSE1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa12 <- ggplot(data=signature1_a, aes(x=pseudotime, y=GTSE1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn12|pa12

pn13 <- ggplot(data=signature1_n, aes(x=pseudotime, y=HMGB1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa13 <- ggplot(data=signature1_a, aes(x=pseudotime, y=HMGB1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn13|pa13

pn14 <- ggplot(data=signature1_n, aes(x=pseudotime, y=HSP90AB1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa14 <- ggplot(data=signature1_a, aes(x=pseudotime, y=HSP90AB1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn14|pa14

pn15 <- ggplot(data=signature1_n, aes(x=pseudotime, y=IQGAP1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa15 <- ggplot(data=signature1_a, aes(x=pseudotime, y=IQGAP1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn15|pa15

pn16 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ITGA5, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa16 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ITGA5, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn16|pa16

pn17 <- ggplot(data=signature1_n, aes(x=pseudotime, y=KLF6, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa17 <- ggplot(data=signature1_a, aes(x=pseudotime, y=KLF6, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn17|pa17

pn18 <- ggplot(data=signature1_n, aes(x=pseudotime, y=NPM1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa18 <- ggplot(data=signature1_a, aes(x=pseudotime, y=NPM1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn18|pa18

pn19 <- ggplot(data=signature1_n, aes(x=pseudotime, y=NRN1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa19 <- ggplot(data=signature1_a, aes(x=pseudotime, y=NRN1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn19|pa19

pn20 <- ggplot(data=signature1_n, aes(x=pseudotime, y=NUSAP1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa20 <- ggplot(data=signature1_a, aes(x=pseudotime, y=NUSAP1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn20|pa20

pn21 <- ggplot(data=signature1_n, aes(x=pseudotime, y=PTMS, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa21 <- ggplot(data=signature1_a, aes(x=pseudotime, y=PTMS, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn21|pa21

pn22 <- ggplot(data=signature1_n, aes(x=pseudotime, y=PTPRF, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa22 <- ggplot(data=signature1_a, aes(x=pseudotime, y=PTPRF, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn22|pa22

pn23 <- ggplot(data=signature1_n, aes(x=pseudotime, y=QSOX1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa23 <- ggplot(data=signature1_a, aes(x=pseudotime, y=QSOX1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn23|pa23

pn24 <- ggplot(data=signature1_n, aes(x=pseudotime, y=RPL29, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa24 <- ggplot(data=signature1_a, aes(x=pseudotime, y=RPL29, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn24|pa24

pn25 <- ggplot(data=signature1_n, aes(x=pseudotime, y=RRM1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa25 <- ggplot(data=signature1_a, aes(x=pseudotime, y=RRM1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn25|pa25

pn26 <- ggplot(data=signature1_n, aes(x=pseudotime, y=SERF2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa26 <- ggplot(data=signature1_a, aes(x=pseudotime, y=SERF2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn26|pa26

pn27 <- ggplot(data=signature1_n, aes(x=pseudotime, y=SPINT1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa27 <- ggplot(data=signature1_a, aes(x=pseudotime, y=SPINT1, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn27|pa27

pn28 <- ggplot(data=signature1_n, aes(x=pseudotime, y=TCF7L2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa28 <- ggplot(data=signature1_a, aes(x=pseudotime, y=TCF7L2, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn28|pa28

pn29 <- ggplot(data=signature1_n, aes(x=pseudotime, y=TOP2A, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa29 <- ggplot(data=signature1_a, aes(x=pseudotime, y=TOP2A, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn29|pa29

pn30 <- ggplot(data=signature1_n, aes(x=pseudotime, y=WWC3, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa30 <- ggplot(data=signature1_a, aes(x=pseudotime, y=WWC3, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn30|pa30


dev.off()
SERPINE1_rescale <- rescale(CTRL.k27.norm["SERPINE1",],to = c(0, 0.5))
signature1_n <- CTRL.k27.norm[c("CD74","ELF3","HTRA4"),]
signature1_n <- as.data.frame(signature1_n)
signature1_n[4,] <- normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(normal.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
signature1_n <- t(signature1_n)
colnames(signature1_n) <- c("TFAP2A","ELF3","HTRA4","pseudotime")

k27_CTB_EVT_usm_501_mat <- as.matrix(k27_CTB_EVT_usm_1174_mat)
signature1_a <- USM.k27.norm[c("TFAP2A","ELF3","HTRA4"),]
signature1_a <- as.data.frame(signature1_a)
signature1_a[4,] <- abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime[order(abortion.evt.cds@principal_graph_aux@listData$UMAP$pseudotime)]
rownames(signature1_a) <- c("TFAP2A","ELF3","HTRA4","pseudotime")
signature1_a <- t(signature1_a)
signature1_a <- as.data.frame(signature1_a)
dim(signature1_a)
write.csv(signature1_n,"signature1_n.csv")
pdf("signature3_k27.pdf", width = 10, height = 5)
library(ggplot2)
pn1 <- ggplot(data=signature1_n, aes(x=pseudotime, y=TFAP2A, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa1 <- ggplot(data=signature1_a, aes(x=pseudotime, y=TFAP2A, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn1|pa1

pn2 <- ggplot(data=signature1_n, aes(x=pseudotime, y=ELF3, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa2 <- ggplot(data=signature1_a, aes(x=pseudotime, y=ELF3, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn2|pa2

pn3 <- ggplot(data=signature1_n, aes(x=pseudotime, y=STAT3, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pa3 <- ggplot(data=signature1_a, aes(x=pseudotime, y=STAT3, color=pseudotime))+
  stat_smooth(method = "loess",span = 0.75)
pn3|pa3
dev.off()

motif_dat <- read.csv("motif_dat.csv",header = T)
motif_dat$logp <- -log10(motif_dat$p.value)
motif_dat$id <- rownames(motif_dat)
motif_dat$x <- "intron" 

pdf("serpine1_motif_var.pdf")
ggplot(motif_dat,aes(x=snp,y=motif_id))+
  geom_point(aes(size=`similarity.score`,
                 color=`logp`))+
  scale_colour_gradient(low="blue",high="yellow")
dev.off()
table(tropho.k27ac$seurat_clusters)
STC_id <- colnames(tropho.k27ac[,which(tropho.k27ac$predict %in% c("STB","STCs"))])
write.table(STC_id,"./bw/STC_id.txt")
FB_id <- colnames(tropho.k27ac[,which(tropho.k27ac$predict %in% c("FB","FBs"))])
write.table(FB_id,"./bw/FB_id.txt")
HC_id <- colnames(tropho.k27ac[,which(tropho.k27ac$predict %in% c("HB","HCs"))])
write.table(HC_id,"./bw/HC_id.txt")
Endo_id <- colnames(tropho.k27ac[,which(tropho.k27ac$predict %in% c("Endo"))])
write.table(Endo_id,"./bw/Endo_id.txt")
Myc_id <- colnames(tropho.k27ac[,which(tropho.k27ac$predict %in% c("Myc","MyCs"))])
write.table(Myc_id,"./bw/Myc_id.txt")
Lyms_id <- colnames(tropho.k27ac[,which(tropho.k27ac$predict %in% c("Lyms","NKs"))])
write.table(Lyms_id,"./bw/Lyms_id.txt")
Ery_id <- colnames(tropho.k27ac[,which(tropho.k27ac$predict %in% c("Ery"))])
write.table(Ery_id,"./bw/Ery_id.txt")
