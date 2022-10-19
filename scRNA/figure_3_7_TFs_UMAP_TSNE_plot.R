options(stringsAsFactors = FALSE)
##Load packages
library(ggplot2)
library(plyr)
library(ggsci)
library(ggrepel)

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor_all)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(55)

## calculating the position of cluster labels
get_label_pos <- function(data, emb = "tSNE", group.by="ClusterID") {
  new.data <- data[, c(paste(emb, 1:2, sep = "_"), group.by)]
  colnames(new.data) <- c("x","y","cluster")
  clusters <- names(table(new.data$cluster))
  new.pos <- lapply(clusters, function(i) {
    tmp.data = subset(new.data, cluster == i)
    data.frame(
      x = median(tmp.data$x),
      y = median(tmp.data$y),
      label = i)
  })
  do.call(rbind, new.pos)
}

#read Seurat object and load initial coldata and matrix
target<- readRDS(file = "/mnt/data/chenwei/gongchen/3.seurat_result/final_GC_seurat_0906.rds")
unique(target[[target_group]][,target_group])
ALL_coldata<-target@meta.data
##plot cells
Cells_col0<-colorRampPalette(colors = rev(Cells_col_raw))(55)
length(unique(Cells_col))
barplot(rep(1,55), col=Cells_col)

DimPlot(target, group.by = "final_major_subgroup_brief",cols =rev(Cells_col0))
DimPlot(target, group.by = "final_major_subgroup_brief",cols =Cells_col)

#DimPlot
cell_name="ALLcell_all_genes"
cell.info_target1<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s4_",cell_name,".info.rds"))
col_used<-rev(Cells_col0)
col_used<-Cells_col

#tSNE
Cell_tsne <-ggplot(cell.info_target1, aes(tSNE_1, tSNE_2, color=final_major_subgroup_brief)) + 
  geom_point(size=.1) + 
  geom_label_repel(inherit.aes = F, data = get_label_pos(cell.info_target1, emb = "tSNE",group.by="final_major_subgroup_brief"), aes(x,y,label=label),size=3, box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))+
  # geom_text_repel(inherit.aes = F, data = get_label_pos(cell.info_target, emb = "tSNE",group.by="sub_cluster_brief"), aes(x,y,label=label),size=3, box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))+
  # geom_text(inherit.aes = F, data = get_label_pos(cell.info_target, emb = "tSNE",group.by="sub_cluster_brief"), aes(x,y,label=label), size=3) + 
  labs(title ="All Cells")+ scale_color_manual(values=col_used)+
  theme_bw(base_size = 15) + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"))

#UMAP
Cell_umap <-ggplot(cell.info_target1, aes(UMAP_1, UMAP_2, color=final_major_subgroup_brief)) + geom_point(size=.1) + 
  geom_label_repel(inherit.aes = F, data = get_label_pos(cell.info_target1, emb = "UMAP",group.by="final_major_subgroup_brief"), aes(x,y,label=label),size=3, box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))+
  # geom_text_repel(inherit.aes = F, data = get_label_pos(cell.info_target, emb = "UMAP",group.by="sub_cluster_brief"), aes(x,y,label=label),size=3, box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))+
  # geom_text(inherit.aes = F, data = get_label_pos(cell.info_target, emb = "UMAP",group.by="sub_cluster_brief"), aes(x,y,label=label), size=3) + 
  labs(title ="All Cells")+scale_color_manual(values=col_used)+
  theme_bw(base_size = 15) + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"))

saveRDS(Cell_umap,paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/",cell_name,"Cell_umap.rds"))
saveRDS(Cell_tsne,paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/",cell_name,"Cell_tsne.rds"))
##

Cell_umap<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/",cell_name,"Cell_umap.rds"))
Cell_tsne<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/",cell_name,"Cell_tsne.rds"))

#Cell_tsne+facet_grid(.~Age_group)
Cell_tsne1<-Cell_tsne+NoLegend()+facet_wrap(.~Treat,ncol =2 )
Cell_umap1<-Cell_umap+NoLegend()+facet_wrap(.~Treat,ncol =2 )
Seurate_plot<-DimPlot(target, group.by = "final_major_subgroup_brief",split.by = "Treat", cols =col_used)
Seurate_plot
Cell_umap1
Cell_tsne1
head(cell.info_target1)
#for target tissue
tSNEPlot.sample<- function(data, sample_use, topn) {
  # data=cell.info_target;tissue.use="Base";topn=4
  clusters.use <- table(subset(data, sample_code == sample_use)$Treat)
  if(topn > length(clusters.use)) {
    topn <- length(clusters.use)
  }
  clusters.use <- names(head(sort(clusters.use, decreasing = T), n=topn))
  data$new.cluster <- ifelse(data$sample_code == sample_use & data$Treat %in% clusters.use, as.character(data$final_major_subgroup_brief), "Others")
  data$new.cluster <- factor(data$new.cluster, levels = c(setdiff(names(table(data$new.cluster)), "Others"), "Others"))
  data$pt.size <- ifelse(data$new.cluster == "Others", 0.05, 0.2)
  n_length<-length(levels(data$new.cluster))
  ggplot(data, aes(tSNE_1, tSNE_2, color=new.cluster)) + 
    geom_point(size=data$pt.size,alpha=0.5) + scale_color_manual(values=c(col_used[1:(n_length-1)],"grey"))+
    #scale_color_manual(values = c(pal_d3("category10")(topn), "grey")) + 
    guides(colour = guide_legend(keyheight = .7, keywidth = .1,override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + ggtitle(sample_use) + 
    theme(#legend.justification = c(0,0),
      legend.title = element_blank(),
      # legend.position = c(0,0),
      # legend.key = element_rect(fill = alpha("white", 0)),
      # legend.background = element_rect(fill=alpha('white', 0)),
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      axis.line = element_line(color="black"),
      plot.title = element_text(hjust = .5, face = "bold"))
}

P1<-tSNEPlot.sample(cell.info_target1, "CTRL_1", 2)
P2<-tSNEPlot.sample(cell.info_target1, "CTRL_2", 2)
P3<-tSNEPlot.sample(cell.info_target1, "CTRL_3", 2)
P4<-tSNEPlot.sample(cell.info_target1, "CTRL_4", 2)
P5<-tSNEPlot.sample(cell.info_target1, "CTRL_5", 2)
P6<-tSNEPlot.sample(cell.info_target1, "Arrest_1", 2)
P7<-tSNEPlot.sample(cell.info_target1, "Arrest_2", 2)
P8<-tSNEPlot.sample(cell.info_target1, "Arrest_3", 2)
P9<-tSNEPlot.sample(cell.info_target1, "Arrest_4", 2)
P1+P2+P3+P4+P5
P6+P7+P8+P9

#for target cell group
tSNEPlot.cellType <- function(data, cell_type, topn) {
  # data=cell.info_target;cell_type="HCs"; topn=12
  tissues_use <- table(as.character(subset(data, final_major_subgroup_brief == cell_type)$sample_code))
  if(topn > length(tissues_use)) {
    topn <- length(tissues_use)
  }
  tissues_use <- names(head(sort(tissues_use, decreasing = T), n=topn))
  data$new_tissue <- ifelse(data$final_major_subgroup_brief == cell_type & data$sample_code %in% tissues_use, as.character(data$sample_code), "Others")
  data$new_tissue <- factor(data$new_tissue, levels = c(setdiff(names(table(data$sample_code)), "Others"), "Others"))
  data$pt_size <- ifelse(data$new_tissue == "Others", 0.05, 0.2)
  
  ggplot(data, aes(tSNE_1, tSNE_2, color=new_tissue)) + 
    geom_point(size=data$pt_size) + 
    scale_color_manual(values = c(pal_d3("category10")(topn), "grey")) + 
    guides(colour = guide_legend(keyheight = .7,keywidth = .1, override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + 
    ggtitle(cell_type) + 
    theme(#legend.justification = c(0,0), legend.position = c(0,0),
      legend.title = element_blank(),
      #legend.key = element_rect(fill = alpha("white", 0)),
      #legend.background = element_rect(fill=alpha('white', 0)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color="black"),
      plot.title = element_text(hjust = .5, face = "bold")
    )
}
table(cell.info_target$final_major_subgroup_brief)
P1<-tSNEPlot.cellType(cell.info_target1, "HCs", 10)
P2<-tSNEPlot.cellType(cell.info_target1, "Masts",10)
P3<-tSNEPlot.cellType(cell.info_target1, "Bs", 10)
P4<-tSNEPlot.cellType(cell.info_target1, "Epi", 10)
P5<-tSNEPlot.cellType(cell.info_target1, "Ery", 10)
P6<-tSNEPlot.cellType(cell.info_target1, "Ts", 10)
P7<-tSNEPlot.cellType(cell.info_target1, "Endo",10)
P8<-tSNEPlot.cellType(cell.info_target1, "NKs", 10)
P9<-tSNEPlot.cellType(cell.info_target1, "Troph_mix", 10)

P1+P2+P3+P4+P5+P6+P7+P8+P9

table(cell.info_target$final_major_subgroup_brief)
P1<-tSNEPlot.cellType(cell.info_target1, "CTBs_0", 10)
P2<-tSNEPlot.cellType(cell.info_target1, "CTBs_1",10)
P3<-tSNEPlot.cellType(cell.info_target1, "CTBs_2", 10)
P4<-tSNEPlot.cellType(cell.info_target1, "CTBs_3", 10)
P5<-tSNEPlot.cellType(cell.info_target1, "EVTs_1", 10)
P6<-tSNEPlot.cellType(cell.info_target1, "EVTs_2", 10)
P7<-tSNEPlot.cellType(cell.info_target1, "EVTs_3",10)
P8<-tSNEPlot.cellType(cell.info_target1, "STBs_1", 10)
P9<-tSNEPlot.cellType(cell.info_target1, "STBs_2", 10)
P10<-tSNEPlot.cellType(cell.info_target1, "STCs", 10)
P11<-tSNEPlot.cellType(cell.info_target1, "FBs", 10)
P12<-tSNEPlot.cellType(cell.info_target1, "MyCs", 10)

P1+P2+P3+P4+P5+P6+P7+P8+P9+P10+P11+P12
##for specific regulon for special cell types
###=========
##待写


##step 3 ：： Quantifying cell-type specificity score
#在通过SCENIC分析得到了regulon的结果后，作者定义了RSS（Regulon Specificity Score）来寻找细胞类型特异的转录调控网络。
##作者首先定义了RAS（Regulon Activity Scores）的概率分布
## Load packages
library(data.table)
library(pbapply)
library(plyr)
# philentropy包用于计算JS散度，关于这个package的详细介绍在此：# http://blog.fens.me/r-entropy/
library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)

## 读入RAS矩阵
##rasMat(rds可直接读取): rows=cellID, cols=regulon, values=Regulon Acticity Score
rasMat <- fread(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".AUCell.txt"), sep = "\t",  header = T, data.table = F)
rownames(rasMat) <- rasMat$V1
colnames(rasMat) <- sub("(+)", "", colnames(rasMat), fixed = T)
rasMat <- rasMat[, -1]
saveRDS(rasMat, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".rasMat.rds"))

## 读入细胞类型矩阵
cell.info_target<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s4_",cell_name,".info.rds"))
#cell.info_target<-ALL_coldata[which(ALL_coldata$sub_cluster_brief == cell_name),]
#cell.info_target <- cbind(cell.info_target, emb.tsne[rownames(cell.info_target), ])
#cell.info_target <- cbind(cell.info_target, emb.umap[rownames(cell.info_target), ])
head(cell.info_target)

#compare_group<-"sub_cluster_brief_Age"
compare_group<-"sub_cluster_brief"
cell_types <- names(table(cell.info_target[compare_group]))
ctMat <- lapply(cell_types, function(i) {
  as.numeric(cell.info_target[compare_group] == i)
})
ctMat <- do.call(cbind, ctMat)
colnames(ctMat) <- cell_types
rownames(ctMat) <- rownames(cell.info_target)

## 计算RSS矩阵(Regulon Specificity Score)(rds可直接读取): 当RSS值越大，Regulon R相对于cell type C就越特异。
#input: ctMat ; rasMat
#/mnt/data/chenwei/gongchen/4.SCENIC/loom_result: rssMat, rows=regulon, cols=cell.type, values=Regulon Specific Score

##rasMat中不能出现某个基因列完全为0,出现这种情况则需要过滤
#rasMat<-rasMat[,colSums(rasMat)>0]
rasMat<-readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".rasMat.rds"))

rssMat <- pblapply(colnames(rasMat), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
    #两个概率分布(RAS和细胞类群)的JS散度（Jensen-Shannon Divergence） #empirical:做归一化
  })
})
rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(rasMat)
colnames(rssMat) <- colnames(ctMat)
#for cell type_specific
saveRDS(rssMat, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s5_",cell_name,"_whole_cell_type.rssMat.rds"))

##plot specific TFs for each cell types
rssMat <- readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s5_",cell_name,"_whole_cell_type.rssMat.rds"))
binMat <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".binary_mtx.txt"), sep = "\t", header = T, row.names = 1, check.names = FALSE)
colnames(binMat) <- sub("(+)", "", colnames(binMat), fixed = T)

## Regulon Rank Plot
source("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/plotRegulonRank.R")
colnames(rssMat)
PlotRegulonRank(rssMat, "HBs")+PlotRegulonRank(rssMat,"F_Nai_CD4_Ts")
#PlotRegulonRank(rssMat, "HBs_old")+PlotRegulonRank(rssMat, "HBs_young")
#PlotRegulonRank(rssMat, "HBs_old",10)+PlotRegulonRank(rssMat, "HBs_young",10)

## DimPlot
### Load data
cell.info_target2 <- cbind(cell.info_target, binMat[rownames(cell.info_target), ])
cell.info_target2[1:4,1:4]
source("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/DimPlot.R")

fig2Plot <- function(group,cell.type, regulon) {
  p.list <- list(
    PlotRegulonRank(rssMat, cell.type),
    DimPlot(cell.info_target2,group=group, cell.type = cell.type),
    DimPlot(cell.info_target2,group=group, regulon = regulon)
  )
  cowplot::plot_grid(plotlist = p.list, ncol = 3,rel_widths = c(3,5,5))
}

plot_group<-"Ery"   
plot_TFs<- names(which(rssMat[,plot_group] == max(rssMat[,plot_group])))
DimPlot(cell.info_target2,group=compare_group,cell.type =plot_group)+DimPlot(cell.info_target2,group=compare_group,regulon = plot_TFs)
#DimPlot(cell.info_target2,group="sub_cluster_brief_Age",cell.type = "HBs_old")+DimPlot(cell.info_target2,group="sub_cluster_brief_Age",regulon = "CEBPB")
DimPlot(cell.info_target2,group="sub_cluster_brief",cell.type = "B_Bs")+DimPlot(cell.info_target2,group="sub_cluster_brief",regulon ="TCF3")

## 
p_tsne_list<-list();cell_target_name<-c()
for ( cell_target in colnames(rssMat)){
  print(cell_target)
  plot_group<-cell_target  
  plot_TFs<- names(which(rssMat[,plot_group] == max(rssMat[,plot_group])))
  p_combin<-fig2Plot(compare_group,plot_group,plot_TFs)
  p_tsne_list<-c(p_tsne_list,list(p_combin))
  cell_target_name<-c(cell_target_name,cell_target)
}
length(p_tsne_list);length(cell_target_name)#55 55
names(p_tsne_list)<-cell_target_name
p_tsne_list[["HBs"]]

saveRDS(p_tsne_list, paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s5_",cell_name,"_cell_TFs_tsne_plot.rds"))
p_tsne_list <- readRDS(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s5_",cell_name,"_cell_TFs_tsne_plot.rds"))

#p_tsne_group <- ggpubr::ggarrange(p_tsne_list[[1]],p_tsne_list[[2]],p_tsne_list[[3]],p_tsne_list[[4]],p_tsne_list[[5]],
#                                 p_tsne_list[[6]],p_tsne_list[[7]],p_tsne_list[[8]],p_tsne_list[[9]],p_tsne_list[[10]],
#                                 #labels = c('A', 'B', 'C', 'D'),  font.label = list(color = 'black'),
#                                 nrow = 5, ncol = 2)
#annotate_figure(p_pie_group, top=text_grob(paste0("In_Vivo-In_vitro_DMRs_myDiff50p5",DMR_type), color = "black",face = "bold", size=12))
#fig2Plot("sub_cluster_brief","Ery","IRF9")

##Functional validation:: 作者对找到的cell specific regulon利用其它的方法进行了验证。

#（1）SEEK analysis，2000+ GEO datasets. 
## 检验  a）regulon genes是否是共表达的，b）regulon基因是否和给定细胞类型有相关性
##（2）CoCiter analysis，文献挖掘，检验一组基因是否和某个term（本文使用的是细胞类型）有相关性。

### 利用SEEK对regulon进行验证# http://seek.princeton.edu/ # using SEEK to analysis.
#SEEK只能最多接收100个基因作为input。如果我们的regulon超过这个数目，我们在需要验证的regulon中随机挑选100个基因进行分析
regulon_file <- read.table(paste0("/mnt/data/chenwei/gongchen/4.SCENIC/loom_result/",cell_name,"/s3_",cell_name,".regulons.txt"), sep = "\t", header = F)
head(regulon_file)
length(unlist(strsplit(regulon_file[grep("STAT4",regulon_file$V1),]$V3,split=',')))
unlist(strsplit(regulon_file[grep("MITF",regulon_file$V1),]$V3,split=','))[1:100]

source("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/SeekPlot.R")
ELF1.seek <- read.csv("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/seek_result/ELF1_seek.txt", sep = "\t", header = T)
MITF.seek <- read.csv("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/seek_result/MITF_seek.txt", sep = "\t", header = T)

p3.1<-SeekPlot(ELF1.seek, key.words_in= "(placent)|(Hofbauer)|( aging)|( ageing)|( age)",key.words_out="(agent)")
p4.1<-SeekPlot(MITF.seek, key.words_in= "(placent)|(Hofbauer)|( aging)|( ageing)|( age)",key.words_out="(agent)")

#p1.2 <- cowplot::plot_grid(p1, p1.1, rel_widths = c(13,3))
#p2.2 <- cowplot::plot_grid(p2, p2.1, rel_widths = c(13,3))
p3.2 <- cowplot::plot_grid(p3, p3.1, rel_widths = c(13,3))
p4.2 <- cowplot::plot_grid(p4, p4.1, rel_widths = c(13,3))
cowplot::plot_grid(p3.2, p4.2, ncol = 1)

### save plot
#if(!dir.exists("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/figure")) {
#  dir.create("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/figure")
#}
#pdf("/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/figures/Figure2.pdf", width = 16, height = 20, units = "in")
#cowplot::plot_grid(p1.2, p2.2, p3.2, p4.2, ncol = 1)
#dev.off()


sessionInfo()

