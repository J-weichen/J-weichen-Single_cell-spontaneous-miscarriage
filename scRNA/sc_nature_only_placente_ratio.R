single_cell_placent_all<-readRDS(file = "/mnt/data/chenwei/qinmen_BR/nature_single_cell_placent_all.rds")
Metadate_plot<-single_cell_placent_all@meta.data
unique(Metadate_plot$seurat_cluster)

Metadate_plot<-single_cell_placent_all@meta.data
Metadate_plot2<-Metadate_plot[which(Metadate_plot$location =="Placenta"),]
table(Metadate_plot2$annotation)
Metadate_plot2$major_group_brief<-as.character(Metadate_plot2$annotation)

Metadate_plot2[which(Metadate_plot2$annotation %in% c("DC1","DC2","dM1","dM2","dM3","HB","MO")),]$major_group_brief<- "Myeloid_Cell"
Metadate_plot2[which(Metadate_plot2$annotation %in% c("dNK p","dNK1","dNK2","NK CD16+","Tcells","Plasma","ILC3","Granulocytes")),]$major_group_brief<- "T_NK_B_Cell"
Metadate_plot2[which(Metadate_plot2$annotation %in% c("dS1","dS2")),]$major_group_brief<- "Decidual_stromal_cells"
Metadate_plot2[which(Metadate_plot2$annotation %in% c("Endo (f)","Endo (m)","Endo L")),]$major_group_brief<- "Endothelial_Cell"
Metadate_plot2[which(Metadate_plot2$annotation %in% c("fFB1","fFB2")),]$major_group_brief<- "Fibroblasts"
Metadate_plot2[which(Metadate_plot2$annotation %in% c("SCT","VCT","EVT","Epi1","Epi2")),]$major_group_brief<- "Trophoblast"
table(Metadate_plot2$major_group_brief)
table(Metadate_plot2$annotation)

Metadate_plot2$major_group_brief<-factor(Metadate_plot2$major_group_brief,levels = c("Trophoblast","Fibroblasts","Decidual_stromal_cells","Endothelial_Cell","Myeloid_Cell","T_NK_B_Cell"))

dim(Metadate_plot2)
type_pieplot<-ggplot(as.data.frame(table(Metadate_plot2$major_group_brief)),aes(x="",y=Freq,fill=Var1))+ 
  geom_bar(stat="identity")+ 
  coord_polar("y",start=1) + 
  # geom_text(aes(x=1,label=percent(Freq/sum(Freq))), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+ labs(title = "nature_placenta")+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(), legend.title = element_blank())+ 
  scale_fill_manual(values=ppCor_all)
type_pieplot


cell.prop1<-as.data.frame(prop.table(table(Metadate_plot2$major_group_brief)))
cell.prop1$group<-"nature_placenta"
colnames(cell.prop1)<-c("Cell_type","proportion","group")
type_barplot<-ggplot(cell.prop1,aes(group,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=ppCor_all2)+coord_flip() 
type_pieplot/type_barplot
type_pieplot+type_barplot
