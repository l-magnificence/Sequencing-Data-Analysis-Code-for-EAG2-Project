setwd("~/scRNA_GBM_integration/Martin_project/")
library(data.table)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(dplyr)

min.max.norm <- function(x){
  ((x-min(x))/(max(x)-min(x)))
}
gmtPathways <- function(gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  pathways
}

### glass datasets
expr<-fread("~/bioinfo_mill/dataset/GLASS/data/gene_tpm_matrix_all_samples.tsv",data.table = F)
rownames(expr)<-expr$Gene_symbol
expr<-expr[,-1]
expr[1:5,1:5]
colnames(expr)<-gsub("\\.","-",colnames(expr))

## select TP and R1
paired<-fread("~/bioinfo_mill/dataset/GLASS/tables/analysis_rnaseq_pairs.csv",data.table = F)
paired$group<-paste(paired$sample_type_a,paired$sample_type_b,sep = "-")
table(paired$group)
paired<-paired[paired$group=="TP-R1",]
table(paired$comparison_type)

paired<-select(paired,c("tumor_pair_barcode","case_barcode","tumor_barcode_a","tumor_barcode_b","sample_type_a","sample_type_b"))
paired_tp<-paired[,c(1:3,5)]
paired_r1<-paired[,c(1:2,4,6)]
colnames(paired_tp)<-c("tumor_pair_barcode","case_barcode","tumor_barcode","sample_type")
colnames(paired_r1)<-c("tumor_pair_barcode","case_barcode","tumor_barcode","sample_type")
paired<-rbind(paired_tp,paired_r1)
paired$sample_barcode<-substr(paired$tumor_barcode,1,15)
paired<-paired[!duplicated(paired$sample_barcode),]
rownames(paired)<-NULL

## merge with clinical information
clinic<-fread("~/bioinfo_mill/dataset/GLASS/tables/clinical_surgeries.csv",data.table = F)
clinic<-merge(paired,clinic,by="sample_barcode",all.x=T)

clinic_idhwt<-clinic[clinic$idh_codel_subtype=="IDHwt",]
table(clinic_idhwt$idh_codel_subtype,clinic_idhwt$histology)

expr_idhwt<-expr[,clinic_idhwt$tumor_barcode]
expr_idhwt<-log2(expr_idhwt+1)
identical(colnames(expr_idhwt),clinic_idhwt$tumor_barcode)

### expression of EAG2/Kvβ2
tmp<-clinic_idhwt
tmp$KCNH5<-t(expr_idhwt["KCNH5",])[,1]
tmp$KCNAB2<-t(expr_idhwt["KCNAB2",])[,1]
tmp[tmp=="TP"]<-"Initial"
tmp[tmp=="R1"]<-"Recurrent"

tmp$eag2_kvbeta2<-(tmp$KCNH5+tmp$KCNAB2)/2
tapply(tmp$eag2_kvbeta2, tmp$sample_type, summary)

p1<-ggplot(tmp ,aes(x = sample_type, y = KCNH5, fill = sample_type))+
  scale_fill_manual(values =c("#F5A200","#A70B00"))+
  geom_violin(alpha=0.25, position = position_dodge(width = .75),size=1,color="NA") +
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.2, alpha = 0.7)+
  theme_classic()+
  geom_line(aes(group = tumor_pair_barcode), color = 'gray', lwd = 0.1)+
  labs(title = "GLASS cohort")+
  ylab("KCNH5 expression") +
  xlab("") +
  stat_compare_means(aes(group=sample_type),method = "wilcox.test",hide.ns = F,label.x = 1.5,#label.y = 2.35,
                     paired = TRUE,label = "p.format")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

p2<-ggplot(tmp ,aes(x = sample_type, y = KCNAB2, fill = sample_type))+
  scale_fill_manual(values =c("#F5A200","#A70B00"))+
  geom_violin(alpha=0.25, position = position_dodge(width = .75),size=1,color="NA") +
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.2, alpha = 0.7)+
  theme_classic()+
  geom_line(aes(group = tumor_pair_barcode), color = 'gray', lwd = 0.1)+
  labs(title = "GLASS cohort")+
  ylab("KCNAB2 expression") +
  xlab("") +
  stat_compare_means(aes(group=sample_type),method = "wilcox.test",hide.ns = F,label.x = 1.5,#label.y = 2.35,
                     paired = TRUE,label = "p.format")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

p3<-ggplot(tmp ,aes(x = sample_type, y = eag2_kvbeta2, fill = sample_type))+
  scale_fill_manual(values =c("#F5A200","#A70B00"))+
  geom_violin(alpha=0.25, position = position_dodge(width = .75),size=1,color="NA") +
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.2, alpha = 0.7)+
  theme_classic()+
  geom_line(aes(group = tumor_pair_barcode), color = 'gray', lwd = 0.1)+
  labs(title = "GLASS cohort")+
  ylab("Average expression of EAG2/Kvβ2") +
  xlab("") +
  stat_compare_means(aes(group=sample_type),method = "wilcox.test",hide.ns = F,label.x = 1.5,#label.y = 2.35,
                     paired = TRUE,label = "p.format")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

tiff('./version2/fig/GLASS_KCNH5_KCNAB2.tiff',units ="in",width = 2.1,height = 3.5,res = 1200,compression ='zip')
p3
dev.off()

### cluster7 and neuron signature
Ctrl_DIP_sample_Markers <- readRDS("~/scRNA_GBM_integration/Martin_project/version2/processed_data/tumor_Markers.rds")
Ctrl_DIP_sample_Markers<-Ctrl_DIP_sample_Markers[Ctrl_DIP_sample_Markers$cluster %in% c("Cluster7"),]
Ctrl_DIP_sample_Markers<-Ctrl_DIP_sample_Markers[Ctrl_DIP_sample_Markers$avg_log2FC>0.5 & Ctrl_DIP_sample_Markers$p_val_adj<0.05,]

GOBP <- gmtPathways("~/software/msigdb_v7.2/c5.go.bp.v7.2.symbols.gmt")
neuron_related<-GOBP[c(grep("AXON",names(GOBP)),grep("NEURON",names(GOBP)),grep("SYNAPSE",names(GOBP)))]
neuron_related<-unique(as.character(unlist(neuron_related)))

gene_list<-list(cluster7=Ctrl_DIP_sample_Markers$gene,NEURON=neuron_related)
ssGSEA<-gsva(as.matrix(expr_idhwt), gene_list,method="ssgsea")
ssGSEA<-as.data.frame(t(ssGSEA))
identical(tmp$tumor_barcode,rownames(ssGSEA))
ssGSEA<-cbind(ssGSEA,tmp)
ssGSEA[,1:2]<-apply(ssGSEA[,1:2],2,min.max.norm)
saveRDS(ssGSEA,file = "./version2/processed_data/glass_ssgsea.rds")

tapply(ssGSEA$cluster7, ssGSEA$sample_type, summary)
tapply(ssGSEA$NEURON, ssGSEA$sample_type, summary)

p4<-ggplot(ssGSEA ,aes(x = sample_type, y = cluster7, fill = sample_type))+
  scale_fill_manual(values =c("#F5A200","#A70B00"))+
  geom_violin(alpha=0.25, position = position_dodge(width = .75),size=1,color="NA") +
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.2, alpha = 0.7)+
  theme_classic()+
  geom_line(aes(group = tumor_pair_barcode), color = 'gray', lwd = 0.1)+
  labs(title = "GLASS cohort")+
  ylab("Cluster7 signatures enrichment") +
  xlab("") +
  stat_compare_means(aes(group=sample_type),method = "wilcox.test",hide.ns = F,label.x = 1.5,#label.y = 2.35,
                     paired = TRUE,label = "p.format")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

p5<-ggplot(ssGSEA ,aes(x = sample_type, y = NEURON, fill = sample_type))+
  scale_fill_manual(values =c("#F5A200","#A70B00"))+
  geom_violin(alpha=0.25, position = position_dodge(width = .75),size=1,color="NA") +
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.2, alpha = 0.7)+
  theme_classic()+
  geom_line(aes(group = tumor_pair_barcode), color = 'gray', lwd = 0.1)+
  labs(title = "GLASS cohort")+
  ylab("Neuronal signatures enrichment") +
  xlab("") +
  stat_compare_means(aes(group=sample_type),method = "wilcox.test",hide.ns = F,label.x = 1.5,#label.y = 2.35,
                     paired = TRUE,label = "p.format")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

p3+p4+p5

tiff('./version2/fig/GLASS_Cluster7.tiff',units ="in",width = 2.1,height = 3.5,res = 1200,compression ='zip')
p4
dev.off()

tiff('./version2/fig/GLASS_neuron.tiff',units ="in",width = 2.1,height = 3.5,res = 1200,compression ='zip')
p5
dev.off()

### heatmap
heatmap_ssGSEA<-ssGSEA
heatmap_ssGSEA[,c(1:2,39:41)]<-scale(heatmap_ssGSEA[,c(1:2,39:41)])

dat_tp<-heatmap_ssGSEA[heatmap_ssGSEA$sample_type=="Initial",]
rownames(dat_tp)<-dat_tp$tumor_pair_barcode
dat_tp2<-dat_tp[,c(1:2,39:41)]
colnames(dat_tp2)<-paste("I",colnames(dat_tp2),sep = "_")

dat_r1<-heatmap_ssGSEA[heatmap_ssGSEA$sample_type=="Recurrent",]
rownames(dat_r1)<-dat_r1$tumor_pair_barcode
identical(rownames(dat_r1),rownames(dat_tp))
dat_r12<-dat_r1[,c(1:2,39:41)]
colnames(dat_r12)<-paste("R",colnames(dat_r12),sep = "_")

dat<-cbind(dat_tp2,dat_r12)
dat<-select(dat,c("I_KCNH5","R_KCNH5","I_KCNAB2","R_KCNAB2","I_eag2_kvbeta2","R_eag2_kvbeta2","I_cluster7","R_cluster7",
                  "I_NEURON","R_NEURON"))
dat<-dat[order(dat$R_eag2_kvbeta2,decreasing = T),]

pheatmap(t(dat[,5:10]),show_colnames = F,cluster_rows = F,cluster_cols = T)

row_anno_dat<-as.data.frame(colnames(dat[,5:10]))
row_anno_dat$Signature<-stringr::str_split(row_anno_dat[,1],"_",simplify = T)[,2]
row_anno_dat$Type<-stringr::str_split(row_anno_dat[,1],"_",simplify = T)[,1]
rownames(row_anno_dat)<-row_anno_dat[,1]
row_anno_dat<-row_anno_dat[,-1]
row_anno_dat[row_anno_dat=="eag2"]<-"EAG2/Kvβ2"
row_anno_dat[row_anno_dat=="NEURON"]<-"Neuronal signature"
row_anno_dat[row_anno_dat=="cluster7"]<-"Cluster7 signature"
row_anno_dat[row_anno_dat=="I"]<-"Initial"
row_anno_dat[row_anno_dat=="R"]<-"Recurrent"

row_ha = rowAnnotation(Signature=row_anno_dat$Signature,Type=row_anno_dat$Type,
  col = list(Signature =structure(names = c("EAG2/Kvβ2","Cluster7 signature","Neuronal signature"), c("#7876B1","#20854E","#EE4C97")),
             Type =structure(names = c("Initial","Recurrent"), c("#F5A200","#A70B00"))),
                       na_col = "white",border = F,
                       gap = unit(3, "points"),
  annotation_name_side = "top")

tiff('./version2/fig/GLASS_heatmap.tiff',units ="in",width = 8.0,height = 3.8,res = 1200,compression ='zip')
Heatmap(t(dat[,5:10]),name="Enrichment",cluster_columns = T,cluster_rows =F,border = T,
        col=colorRampPalette(c(paletteer::paletteer_c("grDevices::Spectral", 30,direction = -1)))(30),
        show_column_names =F,show_row_names = F,left_annotation = row_ha,
        row_names_gp = gpar(fontsize = 11),width  = unit(15, "cm"),
        column_names_gp = gpar(fontsize = 11),height  = unit(6, "cm"),gap = unit(8, "points"),
        #height  = unit(15, "cm"),
        row_split = c(rep(c("A"), 2),rep(c("B"), 2),rep(c("C"), 2)),row_title=NULL,
        column_title = "IDH-wide-type patients in GLASS cohort (n=113)")
dev.off()
