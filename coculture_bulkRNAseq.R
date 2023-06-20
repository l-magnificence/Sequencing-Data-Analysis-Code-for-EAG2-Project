setwd("~/scRNA_GBM_integration/Martin_project/")
library(data.table)
library(ggplot2)
library(ggpubr)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(GseaVis)

###------------tumor cell------------
expr<-fread("/export2/liuhw/Martin_project_bulkRNAseq/human_RNA_output/gene_tpm_matrix.csv",data.table = F)
tmp<-dplyr::select(expr,"gene_id")
tmp$ENSEMBL<-stringr::str_split(tmp$gene_id,"\\|",simplify = T)[,1]
tmp$SYMBOL<-stringr::str_split(tmp$gene_id,"\\|",simplify = T)[,2]
expr$gene_id<-tmp$SYMBOL
expr<-aggregate(x = expr[,2:(ncol(expr))],
                by = list(expr$gene_id),
                FUN = mean)
row.names(expr)<-expr[,1]
expr<-expr[,-1]
##log2(tpm+1)
exprSet<-log2(expr+1)

##raw count
rawcount<-fread("/export2/liuhw/Martin_project_bulkRNAseq/human_RNA_output/gene_count_matrix.csv",data.table = F)
identical(rawcount$gene_id,tmp$gene_id)
rawcount$gene_id<-tmp$SYMBOL
rawcount<-aggregate(x = rawcount[,2:(ncol(rawcount))],
                    by = list(rawcount$gene_id),
                    FUN = mean)
row.names(rawcount)<-rawcount[,1]
rawcount<-rawcount[,-1]

#meta data
phe<-as.data.frame(colnames(exprSet))
colnames(phe)<-"sample_id"
phe$Treatment<-stringr::str_split(phe[,1],"_",simplify = T)[,1]
phe$Replicate<-stringr::str_split(phe[,1],"_",simplify = T)[,3]

saveRDS(phe,file = "./version2/processed_data/coculture_RNAseq/coculture_tumor_phe.rds")
saveRDS(exprSet,file = "./version2/processed_data/coculture_RNAseq/coculture_tumor_expr_logtpm.rds")
saveRDS(rawcount,file = "./version2/processed_data/coculture_RNAseq/coculture_tumor_expr_rawcount.rds")
rm(list=ls())

phe <- readRDS("~/scRNA_GBM_integration/Martin_project/version2/processed_data/coculture_RNAseq/coculture_tumor_phe.rds")
expr <- readRDS("~/scRNA_GBM_integration/Martin_project/version2/processed_data/coculture_RNAseq/coculture_tumor_expr_logtpm.rds")
count <- readRDS("~/scRNA_GBM_integration/Martin_project/version2/processed_data/coculture_RNAseq/coculture_tumor_expr_rawcount.rds")

# PCA
pca <- prcomp(t(expr), scale=F)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
identical(phe$sample_id,rownames(pca.data))
pca.data<-cbind(pca.data,phe)
pca.var <- pca$sdev^2  
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  

tiff('./version2/fig/coculture_RNAseq/tumor_pca.tiff',units ="in",width = 4.2,height = 2.8,res = 1200,compression ='zip')
ggplot(pca.data, aes(x=X,y=Y,color=Treatment),size=5) +
  geom_point(size=3) +#stat_ellipse(level = 0.95, show.legend = F,geom = "polygon",alpha = 1/5, aes(fill = patient))+
  scale_colour_manual(values =c("#F5A200","#A70B00"))+
  xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep="")) +
  ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep="")) +
  theme_classic()+
  labs(title = "Tumor cell in co-culture condition")+
  #ggsci::scale_color_lancet()+
  theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.8, "cm"),legend.text = element_text(size=15),legend.title =element_text(size=15))
dev.off()

# DEG
group<-phe$Treatment
group_list=factor(group,levels = c("K59","K90"))
colData <- data.frame(row.names=colnames(count),group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = round(count),
                              colData = colData,design = ~ group_list)
dds2 <- DESeq(dds)
alldiff <- results(dds2,contrast=c("group_list","K90","K59"))
alldiff <- as.data.frame(alldiff[order(alldiff$log2FoldChange,decreasing = T),])  
alldiff = na.omit(alldiff)
alldiff$significant<-ifelse(alldiff$log2FoldChange>0.5 &alldiff$padj<0.05,"Up in K90",
                                ifelse(alldiff$log2FoldChange< -0.5 &alldiff$padj<0.05,"Down in K90","Usignificant"))
table(alldiff$significant)
saveRDS(alldiff,file = "./version2/processed_data/coculture_RNAseq/coculture_tumor_alldiff.rds")

# GO KEGG
DEG<-alldiff[alldiff$significant %in%c("Up in K90","Down in K90"),]
ids <- bitr(row.names(DEG), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db",drop = TRUE)
DEG$gene<-rownames(DEG)
DEG_list=merge(DEG,ids,by.x='gene',by.y='SYMBOL')
DEG_list<-DEG_list[order(DEG_list$log2FoldChange,decreasing = T),]
table(DEG_list$significant)

##gseGO
geneList<-DEG_list$log2FoldChange
names(geneList)<-DEG_list$ENTREZID
gseGO <- gseGO(geneList = geneList, 
               OrgDb = org.Hs.eg.db, 
               keyType = "ENTREZID", 
               ont="BP",
               pvalueCutoff =0.05)
gseGO <-setReadable(gseGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
saveRDS(gseGO,"./version2/processed_data/coculture_RNAseq/coculture_tumor_gseGO.rds")

gseGO <- readRDS("~/scRNA_GBM_integration/Martin_project/version2/processed_data/coculture_RNAseq/coculture_tumor_gseGO.rds")
tmp<-as.data.frame(gseGO)
write.csv(tmp,"./version2/processed_data/coculture_RNAseq/coculture_tumor_gseGO_all.csv")

go_id<-tmp[tmp$NES>0 & tmp$qvalue<0.05,"ID"]
mat <- GO_similarity(go_id,
                    ont =  "BP",db = 'org.Hs.eg.db',measure = "Rel")
saveRDS(mat,"./version2/processed_data/coculture_RNAseq/coculture_tumor_upk90_gseGO_similarity_mat.rds")

tiff('./version2/fig/coculture_RNAseq/tumor_upk90_gseGO.tiff',units ="in",width = 10.5,height = 6.7,res = 1200,compression ='zip')
df <- simplifyGO(mat)
dev.off()

saveRDS(df,"./version2/processed_data/coculture_RNAseq/coculture_tumor_upk90_gseGO_similarity_cluster.rds")

go_id<-tmp[tmp$NES<0 & tmp$qvalue<0.05,"ID"]
mat2 <- GO_similarity(go_id,
                     ont =  "BP",db = 'org.Hs.eg.db',measure = "Rel")
saveRDS(mat2,"./version2/processed_data/coculture_RNAseq/coculture_tumor_downk90_gseGO_similarity_mat.rds")

tiff('./version2/fig/coculture_RNAseq/tumor_downk90_gseGO.tiff',units ="in",width = 10.5,height = 6.7,res = 1200,compression ='zip')
df2 <- simplifyGO(mat2)
dev.off()

saveRDS(df2,"./version2/processed_data/coculture_RNAseq/coculture_tumor_downk90_gseGO_similarity_cluster.rds")

## select
GO_select<-tmp[c(grep("programmed cell death",tmp$Description),grep("cell death",tmp$Description),grep("neuron death",tmp$Description),
                 grep("apoptotic process",tmp$Description),
                 grep("cell cycle",tmp$Description),grep("division",tmp$Description),grep("mitotic",tmp$Description),
                 grep("repair",tmp$Description),grep("replication",tmp$Description),
                 grep("chromosome",tmp$Description),grep("microtubule",tmp$Description)),]
GO_select<-GO_select[-c(grep("negative",GO_select$Description),
                 grep("positive",GO_select$Description),
                 grep("meio",GO_select$Description),
                 grep("regulation",GO_select$Description),
                 grep("epithelial",GO_select$Description)),]
GO_select<-rbind(rbind(GO_select,tmp[c(grep("negative regulation of MAPK cascade",tmp$Description),grep("p38MAPK cascade",tmp$Description),
                       grep("negative regulation of ERK1 and ERK2 cascade",tmp$Description),
                       grep("Wnt signaling pathway",tmp$Description),
                       grep("peptide",tmp$Description),grep("drug",tmp$Description),
                       grep("stress",tmp$Description),grep("radiation",tmp$Description)),]))
tmp<-GO_select[!duplicated(GO_select$Description),]

GO_select<-tmp[c(grep("peptide",tmp$Description),grep("drug",tmp$Description),
                 grep("programmed cell death",tmp$Description),grep("cell death",tmp$Description),grep("apoptotic process",tmp$Description),
                 grep("cell cycle",tmp$Description),grep("division",tmp$Description),grep("mitotic",tmp$Description),
                 grep("repair",tmp$Description),grep("replication",tmp$Description),
                 grep("microtubule",tmp$Description),
                 grep("p38MAPK cascade",tmp$Description),grep("negative regulation of MAPK cascade",tmp$Description),
                 grep("negative regulation of ERK1 and ERK2 cascade",tmp$Description),
                 grep("oxidative",tmp$Description)),]
GO_select<-GO_select[!duplicated(GO_select$Description),]
write.csv(GO_select,"./version2/processed_data/coculture_RNAseq/coculture_tumor_gseGO_select.csv")

GO_select$Description<-Hmisc::capitalize(GO_select$Description)
go_list<-list("Group1"=data.frame(GO_select[1:5,"Description"], "col"="#BC3C29"),
              "Group2"=data.frame(GO_select[6:8,"Description"], "col"="#80796B"),
              "Group3"=data.frame(GO_select[9:20,"Description"], "col"="#EE4C97"),
              "Group4"=data.frame(GO_select[21:26,"Description"], "col"="#0072B5"),
              "Group5"=data.frame(GO_select[27:32,"Description"], "col"="#E18727"),
              "Group6"=data.frame(GO_select[33:36,"Description"], "col"="#20854E"),
              "Group7"=data.frame(GO_select[37:38,"Description"], "col"="#7876B1"))
              
gene_select<-data.frame()
for (i in 1:nrow(GO_select)) {
  gene_select<-rbind(gene_select,data.frame("gene"=c(stringr::str_split(GO_select[i,"core_enrichment"],"/",simplify = F)[[1]]),
                                            "GO"=GO_select[i,"Description"]))
}
gene_select2<-data.frame()
for (i in 1:length(go_list)) {
  gene_select2<-rbind(gene_select2,data.frame("gene"=unique(gene_select[gene_select$GO %in% go_list[[i]][,1],"gene"]),
                                             "Group"=names(go_list[i])))
}
gene_select2<-gene_select2[!gene_select2$gene %in% c("SMC5","RADX","CBX8","BACH1","REV3L","UBC","ID3"),]
rownames(gene_select2)<-NULL

gap<-as.data.frame(table(gene_select2$Group))
gap$split<-1
for (i in 1:nrow(gap)) {
  gap[i,"split"]<-sum(gap[1:i,"Freq"])
}

annotation_colors = list(Treatment=c("K59"="#F5A200","K90"="#A70B00"))
# genes_to_show = c('ATF3', 'GAS1', 'TNFRSF9', 'BCL3', 'RHOB', 
#                   "WNT9A","WNT11","DUSP10","DUSP1","RGS2",
#                   "IGFBP1","IGFBP4","MAP3K8","NFKBIA","JUND",
#                   "AURKA","AURKB",'BUB1B', "CCNA2","TPX2",
#                   "EME1","BRCA1","AUNIP","EXO5","GLI1",
#                   "MKI67","CDC20","HJURP",
#                   "INCENP","SPDL1","TACC3","KIF14","BORA")
# 
# genes_to_show_idnex<-c(rownames(gene_select2[gene_select2$Group=="Group1",])[match(genes_to_show[1:5],gene_select2[gene_select2$Group=="Group1","gene"])],
#                        rownames(gene_select2[gene_select2$Group=="Group2",])[match(genes_to_show[6:10],gene_select2[gene_select2$Group=="Group2","gene"])],
#                        rownames(gene_select2[gene_select2$Group=="Group3",])[match(genes_to_show[11:15],gene_select2[gene_select2$Group=="Group3","gene"])],
#                        rownames(gene_select2[gene_select2$Group=="Group4",])[match(genes_to_show[16:20],gene_select2[gene_select2$Group=="Group4","gene"])],
#                        rownames(gene_select2[gene_select2$Group=="Group5",])[match(genes_to_show[21:25],gene_select2[gene_select2$Group=="Group5","gene"])],
#                        rownames(gene_select2[gene_select2$Group=="Group6",])[match(genes_to_show[26:28],gene_select2[gene_select2$Group=="Group6","gene"])],
#                        rownames(gene_select2[gene_select2$Group=="Group7",])[match(genes_to_show[29:33],gene_select2[gene_select2$Group=="Group7","gene"])])

tiff('./version2/fig/coculture_RNAseq/tumor_gseGO_gene_heatmap2.tiff',units ="in",width = 8.0,height = 10.0,res = 1200,compression ='zip')
ComplexHeatmap::pheatmap(as.matrix(expr[gene_select2$gene,]),scale = "row",cluster_rows = F,cluster_cols = T,
                         show_colnames = F,show_rownames = F,treeheight_col =20, name = "Expression",
                         gaps_row = gap$split,
                         col=colorRampPalette(c(paletteer::paletteer_c("grDevices::Spectral", 30,direction = -1)))(30),
                         annotation_col = dplyr::select(phe,"Treatment"),
                         annotation_colors = annotation_colors,
                         # left_annotation = rowAnnotation(link = anno_mark(at = as.numeric(genes_to_show_idnex), 
                         #                                                  side="left",
                         #                                                  labels_gp = gpar(fontsize =8),
                         #                                                 labels = genes_to_show)),
                         right_annotation = rowAnnotation(
                           textbox = anno_textbox(gene_select2$Group, go_list,
                                                  background_gp = gpar(fill = "#EEEEEE", col = "white"),
                                                  word_wrap = TRUE, add_new_line = TRUE)))
dev.off()

GO_select$Description<-factor(GO_select$Description,levels=unique(GO_select$Description))
GO_select$Description<-forcats::fct_rev(GO_select$Description)

tiff('./version2/fig/coculture_RNAseq/tumor_gseGO_plot2.tiff',units ="in",width = 6.2,height = 9.0,res = 1200,compression ='zip')
ggplot(data = GO_select, aes(x =NES , y = Description)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_segment(aes(xend=0, yend = Description,color=qvalue))+
  geom_point(aes(color=qvalue, size = abs(NES)))+
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2","#E18727"),
                        guide=guide_colorbar(reverse=TRUE,order=1)) +
  #scale_size_continuous(range=c(2, 10)) +
  ylab("") + xlab("Normalized Enrichment Score") +
  labs(title = "GSEA GO BP in tumor")+
  ggplot2::annotate("rect",
                    ymin = c(0.5,2.5,6.5,12.5,18.5,30.5,33.5),
                    ymax = c(2.5,6.5,12.5,18.5,30.5,33.5,38.5),
                    xmin = -4, xmax = 2, alpha = 0.2,
                    fill = c("#7876B1", "#20854E","#E18727","#0072B5","#EE4C97","#80796B","#BC3C29"))+
  scale_x_continuous(expand = c(0,0))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y.left = (element_text(color=c(go_list[[7]]$col,go_list[[6]]$col,go_list[[5]]$col,go_list[[4]]$col,go_list[[3]]$col,go_list[[2]]$col,go_list[[1]]$col))),
        panel.grid.major=element_line(colour = "white"),
        panel.grid.minor=element_line(colour = "white"))
dev.off()




