setwd("~/scRNA_GBM_integration/Martin_project/")

library(Seurat)
library(ggplot2)
library(forcats)
library(ggpubr)
library(biomaRt)
library(dplyr)
library(dittoSeq)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AUCell)
library(data.table)
library(enrichplot)

min.max.norm <- function(x){
  ((x-min(x))/(max(x)-min(x)))
}
gmtPathways <- function(gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  pathways
}

###-------------------- tumor recluster--------------------
tumor<-readRDS(file = "./processed_data/Ctrl_DIP_tumor_sample.rds")
tumor$RNA_snn_res.1.2<-factor(tumor$RNA_snn_res.1.2,levels = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
DimPlot(tumor, reduction = "umap",group.by = "RNA_snn_res.1.2",label = T,split.by ="orig.ident")
table(tumor$cluster)
##rename cell
tumor@meta.data$cluster <- tumor@meta.data$RNA_snn_res.1.2
tumor@meta.data$cluster <- dplyr::recode(tumor@meta.data$cluster,
                                                     "0"="Cluster5","1"="Cluster1",
                                                     "2"="Cluster4","3"="Cluster1",
                                                     "4"="Cluster3","5"="Cluster2",
                                                     "6"="Cluster2","7"="Cluster6",
                                                     "8"="Cluster7","9"="Cluster3","10"="Cluster8",
                                                     "11"="Cluster4","12"="Cluster3",
                                                     "13"="Cluster9","14"="Cluster7",
                                                     "15"="Cluster9")
tumor$cluster<-factor(tumor$cluster,levels = c("Cluster1","Cluster2","Cluster3","Cluster4",
                                               "Cluster5","Cluster6","Cluster7","Cluster8","Cluster9"))

tiff('./version2/fig/scRNAseq_tumor_cell_cluster.tiff',units ="in",width = 8,height = 4.5,res = 1200,compression ='zip')
DimPlot(tumor, reduction = "umap",group.by = "cluster",label = F,split.by ="orig.ident",
        cols=paletteer::paletteer_d("ggthemes::Tableau_20"))+
  labs(title="Tumor cell")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

tiff('./version2/fig/scRNAseq_tumor_cell_cluster_prop.tiff',units ="in",width = 4,height = 2.5,res = 1200,compression ='zip')
dittoBarPlot(tumor, "orig.ident", group.by = "cluster",xlab=NULL,color.panel =c("#EF6A63","#A1C8DC"),
             retain.factor.levels = T)+labs(title="")
dev.off()

table(tumor$cluster,tumor$orig.ident)

tiff('./version2/fig/scRNAseq_tumor_cell_Phase_prop.tiff',units ="in",width = 2.8,height = 2.5,res = 1200,compression ='zip')
dittoBarPlot(tumor, "orig.ident", group.by = "Phase",xlab=NULL,color.panel =c("#EF6A63","#A1C8DC"),
             retain.factor.levels = T)+labs(title="")
dev.off()


saveRDS(tumor,file = "./version2/processed_data/Ctrl_DIP_tumor_sample.rds")

###-------------------- markers in each cluster--------------------
Idents(tumor) <- "cluster"
tumor_Markers <- FindAllMarkers(tumor, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(tumor_Markers,file = "./version2/processed_data/tumor_Markers.rds")

tumor_cluster_DEG <- FindAllMarkers(tumor, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(tumor_cluster_DEG,file = "./version2/processed_data/tumor_cluster_DEG.rds")

###--------------------specific marker genes of cluster7--------------------
Markers_7<-tumor_Markers[tumor_Markers$cluster %in% c("Cluster7"),]
Markers_7<-Markers_7[Markers_7$avg_log2FC>0.5 & Markers_7$p_val_adj<0.05,]
write.csv(Markers_7,file = "./version2/processed_data/Markers_cluster7.csv")

gene_list <- bitr(Markers_7$gene, fromType = "SYMBOL",
                  toType =  "ENTREZID",
                  OrgDb = org.Hs.eg.db)

go_cluster7 <- enrichGO(gene_list$ENTREZID, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05,keyType = 'ENTREZID',readable =T)
saveRDS(go_cluster7,file = "./version2/processed_data/go_cluster7.rds")

tmp<-go_cluster7@result
tmp<-tmp[tmp$qvalue<0.05,]
write.csv(tmp,file = "./version2/processed_data/go_cluster7_result.csv")

neuron_related_enrichment<-tmp[c(grep("axon",tmp$Description),grep("neuron",tmp$Description),grep("synapse",tmp$Description)),]
write.csv(neuron_related_enrichment,file = "./version2/processed_data/go_cluster7_neuron_related.csv")

tiff('./version2/fig/marker7_GOBP.tiff',units ="in",width = 6.5,height = 4.5,res = 1200,compression ='zip')
ggplot(go_cluster7, showCategory =go_cluster7[go_cluster7@result$ID %in% neuron_related_enrichment[-c(6:7),]$ID ,"Description"],
       aes(Count,fct_reorder(Description, Count))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2","#7e62a3"),trans = "log10",
    guide=guide_colorbar(reverse=TRUE,order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  DOSE::theme_dose(12) +
  xlab("") +ylab("")+
  labs(title = "GO BP enrichment")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###--------------------specific marker genes of cluster9--------------------
Markers_9<-tumor_Markers[tumor_Markers$cluster %in% c("Cluster9"),]
Markers_9<-Markers_9[Markers_9$avg_log2FC>0.5 & Markers_9$p_val_adj<0.05,]

gene_list <- bitr(Markers_9$gene, fromType = "SYMBOL",
                  toType =  "ENTREZID",
                  OrgDb = org.Hs.eg.db)

go_cluster9 <- enrichGO(gene_list$ENTREZID, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05,keyType = 'ENTREZID',readable =T)
saveRDS(go_cluster9,file = "./version2/processed_data/go_cluster9.rds")

tmp<-go_cluster9@result
tmp<-tmp[tmp$qvalue<0.05,]
write.csv(tmp,file = "./version2/processed_data/go_cluster9_result.csv")

neuron_related_enrichment<-tmp[c(grep("axon",tmp$Description),grep("neuron",tmp$Description),grep("synapse",tmp$Description)),]
write.csv(neuron_related_enrichment,file = "./version2/processed_data/go_cluster9_neuron_related.csv")



###--------------------gsea between cluster7 and cluster9--------------------
cluster7and9_markers <- FindMarkers(tumor, ident.1="Cluster7", ident.2 ="Cluster9", test.use ="DESeq2",
                       only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
deg_list<-cluster7and9_markers[order(cluster7and9_markers$avg_log2FC,decreasing = T),]
deg_list<-deg_list[deg_list$p_val_adj<0.05,]
deg_list2 <- bitr(rownames(deg_list), fromType = "SYMBOL",
                  toType =  "ENTREZID",
                  OrgDb = org.Hs.eg.db)
deg_list2$avg_log2FC <- deg_list$avg_log2FC[match(deg_list2$SYMBOL,rownames(deg_list))]
geneList<-deg_list2$avg_log2FC
names(geneList)<-deg_list2$ENTREZID

gseGO <- gseGO(geneList = geneList, 
               OrgDb = org.Hs.eg.db, 
               keyType = "ENTREZID", 
               ont="BP",
               pvalueCutoff =1)
sort_gseGO<-gseGO[order(gseGO$NES,decreasing=T)]
sort_gseGO<-sort_gseGO[sort_gseGO$qvalues<0.05,]

###--------------------specific marker genes of cell cycle cluster--------------------
Markers<-tumor_Markers[tumor_Markers$cluster %in% c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6"),]
Markers<-Markers[Markers$avg_log2FC>0.5 & Markers$p_val_adj<0.05,]

gene_list <- bitr(Markers$gene, fromType = "SYMBOL",
                  toType =  "ENTREZID",
                  OrgDb = org.Hs.eg.db)

go_cluster <- enrichGO(gene_list$ENTREZID, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05,keyType = 'ENTREZID',readable =T)
tmp<-go_cluster@result
tmp<-tmp[tmp$qvalue<0.05,]

###--------------------Hallmark AUC enrichment--------------------
Hallmarker <- gmtPathways("~/bioinfo_mill/transfrom_proj/data/h.all.v7.1.symbols.gmt")
names(Hallmarker)<-stringr::str_split(names(Hallmarker),"HALLMARK_",simplify = T)[,2]

cells_rankings <- AUCell_buildRankings(tumor@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(Hallmarker, cells_rankings)

res<-tumor@meta.data
res<-dplyr::select(res,c("orig.ident","cluster"))
for (i in cells_AUC@NAMES) {
  aucs <- as.numeric(getAUC(cells_AUC)[i,])
  res$AUC<-aucs
  colnames(res)[colnames(res)=="AUC"]<-i
}
saveRDS(res,file = "./version2/processed_data/Hallmark_AUC_tumor.rds")

res2<-aggregate(x = res[,3:ncol(res)],
                by = list(res$cluster),
                FUN = mean)
rownames(res2)<-res2$Group.1
res2<-res2[,-1]
res2<-as.matrix(t(res2))

sort_row <- c()
sort_column <- c()

for(i in colnames(res2)){
  select_row <- which(matrixStats::rowMaxs(res2,na.rm = T) == res2[,i])
  tmp <- rownames(res2)[select_row][order(res2[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(res2[sort_row,],2,function(x) order(x)[nrow(res2)])
sort_column <- names(sort_column)

tiff('./version2/fig/scRNAseq_tumor_cell_cluster_Hallmarkpathway.tiff',units ="in",width = 5.5,height = 8,res = 1200,compression ='zip')
ComplexHeatmap::pheatmap(res2[sort_row,sort_column],cluster_cols = F,cluster_rows = F,scale = "row",
                         fontsize_row=9,fontsize_col=9,border=FALSE,
                         color=colorRampPalette(c("#1E3163", "#00C1D4", "#FFED99","#FF7600"))(30))
dev.off()

###--------------------gene sets AUC enrichment--------------------
tumor<-readRDS(file = "./version2/processed_data/Ctrl_DIP_tumor_sample.rds")

Hallmarker <- gmtPathways("~/bioinfo_mill/transfrom_proj/data/h.all.v7.1.symbols.gmt")
GOBP <- gmtPathways("~/software/msigdb_v7.2/c5.go.bp.v7.2.symbols.gmt")

neuron_related<-GOBP[c(grep("AXON",names(GOBP)),grep("NEURON",names(GOBP)),grep("SYNAPSE",names(GOBP)))]
neuron_related<-unique(as.character(unlist(neuron_related)))
EMT_related<-GOBP[c(grep("MESENCHY",names(GOBP)))]
EMT_related<-unique(as.character(unlist(EMT_related)))
ECM_related<-GOBP[c(grep("EXTRACELLULAR_MATRIX",names(GOBP)))]
ECM_related<-unique(as.character(unlist(ECM_related)))
synape_gene<-c("PTPRS", "ARHGEF2","GRIK2", "DNM3", "LRRTM2", "GRIK5", "NLGN4X", "NRCAM","MAP2", "INA" ,"TMPRSS9")

gene_list<-list("Neuron"=neuron_related,
                "Hypoxia"=Hallmarker$HALLMARK_HYPOXIA,
                "EMT"=EMT_related,
                "ECM"=ECM_related,
                "Synape"=synape_gene)

cells_rankings <- AUCell_buildRankings(tumor@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(gene_list, cells_rankings)
saveRDS(cells_AUC,file = "./version2/processed_data/genesets_AUC_tumor.rds")

res<-tumor@meta.data
res<-dplyr::select(res,c("orig.ident","cluster"))
for (i in cells_AUC@NAMES) {
  aucs <- as.numeric(getAUC(cells_AUC)[i,])
  res$AUC<-aucs
  colnames(res)[colnames(res)=="AUC"]<-i
}

res[,3:7]<-apply(res[,3:7],2,min.max.norm)

identical(rownames(res),colnames(tumor))
tumor <- AddMetaData(object = tumor, metadata = res[,3:7], col.name = colnames(res)[3:7])

tiff('./version2/fig/cluster_genesets_enrich.tiff',units ="in",width = 8.6,height = 8,res = 1200,compression ='zip')
FeaturePlot(tumor,features = c("MKI67","CDK1","TOP2A","Neuron","EMT",
                               "ECM","Hypoxia"),reduction = "umap",ncol =3,
            cols =paletteer::paletteer_c("grDevices::Purple-Yellow", 5 ,direction = -1)) 
dev.off()

tiff('./version2/fig/cluster_genesets_enrich2.tiff',units ="in",width = 8.6,height = 2.7,res = 1200,compression ='zip')
FeaturePlot(tumor,features = c("Neuron","EMT","ECM"),reduction = "umap",ncol =3,min.cutoff = 0.35,
            cols =paletteer::paletteer_c("grDevices::Purple-Yellow", 5 ,direction = -1)) 
dev.off()

##split Neuronal projection, Axonogenesis(Axon regeneration)
neuron_related<-GOBP[c(grep("AXON",names(GOBP)),grep("NEURON",names(GOBP)),grep("SYNAPSE",names(GOBP)))]
AXON <-GOBP[c(grep("AXON",names(GOBP)))]
NEURON <-GOBP[c(grep("NEURON",names(GOBP)))]
SYNAPSE <-GOBP[c(grep("SYNAPSE",names(GOBP)))]
Neuronal_projection <-neuron_related[c(grep("PROJECTION",names(neuron_related)))]
Axonogenesis <- neuron_related[c(grep("REGENERATION",names(neuron_related)),grep("AXONOGENESIS",names(neuron_related)))]

gene_list<-list("Neuronal_projection"=unique(as.character(unlist(Neuronal_projection))),
                "Axonogenesis"=unique(as.character(unlist(Axonogenesis))),
                "NEURON"=unique(as.character(unlist(NEURON))),
                "AXON"=unique(as.character(unlist(AXON))),
                "SYNAPSE"=unique(as.character(unlist(SYNAPSE))))
cells_rankings <- AUCell_buildRankings(tumor@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(gene_list, cells_rankings)
saveRDS(cells_AUC,file = "./version2/processed_data/genesets_AUC_tumor2.rds")

res<-tumor@meta.data
res<-dplyr::select(res,c("orig.ident","cluster"))
for (i in cells_AUC@NAMES) {
  aucs <- as.numeric(getAUC(cells_AUC)[i,])
  res$AUC<-aucs
  colnames(res)[colnames(res)=="AUC"]<-i
}
res[,3:ncol(res)]<-apply(res[,3:ncol(res)],2,min.max.norm)
identical(rownames(res),colnames(tumor))
tumor <- AddMetaData(object = tumor, metadata = res[,3:ncol(res)], col.name = colnames(res)[3:ncol(res)])

FeaturePlot(tumor,features = c("Neuronal_projection","Axonogenesis"),reduction = "umap",ncol =3,min.cutoff = 0.35,
            cols =paletteer::paletteer_c("grDevices::Purple-Yellow", 5 ,direction = -1)) 

tiff('./version2/fig/cluster_genesets_enrich3.tiff',units ="in",width = 5.8,height = 5.3,res = 1200,compression ='zip')
FeaturePlot(tumor,features = c("NEURON","AXON","Neuronal_projection","Axonogenesis"),reduction = "umap",ncol =2,min.cutoff = 0.4,
            cols =paletteer::paletteer_c("grDevices::Purple-Yellow", 5 ,direction = -1)) 
dev.off()

tiff('./version2/fig/cluster_genesets_enrich4.tiff',units ="in",width = 2.9,height = 2.6,res = 1200,compression ='zip')
FeaturePlot(tumor,features = c("SYNAPSE"),reduction = "umap",ncol =1,min.cutoff = 0.55,
            cols =paletteer::paletteer_c("grDevices::Purple-Yellow", 5 ,direction = -1)) 
dev.off()








