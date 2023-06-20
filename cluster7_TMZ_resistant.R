setwd("~/scRNA_GBM_integration/Martin_project/")
library(GSVA)
library(ggplot2)
library(ggpubr)
library(limma)
library(ggnewscale)
library(clusterProfiler)

min.max.norm <- function(x){
  ((x-min(x))/(max(x)-min(x)))
}

# Ctrl_DIP_sample_Markers <- readRDS("~/scRNA_GBM_integration/Martin_project/version2/processed_data/tumor_Markers.rds")
# Ctrl_DIP_sample_Markers<-Ctrl_DIP_sample_Markers[Ctrl_DIP_sample_Markers$cluster %in% c("Cluster7"),]
# Ctrl_DIP_sample_Markers<-Ctrl_DIP_sample_Markers[Ctrl_DIP_sample_Markers$avg_log2FC>0.5 & Ctrl_DIP_sample_Markers$p_val_adj<0.05,]

tumor_cluster_DEG<-readRDS(file = "./version2/processed_data/tumor_cluster_DEG.rds")
tumor_cluster_DEG<-tumor_cluster_DEG[tumor_cluster_DEG$cluster %in% c("Cluster7"),]
tumor_cluster_DEG_up<-tumor_cluster_DEG[tumor_cluster_DEG$avg_log2FC>0.5 & tumor_cluster_DEG$p_val_adj<0.05,]
tumor_cluster_DEG_down<-tumor_cluster_DEG[tumor_cluster_DEG$avg_log2FC< -0.5 & tumor_cluster_DEG$p_val_adj<0.05,]
# identical(Ctrl_DIP_sample_Markers$gene,Ctrl_DIP_sample_Markers$gene)

gene_list<-list(cluster7_up=tumor_cluster_DEG_up$gene,
                cluster7_down=tumor_cluster_DEG_down$gene)

EMTAB_2693_expr <- readRDS("~/bioinfo_mill/dataset/glioma/processed_data/EMTAB_2693_expr.rds")
EMTAB_2693_phe <- readRDS("~/bioinfo_mill/dataset/glioma/processed_data/EMTAB_2693_phe.rds")
EMTAB_2693_phe$patient<-stringr::str_split(EMTAB_2693_phe$clone,"_",simplify = T)[,1]

#PCA
pca <- prcomp(t(EMTAB_2693_expr), scale=F)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
identical(rownames(EMTAB_2693_phe),rownames(pca.data))
pca.data<-cbind(pca.data,EMTAB_2693_phe)

pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation

tiff('./version2/fig/cluster7_TMZtherapy_GBM_PCA.tiff',units ="in",width = 5.8,height = 3.8,res = 1200,compression ='zip')
ggplot(pca.data, aes(x=X,y=Y,color=type,shape=patient),size=5) +
  geom_point(size=3) +#stat_ellipse(level = 0.95, show.legend = F,geom = "polygon",alpha = 1/5, aes(fill = patient))+
  scale_colour_manual(values =c("#F5A200","#A70B00"))+
  xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep="")) +
  ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep="")) + 
  theme_classic()+
  theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
                      legend.key.size = unit(0.8, "cm"),legend.text = element_text(size=15),legend.title =element_text(size=15))
dev.off()

##cluster7 signature enrichment
ssGSEA<-gsva(as.matrix(EMTAB_2693_expr), gene_list,method="ssgsea")
ssGSEA<-as.data.frame(t(ssGSEA))
identical(rownames(EMTAB_2693_phe),rownames(ssGSEA))
ssGSEA<-cbind(ssGSEA,EMTAB_2693_phe)

ggstripchart(ssGSEA[ssGSEA$patient %in% c("482","472"),], "type", "cluster7",
             color = "type", palette = c("#F5A200","#A70B00"),
             add = "mean_sd")+
  theme_classic()+
  labs(title = "TMZ-therapy GBM clones")+
  ylab("Cluster7 signatures") +
  xlab("") +
  stat_compare_means(aes(group=type),method = "wilcox.test",hide.ns = F,label.x = 1.5,label.y = 4,
                     label = "p.format")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

tmp<-ssGSEA[ssGSEA$patient %in% c("482","472"),]
tmp$type<-factor(tmp$type,levels = c("TMZ sensitive","TMZ resistant"))
tmp[,1:2]<-apply(tmp[,1:2],2,min.max.norm)

tiff('./version2/fig/cluster7_up_TMZtherapy_GBM_distribution.tiff',units ="in",width = 2.7,height = 3.5,res = 1200,compression ='zip')
ggplot(tmp, aes(type, cluster7_up,fill =type ))+
  scale_fill_manual(values =c("#F5A200","#A70B00"))+
  geom_violin(alpha=0.25, size=1,color="NA",trim=FALSE)+
  geom_boxplot( outlier.size = -1, color="black",lwd=0.2, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 0.6) +
  theme_classic()+
  labs(title = "TMZ-therapy GBM clones")+
  ylab("Cluster7 upregulated signatures") +
  xlab("") +
  stat_compare_means(aes(group=type),method = "wilcox.test",hide.ns = F,label.x = 1.2,label.y = 1.2,
                     label = "p.format")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
dev.off()

tiff('./version2/fig/cluster7_down_TMZtherapy_GBM_distribution.tiff',units ="in",width = 2.7,height = 3.5,res = 1200,compression ='zip')
ggplot(tmp, aes(type, cluster7_down,fill =type ))+
  scale_fill_manual(values =c("#F5A200","#A70B00"))+
  geom_violin(alpha=0.25, size=1,color="NA",trim=FALSE)+
  geom_boxplot( outlier.size = -1, color="black",lwd=0.2, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 0.6) +
  theme_classic()+
  labs(title = "TMZ-therapy GBM clones")+
  ylab("Cluster7 downregulated signatures") +
  xlab("") +
  stat_compare_means(aes(group=type),method = "wilcox.test",hide.ns = F,label.x = 1.4,label.y = 1.2,
                     label = "p.format")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
dev.off()

## DEG between TMZ-sensitive clones and TMZ-resistant clones
EMTAB_2693_phe_select<-EMTAB_2693_phe[EMTAB_2693_phe$patient %in% c("482","472"),]
EMTAB_2693_expr_select<-EMTAB_2693_expr[,rownames(EMTAB_2693_phe_select)]

group <- EMTAB_2693_phe_select$type
group <- factor(group,levels = c("TMZ sensitive","TMZ resistant"))
design <- model.matrix(~group)
fit <- lmFit(EMTAB_2693_expr_select,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 

allDiff$significant<-ifelse(allDiff$logFC>0.5&allDiff$adj.P.Val<0.05,"TMZ resistant",
                            ifelse(allDiff$logFC< -0.5&allDiff$adj.P.Val<0.05,"TMZ sensitive","unsignificant"))
table(allDiff$significant)

TMZ_resistant_inter<-intersect(rownames(allDiff[allDiff$significant=="TMZ resistant",]),tumor_cluster_DEG_up$gene)
TMZ_sensitive_inter<-intersect(rownames(allDiff[allDiff$significant=="TMZ sensitive",]),tumor_cluster_DEG_down$gene)
length(TMZ_resistant_inter)
length(TMZ_sensitive_inter)

rownames(tumor_cluster_DEG_up)<-tumor_cluster_DEG_up$gene
tmp1<-tumor_cluster_DEG_up[TMZ_resistant_inter,]
tmp2<-allDiff[TMZ_resistant_inter,]

comb<-tmp2
comb$log2FC_scRNAseq<-tmp1$avg_log2FC
comb<-comb[,c(1,7:8)]
colnames(comb)[1:2]<-c("log2FC_clones","type")

rownames(tumor_cluster_DEG_down)<-tumor_cluster_DEG_down$gene
tmp1<-tumor_cluster_DEG_down[TMZ_sensitive_inter,]
tmp2<-allDiff[TMZ_sensitive_inter,]
tmp2$log2FC_scRNAseq<-tmp1$avg_log2FC
tmp2<-tmp2[,c(1,7:8)]
colnames(tmp2)[1:2]<-c("log2FC_clones","type")
comb<-rbind(comb,tmp2)
comb<-comb[order(comb$log2FC_clones,decreasing = T),]
comb$gene<-rownames(comb)
saveRDS(comb,file = "./version2/processed_data/cluster7_TMZtherapy_gene.rds")

# comb<-comb[comb$type=="TMZ resistant",]
tmp<-reshape2::melt(comb)
tmp$var <- rep(1:nrow(comb),2)

res<-tmp[tmp$variable=="log2FC_scRNAseq",]
rownames(res)<-1:nrow(res)
res$ang <- seq(from = (360/nrow(res)) / 1.5,
               to = (1.5* (360/nrow(res))) - 360,
               length.out = nrow(res)) + 80

res$hjust <- 0
res$hjust[which(res$ang < -90)] <- 1
res$ang[which(res$ang < -90)] <- (180+res$ang)[which(res$ang < -90)]
res$color<-ifelse(res$type=="TMZ resistant","#A70B00","#F5A200")

colnames(tmp)[4]<-"Log2FC"

# ggplot(tmp,aes(x = var,y = variable)) +
#   geom_tile(aes(fill = Log2FC),color = 'white') +
#   scale_fill_gradient2(midpoint = 0,
#                        low = '#3C8DAD',
#                        mid="white",
#                        high = '#A71B4B') +
#   scale_y_discrete(expand = expansion(mult = c(3,0))) +
#   scale_x_discrete(expand = expansion(mult = c(0,0.05))) +
#   coord_polar(theta = 'x') +
#   theme_void() +
#   geom_text(data = res,
#             aes(x = as.numeric(rownames(res)),
#                 y = 2.6,
#                 label = gene, angle = ang, hjust = hjust),
#             color = res$color,
#             size = 3.5)+
#   theme(legend.position = c(0.5,0.5))
# 
# tiff('./version2/fig/cluster7_TMZtherapy_gene.tiff',units ="in",width = 6.5,height = 6.5,res = 1200,compression ='zip')
# ggplot(tmp,aes(x = var,y = variable)) +
#   geom_tile(aes(fill = Log2FC),color = 'white') +
#   scale_fill_gradientn(colours =paletteer::paletteer_c("ggthemes::Classic Red", 30)) +
#   scale_y_discrete(expand = expansion(mult = c(3,0))) +
#   scale_x_discrete(expand = expansion(mult = c(0,0.05))) +
#   coord_polar(theta = 'x') +
#   theme_void() +
#   geom_text(data = res,
#             aes(x = as.numeric(rownames(res)),
#                 y = 2.6,
#                 label = gene, angle = ang, hjust = hjust),
#             color = res$color,
#             size = 3.5)+
#   theme(legend.position = c(0.5,0.5))
# dev.off()
# 
# tiff('./version2/fig/cluster7_TMZtherapy_gene.tiff',units ="in",width = 6.5,height = 4.8,res = 1200,compression ='zip')
# ggplot() +
#   geom_tile(data = tmp[which(tmp$variable == 'log2FC_scRNAseq'),],
#             aes(x = 1:nrow(res),y = 1,fill = Log2FC),
#             color = 'white') +
#   scale_fill_gradientn(colours =paletteer::paletteer_c("ggthemes::Classic Red", 30),name="log2FC in sc")+
#   new_scale("fill") +
#   geom_tile(data = tmp[which(tmp$variable == 'log2FC_clones'),],
#             aes(x = 1:nrow(res),y = 2,fill = Log2FC),
#             color = 'white') +
#   scale_fill_gradient2(midpoint = 0,low = '#3C8DAD', mid="white", high = '#A71B4B',name="log2FC in clones") +
#   ylim(-2,3) +
#   coord_polar(theta = 'x') +
#   theme_void()+
#   geom_text(data = res,
#             aes(x = as.numeric(rownames(res)),
#                 y = 2.6,
#                 label = gene, angle = ang, hjust = hjust),
#             color = res$color,
#             size = 3.5,) +
#   xlim(0,55)
# dev.off()

tiff('./version2/fig/cluster7_TMZtherapy_gene.tiff',units ="in",width = 9.5,height = 6.5,res = 1200,compression ='zip')
ggplot() +
  geom_tile(data = tmp[which(tmp$variable == 'log2FC_scRNAseq'),],
            aes(x = 1:nrow(res),y = 1,fill = Log2FC),
            color = 'white') +
  scale_fill_gradientn(colours =paletteer::paletteer_c("ggthemes::Red-Green-Gold Diverging", 30,direction = -1),name="log2FC in sc")+
  new_scale("fill") +
  geom_tile(data = tmp[which(tmp$variable == 'log2FC_clones'),],
            aes(x = 1:nrow(res),y = 2,fill = Log2FC),
            color = 'white') +
  scale_fill_gradient2(midpoint = 0,low = '#3C8DAD', mid="white", high = '#A71B4B',name="log2FC in clones") +
  ylim(-2,3) +
  coord_polar(theta = 'x') +
  theme_void()+
  geom_text(data = res,
            aes(x = as.numeric(rownames(res)),
                y = 2.6,
                label = gene, angle = ang, hjust = hjust),
            color = res$color,
            size = 3.5,)+
  xlim(0,138)
dev.off()

## enrichment
ids <- bitr(comb$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db",drop = TRUE)
comb_id=merge(comb,ids,by.x='gene',by.y='SYMBOL',sort=F)

enrichGO <- compareCluster(ENTREZID~type, data=comb_id, fun ="enrichGO",ont="BP",OrgDb='org.Hs.eg.db',readable =T) 
saveRDS(enrichGO,file = "./version2/processed_data/cluster7_TMZtherapy_gene_enrichGO.rds")

comb_id<-as.data.frame(enrichGO)
comb_id<-comb_id[comb_id$qvalue<0.05,]
write.csv(comb_id,file = "./version2/processed_data/cluster7_TMZtherapy_gene_enrichGO_result.csv")

select_enrichment<-comb_id[c(grep("hypoxia",comb_id$Description),grep("decreased oxygen levels",comb_id$Description),
                         grep("activin receptor signaling pathway",comb_id$Description),grep("osmotic stress",comb_id$Description),
                         grep("mitotic G1 DNA damage checkpoint",comb_id$Description),
                         grep("GO:0000082",comb_id$ID),
                         grep("nuclear DNA replication",comb_id$Description),grep("GO:0006282",comb_id$ID),
                         grep("CD40 signaling pathway",comb_id$Description),
                         grep("response to radiation",comb_id$Description)),]
select_enrichment$Description<-factor(select_enrichment$Description,levels = c(unique(select_enrichment$Description)))

tiff('./version2/fig/cluster7_TMZtherapy_gene_GOenrich.tiff',units ="in",width = 4.1,height = 3.8,res = 1200,compression ='zip')
ggplot(select_enrichment, aes(x=type,y=Description,size=Count,color=p.adjust))  + 
  geom_point() +
  theme(panel.background = element_blank(),
        panel.grid.major=element_line(colour = "grey"),
        panel.border = element_rect(colour = "black", fill=NA,size = 1))+
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2","#7e62a3"),
                        #trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  labs(title="GO enrichment", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#plot enrichment and gene logfc together 
unlist(stringr::str_split(select_enrichment[1,"geneID"],"/"))
go_res<-data.frame()
for (i in 1:nrow(select_enrichment)) {
  go_tmp<-data.frame(gene=unlist(stringr::str_split(select_enrichment[i,"geneID"],"/")),
                     GOenrichment=rep(select_enrichment[i,"Description"],length(unlist(stringr::str_split(select_enrichment[i,"geneID"],"/")))),
                     padjust=select_enrichment[i,"p.adjust"])
  go_res<-rbind(go_res,go_tmp)
}
go_res$count<- -0.5

tmp2<-tmp[which(tmp$variable == 'log2FC_scRNAseq'),]
tmp2<-merge(tmp2,go_res,by="gene",all=T,sort=F)

ggplot() +
  geom_tile(data = tmp[which(tmp$variable == 'log2FC_scRNAseq'),],
            aes(x = 1:nrow(res),y = 0.5,fill = Log2FC),
            color = 'white') +
  scale_fill_gradientn(colours =paletteer::paletteer_c("ggthemes::Red-Green-Gold Diverging", 30,direction = -1),name="log2FC in sc")+
  new_scale("fill") +
  geom_tile(data = tmp[which(tmp$variable == 'log2FC_clones'),],
            aes(x = 1:nrow(res),y = 1.5,fill = Log2FC),
            color = 'white') +
  scale_fill_gradient2(midpoint = 0,low = '#3C8DAD', mid="white", high = '#A71B4B',name="log2FC in clones")  +
  new_scale("fill")+
  geom_bar(data=tmp2,aes(x=var, y=count, fill=GOenrichment), stat="identity", alpha=1)+
  scale_fill_manual(values=paletteer::paletteer_d("ggsci::category20_d3"),name="GO enrichment")+
  ylim(-3,3) +
  coord_polar(theta = 'x') +
  theme_void()+
  geom_text(data = res,
            aes(x = as.numeric(rownames(res)),
                y = 2.1,
                label = gene, angle = ang, hjust = hjust),
            color = res$color,
            size = 3.5)
