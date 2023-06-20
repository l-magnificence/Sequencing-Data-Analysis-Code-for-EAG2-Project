setwd("~/scRNA_GBM_integration/Martin_project/")

library(Seurat)
library(ggplot2)
library(biomaRt)
library(dplyr)

###--------------------1.split matrix to human and mouse--------------------
sample_all <- Read10X(data.dir = "/export2/liuhw/Martin_project/map2hg38_mm10_cellranger_aggr/aggr/outs/filtered_feature_bc_matrix/")
sample_all[1:5,1:5]

sample_info<-as.data.frame(colnames(sample_all))
sample_info$type<-stringr::str_split(sample_info[,1],"-",simplify = T)[,2]
sample_info$type<-ifelse(sample_info$type==1,"Ctrl","DIP")
table(sample_info$type)

gene<-as.data.frame(sample_all@Dimnames[[1]])
gene$species<-stringr::str_split(gene[,1],"_",simplify = T)[,1]
table(gene$species)

##ctrl sample
sample<-sample_all[,sample_info[sample_info$type=="Ctrl",1]]
map2human<-sample[gene[gene$species=="GRCh38",1],]
map2mouse<-sample[gene[gene$species=="mm10",1],]

map2human_stat<-as.data.frame(apply(map2human, 2, sum))
map2mouse_stat<-as.data.frame(apply(map2mouse, 2, sum))
identical(rownames(map2human_stat),row.names(map2mouse_stat))
map_stat<-cbind(map2human_stat,map2mouse_stat)
colnames(map_stat)<-c("map2human","map2mouse")
map_stat$total<-apply(map_stat,1,sum)
map_stat$graft_percentage<-map_stat$map2human / map_stat$total
map_stat$species<-ifelse(map_stat$graft_percentage>0.9,"graft",ifelse(map_stat$graft_percentage<0.1,"host","multiplets"))
table(map_stat$species)

p1<-ggplot(data=map_stat, aes(x=total, y=graft_percentage,fill =species))+
  geom_point(shape = 21, size=3, 
             color="black", alpha = 0.6)+
  theme_classic()+
  labs(title="Ctrl sample")+
  ylab("Percentage of reads from graft") +
  xlab("Number of reads") +
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 12))+
  geom_hline(yintercept=c(0.1,0.9),linetype=4)

pdf(file="./fig/split_huamn_mouse_ctrl.pdf",bg="white",width=5,height=4,pointsize=12)
p1
dev.off()

tiff('./20211026_paper_use/split_huamn_mouse_ctrl.tiff',units ="in",width = 4.5,height = 3.7,res = 1200,compression ='zip') 
p1
dev.off()

map2human<-map2human[,rownames(map_stat[map_stat$species=="graft",])]
map2mouse<-map2mouse[,rownames(map_stat[map_stat$species=="host",])]
rownames(map2human)<- stringr::str_split(rownames(map2human),"_",simplify = T)[,2]
rownames(map2mouse)<- stringr::str_split(rownames(map2mouse),"___",simplify = T)[,2]

#Transfer mouse gene symbol into human gene symbol
Mouse_gene_id <- rownames(map2mouse)
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org") #sometimes, there will be some connection problem with biomaRt server.
class(human)
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl",host = "useast.ensembl.org")
class(mouse)
human_gene_id <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                        values = Mouse_gene_id, mart = mouse,
                        attributesL = c("hgnc_symbol"), martL = human,
                        uniqueRows = T)
colnames(human_gene_id) <- c("gene_short_name","human_symbol")
#Replace rownames of the matrix in mouse gene name into human gene symbol
Human_gene_symbol = human_gene_id$human_symbol[match(rownames(map2mouse),human_gene_id$gene_short_name)]
rownames(map2mouse) <- Human_gene_symbol
map2mouse<-map2mouse[which(!is.na(rownames(map2mouse))),] #still have more mouse genes match to one human gene

##create a merged seurat
huamn_cell <- CreateSeuratObject(counts = map2human, project = "Ctrl", 
                             min.cells =0, min.features = 0)
huamn_cell$origin<-"PDX Human cells"
mouse_cell <- CreateSeuratObject(counts = map2mouse, project = "Ctrl", 
                                 min.cells =0, min.features = 0)
mouse_cell$origin<-"PDX mouse cells"

ctrl_sample<-merge(huamn_cell,mouse_cell)

##DIP sample
sample<-sample_all[,sample_info[sample_info$type=="DIP",1]]
map2human<-sample[gene[gene$species=="GRCh38",1],]
map2mouse<-sample[gene[gene$species=="mm10",1],]

map2human_stat<-as.data.frame(apply(map2human, 2, sum))
map2mouse_stat<-as.data.frame(apply(map2mouse, 2, sum))
identical(rownames(map2human_stat),row.names(map2mouse_stat))
map_stat<-cbind(map2human_stat,map2mouse_stat)
colnames(map_stat)<-c("map2human","map2mouse")
map_stat$total<-apply(map_stat,1,sum)
map_stat$graft_percentage<-map_stat$map2human / map_stat$total
map_stat$species<-ifelse(map_stat$graft_percentage>0.9,"graft",ifelse(map_stat$graft_percentage<0.1,"host","multiplets"))
table(map_stat$species)

p2<-ggplot(data=map_stat, aes(x=total, y=graft_percentage,fill =species))+
  geom_point(shape = 21, size=3, 
             color="black", alpha = 0.6)+
  theme_classic()+
  labs(title="DIP sample")+
  ylab("Percentage of reads from graft") +
  xlab("Number of reads") +
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 12))+
  geom_hline(yintercept=c(0.1,0.9),linetype=4)

pdf(file="./fig/split_huamn_mouse_dip.pdf",bg="white",width=5,height=4,pointsize=12)
p2
dev.off()

tiff('./20211026_paper_use/split_huamn_mouse_dip.tiff',units ="in",width = 4.5,height = 3.7,res = 1200,compression ='zip') 
p2
dev.off()

map2human<-map2human[,rownames(map_stat[map_stat$species=="graft",])]
map2mouse<-map2mouse[,rownames(map_stat[map_stat$species=="host",])]
rownames(map2human)<- stringr::str_split(rownames(map2human),"_",simplify = T)[,2]
rownames(map2mouse)<- stringr::str_split(rownames(map2mouse),"___",simplify = T)[,2]

#Transfer mouse gene symbol into human gene symbol
Mouse_gene_id <- rownames(map2mouse)
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org") #sometimes, there will be some connection problem with biomaRt server.
class(human)
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl",host = "useast.ensembl.org")
class(mouse)
human_gene_id <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                        values = Mouse_gene_id, mart = mouse,
                        attributesL = c("hgnc_symbol"), martL = human,
                        uniqueRows = T)
colnames(human_gene_id) <- c("gene_short_name","human_symbol")
#Replace rownames of the matrix in mouse gene name into human gene symbol
Human_gene_symbol = human_gene_id$human_symbol[match(rownames(map2mouse),human_gene_id$gene_short_name)]
rownames(map2mouse) <- Human_gene_symbol
map2mouse<-map2mouse[which(!is.na(rownames(map2mouse))),] #still have more mouse genes match to one human gene

##create a merged seurat
huamn_cell <- CreateSeuratObject(counts = map2human, project = "DIP", 
                                 min.cells =0, min.features = 0)
huamn_cell$origin<-"PDX Human cells"
mouse_cell <- CreateSeuratObject(counts = map2mouse, project = "DIP", 
                                 min.cells =0, min.features = 0)
mouse_cell$origin<-"PDX mouse cells"

DIP_sample<-merge(huamn_cell,mouse_cell)

###--------------------2.merge ctrl and dip sample into one object--------------------
Ctrl_DIP_sample<-merge(ctrl_sample,DIP_sample)
Ctrl_DIP_sample <- PercentageFeatureSet(Ctrl_DIP_sample, pattern = "^MT-", col.name = "percent.mt")

pdf(file="./fig/QC.pdf",bg="white",width=9,height=4,pointsize=12)
VlnPlot(Ctrl_DIP_sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

Ctrl_DIP_sample <- subset(Ctrl_DIP_sample, subset = nFeature_RNA > 500 & nFeature_RNA <8000 & percent.mt < 10 )
Ctrl_DIP_sample <- NormalizeData(Ctrl_DIP_sample, normalization.method = "LogNormalize", scale.factor = 100000, verbose=T)
Ctrl_DIP_sample <- FindVariableFeatures(Ctrl_DIP_sample, selection.method = "vst", nfeatures = 3000)
Ctrl_DIP_sample <- ScaleData(Ctrl_DIP_sample,features = rownames(Ctrl_DIP_sample))

Ctrl_DIP_sample<-RunPCA(Ctrl_DIP_sample, npcs = 50, verbose = FALSE)
ElbowPlot(Ctrl_DIP_sample, ndims = 50, reduction = "pca")

Ctrl_DIP_sample <- RunUMAP(Ctrl_DIP_sample, reduction = "pca",dims = 1:25,
                       n.neighbors = 30L, n.components = 2L,
                       min.dist = 0.1, spread = 1)
#Ctrl_DIP_sample <- RunTSNE(Ctrl_DIP_sample, reduction = "pca", dims = 1:30)
Ctrl_DIP_sample <- FindNeighbors(Ctrl_DIP_sample, dims = 1:25)
Ctrl_DIP_sample <- FindClusters(Ctrl_DIP_sample, resolution = c(0.5,0.8,1,1.2))

DimPlot(Ctrl_DIP_sample, reduction = "umap",group.by = "RNA_snn_res.1.2",label = T)
DimPlot(Ctrl_DIP_sample, reduction = "umap",group.by = "origin")
DimPlot(Ctrl_DIP_sample, reduction = "umap",group.by = "orig.ident")
DimPlot(Ctrl_DIP_sample, reduction = "umap",group.by = "origin",split.by = "orig.ident")

Idents(Ctrl_DIP_sample) <- "RNA_snn_res.1.2"
FeaturePlot(Ctrl_DIP_sample,features = c("MKI67","FUT4","PROM1","CD44","NES"),reduction = "umap")
FeaturePlot(Ctrl_DIP_sample,features = c("KCNH5","PIEZO1","PIEZO2","CLIC1"),reduction = "umap")
FeaturePlot(Ctrl_DIP_sample,features = c("KCNAB2","KCNAB1","KCNH5"),reduction = "umap")

###--------------------3.define cell type--------------------
Ctrl_DIP_sample_Markers <- FindAllMarkers(Ctrl_DIP_sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(Ctrl_DIP_sample_Markers,file = "./processed_data/Ctrl_DIP_sample_Markers.rds")

top15 <- Ctrl_DIP_sample_Markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)

Gene_list<-c('MOG',"MAG","OLIG2", #Oligodendrocyte
          "GFAP","AQP4","CLU", #Astrocyte
          'PTPRC', 'P2RY12',"TMEM119","CTSS","CX3CR1",#Microglia
          "COL3A1","COL1A2","DCN", #Fibroblast
          "CD34","PECAM1","VWF", #Endothelial
          #"SOSTDC1","TRPM3","PCSK2", #neuron
          #"MKI67","SOX2",
          "PDGFRA","EGFR" ,"DLL3") #TUMOR

DotPlot(Ctrl_DIP_sample, features =Gene_list) + RotatedAxis() + #coord_flip()+  
  scale_colour_gradientn(
    limits = c(-1,3), # here set -0.2 as a cutoff for visulization
    colours=c("white", "darkmagenta", "darkorange1"),na.value="#FFFFFF1a") +  ## values less than 0.001 will be transparent white (to reduce noise)
  theme_linedraw() + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab("") + ylab("Cluster identity")

Ctrl_DIP_sample@meta.data$cell_type <- plyr::mapvalues(Ctrl_DIP_sample@meta.data$RNA_snn_res.1.2,
                                                                from=0:23,
                                                                to=c("Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                                                                     "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                                                                     "Tumor","Tumor","Tumor","Tumor",
                                                                     "Microglia",
                                                                     "Low quality",
                                                                     "Microglia",
                                                                     "Endothelial",
                                                                     "Astrocyte",
                                                                     "Fibroblast",
                                                                     "Unknown",
                                                                     "Oligodendrocyte"))

p1<-DotPlot(Ctrl_DIP_sample, features =Gene_list,group.by = "cell_type") + RotatedAxis() + #coord_flip()+  
  scale_colour_gradientn(
    limits = c(-0.2,3), # here set -0.2 as a cutoff for visulization
    colours=c("white", "darkmagenta", "darkorange1"),na.value="#FFFFFF1a") +  ## values less than 0.001 will be transparent white (to reduce noise)
  theme_linedraw() + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab("") + ylab("Cell type")

pdf(file="./fig/cellmarker_dotplot.pdf",bg="white",width=7.5,height=4,pointsize=12)
p1
dev.off()

tiff('./20211026_paper_use/cellmarker_dotplot.tiff',units ="in",width = 7.5,height = 4,res = 1200,compression ='zip')
p1
dev.off()


color_cell_type <- c('#F3B1A0', '#57C3F3', '#F1BB72', '#53A85F', '#D6E7A3', '#E95C59', '#476D87','#E5D2DD')
names(color_cell_type) <- unique(Ctrl_DIP_sample$cell_type)

color_condition <- c("#C7E9B4","#41B6C4")
names(color_condition) <- unique(Ctrl_DIP_sample$origin)

pdf(file="./fig/umap.pdf",bg="white",width=7,height=5,pointsize=12)
DimPlot(Ctrl_DIP_sample, reduction = "umap",group.by = "cell_type",label = F,repel = F,cols=color_cell_type)
dev.off()

pdf(file="./fig/expr_umaps.pdf",bg="white",width=7,height=5,pointsize=12)
DimPlot(Ctrl_DIP_sample, reduction = "umap",group.by = "orig.ident",label = F,repel = F)
dev.off()

pdf(file="./fig/genome_umaps.pdf",bg="white",width=7,height=5,pointsize=12)
DimPlot(Ctrl_DIP_sample, reduction = "umap",group.by = "origin",label = F,repel = F)
dev.off()

pdf(file="./fig/nFeature_RNA.pdf",bg="white",width=7,height=5,pointsize=12)
FeaturePlot(Ctrl_DIP_sample,features = "nFeature_RNA",reduction = "umap")
dev.off()

##
tiff('./20211026_paper_use/scRNAseq_celltype.tiff',units ="in",width = 6,height = 4.5,res = 1200,compression ='zip') 
DimPlot(Ctrl_DIP_sample, reduction = "umap",group.by = "cell_type",label = F,repel = F,cols=color_cell_type)+
  labs(title="")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

tiff('./20211026_paper_use/scRNAseq_origin.tiff',units ="in",width = 6,height = 4.5,res = 1200,compression ='zip') 
DimPlot(Ctrl_DIP_sample, reduction = "umap",group.by = "origin",label = F,repel = F,cols=color_condition)+
  labs(title="")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


dittoBarPlot(Ctrl_DIP_sample, "orig.ident", group.by = "cell_type",xlab=NULL,color.panel =color_condition)


saveRDS(Ctrl_DIP_sample,file = "./processed_data/Ctrl_DIP_sample_aggr_merge.rds")




