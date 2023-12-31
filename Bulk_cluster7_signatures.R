setwd("~/scRNA_GBM_integration/Martin_project/")

library(GSVA)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(RColorBrewer)

min.max.norm <- function(x){
  ((x-min(x))/(max(x)-min(x)))
}
gmtPathways <- function(gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  pathways
}

tcga_expr<-readRDS("~/bioinfo_mill/dataset/tcga_glioma/processed/TCGA_glioma_expr_log2tpm_2016cell.rds")
tcga_meta<-readRDS("~/bioinfo_mill/dataset/tcga_glioma/processed/TCGA_glioma_meta_2016cell.rds")
tcga_meta<-as.data.frame(tcga_meta)
tcga_expr[1:5,1:5]
rownames(tcga_meta)<-tcga_meta$Case
identical(rownames(tcga_meta),colnames(tcga_expr))

Ctrl_DIP_sample_Markers <- readRDS("~/scRNA_GBM_integration/Martin_project/version2/processed_data/tumor_Markers.rds")
Ctrl_DIP_sample_Markers<-Ctrl_DIP_sample_Markers[Ctrl_DIP_sample_Markers$cluster %in% c("Cluster7"),]
Ctrl_DIP_sample_Markers<-Ctrl_DIP_sample_Markers[Ctrl_DIP_sample_Markers$avg_log2FC>0.5 & Ctrl_DIP_sample_Markers$p_val_adj<0.05,]

gene_list<-list(cluster7=Ctrl_DIP_sample_Markers$gene)
ssGSEA<-gsva(as.matrix(tcga_expr), gene_list,method="ssgsea")
ssGSEA<-as.data.frame(t(ssGSEA))
identical(rownames(tcga_meta),rownames(ssGSEA))
ssGSEA<-cbind(ssGSEA,tcga_meta)
ssGSEA$`Survival (months)`<-as.numeric(ssGSEA$`Survival (months)`)
ssGSEA$`Vital status (1=dead)`<-as.numeric(ssGSEA$`Vital status (1=dead)`)
saveRDS(ssGSEA,file = "./version2/processed_data/tcga_bulk_cluster7.rds")
rm(list=ls())

###---------------------- cluster7 siganture distribution in tcga----------------------
ssGSEA <- readRDS("~/scRNA_GBM_integration/Martin_project/version2/processed_data/tcga_bulk_cluster7.rds")
ssGSEA[ssGSEA=="G2"]<-"LGG"
ssGSEA[ssGSEA=="G3"]<-"LGG"
ssGSEA[ssGSEA=="G4"]<-"GBM"
ssGSEA$Grade<-factor(ssGSEA$Grade,levels = c("LGG","GBM"))
#my_comparisons <- list(c("G2","G3"),c("G3","G4"),c("G2","G4"))

tiff('./version2/fig/cluster7_tcga_glioma_distribution.tiff',units ="in",width = 2,height = 3.8,res = 1200,compression ='zip')
ggplot(ssGSEA ,aes(x = Grade, y = cluster7, fill = Grade))+
  scale_fill_manual(values =c("#F5A200","#A70B00"))+
  geom_violin(alpha=0.25, position = position_dodge(width = .75),size=1,color="NA") +
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.2, alpha = 0.7)+
  theme_classic()+
  labs(title = "TCGA glioma")+
  ylab("Cluster7 signatures") +
  xlab("") +
  stat_compare_means(aes(group=Grade),method = "wilcox.test",hide.ns = F,label.x = 1.5,label.y = 2.35,
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

###---------------------- cluster7 signature and clinical feature----------------------
Hallmarker <- gmtPathways("~/bioinfo_mill/transfrom_proj/data/h.all.v7.1.symbols.gmt")
names(Hallmarker)<-stringr::str_split(names(Hallmarker),"HALLMARK_",simplify = T)[,2]
Hallmarker<-Hallmarker[c("HYPOXIA","GLYCOLYSIS","MTORC1_SIGNALING","TNFA_SIGNALING_VIA_NFKB","HEDGEHOG_SIGNALING",
                         "KRAS_SIGNALING_UP","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE",
                         "WNT_BETA_CATENIN_SIGNALING",'IL2_STAT5_SIGNALING',"NOTCH_SIGNALING","P53_PATHWAY")]
ssGSEA_pathway<-gsva(as.matrix(tcga_expr), Hallmarker,method="ssgsea")
ssGSEA_pathway<-as.data.frame(t(ssGSEA_pathway))

ssGSEA <- readRDS("./version2/processed_data/tcga_bulk_cluster7.rds")
identical(rownames(ssGSEA_pathway),rownames(ssGSEA))
ssGSEA<-cbind(ssGSEA_pathway,ssGSEA)

ssGSEA<-ssGSEA[order(ssGSEA$cluster7,decreasing = F),]

heatmap_data<-ssGSEA[,c("HYPOXIA","GLYCOLYSIS","MTORC1_SIGNALING","TNFA_SIGNALING_VIA_NFKB",
                        "KRAS_SIGNALING_UP","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE",
                        'IL2_STAT5_SIGNALING',"P53_PATHWAY")]
heatmap_data<-apply(heatmap_data, 2, min.max.norm)
meta<-ssGSEA
colnames(ssGSEA)

ha1 = HeatmapAnnotation(Cluster7_siganture=meta$cluster7,
                        Age=ifelse(meta$`Age (years at diagnosis)`>60,">60","<60"),
                        Gender=meta$Gender,
                        Grade=meta$Grade,
                        Histology=meta$Histology,
                        IDH=meta$`IDH status`,
                        Chr1p19q=meta$`1p/19q codeletion`,
                        ATRX=meta$`TERT promoter status`,
                        MGMT=meta$`MGMT promoter status`,
                        Genome_subtype=meta$`IDH/codel subtype`,
                        Transcriptome_subtype =meta$`Transcriptome Subtype`,
                        Methylation_cluster=meta$`Pan-Glioma DNA Methylation Cluster`,
                        Molecular_subtype=meta$`Supervised DNA Methylation Cluster`,
                        col = list(Cluster7_siganture=circlize::colorRamp2(c(seq(1, 2.5, 0.05)[-30]), paletteer::paletteer_c("grDevices::Inferno", 30,direction = -1)),
                                   Age =structure(names = c(">60", "<60"), brewer.pal(3, "Dark2")[1:2]),
                                   Gender =structure(names = c("female", "male"), c("#98D9E4","#FFB977")),
                                   Grade=structure(names = c("G2", "G3","G4"), c("#FFF7BD","#F5A200","#A70B00")),
                                   Histology=structure(names = c("astrocytoma", "oligodendroglioma","oligoastrocytoma","glioblastoma"),brewer.pal(8, "Set1")[1:4]),
                                   IDH =structure(names = c("Mutant", "WT","NA"), c("red", "blue","white")),
                                   Chr1p19q=structure(names = c("codel", "non-codel","NA"), c("red", "blue","white")),
                                   ATRX =structure(names = c("Mutant", "WT","NA"), c("red", "blue","white")),
                                   MGMT =structure(names = c("Methylated", "Unmethylated","NA"), c("red", "blue","white")),
                                   Genome_subtype =structure(names = c("IDHmut-codel", "IDHmut-non-codel","IDHwt","NA"), c("red", "green","blue","white")),
                                   Transcriptome_subtype =structure(names = c("ME", "NE","PN","CL","NA"), c(brewer.pal(4, "Set2"),"white")),
                                   Methylation_cluster =structure(names = c("LGm1", "LGm2","LGm3","LGm4","LGm5","LGm6","NA"), c(brewer.pal(8, "Set1")[3:8],"white")),
                                   Molecular_subtype =structure(names = c("Classic-like", "Codel","G-CIMP-high","G-CIMP-low","LGm6-GBM","Mesenchymal-like","PA-like","NA"), c("#AEC7E8","#FFBB78","#98DF8A","#FF9896","#C5B0D5","#A4BEB8","#DBDB8D","white"))),
                        na_col = "white",border = F,
                        gap = unit(3, "points"),
                        height = unit(4, "cm")
)

tiff('./version2/fig/cluster7_clinic_tcga_glioma.tiff',units ="in",width = 13,height = 6,res = 1200,compression ='zip')
Heatmap(as.matrix(t(heatmap_data)),name="score",cluster_columns = F,cluster_rows =F,
        col=colorRampPalette(c("#1E3163", "#00C1D4", "#FFED99","#FF7600"))(30),
        show_column_names =F,row_names_gp = gpar(fontsize = 10),width  = unit(10, "cm"),
        height  = unit(6, "cm"),
        show_row_names =T,top_annotation =ha1)

dev.off()

###---------------------- survival by median group----------------------
splot<-list()
## lgg+gbm sample
dat_select<-ssGSEA
mySurv<-Surv(dat_select$`Survival (months)`, dat_select$`Vital status (1=dead)`)
dat_select$group<-ifelse(dat_select[,"cluster7"]>median(dat_select[,"cluster7"]),"High","Low")
group <-dat_select$group
group <- factor(group, levels = c("Low", "High"))
survival_dat <- data.frame(group = group)
fit <- survfit(mySurv ~ group)

data.survdiff <- survdiff(mySurv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- dat_select[order(dat_select[,"cluster7"]),]

splot[[1]]<-ggsurvplot(fit, data = survival_dat ,
                       #ggtheme = theme_bw(), 
                       conf.int = T, 
                       #conf.int.style = "step",
                       censor = F, 
                       palette = c("#1B9E77","#D95F02"), 
                       legend.title = paste0("Cluster7 signature"," expression"),
                       font.legend = 10,
                       font.title = 12,font.x = 11,font.y = 11,
                       legend.labs=c(paste0("Low","(",fit$n[1],")"),
                                     paste0("High","(",fit$n[2],")")),
                       pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                  paste("p = ",round(p.val,3), sep = "")),
                                    HR, CI, sep = "\n"),
                       pval.coord=c(110,0.9),
                       pval.size=4,
                       xlab="Month",
                       title="TCGA glioma")
tiff('./version2/fig/tcga_glioma_sur.tiff',units ="in",width = 5,height = 5,res = 1200,compression ='zip') 
splot[[1]]$plot+theme(plot.title = element_text(hjust = 0.5))
dev.off()

###---------------------- cluster7 siganture in CGGA----------------------
cgga_expr <- readRDS("~/bioinfo_mill/dataset//CGGA_data/processed_data/cgga_expr.rds")
cgga_clinic <- readRDS("~/bioinfo_mill/dataset//CGGA_data/processed_data/cgga_clinic.rds")

cgga_clinic_select<-cgga_clinic[cgga_clinic$batch=="mRNA_693" & cgga_clinic$PRS_type %in% c("Primary" ,"Recurrent"),]
table(cgga_clinic_select$PRS_type)
#cgga_clinic_select<-cgga_clinic[cgga_clinic$batch=="mRNA_325" & cgga_clinic$PRS_type %in% c("Primary" ,"Recurrent"),]
cgga_expr_selcet<-cgga_expr[,cgga_clinic_select$CGGA_ID]

ssGSEA<-gsva(as.matrix(cgga_expr_selcet), gene_list,method="ssgsea")
ssGSEA<-as.data.frame(t(ssGSEA))
identical(cgga_clinic_select$CGGA_ID,rownames(ssGSEA))
ssGSEA<-cbind(ssGSEA,cgga_clinic_select)
saveRDS(ssGSEA,file = "./version2/processed_data/cgga_bulk_cluster7.rds")

splot<-list()
## TMZ therapy only
ssGSEA_select<-ssGSEA[ssGSEA$`Chemo_status (TMZ treated=1;un-treated=0)`==1,]
ssGSEA_select<-na.omit(ssGSEA_select)
ssGSEA_select<-ssGSEA_select[ssGSEA_select$`Radio_status (treated=1;un-treated=0)`==0,]

# survival in glioma
dat_select<-ssGSEA_select
res.cut <- surv_cutpoint(dat_select, time = "OS", 
                         event = "Censor (alive=0; dead=1)", 
                         variables = "cluster7", 
                         minprop = 0.5) #sample number don't less than 30%

summary(res.cut)
res.cat <- surv_categorize(res.cut)
mySurv <- Surv(res.cat$OS, res.cat$`Censor (alive=0; dead=1)`)
group <- res.cat[,3] 
group <- factor(group, levels = c("low", "high"))
survival_dat <- data.frame(group = group)
fit <- survfit(mySurv ~ group)

data.survdiff <- survdiff(mySurv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- dat_select[order(dat_select[,"cluster7"]),]

splot[[1]]<-ggsurvplot(fit, data = survival_dat ,
                       #ggtheme = theme_bw(), 
                       conf.int = T, 
                       #conf.int.style = "step",
                       censor = F, 
                       palette = c("#1B9E77","#D95F02"), 
                       legend.title = paste0("Cluster7 signature"," expression"),
                       font.legend = 10,
                       font.title = 12,font.x = 11,font.y = 11,
                       legend.labs=c(paste0("Low","(",fit$n[1],")"),
                                     paste0("High","(",fit$n[2],")")),
                       pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                  paste("p = ",round(p.val,3), sep = "")),
                                    HR, CI, sep = "\n"),
                       pval.coord=c(4000,0.9),
                       pval.size=2,
                       xlab="Day",
                       title="CGGA glioma with TMZ therapy only")
tiff('./version2/fig/cgga_glioma_TMZ_only_sur.tiff',units ="in",width = 5,height = 5,res = 1200,compression ='zip') 
splot[[1]]$plot+theme(plot.title = element_text(hjust = 0.5))
dev.off()

## TMZ/radiation therapy 
# survival in glioma
ssGSEA_select<-ssGSEA[ssGSEA$`Chemo_status (TMZ treated=1;un-treated=0)`==1,]
ssGSEA_select<-na.omit(ssGSEA_select)
ssGSEA_select<-ssGSEA_select[ssGSEA_select$`Radio_status (treated=1;un-treated=0)`==1,]

dat_select<-ssGSEA_select
res.cut <- surv_cutpoint(dat_select, time = "OS", 
                         event = "Censor (alive=0; dead=1)", 
                         variables = "cluster7", 
                         minprop = 0.5) #sample number don't less than 30%

summary(res.cut)
res.cat <- surv_categorize(res.cut)
mySurv <- Surv(res.cat$OS, res.cat$`Censor (alive=0; dead=1)`)
group <- res.cat[,3] 
group <- factor(group, levels = c("low", "high"))
survival_dat <- data.frame(group = group)
fit <- survfit(mySurv ~ group)

data.survdiff <- survdiff(mySurv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- dat_select[order(dat_select[,"cluster7"]),]

splot[[2]]<-ggsurvplot(fit, data = survival_dat ,
                       #ggtheme = theme_bw(), 
                       conf.int = T, 
                       #conf.int.style = "step",
                       censor = F, 
                       palette = c("#1B9E77","#D95F02"),
                       legend.title = paste0("Cluster7 signature"," expression"),
                       font.legend = 10,
                       font.title = 12,font.x = 11,font.y = 11,
                       legend.labs=c(paste0("Low","(",fit$n[1],")"),
                                     paste0("High","(",fit$n[2],")")),
                       pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                  paste("p = ",round(p.val,3), sep = "")),
                                    HR, CI, sep = "\n"),
                       pval.coord=c(2500,0.9),
                       pval.size=4,
                       xlab="Day",
                       title="CGGA glioma with TMZ/radiotherapy")
tiff('./version2/fig/cgga_glioma_TMZ_radiotherapy_sur.tiff',units ="in",width = 5,height = 5,res = 1200,compression ='zip') 
splot[[2]]$plot+theme(plot.title = element_text(hjust = 0.5))
dev.off()

###---------------------- cluster7 siganture distribution in cgga----------------------
ssGSEA <- readRDS("~/scRNA_GBM_integration/Martin_project/version2/processed_data/cgga_bulk_cluster7.rds")
ssGSEA[ssGSEA=="WHO II"]<-"LGG"
ssGSEA[ssGSEA=="WHO III"]<-"LGG"
ssGSEA[ssGSEA=="WHO IV"]<-"GBM"
ssGSEA$Grade<-factor(ssGSEA$Grade,levels = c("LGG","GBM"))
#my_comparisons <- list(c("G2","G3"),c("G3","G4"),c("G2","G4"))

tiff('./version2/fig/cluster7_cgga_glioma_distribution.tiff',units ="in",width = 2,height = 3.8,res = 1200,compression ='zip')
ggplot(ssGSEA[-which(is.na(ssGSEA$Grade)),] ,aes(x = Grade, y = cluster7, fill = Grade))+
  scale_fill_manual(values =c("#F5A200","#A70B00"))+
  geom_violin(alpha=0.25, position = position_dodge(width = .75),size=1,color="NA") +
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.2, alpha = 0.7)+
  theme_classic()+
  labs(title = "CGGA glioma")+
  ylab("Cluster7 signatures") +
  xlab("") +
  stat_compare_means(aes(group=Grade),method = "wilcox.test",hide.ns = F,label.x = 1.5,label.y = 2.25,
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



