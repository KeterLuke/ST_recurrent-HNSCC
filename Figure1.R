####load pcakges####
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(data.table)
library(hdf5r)
library(STutility)
library(future)
library(SCpubr)
library(spdep)
library(harmony)
library(clustree)
library(purrr)
library(ComplexHeatmap) 
library(ggpubr)
library(magrittr)
plan('multisession',workers = 1)
options(future.globals.maxSize = 128000 * 1024^2)
load('./T.I.N.rData')
p1$labels <- factor(p1$labels,levels = c('N','I','T'))
r1$labels <- factor(r1$labels,levels = c('N','I','T'))
p2$labels <- factor(p2$labels,levels = c('N','I','T'))
r2$labels <- factor(r2$labels,levels = c('N','I','T'))

set.seed(123)
dir.create('./Figures')
dir.create('./Figures/Figure1')
cols1 <-c('N'="#e7a801",'I'="#13274f",'T'="#ce1141")
cols2 <- c('Malig'="#ce1141",'Non-malig'="#e7a801")
cols3 <- c('cluster1'="#13274f",'cluster2'="#ce1141",'cluster3'="#e7a801")
cols4 <- c('primary I' = "#fff89a",'primary T' = "#ffc920",'recurrent I' = "#086e7d",'recurrent T' = "#1a5f7a")
####Figure 1B####
pdf('./Figures/Figure1/T.N.stdimplot.r1.pdf',width = 6,height = 6)
FeatureOverlay(r1,features = 'malig_spot',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols2) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()
pdf('./Figures/Figure1/T.N.stdimplot.p1.p2.r2.pdf',width = 6,height = 6)
FeatureOverlay(p1,features = 'malig_spot',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols2) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
FeatureOverlay(p2,features = 'malig_spot',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols2) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
FeatureOverlay(r2,features = 'malig_spot',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols2) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))


dev.off()
####Figure 1C####
pdf('./Figures/Figure1/cancer.deconv.r1.pdf',width = 6,height = 6)
FeatureOverlay(r1,features = 'Cancer',add.alpha = T,type = 'raw',pt.size = 2.5) + 
  ggtitle(label = '',subtitle = '') + 
  scale_fill_viridis(option = 'H') + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'Malig score')
dev.off()
pdf('./Figures/Figure1/cancer.deconv.p1.p2.r2.pdf',width = 6,height = 6)
FeatureOverlay(p1,features = 'Cancer',add.alpha = T,type = 'raw',pt.size = 2.5) + 
  ggtitle(label = '',subtitle = '') + 
  scale_fill_viridis(option = 'H') + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'Malig score')
FeatureOverlay(p2,features = 'Cancer',add.alpha = T,type = 'raw',pt.size = 2.5) + 
  ggtitle(label = '',subtitle = '') + 
  scale_fill_viridis(option = 'H') + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'Malig score')
FeatureOverlay(r2,features = 'Cancer',add.alpha = T,type = 'raw',pt.size = 2.5) + 
  ggtitle(label = '',subtitle = '') + 
  scale_fill_viridis(option = 'H') + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'Malig score')
dev.off()
####Figure 1D####
r1$Node <- ifelse(r1$Node=='Clone1','cluster1',
                  ifelse(r1$Node=='Clone2','cluster2','cluster3'))
p1$Node <- ifelse(p1$Node=='Clone1','cluster1',
                  ifelse(p1$Node=='Clone2','cluster2','cluster3'))
p2$Node <- ifelse(p2$Node=='Clone1','cluster1',
                  ifelse(p2$Node=='Clone2','cluster2','cluster3'))
r2$Node <- ifelse(r2$Node=='Clone1','cluster1',
                  ifelse(r2$Node=='Clone2','cluster2','cluster3'))
table(r1$Node)
pdf('./Figures/Figure1/infercnv.cluster.r1.pdf',6,6)
FeatureOverlay(r1,features = 'Node',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols3) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_text(size = 15,face = 'bold',hjust = .9),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'top') + 
  labs(fill = 'CNV cluster') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()
pdf('./Figures/Figure1/infercnv.cluster.p1,p2,r2.pdf',6,6)
FeatureOverlay(p1,features = 'Node',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols3) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_text(size = 15,face = 'bold',hjust = .9),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'top') + 
  labs(fill = 'CNV cluster') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
FeatureOverlay(p2,features = 'Node',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols3) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_text(size = 15,face = 'bold',hjust = .9),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'top') + 
  labs(fill = 'CNV cluster') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
FeatureOverlay(r2,features = 'Node',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols3) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_text(size = 15,face = 'bold',hjust = .9),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'top') + 
  labs(fill = 'CNV cluster') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()
####Figure 1E####
r1_cnv_table <- read.table('./InferCNVrun_outputs/r1/infercnv.observations.txt',header = T)
colnames(r1_cnv_table) <- gsub("\\.", "-", colnames(r1_cnv_table))
cnvScore <- function(data){
  data <- data %>% as.matrix()
  cnv_score <- as.data.frame(colSums(data))
  return(cnv_score)
}
r1_cnv_score <- cnvScore(r1_cnv_table) %>% as.matrix()
r1_cnv_score <- scales::rescale(r1_cnv_score,to=c(0,1)) %>% as.data.frame()
colnames(r1_cnv_score) <- 'CNV_score'
r1 <- AddMetaData(r1,metadata = r1_cnv_score)
summary(r1$CNV_score)
pdf('./Figures/Figure1/cnv.score.vlnplot.r1.pdf')
SCpubr::do_ViolinPlot(r1,features = 'CNV_score',group.by = 'Node',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = r1$Node)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .3,varwidth = T) + 
  stat_compare_means(comparisons = list(c('cluster1','cluster2'),
                                        c('cluster1','cluster3')),label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = .05) + 
  labs(x = '') + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols3) + 
  scale_color_manual(values = cols3)
dev.off()

p1_cnv_table <- read.table('./InferCNVrun_outputs/p1/infercnv.observations.txt',header = T)
colnames(p1_cnv_table) <- gsub("\\.", "-", colnames(p1_cnv_table))
p1_cnv_score <- cnvScore(p1_cnv_table) %>% as.matrix()
p1_cnv_score <- scales::rescale(p1_cnv_score,to=c(0,1)) %>% as.data.frame()
colnames(p1_cnv_score) <- 'CNV_score'
p1 <- AddMetaData(p1,metadata = p1_cnv_score)
p2_cnv_table <- read.table('./InferCNVrun_outputs/p2/infercnv.observations.txt',header = T)
colnames(p2_cnv_table) <- gsub("\\.", "-", colnames(p2_cnv_table))
p2_cnv_table <- cnvScore(p2_cnv_table) %>% as.matrix()
p2_cnv_table <- scales::rescale(p2_cnv_table,to=c(0,1)) %>% as.data.frame()
colnames(p2_cnv_table) <- 'CNV_score'
p2 <- AddMetaData(p2,metadata = p2_cnv_table)
r2_cnv_table <- read.table('./InferCNVrun_outputs/r2/infercnv.observations.txt',header = T)
colnames(r2_cnv_table) <- gsub("\\.", "-", colnames(r2_cnv_table))
r2_cnv_score <- cnvScore(r2_cnv_table) %>% as.matrix()
r2_cnv_score <- scales::rescale(r2_cnv_score,to=c(0,1)) %>% as.data.frame()
colnames(r2_cnv_score) <- 'CNV_score'
r2 <- AddMetaData(r2,metadata = r2_cnv_score)
pdf('./Figures/Figure1/cnv.score.vlnplot.p1.p2.r2.pdf')
SCpubr::do_ViolinPlot(p1,features = 'CNV_score',group.by = 'Node',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = p1$Node)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .3,varwidth = T) + 
  stat_compare_means(comparisons = list(c('cluster1','cluster3'),
                                        c('cluster2','cluster3')),label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = .05) + 
  labs(x = '',title = 'p1') + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols3) + 
  scale_color_manual(values = cols3)
SCpubr::do_ViolinPlot(p2,features = 'CNV_score',group.by = 'Node',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = p2$Node)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .3,varwidth = T) + 
  stat_compare_means(comparisons = list(c('cluster1','cluster3'),
                                        c('cluster2','cluster3')),label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = .05) + 
  labs(x = '',title = 'p2') + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols3) + 
  scale_color_manual(values = cols3)
SCpubr::do_ViolinPlot(r2,features = 'CNV_score',group.by = 'Node',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = r2$Node)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .3,varwidth = T) + 
  stat_compare_means(comparisons = list(c('cluster1','cluster3'),
                                        c('cluster2','cluster3')),label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = .05) + 
  labs(x = '',title = 'r2') + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols3) + 
  scale_color_manual(values = cols3)
dev.off()
####Figure 1F####
pdf('./Figures/Figure1/T.I.N.stdimplot.r1.pdf',width = 6,height = 6)
FeatureOverlay(r1,features = 'labels',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols1) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()
pdf('./Figures/Figure1/T.I.N.stdimplot.p1.p2.r2.pdf',width = 6,height = 6)
FeatureOverlay(p1,features = 'labels',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols1) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
FeatureOverlay(p2,features = 'labels',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols1) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
FeatureOverlay(r2,features = 'labels',pt.alpha = 1,type = 'raw',pt.size = 2.5,cols = cols1) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))

dev.off()
####Figure 1G####
s2 <- readRDS('./Boundary/GSE208253/s2.addlocation.rds')
s2$labels <- ifelse(s2$Location=='Bdy','I',
                    ifelse(s2$Location=='Mal','T','N'))
table(s2$labels)
s2$labels <- factor(s2$labels,levels = c('N','I','T'))
pdf('./Figures/Figure1/T.I.N.GSE208253.s2.pdf',width = 6,height = 6)
Idents(s2) <- s2$labels
SpatialDimPlot(s2,pt.size = 1.8,cols = cols1) + 
  ggtitle(label = '',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()
s1 <- readRDS('./Boundary/GSE208253/s1.addlocation.rds')
s1$labels <- ifelse(s1$Location=='Bdy','I',
                    ifelse(s1$Location=='Mal','T','N'))
table(s1$labels)
s1$labels <- factor(s1$labels,levels = c('N','I','T'))
Idents(s1) <- s1$labels
plot1 <- SpatialDimPlot(s1,pt.size = 1.8,cols = cols1,image.alpha = 0) + 
  ggtitle(label = 'GSE208253-S1',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'top') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
s5 <- readRDS('./Boundary/GSE208253/s5.addlocation.rds')
s5$labels <- ifelse(s5$Location=='Bdy','I',
                    ifelse(s5$Location=='Mal','T','N'))
table(s5$labels)
s5$labels <- factor(s5$labels,levels = c('N','I','T'))
Idents(s5) <- s5$labels
plot2 <- SpatialDimPlot(s5,pt.size = 1.8,cols = cols1,image.alpha = 0) + 
  ggtitle(label = 'GSE208253-S5',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'top') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
s11 <- readRDS('./Boundary/GSE208253/s11.addlocation.rds')
s11$labels <- ifelse(s11$Location=='Bdy','I',
                     ifelse(s11$Location=='Mal','T','N'))
table(s11$labels)
s11$labels <- factor(s11$labels,levels = c('N','I','T'))
Idents(s11) <- s11$labels
plot3 <- SpatialDimPlot(s11,pt.size = 1.8,cols = cols1,image.alpha = 0) + 
  ggtitle(label = 'GSE208253-S11',subtitle = '') + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15,face = 'bold'),
        legend.position = 'top') + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
pdf('./Figures/Figure1/T.I.N.GSE208253.s1.s5.s11.pdf',width = 12,height = 5)
plot1 + plot2 + plot3 + plot_layout(guides = 'collect',ncol = 3) & 
  theme(legend.position = 'bottom')
dev.off()
####Figure 1H####
gseares <- readRDS('./IvsT.GSEAhallmark.rds')
gseadf <- gseares@result
library(GseaVis)
id <- c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
        'HALLMARK_E2F_TARGETS',
        'HALLMARK_MYC_TARGETS_V1',
        'HALLMARK_TGF_BETA_SIGNALING')
pdf('./Figures/Figure1/IvsT.gsea.hallmark.pdf',width = 16,height = 4)
plot_list <- lapply(id,function(i){
  plot <- gseaNb(object = gseares,
                 geneSetID = i,
                 subPlot = 2,
                 addPval = T,
                 rmHt = F,
                 rmSegment = T,
                 pvalX = .9,
                 pvalY = .7,
                 rank.gene.nudgey = .5,termWidth = 20,arrowAngle = 20,lineSize = 2,htHeight = .5,pDigit = 4)
  return(plot)
})
plot_list[[1]][[1]]<- plot_list[[1]][[1]] + ggtitle(label = 'EMT Transition')
plot_grid(plotlist = plot_list,ncol = 4)
dev.off()
####Figure 1I####
load('./add.modulescore.rData')
sce$labels <- factor(sce$labels,levels = c('N','I','T'))
sce.sub$labels <- factor(sce.sub$labels,levels = c('primary I','primary T','recurrent I','recurrent T'))
signature <- list(
  Hypoxia=c('KRTDAP','ERO1A','HK2','HIF1A','HILPDA','VEGFA'),
  EMT=c('CTNNB1','TGFBR1','CDH1','CDH2','TGFB1','TGFB2','TGFB3','MMP10'),
  Exhaustion=c('PDCD1','CTLA4','TIGIT','HAVCR2','LAG3','CXCL13','ENTPD1',
               'TOX','DSC2','PLXND1','MYO10','HTRA3'))
plot1 <- do_EnrichmentHeatmap(sample = sce,
                              input_gene_list = signature,ctrl = 15,nbin = 15,legend.width = 1,legend.length = 10,
                              enforce_symmetry = T,use_viridis = F,scale_scores = T,legend.title = 'Score',
                              group.by = 'labels',legend.position = 'right') + scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai1',type = 'continuous',direction = -1)) 
df <- plot1[[1]]$data
df$mean <- round(df$mean,3)
plot1 <- plot1 + geom_text(data = df,aes(label = mean),color='white',size = 4)
plot1

comparisons <- list(c('N','I'),
                    c('N','T'),
                    c('I','T'))
vlnplot1 <- SCpubr::do_ViolinPlot(sce,features = 'Hypoxia1',group.by = 'labels',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = sce$labels)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .3,varwidth = T) + 
  stat_compare_means(comparisons = comparisons,label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = 0.05) + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols1) + 
  labs(x = '') + 
  scale_color_manual(values = cols1)
vlnplot2 <- SCpubr::do_ViolinPlot(sce,features = 'EMT2',group.by = 'labels',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = sce$labels)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .3,varwidth = T) + 
  stat_compare_means(comparisons = comparisons,label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = 0.05) + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols1) + 
  labs(x = '') + 
  scale_color_manual(values = cols1)
vlnplot3 <- SCpubr::do_ViolinPlot(sce,features = 'Exhaustion3',group.by = 'labels',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = sce$labels)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .3,varwidth = T) + 
  stat_compare_means(comparisons = comparisons,label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = 0.05) + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols1) + 
  scale_color_manual(values = cols1) + 
  labs(x = '')
pdf('./Figures/Figure1/T.I.N.signature.pdf',width = 16,height = 4)
plot_grid(plot1,vlnplot1, vlnplot2,vlnplot3,axis = 'b',ncol = 4,align = 'h',rel_widths = c(4,4,4,4),rel_heights = c(4,4,4,4))
dev.off()
####Figure 1J####
table(sce.sub$region)
signature <- list(
  Hypoxia=c('KRTDAP','ERO1A','HK2','HIF1A','HILPDA','VEGFA'),
  EMT=c('CTNNB1','TGFBR1','CDH1','CDH2','TGFB1','TGFB2','TGFB3','MMP10'),
  Exhaustion=c('PDCD1','CTLA4','TIGIT','HAVCR2','LAG3','CXCL13','ENTPD1',
               'TOX','DSC2','PLXND1','MYO10','HTRA3'))

plot2 <- SCpubr::do_EnrichmentHeatmap(sample = sce.sub,
                                      input_gene_list = signature,ctrl = 15,nbin = 25,legend.width = 1,legend.length = 10,
                                      enforce_symmetry = T,use_viridis = F,scale_scores = T,legend.title = 'Score',
                                      group.by = 'region',legend.position = 'right') + scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai1',type = 'continuous',direction = -1))
df <- plot2[[1]]$data
df$mean <- round(df$mean,3)
plot2 <- plot2 + geom_text(data = df,aes(label = mean),color='white',size = 4)
print(plot2)
sce.sub <- AddModuleScore(sce.sub,features = signature,name = names(signature),ctrl = 15,nbin = 25)
sce.sub$Hypoxia1 <- vegan::decostand(sce.sub$Hypoxia1,MARGIN = 2,method = 'range')
sce.sub$EMT2 <- vegan::decostand(sce.sub$EMT2,MARGIN = 2,method = 'range')
sce.sub$Exhaustion3 <- vegan::decostand(sce.sub$Exhaustion3,MARGIN = 2,method = 'range')
comparisons <- list(c('primary I','primary T'),
                    c('recurrent I','recurrent T'),
                    c('primary I','recurrent I'),
                    c('primary T','recurrent T'))
vlnplot4 <- SCpubr::do_ViolinPlot(sce.sub,features = 'Hypoxia1',group.by = 'region',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = sce.sub$region)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .2,varwidth = F) + 
  stat_compare_means(comparisons = comparisons,label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = 0.05) + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols4) + 
  scale_color_manual(values = cols4) + 
  labs(x = '')
vlnplot5 <- SCpubr::do_ViolinPlot(sce.sub,features = 'EMT2',group.by = 'region',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = sce.sub$region)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .2,varwidth = F) + 
  stat_compare_means(comparisons = comparisons,label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = 0.05) + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols4) + 
  scale_color_manual(values = cols4) + 
  labs(x = '')
vlnplot6 <- SCpubr::do_ViolinPlot(sce.sub,features = 'Exhaustion3',group.by = 'region',pt.size = 0,axis.text.face = 'bold') + 
  geom_violin(aes(color = sce.sub$region)) + 
  geom_boxplot(color = 'black',fill = 'white',outliers = F,width = .2,varwidth = F) + 
  stat_compare_means(comparisons = comparisons,label = 'p.signif',size = 8,tip.length = 0,vjust = .5,label.y.npc = 0.05) + 
  theme(axis.text = element_text(face = 'bold',size = 15)) + 
  scale_fill_manual(values = cols4) + 
  scale_color_manual(values = cols4) + 
  labs(x = '')
pdf('./Figures/Figure1/pr.re.I.N.signature.pdf',width = 16,height = 4)
plot_grid(plot2,vlnplot4,vlnplot5,vlnplot6,ncol = 4,align = 'h',axis = 'a',rel_widths = c(4,4,4,4),rel_heights = c(4,4,4,4))
dev.off()
####Figure 1K####
load('./IvsT.deg.rData')
load('./I_deg.rData')
load('./T_deg.rData')
I.sig.deg <- I_deg %>% dplyr::filter(avg_log2FC>0.1 & p_val_adj < 0.05)
idents <- (I.sig.deg$pct.1 - I.sig.deg$pct.2)>0
I.sig.deg <- I.sig.deg[idents,]
T.sig.deg <- T_deg %>% dplyr::filter(avg_log2FC>0.1 & p_val_adj < 0.05)
idents <- (T.sig.deg$pct.1 - T.sig.deg$pct.2)>0
T.sig.deg <- T.sig.deg[idents,]
expr <- data.table::fread('./8.Bulk/bulk_raw/TCGA-HNSC.htseq_counts.tsv.gz',data.table = F)
probe <- data.table::fread('./8.Bulk/bulk_raw/gencode.v22.annotation.gene.probeMap.txt',data.table = F)
expr$Symbol <- probe$gene[match(expr$Ensembl_ID,probe$id)]
sum(is.na(expr$Symbol))
expr <- expr[!is.na(expr$Symbol),]
expr <- expr[!duplicated(expr$Symbol),]
rownames(expr) <- NULL
expr <- expr %>% tibble::column_to_rownames(var = 'Symbol')
expr <- expr[,-1]
group=sapply(strsplit(colnames(expr),"\\-"),"[",4)
table(group)
group <- ifelse(grepl(group,pattern = '11'),'Normal','Tumor')
table(group)
expr <- expr[,group=='Tumor']
dim(expr)
max(expr)
###log2(count+1)
expr <- 2^expr - 1
max(expr)
expr <- as.matrix(expr)
surv <- data.table::fread('./8.Bulk/bulk_raw/TCGA-HNSC.survival.tsv',data.table = F,header = T,)
library(stringr)
k = colnames(expr)%in%surv$sample 
table(k)
expr <- expr[,k]
k = surv$sample%in%colnames(expr)
table(k)
surv = surv[k,]
k <- match(surv$sample,colnames(expr))
sum(is.na(k))
k2 = !(is.na(surv$OS.time)|is.na(surv$OS))
table(k2)
surv = surv[k2,]
expr = expr[,sort(colnames(expr))]
colnames(expr) <- str_sub(colnames(expr),1,12)
colMeans <- function(x){
  exp_m <- as.matrix(x)
  exp_t <- t(exp_m)
  expr_t <- limma::avereps(exp_t)
  return(t(expr_t))
}
expr <- colMeans(expr)
surv <- limma::avereps(surv,surv$'_PATIENT')
surv <- as.data.frame(surv)
expr <- expr[,colnames(expr)%in%surv$'_PATIENT']
surv <- surv[surv$'_PATIENT'%in%colnames(expr),]
dim(surv)
dim(expr)
table(surv$OS)
surv$OS.time <- as.numeric(surv$OS.time)
surv$OS.time = surv$OS.time/30
range(surv$OS.time)
dim(surv)
rownames(surv) <- NULL
surv <- tibble::column_to_rownames(surv,var = 'sample')
expr <- expr[,colnames(expr)%in%surv$'_PATIENT']
dim(expr)
colnames(surv)[1:3]=c('event','patient','time')

range(surv$time)
expr <- as.matrix(expr)
library(survival)
library(survminer)
surv$event <- as.numeric(surv$event)

ct <- log2(edgeR::cpm(expr)+1)
max(ct)
library(GSVA)

signature <- list(Interface=I.sig.deg
)
ssGSEA_matrix <- gsva(expr = ct,
                      gset.idx.list = signature,
                      method = 'ssgsea',
                      kcdf = "Gaussian",
                      abs.ranking=T
)
library(vegan)
ssGSEA.norm <- decostand(ssGSEA_matrix,'range',1) %>% t() %>% as.data.frame()
ssGSEA.norm$patient <- rownames(ssGSEA.norm)


surv <- full_join(surv,ssGSEA.norm,by='patient')
surv$group = ifelse(surv$Interface>median(surv$Interface),'high','low')
table(surv$group)
sfit=survfit(Surv(time, event)~group, data=surv)
data <- summary(coxph(Surv(time, event)~group, data=surv))
data$coefficients[5]
pdf('./Figures/Figure1/I_ssgsea.survplot.pdf')
ggsurvplot(sfit,pval =TRUE, data = surv,palette = 'jco',surv.median.line = 'hv',pval.method = T,conf.int = T,legend.labs = c('High','Low'),legend.title='I_score')
dev.off()
####Figure S1A,S1B####
plot1 <- FeatureOverlay(p1,features = 'nCount_RNA',add.alpha = F,type = 'raw',pt.size = 2.5,max.cutoff = 40000) + 
  ggtitle(label = 'p1',subtitle = '') + 
  scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai3',type = 'continuous',direction = -1)) + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'nCount')
plot2 <- FeatureOverlay(r1,features = 'nCount_RNA',add.alpha = F,type = 'raw',pt.size = 2.5,max.cutoff = 40000) + 
  ggtitle(label = 'r1',subtitle = '') + 
  scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai3',type = 'continuous',direction = -1)) + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'nCount')
plot3 <- FeatureOverlay(p2,features = 'nCount_RNA',add.alpha = F,type = 'raw',pt.size = 2.5,max.cutoff = 40000) + 
  ggtitle(label = 'p2',subtitle = '') + 
  scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai3',type = 'continuous',direction = -1)) + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'nCount')
plot4 <- FeatureOverlay(r2,features = 'nCount_RNA',add.alpha = F,type = 'raw',pt.size = 2.5,max.cutoff = 40000) + 
  ggtitle(label = 'r2',subtitle = '') + 
  scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai3',type = 'continuous',direction = -1)) + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'nCount')
pdf('./Figures/Figure1/nCount.stfeatureplot.pdf',width = 16,height = 5)
plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 4,guides = 'collect')  &
  theme(legend.position='bottom')
dev.off()

plot1 <- FeatureOverlay(p1,features = 'nFeature_RNA',add.alpha = F,type = 'raw',pt.size = 2.5,max.cutoff = 8000) + 
  ggtitle(label = 'p1',subtitle = '') + 
  scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai1',type = 'continuous',direction = -1)) + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'nFeature')
plot2 <- FeatureOverlay(r1,features = 'nFeature_RNA',add.alpha = F,type = 'raw',pt.size = 2.5,max.cutoff = 8000) + 
  ggtitle(label = 'r1',subtitle = '') + 
  scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai1',type = 'continuous',direction = -1)) + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'nFeature')
plot3 <- FeatureOverlay(p2,features = 'nFeature_RNA',add.alpha = F,type = 'raw',pt.size = 2.5,max.cutoff = 8000) + 
  ggtitle(label = 'p2',subtitle = '') + 
  scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai1',type = 'continuous',direction = -1)) + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'nFeature')
plot4 <- FeatureOverlay(r2,features = 'nFeature_RNA',add.alpha = F,type = 'raw',pt.size = 2.5,max.cutoff = 8000) + 
  ggtitle(label = 'r2',subtitle = '') + 
  scale_fill_gradientn(colours = MetBrewer::met.brewer('Hokusai1',type = 'continuous',direction = -1)) + 
  theme(legend.title = element_text(size = 12,face = 'bold',vjust = .9),
        legend.text = element_text(size = 12,face = 'bold'),
        legend.position = 'top',
        legend.key.height = unit(1,'lines'),
        legend.key.width = unit(3,'lines')) + 
  labs(fill = 'nFeature')
pdf('./Figures/Figure1/nFeature.stfeatureplot.pdf',width = 16,height = 5)
plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 4,guides = 'collect')  &
  theme(legend.position='bottom')
dev.off()
