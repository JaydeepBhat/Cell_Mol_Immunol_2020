# Author: Jaydeep Bhat
# Affiliations: 1) Institute of Immunology, Christian-Albrechts-University Kiel & University Hospital Schleswig-Holstein, Campus Kiel, Kiel, Germany; 
# 2) Metabolic Programming, School of Life Sciences Weihenstephan, Technical University Munich (TUM), 85354 Freising, Germany
# Purpose: this script contains codes used in the manuscript
# Please cite this article, if you use this code:
# Bhat, J., Bergmann, A.K., Waschina, S. et al. 
# DNA methylation profile of a hepatosplenic gamma/delta T-cell lymphoma patient associated with response to interferon-α therapy. 
# Cell Mol Immunol (2020). https://doi.org/10.1038/s41423-020-0518-4
# GitHub Repo: https://github.com/JaydeepBhat/Cell_Mol_Immunol_2020

# basic setup
setwd("../HSTL_analysis/")
#try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("GO.db")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ggplot2")
#install.packages("ggplot2") # only if not installed already!
source("../cSplit.R") # depends on dir
library(proxy)
library(ggplot2)
library(cluster)
library(gclus)
library(reshape)
library(reshape2)
library(gplots)
library(splitstackshape)
library(RColorBrewer)
library(pheatmap)
library(topGO)
library(dplyr)
library(GO.db)
library(grid)
library(genomation)
library(GenomicRanges)
library(regioneR)
library(Homo.sapiens)
library(GenomicFeatures)
library(genefilter)
library(Repitools)
library(FactoMineR)
library(factoextra)
library(data.table)

###############
# Import data #
###############
May2019CpGlist <- read.delim(file = "May2019List.txt", header = T)
dim(May2019CpGlist)
AnnotFile <- read.csv("~/Documents/Scripts/HumanMethylation450_15017482_v1-2.csv")
dim(AnnotFile)

# Merge files and rearrange columns
CpGListAnnotated0 <- merge(May2019CpGlist, AnnotFile, "TargetID")
dim(CpGListAnnotated0)
colnames(CpGListAnnotated0)
write.table(CpGListAnnotated0, "CpGListAnnotated_HSTL_full.txt", quote = F, row.names = F, sep = "\t")

# healthy T-cell CpG list
HealthyCpG <- read.delim("../p=0.001.txt")
colnames(HealthyCpG)
HSTL_Healthy0 <- merge(May2019CpGlist,HealthyCpG, "TargetID")
colnames(HSTL_Healthy0)
HSTL_Healthy <- HSTL_Healthy0[,c(1,5:25,39)]
colnames(HSTL_Healthy)
write.table(HSTL_Healthy, "CpGList_HSTL_healthy.txt", quote = F, row.names = F, sep = " ")

# calculate mean of biological replicates and used it further
HSTL_Healthy1 <- cbind(HSTL_Healthy,apply(HSTL_Healthy[,7:22],1,mean))

######################
# Data visualization #
######################
# simple boxplot (just to check!)
CpGListAnnotatedClean <- read.delim("CpGListAnnotated_HSTL.txt", header = T, sep = "/")
boxplot(CpGListAnnotatedClean[c(2:6)], main = "Methylation levels", ylab = "beta values", xlab = "Follow-up", col = c("red","yellow","green","pink","blue"), outline = FALSE)

# violin plot used in the Figure 1 #
P <- subset(CpGListAnnotatedClean[c(2:6)])
Pm <- melt(P)
ggplot(Pm, aes(variable, value), fill=variable) + geom_boxplot()
pg <- ggplot(Pm, aes(variable, value), fill=variable) + geom_violin(fill="skyblue")
pg +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(width=0.15, fill=c("gray")) + scale_y_continuous(breaks=seq(0,1,0.2)) + theme(axis.text.x = element_text(size=14)) + theme(axis.text.y = element_text(size=14))
ggsave("violin_plot_pretransform.pdf",dpi = 300)

# heatmap used in the figure 1
pdf("pearson_heatmap.pdf")
pearson_heatmap <- pheatmap(cor(CpGListAnnotatedClean[,2:6], method=c("pearson")),main = "Pearson", cluster_rows = F, cluster_cols = F, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", legend = TRUE)
print(pearson_heatmap)
dev.off()

# simple scatter plot (to check differences)
Pat_5_1 <- plot(CpGListAnnotatedClean$GD7753, CpGListAnnotatedClean$GD5060, col = adjustcolor("blue", alpha.f=.5), pch = 16, xlim = c(0,1.0), ylim = c(0,1.0), xlab = "Visit 5", ylab = "Visit 1")
abline(v=c(.2, .8), lty = 1, col = "blue")
abline(h=c(.2, .8), lty = 1, col = "blue")

# Principal component analysis and PCA plots
pc_Patient <- princomp(CpGListAnnotatedClean[c(2:6)])
p <- pc_Patient
loadings <- pc_Patient$loadings[]
p.variance.explained = p$sdev^2 / sum(p$sdev^2)
# plot percentage of variance explained for each principal component    
barplot(100*p.variance.explained, las=2, ylim = c(0,100), xlab='', ylab='% Variance', col = "gray25")
# plot by components 1 and 2 (for figure 1)
rownames(CpGListAnnotatedClean) <- CpGListAnnotatedClean[,1]
pc_Patient <- princomp(CpGListAnnotatedClean[c(2:6)])
pdf("PCA_comp1_2.pdf")
PCA_comp1_2 <- fviz_pca_var(pc_Patient, geom.var = c("point", "text"))
print(PCA_comp1_2)
dev.off()
# check scorings (analysis for revision)
dt <- data.table(cpg = CpGListAnnotatedClean[1],
                 gene = CpGListAnnotatedClean[8],
                 pc1.score = pc_Patient$scores[,1],
                 pc2.score = pc_Patient$scores[,2])
dt1 <- dt[order(dt[,abs(dt$pc1.score)],decreasing=T),]
fwrite(dt1, "pc_genes_pc1.txt",quote = F,row.names = F,sep = "\t")
dt2 <- dt[order(dt[,abs(dt$pc2.score)],decreasing=T),]
fwrite(dt2, "pc_genes_pc2.txt",quote = F,row.names = F,sep = "\t")
# plot by components 2 and 3
pdf("PCA_comp2_3.pdf")
PCA_comp2_3 <- fviz_pca_var(pc_Patient, axes = c(2,3), geom.var = c("point", "text"))
print(PCA_comp2_3)
dev.off()

# heatmap representation in Figure 1
# row scale normalized
HSTL <- data.matrix(CpGListAnnotatedClean[,c(2:6)])
rownames(HSTL) <- CpGListAnnotatedClean[,1]
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
pdf("heatmap_normalized_HSTL.pdf")
heatmap_normalized_HSTL <- heatmap.2(HSTL,
                                     scale = "row", #scale = c("none","row", "column")
                                     col = color,
                                     dendrogram = "row", #c("both","row","column","none")
                                     trace = 'none',
                                     density.info = 'none',
                                     Rowv = T,
                                     Colv = F)
print(heatmap_normalized_HSTL)
dev.off()
# none normalized
HSTL <- data.matrix(CpGListAnnotatedClean[,c(2:6)])
rownames(HSTL) <- CpGListAnnotatedClean[,1]
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
pdf("heatmap_nonormalized_HSTL.pdf")
heatmap_nonormalized_HSTL <- heatmap.2(HSTL,
                                       scale = "none", #scale = c("none","row", "column")
                                       col = color,
                                       dendrogram = "row", #c("both","row","column","none")
                                       trace = 'none',
                                       density.info = 'none',
                                       Rowv = T,
                                       Colv = F)
print(heatmap_nonormalized_HSTL)
dev.off()

# genomic features distribution (supplementary figure 2)
fc1 <- data.frame(feature=unlist(strsplit(as.character(CpGListAnnotatedClean1$UCSC_RefGene_Group),";")),test=NA)
fc1$feature <- factor(fc1$feature, levels = c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR"))
ggplot(fc1,aes(feature,fill=feature)) + geom_bar() + 
  theme_bw() + ylim(0,150)+ ylab("No. of CpG sites") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1)
dev.copy(pdf, "UCSC_RefGene_Group_CpG_distribution.pdf")
dev.off()

# CpG distribution analysis (supplementary figure 2)
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
plot(CpGListAnnotatedClean$Relation_to_UCSC_CpG_Island,
     ylim=c(0,160),
     col = ggplotColours(n = 6))
dev.copy(pdf, "Relation_to_UCSC_CpG_Island_all_CpG_distribut.pdf")
dev.off()

## data wrangling to visualize data according to the genomic features
# wide- to - long transformation
source("~/Desktop/NIH_USB_Drive/Scripts/cSplit.R") # depends on location of file
fc2 <- cSplit(CpGListAnnotatedClean1,splitCols = 'UCSC_RefGene_Group',sep = ';') # wide-column
fc3 <- reshape2::melt(fc2,id.vars = colnames(fc2)[1:6],measure.vars = c("UCSC_RefGene_Group_01","UCSC_RefGene_Group_02","UCSC_RefGene_Group_03","UCSC_RefGene_Group_04","UCSC_RefGene_Group_05",
                                                                        "UCSC_RefGene_Group_06","UCSC_RefGene_Group_07","UCSC_RefGene_Group_08","UCSC_RefGene_Group_09","UCSC_RefGene_Group_10",
                                                                        "UCSC_RefGene_Group_11","UCSC_RefGene_Group_12","UCSC_RefGene_Group_13","UCSC_RefGene_Group_14","UCSC_RefGene_Group_15",
                                                                        "UCSC_RefGene_Group_16","UCSC_RefGene_Group_17","UCSC_RefGene_Group_18","UCSC_RefGene_Group_19","UCSC_RefGene_Group_20",
                                                                        "UCSC_RefGene_Group_21","UCSC_RefGene_Group_22","UCSC_RefGene_Group_23","UCSC_RefGene_Group_24","UCSC_RefGene_Group_25",
                                                                        "UCSC_RefGene_Group_26","UCSC_RefGene_Group_27","UCSC_RefGene_Group_28","UCSC_RefGene_Group_29","UCSC_RefGene_Group_30")) # long column
head(fc3) # check!
fc3$value <- factor(fc3$value, levels = c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR"))

# Plots for genomic features analysis for all visits (supplementary figure 2)
# TSS1500 #
fc41 <- subset(fc3, fc3$value == c("TSS1500"), select = c("TargetID","GD5060","GD5239","GD5681","GD6393","GD7753","value"))
head(fc41)
fc410 <- fc41[,c(2:6)]
Pm <- melt(fc410)
pg <- ggplot(Pm, aes(variable, value), fill=variable) + geom_violin(fill="skyblue")
pg +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(width=0.15, fill=c("gray")) + scale_y_continuous(breaks=seq(0.0,1.0,0.25)) + theme(axis.text.x = element_text(size=14)) + theme(axis.text.y = element_text(size=14))
ggsave("violin_plot_TSS1500.pdf",dpi = 300)
# TSS200 #
fc42 <- subset(fc3, fc3$value == c("TSS200"), select = c("TargetID","GD5060","GD5239","GD5681","GD6393","GD7753","value"))
head(fc42)
fc420 <- fc42[,c(2:6)]
Pm <- melt(fc420)
pg <- ggplot(Pm, aes(variable, value), fill=variable) + geom_violin(fill="skyblue")
pg +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(width=0.15, fill=c("gray")) + scale_y_continuous(breaks=seq(0.0,1.0,0.25)) + theme(axis.text.x = element_text(size=14)) + theme(axis.text.y = element_text(size=14))
ggsave("violin_plot_TSS200.pdf",dpi = 300)
# 5'UTR #
fc43 <- subset(fc3, fc3$value == c("5'UTR"), select = c("TargetID","GD5060","GD5239","GD5681","GD6393","GD7753","value"))
head(fc43)
fc430 <- fc43[,c(2:6)]
Pm <- melt(fc430)
pg <- ggplot(Pm, aes(variable, value), fill=variable) + geom_violin(fill="skyblue")
pg +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(width=0.15, fill=c("gray")) + scale_y_continuous(breaks=seq(0.0,1.0,0.25)) + theme(axis.text.x = element_text(size=14)) + theme(axis.text.y = element_text(size=14))
ggsave("violin_plot_5UTR.pdf",dpi = 300)
# 1stExon #
fc44 <- subset(fc3, fc3$value == c("1stExon"), select = c("TargetID","GD5060","GD5239","GD5681","GD6393","GD7753","value"))
head(fc44)
fc440 <- fc44[,c(2:6)]
Pm <- melt(fc440)
pg <- ggplot(Pm, aes(variable, value), fill=variable) + geom_violin(fill="skyblue")
pg +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(width=0.15, fill=c("gray")) + scale_y_continuous(breaks=seq(0.0,1.0,0.25)) + theme(axis.text.x = element_text(size=14)) + theme(axis.text.y = element_text(size=14))
ggsave("violin_plot_1stExon.pdf",dpi = 300)
# Body #
fc45 <- subset(fc3, fc3$value == c("Body"), select = c("TargetID","GD5060","GD5239","GD5681","GD6393","GD7753","value"))
head(fc45)
fc450 <- fc45[,c(2:6)]
Pm <- melt(fc450)
pg <- ggplot(Pm, aes(variable, value), fill=variable) + geom_violin(fill="skyblue")
pg +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(width=0.15, fill=c("gray")) + scale_y_continuous(breaks=seq(0.0,1.0,0.25)) + theme(axis.text.x = element_text(size=14)) + theme(axis.text.y = element_text(size=14))
ggsave("violin_plot_Body.pdf",dpi = 300)
# 3'UTR
fc46 <- subset(fc3, fc3$value == c("3'UTR"), select = c("TargetID","GD5060","GD5239","GD5681","GD6393","GD7753","value"))
head(fc46)
fc460 <- fc46[,c(2:6)]
Pm <- melt(fc460)
pg <- ggplot(Pm, aes(variable, value), fill=variable) + geom_violin(fill="skyblue")
pg +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(width=0.15, fill=c("gray")) + scale_y_continuous(breaks=seq(0.0,1.0,0.25)) + theme(axis.text.x = element_text(size=14)) + theme(axis.text.y = element_text(size=14))
ggsave("violin_plot_3UTR.pdf",dpi = 300)

# Plots for genomic features analysis for all visits
fc30 <- fc3[!is.na(fc3$value),] # remove NA
ggplot(fc30,aes(x=value,y=GD5060)) +
  stat_summary(aes(group=1),fun.y="median", geom="line", col='red') + theme_bw() + ylim(0.4,0.8) + ylab("Beta value") + xlab("features") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_summary(aes(x=value,y=GD5239,group=1),fun.y = "median",geom = "line",col='green') +
  stat_summary(aes(x=value,y=GD5681,group=1),fun.y = "median",geom = "line",col='dimgray') +
  stat_summary(aes(x=value,y=GD6393,group=1),fun.y = "median",geom = "line",col='magenta') +
  stat_summary(aes(x=value,y=GD7753,group=1),fun.y = "median",geom = "line",col='blue')
ggsave("genomic_features_all-vis.pdf",dpi = 300)

# Plots for genomic features analysis for visit 1 and 5 (supplementary figure 2)
fc31 <- fc3[!is.na(fc3$value),] # remove NA
ggplot(fc31,aes(x=value,y=GD5060)) +
  stat_summary(aes(group=1),fun.y="median", geom="line", col='red') + theme_bw() + ylim(0.4,0.8) + ylab("Beta value") + xlab("features") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_summary(aes(x=value,y=GD7753,group=1),fun.y = "median",geom = "line",col='blue') + theme(axis.text.x = element_text(size=14)) + theme(axis.text.y = element_text(size=14)) + theme(aspect.ratio = 1)
#ggsave("genomic_features_vis15.pdf",dpi = 300)

# use merged healthy and HSTL data
# heatmaps representation (supplementary figure 2)
HSTL_Healthy2 <- HSTL_Healthy1[c(1:6,24,23)]
HSTL_Healthy3 <- data.matrix(HSTL_Healthy2[,c(2:7)])
rownames(HSTL_Healthy3) <- HSTL_Healthy2[,1]
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
pdf("heatmap_normalized_HSTL_HealthyMean.pdf")
heatmap_normalized_HSTL_HealthyMean <- heatmap.2(HSTL_Healthy3,
                                                 scale = "row", #scale = c("none","row", "column")
                                                 col = color,
                                                 dendrogram = "both", #c("both","row","column","none")
                                                 trace = 'none',
                                                 density.info = 'none',
                                                 Rowv = T,
                                                 Colv = T)
print(heatmap_normalized_HSTL_HealthyMean)
dev.off()

########################
# Statistical analysis #
########################
# Data Summary
Sum_Data <- summary(CpGListAnnotatedClean)
capture.output(Sum_Data, file = "Sum_Data_HSTL.txt")
Sum_Data
# significance test between visit 1 and 5
Wilcox_Signi_pat_1_5 <- wilcox.test(CpGListAnnotatedClean$GD5060,CpGListAnnotatedClean$GD7753, paired = T)
capture.output(Wilcox_Signi_pat_1_5, file = "Wilcox_Signi_HSTL_1_5.txt")
Wilcox_Signi_pat_1_5
pvals15 <- c()
pvals15 <- c(pvals15,Wilcox_Signi_pat_1_5$p.value)
pvals15 <- p.adjust(pvals15,method = 'fdr')
capture.output(Wilcox_Signi_pat_1_5, file = "Wilcox_Signi_HSTL_1_5_with_adj_pval.txt")
pvals15
# significance test between genomic features of visit 1 and 5
pvals1 <- c()
ind <- fc3$value=="TSS1500"
wcp <- wilcox.test(fc3$GD5060[ind], fc3$GD7753[ind], paired = T)
pvals1 <- c(pvals1,wcp$p.value)
ind <- fc3$value=="TSS200"
wcp <- wilcox.test(fc3$GD5060[ind], fc3$GD7753[ind], paired = T)
pvals1 <- c(pvals1,wcp$p.value)
ind <- fc3$value=="5'UTR"
wcp <- wilcox.test(fc3$GD5060[ind], fc3$GD7753[ind], paired = T)
pvals1 <- c(pvals1,wcp$p.value)
ind <- fc3$value=="1stExon"
wcp <- wilcox.test(fc3$GD5060[ind], fc3$GD7753[ind], paired = T)
pvals1 <- c(pvals1,wcp$p.value)
ind <- fc3$value=="Body"
wcp <- wilcox.test(fc3$GD5060[ind], fc3$GD7753[ind], paired = T)
pvals1 <- c(pvals1,wcp$p.value)
ind <- fc3$value=="3'UTR"
wcp <- wilcox.test(fc3$GD5060[ind], fc3$GD7753[ind], paired = T)
pvals1 <- c(pvals1,wcp$p.value)
pvals1 <- p.adjust(pvals1,method = 'fdr')
capture.output(pvals1, file = "fdr_pvals_HSTL_1_5.txt")
pvals1
# Correlation analysis
cor_all <- cor(CpGListAnnotatedClean[,2:6], method = c("pearson"))
capture.output(cor_all, file = "cor_all_visits_HSTL.txt")
cor_all

session_info()
