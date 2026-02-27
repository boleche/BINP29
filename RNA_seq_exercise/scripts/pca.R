#!/usr/bin/env Rscript

## R script to create pca plot

# only applicable when running on rackham cluster
# run before running R
# module load R/4.0.0
# module load R_packages/4.0.0
## .libPaths("/sw/apps/R_packages/4.0.0/rackham")

library(ggplot2)

# load deg data
table_res <- readRDS("../5_dge/dge_results_full.Rds")
table_res <- table_res[!is.na(table_res$padj),]
table_res$sig <- ifelse(table_res$padj<0.05,"Sig","NotSig")

# load vst count data
table_vst <- readRDS("../5_dge/counts_vst_full.Rds")

# create metadata
table_meta <- data.frame(accession=c("SRR3222409","SRR3222410","SRR3222411","SRR3222412","SRR3222413","SRR3222414"),
condition=c(rep(c("KO","Wt"),each=3)),replicate=rep(1:3,2),stringsAsFactors=F)
table_meta$condition <- factor(table_meta$condition,levels=c("Wt","KO"))
table_meta$replicate <- as.factor(table_meta$replicate)
rownames(table_meta) <- table_meta$accession

# match order of counts and metadata
mth <- match(colnames(table_vst),rownames(table_meta))
table_vst <- table_vst[,mth]
all.equal(rownames(table_meta),colnames(table_vst))

# pca
pcaobj <- prcomp(x=t(table_vst))
pcs <- round(pcaobj$sdev^2/sum(pcaobj$sdev^2)*100,2)

pcamat1 <- as.data.frame(pcaobj$x)
pcamat2 <- merge(pcamat1,table_meta,by=0)

p <- ggplot(pcamat2,aes(PC1,PC2,colour=condition))+
	geom_point()+
	geom_text(aes(label=accession),size=3,nudge_x=1,hjust="inward")+
	theme_bw()+
	theme(legend.title=element_blank(),
	      legend.position="top",
	      legend.justification="right")

ggsave("pca.png",p,height=12,width=12,units="cm",dpi=250)

