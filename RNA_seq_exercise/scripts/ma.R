#usr/bin/env Rscript

## R script to create ma plot

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

p <- ggplot(table_res,aes(x=baseMean,y=log2FoldChange,colour=sig))+
	geom_point()+
	scale_x_log10()+
	scale_colour_manual(values=c("grey40","#80b1d3"))+
	labs(x="Log10 Mean expression",y="Log2 Fold Change")+
	theme_bw()+
	theme(legend.title=element_blank(),
	      legend.position="top",
	      legend.justification="right")

ggsave("ma.png",p,height=12,width=12,units="cm",dpi=250)


