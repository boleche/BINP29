# loading in required packages
library(magrittr)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# load the count matrix and save gene names
rawData <- read.csv("/home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/results/star_salmon/salmon.merged.gene_counts.tsv",
    sep = "\t", row.names = "gene_id")

# save gene names for later
geneNames <- data.frame(gene_id = rownames(rawData),
                        gene_name = rawData$gene_name)

# create count matrix
countMatrix <- rawData %>%
    dplyr::select(-gene_name) %>%    # remove gene_name column
    apply(2, as.numeric) %>%         # convert strings to numbers
    round()

# restore row names (apply drops them)
rownames(countMatrix) <- geneNames$gene_id

# check the column names of the count matrix
print(colnames(countMatrix))
print(head(countMatrix))
print(dim(countMatrix))
print(class(countMatrix))

# define the sample names (must match column names in count matrix exactly)
sampleNames <- c("SRR3222409_KO",
                 "SRR3222410_KO",
                 "SRR3222411_KO",
                 "SRR3222412_WT",
                 "SRR3222413_WT",
                 "SRR3222414_WT")
sampleNames # check


# define the sample conditions
sampleConditions <- c(rep("KO", 3),
                      rep("WT", 3))
sampleConditions # check

# create a sample table from the sample conditions
sampleTable <- data.frame(condition = as.factor(sampleConditions))
row.names(sampleTable) <- sampleNames


# create the DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(countMatrix),
                              colData = sampleTable,
                              design = ~ condition)

# running the DESeq analysis
dds <- DESeq(dds)

# check which contrasts are availale
DESeq2::resultsNames(dds)
# you can see that condition_Urtica_vs_Ribes is available

# set up the results table (KO vs WT)
res_table <- results(dds, name="condition_WT_vs_KO")

# add gene names to results
res_table$gene_id <- rownames(res_table)
res_table <- merge(as.data.frame(res_table), geneNames, by = "gene_id")

# reorder so gene_name is next to gene_id
res_table <- res_table[, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# re-sort by padj after merge scrambled the order
res_table <- res_table[order(res_table$padj),]

# creating a subset with out log2fc max of 1 and a padj value of 0.01
resSig <- subset(res_table, padj < 0.05 & abs(log2FoldChange) > 1)

# write the table to output
write.table(resSig, file = "/home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/R_analysis/diffExpr.sig01.tab", sep = "\t", quote = FALSE, row.names = TRUE)

# read in the subset table ( and treat the first line as a header )
subset_table <- read.csv("diffExpr.sig01.tab", header = TRUE)

# count the rows
print(nrow(subset_table))

#############################################################

# ---- VST TRANSFORMATION ----
# apply VST transformation
vst_counts <- varianceStabilizingTransformation(dds, blind = FALSE)

# extract VST matrix
vst_matrix <- assay(vst_counts)

# add gene names to vst matrix
vst_output <- as.data.frame(vst_matrix)
vst_output$gene_id <- rownames(vst_output)
vst_output <- merge(data.frame(gene_id = rownames(vst_output), vst_output), 
                    geneNames, by = "gene_id")

# reorder so gene_name is next to gene_id
vst_output <- vst_output[, c("gene_id", "gene_name", sampleNames)]

# write VST matrix to file
write.table(vst_output, 
            file = "/home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/R_analysis/vst_counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("VST counts saved!\n")

#############################################################

# get all results (not just significant) for the volcano plot
all_results <- results(dds, name="condition_WT_vs_KO")
all_results$gene_id <- rownames(all_results)
all_results <- merge(as.data.frame(all_results), geneNames, by = "gene_id")

# remove NAs
all_results <- all_results[!is.na(all_results$padj) & !is.na(all_results$log2FoldChange), ]

# add significance column
all_results$sig <- ifelse(all_results$padj < 0.05, "Sig", "NotSig")

# volcano plot
p <- ggplot(all_results, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) +
    geom_point() +
    scale_colour_manual(values = c("grey40", "#80b1d3")) +
    labs(x = "Log2 Fold Change",
         y = "-Log10 BH adjusted p-value",
         title = "KO vs WT Volcano Plot") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.justification = "right")

# save the plot
ggsave("/home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/R_analysis/volcano_plot.png",
       p, height = 12, width = 12, units = "cm", dpi = 250)

cat("Volcano plot saved!\n")
#############################################################

# create metadata for heatmap
table_meta <- data.frame(
    accession = c("SRR3222409_KO", "SRR3222410_KO", "SRR3222411_KO",
                  "SRR3222412_WT", "SRR3222413_WT", "SRR3222414_WT"),
    condition = c(rep("KO", 3), rep("WT", 3)),
    replicate = rep(1:3, 2),
    stringsAsFactors = FALSE)
table_meta$condition <- factor(table_meta$condition, levels = c("WT", "KO"))
table_meta$replicate <- as.factor(table_meta$replicate)
rownames(table_meta) <- table_meta$accession

# define colours
col_cond <- c("#fb8072", "#8dd3c7")
names(col_cond) <- levels(table_meta$condition)

col_rep <- c("#80b1d3", "#fdb462", "#bebada")
names(col_rep) <- levels(table_meta$replicate)

# get top 50 genes with lowest pvalue (from DGE)
res_pval <- results(dds, name="condition_WT_vs_KO")
res_pval <- res_pval[order(res_pval$pvalue),]
top50_genes <- head(rownames(res_pval), 50)

# get gene names for top 50 IN THE SAME ORDER as top50_genes
top50_geneNames <- geneNames[match(top50_genes, geneNames$gene_id), ]

# extract VST counts for top 50 genes
top50_counts <- vst_matrix[top50_genes, ]

# replace ensembl ids with gene names as row labels
rownames(top50_counts) <- top50_geneNames$gene_name
# verify gene order is correct
all.equal(rownames(top50_counts), top50_genes)

# replace ensembl ids with gene names as row labels
rownames(top50_counts) <- top50_geneNames$gene_name

# match column order of counts to metadata
mth <- match(colnames(top50_counts), rownames(table_meta))
top50_counts <- top50_counts[, mth]

# verify columns match metadata
all.equal(colnames(top50_counts), rownames(table_meta))

# plot and save heatmap
png("/home/inf-25-2025/Desktop/binp29/RNA_seq_exercise/R_analysis/heatmap.png",
    height = 22, width = 12, units = "cm", res = 300)
pheatmap::pheatmap(top50_counts,
                   annotation_col = table_meta[, c("condition", "replicate")],
                   annotation_colors = list(condition = col_cond,
                                            replicate = col_rep),
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   scale = "row",
                   border_color = NA,
                   fontsize_row = 8,
                   main = "Top 50 DE Genes by P-value (KO vs WT)")
dev.off()
cat("Heatmap saved!\n")