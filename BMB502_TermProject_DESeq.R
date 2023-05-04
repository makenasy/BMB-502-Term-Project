# BMB 502 Term Project
# DESeq2

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
library(tidyverse)

# Read in Galaxy hiseq counts table
htsq_counts <- read.table("htseq_counts_SRR_labeled.txt",header=TRUE,sep="\t",row.names = 1)


# Read in sample info
colData <- read.csv("FSHD_col_data.csv",row.names = 1)


# make sure row names in colData match to column names in htsq_counts
all(colnames(htsq_counts) %in% rownames(colData))

# are they in the same order?
all(colnames(htsq_counts) == rownames(colData))

# Construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = htsq_counts,
                       colData = colData,
                       design = ~ doxycycline)

# removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set factor level
# dds$doxycycline <- relevel(dds$doxycycline, ref = "empty_no_treatment")

# NOTE: collapse technical replicates

# Run DESeq
dds <- DESeq(dds)
res <- results(dds)

# Explore Results
summary(res)

# Different p-value
# res0.01 <- results(dds, alpha = 0.01)
# summary(res0.01)

# MA plot
plotMA(res)

# PCA plot
plotPCA(rlog(dds), intgroup = "Sample_Group")
install.packages("ggplot2")
library(ggplot2)

# May or may not need
# install.packages("ggfortify")
# library("ggfortify")

# Save the results as a data frame
results <- data.frame(res)

write.table(results,"FSHD_DESeq2.txt",sep="\t",row.names=TRUE)

BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
# results <- as.dataframe(res)
results$symbol <- mapIds(org.Mm.eg.db, keys = rownames(results), keytype = "ENSEMBL", column = "SYMBOL")

# Volcano plot
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

EnhancedVolcano(results, x = "log2FoldChange", y = "padj", lab = results$symbol)
                # pCutoff = 1e-4, FCcutoff = 1)

#------------------------------------------------------------------------------#

