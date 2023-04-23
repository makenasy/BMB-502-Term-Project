# BMB 502 Term Project
# DESeq2

# install BiocManager package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# install DESeq2 package
BiocManager::install("DESeq2")

# load libraries
library("DESeq2")
#library(tidyverse)
library(airway)

# Read in Galaxy hiseq counts table
htsq_counts <- read.table("htseq_counts_SRR_labeled.txt",header=TRUE,sep="\t",row.names = 1)
#rownames <- htsq_counts[,1]

# Read in sample info
colData <- read.csv("FSHD_col_data.csv",row.names = 1)
#colnames <- colData[0:1]

# make sure row names in colData match to column names in htsq_counts
# all(colnames(htsq_counts) %in% rownames(colData))

# are they in the same order?
# all(colnames(htsq_counts) == rownames(colData))

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

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

# MA plot
plotMA(res)

plotPCA(res)
# Install ggfortify for PCA plot
# library(pheatmap)
# library(apeglm)
install.packages("ggplot2")
library(ggplot2)
install.packages("ggfortify")
library("ggfortify")


#----------------------------------------------------------------------------------------------------#

