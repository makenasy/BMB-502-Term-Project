# BMB 502 Term Project
# DESeq2

#------------------------------------------------------------------------------#
# Ian's DESeq2 Workflow
library("DESeq2")

# Read in Galaxy hiseq counts table
htseq_counts <- read.table("htseq_counts_SRR_labeled.txt",header=TRUE,sep="\t",row.names = 1)

# Change this to match the number of samples
# rnaseqMatrix <- htseq_counts[,c(1:12)]
# rownames(rnaseqMatrix) <- htseq_counts[,1]
# head(rnaseqMatrix)

rnaseqMatrix <- htseq_counts[1:12]
rownames(rnaseqMatrix) <- htseq_counts[,1]
colnames(rnaseqMatrix) <- c("Emp_NoDox_1","Dux_NoDox_1","Emp_Dox_1","Dux_Dox_1","Emp_NoDox_2","Dux_NoDox_2","Emp_Dox_2","Dux_Dox_2","Emp_NoDox_3","Dux_NoDox_3","Emp_Dox_3","Dux_Dox_3")
head(rnaseqMatrix)

# Define the sample mappings
samples <- data.frame(matrix(c("Emp_NoDox_1","Dux_NoDox_1","Emp_Dox_1","Dux_Dox_1",
                               "Emp_NoDox_2","Dux_NoDox_2","Emp_Dox_2","Dux_Dox_2",
                               "Emp_NoDox_3","Dux_NoDox_3","Emp_Dox_3","Dux_Dox_3",
                               "Emp_NoDox","Dux_NoDox","Emp_Dox","Dux_Dox",
                               "Emp_NoDox","Dux_NoDox","Emp_Dox","Dux_Dox",
                               "Emp_NoDox","Dux_NoDox","Emp_Dox","Dux_Dox"),ncol=2))

names(samples) <- c("ID","Treatment")
rownames(samples) <- samples[,1]
samples$Genotype <- factor(c("Emp_NoDox","Dux_NoDox","Emp_Dox","Dux_Dox",
                               "Emp_NoDox","Dux_NoDox","Emp_Dox","Dux_Dox",
                               "Emp_NoDox","Dux_NoDox","Emp_Dox","Dux_Dox"))

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData = rnaseqMatrix,
                                     colData = samples,
                                     design = ~ Treatment)

dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 10, ])

deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 10, ]

deseq2Data <- DESeq(deseq2Data)


rld <- rlog(deseq2Data, blind=FALSE)
rlogcounts <- data.frame(assay(deseq2Data))
rownames(rlogcounts) <- rownames(deseq2Data)

res_Dux_On_Off <- results(deseq2Data, contrast=c("Treatment", "Dux_Dox", "Dux_NoDox"))
resOrdered_Dux_On_Off  <- res_Dux_On_Off [order(res_Dux_On_Off$pvalue),]

results <- data.frame(resOrdered_Dux_On_Off)
head(results)

write.table(results,"Dux_On_Off_DESeq2_paper.txt",sep="\t",row.names=TRUE)

#------------------------------------------------------------------------------#
# Makena's DESeq2 Workflow

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
plotMA(results)

# PCA plot
plotPCA(rlog(dds), intgroup = "Sample_Group")
install.packages("ggplot2")
library(ggplot2)

# May or may not need
# install.packages("ggfortify")
# library("ggfortify")

# Save the results as a data frame
results <- data.frame(res)
write.table(results,"FSHD_DESeq2_gene_symbols.txt",sep="\t",row.names=TRUE)
#------------------------------------------------------------------------------#
# VOLCANO PLOT

# OPTION 1: Easy Workflow, difficulty level = 5
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
# results <- as.dataframe(res)
results$gene <- mapIds(org.Mm.eg.db, keys = rownames(results), keytype = "ENSEMBL", column = "SYMBOL")

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano) 

library(gridExtra)
library(grid)

mypval=0.001
myfc=2 
mypadj=0.001
mylog2fc=log2(myfc)

# Labeling top 10 Up and Down regulated genes
top <- 10
top_genes <- bind_rows(
  results %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top),
  results %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top)
)
top_genes %>% 
  knitr_table()

write.table(top_genes,"FSHD_TP_top_10_gene_expression.txt",sep="\t",row.names=TRUE)

p1 <- EnhancedVolcano(results,
                lab = results$gene, 
                x = "log2FoldChange", 
                y = "padj",
                pCutoff = mypadj,
                FCcutoff = mylog2fc,
                labSize = 3.0,
                colAlpha = 1,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                legendPosition = 'top',
                legendLabSize = 10,
                legendIconSize = 5.0,
                drawConnectors = F, # Add lines to show more labels
                widthConnectors = 0.2,
                colConnectors = 'grey30') 
p1

p2 <- EnhancedVolcano(top_genes,
                lab = top_genes$gene, 
                x = "log2FoldChange", 
                y = "padj",
                pCutoff = mypadj,
                FCcutoff = mylog2fc,
                labSize = 3.0,
                colAlpha = 1,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                legendPosition = 'top',
                legendLabSize = 10,
                legendIconSize = 5.0,
                drawConnectors = F, # Add lines to show more labels
                widthConnectors = 0.2,
                colConnectors = 'grey30') 

p2

grid.arrange(p1, p2, ncol=2)

#------------------------------------------------------------------------------#
# OPTION 2: Painful Workflow, difficulty level = no sleep
# Load libraires
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(tidyverse)
library(ggrepel)

# Add gene symbols to results
results$gene <- mapIds(org.Mm.eg.db, keys = rownames(results), keytype = "ENSEMBL", column = "SYMBOL")

# A short function for outputting the tables
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

# A simple volcano plot
# p1 <- ggplot(results, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion  
#   geom_point(size = 2/5) +
#   xlab(expression("log"[2]*"log2FoldChange")) + 
#   ylab(expression("-log"[10]*"padj"))
# p1

# Adding color to differentially expressed genes (DEGs)
results <- results %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(results) %>% 
  knitr_table()

# Conversion of the FDR values to their -log10 can be done at this step
p1 <- ggplot(results, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"Fold Change")) + 
  ylab(expression("-log"[10]*"FDR"))

p1

# Adding color to DEGs
results <- results %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(results) %>% 
  knitr_table()

# Map the column ´Expression’ to the color aesthetic of geom_point() and color the points according to their expression classification
p2 <- ggplot(results, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"Fold Change")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("hotpink", "gray50", "cyan")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))

p2

# To know how many genes are up- or down-regulated, or unchanged, we can use dplyr’s count() function
results %>% 
  count(Expression) %>% 
  knitr_table()

# Since we already know that the genes towards the right are up-regulated 
# and the genes towards the left are down-regulated
# it would be more informative if we colored the points according to their significance level instead. 
# Let’s create another column, named ‘Significance'
# classify the genes according to significance thresholds (0.05, 0.01, and 0.001)
results <- results %>% 
  mutate(
    Significance = case_when(
      abs(log2FoldChange) >= log(2) & padj <= 0.05 & padj > 0.01 ~ "FDR 0.05", 
      abs(log2FoldChange) >= log(2) & padj <= 0.01 & padj > 0.001 ~ "FDR 0.01",
      abs(log2FoldChange) >= log(2) & padj <= 0.001 ~ "FDR 0.001", 
      TRUE ~ "Unchanged")
  )
head(results) %>% 
  knitr_table()

# Use the color aesthetic to map the color of the points to their corresponding significance thresholds
p3 <- ggplot(results, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"Fold Change")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3

# Count how many genes are up- or down-regulated according to the different significance thresholds with count()
results %>% 
  count(Expression, Significance) %>% 
  knitr_table()

# Labeling top 10 Up and Down regulated genes
top <- 10
top_genes <- bind_rows(
  results %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top),
  results %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top)
)
top_genes %>% 
  knitr_table()

# Write top and bottom expressed genes to table 
write.table(top_genes,"FSHD_TP_top_10_gene_expression.csv",sep="\t",row.names=TRUE)


# Plot with labels and Significance
p3 <-  p3 +
  geom_label_repel(results = top_genes,
                   mapping = aes(log2FoldChange, -log(padj,10), label = symbol),
                   size = 2)
p3

# Fixing "too many overlaps" error, for some reason you need to use this after running code above (not before)
options(geomrepel.max.overlaps = Inf)

# Plot with labels and Expression
p2 <-  p2 +
  geom_label_repel(results = top_genes,
                   mapping = aes(log2FoldChange, -log(padj,10), label = gene),
                   size = 2)
p2

# Fixing "too many overlaps" error, for some reason you need to use this after running code above (not before)
options(geomrepel.max.overlaps = Inf)

#------------------------------------------------------------------------------#
# OPTION 3: Kinda works...
library(ggplot2)
library(dplyr)

top_genes <- results %>% 
  group_by(Expression) %>% 
  dplyr::top_n(10, wt = padj)

ggplot(results, aes(x = log2FoldChange, y = padj, col = Expression, label = gene)) + 
  geom_point() + 
  ggrepel::geom_text_repel(show.legend = FALSE) + 
  #geom_text(vjust = -.1, show.legend = FALSE) + 
  geom_point(size = 3, show.legend = FALSE) + 
  theme_minimal() + 
  scale_color_manual(values = c("hotpink", "gray", "cyan")) + 
  geom_vline(xintercept=c(-1.6, 1.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") 

options(geomrepel.max.overlaps = Inf)

# end of code
#------------------------------------------------------------------------------#

