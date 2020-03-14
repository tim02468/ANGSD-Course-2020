library(ggplot2)
library(DESeq2)
library(magrittr)

df <- read.table("../data/featureCounts_human.txt", header = TRUE, stringsAsFactors = FALSE)
head(df)

names(df) <- c(names(df)[1:6], "MSDP075", "MSDP080", "SDPC082", "SDPC087")
row.names(df) <- make.names(df$Geneid)

df <- df[ ,-c(1:6)]
col_data <- DataFrame(condition = gsub("[0-9]+", "", names(df)), row.names = names(df))

# Create DESeq object
DESeq.df <- DESeqDataSetFromMatrix(countData = df, colData = col_data, design = ~ condition)
dim(DESeq.df)

# remove genes that have no counts
DESeq.df <- DESeq.df[rowSums(counts(DESeq.df)) > 0, ]
dim(DESeq.df)

# calculate SFs and add them to the object
DESeq.df <- estimateSizeFactors(DESeq.df)

plot(sizeFactors(DESeq.df), colSums(counts(DESeq.df)), ylab = "library size", xlab = "size factors")

# par(mfrow = c(1,1))
# boxplot(log2(counts(DESeq.df) + 1), notch = TRUE, main = "Non-normalized read counts", ylab = "log2(read counts)")
# boxplot(log2(counts(DESeq.df, normalize = TRUE) + 1), notch = TRUE, main = "Size-Factor normalized read counts", ylab = "log2(read counts)")

# assign log normalized counts to a matrix
assay(DESeq.df, "log.norm.counts") <- log2(counts(DESeq.df, normalized = TRUE) + 1)

DESeq.rlog <- rlog(DESeq.df, blind = TRUE)

# plot dendrogram
corr_coeff <- cor(assay(DESeq.df, "log.norm.counts"), method = "pearson")
rlog.norm.counts <- assay(DESeq.df, "log.norm.counts")
as.dist(1-corr_coeff) %>% hclust %>% plot( ., labels = colnames(rlog.norm.counts),main = "rlog transformed read counts")

# run pca and plot the result
pca_data <- plotPCA(DESeq.rlog, returnData=TRUE)

p <- ggplot(pca_data,aes(x=PC1,y=PC2,color=group))
p <- p+geom_point()
p

# run pcaExplorer
pcaExplorer::pcaExplorer(dds = DESeq.df, dst = DESeq.rlog)

