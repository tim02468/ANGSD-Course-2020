library(ggplot2)
library(DESeq2)
library(magrittr)

df <- read.table("../data/featureCounts_human.txt", header = TRUE, stringsAsFactors = FALSE)
head(df)

names(df) <- c(names(df)[1:6], sprintf("MSDP%s", seq(1:10)), sprintf("SDP%s", seq(1:10)), sprintf("SDPC%s", seq(1:10)))
# names(df) <- c(names(df)[1:6], sprintf("MSDP%s", seq(1:20)), sprintf("SDPC%s", seq(1:10)))
row.names(df) <- make.names(df$Geneid)

df <- df[c(names(df)[1:6], sprintf("SDPC%s", seq(1:10)), sprintf("MSDP%s", seq(1:10)))]
# df <- df[c(names(df)[1:6], sprintf("SDPC%s", seq(1:10)), sprintf("SDP%s", seq(1:10)))]

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
# pcaExplorer::pcaExplorer(dds = DESeq.df, dst = DESeq.rlog)


# run DESeq test
DESeq.df <- DESeq(DESeq.df)

DESeq.df$condition <- relevel(DESeq.df$condition, ref="SDPC")


DGE.results <- results(DESeq.df, independentFiltering = TRUE, alpha = 0.05)
# the first line will tell you which comparison was done to achieve the log2FC
head(DGE.results)


hist(DGE.results$padj,
     col="grey", border="white", xlab="", ylab="", main="frequencies of adj. p-values\n(all genes)", cex = 0.4)


DGEgenes <- rownames(subset(DGE.results, padj < 0.05))



# extract rlog-transformed values into a matrix
rlog.dge <- DESeq.rlog[DGEgenes,] %>% assay

# show dendrogram and pca using filtered gene
corr_coeff <- cor(rlog.dge, method = "pearson")
rlog.norm.counts <- assay(DESeq.df[DGEgenes,], "log.norm.counts")
as.dist(1-corr_coeff) %>% hclust %>% plot( ., labels = colnames(rlog.norm.counts),main = "rlog transformed read counts")

pca_data <- plotPCA(DESeq.rlog[DGEgenes,], returnData=TRUE)

p <- ggplot(pca_data,aes(x=PC1,y=PC2,color=group))
p <- p+geom_point()
p

pcaExplorer::pcaExplorer(dds = DESeq.df[DGEgenes,], dst = DESeq.rlog[DGEgenes,])


library(pheatmap)

pheatmap(rlog.dge, scale="none",
         show_rownames = FALSE, main = "DGE (no scaling)")

pheatmap(rlog.dge, scale="row",
         show_rownames = FALSE, main = "DGE (row-based z-score)")

plotMA(DGE.results, alpha = 0.05,
       main = "Test: p.adj.value < 0.05", ylim = c(-4,4))

library(EnhancedVolcano)
vp1 <- EnhancedVolcano(DGE.results,
                       lab = rownames(DGE.results), x = 'log2FoldChange',
                       y = 'padj', pCutoff = 0.05,
                       title = "SDPC / MSDP") 
print(vp1)




DGE.results.shrnk <- lfcShrink(DESeq.df, coef = 2, type = "apeglm")
par(mfrow = c(1,1))

#plotMA(DGE.results, alpha = 0.05, main = "no shrinkage", ylim = c(-4,4))

plotMA(DGE.results.shrnk,
       alpha = 0.05,
       main = "with logFC shrinkage", ylim = c(-4,4))


vp2 <- EnhancedVolcano(DGE.results.shrnk, lab = rownames(DGE.results.shrnk),
                       x = 'log2FoldChange',
                       y = 'padj', pCutoff = 0.05,
                       title = "with logFC shrinkage")
library(patchwork)
vp2


library(TxDb.Hsapiens.UCSC.hg38.knownGene)

keytypes(TxDb.Hsapiens.UCSC.hg38.knownGene)
# list columns that can be retrieved from the annotation data base
columns(TxDb.Hsapiens.UCSC.hg38.knownGene)



anno.DGE <- select(TxDb.Hsapiens.UCSC.hg38.knownGene,
                   keys = DGEgenes, 
                   keytype="GENEID", 
                   columns=c("EXONNAME")) 




keys(TxDb.Hsapiens.UCSC.hg38.knownGene)


