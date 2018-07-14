source("functions.R")

# Using the public ERCC data at:
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/ercc

library(Matrix)
counts = readMM("ercc/ercc92/matrix.mtx")
genes = read.table("ercc/ercc92/genes.tsv")
rownames(counts) = genes[,1]
keep = rowMeans(counts) > 0.1
counts = counts[keep,]

MGC = run_MAGIC(as.matrix(counts), t=10)

# Performing t-SNE on the result, compared to the original.
set.seed(2000)
library(Rtsne)
mgc_TSNE = Rtsne(MGC, pca=FALSE, perplexity=90)

lcounts = lognormalize(as.matrix(counts))
normal_TSNE = Rtsne(lcounts, pca=FALSE, perplexity=90)

# Creating plots.
pdf("pics/ercc_results.pdf")
col = ifelse(mgc_TSNE$Y[,1] > -5, "forestgreen", "goldenrod")
plot(mgc_TSNE$Y[,1], mgc_TSNE$Y[,2], xlab="TSNE1", ylab="TSNE2", pch=16, col=col)
plot(normal_TSNE$Y[,1], normal_TSNE$Y[,2], xlab="TSNE1", ylab="TSNE2", pch=16, col=col)
dev.off()
