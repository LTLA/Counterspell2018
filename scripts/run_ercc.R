# This script tests the behaviour of MAGIC on the public ERCC data from 10X Genomics.

source("functions.R")

##############################################

library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
fname <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/ercc/ercc_raw_gene_bc_matrices.tar.gz")
tempdir <- tempfile()
dir.create(tempdir)
untar(fname, exdir=tempdir)

# Reading in the libraries.
library(DropletUtils)
sce <- read10xCounts(file.path(tempdir, "matrices_mex/ercc92"))

##############################################

library(Matrix)
counts <- counts(sce)
counts <- counts[,colSums(counts) > 100] # retaining some empty droplets for testing.
counts <- counts[rowMeans(counts) > 0.01,]

MGC <- run_MAGIC(as.matrix(counts), t=10)
lcounts <- as.matrix(lognormalize(counts))

##############################################

# Performing t-SNE on the result, compared to the original.
set.seed(2000)
library(Rtsne)
mgc_TSNE <- Rtsne(MGC, pca=FALSE, perplexity=30)
normal_TSNE <- Rtsne(lcounts, pca=FALSE, perplexity=30)

# Creating plots.
pdf("pics/ercc_results.pdf")
col = viridis::viridis(100)[cut(log(colSums(counts)), 100)]
plot(mgc_TSNE$Y[,1], mgc_TSNE$Y[,2], xlab="TSNE1", ylab="TSNE2", pch=16, col=col)
plot(normal_TSNE$Y[,1], normal_TSNE$Y[,2], xlab="TSNE1", ylab="TSNE2", pch=16, col=col)
dev.off()
