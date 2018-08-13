# This script tests the behaviour of MAGIC on the public 293T data from 10X Genomics.

source("functions.R")

##############################################

library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
fname <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/293t_filtered_gene_bc_matrices.tar.gz")
tempdir <- tempfile()
dir.create(tempdir)
untar(fname, exdir=tempdir)

# Reading in the libraries.
library(DropletUtils)
sce <- read10xCounts(file.path(tempdir, "filtered_matrices_mex/hg19"))

##############################################

library(Matrix)
counts <- counts(sce)
counts <- counts[rowMeans(counts) > 0.01,]

MGC <- run_MAGIC(counts)
lcounts <- lognormalize(counts)

##############################################

# Performing PCA on the result, compared to the original.
set.seed(1000)    
library(irlba)
mgc_PCA <- prcomp_irlba(MGC, n=20)
total_var_mgc <- sum(apply(MGC, 2, var))

normal_PCA <- prcomp_irlba(lcounts, n=20)
total_var_normal <- sum(apply(lcounts, 2, var))

##############################################

# Creating plots.
pdf("pics/293t_results.pdf")
par(cex.lab=1.4)
discolored = mgc_PCA$x[,2] < 0 & mgc_PCA$x[,1] < -5
col = ifelse(discolored, "forestgreen", "goldenrod")

plot(mgc_PCA$x[,1], mgc_PCA$x[,2], 
    xlab=sprintf("PC1 (%.1f%%)", mgc_PCA$sdev[1]^2/total_var_mgc * 100),
    ylab=sprintf("PC2 (%.1f%%)", mgc_PCA$sdev[2]^2/total_var_mgc * 100),
    pch=16, col=col)

plot(normal_PCA$x[,1], normal_PCA$x[,2], 
    xlab=sprintf("PC1 (%.1f%%)", normal_PCA$sdev[1]^2/total_var_normal * 100),
    ylab=sprintf("PC2 (%.1f%%)", normal_PCA$sdev[2]^2/total_var_normal * 100),
    pch=16, col=col)
dev.off()
