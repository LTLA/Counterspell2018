source("functions.R")

# Using the public Jurkat data at:
# http://cf.10xgenomics.com/samples/cell-exp/1.1.0/jurkat/jurkat_filtered_gene_bc_matrices.tar.gz

library(Matrix)
counts = readMM("jurkat/hg19/matrix.mtx")
genes = read.table("jurkat/hg19/genes.tsv")
rownames(counts) = genes[,1]
keep = rowMeans(counts) > 0.1
counts = counts[keep,]

library(Rmagic)
MGC = run_magic(as.matrix(t(counts)), t_diffusion=10)

# Performing PCA on the result, compared to the original.
set.seed(1000)    
library(irlba)
lmgc = lognormalize(t(MGC))
mgc_PCA = prcomp_irlba(lmgc, n=20)
total_var_mgc = sum(apply(lmgc, 2, var))

lcounts = lognormalize(as.matrix(counts))
normal_PCA = prcomp_irlba(lcounts, n=20)
total_var_normal = sum(apply(lcounts, 2, var))

# Creating plots.
pdf("pics/jurkat_results.pdf")
col = viridis::viridis(100)[cut(mgc_PCA$x[,1], 100)]

plot(mgc_PCA$x[,1], mgc_PCA$x[,2], 
    xlab=sprintf("PC1 (%.1f%%)", mgc_PCA$sdev[1]^2/total_var_mgc * 100),
    ylab=sprintf("PC2 (%.1f%%)", mgc_PCA$sdev[2]^2/total_var_mgc * 100),
    pch=16, col=col)

plot(normal_PCA$x[,1], normal_PCA$x[,2], 
    xlab=sprintf("PC1 (%.1f%%)", normal_PCA$sdev[1]^2/total_var_normal * 100),
    ylab=sprintf("PC2 (%.1f%%)", normal_PCA$sdev[2]^2/total_var_normal * 100),
    pch=16, col=col)
dev.off()
