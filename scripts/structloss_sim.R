source("functions.R")

# Exploring suppression of minor structure.

set.seed(600)
ngenes = 800
ncells = 500
lambda = 2^rexp(ngenes)
counts = rbind(
    matrix(rpois(ncells * ngenes, lambda=lambda), ncol=ncells),
    matrix(rpois(ncells * 50, lambda=rep(c(0, 5), ncells/2)), ncol=ncells, byrow=TRUE),
    matrix(rpois(ncells * 50, lambda=rep(c(5, 0), ncells/2)), ncol=ncells, byrow=TRUE),
    matrix(rpois(ncells * 50, lambda=rep(c(0, 2), each=ncells/2)), ncol=ncells, byrow=TRUE),
    matrix(rpois(ncells * 50, lambda=rep(c(2, 0), each=ncells/2)), ncol=ncells, byrow=TRUE)
)

color = rep(c("forestgreen", "goldenrod"), each=ncells/2)
pch = rep(c(1,4), ncells/2)

# Creating a PCA plot on the log-normalized counts from MAGIC.
library(Rmagic) 
MGC = run_magic(t(counts), t_diffusion=10)

pdf("pics/minor_with_magic.pdf")
lmgc = lognormalize(t(MGC))
run_PCA(lmgc, xlim=c(-15, 15), ylim=c(-15, 15), col=color, pch=pch)
dev.off()

# Creating a PCA plot on the reference.
pdf("pics/minor_without_magic.pdf")
lref = lognormalize(counts)
run_PCA(lref, xlim=c(-15, 15), ylim=c(-15, 15), col=color, pch=pch)
dev.off()
