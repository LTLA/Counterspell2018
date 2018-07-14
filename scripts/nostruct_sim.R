source("functions.R")

# Generating Poisson-distributed counts.
set.seed(200)
ngenes = 1000
ncells = 500
lambda = runif(ngenes, 0, 10)
size_fac = runif(ncells, 1, 2)
counts = matrix(rpois(ncells * ngenes, lambda=outer(lambda, size_fac)), ncol=ncells)

library(viridis)
col = viridis(100)[cut(size_fac, 100)]

# Creating a PCA plot on the log-normalized counts from MAGIC.
MGC = run_MAGIC(counts, t=10)

pdf("pics/nostructure_with_magic.pdf")
run_PCA(MGC, same_xy=TRUE, col=col)
plot_correlated(MGC)
dev.off()

# Creating a PCA plot on the reference.
pdf("pics/nostructure_without_magic.pdf")
lref = lognormalize(counts)
run_PCA(lref, same_xy=TRUE, col=col)
plot_correlated(lref, position="topleft")
dev.off()
