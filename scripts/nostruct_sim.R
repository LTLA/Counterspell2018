source("functions.R")

# Generating Poisson-distributed counts.
set.seed(200)
ngenes <- 1000
ncells <- 500
lambda <- runif(ngenes, 0, 10)
size_fac <- runif(ncells, 1, 2)
counts <- matrix(rpois(ncells * ngenes, lambda=outer(lambda, size_fac)), ncol=ncells)

library(viridis)
col <- viridis(100)[cut(size_fac, 100)]

MGC <- run_MAGIC(counts, t=10)
lref <- lognormalize(counts)

# Creating a PCA plot on the log-normalized counts from MAGIC.
pdf("pics/nostructure.pdf")
par(cex.lab=1.4)
run_PCA(MGC, same_xy=TRUE, col=col)
plot_correlated(MGC)

run_PCA(lref, same_xy=TRUE, col=col)
plot_correlated(lref, position="topleft")
dev.off()
