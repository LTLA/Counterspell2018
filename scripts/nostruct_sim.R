source("functions.R")

# Generating Poisson-distributed counts.
set.seed(1000)
ngenes = 1000
ncells = 500
lambda = runif(ngenes, 0, 10)
counts = matrix(rpois(ncells * ngenes, lambda=lambda), ncol=ncells)

# Creating a PCA plot on the log-normalized counts from MAGIC.
library(Rmagic) 
MGC = run_magic(t(counts), t_diffusion=10)

pdf("pics/nostructure_with_magic.pdf")
lmgc = lognormalize(t(MGC))
run_PCA(lmgc, same_xy=TRUE)
plot_correlated(lmgc)
dev.off()

# Creating a PCA plot on the reference.
pdf("pics/nostructure_without_magic.pdf")
lref = lognormalize(counts)
run_PCA(lref, same_xy=TRUE)
plot_correlated(lref)
dev.off()

# Varying the diffusion time.
output = numeric(20)
for (i in seq_along(output)) {
    MGC = run_magic(t(counts), t_diffusion=i)
    lmgc = log2(MGC/rowSums(MGC) + 1)
    var_exp = run_PCA(lmgc, plot=FALSE)
    output[i] = var_exp[1]/sum(var_exp) * 100
}

pdf("pics/nostructure_time.pdf", width=10, height=5)
plot(output, xlab="Diffusion time", ylab="Variance explained by PC1 (%)", pch=16)
lines(output)
dev.off()
