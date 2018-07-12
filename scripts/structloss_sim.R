source("functions.R")

# Exploring suppression of minor structure.

set.seed(600)
ngenes = 800
ncells = 500
lambda = runif(ngenes, 0, 10)
counts = rbind(
    matrix(rpois(ncells * ngenes, lambda=lambda), ncol=ncells),
    matrix(rpois(ncells * 50, lambda=rep(c(1, 5), ncells/2)), ncol=ncells, byrow=TRUE),
    matrix(rpois(ncells * 50, lambda=rep(c(5, 1), ncells/2)), ncol=ncells, byrow=TRUE),
    matrix(rpois(ncells * 50, lambda=rep(c(1, 2), each=ncells/2)), ncol=ncells, byrow=TRUE),
    matrix(rpois(ncells * 50, lambda=rep(c(2, 1), each=ncells/2)), ncol=ncells, byrow=TRUE)
)

primary_sep = rep(c(TRUE, FALSE), ncells/2)
pch = ifelse(primary_sep, 1, 4)
secondary_sep = rep(c(TRUE, FALSE), each=ncells/2)
color = ifelse(secondary_sep, "forestgreen", "goldenrod")
de_cols = rep(c("grey60", "dodgerblue"), c(ngenes, nrow(counts) - ngenes))

# Creating a PCA plot on the log-normalized counts from MAGIC.
library(Rmagic)
MGC = run_MAGIC(counts, t=10)

pdf("pics/minor_with_magic.pdf")
run_PCA(MGC, same_xy=TRUE, col=color, pch=pch)

plot(compute_logFC(MGC, primary_sep),
    compute_logFC(MGC, secondary_sep),
    xlab="Difference (primary)",
    ylab="Difference (secondary)",
    xlim=c(-2, 2), ylim=c(-2, 2),
    col=de_cols, pch=16)
dev.off()

# Creating a PCA plot on the reference.
pdf("pics/minor_without_magic.pdf")
lref = lognormalize(counts)
run_PCA(lref, same_xy=TRUE, col=color, pch=pch)

plot(compute_logFC(lref, primary_sep),
    compute_logFC(lref, secondary_sep),
    xlab=expression(Log[2]~"fold change (primary)"),
    ylab=expression(Log[2]~"fold change (secondary)"),
    xlim=c(-2, 2), ylim=c(-2, 2),
    col=de_cols, pch=16)
dev.off()
