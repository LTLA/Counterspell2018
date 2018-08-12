source("functions.R")

# Scenario 2: DE between two clusters

set.seed(400)
ngenes <- 9900
ncells <- 500
lambda <- runif(ngenes, 0, 10)
counts <- rbind(
    matrix(rpois(ncells * ngenes, lambda=lambda), ncol=ncells),
    matrix(rpois(ncells * 50, lambda=rep(c(0, 5), each=ncells/2)), ncol=ncells, byrow=TRUE),
    matrix(rpois(ncells * 50, lambda=rep(c(5, 0), each=ncells/2)), ncol=ncells, byrow=TRUE)
)

cluster <- rep(1:2, each=ncells/2)
cluster_cols <- c(`1`="forestgreen", `2`="goldenrod")
is_de <- c(logical(ngenes), !logical(nrow(counts)-ngenes))
de_cols <- c(`TRUE`="dodgerblue", `FALSE`="grey60")    

# Running MAGIC.
MGC <- run_MAGIC(counts, t=10)
lref <- lognormalize(counts)

# Applying a series of Wilcoxon rank-sum tests.

p.mgc <- numeric(nrow(counts))
for (i in seq_along(p.mgc)) {
    p.mgc[i] <- wilcox.test(MGC[cluster==1,i][,1], MGC[cluster==2,i][,1])$p.value
}

p.log <- numeric(nrow(counts))
for (i in seq_along(p.log)) {
    p.log[i] <- wilcox.test(lref[cluster==1,i], lref[cluster==2,i])$p.value
}

pdf("pics/fakede.pdf")
breaks <- seq(0, 1, length.out=50)
hist(p.mgc[!is_de], breaks=breaks, col="violet")
hist(p.log[!is_de], breaks=breaks, col="lavender")
dev.off()
