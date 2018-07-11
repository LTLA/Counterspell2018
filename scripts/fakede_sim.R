source("functions.R")

# Scenario 2: DE between two clusters

set.seed(400)
ngenes = 900
ncells = 500
lambda = 2^rexp(ngenes)
counts = rbind(
    matrix(rpois(ncells * ngenes, lambda=lambda), ncol=ncells),
    matrix(rpois(ncells * 50, lambda=rep(c(0, 5), each=ncells/2)), ncol=ncells, byrow=TRUE),
    matrix(rpois(ncells * 50, lambda=rep(c(5, 0), each=ncells/2)), ncol=ncells, byrow=TRUE)
)

cluster = rep(1:2, each=ncells/2)
cluster_cols = c(`1`="forestgreen", `2`="goldenrod")
is_de = c(logical(ngenes), !logical(nrow(counts)-ngenes))
de_cols = c(`TRUE`="lavender", `FALSE`="violet")    

# Creating a PCA plot on the log-normalized counts from MAGIC.
library(Rmagic) 
MGC = run_magic(t(counts), t_diffusion=10)

png("pics/clusters_with_magic.png", units="in", width=7, height=7, res=150, pointsize=12)
lmgc = lognormalize(t(MGC))
create_heatmap(t(lmgc), cluster=cluster, is_de=is_de, limit=0.2, cluster_cols=cluster_cols, de_cols=de_cols)
dev.off()

# Creating a PCA plot on the reference.
png("pics/clusters_without_magic.png", units="in", width=7, height=7, res=150, pointsize=12)
lref = lognormalize(counts)
create_heatmap(t(lref), cluster=cluster, is_de=is_de, limit=0.2, cluster_cols=cluster_cols, de_cols=de_cols)
dev.off()
