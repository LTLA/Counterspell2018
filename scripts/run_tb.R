# This tests for induced correlations between irrelevant genes upon using MAGIC.

source("functions.R")

#########################################

library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)

cd4.fname <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd4_t_helper/cd4_t_helper_filtered_gene_bc_matrices.tar.gz")
cd4.tempdir <- tempfile()
dir.create(cd4.tempdir)
untar(cd4.fname, exdir=cd4.tempdir)

b.fname <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz")
b.tempdir <- tempfile()
dir.create(b.tempdir)
untar(b.fname, exdir=b.tempdir)

#########################################

# Reading in the libraries.
library(DropletUtils)
sce.cd4 <- read10xCounts(file.path(cd4.tempdir, "filtered_matrices_mex/hg19"))
sce.b <- read10xCounts(file.path(b.tempdir, "filtered_matrices_mex/hg19"))

library(Matrix)
SPAWNER <- function(nchosen) { 
    sub.b <- sce.b[,sample(ncol(sce.b), nchosen)]
    sub.cd4 <- sce.cd4[,sample(ncol(sce.cd4), nchosen)]
    sce <- cbind(sub.b, sub.cd4)
    counts <- counts(sce)
    counts <- counts[rowMeans(counts) > 0.01,]
    counts
}

#########################################

# Running multiple times to ensure results are representative.
set.seed(10000)
pdf("pics/tb_results.pdf")
for (it in 1:10) { # LOOP START

#########################################

library(Rmagic)

nlow <- 200
counts.low <- SPAWNER(nlow)
MGC.low <- run_MAGIC(counts.low)

nhigh <- 5000
counts.high <- SPAWNER(nhigh)
MGC.high <- run_MAGIC(counts.high)

#########################################

gene1 <- "ENSG00000177455" # CD19
gene2 <- "ENSG00000010610" # CD4
symb1 <- "CD19"
symb2 <- "CD4"

cols <- c("violet", "dodgerblue")
pchs <- c(16, 4)

par(cex.lab=1.4, cex.main=1.5)
plot(MGC.low[,gene1], MGC.low[,gene2], xlab=symb1, ylab=symb2, 
    main=paste(nlow, "cells each"), 
	col=rep(cols, each=nlow), 
	pch=rep(pchs, each=nlow)) 
legend("topright", pch=pchs, col=cols, legend=c("B cell", "CD4+ T cell"), cex=1.4)
plot(MGC.high[,gene1], MGC.high[,gene2], xlab=symb1, ylab=symb2, 
    main=paste(nhigh, "cells each"), 
    col=rep(cols, each=nhigh), 
    pch=rep(pchs, each=nhigh))

#########################################

gc()
} # LOOP END
dev.off()
