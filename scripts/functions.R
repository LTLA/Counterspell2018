#' Runs PCA on a matrix of log-counts.
run_PCA = function(log_counts, plot=TRUE, pch=16, ..., xlim=NULL, ylim=NULL, same_xy=FALSE) {
    PCA = prcomp(log_counts)
    var_exp = PCA$sdev^2 

    if (plot) { 
        X = PCA$x[,1]
        Y = PCA$x[,2]

        if (same_xy) {
            xlim = range(X)
            ylim = range(Y)
            xlim = ylim = range(c(xlim, ylim))
        } 

        plot(PCA$x[,1], PCA$x[,2], ..., pch=pch, xlim=xlim, ylim=ylim,
	        xlab=sprintf("PC1 (%.1f%%)", var_exp[1]/sum(var_exp) * 100),
        	ylab=sprintf("PC2 (%.1f%%)", var_exp[2]/sum(var_exp) * 100))
    }
    return(invisible(var_exp))
}

#' Computes log-transformed expression values; returns a transposed matrix.
lognormalize = function(counts) {
    lib_size = colSums(counts)
    lib_size = lib_size/mean(lib_size)
    log2(t(counts)/lib_size + 1)
}

#' Identifies and plots the most correlated gene pair.
plot_correlated = function(log_counts, ...) {
    COR = cor(log_counts)
    diag(COR) = 0
    max_COR = max(COR, na.rm=TRUE)
    i = which(COR == max_COR, arr.ind=TRUE)[1,]
    i1 = i[1]
    i2 = i[2]

    X = log_counts[,i1]
    Y = log_counts[,i2]
    plot(X, Y, ...,  
        xlab=sprintf("Gene %i", i1),  
        ylab=sprintf("Gene %i", i2))

    fit <- lm(Y ~ X)
    abline(a=coef(fit)[1], b=coef(fit)[2], col="red")
    legend("bottomright", bty="n", legend=substitute(R^2~"="~r2, 
        list(r2=round(max_COR^2, 2))), text.col="red")
}

#' Calculates log-fold changes for all genes.
compute_logFC = function(log_counts, groupings) {
    all_groups = unique(groupings)
    colMeans(log_counts[groupings==all_groups[1],]) - colMeans(log_counts[groupings==all_groups[2],])
}

#' Creates a heatmap with mean-centering and limits.
create_heatmap = function(exprs, cluster, is_de, limit=0.5, cluster_cols, de_cols) {
    uniq_clusters = unique(cluster)
    reorder = order(rowMeans(exprs[,cluster==uniq_clusters[1]]) - rowMeans(exprs[,cluster==uniq_clusters[2]]))
    exprs = exprs[reorder,]
    is_de = is_de[reorder]

    exprs = exprs - rowMeans(exprs)
    exprs[exprs > limit] = limit
    exprs[exprs < -limit] = -limit

    colnames(exprs) = seq_len(ncol(exprs))
    rownames(exprs) = seq_len(nrow(exprs))

    require(pheatmap)
    pheatmap(exprs, cluster_col=FALSE, cluster_row=FALSE,
        annotation_col=data.frame(Cluster=as.character(cluster), row.names=colnames(exprs)),
        annotation_row=data.frame(DE=as.character(is_de), row.names=rownames(exprs)),
        annotation_colors=list(Cluster=cluster_cols, DE=de_cols),
        show_rownames=FALSE, show_colnames=FALSE)
}

#' Also creating a directory.
dir.create("pics", showWarnings=FALSE)
