fgsea.mod <- function(index, pr.stat, nperm, minsize = 10, maxsize = 300, mccores = 1, ...){
    res <- fgsea(index, stats = pr.stat, nperm = nperm, minSize = minsize, maxSize = maxsize, nproc = mccores, ...)
    res <- as.data.frame(res)
    rownames(res) <- res$pathway
    res <- res[,c(7,4,2,3,5,8)]
    names(res) <- c("ngenes", "est", "pval","adj.pvalue","nes", "leadingEdge")
    res <- res[order(res$pval),]
    res <- list(res = res, stats = pr.stat)
    class(res) <- c("fgsea.mod")
    return(res)
}

print.fgsea.mod <- function(x)  print(as.data.frame(x$res)[,-ncol(x$res)])
