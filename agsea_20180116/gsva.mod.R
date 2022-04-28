gsvacompetitive <- function(y,  design, covar, lgroups, index, nperm = 1000, mccores = 1, comp = FALSE, sizeId = 50, ...){

    y <- t(scale(t(y)))
    rownames(y) <- 1:nrow(y)
    gsvaD <- gsva(y, index, ...)
    lmf <- lmFit(gsvaD, design)
    cont.mat <- contMatFunction(covar, form = form, lgroups = lgroups)

    lmc <- contrasts.fit(lmf, cont.mat)
    eb <- eBayes(lmc)
    dsgs <- topTable(eb, adjust="BH", number=nrow(gsvaD), coef = ncol(eb), confint=T)
    dsgs <- dsgs[order(-abs(dsgs$t)),]

    dsgs$ngenes <- sapply(index[rownames(dsgs)],length)
    dsgs <- dsgs[,c("ngenes","logFC","P.Value","adj.P.Val","t")]
    names(dsgs) <- c("ngenes","est","pval", "adj.pval","t")

    if(comp){
     sam <- mclapply(1:nperm, function(j){ sample(rownames(y),sizeId)}, mc.cores = mccores)
     gsvaR <- gsva(y, sam, ...)
     lmfr <- lmFit(gsvaR, design)
     lmcr <- contrasts.fit(lmfr, cont.mat)
     ebr <- eBayes(lmc)
     dsgs.rand <- topTable(ebr, adjust="BH", number=nrow(gsvaR), coef = ncol(ebr), confint=T)
     dsgs.rand <- as.data.frame(dsgs.rand)

     pval.comp <- unlist(mclapply(dsgs$t, function(o){ sum(o <= dsgs.rand$t)}, mc.cores = mccores))
     pval.comp <- 2 * pmin(pval.comp, nperm - pval.comp)
     dsgs$pval.comp <- (pval.comp + 1) /(nperm + 1)
     dsgs$apval.comp <- p.adjust(dsgs$pval.comp,method="BH")
 }

    return(list(res = dsgs, stats = gsvaD))

}


