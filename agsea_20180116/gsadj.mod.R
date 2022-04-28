batchcorlm <- function(y, design, contrast){
    lmf1 <- lmFit(y, design[,-contrast])
    y <- y - lmf1$coef[,-contrast]%*%t(design[,-contrast])
    return(y)
}


gsadjcompetitive <- function(y, design, contrast, index, nperm = 1000, mccores = 1, comp = FALSE, sizeId = 50, addGS =TRUE, geneCorBatch = FALSE){
    if(geneCorBatch)
        y <- batchcorlm(y, design, contrast)

    y <- t(scale(t(y)))
    GS <- as.numeric( (rep(1,dim(y)[1]) %*% y)/dim(y)[1])
    gsvaD <- t(sapply(index, function(o) (as.numeric(rep(1,length(o))%*%y[o,])/length(o))))
    colnames(gsvaD) <- colnames(y)
    if(addGS){
        design <- cbind(GS,design)
        contrast <- contrast + 1
    }

#    lmcr <- contrasts.fit(lmf2, contrast)
#     ebr <- eBayes(lmc)
    dsgs <- t(apply(gsvaD, 1, function(x){
        lm1 <- summary(lm(x ~ design))$coeff
        as.numeric(lm1[contrast,])
    }))
    if(addGS){
        lmf2 <- lmFit(gsvaD,design)
        gsvaD <- gsvaD - lmf2$coef[,1]%*%t(design[,1])
    }
    colnames(dsgs) <- c("est","sd","t","pval")
    dsgs <- as.data.frame(dsgs)
    dsgs <- dsgs[order(-abs(dsgs$t)),]
    dsgs$adj.pval <- p.adjust(dsgs$pval,method="BH")

    dsgs$ngenes <- sapply(index[rownames(dsgs)],length)
    dsgs <- dsgs[,c("ngenes","est","pval","adj.pval","t")]

    if(comp){
     sam <- mclapply(1:nperm, function(j){ sample(rownames(y),sizeId)}, mc.cores = mccores)
     gsva.rand <- mclapply(sam, function(o){ (as.numeric(rep(1,length(o))%*%y[o,])/length(o))},mc.cores = mccores)

     dsgs.rand <- do.call(rbind, mclapply(gsva.rand, function(x){
        lm1 <- summary(lm(x ~ design))$coeff
        as.numeric(lm1[contrast,])},mc.cores = mccores))
     colnames(dsgs.rand) <- c("est","sd","t","pval")

     dsgs.rand <- as.data.frame(dsgs.rand)
     pval.comp <- unlist(mclapply(dsgs$t, function(o){ sum(o <= dsgs.rand$t)}, mc.cores = mccores))
     pval.comp <- 2 * pmin(pval.comp, nperm - pval.comp)
     dsgs$pval.comp <- (pval.comp + 1) /(nperm + 1)
     dsgs$apval.comp <- p.adjust(dsgs$pval.comp,method="BH")
    }
    return(list(res = dsgs, stats = gsvaD))
}





gsadjrepresentation <- function(y, design, geneset, testvar = colnames(design)[ncol(design)], PLOT= FALSE,
                                addGS = TRUE,  coef2cor = NULL, alpha = 0.90, maintitle = "", subtitle = "",  ...){
    y <- t(scale(t(y)))
    GS <- as.numeric( (rep(1,dim(y)[1]) %*% y)/dim(y)[1])
    geneset <- intersect(rownames(y),geneset)
    gsvaD <- (as.numeric(rep(1,length(geneset))%*%y[geneset,])/length(geneset))

    if(addGS){
        design <- cbind(GS,design)
        lmf2 <- lm(gsvaD ~ -1 + design)
        coef2cor <- c(coef2cor, "GS")
    }
    else gsvaD2 <- gsvaD

    lmf2 <- lm(gsvaD ~ -1 + design)
    gsvaD2 <- as.numeric(gsvaD - lmf2$coef[-match(coef2cor,colnames(design))]%*%t(design[,-match(coef2cor,colnames(design))]))
    pr <- predict(lm(gsvaD2 ~ -1 + design[,match(coef2cor,colnames(design))]), se.fit=TRUE)

    Qt <- c(-1, 1) * qt((1 - alpha) / 2, lmf2$df, lower.tail = FALSE)

    ## 90% confidence interval
    puntpred <- as.numeric(gsvaD - lmf2$coef[paste0("design",coef2cor)]%*%t(design[,coef2cor]))
    CI <-  puntpred + outer(pr$se.fit, Qt)
    colnames(CI) <- c("lwr", "upr")
    if(PLOT){
        if (testvar[1] != ""){
            plot(design[,testvar], CI[,1],  xlab =  testvar, ylab = paste0("Predicted expression geneset"), pch = "-", ylim = c(min(CI),max(CI)), main = maintitle, ...)
            mtext(subtitle, 3, line=0.5, col= "darkgray")
            points(design[,testvar], CI[,2],  xlab =  testvar, ylab = paste0("Predicted expression geneset"), pch = "-")
            for(o in 1:nrow(CI)) lines(rep(design[o,testvar],2), CI[o,], col= "grey")
            text(design[,testvar], puntpred, colnames(y))
        }
        else{
            plot(CI[,1],  xlab =  "", xaxt = "n", ylab = paste0("Predicted expression geneset"), pch = "-", ylim = c(min(CI),max(CI)))
            points(CI[,2],  xlab =  testvar, ylab =  paste0("Predicted expression geneset"), pch =  "-")
            for(o in 1:nrow(CI)) lines(rep(o,2), CI[o,], col= "grey")
            text(1:length(puntpred), puntpred, colnames(y))
        }
    }
    return(list(puntpred= puntpred, CI = CI))
}

#    plot(CI[,1],  xlab =  testvar, ylab = paste0("expression geneset ", genesetname), pch = "-", ylim = c(min(CI),max(CI)))
#    points(CI[,2],  xlab =  testvar, ylab = paste0("expression geneset ", genesetname), pch =  "-")
#    for(o in 1:nrow(CI)) lines(rep(o,2), CI[o,], col= "grey")
#    text(1:length(gsvaD2), as.numeric(gsvaD - lmf2$coef[1]%*%t(design[,1])), colnames(y))




