
enrichmentScore <- function(fchr,sign) {
  nfchr <- length(fchr)
  p <- numeric(nfchr)
  p[] <- -1 / (nfchr - length(sign))
  absfchr <- abs(fchr[sign])
  p[sign] <- absfchr / sum(absfchr)
  cumsum(p)
}

enrichmentScore.mat <- function(fchr, indexmat) {
    nfchr <- length(fchr)
    nsign <- Matrix::colSums(indexmat)
    p <- array( rep(-1 / (nfchr - nsign),each=nfchr), dim=c(nfchr,length(sign)))
    absfchr <- abs(fchr * indexmat)
    p <- t(t(absfchr)/Matrix::colSums(absfchr)) + (1-indexmat) * p
    apply(p,2,cumsum)
}



romer.mod <- function (y, index, design = NULL, contrast = ncol(design), array.weights = NULL,
                       block = NULL, correlation = NULL, set.statistic = "mean", psel = NULL,
                       nrot = 9999, shrink.resid = TRUE, restand = TRUE, mccores = 1, executation.info = TRUE, ...)
{
   dots <- names(list(...))
   if (length(dots))
    warning("Extra arguments disregarded: ", sQuote(dots))
   y <- as.matrix(y)
   ngenes <- nrow(y)
   n <- ncol(y)
   if (!is.list(index))
    index <- list(set = index)
   nsets <- length(index)
   if (nsets == 0)
    stop("index is empty")
   SetSizes <- unlist(lapply(index, length))
   if (is.null(design))
   {
    stop("design matrix not specified")
   }else
   {
    design <- as.matrix(design)
    if (mode(design) != "numeric")
      stop("design must be a numeric matrix")
   }
   if (nrow(design) != n)
    stop("row dimension of design matrix must match column dimension of data")
   ne <- nonEstimable(design)
   if (!is.null(ne))
    cat("Coefficients not estimable:", paste(ne, collapse = " "),
        "\n")
   p <- ncol(design)
   if (p < 2L)
    stop("design needs at least two columns")
   p0 <- p - 1L
   d <- n - p
   if (length(contrast) == 1L)
   {
    contrast <- round(contrast)
    if (contrast < p)
    {
      i <- 1L:p
      design <- design[, c(i[-contrast], contrast)]
    }
   }else
   {
    design <- contrastAsCoef(design = design, contrast = contrast,
                             first = FALSE)$design
   }

   if (!is.null(array.weights))
   {
     if (any(array.weights <= 0))
      stop("array.weights must be positive")
     if (length(array.weights) != n)
      stop("Length of array.weights doesn't match number of array")
     design <- design * sqrt(array.weights)
     y <- t(t(y) * sqrt(array.weights))
   }
   if (!is.null(block))
   {
     block <- as.vector(block)
     if (length(block) != n)
      stop("Length of block does not match number of arrays")
     ub <- unique(block)
     nblocks <- length(ub)
     Z <- matrix(block, n, nblocks) == matrix(ub, n, nblocks,
                 byrow = TRUE)
     cormatrix <- Z %*% (correlation * t(Z))
     diag(cormatrix) <- 1
     R <- chol(cormatrix)
     y <- t(backsolve(R, t(y), transpose = TRUE))
      design <- backsolve(R, design, transpose = TRUE)
   }
   qr <- qr(design)
   signc <- sign(qr$qr[p, p])
   effects <- qr.qty(qr, t(y))
   s2 <- colMeans(effects[-(1:p), , drop = FALSE]^2)
   sv <- squeezeVar(s2, df = d)
   d0 <- sv$df.prior
   s02 <- sv$var.prior
   sd.post <- sqrt(sv$var.post)
   Y <- effects[-(1:p0), , drop = FALSE]
   YY <- colSums(Y^2)
   B <- Y[1, ]
   modt <- signc * B/sd.post



  if (shrink.resid)
  {
    p.value <- 2 * pt(abs(modt), df = d0 + d, lower.tail = FALSE)
    proportion <- 1 - propTrueNull(p.value)
    stdev.unscaled <- rep_len(1/abs(qr$qr[qr$rank, qr$rank]), ngenes)
    var.unscaled <- stdev.unscaled^2
    df.total <- rep_len(d, ngenes) + sv$df.prior
    stdev.coef.lim <- c(0.1, 4)
    var.prior.lim <- stdev.coef.lim^2/sv$var.prior
    var.prior <- tmixture.vector(modt, stdev.unscaled, df.total,
                                 proportion, var.prior.lim)
    if (any(is.na(var.prior)))
    {
      var.prior[is.na(var.prior)] <- 1/sv$var.prior
      warning("Estimation of var.prior failed - set to default value")
    }
    r <- (var.unscaled + var.prior)/var.unscaled
    if (sv$df.prior > 10^6)
    {
      kernel <- modt^2 * (1 - 1/r)/2
    }else
    {
      kernel <- (1 + df.total)/2 * log((modt^2 + df.total)/(modt^2/r + df.total))
    }
      lods <- log(proportion/(1 - proportion)) - log(r)/2 + kernel
    ProbDE <- exp(lods)/(1 + exp(lods))
    ProbDE[is.na(ProbDE)] <- 1
    Y[1, ] <- Y[1, ] * sqrt(var.unscaled/(var.unscaled + var.prior * ProbDE))
  }

   if(!is.null(psel)){
      modt <- modt[psel]
      sd.post <- sd.post[psel]
      YY <- YY[psel]
      Y <- Y[,psel]
      s2 <- s2[psel]

      names(sd.post) <- colnames(Y) <- names(s2) <- names(modt) <- names(YY) <- names(psel)
      ngenes <- length(modt)
   }

  set.statistic <- match.arg(set.statistic, c("mean", "mean.rank", "absmean", "median", "maxmean", "GSEA", "GSVA"))
  if(executation.info)   pb <- txtProgressBar(min = 0, max = nrot, style = 3)

  if (set.statistic == "maxmean")
  {
    estpos <- sapply(index, function(o, modt) sum((modt[o] + abs(modt[o]))/2)/length(o), modt)
    estneg <- sapply(index, function(o, modt) sum((abs(modt[o]) - modt[o])/2)/length(o), modt)
    nuapos <- mean(modt*(modt > 0))
    nuaneg <- mean(-modt*(modt < 0))
    sdapos <- sd(modt*(modt > 0))
    sdaneg <- sd(modt*(modt < 0))

    #print(estpos)
    #print(estneg)
    
    
    if (restand)
    {
      estpos <- (estpos - nuapos)/sdapos
      estneg <- (estneg - nuaneg)/sdaneg
    }
    est <- pmax(estpos, estneg)
    #print(est)
    #print(estpos)
    #print(estneg)
    est[estneg > estpos] = -1 * est[estneg > estpos]
    p.value <- matrix(0, nrow = nsets, ncol = 2)
    est0.all <- do.call(cbind,mclapply(seq_len(nrot),function(i)
    {
      if(executation.info ) setTxtProgressBar(pb, i)
      R <- matrix(rnorm((d + 1)), 1, d + 1)
      R <- R/sqrt(rowSums(R^2))
      Br <- R %*% Y
      s2r <- (YY - Br^2)/d
      if (is.finite(d0))
      {
        sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
      }else
      {
        sdr.post <- sqrt(s02)
      }
      modtr <- signc * Br/sdr.post
      estpos0 <- sapply(index, function(o, modtr) sum((modtr[o] + abs(modtr[o]))/2)/length(o), modtr)
      estneg0 <- sapply(index, function(o, modtr) sum((abs(modtr[o]) - modtr[o])/2)/length(o), modtr)
      nuapos0 <- mean(modtr*(modtr > 0))
      nuaneg0 <- mean(-modtr*(modtr < 0))
      sdapos0 <- sd(modtr*(modtr > 0))
      sdaneg0 <- sd(modtr*(modtr < 0))
      if (restand)
      {
        estpos0 <- (estpos0 - nuapos0)/sdapos0
        estneg0 <- (estneg0 - nuaneg0)/sdaneg0
      }
      est0 <- pmax(estpos0, estneg0)
      est0[estneg0 > estpos0] = -1 * est0[estneg0 > estpos0]
     # p.value <- est0 <= est
      est0
    }, mc.cores = mccores))
   p.value <- sapply(seq_len(dim(est0.all)[1]), function(o) sum(est[o]>=est0.all[o,]))
   pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
   nes <- est/ apply(est0.all,1,sd)
  }

  if (set.statistic == "mean")
  {
    est <- sapply(index, function(o, modt) mean(modt[o]), modt)
    nua <- mean(modt)
    sda <- sd(modt)
    if (restand)
    {
      est <- (est - nua)/sda
    }
    p.value <- matrix(0, nrow = nsets, ncol = 2)
    est0.all <- do.call(cbind,mclapply(seq_len(nrot),function(i)
    {
      if(executation.info ) setTxtProgressBar(pb, i)
      R <- matrix(rnorm((d + 1)), 1, d + 1)
      R <- R/sqrt(rowSums(R^2))
      Br <- R %*% Y
      which(is.na(Br))
      s2r <- (YY - Br^2)/d
      if (is.finite(d0))
      {
        sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
      }else
      {
        sdr.post <- sqrt(s02)
      }
      modtr <- signc * Br/sdr.post
      est0 <- sapply(index, function(o, modtr) mean(modtr[o]), modtr)
      nua0 <- mean(modtr)
      sda0 <- sd(modtr)
      if (restand)
      {
        est0 <- (est0 - nua0)/sda0
      }
      est0 #p.value <- est0 <= est #
  }, mc.cores = mccores))
   p.value <- sapply(seq_len(dim(est0.all)[1]), function(o) sum(est[o]>=est0.all[o,]))
   pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
   nes <- est/ apply(est0.all,1,sd)
  }

  if (set.statistic == "mean.rank")
  {
    obs.ranks <- matrix(0, ngenes, 3)
    obs.ranks[, 1] <- rank(modt)
    est <- sapply(index, function(o, modt) mean(modt[o]), modt)
    nua <- mean(modt)
    sda <- sd(modt)
    if (restand)
    {
      est <- (est - nua)/sda
    }
    obs.ranks[, 2] <- ngenes - obs.ranks[, 1] + 1
    obs.ranks[, 3] <- rank(abs(modt))
    AllIndices <- unlist(index)
    Set <- rep(1:nsets, SetSizes)
    obs.set.ranks <- rowsum(obs.ranks[AllIndices, ], group = Set,
                            reorder = FALSE)
    obs.set.ranks <- obs.set.ranks/SetSizes
    rot.ranks <- obs.ranks
    p.value <- matrix(0, nrow = nsets, ncol = 3)
    p.value <- mclapply(seq_len(nrot),function(i)
    {
      if(executation.info ) setTxtProgressBar(pb, i)
      R <- matrix(rnorm((d + 1)), 1, d + 1)
      R <- R/sqrt(rowSums(R^2))
      Br <- R %*% Y
      s2r <- (YY - Br^2)/d
      if (is.finite(d0))
        sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
      else sdr.post <- sqrt(s02)
      modtr <- signc * Br/sdr.post
      rot.ranks[, 1] <- rank(modtr)
      rot.ranks[, 2] <- ngenes - rot.ranks[, 1] + 1
      rot.ranks[, 3] <- rank(abs(modtr))
      rot.set.ranks <- rowsum(rot.ranks[AllIndices, ],
                              group = Set, reorder = FALSE)
      rot.set.ranks <- rot.set.ranks/SetSizes
      (rot.set.ranks <= obs.set.ranks)
   }, mc.cores = mccores)
   pval <- sapply(c(1,3),function(o) apply(sapply(p.value,function(x) x[,o]),1,sum))
   pval[,2] <- nrot - pval[,2]
   pval[,1] <- 2 * pmin(pval[,1],nrot - pval[,1])
   colnames(pval) <- c("pval.diff","pval.abs")
   rownames(pval) <- names(index)
   modt <- rank(modt)
  }

  if (set.statistic == "median")
  {
    est <- sapply(index, function(o, modt) median(modt[o]), modt)
    nua <- median(modt)
    sda <- mad(modt)
    if (restand)
    {
      est <- (est - nua)/sda
    }
    p.value <- matrix(0, nrow = nsets, ncol = 2)
    p.value <- apply(do.call(cbind,mclapply(seq_len(nrot),function(i)
    {
      if(executation.info ) setTxtProgressBar(pb, i)
      R <- matrix(rnorm((d + 1)), 1, d + 1)
      R <- R/sqrt(rowSums(R^2))
      Br <- R %*% Y
      s2r <- (YY - Br^2)/d
      if (is.finite(d0))
      {
        sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
      }else
      {
        sdr.post <- sqrt(s02)
      }
      modtr <- signc * Br/sdr.post
      est0 <- sapply(index, function(o, modtr) median(modtr[o]), modtr)
      nua0 <- median(modtr)
      sda0 <- mad(modtr)
      if (restand)
      {
        est0 <- (est0 - nua0)/sda0
      }
         p.value <- est0 <= est
    }, mc.cores = mccores)),1,sum)
    pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
  }

  if (set.statistic == "absmean")
  {
    modt <- abs(modt)
    est <- sapply(index, function(o, modt) mean(modt[o]), modt)
    nua <- mean(modt)
    sda <- sd(modt)
    if (restand)
    {
      est <- (est - nua)/sda
    }
    p.value <- matrix(0, nrow = nsets, ncol = 2)
    p.value <- apply(do.call(cbind,mclapply(seq_len(nrot),function(i)
    {
      if(executation.info ) setTxtProgressBar(pb, i)
      R <- matrix(rnorm((d + 1)), 1, d + 1)
      R <- R/sqrt(rowSums(R^2))
      Br <- R %*% Y
      s2r <- (YY - Br^2)/d
      if (is.finite(d0))
      {
        sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
      }else
      {
        sdr.post <- sqrt(s02)
      }
      modtr <- signc * Br/sdr.post
      est0 <- sapply(index, function(o, modtr) mean(abs(modtr[o])), modtr)
      nua0 <- mean(abs(modtr))
      sda0 <- sd(abs(modtr))
      if (restand)
      {
        est0 <- (est0 - nua0)/sda0
      }
          p.value <- est0 >= est
    }, mc.cores = mccores)),1,sum)
    pval <- as.matrix(p.value)
  }

  if (set.statistic == "GSEA")
  {
      ord <- order(-modt)
      modt2 <- modt[ord]
#      index2 <- lapply(index, function(o) which(ord%in%o))
      rankmod = rank(-modt, ties.method = "first")
      index2 <- lapply(index, function(o) rankmod[o])

      es <- sapply(index2, function(o, modt2) enrichmentScore(modt2, o), modt2)
      est <- apply(es, 2, function(o)
           {
             es.range <- range(o)
             es.range[abs(es.range)==max(abs(es.range))][1]
         })
    p.value <- matrix(0, nrow = nsets, ncol = 2)
    est0.all <- do.call(cbind,mclapply(seq_len(nrot),function(i)
    {
      if(executation.info ) setTxtProgressBar(pb, i)
      R <- matrix(rnorm((d + 1)), 1, d + 1)
      R <- R/sqrt(rowSums(R^2))
      Br <- R %*% Y
      (s2r <- (YY - Br^2)/d)

      if (is.finite(d0))
      {
        sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
      }else
      {
        sdr.post <- sqrt(s02)
      }
      modtr <- signc * Br/sdr.post
      ordr <- order(-modtr)
      modtr2 <- modtr[ordr]
      rankmod = rank(-modtr, ties.method = "first")
      index2 <- lapply(index, function(o) rankmod[o])
#      index2 <- lapply(index, function(o) which(ordr%in%o))
      es0 <- sapply(index2, function(o, modtr2) enrichmentScore(modtr2, o), modtr2)
      est0 <- apply(es0, 2, function(o)
              {
                es.range <- range(o)
                es.range[abs(es.range)==max(abs(es.range))][1]
              })
      est0 #  p.value <- est0 >= est #
   }, mc.cores = mccores))
   p.value <- sapply(seq_len(dim(est0.all)[1]), function(o) sum(est[o] >= as.numeric(est0.all[o,])))
   pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
   nes <- sapply(seq_len(dim(est0.all)[1]), function(o){
     est0.a <- est0.all[o,]
     est[o] / abs(mean(est0.a[sign(est0.a) == sign(est[o])]))
   })
  }
  if (set.statistic == "GSVA")
  {
    ord <- order(-modt)
    modt2 <- modt[ord]
#    index2 <- lapply(index, function(o) which(ord%in%o))
    rankmod = rank(-modt, ties.method = "first")
    index2 <- lapply(index, function(o) rankmod[o])
    es <- sapply(index2, function(o, modt2) enrichmentScore(modt2, o), modt2)
    estpos <- apply(es, 2, function(o) if (any(o > 0)) o[o==max(o)] else 0)
    estneg <- apply(es, 2, function(o) if (any(o < 0)) o[o==min(o)] else 0)
    est <- abs(estpos) - abs(estneg)
    p.value <- matrix(0, nrow = nsets, ncol = 2)
    p.value <- apply(do.call(cbind,mclapply(seq_len(nrot),function(i)
    {
      if(executation.info ) setTxtProgressBar(pb, i)
      R <- matrix(rnorm((d + 1)), 1, d + 1)
      R <- R/sqrt(rowSums(R^2))
      Br <- R %*% Y
      s2r <- (YY - Br^2)/d
      if (is.finite(d0))
      {
        sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
      }else
      {
        sdr.post <- sqrt(s02)
      }
      modtr <- signc * Br/sdr.post
      ordr <- order(-modtr)
      modtr2 <- modtr[ordr]
#      index2 <- lapply(index, function(o) which(ordr %in% o))
      rankmod = rank(-modt, ties.method = "first")
      index2 <- lapply(index, function(o) rankmod[o])
      es0 <- sapply(index2, function(o, modtr2) enrichmentScore(modtr2, o), modtr2)
      estpos0 <- apply(es0, 2, function(o) if (any(o > 0)) o[o==max(o)] else 0)
      estneg0 <- apply(es0, 2, function(o) if (any(o < 0)) o[o==min(o)] else 0)
      est0 <- abs(estpos0) - abs(estneg0)
         p.value <- est0 >= est
    }, mc.cores = mccores)),1,sum)
    pval <- as.matrix(2 * c(pmin(p.value, nrot-p.value)))
  }

  r <- (pval + 1)/(nrot + 1)
  rownames(r) <- names(index)
  if (set.statistic %in% c('GSEA',"mean", "maxmean"))
  {
    r <- cbind(NGenes = SetSizes, est=est, nes = nes, r)
    r <- data.frame(r)
    colnames(r) <- c("ngenes", "est", "nes", "pval")
    r$adj.pval <- p.adjust(r$pval, method = "BH")
    r$sign <- sign(r$est)
    r <- r[order(-abs(r$nes)), c("ngenes", "est", "nes", "pval", "adj.pval")]
  }
  else{
      r <- cbind(NGenes = SetSizes, est=est, r)
      r <- data.frame(r)
      if (set.statistic!='mean.rank')
      {
          colnames(r) <- c("ngenes", "est", "pval")
          r$adj.pval <- p.adjust(r$pval, method = "BH")
          r$sign <- sign(r$est)
          r <- r[order(r$pval), c("ngenes", "est", "pval", "adj.pval")]
      }else
      {
          colnames(r)[1] <- "ngenes"
          colnames(r)[2] <- "est"
          colnames(r)[3] <- "pval.diff"
          colnames(r)[4] <- "pval.mixed"
          r$adj.pval.diff <- p.adjust(r$pval.diff, method = "BH")
          r$adj.pval.mixed <- p.adjust(r$pval.mixed, method = "BH")
          r$sign <- sign(r$est)
          r <- r[order(r$pval.diff),
                 c("ngenes", "est", "pval.diff", "adj.pval.diff","pval.mixed",
                   "adj.pval.mixed")]
      }
  }
  return(list(res = r, stats = modt))
}
