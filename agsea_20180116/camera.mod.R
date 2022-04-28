rankSumTestWithCorrelation.mod <- function (index, statistics, correlation = 0, df = Inf)
{
    n <- length(statistics)
    r <- rank(statistics)
    r1 <- r[index]
    n1 <- length(r1)
    n2 <- n - n1
    U <- n1 * n2 + n1 * (n1 + 1)/2 - sum(r1)
    mu <- n1 * n2/2
    if (correlation == 0 || n1 == 1) {
        sigma2 <- n1 * n2 * (n + 1)/12
    }
    else {
        sigma2 <- asin(1) * n1 * n2 + asin(0.5) * n1 * n2 * (n2 -
            1) + asin(correlation/2) * n1 * (n1 - 1) * n2 * (n2 -
            1) + asin((correlation + 1)/2) * n1 * (n1 - 1) *
            n2
        sigma2 <- sigma2/2/pi
    }
    TIES <- (length(r) != length(unique(r)))
    if (TIES) {
        NTIES <- table(r)
        adjustment <- sum(NTIES * (NTIES + 1) * (NTIES - 1))/(n *
            (n + 1) * (n - 1))
        sigma2 <- sigma2 * (1 - adjustment)
    }
    zlowertail <- (U + 0.5 - mu)/sqrt(sigma2)
    zuppertail <- (U - 0.5 - mu)/sqrt(sigma2)
    pvalues <- c(less = pt(zuppertail, df = df, lower.tail = FALSE), stat = (U - 0.5 - mu))
    pvalues
}

camera.mod <- function(y, index, design, contrast = ncol(design), weights = NULL, set.statistic = "mean", allow.neg.cor=FALSE,
                        inter.gene.cor = NULL, trend.var = FALSE, sort = TRUE, df.camera = NULL, ...)
{
    dots <- names(list(...))
    if (length(dots))
        warning("Extra arguments disregarded: ", sQuote(dots))
    y <- getEAWP(y)
    G <- nrow(y$exprs)
    n <- ncol(y$exprs)
    ID <- rownames(y$exprs)
    use.ranks <- set.statistic == "mean.ranks"

    if (G < 3)
        stop("Two few genes in dataset: need at least 3")
    if (!is.list(index))
        index <- list(set1 = index)
    nsets <- length(index)
    if (nsets == 0L)
        stop("index is empty")
    if (is.null(design))
        design <- y$design
    if (is.null(design))
        stop("design matrix not specified")
    else {
        design <- as.matrix(design)
        if (mode(design) != "numeric")
            stop("design must be a numeric matrix")
    }
    if (nrow(design) != n)
        stop("row dimension of design matrix must match column dimension of data")
    p <- ncol(design)
    df.residual <- n - p
    if (df.residual < 1)
        stop("No residual df: cannot compute t-tests")
    if (is.null(weights))
        weights <- y$weights
    fixed.cor <- !(is.na(inter.gene.cor) || is.null(inter.gene.cor))
    if (fixed.cor) {
        if(is.null(df.camera)){
            if (use.ranks)
             df.camera <- Inf
            else df.camera <- G - 2
        }
    }
    else {
        if(is.null(df.camera)) df.camera <- min(df.residual, G - 2)
    }
    y <- y$exprs
    if (!is.null(weights)) {
        if (any(weights <= 0))
            stop("weights must be positive")
        if (length(weights) == n) {
            sw <- sqrt(weights)
            y <- t(t(y) * sw)
            design <- design * sw
            weights <- NULL
        }
    }
    if (!is.null(weights)) {
        if (length(weights) == G)
            weights <- matrix(weights, G, n)
        weights <- as.matrix(weights)
        if (any(dim(weights) != dim(y)))
            stop("weights not conformal with y")
    }
    if (is.character(contrast)) {
        contrast <- which(contrast == colnames(design))
        if (length(contrast) == 0)
            stop("coef ", contrast, " not found")
    }
    if (length(contrast) == 1) {
        j <- c((1:p)[-contrast], contrast)
        if (contrast < p)
            design <- design[, j]
    }
    else {
        QR <- qr(contrast)
        design <- t(qr.qty(QR, t(design)))
        if (sign(QR$qr[1, 1] < 0))
            design[, 1] <- -design[, 1]
        design <- design[, c(2:p, 1)]
    }
    if (is.null(weights)) {
        QR <- qr(design)
        if (QR$rank < p)
            stop("design matrix is not of full rank")
        effects <- qr.qty(QR, t(y))
        unscaledt <- effects[p, ]
        if (QR$qr[p, p] < 0)
            unscaledt <- -unscaledt
    }
    else {
        effects <- matrix(0, n, G)
        colnames(effects) <- ID
        unscaledt <- rep.int(0, G)
        names(unscaledt) <- ID
        sw <- sqrt(weights)
        yw <- y * sw
        for (g in 1:G) {
            xw <- design * sw[g, ]
            QR <- qr(xw)
            if (QR$rank < p)
                stop("weighted design matrix not of full rank for gene ",
                  g)
            effects[, g] <- qr.qty(QR, yw[g, ])
            unscaledt[g] <- effects[p, g]
            if (QR$qr[p, p] < 0)
                unscaledt[g] <- -unscaledt[g]
        }
    }
    U <- effects[-(1:p), , drop = FALSE]
    sigma2 <- colMeans(U^2)
    U <- t(U)/sqrt(pmax(sigma2, 1e-08))
    if (trend.var)
        A <- rowMeans(y)
    else A <- NULL
    sv <- squeezeVar(sigma2, df = df.residual, covariate = A)
    modt <- unscaledt/sqrt(sv$var.post)
    if (use.ranks)
        Stat <- modt
    else {
        df.total <- min(df.residual + sv$df.prior, G * df.residual)
        Stat <- zscoreT(modt, df = df.total, approx = TRUE)
    }
    meanStat <- mean(Stat)
    varStat <- var(Stat)
    tab <- matrix(0, nsets, 5)
    rownames(tab) <- names(index)
    colnames(tab) <- c("ngenes", "est", "pval","adj.pvalue","correlation")

    if(fixed.cor){
        if(length(inter.gene.cor)==1) inter.gene.cor <- rep(inter.gene.cor, length(index))
        else inter.gene.cor <- inter.gene.cor[names(index)]
    }

    for (i in 1:nsets) {
        iset <- index[[i]]
        if (is.character(iset))
            iset <- which(ID %in% iset)
        StatInSet <- Stat[iset]
        m <- length(StatInSet)
        m2 <- G - m
        if (fixed.cor) {
            correlation <- inter.gene.cor[i]
            vif <- 1 + (m - 1) * correlation
        }
        else {
            if (m > 1) {
                Uset <- U[iset, , drop = FALSE]
                vif <- m * mean(colMeans(Uset)^2)
                correlation <- (vif - 1)/(m - 1)
            }
            else {
                vif <- 1
                correlation <- NA
            }
        }
        tab[i, 1] <- m
        tab[i, 5] <- correlation
        if (use.ranks) {
            if (!allow.neg.cor)
                correlation <- max(0, correlation)
            statandpval <- rankSumTestWithCorrelation.mod(iset, statistics = Stat,
                                            correlation = correlation, df = df.camera)
            tab[i, 3] <- 2 * pmin(statandpval[1], 1- statandpval[1])
            tab[i, 2] <- statandpval[2]
        }
        else {
            if (!allow.neg.cor)
                vif <- max(1, vif)
            meanStatInSet <- mean(StatInSet)
            delta <- G/m2 * (meanStatInSet - meanStat)
            varStatPooled <- ((G - 1) * varStat - delta^2 * m *
                m2/G)/(G - 2)
            two.sample.t <- delta/sqrt(varStatPooled * (vif/m +
                1/m2))
            tab[i, 3] <- pt(two.sample.t, df = df.camera)
            tab[i, 3] <- 2 * pmin(tab[i, 3], 1 - tab[i, 3])
            tab[i, 2] <- delta
        }
    }
    tab <- data.frame(tab, stringsAsFactors = FALSE)
    if (fixed.cor)
        tab$correlation <- NULL
    if (nsets > 1)
        tab$adj.pvalue <- p.adjust(tab$pval, method = "BH")
    if (sort && nsets > 1) {
        o <- order(tab$pval)
        tab <- tab[o, ]
    }
    return(list(res = tab, stats = Stat))
}
