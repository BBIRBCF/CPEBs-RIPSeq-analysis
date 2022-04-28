library(data.table)
library(Matrix)

gsa.mod <- function (eset, fac, gs, nperm = NULL, tests = c("unpaired",
    "paired"), method = c("maxmean", "mean", "absmean"), minsize = 1,
    maxsize = Inf, restandardize = TRUE, npermBreaks = 2000,
    rmGSGenes = c("gene", "stop", "gs"), verbose = TRUE)
{
    fac <- (2-as.numeric(fac))
    res1 <- GSALight.mod(eset, fac, gs, nperm = nperm, tests = tests, method = method, minsize =minsize,
                    maxsize = maxsize, restandardize = restandardize, npermBreaks = npermBreaks,
                    rmGSGenes = rmGSGenes, verbose = verbose)
    res <- res1$res
    if(method=="absmean")
    {
        res <- as.data.frame(res[,c(4,3,1)])
        names(res) <- c("ngenes","est","pval")
        res$pval <-  c(res$pval*nperm+1)/(nperm+1)
        res$pval.adj <- p.adjust(res$pval, method = "BH")
    }
    else{
        res <- as.data.frame(res[,c(6,5,1)])
        names(res) <- c("ngenes","est","pval")
        res[,3] <- 2 * pmin(res$pval,1-res$pval)
        res$est <- -res$est
        res$pval <-  c(res$pval*nperm+1)/(nperm+1)
        res$pval.adj <- p.adjust(res$pval, method = "BH")
    }
    res <- res[order(res$pval),]
    return(list(res = res, stats = res1$stats))
}
 pvalFromPermMat.mod <- function (obs, perms)
{
    tempObs <- rep(obs, ncol(perms))
    dim(tempObs) <- dim(perms)
    Matrix::rowSums(perms <= tempObs)
}
data.table.mod <- function (..., keep.rownames = FALSE, check.names = FALSE, key = NULL,
    stringsAsFactors = FALSE)
{
    x <- list(...)
    .Call(CcopyNamedInList, x)
    if (identical(x, list(NULL)) || identical(x, list(list())) ||
        identical(x, list(data.frame(NULL))) || identical(x,
        list(data.table(NULL))))
        return(null.data.table())
    tt <- as.list(substitute(list(...)))[-1L]
    vnames = names(tt)
    if (is.null(vnames))
        vnames = rep.int("", length(x))
    vnames[is.na(vnames)] = ""
    novname = vnames == ""
    if (any(!novname)) {
        if (any(vnames[!novname] == ".SD"))
            stop("A column may not be called .SD. That has special meaning.")
    }
    for (i in which(novname)) {
        if (is.null(ncol(x[[i]]))) {
            if ((tmp <- deparse(tt[[i]])[1L]) == make.names(tmp))
                vnames[i] <- tmp
        }
    }
    tt = vnames == ""
    if (any(tt))
        vnames[tt] = paste0("V", which(tt))
    n <- length(x)
    if (n < 1L)
        return(null.data.table())
    if (length(vnames) != n)
        stop("logical error in vnames")
    vnames <- as.list.default(vnames)
    nrows = integer(n)
    numcols = integer(n)
    for (i in seq_len(n)) {
        xi = x[[i]]
        if (is.null(xi))
            stop("column or argument ", i, " is NULL")
        if ("POSIXlt" %chin% class(xi)) {
            warning("POSIXlt column type detected and converted to POSIXct. We do not recommend use of POSIXlt at all because it uses 40 bytes to store one date.")
            x[[i]] = as.POSIXct(xi)
        }
        else if (is.matrix(xi) || is.data.frame(xi)) {
            xi = as.data.table(xi, keep.rownames = keep.rownames)
            x[[i]] = xi
            numcols[i] = length(xi)
        }
        else if (is.table(xi)) {
            x[[i]] = xi = as.data.table.table(xi, keep.rownames = keep.rownames)
            numcols[i] = length(xi)
        }
        else if (is.function(xi)) {
            x[[i]] = xi = list(xi)
        }
        nrows[i] <- NROW(xi)
        if (numcols[i] > 0L) {
            namesi <- names(xi)
            if (length(namesi) == 0L)
                namesi = rep.int("", ncol(xi))
            namesi[is.na(namesi)] = ""
            tt = namesi == ""
            if (any(tt))
                namesi[tt] = paste0("V", which(tt))
            if (novname[i])
                vnames[[i]] = namesi
            else vnames[[i]] = paste(vnames[[i]], namesi, sep = ".")
        }
    }
    nr <- max(nrows)
    ckey = NULL
    recycledkey = FALSE
    for (i in seq_len(n)) {
        xi = x[[i]]
        if (is.data.table(xi) && haskey(xi)) {
            if (nrows[i] < nr)
                recycledkey = TRUE
            else ckey = c(ckey, key(xi))
        }
    }
    for (i in which(nrows < nr)) {
        xi <- x[[i]]
        if (identical(xi, list())) {
            x[[i]] = vector("list", nr)
            next
        }
        if (nrows[i] == 0L)
            stop("Item ", i, " has no length. Provide at least one item (such as NA, NA_integer_ etc) to be repeated to match the ",
                nr, " row", if (nr > 1L)
                  "s", " in the longest column. Or, all columns can be 0 length, for insert()ing rows into.")
        if (nr%%nrows[i] != 0L)
            warning("Item ", i, " is of size ", nrows[i], " but maximum size is ",
                nr, " (recycled leaving remainder of ", nr%%nrows[i],
                " items)")
        if (is.data.frame(xi)) {
            ..i = rep(seq_len(nrow(xi)), length.out = nr)
            x[[i]] = xi[..i, , drop = FALSE]
            next
        }
        if (is.atomic(xi) || is.list(xi)) {
            x[[i]] = rep(xi, length.out = nr)
            next
        }
        stop("problem recycling column ", i, ", try a simpler type")
    }
    if (any(numcols > 0L)) {
        value = vector("list", sum(pmax(numcols, 1L)))
        k = 1L
        for (i in seq_len(n)) {
            if (is.list(x[[i]]) && !is.ff(x[[i]])) {
                for (j in seq_len(length(x[[i]]))) {
                  value[[k]] = x[[i]][[j]]
                  k = k + 1L
                }
            }
            else {
                value[[k]] = x[[i]]
                k = k + 1L
            }
        }
    }
    else {
        value = x
    }
    vnames <- unlist(vnames)
    if (check.names)
        vnames <- make.names(vnames, unique = TRUE)
    setattr(value, "names", vnames)
    setattr(value, "row.names", .set_row_names(nr))
    setattr(value, "class", c("data.table", "data.frame"))
    if (!is.null(key)) {
        if (!is.character(key))
            stop("key argument of data.table() must be character")
        if (length(key) == 1L) {
            key = strsplit(key, split = ",")[[1L]]
        }
        setkeyv(value, key)
    }
    else {
        if (length(ckey) && !recycledkey && !any(duplicated(ckey)) &&
            all(ckey %in% names(value)) && !any(duplicated(names(value)[names(value) %in%
            ckey])))
            setattr(value, "sorted", ckey)
    }
    if (isTRUE(stringsAsFactors))
        setfactor(value, which(vapply(value, is.character, TRUE)),
            FALSE)
    alloc.col(value)
}

 dataTable2Mat.mod <- function (gsTable)
{
    tmp <- factor(gsTable$gene)
    tmpGenes <- levels(tmp)
    geneFactor <- as.numeric(tmp)
    tmp <- factor(gsTable$geneSet)
    tmpGS <- levels(tmp)
    gs <- as.numeric(tmp)
    mat <- sparseMatrix(gs, geneFactor, x = 1)
    if (any(mat > 1)) {
        warning("There appears to be duplicated genes in some gene sets. The duplicated genes are removed.")
        mat[mat > 1] <- 1
    }
    rownames(mat) <- tmpGS
    colnames(mat) <- tmpGenes
    mat
}

list2DataTable.mod <- function (geneList)
{
    gsLength <- sapply(geneList, length)
    data.table(geneSet = rep(names(geneList), times = gsLength),
        gene = unlist(geneList))
}

rowPairedTtests.mod <- function (eset, fac, method = c("maxmean", "mean", "absmean"))
{
    numSubject <- length(fac)/2
    mat <- Matrix(0, ncol(eset), numSubject)
    for (i in 1:numSubject) mat[which(abs(fac) == i), i] <- sign(fac[which(abs(fac) ==
        i)])
    d <- eset %*% mat
    meand <- rowMeans(d)
    vard <- Matrix::rowSums((d - meand)^2)/(numSubject - 1)
    results <- meand/(sqrt(vard)/sqrt(numSubject))
    if (method == "mean")
        return(results)
    else if (method == "absmean")
        return(abs(results))
    else if (method == "maxmean") {
        results1 <- pmax(results, 0)
        results2 <- -pmin(results, 0)
        return(ls = list(results1 = results1, results2 = results2))
    }
}
 rowtests.mod <- function (eset, fac, method = c("maxmean", "mean", "absmean"))
{
    fac <- as.numeric(as.factor(fac)) - 1
    if (sum(fac == 1) > sum(fac == 0))
        fac <- abs(fac - 1)
    numx <- sum(fac == 1)
    numy <- sum(fac == 0)
    x1 <- eset[, fac == 1]
    x2 <- eset[, fac == 0]
    mean1 <- rowMeans(x1)
    mean2 <- rowMeans(x2)
    var1 <- Matrix::rowSums(x1^2) - numx * mean1^2
    var2 <- Matrix::rowSums(x2^2) - numy * mean2^2
    results <- {
        mean1 - mean2
    }/sqrt({
        {
            var1 + var2
        }/{
            length(fac) - 2
        }
    } * {
        1/numx + 1/numy
    })
    if (method == "mean")
        return(results)
    else if (method == "absmean")
        return(abs(results))
    else if (method == "maxmean") {
        results1 <- pmax(results, 0)
        results2 <- -pmin(results, 0)
        return(ls = list(results1 = results1, results2 = results2))
    }
}

GSApairedfunc.mod <- function (eset, fac, nperm, method = c("maxmean", "mean", "absmean"))
{
    numSubject <- length(fac)/2
    mat <- Matrix(0, ncol(eset), numSubject)
    for (i in 1:numSubject) mat[which(abs(fac) == i), i] <- sign(fac[which(abs(fac) ==
        i)])
    d <- eset %*% mat
    permMat <- matrix(0, numSubject, nperm)
    for (i in 1:nperm) permMat[, i] <- rbinom(numSubject, 1,
        prob = 0.5)
    permMat <- Matrix(permMat)
    sumd <- Matrix::rowSums(d)
    sumpermd <- 2 * d %*% permMat
    meanx <- (sumd - sumpermd)/numSubject
    sumdsq <- Matrix::rowSums(d^2)
    sdx <- sqrt((sumdsq - numSubject * meanx^2)/(numSubject -
        1))
    resultsMat <- meanx/(sdx/sqrt(numSubject))
    if (method == "mean")
        return(resultsMat)
    else if (method == "absmean")
        return(abs(resultsMat))
    else if (method == "maxmean") {
        resultsMat1 <- pmax(as.matrix(resultsMat), 0)
        resultsMat2 <- -pmin(as.matrix(resultsMat), 0)
        return(ls = list(resultsMat1 = resultsMat1, resultsMat2 = resultsMat2))
    }
}

GSAfunc.mod <- function (eset, fac, nperm, method = c("maxmean", "mean", "absmean"))
{
    fac <- as.numeric(as.factor(fac)) - 1
    if (sum(fac == 1) > sum(fac == 0))
        fac <- abs(fac - 1)
    numX <- sum(fac == 1)
    numY <- sum(fac == 0)
    permMat <- fac %*% t(rep(1, nperm))
    for (i in 1:nperm) permMat[, i] <- sample(permMat[, i])
    permMat <- Matrix(permMat)
    sumAll <- Matrix::rowSums(eset)
    sumsqAll <- Matrix::rowSums(eset^2)
    SumX <- eset %*% permMat
    SumXsq <- eset^2 %*% permMat
    SumY <- sumAll - SumX
    SumYsq <- sumsqAll - SumXsq
    MeanX <- SumX/numX
    MeanY <- SumY/numY
    resultsMat <- {
        (MeanX - MeanY)/sqrt(SumXsq - numX * MeanX^2 + SumYsq -
            numY * MeanY^2)
    } * sqrt((length(fac) - 2)/(1/numX + 1/numY))
    if (method == "mean")
        return(resultsMat)
    else if (method == "absmean")
        return(abs(resultsMat))
    else if (method == "maxmean") {
        resultsMat1 <- pmax(as.matrix(resultsMat), 0)
        resultsMat2 <- -pmin(as.matrix(resultsMat), 0)
        return(ls = list(resultsMat1 = resultsMat1, resultsMat2 = resultsMat2))
    }
}
GSALight.mod <- function (eset, fac, gs, nperm = NULL, tests = c("unpaired",
    "paired"), method = c("maxmean", "mean", "absmean"), minsize = 1,
    maxsize = Inf, restandardize = TRUE, npermBreaks = 2000,
    rmGSGenes = c("stop", "gene", "gs"), verbose = TRUE)
{
    restandardize <- isTRUE(restandardize)
    verbose <- isTRUE(verbose)
    tests <- match.arg(tests)
    method <- match.arg(method)
    rmGSGenes <- match.arg(rmGSGenes)
    if (tests == "paired") {
        if (!is.integer(fac))
            stop(" For paired t-test, fac must be an integer vector.")
        if (any(fac == 0))
            stop(" For paired t-test, 0 is not allowed for subject pair index.")
        tab <- table(abs(fac))
        if (any(tab != 2) | sum(fac > 0) != sum(fac < 0))
            stop(" Some values in fac are not paired. For paired t-tests, fac must be a vector of 1,-1,2,-2,..., where each number represents a pair, and the sign represents the conditions. ")
        sortedfac <- sort(unique(abs(fac)))
        if (!all(sortedfac == c(1:(length(fac)/2))))
            stop(" Some values in fac are skipped. For paired t-tests, the numbering of the subject pairs must be 1,2,3,..., with no skipped integers. ")
    }
    if (tests == "unpaired") {
        if (!is.factor(fac))
            fac <- as.factor(fac)
        if (length(unique(fac)) > 2)
            stop("more than two classes detected. GSALightning only supports two-sample t-tests.")
    }
    if (is.data.table(gs))
        mat <- dataTable2Mat.mod(gs)
    else if (is.list(gs)) {
        if (is.null(names(gs)))
            stop("Gene set names missing: each element in the gene set list must be named.")
        gs <- list2DataTable.mod(gs)
        mat <- dataTable2Mat.mod(gs)
    }
    else mat <- gs
    rm(gs)
    if (is.null(rownames(eset)))
        stop("Gene names missing: each row of the expression data must be named after the gene.")
    if (any(is.na(eset)))
        stop("The Expression data contain missing values.")
    if (any(apply(eset, 1, var) == 0))
        stop("Some genes has 0 sample variance, please remove those genes prior to running GSALightning. Also consider removing genes with small sample variance.")
    if (!all(colnames(mat) %in% rownames(eset))) {
        if (rmGSGenes == "gene") {
            if (verbose)
                message("Some genes within the gene sets are not contained in the expression data set.\n These genes are removed from the gene sets since rmGSGenes == 'gene'.")
            mat <- mat[, colnames(mat) %in% rownames(eset)]
        }
        else if (rmGSGenes == "gs") {
            if (verbose)
                message("Some genes within the gene sets are not contained in the expression data set.\n Gene sets with missing genes are removed since rmGSGenes == 'gs'.")
            numGenes <- Matrix::rowSums(mat)
            newNumGenes <- Matrix::rowSums(mat[, colnames(mat) %in% rownames(eset)])
            mat <- mat[numGenes == newNumGenes, ]
        }
        else stop("Some genes within the gene sets are not contained in the expression data set.\n Set rmGSGenes = 'gene' or 'gs' to remove respectively the missing genes or gene sets.")
    }
    setsize <- Matrix::rowSums(mat)
    mat <- mat[setsize >= minsize & setsize <= maxsize, ]
    mat <- mat[, Matrix::colSums(mat) >= 1]
    eset <- eset[colnames(mat), ]
    numGenes <- Matrix::rowSums(mat)
    if (is.null(nperm)) {
        nperm <- nrow(mat)/0.05 * 2
        message("Number of permutations is not specified. Automatically set to ",
            nperm, ".")
    }
    if (verbose)
        message("After gene set size filtering, there are ",
            nrow(mat), " gene sets,\n containing a total of ",
            nrow(eset), " genes for analysis.")
    if (verbose)
        message("Obtaining observed gene set statistics.")
    if (tests == "paired")
        obs <- rowPairedTtests.mod(as.matrix(eset), fac, method)
    else obs <- rowtests.mod(as.matrix(eset), fac, method)
    stats <- obs
    if (restandardize) {
        if (method == "maxmean") {
            numCatGenes <- Matrix::colSums(mat)
            totCatGenes <- sum(numCatGenes)
            meanobs1 <- sum(obs$results1 * numCatGenes)/totCatGenes
            sdobs1 <- sqrt({
                sum(numCatGenes * {
                  obs$results1^2
                }) - totCatGenes * meanobs1^2
            }/(totCatGenes - 1))
            meanobs2 <- sum(obs$results2 * numCatGenes)/totCatGenes
            sdobs2 <- sqrt({
                sum(numCatGenes * {
                  obs$results2^2
                }) - totCatGenes * meanobs2^2
            }/(totCatGenes - 1))
        }
        else {
            numCatGenes <- Matrix::colSums(mat)
            totCatGenes <- sum(numCatGenes)
            meanobs <- sum(obs * numCatGenes)/totCatGenes
            sdobs <- sqrt({
                sum(numCatGenes * {
                  obs^2
                }) - totCatGenes * meanobs^2
            }/(totCatGenes - 1))
        }
    }
    if (restandardize) {
        if (method == "maxmean") {
            obs1 <- as.vector(mat %*% obs$results1)
            obs2 <- as.vector(mat %*% obs$results2)
            obs1 <- obsOrig1 <- (obs1/numGenes - meanobs1)/sdobs1
            obs2 <- obsOrig2 <- (obs2/numGenes - meanobs2)/sdobs2
            obs <- pmax(obs1, obs2)
            obs[obs2 > obs1] <- -1 * obs[obs2 > obs1]
            obsOrig <- pmax(obsOrig1, obsOrig2)
            obsOrig[obsOrig2 > obsOrig1] <- -1 * obsOrig[obsOrig2 >
                obsOrig1]
        }
        else {
            obs <- mat %*% obs
            obs <- as.vector(obs)
            obs <- obsOrig <- (obs/numGenes - meanobs)/sdobs *
                sqrt(numGenes)
        }
    }
    else {
        if (method == "maxmean") {
            obs1 <- as.vector(mat %*% obs$results1)
            obs2 <- as.vector(mat %*% obs$results2)
            obsOrig1 <- obs1/numGenes
            obsOrig2 <- obs2/numGenes
            obs <- pmax(obs1, obs2)
            obs[obs2 > obs1] <- -1 * obs[obs2 > obs1]
            obsOrig <- pmax(obsOrig1, obsOrig2)
            obsOrig[obsOrig2 > obsOrig1] <- -1 * obsOrig[obsOrig2 >
                obsOrig1]
        }
        else {
            obs <- mat %*% obs
            obs <- as.vector(obs)
            obsOrig <- obs/numGenes
        }
    }
    if (nperm <= npermBreaks) {
        if (tests == "paired")
            permMat <- GSApairedfunc.mod(as.matrix(eset), fac, nperm,
                method)
        else permMat <- GSAfunc.mod(as.matrix(eset), fac, nperm,
            method)
        if (restandardize) {
            if (method == "maxmean") {
                meanStar1 <- Matrix::colSums(permMat$resultsMat1 * numCatGenes)/totCatGenes
                sdStar1 <- sqrt((Matrix::colSums({
                  permMat$resultsMat1^2
                } * numCatGenes) - totCatGenes * meanStar1^2)/(totCatGenes -
                  1))
                permMat1 <- as.matrix(mat %*% permMat$resultsMat1)
                permMat1 <- t((t(permMat1/numGenes) - meanStar1)/sdStar1)
                meanStar2 <- Matrix::colSums(permMat$resultsMat2 * numCatGenes)/totCatGenes
                sdStar2 <- sqrt((Matrix::colSums({
                  permMat$resultsMat2^2
                } * numCatGenes) - totCatGenes * meanStar2^2)/(totCatGenes -
                  1))
                permMat2 <- as.matrix(mat %*% permMat$resultsMat2)
                permMat2 <- t((t(permMat2/numGenes) - meanStar2)/sdStar2)
                permMat <- pmax(permMat1, permMat2)
                permMat[permMat2 > permMat1] <- -1 * permMat[permMat2 >
                  permMat1]
            }
            else {
                meanStar <- Matrix::colSums(permMat * numCatGenes)/totCatGenes
                sdStar <- sqrt((Matrix::colSums({
                  permMat^2
                } * numCatGenes) - totCatGenes * meanStar^2)/(totCatGenes -
                  1))
                permMat <- as.matrix(mat %*% permMat)
                permMat <- t((t(permMat/numGenes) - meanStar)/sdStar) *
                  sqrt(numGenes)
            }
        }
        else {
            if (method == "maxmean") {
                permMat1 <- as.matrix(mat %*% permMat$resultsMat1)
                permMat2 <- as.matrix(mat %*% permMat$resultsMat2)
                permMat <- pmax(permMat1, permMat2)
                permMat[permMat2 > permMat1] <- -1 * permMat[permMat2 >
                  permMat1]
            }
            else permMat <- as.matrix(mat %*% permMat)
        }
        pvalSums <- pvalFromPermMat.mod(obs, permMat)
        pval <- pvalSums/nperm
    }
    else {
        if (verbose)
            message("Running batch-mode permutation.")
        permBreaks <- ceiling(nperm/npermBreaks)
        pvalSums <- rep(0, nrow(mat))
        for (i in 1:permBreaks) {
            if (verbose)
                message("Running batch ", i, " of ", permBreaks,
                  " batches.")
            if (tests == "paired")
                permMat <- GSApairedfunc.mod(as.matrix(eset), fac,
                  ifelse(i != permBreaks, npermBreaks, nperm -
                    {
                      i - 1
                    } * npermBreaks), method)
            else permMat <- GSAfunc.mod(as.matrix(eset), fac, ifelse(i !=
                permBreaks, npermBreaks, nperm - {
                i - 1
            } * npermBreaks), method)
            if (restandardize) {
                if (method == "maxmean") {
                  meanStar1 <- Matrix::colSums(permMat$resultsMat1 *
                    numCatGenes)/totCatGenes
                  sdStar1 <- sqrt((Matrix::colSums({
                    permMat$resultsMat1^2
                  } * numCatGenes) - totCatGenes * meanStar1^2)/(totCatGenes -
                    1))
                  permMat1 <- as.matrix(mat %*% permMat$resultsMat1)
                  permMat1 <- t((t(permMat1/numGenes) - meanStar1)/sdStar1)
                  meanStar2 <- Matrix::colSums(permMat$resultsMat2 *
                    numCatGenes)/totCatGenes
                  sdStar2 <- sqrt((Matrix::colSums({
                    permMat$resultsMat2^2
                  } * numCatGenes) - totCatGenes * meanStar2^2)/(totCatGenes -
                    1))
                  permMat2 <- as.matrix(mat %*% permMat$resultsMat2)
                  permMat2 <- t((t(permMat2/numGenes) - meanStar2)/sdStar2)
                  permMat <- pmax(permMat1, permMat2)
                  permMat[permMat2 > permMat1] <- -1 * permMat[permMat2 >
                    permMat1]
                }
                else {
                  meanStar <- Matrix::colSums(permMat * numCatGenes)/totCatGenes
                  sdStar <- sqrt((Matrix::colSums({
                    permMat^2
                  } * numCatGenes) - totCatGenes * meanStar^2)/(totCatGenes -
                    1))
                  permMat <- as.matrix(mat %*% permMat)
                  permMat <- t((t(permMat/numGenes) - meanStar)/sdStar) *
                    sqrt(numGenes)
                }
            }
            else {
                if (method == "maxmean") {
                  permMat1 <- as.matrix(mat %*% permMat$resultsMat1)
                  permMat2 <- as.matrix(mat %*% permMat$resultsMat2)
                  permMat <- pmax(permMat1, permMat2)
                  permMat[permMat2 > permMat1] <- -1 * permMat[permMat2 >
                    permMat1]
                }
                else permMat <- as.matrix(mat %*% permMat)
            }
            pvalSums <- pvalSums + pvalFromPermMat.mod(obs, permMat)
        }
        pval <- pvalSums/nperm
    }
    if (verbose)
        message("Permutation done. Evaluating P-values.")
    if (method == "absmean") {
        pvals <- 1 - pval
        qvals <- p.adjust(pvals, method = "BH")
        results <- cbind(pvals, qvals, obsOrig, Matrix::rowSums(mat))
        colnames(results) <- c("p-value", "q-value", "statistics",
            "# genes")
    }
    else {
        if (tests == "unpaired") {
            lvls <- levels(as.factor(fac))
            fac <- as.numeric(as.factor(fac)) - 1
            if (sum(fac == 1) > sum(fac == 0))
                lvls <- rev(lvls)
            pvals <- cbind(pval, 1 - pval)
            rownames(pvals) <- rownames(mat)
            qvals <- cbind(p.adjust(pvals[, 1], "BH"), p.adjust(pvals[,
                2], "BH"))
            results <- cbind(pvals, qvals, obsOrig, Matrix::rowSums(mat))
            colnames(results) <- c(paste("p-value:up-regulated in",
                lvls[1]), paste("p-value:up-regulated in", lvls[2]),
                paste("q-value:up-regulated in", lvls[1]), paste("q-value:up-regulated in",
                  lvls[2]), paste("statistics: up-regulated in",
                  lvls[2]), "# genes")
        }
        else {
            pvals <- cbind(1 - pval, pval)
            qvals <- cbind(p.adjust(pvals[, 1], "BH"), p.adjust(pvals[,
                2], "BH"))
            results <- cbind(pvals, qvals, obsOrig, Matrix::rowSums(mat))
            colnames(results) <- c(paste("p-value:up-regulated in positives"),
                paste("p-value:up-regulated in negatives"), paste("q-value:up-regulated in positives"),
                paste("q-value:up-regulated in negatives"), "statistics (up-regulated in positives)",
                "# genes")
        }
    }
    return(list(res = results, stats =stats))
}
