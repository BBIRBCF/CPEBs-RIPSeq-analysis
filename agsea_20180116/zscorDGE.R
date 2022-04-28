zscoreDGE <- function (y, covar, testvar, design, contrast = ncol(design))
{
    y <- DGEList(y,  group = covar[[testvar]])
    y <- estimateDisp(y, design)
    allzero <- rowSums(y$counts > 1e-08) == 0

    if (any(allzero))  warning(sum(allzero), "rows with all zero counts")
    dispersion <- edgeR::getDispersion(y)
    if (is.null(dispersion)) stop("Dispersion estimate not found. Please estimate the dispersion(s) before you proceed.")
    if (is.null(design))
        design <- y$design
    if (is.null(design)) {
        if (nlevels(y$samples$group) < 2)
            stop("design not supplied and samples all belong to the same group")
        design <- model.matrix(~y$samples$group)
        rownames(design) <- colnames(y)
    }
    nbeta <- ncol(design)
    if (nbeta < 2)
        stop("design matrix must have at least two columns")
    if (is.character(contrast)) {
        if (length(contrast) > 1)
            stop("contrast should specify only one column of design")
        contrast <- which(contrast == colnames(design))
        if (!length(contrast))
            stop("contrast doesn't match any column of design")
    }
    if (length(contrast) == 1) {
        design0 <- design[, -contrast, drop = FALSE]
    }
    else {
        design <- contrastAsCoef(design, contrast = contrast, first = FALSE)$design
        design0 <- design[, -nbeta, drop = FALSE]
    }
    fit.null <- edgeR::glmFit(y, design0, prior.count = 0)
    y <- edgeR::zscoreNBinom(y$counts, mu = pmax(fit.null$fitted.values, 1e-17), size = 1/dispersion)
    y
}
