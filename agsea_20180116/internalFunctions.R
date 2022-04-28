
gpcor <- function(nds, nins = 500, mccores = 4, sigsize = 20, form = NULL, name.coef = "treat")
{
 y  <- t(scale(t(exprs(nds))))
 GS <- apply(y,2,mean)
 nds$GS <- GS

 if (is.null(form)) form <- formula("tr ~ treat + GS")

 coefs <- mclapply(seq_len(nins), function(o){
    sam2 <- sample(rownames(y), sigsize)
    tr <- apply(y[sam2,],2,mean)
    nds$tr <- tr

    coefs <- summary(lm(form, data = pData(nds)))$coef
    coefs <- coefs[regexpr(name.coef,rownames(coefs))>0,3:4]

    c(qnorm(coefs[2]/2) * sign(coefs[1]),coefs[2])
 }, mc.cores = mccores)
 coefs <- do.call(rbind,coefs)

 return(coefs)
}


## average correlation
avgcorfun <- function(x){
    x <- t(scale(t(x)))
    p <- dim(x)[1]
    n <- dim(x)[2]
    c1 <- ((rep(1,p) %*% x) %*% (t(x) %*% as.matrix(rep(1,p)))) / ((n-1)*p^2)
    return( (c1 - 1/p) * p/(p-1))
}


formatApval <- function(x){
    st <- sapply(x, function(x){
        st = ""
        st <- ifelse(x < 0.25, "+", st)
        st <- ifelse(x < 0.10, "*", st)
        st <- ifelse(x < 0.05, "**", st)
        st <- ifelse(x < 0.01, "***", st)
        if(is.na(st)) st <- ""
        return(st)})
    return(st)
}


## design from formula and expressionSet
designFunction <- function(nd, form = NULL){
    if(is.null(form)) form <-  paste0("~ -1 + ev", sep='')
    design <- model.matrix(as.formula(form), data=nd)
    return(design)
}

## contrast matrix from formula and groups
contMatFunction <- function(pd, form = NULL, lgroups = NULL){
    lsigns <- c(sapply(1:length(lgroups), function(j) c(1, -1), simplify = F))
    conts <- sapply(1:length(lsigns), function(j)
                make.contrasts(d=pd, form=as.formula(form),
                               lsel=lgroups[[j]],
                               signs=lsigns[[j]]), simplify=F)
    cont.mat <- sapply(conts, function(o) o$cont)
    colnames(cont.mat) <- names(lgroups)
    return(cont.mat)
}

## differential expression using Limma
dsfunction <- function(nd, form = NULL, lgroups = NULL){
 design <- designFunction(nd, form = form)
 lmf <- lmFit(nd, design)
 betas <- coef(lmf)
 cont.mat <- contMatFunction(pData(nd), form = form, lgroups = lgroups)

 lmc <- contrasts.fit(lmf, cont.mat)
 eb <- eBayes(lmc)
 ds <- topTable(eb, adjust="BH", number=nrow(exprs(nd)), coef=colnames(cont.mat)[1], confint=T)

 return(ds)
}
