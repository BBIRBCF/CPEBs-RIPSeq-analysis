## contrast matrix from formula and groups
contMatFunction <- function(covar, form, lgroups){
    lsigns <- c(sapply(1:length(lgroups), function(j) c(1, -1), simplify = F))
    conts <- sapply(1:length(lsigns), function(j)
                make.contrasts(d=covar, form=as.formula(form),
                               lsel=lgroups[[j]],
                               signs=lsigns[[j]]), simplify=F)
    cont.mat <- sapply(conts, function(o) o$cont)
    colnames(cont.mat) <- names(lgroups)
    return(cont.mat)
}
