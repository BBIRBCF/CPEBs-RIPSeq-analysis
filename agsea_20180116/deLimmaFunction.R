## differential expression using Limma
dsfunction <- function(y, covar, form, lgroups){

 design <- designFunction(nd, form = form)
 lmf <- lmFit(nd, design)
 betas <- coef(lmf)
 cont.mat <- contMatFunction(nd, form = form, lgroups = lgroups)

 lmc <- contrasts.fit(lmf, cont.mat)
 eb <- eBayes(lmc)
 ds <- topTable(eb, adjust="BH", number=nrow(exprs(nd)), coef=colnames(cont.mat)[1], confint=T)

 return(ds)
}
