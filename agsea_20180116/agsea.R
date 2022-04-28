################################################################################################################
#################################################################################################################
### Main functions
#################################################################################################################
#################################################################################################################

#source("/Volumes/biostats/acaballe/functions_for_all_analysis/GSEA_mod.R")
#source("/Volumes/biostats/acaballe/functions_for_all_analysis/utilsACM.R")
#source("/Volumes/biostats/acaballe/functions_for_all_analysis/make.contrasts.R")

#source("/Volumes/biostats/acaballe/GSEA_significance/code/romer.mod.R")
#source("/Volumes/biostats/acaballe/GSEA_significance/code/camera.mod.R")
#source("/Volumes/biostats/acaballe/GSEA_significance/code/gsa.mod.R")
#source("/Volumes/biostats/acaballe/GSEA_significance/code/fgsea.mod.R")
#source("/Volumes/biostats/acaballe/GSEA_significance/code/gsva.mod.R")
#source("/Volumes/biostats/acaballe/GSEA_significance/code/gsadj.mod.R")


#########
#statistics \ methods | romer          camera         roast        prgsea       gsa         gsva           gsadj
#-----------------------------------------------------------------------------------------------------------------
#                       mean           mean           mean                      mean
#                       mean.rank      mean.rank
#                       abs.mean                      abs.mean                  abs.mean
#                       median                        median
#                       max.mean                                                max.mean
#                       GSEA                                        V
#                       GSVA                                                                V
#                                                                                                           V
#-----------------------------------------------------------------------------------------------------------------

##############################################
######### main function:  agsea ##############
##############################################
# arguments:
# y: expression matrix;  covar: covariates for desgin; form: formula for DE;
# contrast: column or name of the contrast (if NA take last column of design matrix)
# counts: are y counts?; method: to chose from "romer", "gsadj", "roast", "camera", "gsea", "gsva", "gsa", "prgsea";
# set.statistic: to chose from  "mean", "mean.rank", "abs.mean", "max.mean", "GSEA", "GSVA"
# pr.stat: ranked score if method == GSEA; testvar = name of testing variable in covar
# rmGSGenes: (from GSA function); addGS: GS adjustment if method= gsadj?
# geneSymbol: genesymbol if counts data in entrez;  normalizationMethod: to chose from  "voom", "rlog", "zscoreNBinom",
# geneCorBatch: correct by adjusting variables in gsadj genewise a priori;
# inter.gene.cor: preset correlation, if null estimated (for camera)
# index: rows in y of genesets of interest (if NULL use info of gsetsel and gspath)
# gsetsel: name of pathway; gspath: dir of pathway,
# minsize: minimum size of geneset; maxsize: maximum size of geneset; nrot: number of rotations/permutations;
# mccores: number of cores for parallel computations;  geneRandom: randomization of gene label (just for diagnostic)
# comp: gsadj/gsva with competitive HT; romerMeanOpt: optimal romer/mean computation when total number of genes and number of pathways are not large.

agsea <- function(y, covar, form, contrast = NA, counts =FALSE, testvar = NA, method = "romer", set.statistic = "mean", inter.gene.cor = NULL, psel = NULL,
                  pr.stat = NA, lgroups = NA,  rmGSGenes ="gene", addGS = TRUE, geneSymbol = NULL, normalizationMethod = "rlog", geneCorBatch = FALSE,
                  index = NULL, gsetsel = "Hallmarks", gspath = "/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/", executation.info = TRUE,
                  minsize = 10, maxsize = 500, nrot = 9999, mccores = 1,  geneRandom = FALSE, comp = FALSE, romerMeanOpt = FALSE, ...){

    method <- match.arg(method, c("romer", "camera", "roast", "prgsea", "gsa", "gsva", "gsadj"))

    # Matrix design
    design <- model.matrix(as.formula(form), data = covar)
    terms <- colnames(design)
    if(is.na(contrast[1])) contrast = ncol(design)
    if(is.character(contrast)) contrast <- which(colnames(design) == contrast)

    # probset level?
    if(is.null(psel)) y2 <- y
    else{
        y2 <- y[psel,]
        rownames(y2) <- names(psel)
    }

    # Indeces
    if(is.null(index)){
        index.o <- indexFunction(y2, gsetsel, gspath = gspath, geneRandom = geneRandom, minsize = minsize, maxsize = maxsize, forPRgsea = TRUE)

        if(method %in% c("prgsea", "gsa")) index <- index.o
        else   index <- indexFunction(y2, gsetsel, gspath = gspath, geneRandom = geneRandom, minsize = minsize, maxsize = maxsize, forPRgsea = FALSE)

        index.o <- index.o[names(index)]
        index2 <- sapply(index.o, function(x) x[x%in%rownames(y2)])
    }
    else index2 <- index.o <- index

    # counts?
    if(counts){
        if(normalizationMethod == "voom") ynorm <- voom(y, design = design)
        if(normalizationMethod == "zscoreNBinom") ynorm <- zscoreDGE(y, covar, testvar, design, contrast = contrast)
        if(normalizationMethod == "rlog"){
            dds <- DESeqDataSetFromMatrix(countData = y[,rownames(covar)], colData=covar,design=~1)
            rlog.ass <- rlog(dds)
            ynorm <- assays(rlog.ass)[[1]]
        }
        ynormgene <- ynorm[(!is.na(geneSymbol)),]
        genes <- geneSymbol[(!is.na(geneSymbol))]

        nd <- ynormgene
        nds <- nd
        mvg <- sapply(unique(genes), function(o){
            wh <- which(genes == o)
            if(length(wh) == 1)
                return(rownames(nds)[wh])
            else
                return(rownames(nds)[wh[which.max(apply(nds[wh,],1,mad))]])
        })
        nd <- ynormgene
        nds <- nd[mvg,]
        rownames(nds) <- names(mvg)#tolower(names(mvg))
        y <- nds[apply(nds,1,var) > 0, ]
    }


    ########################### Results ######################
    if(method == "romer"){
         if(romerMeanOpt & set.statistic == "mean")
           res <- romer.modMeanOpt(y = y, index = index, design = design, contrast = contrast, set.statistic = set.statistic, nrot = nrot, mccores = mccores, ...)
         else
           res <- romer.mod(y = y, index = index, design = design, contrast = contrast, set.statistic = set.statistic, nrot = nrot, mccores = mccores, psel = psel, executation.info = executation.info, ...)
     }
    if(method == "camera")
        res <- camera.mod(y = y, index = index, design = design, contrast = contrast, inter.gene.cor = inter.gene.cor, set.statistic = set.statistic, ...)

    if(method == "roast"){
       set.statistic <-  match.arg(set.statistic, c("mean", "absmean", "median"))

       if(romerMeanOpt & set.statistic == "mean")
        res <- romer.modMeanOpt(y = y, index = index, design = design, contrast = contrast, set.statistic = set.statistic, nrot = nrot, mccores = mccores, restand = FALSE, ...)
       else
        res <- romer.mod(y = y, index = index, design = design, contrast = contrast, set.statistic = set.statistic, nrot = nrot, mccores = mccores, restand = FALSE, psel = psel, executation.info = executation.info, ...)
   }

    if(method == "prgsea"){
       if(any(is.na(pr.stat))) stop("pr.stat not properly defined")

       res <- fgsea.mod(index = index, pr.stat = pr.stat, nperm = nrot, minsize = minsize, maxsize = maxsize, mccores = mccores, ...)

    }

    if(method == "gsa"){
#        rownames(y) <- tolower(rownames(y))
        nd <- ExpressionSet(y, AnnotatedDataFrame(covar))
        res <- gsa.mod(nd, nd[[testvar]], gs = index,  nperm = nrot, method = set.statistic, minsize = minsize, maxsize = maxsize, rmGSGenes = rmGSGenes, ...)
    }

    if(method == "gsva")
        res <- gsvacompetitive(y = y,  design = design, covar = covar, lgroups = lgroups, index = index, nperm = nrot, mccores = mccores, ...)

    if(method == "gsadj")
       res <- gsadjcompetitive(y = y, design = design, contrast = contrast, index = index, nperm = nrot, mccores = mccores, addGS = addGS, geneCorBatch =geneCorBatch, ...)

    res$res <- data.frame(sapply(index.o,length)[rownames(res$res)],res$res)
    names(res$res)[c(1,2)] <- c("total_genes","overlap_genes")
    res$index <- index2
    res$covar <- covar
    res$testvar <- testvar
    res$statistic <- set.statistic
    res$design <- design
    res$contrast <- contrast
    class(res) <- c("agsea",method)
    return(res)
}


##############################################
#########  function:  print agsea ############
##############################################

print.agsea <- function(x, ...){
   if(inherits(x,"prgsea"))
       print(x$res[,!names(x$res)%in%"leadingEdge"])
   else
       print(x$res)
}

##############################################
#########  function:  html agsea #############
##############################################

htmlagsea <- function(res, htmlpath = "", htmlname = "file.html", plotpath ="", PLOT = TRUE, heatmap = TRUE, y, whplots = NULL, tit = "", topgenes = 10,
                      margins = c(5,16), dendrogram = "n", Rowv =TRUE, Colv = FALSE, mycol = c("black","red"), orderby = 1, scaling = FALSE, ... ){
    if(!inherits(res,"agsea")) stop("not a agsea object")

    x <- data.frame(geneset  = rownames(res$res),res$res)

    if(PLOT) dir.create(paste0(htmlpath,plotpath))
    if(PLOT & !inherits(res,c("gsadj","gsva"))){
        if(is.null(whplots)) whplots <- names(index)
        if(!is.na(whplots[1])){
            stats <- sort(res$stats)
            index <- sapply(res$index, function(x) which(names(stats)%in%x))
            for(k in whplots){
                png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),".png"))
                if(res$statistic %in% c("GSEA","GSVA")) plotGSEA(res, whplot = k, ...)
                if(res$statistic %in% c("absmean","mean", "mean.rank", "median","maxmean")) plotMean(res, whplot = k, ...)
                dev.off()
            }
         }
    }
    if(heatmap & PLOT){
        if(is.null(whplots)) whplots <- names(index)
        if(!is.na(whplots[1])){
            for(k in whplots){
                png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),"_heatmap.png"))
                plotindheatmap(y, res$index[[k]], res$covar,  res$testvar, margins = margins,
                               dendrogram = dendrogram, Rowv= Rowv, Colv = Colv, mycol = mycol,
                               orderby = orderby, scaling = scaling, design = res$design, contrast = res$contrast, ...)
                dev.off()
            }
         }


    }
    if(topgenes > 0){
        genesord <- sapply(res$index, function(ind) sort(res$stats[ind]))
        genestop <- sapply(rownames(res$res), function(i) if(res$res[i,"est"]<0) return(head(names(genesord[[i]]),min(topgenes, length(genesord[[i]])))) else  return(tail(names(genesord[[i]]),min(topgenes, length(genesord[[i]])))))
        if(topgenes>1) x$topgenes <- apply(genestop,2, function(y) paste(y,collapse =", "))
    }
    if(PLOT){
        x$plot <- NA
        if(heatmap) x$heatmap <- NA
    }
    links <- vector('list',length=ncol(x));
    names(links) <- colnames(x);
    plots <- links;
    if(PLOT)   links$plot <- plots$plot <- paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),".png")
    if(PLOT & heatmap)   links$heatmap <- plots$heatmap <- paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),"_heatmap.png")
    write.html.mod(x, file=paste0(htmlpath,htmlname), links=links,tiny.pic=plots,  title=tit)

}

##############################################
#########  function:  plot heatmap  ##########
##############################################

library(gplots)
## plotindheatmap <- function(y, geneset, covar,  testvar, design, contrast, margins = c(5,16), dendrogram = "n", Rowv =TRUE, Colv = FALSE, mycol = c("black","red"), orderby = 1, scaling = FALSE, filt = 3, ...){
##          resstats <- y
##          resstats <- resstats[,order(covar[[orderby]])]
##          resstats <- resstats[geneset,]
##          resstats <- resstats[apply(resstats,1,mean)>filt,]
##          design <- design[colnames(resstats),]
##          resstats <- batchcorlm(resstats, design, contrast)
##          if(scaling) resstats <- t(scale(t(resstats)))
##          covarord <- covar[order(covar[[orderby]]),]
##          mycol <- c("black","red")
##          heatmap.2(resstats, dendrogram = dendrogram, col=bluered(100),trace='none',margins=margins,notecol='black',notecex=1,keysize=.9,cexCol=1.5,
##                    ColSideColors = mycol[covarord[[testvar]]], Rowv= Rowv, Colv = Colv, ... )
## }

# new plotindheatmap by default does not adjust by GS
plotindheatmap <- function(y, geneset, covar,  testvar, design, contrast, margins = c(12,8), dendrogram = "n", Rowv =TRUE, Colv = FALSE, mycol = c("black","red"),
                           adjusting = FALSE, orderby = 1, scaling = FALSE, filt = 3, ...){
    resstats <- y
    resstats <- resstats[,order(covar[[orderby]])]
    resstats <- resstats[geneset,]
    resstats <- resstats[apply(resstats,1,mean)>filt,]
    #print(colnames(resstats))
    design <- design[colnames(resstats),]
    if(adjusting) resstats <- batchcorlm(resstats, design, contrast)
    if(scaling) resstats <- t(scale(t(resstats)))
    covarord <- covar[order(covar[[orderby]]),]
    mycol <- c("black","red")
    heatmap.2(resstats, dendrogram = dendrogram, col=bluered(100),trace='none',margins=margins,notecol='black',notecex=1,keysize=.9,cexCol=1.5,
              ColSideColors = mycol[covarord[[testvar]]], Rowv= Rowv, Colv = Colv, ... )
}


##############################################
#########  function:  plot gsadj #############
##############################################

library(gplots)
plotgsadj <- function(res,  margins = c(5,16), dendrogram = "n", Rowv =FALSE, Colv = FALSE, mycol = c("black","red"), orderby = 1, whplot = rownames(res$res), scaling = FALSE, ...){
   if(inherits(res,c("gsadj","gsva"))){
         covar <- res$covar
         resstats <- res$stats
         resstats <- resstats[,order(covar[[orderby]])]
         resstats <- resstats[whplot,]
         if(scaling) resstats <- t(scale(t(resstats)))
         covarord <- covar[order(covar[[orderby]]),]
         nbar <- as.numeric(as.factor(covarord[[res$testvar]]))
         heatmap.2(resstats, dendrogram = dendrogram, col=bluered(100),trace='none',margins=margins,notecol='black',notecex=1,keysize=.9,cexCol=1.5,
                   ColSideColors = mycol[nbar], Rowv= Rowv, Colv = Colv, ... )
     }
}

##############################################
#########  function:  plot GSEA  #############
##############################################

library(RColorBrewer)
plotGSEA <- function(obj, whplot = 1, maintitle = "",  ...) {
  res <- obj$res
  index <- obj$index[rownames(res[whplot,])]
  stats <- obj$stats

  modt2 <- stats[order(-stats)]
  if(is.list(index)) index2 <- lapply(index, function(o) which(names(modt2)%in%o))
  else index2 <-  list(which(names(modt2)%in%index))
  def.par <- par(no.readonly = TRUE)

  for(K in names(index)){
      if(maintitle == ""){
          maintitle <- K
          subtitle  <- paste(names(res[K,]),round( as.numeric(res[K,]),3),sep=" = ",collapse="; " )
      }
      else subtitle <- ""
      par(mar=c(1,4,5,1))
      layout(c(1,2),heights=c(4,2))

      es <- sapply(index2, function(o, modt2) enrichmentScore(modt2, o), modt2)
      plot(es[,K],type='l',axes=FALSE, ylab="ES", xlab='', oma=c(0,0,0,0), col='darkgreen',lwd=2, main = maintitle)
      mtext(subtitle, 3, line=0.5, col= "darkgray")
      abline(h=0)
      axis(2)
      abline(h=0)
      abline(v= which.max(abs(es[,K])), col=2,lty=3)

      par(mar=c(5,4,1,1))
      plot(NA,NA,xlim=c(1,length(es[,K])),ylim=c(0,1),axes=FALSE,ylab='',xlab='Gene list rank',sub="",col='green',lwd=2)
      myTicks <- axTicks(1); myTicks[myTicks==0] <- 1; axis(1,lty=0,at=myTicks)
      s <- rep(FALSE, length(es[,K])); s[index2[[K]]] <- TRUE
      if (sum(s)<50) myCols <- "darkgreen" else myCols <- densCols(which(s[1:length(s)]),colramp=colorRampPalette(brewer.pal(9, "Greens")[-c(1:3)]))
      abline(v=which(s[1:length(s)]),col=myCols)
      polygon(c(1,length(s)/15,length(s)/15,1), c(0-0.05,0-0.05,1+0.05,1+0.05), col= rgb(0.5, 0, 0,0.3), border = NA)
      polygon(c(length(s), length(s) - length(s)/15, length(s) - length(s)/15, length(s)), c(0-0.05,0-0.05,1+0.05,1+0.05), col= rgb(0,0,.5,0.3), border = NA)
      text(length(s)/30, 0.5, "+", cex=2)
      text(length(s)-length(s)/30, 0.5, "-", cex=1.5)
  }
  par(def.par)
}

##############################################################################
#########  function:  plot Mean  (assume restand == TRUE by now) #############
##############################################################################

plotMean <- function(obj, whplot = 1, maintitle = "", ylimAll = TRUE, ylim = NULL,  maxDensPlot = 20, ...) {
  res <- obj$res
  index <- obj$index[rownames(res[whplot,])]
  stats <- obj$stats
  statistic <- obj$statistic

  modt2 <- stats[order(-stats)]
  if(is.list(index)) index2 <- lapply(index, function(o) which(names(modt2)%in%o))
  else index2 <-  list(which(names(modt2)%in%index))
  def.par <- par(no.readonly = TRUE)

  if(inherits(obj, "romer")){
      if(statistic == "median") stt <- (modt2 - median(stats))/mad(stats)
      else stt <- (modt2 - mean(stats))/sd(stats)
  }
  else{
      stt <- modt2
  }
  es <- sapply(index2, function(o) stt[o])

  for(K in names(index2)){
      if(maintitle == ""){
          maintitle <- K
          subtitle  <- paste(names(res[K,]),round( as.numeric(res[K,]),3),sep=" = ",collapse="; " )
      }
      else subtitle <- ""

      le <- length(es[,K])
      if(is.null(ylim)){
          if(ylimAll) ylim <- c(min(stt),max(stt))
          else ylim <- c(min(es[,K]), max(es[,K]))
      }
      s <- rep(FALSE, length(stt)); s[index2[[K]]] <- TRUE
      if(sum(s) > maxDensPlot)
          layout(t(as.matrix(c(1,2,3))),heights=c(4),widths=c(4,.5,1))
      else layout(t(as.matrix(c(1,2))),heights=c(4),widths=c(4,1))

      par(mar=c(3,4,5,1))
      plot(1:le,es[,K], type='l',axes=FALSE, ylab="stat", xlab="ordered genes", oma=c(0,0,0,0),
           col='darkgreen',lwd=2, main = maintitle, ylim=ylim)
      polygon( c(1:le, le,0),c(es[,K],0,0),col=4, border = 4)
      mtext(subtitle, 3, line=0.5, col= "darkgray")
      axis(2)
      whchange <- ifelse(sum(es[,K]<0)>0, which(es[,K]<0)[1], le)
      axis(1, at = c(1,whchange,le))
      abline(v= whchange, col=1, lty=3)

      par(mar=c(3,1,5,1))
#      myTicks <- axTicks(2);
      s <- rep(FALSE, length(stt)); s[index2[[K]]] <- TRUE
      if (sum(s)<50) myCols <- "darkgreen" else myCols <- densCols(which(s[1:length(s)]),colramp=colorRampPalette(brewer.pal(9, "Greens")[-c(1:3)]))

      plot(NA,NA, ylim=ylim, xlim=c(0,1),axes=FALSE,xlab='',ylab='',sub="",col='green',lwd=2)
      abline(h=es[,K],col=myCols)
      maxmin <- ylim[2] - ylim[1]
      polygon( c(0-0.05,0-0.05,1+0.05,1+0.05), c(ylim[1], ylim[1] + maxmin/15, ylim[1] + maxmin/15, ylim[1]), col= rgb(0, 0, 0.5,0.3), border = NA)
      polygon(c(0-0.05,0-0.05,1+0.05,1+0.05), c(ylim[2], ylim[2] - maxmin/15, ylim[2] - maxmin/15, ylim[2]),  col= rgb(0.5, 0, 0,0.3), border = NA)

      text( 0.5,ylim[1] + maxmin/30, "-", cex=1.5)
      text( 0.5,ylim[2] - maxmin/30, "+", cex=1.5)

      if(sum(s)>20){
          par(mar=c(3,1,5,2))
          ds <- density(es[,K])
          plot(ds$y,ds$x, ylim=ylim, type="l", axes=FALSE,oma=c(0,0,0,0), main="", ylab="", col="darkblue")
          abline(h=0, col=1,lty=3)
      }
  #    axis(1)
  }
  par(def.par)
}


