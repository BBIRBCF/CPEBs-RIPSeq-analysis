bestextpos <- function(x, y, limx, limy, down = TRUE, lsep = 5){
   xsep <- seq(limx[1],limx[2], length.out=lsep)
   wh <- list()
   for(i in 1:(length(xsep)-1))
       wh[[i]] <- which(x>xsep[i] & x<=xsep[i+1])

   ys <- list()
   for(o in 1:length(wh)){
       ys[[o]] <- numeric(length(wh[[o]]))
       if(length(wh[[o]])>0){
           if(down) ys[[o]][order(y[wh[[o]]])] <- seq(limy[1] + 0.05*diff(limy), min(y[wh[[o]]]), length.out = length(y[wh[[o]]]) + 1 )[1:length(y[wh[[o]]])]
       }
   }
   yresp <- unlist(ys)
   yresp[unlist(wh)] <- yresp
   return(yresp)
}



formatApval <- function(x){
    st = "ns"
    st <- ifelse(x < 0.25, "+", st)
    st <- ifelse(x < 0.10, "*", st)
    st <- ifelse(x < 0.05, "**", st)
    st <- ifelse(x < 0.01, "***", st)
    if(is.na(st)) st <- "nc"
    return(st)
}
densPlot <- function(x, y,...){
    cols <- densCols(x,y,   colramp=colorRampPalette(c("black", "white")))
    dens <- col2rgb(cols)[1,] + 1L
    cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
    col <- cols[dens]
    plot(x,y, col=col, ...)
    cc <- ic.ccc(x,y, alfa=0.05)

    legend("bottomright", legend=c(paste0("r = ", round(cor(x,y),3)), paste0("cb = ", round(cc[8],3)),  paste0("ccc = ", round(cc[1],3)) ), col=2,cex=1.5)

}

getGseaNES.ACM <- function(gseadirs)
    {
        nes <- lapply(gseadirs,function(x)
                      {
			  gsea.pos.lf <- list.files(x,pattern='gsea_report_for_na_pos_')
                          gsea.pos <- read.delim(file.path(x,gsea.pos.lf[regexpr("xls",gsea.pos.lf)>0]),as.is=TRUE,dec='.')
			  gsea.neg.lf <- list.files(x,pattern='gsea_report_for_na_neg_')
                          gsea.neg <- read.delim(file.path(x,gsea.neg.lf[regexpr("xls",gsea.neg.lf)>0]),as.is=TRUE,dec='.')

                          df <- rbind(gsea.pos[,c('NAME','NES','FDR.q.val')],gsea.neg[,c('NAME','NES','FDR.q.val')])
                          df <- df[order(df$NAME,decreasing=FALSE),]
                                        #ans <- df$NES; names(ans) <- df$NAME
                          rownames(df) <- df$NAME
                          return(df)
                      })
        ## Mod to fix the scarce but possible case in which one of the ranked lists doesn't have values for all genes in the other (i.e. zeros in RNASeq done with DESeq2), resulting in geneset being excluded for small size
        ## In that case, cells with no value for that geneset will report NA
        gsets <- sort(unique(unlist(lapply(nes,rownames))))
        ans <- data.frame(ID=gsets,do.call(cbind,lapply(nes,function(x) x[gsets,2:3])))
        rownames(ans) <- gsets
        ans
    }

getGseaNES.ACM.2 <- function(gseadirs, vars = c('NAME','NES','FDR.q.val','NOM.p.val'))
    {
        nes <- lapply(gseadirs,function(x)
                      {
                          gsea.pos.lf <- list.files(x,pattern='gsea_report_for_na_pos_')
                          gsea.pos <- read.delim(file.path(x,gsea.pos.lf[regexpr("xls",gsea.pos.lf)>0]),as.is=TRUE,dec='.')
                          gsea.neg.lf <- list.files(x,pattern='gsea_report_for_na_neg_')
                          gsea.neg <- read.delim(file.path(x,gsea.neg.lf[regexpr("xls",gsea.neg.lf)>0]),as.is=TRUE,dec='.')
                                         #print(colnames(gsea.pos))
                          df <- rbind(gsea.pos[,vars],gsea.neg[,vars])
                          df <- df[order(df$NAME,decreasing=FALSE),]
                                        #ans <- df$NES; names(ans) <- df$NAME
                          rownames(df) <- df$NAME
                          return(df)
                      })
        ## Mod to fix the scarce but possible case in which one of the ranked lists doesn't have values for all genes in the other (i.e. zeros in RNASeq done with DESeq2), resulting in geneset being excluded for small size
        ## In that case, cells with no value for that geneset will report NA
        gsets <- sort(unique(unlist(lapply(nes,rownames))))
        ans <- data.frame(ID=gsets,do.call(cbind,lapply(nes,function(x) x[gsets,2:dim(x)[2]])))
        rownames(ans) <- gsets
        ans
    }

ic.icc <- function(varind, vares, n, k, alfa)
{
### ICC y varianza

    icc <- varind/(varind + vares);
    var.icc <- (2*(1-icc)^2*(1+(k-1)*icc)^2)/(n*k*(k-1));


### Intervalos de confianza sin transformacion de Fisher

    delta1 <- qnorm(1-alfa/2)*sqrt(var.icc);
    ic1 <- icc + c(-delta1, delta1);

### Intervalos de confianza con transformacion de Fisher

    z <-  1/2*log((1+icc)/(1-icc));
    var.z <- var.icc/(1-icc^2)^2;
    delta2 <- qnorm(1-alfa/2)*sqrt(var.z);
    ic2.t <- z + c(-delta2, delta2);
    ic2 <- (exp(2*ic2.t)-1)/(exp(2*ic2.t)+1)

### Resultado

    list(icc=icc, var.icc=var.icc, z=z, var.z=var.z, ic1=ic1, ic2=ic2);
}


###############################################################################
    # Funcion para contrastes de hipotesis de icc

cont.icc <- function(icc, icc0, transform=T)
{
    if (transform==T)
    {
        z0 <- 1/2*log((1+icc0)/(1-icc0));
        est <- (icc$z -  z0)/(sqrt(icc$var.z));
        pval <- 1- pnorm(est);
    }
    else
    {
        est <- (icc$icc -  icc0)/(sqrt(icc$var.icc));
        pval <- 1- pnorm(est);
    }

    r <- c(est=est, pval=pval);
    r;

}

###############################################################################
    # ICC e intervalos con paquete nlme

ic.icc2 <- function(y, grp, alfa)
{

### Modelo de efectos mixtos

    m <- lme(y ~ 1, random=~1|grp);
    varind <- as.numeric(VarCorr(m)[1, 1]);
    vares <- as.numeric(VarCorr(m)[2, 1]);
    n <- length(unique(grp));
    k <- length(y)/n;


### ICC y varianza
    icc <- varind/(varind + vares);
    var.icc <- (2*(1-icc)^2*(1+(k-1)*icc)^2)/(n*k*(k-1));


### Intervalos de confianza sin transformacion de Fisher

    delta1 <- qnorm(1-alfa/2)*sqrt(var.icc);
    ic1 <- icc + c(-delta1, delta1);


### Intervalos de confianza con transformacion de Fisher

    z <-  1/2*log((1+icc)/(1-icc));
    var.z <- var.icc/(1-icc^2)^2;
    delta2 <- qnorm(1-alfa/2)*sqrt(var.z);
    ic2.t <- z + c(-delta2, delta2);
    ic2 <- (exp(2*ic2.t)-1)/(exp(2*ic2.t)+1)


### Resultado

    list(icc=icc, var.icc=var.icc, z=z, var.z=var.z, ic1=ic1, ic2=ic2);
}

###############################################################################
    # Funcion para estimacion de ccc

ic.ccc <- function(x, y, alfa=0.05)
{
    n <- length(x);
    ccc <- 2*cov(x, y)/((mean(x)-mean(y))^2 + var(x) + var(y));
    v <- min(c(sd(y)/sd(x), sd(x)/sd(y)));
    u <- (mean(y) - mean(x))/sqrt((sd(x)*sd(y)));
    cb <- 1/((v+1/v + u^2)/2);
    cp <- cor(x, y);
    se <- sqrt(1/(n-2)*((1-cp^2)*(1-ccc^2)*ccc^2/cp^2 + 2*u^2*(1-ccc)*ccc^3/cp -
                                 u^4*ccc^4/(2*cp^2)));

### Intervalo de confianza con transformacion de Fisher

    z <- 1/2*log((1+ccc)/(1-ccc));
    se.z <- sqrt(se^2/(1-ccc^2)^2);
    delta <- qnorm(1-alfa/2)*se.z;
    ciz <- z + c(-delta, delta);
    ci <- (exp(2*ciz) - 1)/(exp(2*ciz) + 1);

    c(ccc=ccc, se=se, alfa=alfa, llci=ci[1], upci=ci[2], u=u, v=v, cb=cb, r=cp,
      z=z, se.z=se.z);

}

###############################################################################
    # Funcion para contrastes de hipotesis de ccc

cont.ccc <- function(ccc, ccc0)
{
    z0 <- 1/2*log((1+ccc0)/(1-ccc0));
    est <- (ccc["z"] -  z0)/(ccc["se.z"]);
    pval <- 1- pnorm(est);

    r <- c(est=est, pval=pval);
    r;
}

###############################################################################
###############################################################################
###############################################################################




write.html.mod.sublinks <- function (x, links, tiny.pic, tiny.pic.size = 100, title = "",
    file, digits = 3, col.align='center', cellpadding=10, colspecial.links = 4, links2, topgenes)
{
    stopifnot(class(x) == "data.frame")
    if (missing(links))
        links <- vector("list", ncol(x))
    if (missing(tiny.pic))
        tiny.pic <- vector("list", ncol(x))
    stopifnot(class(links) == "list")
    stopifnot(class(tiny.pic) == "list")
    stopifnot(length(links) == ncol(x))
    stopifnot(length(tiny.pic) == ncol(x))
    stopifnot(!missing(file))
    column.class <- unlist(lapply(x, class))
    for (j in 1:ncol(x)) {
        if (column.class[j] == "factor")
            x[, j] <- as.character(x[, j])
        if (column.class[j] == "numeric")
            x[, j] <- round(x[, j], digits = digits)
    }
    cat("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n",
        sep = "", file = file)
    cat("<html>\n", file = file, append = T)
    cat("<body>\n", file = file, append = T)
    cat(paste("<CAPTION ALIGN=\"top\"><center><B>", title, "</B></center></CAPTION><BR>\n"),
        sep = "", file = file, append = T)
    cat(paste("<TABLE border=1 cellpadding=", cellpadding,">\n", sep=''), file = file, append = T)
    cat("<TR>\n", file = file, append = T)
    for (j in 1:ncol(x)) {
        cat("<TH>", file = file, append = T)
        cat(colnames(x)[j], file = file, append = T)
        cat("</TH>\n", file = file, append = T)
    }
    cat("</TR>\n", file = file, append = T)
    for (i in 1:nrow(x)) {
        cat("<TR>\n", file = file, append = T)
        for (j in seq_len(ncol(x))) {
            cat(paste("<TD align=", col.align, ">", sep=''), file = file, append = T)
            if(j != colspecial.links){
                if (is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
                    cat(x[i, j], file = file, append = T)
                }
                else if (is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
                    cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                              tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                              "\" width=\"", tiny.pic.size, "\" /></A>",
                              sep = ""), file = file, append = T)
                }
                else if (!is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
                    cat(paste("<A HREF=\"", links[[j]][[i]], "\">",
                              x[i, j], "</A>", sep = ""), file = file, append = T)
                }
                else if (!is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
                    cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                              tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                              "\" width=\"", tiny.pic.size, "\" /></A>",
                              sep = ""), file = file, append = T)
                }

            }
            else{
                for(k in 1:length(links2[i,])){
                    if (!is.na(links2[i,k])) {
                        cat(paste("<A HREF=\"", links2[i,k], "\">",
                              topgenes[i, k], " </A>", sep = ""), file = file, append = T)

                    }
                    else{
                        cat(paste0(topgenes[i, k], " "), file = file, append = T)
                    }
                }
            }
               cat("</TD>\n", file = file, append = T)
        }
        cat("</TR>\n", file = file, append = T)
    }

    cat("</TABLE>\n", file = file, append = T)
    cat("</body>\n", file = file, append = T)
    cat("</html>\n", file = file, append = T)
    sortDragHtmlTable(filename = file)
}


sortDragHtmlTable <- function (filename)
{
    lastSlashPos <- gregexpr(.Platform$file.sep, filename)[[1]][length(gregexpr(.Platform$file.sep,
        filename)[[1]])]
    outputdir <- ifelse(lastSlashPos == -1, getwd(), substr(filename,
        0, lastSlashPos))
    fileCode <- ""
    data(sorttable)
    writeLines(sorttable, con = paste(outputdir, "sorttable.js",
        sep = .Platform$file.sep))
    data(dragtable)
    writeLines(dragtable, con = paste(outputdir, "dragtable.js",
        sep = .Platform$file.sep))
    tmpTxt <- readLines(filename)
    tmpTxt[1] <- paste("<script src=\"sorttable.js\"></script>\n",
        tmpTxt[1])
    tmpTxt[1] <- paste("<script src=\"dragtable.js\"></script>\n",
        tmpTxt[1])
    tmpTxt <- sub("TABLE", "TABLE class=\"draggable sortable\"",
        tmpTxt)
    writeLines(tmpTxt, con = filename)
}
