###########################################################################
########### GSEA roast for ripseq: bduran_201903_ripseq  ##################
###########################################################################

library(Biobase)
library(MDA)
library(limma)
library(phenoTest)
library(ggplot2)
library(quantreg)
library(affy)
library(affyPLM)
library(genefilter)
library(hwriter)
library(openxlsx)
library(pheatmap)

basedir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq'
datadir <- file.path(basedir,'data/')
reportsdir <- file.path(basedir,'reports/')
tablesdir <- file.path(basedir,'reports/tables/')
figsdir <- file.path(basedir,'reports/figs/')
res <- paste0(reportsdir,"/RoastGSEA/")
dir.create(res)

###########################################################################
# 1. Functions to source
###########################################################################

s1 <- list.files(file.path(basedir,'routines/agsea_20180116/'),'*.R',full.names=TRUE)
for (s in s1) source(s)
s2 <- list.files(file.path(basedir,'routines/functions_acaballe_20180116/'),'*.R',full.names=TRUE)
for (s in s2) source(s)

htmlagsea2 <- function(res, htmlpath = "", htmlname = "file.html", plotpath ="", PLOT = TRUE, heatmap = TRUE, y, whplots = NULL, tit = "", topgenes = 10, links2, genestopaux,
                      margins = c(5,16), dendrogram = "n", Rowv =TRUE, Colv = FALSE, mycol = c("black","red"), orderby = 1, scaling = FALSE, filt = 5, ... ){
    if(!inherits(res,"agsea")) stop("not a agsea object")
    ##
    x <- data.frame(geneset  = rownames(res$res),res$res)
    ##
    if(PLOT) dir.create(paste0(htmlpath,plotpath))
    if(PLOT & !inherits(res,c("gsadj","gsva"))) {
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
                pdf(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),"_heatmap.pdf"), width = 12, height = 15)
                plotindheatmap(y, res$index[[k]], res$covar,  res$testvar, margins = margins,
                               dendrogram = dendrogram, Rowv= Rowv, Colv = Colv, mycol = mycol,
                               orderby = orderby, scaling = scaling, filt = filt, design = res$design, contrast = res$contrast, ...)
                dev.off()
            }
         }
    }
    ##
    x$geneDEinfo <- rep("view", dim(x)[1])
    ##
    if(topgenes > 0){
        genesord <- sapply(res$index, function(ind) sort(res$stats[ind]))
        genestop <- sapply(rownames(res$res), function(i) if(res$res[i,"est"]<0) return(head(names(genesord[[i]]),min(topgenes, length(genesord[[i]])))) else  return(tail(names(genesord[[i]]),min(topgenes, length(genesord[[i]])))))
        if(topgenes>1) x$topgenes <- apply(genestop,2, function(y) paste(y,collapse =", "))
    }
    ##
    if(PLOT){
        x$plot <- NA
        if(heatmap) x$heatmap <- NA
    }
    links <- vector('list',length=ncol(x));
    names(links) <- colnames(x);
    plots <- links;
    if(PLOT)   links$plot <- plots$plot <- paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),".png")
    if(PLOT & heatmap)   links$heatmap <- plots$heatmap <- paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),"_heatmap.png")
#    links$geneDEinfo <- paste0(names(res), "indpath/",gsub("", "", gsub("[[:punct:]]", "_", rowanmes(x))), "_genes.html")
#    if(dim(x)[1]>200) links$geneDEinfo[201:dim(x)[1]] <- NA
    write.html.mod.sublinks(x, file=paste0(htmlpath,htmlname), links=links,tiny.pic=plots,  title=tit,  colspecial.links = 8, links2 = links2, topgenes = genestopaux)
#    write.html.mod(x, file=paste0(htmlpath,htmlname), links=links,tiny.pic=plots,  title=tit)
}

htmlagsea3 <- function(res, htmlpath = "", htmlname = "file.html", plotpath ="", PLOT = TRUE, heatmap = TRUE, y, whplots = NULL, tit = "", topgenes = 10, links2, genestopaux,
                      margins = c(10,10), dendrogram = "n", Rowv =TRUE, Colv = FALSE, mycol = c("black","red"), orderby = 1, scaling = FALSE, filt = 5, namfile = "Hallmarks",... ){
    if(!inherits(res,"agsea")) stop("not a agsea object")
    ##
    x <- data.frame(geneset  = rownames(res$res),res$res)
    ##    rownames(y) <- toupper(rownames(y))
    if(PLOT) dir.create(paste0(htmlpath,plotpath))
    if(PLOT & !inherits(res,c("gsadj","gsva"))){
        if(is.null(whplots)) whplots <- names(index)
        if(!is.na(whplots[1])){
            stats <- sort(res$stats)
            index <- sapply(res$index, function(x) which(names(stats)%in%x))
            for(k in whplots){
                png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),"_mean.png"))
                plotMean(res, whplot = k, ...)
                dev.off()
                png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),"_gsea.png"))
                plotGSEA(res, whplot = k, ...)
                dev.off()
            }
         }
    }
    if(heatmap & PLOT){
        if(is.null(whplots)) whplots <- names(index)
        if(!is.na(whplots[1])){
            for(k in whplots){
                ## PDF heatmap
                pdf(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),"_heatmap.pdf"), width = 12, height = 15)
                ## Scaled
                plotindheatmap(y, (res$index[[k]]), res$covar,  res$testvar, margins = margins,
                               dendrogram = dendrogram, Rowv= Rowv, Colv = Colv, mycol = mycol,
                               orderby = orderby, scaling = scaling, filt = filt, design = res$design, contrast = res$contrast, main='Scaled',...)
                ## Scaled
                plotindheatmap(y, (res$index[[k]]), res$covar,  res$testvar, margins = margins,
                               dendrogram = dendrogram, Rowv= Rowv, Colv = Colv, mycol = mycol,
                               orderby = orderby, scaling = FALSE, filt = filt, design = res$design, contrast = res$contrast, main='Not scaled',...)
                dev.off()
                ## PNG heatmap
                png(paste0(htmlpath,plotpath, gsub("[[:punct:]]"," ",k),"_heatmap.png"), width = 600, height = 800)
                plotindheatmap(y, (res$index[[k]]), res$covar,  res$testvar, margins = margins,
                               dendrogram = dendrogram, Rowv= Rowv, Colv = Colv, mycol = mycol,
                               orderby = orderby, scaling = scaling, filt = filt, design = res$design, contrast = res$contrast, ...)
                dev.off()
            }
         }
    }
    ##
    x$geneDEinfo <- rep("view", dim(x)[1])
    ##
    ##    if(topgenes > 0){
    ##        genesord <- sapply(res$index, function(ind) sort(res$stats[ind]))
    ##        genestop <- sapply(rownames(res$res), function(i) if(res$res[i,"est"]<0) return(head(names(genesord[[i]]),min(topgenes, length(genesord[[i]])))) else  return(tail(names(genesord[[i]]),min(topgenes, length(genesord[[i]])))))
    ##        if(topgenes>1) x$topgenes <- apply(genestop,2, function(y) paste(y,collapse =", "))
    ##    }
    ##
    if(PLOT){
        x$plot_mean <- NA
        x$plot_gsea <- NA
        if(heatmap) x$heatmap <- NA
    }
    links <- vector('list',length=ncol(x));
    names(links) <- colnames(x);
    plots <- links;
    if(PLOT)   links$plot_mean <- plots$plot_mean <- paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),"_mean.png")
    if(PLOT)   links$plot_gsea <- plots$plot_gsea <- paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),"_gsea.png")
    ##
    if(PLOT & heatmap)   links$heatmap <- plots$heatmap <- paste0(plotpath, gsub("[[:punct:]]"," ",rownames(x)),"_heatmap.png")
    links$heatmap <- gsub('.png','.pdf',links$heatmap,fixed=TRUE) # Show png thumbnail but point to PDF
    links$geneDEinfo <- paste0(namfile, "indpath/",gsub("", "", gsub("[[:punct:]]", "_", rownames(x))),
                               "_genes.html")
    if(dim(x)[1]>200) links$geneDEinfo[201:dim(x)[1]] <- NA
    ##    write.html.mod.sublinks(x, file=paste0(htmlpath,htmlname), links=links,tiny.pic=plots,  title=tit,  colspecial.links = 8, links2 = links2, topgenes = genestopaux)
    write.html.mod(x, file=paste0(htmlpath,htmlname), links=links,tiny.pic=plots,  title=tit)
    ##
}

###########################################################################
# 1. GENERATE RLOG MATRIX AS USED IN DESEQ2 ANALYSIS
###########################################################################

library(DESeq2)
library(edgeR)

# Load counts, sampleinfo and normalized counts
load(file=file.path(datadir,'featureCounts_bowtie2_genes_allbatches_multioverlap_minq1_notmultimap.RData'))
load(file=file.path(datadir,'sampleinfo_allbatches_forDESeq2_may19.RData'))
load(file=file.path(datadir,'bowtie2_genes_allbatches_dds_DESeq2_fcounts_minq1_may19.RData'))

# Count matrix
y <- fcounts.mapq1$counts
substr(colnames(y),84,100)
colnames(y) <- substr(colnames(y),84,100)
colnames(y) %in% sampleinfo$ID
y <- y[,as.character(sampleinfo$ID)]
colnames(y) == sampleinfo$ID
covar <- sampleinfo
rownames(covar) <- colnames(y) <- covar$Name
testvar <- "Group"

###########################################################################
# 2.2 Normalization
###########################################################################

# Get rlog matrix
dds      <- DESeqDataSetFromMatrix(countData = y, colData = covar, design = ~ 1)
dds      <- estimateSizeFactors(dds)
rlog.ass <- rlog(dds)
ynorm    <- assays(rlog.ass)[[1]]

# Remove non variable genes and low counts
y.bak <- ynorm
y <- ynorm
y <- y[!apply(y,1,var)==0,]
y <- y[apply(y,1,sum)>10,]

covar$nameid2 <- covar$Name
covar$nameid2 <- factor(covar$nameid2, levels = covar$nameid2)
covar$Replicate <- factor(covar$Replicate)

###########################################################################
# 3.  Romer and zscore + (with or without gsadj) signature activation
###########################################################################

## Geneset paths
gspath <- '/Volumes/biostats/databases/BroadGSEA/genesets/xlaevis/'
topgenes <- 10
## Geneset to be tested
gsetselALL <- c("GOSLIM", "KEGG","GOCC","GOMF","GOBP")

###########################################################################
########## 3.1. Romer with maxmean statistic, IP NI/CPEBs VS INPUT NI
###########################################################################

##########
## IP NI
##########

covar.sub <- covar[covar$Group %in% c('IP_NI','Input_NI'),]
for (i in names(covar.sub)) covar.sub[,i] <- as.factor(as.character(covar.sub[,i]))
covar.sub$Group <- relevel(covar.sub$Group,ref='Input_NI')
form <- "~ Group"
design <- model.matrix(as.formula(form), data = covar.sub)
colnames(design) <- gsub('Group','',colnames(design))
design

## report directories
dir.create(paste0(res,'IP_NI_vs_INPUT_NI'))
res4 <- paste0(res, "IP_NI_vs_INPUT_NI/romer_Rdataall/")
res3 <- paste0(res, "IP_NI_vs_INPUT_NI/romer_xlsall/")
res2 <- paste0(res, "IP_NI_vs_INPUT_NI/romer_htmlall/")
dir.create(res2)
dir.create(res3)
dir.create(res4)

# Subset counts
y.sub <- y[,as.character(covar.sub$Name)]
# Remove genes from short arm, as they are the same homolog
sel <- grep('.S',rownames(y.sub),fixed=TRUE)
y.sub <- y.sub[-sel,]
rownames(y.sub) <- gsub('.L','',rownames(y.sub),fixed=TRUE) # as they appear in the genesets without suffix

## agsea
romer.maxmean <- list()
for(gsetsel in gsetselALL) {
    romer.maxmean[[gsetsel]] <- agsea(y.sub, covar.sub, form, testvar = testvar, method = "romer", set.statistic = "maxmean", gsetsel = gsetsel, gspath = gspath, minsize = 10, maxsize = 500, nrot = 2000, mccores = 12)
}
save(romer.maxmean, file= paste0(res4, "filemaxmean_NI.Rdata"))

##########
## CPEB1
##########

covar.sub <- covar[covar$Group %in% c('IP_CPEB1','Input_NI'),]
for (i in names(covar.sub)) covar.sub[,i] <- as.factor(as.character(covar.sub[,i]))
covar.sub$Group <- relevel(covar.sub$Group,ref='Input_NI')
form <- "~ Group"
design <- model.matrix(as.formula(form), data = covar.sub)
colnames(design) <- gsub('Group','',colnames(design))
design

## report directories
dir.create(paste0(res,'IP_CPEB1_vs_INPUT_NI'))
res4 <- paste0(res, "IP_CPEB1_vs_INPUT_NI/romer_Rdataall/")
res3 <- paste0(res, "IP_CPEB1_vs_INPUT_NI/romer_xlsall/")
res2 <- paste0(res, "IP_CPEB1_vs_INPUT_NI/romer_htmlall/")
dir.create(res2)
dir.create(res3)
dir.create(res4)

# Subset counts
y.sub <- y[,as.character(covar.sub$Name)]
# Remove genes from short arm, as they are the same homolog
sel <- grep('.S',rownames(y.sub),fixed=TRUE)
y.sub <- y.sub[-sel,]
rownames(y.sub) <- gsub('.L','',rownames(y.sub),fixed=TRUE) # as they appear in the genesets without suffix

## agsea
romer.maxmean <- list()
for(gsetsel in gsetselALL) {
    romer.maxmean[[gsetsel]] <- agsea(y.sub, covar.sub, form, testvar = testvar, method = "romer", set.statistic = "maxmean", gsetsel = gsetsel, gspath = gspath, minsize = 10, maxsize = 500, nrot = 2000, mccores = 12)
}
save(romer.maxmean, file= paste0(res4, "filemaxmean_CPEB1.Rdata"))

##########
## CPEB2
##########

covar.sub <- covar[covar$Group %in% c('IP_CPEB2','Input_NI'),]
for (i in names(covar.sub)) covar.sub[,i] <- as.factor(as.character(covar.sub[,i]))
covar.sub$Group <- relevel(covar.sub$Group,ref='Input_NI')
form <- "~ Group"
design <- model.matrix(as.formula(form), data = covar.sub)
colnames(design) <- gsub('Group','',colnames(design))
design

## report directories
dir.create(paste0(res,'IP_CPEB2_vs_INPUT_NI'))
res4 <- paste0(res, "IP_CPEB2_vs_INPUT_NI/romer_Rdataall/")
res3 <- paste0(res, "IP_CPEB2_vs_INPUT_NI/romer_xlsall/")
res2 <- paste0(res, "IP_CPEB2_vs_INPUT_NI/romer_htmlall/")
dir.create(res2)
dir.create(res3)
dir.create(res4)

# Subset counts
y.sub <- y[,as.character(covar.sub$Name)]
# Remove genes from short arm, as they are the same homolog
sel <- grep('.S',rownames(y.sub),fixed=TRUE)
y.sub <- y.sub[-sel,]
rownames(y.sub) <- gsub('.L','',rownames(y.sub),fixed=TRUE) # as they appear in the genesets without suffix

## agsea
romer.maxmean <- list()
for(gsetsel in gsetselALL) {
    romer.maxmean[[gsetsel]] <- agsea(y.sub, covar.sub, form, testvar = testvar, method = "romer", set.statistic = "maxmean", gsetsel = gsetsel, gspath = gspath, minsize = 10, maxsize = 500, nrot = 2000, mccores = 12)
}
save(romer.maxmean, file= paste0(res4, "filemaxmean_CPEB2.Rdata"))

##########
## CPEB3
##########

covar.sub <- covar[covar$Group %in% c('IP_CPEB3','Input_NI'),]
for (i in names(covar.sub)) covar.sub[,i] <- as.factor(as.character(covar.sub[,i]))
covar.sub$Group <- relevel(covar.sub$Group,ref='Input_NI')
form <- "~ Group"
design <- model.matrix(as.formula(form), data = covar.sub)
colnames(design) <- gsub('Group','',colnames(design))
design

## report directories
dir.create(paste0(res,'IP_CPEB3_vs_INPUT_NI'))
res4 <- paste0(res, "IP_CPEB3_vs_INPUT_NI/romer_Rdataall/")
res3 <- paste0(res, "IP_CPEB3_vs_INPUT_NI/romer_xlsall/")
res2 <- paste0(res, "IP_CPEB3_vs_INPUT_NI/romer_htmlall/")
dir.create(res2)
dir.create(res3)
dir.create(res4)

# Subset counts
y.sub <- y[,as.character(covar.sub$Name)]
# Remove genes from short arm, as they are the same homolog
sel <- grep('.S',rownames(y.sub),fixed=TRUE)
y.sub <- y.sub[-sel,]
rownames(y.sub) <- gsub('.L','',rownames(y.sub),fixed=TRUE) # as they appear in the genesets without suffix

## agsea
romer.maxmean <- list()
for(gsetsel in gsetselALL) {
    romer.maxmean[[gsetsel]] <- agsea(y.sub, covar.sub, form, testvar = testvar, method = "romer", set.statistic = "maxmean", gsetsel = gsetsel, gspath = gspath, minsize = 10, maxsize = 500, nrot = 2000, mccores = 12)
}
save(romer.maxmean, file= paste0(res4, "filemaxmean_CPEB3.Rdata"))

##########
## CPEB4
##########

covar.sub <- covar[covar$Group %in% c('IP_CPEB4','Input_NI'),]
for (i in names(covar.sub)) covar.sub[,i] <- as.factor(as.character(covar.sub[,i]))
covar.sub$Group <- relevel(covar.sub$Group,ref='Input_NI')
form <- "~ Group"
design <- model.matrix(as.formula(form), data = covar.sub)
colnames(design) <- gsub('Group','',colnames(design))
design

## report directories
dir.create(paste0(res,'IP_CPEB4_vs_INPUT_NI'))
res4 <- paste0(res, "IP_CPEB4_vs_INPUT_NI/romer_Rdataall/")
res3 <- paste0(res, "IP_CPEB4_vs_INPUT_NI/romer_xlsall/")
res2 <- paste0(res, "IP_CPEB4_vs_INPUT_NI/romer_htmlall/")
dir.create(res2)
dir.create(res3)
dir.create(res4)

# Subset counts
y.sub <- y[,as.character(covar.sub$Name)]
# Remove genes from short arm, as they are the same homolog
sel <- grep('.S',rownames(y.sub),fixed=TRUE)
y.sub <- y.sub[-sel,]
rownames(y.sub) <- gsub('.L','',rownames(y.sub),fixed=TRUE) # as they appear in the genesets without suffix

## agsea
romer.maxmean <- list()
for(gsetsel in gsetselALL) {
    romer.maxmean[[gsetsel]] <- agsea(y.sub, covar.sub, form, testvar = testvar, method = "romer", set.statistic = "maxmean", gsetsel = gsetsel, gspath = gspath, minsize = 10, maxsize = 500, nrot = 2000, mccores = 12)
}
save(romer.maxmean, file= paste0(res4, "filemaxmean_CPEB4.Rdata"))

###########################################################################
########## 3.2. Romer with maxmean statistic, XLSX REPORTS + HEATMAP
###########################################################################

# Write xlsx
for (cpeb in 1:4)
    {
        ## dirs
        res4 <- paste0(res, sprintf("IP_CPEB%d_vs_INPUT_NI/romer_Rdataall/",cpeb))
        res3 <- paste0(res, sprintf("IP_CPEB%d_vs_INPUT_NI/romer_xlsall/",cpeb))
        res2 <- paste0(res, sprintf("IP_CPEB%d_vs_INPUT_NI/romer_htmlall/",cpeb))
        ## load romer results
        load(paste0(res4, sprintf("filemaxmean_CPEB%d.Rdata",cpeb)))
        for(gsetsel in gsetselALL) {
            x <- romer.maxmean[[gsetsel]]$res
            x <- data.frame(rownames(x), x)
            colnames(x)[1] <- "gene set"
            write.xlsx(x, file=paste0(res3, gsetsel,"_filemaxmean.xlsx"), sep="\t", row.names = FALSE)
        }
    }

# Write xlsx IP NI
for (cpeb in 'NI')
    {
        ## dirs
        res4 <- paste0(res, sprintf("IP_%s_vs_INPUT_NI/romer_Rdataall/",cpeb))
        res3 <- paste0(res, sprintf("IP_%s_vs_INPUT_NI/romer_xlsall/",cpeb))
        res2 <- paste0(res, sprintf("IP_%s_vs_INPUT_NI/romer_htmlall/",cpeb))
        ## load romer results
        load(paste0(res4, sprintf("filemaxmean_%s.Rdata",cpeb)))
        for(gsetsel in gsetselALL) {
            x <- romer.maxmean[[gsetsel]]$res
            x <- data.frame(rownames(x), x)
            colnames(x)[1] <- "gene set"
            write.xlsx(x, file=paste0(res3, gsetsel,"_filemaxmean.xlsx"), sep="\t", row.names = FALSE)
        }
    }

##################
# Heatmaps GOSLIM
##################

res.goslim <- lapply(c('NI',paste0('CPEB',1:4)),function(cpeb)
                     {
                         ## dirs
                         res4 <- paste0(res, sprintf("IP_%s_vs_INPUT_NI/romer_Rdataall/",cpeb))
                         ## load romer results
                         load(paste0(res4, sprintf("filemaxmean_%s.Rdata",cpeb)))
                         x <- romer.maxmean[['GOSLIM']]$res
                         x
                     })
res.goslim <- lapply(res.goslim,function(x) { x <- x[rownames(res.goslim[[1]]),]; x })
## nes
nes <- as.data.frame(do.call(cbind,lapply(res.goslim,function(x) x[,'nes'])))

## fdr and symbol
fdr <- as.data.frame(do.call(cbind,lapply(res.goslim,function(x) x[,'adj.pval'])))
rownames(nes) <- rownames(fdr) <- rownames(res.goslim[[1]])
colnames(nes) <- colnames(fdr) <- c('NI',paste0('CPEB',1:4))
fdrsym <- matrix(data='',nrow=nrow(nes),ncol=5)
fdrsym[fdr<0.1] <- '+'
fdrsym[fdr<0.05] <- '*'
fdrsym[fdr<0.01] <- '**'
fdrsym[fdr<0.001] <- '***'

## heatmap pdf
pdf(file.path(figsdir,'RoastGSEA_Heatmap_GOSLIM_withNI.pdf'),width=10,height=18)
pheatmap(nes,display_numbers=round(fdr,4))
pheatmap(nes,display_numbers=fdrsym)
dev.off()

##################
# Heatmaps KEGG
##################

res.kegg <- lapply(c('NI',paste0('CPEB',1:4)),function(cpeb)
                   {
                       ## dirs
                       res4 <- paste0(res, sprintf("IP_%s_vs_INPUT_NI/romer_Rdataall/",cpeb))
                       ## load romer results
                       load(paste0(res4, sprintf("filemaxmean_%s.Rdata",cpeb)))
                       x <- romer.maxmean[['KEGG']]$res
                       x
                   })

res.kegg <- lapply(res.kegg,function(x) { x <- x[rownames(res.kegg[[1]]),]; x })
## nes
nes <- as.data.frame(do.call(cbind,lapply(res.kegg,function(x) x[,'nes'])))
## fdr and symbol
fdr <- as.data.frame(do.call(cbind,lapply(res.kegg,function(x) x[,'adj.pval'])))
rownames(nes) <- rownames(fdr) <- rownames(res.kegg[[1]])
colnames(nes) <- colnames(fdr) <- c('NI',paste0('CPEB',1:4))
fdrsym <- matrix(data='',nrow=nrow(nes),ncol=5)
fdrsym[fdr<0.1] <- '+'
fdrsym[fdr<0.05] <- '*'
fdrsym[fdr<0.01] <- '**'
fdrsym[fdr<0.001] <- '***'

## heatmap pdf
pdf(file.path(figsdir,'RoastGSEA_Heatmap_KEGG_withNI.pdf'),width=10,height=18)
pheatmap(nes,display_numbers=round(fdr,4))
pheatmap(nes,display_numbers=fdrsym)
dev.off()

###########################################################################
########## 3.3. Romer with maxmean statistic, HTML REPORTS
###########################################################################

for (cpeb in 1:4)
    {
        ## dirs
        res4 <- paste0(res, sprintf("IP_CPEB%d_vs_INPUT_NI/romer_Rdataall/",cpeb))
        res3 <- paste0(res, sprintf("IP_CPEB%d_vs_INPUT_NI/romer_xlsall/",cpeb))
        res2 <- paste0(res, sprintf("IP_CPEB%d_vs_INPUT_NI/romer_htmlall/",cpeb))
        ## load romer results
        load(paste0(res4, sprintf("filemaxmean_CPEB%d.Rdata",cpeb)))
        ##
        for(gsetsel in gsetselALL[-2]) {
            print(gsetsel)
            aa <- romer.maxmean[[gsetsel]]
            genesord <- sapply(aa$index, function(ind) sort(aa$stats[ind]))
            genestop <-sapply(rownames(aa$res), function(i) if(aa$res[i,"est"]<0) return(head(names(genesord[[i]]),min(topgenes,length(genesord[[i]]))))else  return(tail(names(genesord[[i]]),min(topgenes, length(genesord[[i]])))))
            ##
            linksbase <- paste0("../../tables/Stripcharts/")
            topww <- dir(paste0(res2,linksbase))
            genesha <- (gsub(".png","",sapply(strsplit(topww,"...",fixed=TRUE), "[",2)))
            ##
            AL <- array(NA, dim=dim(genestop))
            for(i in 1:dim(genestop)[2]){
                for(l in 1:length(genestop[,i])){
                    aux <- which(genesha%in%genestop[l,i])
                    if(length(aux)>=1) AL[l,i] <- paste0(linksbase, topww[aux[1]])
                }
            }
            ##
            nn <- romer.maxmean[[gsetsel]]$covar[,"nameid2"][order(romer.maxmean[[gsetsel]]$covar[,"nameid2"])]
            romer.maxmean[[gsetsel]]$covar[,"nameid2"] <- factor(romer.maxmean[[gsetsel]]$covar[,"nameid2"], levels = c(nn[1:4],nn[5:7]))
            ##
            if (nrow(romer.maxmean[[gsetsel]]$res<200)) selplots <- 1:nrow(romer.maxmean[[gsetsel]]$res) else selplots <- 1:100
            htmlagsea3(romer.maxmean[[gsetsel]], htmlpath = res2, htmlname = paste0(gsetsel,"_filemaxmean.html"), plotpath = paste0(gsetsel,"_barcodesmaxmean/"), PLOT = TRUE, heatmap=TRUE,
                       whplots = rownames(romer.maxmean[[gsetsel]]$res)[selplots], tit = "", topgenes = 10, links2 = t(AL), genestopaux = t(genestop ),namfile = gsetsel,
                       y = y[,as.character(covar.sub$Name)], margins = c(15,8), dendrogram = "n", Rowv =TRUE, Colv = FALSE, mycol = c("black","red"), orderby = "nameid2", scaling = TRUE, filt =2)
            ##
            x <- romer.maxmean[[gsetsel]]$res
            x <- data.frame(rownames(x), x)
            colnames(x)[1] <- "gene set"
            write.table(x, file=paste0(res3, gsetsel,"_filemaxmean.xls"), sep="\t", row.names = FALSE)
        }
    }

# add NI results
for (cpeb in 'NI')
    {
        ## dirs
        res4 <- paste0(res, sprintf("IP_%s_vs_INPUT_NI/romer_Rdataall/",cpeb))
        res3 <- paste0(res, sprintf("IP_%s_vs_INPUT_NI/romer_xlsall/",cpeb))
        res2 <- paste0(res, sprintf("IP_%s_vs_INPUT_NI/romer_htmlall/",cpeb))
        ## load romer results
        load(paste0(res4, sprintf("filemaxmean_%s.Rdata",cpeb)))
        ##
        for(gsetsel in gsetselALL[-2]) { # kegg crashes..
            print(gsetsel)
            aa <- romer.maxmean[[gsetsel]]
            genesord <- sapply(aa$index, function(ind) sort(aa$stats[ind]))
            genestop <-sapply(rownames(aa$res), function(i) if(aa$res[i,"est"]<0) return(head(names(genesord[[i]]),min(topgenes,length(genesord[[i]]))))else  return(tail(names(genesord[[i]]),min(topgenes, length(genesord[[i]])))))
            ##
            linksbase <- paste0("../../tables/Stripcharts/")
            topww <- dir(paste0(res2,linksbase))
            genesha <- (gsub(".png","",sapply(strsplit(topww,"...",fixed=TRUE), "[",2)))
            ##
            AL <- array(NA, dim=dim(genestop))
            for(i in 1:dim(genestop)[2]){
                for(l in 1:length(genestop[,i])){
                    aux <- which(genesha%in%genestop[l,i])
                    if(length(aux)>=1) AL[l,i] <- paste0(linksbase, topww[aux[1]])
                }
            }
            ##
            nn <- romer.maxmean[[gsetsel]]$covar[,"nameid2"][order(romer.maxmean[[gsetsel]]$covar[,"nameid2"])]
            romer.maxmean[[gsetsel]]$covar[,"nameid2"] <- factor(romer.maxmean[[gsetsel]]$covar[,"nameid2"], levels = c(nn[1:4],nn[5:7]))
            ##
            if (nrow(romer.maxmean[[gsetsel]]$res<200)) selplots <- 1:nrow(romer.maxmean[[gsetsel]]$res) else selplots <- 1:100
            htmlagsea3(romer.maxmean[[gsetsel]], htmlpath = res2, htmlname = paste0(gsetsel,"_filemaxmean.html"), plotpath = paste0(gsetsel,"_barcodesmaxmean/"), PLOT = TRUE, heatmap=TRUE,
                       whplots = rownames(romer.maxmean[[gsetsel]]$res)[selplots], tit = "", topgenes = 10, links2 = t(AL), genestopaux = t(genestop ),namfile = gsetsel,
                       y = y[,as.character(covar.sub$Name)], margins = c(15,8), dendrogram = "n", Rowv =TRUE, Colv = FALSE, mycol = c("black","red"), orderby = "nameid2", scaling = TRUE, filt =2)
            ##
            x <- romer.maxmean[[gsetsel]]$res
            x <- data.frame(rownames(x), x)
            colnames(x)[1] <- "gene set"
            write.table(x, file=paste0(res3, gsetsel,"_filemaxmean.xls"), sep="\t", row.names = FALSE)
        }
    }


### THE END
