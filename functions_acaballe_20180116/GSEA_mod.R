#################################################################################################################
#################################################################################################################
### GSEA_mod.R: functions to perform and GSEA and leading edge analysis using the java tool installed in pipeline
### Author: Toni Berenguer
#################################################################################################################
#################################################################################################################
    # Oscar R functions for GSEA pre-ranked execution in mulan via ssh to pipeline, modified to be used
    # in a R interactive R session

runGSEApreRanked.mod <- 
	function(rnk, 
                 customgeneset=NULL, score='weighted',summarize='Median_of_probes', 
                 minSize=15, maxSize=500, numplots=25, 
                 permutations=1000,
                 label='my_analysis', seed=149, 
                 outdir=getwd(), 
                 chip='gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip',
                 jarpath='/Volumes/biostats/soft/BroadGSEA/gsea2-2.0.12.jar',
                 userid='aberenguer') 
{
if (FALSE)
{
### Types of scores and gene summarizes available 
score <- c('weighted','weighted_p2','weighted_p1.5','classic')
list.summarize <- c('Median_of_probes','Max_probe');
}

  #if (!file.exists(rankedfile)) stop('Ranked File file does not exist!')
  # preProcess arguments
  # Broad GSEA Parameters

  memory <- '-Xmx1g'
  tool <- 'xtools.gsea.GseaPreranked'
  gmx <- customgeneset
  if (!file.exists(rnk)) stop('RNK File does not exist!')

  rungsea <-
    paste(sprintf('ssh %s@pipeline java -cp %s', userid, jarpath),
          sprintf('%s',memory),
          sprintf('%s',tool),
          sprintf('-gmx %s',gmx),
          sprintf('-collapse false -mode %s -norm meandiv -nperm %d', summarize, permutations),
          # Does not matter, the rnk file comes collapsed in advanced
          sprintf('-rnk %s',rnk),
          sprintf('-scoring_scheme %s', score),
          sprintf('-rpt_label %s',label),
          sprintf('-chip %s',chip),
          sprintf('-include_only_symbols true -make_sets true -plot_top_x %d -rnd_seed %d',numplots,seed),
          sprintf('-set_max %d',maxSize),
          sprintf('-set_min %d',minSize),
          sprintf('-zip_report true -out %s -gui false',outdir),sep=' ')

cat(paste("Running", rungsea))

  res <<- system(rungsea,wait=TRUE,intern=TRUE)
  # system(sprintf('ssh pipeline rm %s',rnk),wait=TRUE,intern=TRUE)
  # And this is how an example call should look like
  #res <<- system("ssh pipeline java -cp /Volumes/biostats/soft/BroadGSEA/gsea2-2.0.12.jar -Xmx512m xtools.gsea.GseaPreranked -gmx /Volumes/biostats/databases/BroadGSEA/genesets/dmelanogaster/KEGG.gmt -collapse false -mode Median_of_probes -norm meandiv -nperm 100 -rnk /Volumes/biostats/consulting/cayetano_gonzalez/jpetrovic_array/reports/tables/fc1.rnk -scoring_scheme weighted -rpt_label my_analysis -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 50 -rnd_seed 149 -set_max 500 -set_min 15 -zip_report true -out /Volumes/biostats/soft/BroadGSEA/GenePattern -gui false",wait=TRUE,intern=TRUE)
  }


#################################################################################################################   # Extract genesets that pass a FDR threshold

extract.gsets <- function(fs, th, zone='both')
{
    ff <- list.files(fs);
    if (zone == "top")
    {   
        fpos <- ff[regexpr("gsea_report_for_", ff) > 0 & regexpr("_pos_", ff) > 0 & regexpr("\\.xls$", ff) > 0];
        dp <- read.table(file=paste(fs, "/", fpos, sep=''), sep='\t', header=T, as.is=T, quote='');
        cn <- unlist(read.table(file=paste(fs, "/", fpos, sep=''), sep='\t', header=F, as.is=T, nrows=1));
        colnames(dp) <- cn;
        
    }else if (zone == 'bottom')
    {
        fneg  <- ff[regexpr("gsea_report_for_", ff) > 0 & regexpr("_neg_", ff) > 0 & regexpr("\\.xls$", ff) > 0];
        dp <- read.table(file=paste(fs, "/", fneg, sep=''), sep='\t', header=T, as.is=T, quote='');
        cn <- unlist(read.table(file=paste(fs, "/", fneg, sep=''), sep='\t', header=F, as.is=T, nrows=1));
        colnames(dp) <- cn;
        
    }else if (zone == 'both')
    {
        fpos <- ff[regexpr("gsea_report_for_", ff) > 0 & regexpr("_pos_", ff) > 0 & regexpr("\\.xls$", ff) > 0];
        fneg <- ff[regexpr("gsea_report_for_", ff) > 0 & regexpr("_neg_", ff) > 0 & regexpr("\\.xls$", ff) > 0];
        dpos <- read.table(file=paste(fs, "/", fpos, sep=''), sep='\t', header=T, as.is=T, quote='');
        dneg <- read.table(file=paste(fs, "/", fneg, sep=''), sep='\t', header=T, as.is=T, quote='');
        dp <- rbind(dpos, dneg);
        cn <- unlist(read.table(file=paste(fs, "/", fpos, sep=''), sep='\t', header=F, as.is=T, nrows=1));
        colnames(dp) <- cn;
    }

    dp[dp[, "FDR q-val"] < th, 1];
}

#################################################################################################################
   # Function to execute leading edge analysis within R

runLeadingEdge.mod <- function(gseadir, gsetnames, label='my_analysis',
                               jarpath="/Volumes/biostats/soft/BroadGSEA/gsea2-2.0.12.jar", outdir=getwd(),
                               mar.gsets=35, mar.genes=10)
{

### Required packages

    require("cluster");
    require("gplots");    
    require("hwriter");

### Directory for results

    dir.create(outdir, F)
    fdir <- paste(outdir, "/", label, "/", sep='');
    dir.create(fdir, F);
    ngs <- 0;

    
### If number of gene sets < 2

    if (length(gsetnames) < 2 | is.null(gsetnames))
    {
        warning("Warning! Number of gene sets < 2!");
        sink(paste(fdir, "Warning_no_gene_sets.txt", sep=''));
        cat("Warning! No gene sets given as input!");
        sink();

        p <- openPage(paste(fdir, "leading_edge_index.html", sep=''));
        hwrite(paste("<br><h1>LEADING EDGE ANALYSIS not performed", "</h1>", sep=''), p);
        hwrite(paste("<br><h3>Less than 2 gene sets selected", "</h3>", sep=''), p);        
        closePage(p);        

    }else
    {
    
    ### Extract results for selected gene sets

        ff <- list.files(gseadir);
        fpos <- ff[regexpr("gsea_report_for_", ff) > 0 & regexpr("_pos_", ff) > 0 & regexpr("\\.xls$", ff) > 0];
        fneg <- ff[regexpr("gsea_report_for_", ff) > 0 & regexpr("_neg_", ff) > 0 & regexpr("\\.xls$", ff) > 0];
        dpos <- read.table(file=paste(gseadir, "/", fpos, sep=''), sep='\t', header=T, as.is=T, quote='');
        dneg <- read.table(file=paste(gseadir, "/", fneg, sep=''), sep='\t', header=T, as.is=T, quote='');
        if (nrow(dpos) > 0)
        {
            dpos$"TOP-BOTTOM" <- "Top";
        }
        if (nrow(dneg) > 0)
        {
            dneg$"TOP-BOTTOM" <- "Bottom";
        }
        if (nrow(dpos) > 0 & nrow(dneg) > 0)
        {
            dp <- rbind(dpos, dneg);
            
        }else if (nrow(dpos)>0 & nrow(dneg) == 0)
        {
            dp <- dpos;
            
        }else if (nrow(dpos)==0 & nrow(dneg) > 0)
        {
            dp <- dneg;
            
        }else if (nrow(dpos)==0 & nrow(dneg) == 0)
        {
            warning("Warning! No gene set selected");
            ngs <- 1;
        }

        if (ngs == 0)
        {
            
            ### GSEA results for selected genesets

            rownames(dp) <- dp[, 1];
            cn <- unlist(read.table(file=paste(gseadir, "/", fpos, sep=''), sep='\t', header=F,
                                    as.is=T, nrows=1));
            colnames(dp)[1:length(cn)] <- cn;
            dp$"NA" <- NULL;
            dgs <- dp[dp[, 1]%in%gsetnames, ];
            rownames(dgs) <- dgs[, 1];        
            write.table(dgs, file=paste(fdir, "GSEA_results.txt", sep=''), sep='\t', row.names=F, quote=F);
            

            ### Gene sets used for GSEA analysis

            d <- read.table(file=paste(gseadir, "/edb/gene_sets.gmt", sep=''), sep='\t', header=F,
                            as.is=T, fill=T,
                            quote='', comment.char='~', col.names=paste("V", 1:5000, sep=''));
            rownames(d) <- d[, 1];
            d <- d[, -1];
            d <- d[dgs[, 1], ];
            d <- d[, apply(d, 2, function(o) any(!is.na(o) & o!=''))];

            
            ### Ranked list of genes

            fr <- list.files(paste(gseadir, "/edb/", sep=''));
            frnk <- fr[regexpr("\\.rnk", fr) > 0];
            dr <- read.table(file=paste(gseadir, "/edb/", frnk, sep=''), sep='\t', header=F, as.is=T );
            dr <- dr[order(dr[, 2], decreasing=T), ];
            rownames(dr) <- dr[, 1];

            
            ### Leading edge genes in each geneset       

            ld <- apply(d, 1, function(o) as.character(o[!is.na(o) & o!='']));
            if (is.matrix(ld) && ncol(ld) > 1) ld <- sapply(1:ncol(ld), function(j, ld) ld[, j], ld, simplify=F)
            else if (is.matrix(ld) && ncol(ld) == 1)
            {
                ld <- list(ld[, 1]);
            }
            names(ld) <- rownames(d);
            
            lle <- list();
            for (j in 1:length(ld))
            {
                gs <- names(ld)[j];
                g <- dr[dr[, 1]%in%ld[[j]], 1];
                
                tag <- as.numeric(gsub("%", "", gsub("tags=", "",
                                                     strsplit(dp[gs, "LEADING EDGE"], split=", ")[[1]][1])))/100;
                topb <- dgs[gs, "TOP-BOTTOM"];
                if (topb == 'Top')
                {
                    g <- g[1:round(length(g)*tag)];
                    
                }else if (topb == 'Bottom')
                {
                    g <- rev(g)[1:round(length(g)*tag)];
                }
                lle[[length(lle) + 1]] <- g;
            }
            names(lle) <- names(ld);
            ld <- lle;


            ### Matrix indicating to which gene sets belong each leading edge gene
            
            gsu <- unique(unlist(ld));
            gsu <- gsu[gsu!='' & !is.na(gsu)];

            dggs<- t(sapply(gsu, function(g)
                        {
                            sapply(ld, function(l, g)
                               {
                                   as.numeric(g%in%l);
                               }, g)
                        }));

            clg <- mona(dggs);
            clgs <- mona(t(dggs));    
            dggs <- dggs[clg$order, clgs$order];
            gind <- 1:nrow(dggs);
            sgind <- apply(dggs, 2, function(o) mean(gind[o==1]));
            dggs <- dggs[, order(sgind)];
            
            
            ### Set to set plot
            
            comb <- combn(colnames(dggs), 2);
            dst <- apply(comb, 2, function(o)
                     {
                         lds <- ld[o];
                         length(intersect(lds[[1]], lds[[2]]))/length(union(lds[[1]], lds[[2]]))
                     });
            dhm <- matrix(NA, ncol=length(gsetnames), nrow=length(gsetnames));
            rownames(dhm) <- colnames(dhm) <- colnames(dggs);
            for (j in 1:length(dst))
            {
                dhm[comb[1, j], comb[2, j]] <- dst[j];
                dhm[comb[2, j], comb[1, j]] <- dst[j];                
            }
            diag(dhm) <- 1;

            bk <- seq(0, 1, 0.1);
            colp <- colorpanel(n=length(bk)-1, low="gray90", mid="yellow", high='red3');        
            pdf(paste(fdir, "set_to_set_plot.pdf", sep=''), width=18, height=16);
            heatmap.2(x=dhm,
                      Rowv=NA,
                      Colv=NA,
                      dendrogram='none',
                      scale='none',
                      breaks=bk,
                      col=colp, 
                      main="Set to set plot (Jacquard coefficients)", cex.main=2,
                      key=T, keysize=0.2, density.info='none', lhei=c(1, 5), lwid=c(1, 3),
                      mar=c(mar.gsets, mar.gsets),
                      cexRow = 0.2 + 1/log10(nrow(dhm))/2, cexCol = 0.2 + 1/log10(ncol(dhm))/2,
                      trace='none', na.color='gray50',
                      );
            dev.off();

            dhmc <- cbind(rownames(dhm), dhm);
            write.table(dhmc, file=paste(fdir, "Jacquard_indexes.txt", sep=''), row.names=F, quote=F, sep='\t');
            
            #        heatmap.gsea(t(dhm), bk=seq(0, 1, 0.1), fn=paste(fdir, "set_to_set_plot", sep=''),
            #                     width=1350, height=1200, 
            #                     #vals=format(t(dhm), format='f', dig=2),
            #                     mar=c(mar.gsets, mar.gsets, 5, 5),
            #                     main='Set to set plot (Jacquard coefficients)',  type='square');

            
            ### Leading edge heatmap clustered -> redone

            dcl <- t(dggs);
            dcl[dcl == 0] <- NA;
            lfs <- list.files(paste(gseadir, "/edb", sep=''));
            fs <- lfs[regexpr("\\.rnk$", lfs) > 0];
            dval <-read.table(file=paste(gseadir, "/edb/", fs, sep=''), sep='\t', header=F, as.is=T, quote='');
            rownames(dval) <- dval[, 1];

            genes <- colnames(dcl);
            for(g in genes)
            {
                gv <- dval[g, 2];
                dcl[dcl[, g] == 1, g] <- gv
            }

            
            bk <- c(min(dval[, 2]), seq(quantile(dval[, 2], 0.05), quantile(dval[, 2], 0.95), length=20),
                    max(dval[, 2]));
            sm <- summary(factor(dp$"TOP-BOTTOM", levels=c("Top", "Bottom")));
            if (sm["Bottom"] == 0)
            {
                colp <- colorpanel(n=length(bk)-1, low="ivory", mid="yellow", high='red3');
                
            }else if (sm["Top"] == 0)
            {

                colp <- colorpanel(n=length(bk)-1, low="gray90", mid="violet", high='blue');            
                
            }else
        {
            colp <- colorpanel(n=length(bk)-1, low="blue", mid="gray90", high='red3');
        }
            

            pdf(paste(fdir, "heatmap_v1.pdf", sep=''), width=18, height=9);
            heatmap.2(x=dcl,
                      Rowv=NA,
                      Colv=NA,
                      dendrogram='none',
                      scale='none',
                      breaks=bk,
                      col=colp, 
                      main='Leading edge heatmap', cex.main=2,
                      key=T, keysize=0.2, density.info='none', lhei=c(1, 5), lwid=c(1, 2),
                      mar=c(mar.genes, mar.gsets),
                      cexRow = 0.2 + 1/log10(nrow(dcl))/2, cexCol = 0.2 + 1/log10(ncol(dcl))/2,
                      trace='none', na.color='gray50',
                      );
            dev.off();


            pdf(paste(fdir, "heatmap_v2.pdf", sep=''), width=18, height=9);
            heatmap.2(x=t(dcl),
                      Rowv=NA,
                      Colv=NA,
                      dendrogram='none',
                      scale='none',
                      breaks=bk,
                      col=colp, 
                      main='Leading edge heatmap', cex.main=2,
                      key=T, keysize=0.2, density.info='none', lhei=c(1, 5), lwid=c(1, 2),
                      mar=c(mar.gsets, mar.genes),
                      cexRow = 0.2 + 1/log10(nrow(t(dcl)))/2, cexCol = 0.2 + 1/log10(ncol(t(dcl)))/2,
                      trace='none', na.color='gray50',
                      );
            dev.off();


            #        heatmap.gsea(dcl, bk=bk, fn=paste(fdir, "heatmap_v1", sep=''), col=colp, 
            #                  mar=c(mar.gsets, mar.genes, 5, 5),
            #                  main='Leading edge heatmap', width=nrow(dcl)*200, height=ncol(dcl)*100,
            #                  type='portrait');
            #        heatmap.gsea(t(dcl), bk=bk, fn=paste(fdir, "heatmap_v2", sep=''), col=colp, 
            #                  mar=c(mar.genes, mar.gsets, 5, 5),      
            #                  main='Leading edge heatmap', height=nrow(dcl)*200, width=ncol(dcl)*100,
            #                  type='landscape');
            

            ### Gene in gene sets plot
            
            gu <- unlist(ld);
            tab <- sort(table(gu), decreasing=T);
            
            pdf(paste(fdir, "/gene_in_genesets.pdf", sep=''), width=20,  height=8);
            par(mar=c(5, mar.genes, 5, 5));
            plot(tab, type='h', axes=F, col='red3', ylim=c(0, max(tab)), lwd=2, 
                 main='Gene in subsets', cex.main=1.5, xlab="", ylab='Number of genesets', cex.lab=1.2);
            axis(1, at=1:length(tab), labels=names(tab), las=2);
            axis(2, at=1:length(tab));
            box();
            dev.off();

            dtab <- cbind(Gene=names(tab), "N.gene.sets"=as.numeric(tab));
            write.table(dtab, file=paste(fdir, "gene_in_genesets.txt", sep=''), sep='\t',
                        row.names=F, quote=F
                        );

            ### Jacquard coefficients plot
            
            dstt <- dst;
            dstt[dstt==0] <- 0.0001;
            pdf(paste(fdir, "/jacquard_freq.pdf", sep=''), width=10, height=8);
            plot(sort(dstt, decreasing=T), type='h', col='blue', lwd=3, axes=F, 
                 main='Ordered Jacquard coefficients', cex.main=1.5,
                 xlab="Indexed pairs of gene sets", ylab='Jacquard coefficients',
                 cex.lab=1.2);
            axis(1, at=1:length(dstt), labels=rep("", length(dstt)));
            axis(2);
            box();
            dev.off()


      ### Copy files needed from GSEA results

            gs <- rownames(dcl);
            gs2 <- gsub(" ", "_", gs);
            lfs <- list.files(gseadir);
            lfs <- lfs[unlist(sapply(1:length(gs), function(j, gs, gs2)
                                     {
                                         which(regexpr(gs[j], lfs) > 0 | regexpr(gs2[j], lfs) > 0);
                                     }, gs, gs2))];
            for (f in lfs) file.copy(paste(gseadir, "/", f, sep=''), paste(fdir, f, sep=''), overwrite=T)
            file.copy(paste(gseadir, "/xtools.css", sep=''), paste(fdir, "xtools.css", sep=''), overwrite=T);
            
      ### Table with leading-edge genes in each gene set
            
            dle <- data.frame(Genes=colnames(dcl), t(dcl), stringsAsFactors=F);
            write.table(dle, file=paste(fdir, "Leading_edge_genes_genesets.txt", sep=''), sep='\t',
                        row.names=F, quote=F);
            bgc <- apply(dle[, -1], 2, function(o, bk, colp)
                     {
                         r <- cut(o, breaks=bk, levels=colp, include.lowest=T, right=F);
                         r <- factor(r, level=levels(r), labels=colp);
                         
                     }, bk, colp)

            cl <- col2rgb("white");
            white.rgb <- rgb(cl[1], cl[2], cl[3], maxColorValue=255);
            bgc[is.na(bgc)] <- white.rgb
            bgc <- cbind(rep(white.rgb, nrow(bgc)), bgc);
            colnames(bgc) <- colnames(dle);
            rownames(bgc) <- rownames(dle);
            dle[,-1] <- apply(dle[, -1], 2, function(o)
                          {
                              r <- formatC(o, format='f', dig=3);
                              r[regexpr("NA", r) > 0] <- '';
                              r;
                          });
            dle <- rbind(c("", rownames(dcl)), dle);
            bgc <- rbind(rep(white.rgb, ncol(bgc)), bgc);

            col.links <- c(NA, paste(rownames(dcl), ".html", sep=''));
            col.links[!col.links%in%list.files(gseadir)] <- NA;
            col.links <- sapply(col.links, function(o, n) c(o, rep(NA, n-1)), nrow(dle), simplify=F);
            col.links[[1]] <- paste("http://www.ncbi.nlm.nih.gov/gene/?term=", c("", colnames(dcl)), sep='');
            names(col.links) <- colnames(dle) <- paste("v", 1:ncol(dle))
            
            p <- openPage(paste(fdir, "Leading_edge_genes_genesets.html", sep=''));
            hwrite('<br><h1>LEADING-EDGE GENES BY GENE SETS</h1><br><br>', p);
            hwrite(dle, page=p, center=F, row.names=F, col.names=F, col.width=rep('5px', ncol(dle)),
                   col.style=rep(c("text-align:center"), ncol(dle)), bgcolor=bgc,
                   col.link=col.links);
            closePage(p);

            
            ### List of gene sets included in the analysis


            dlgs <- rbind(colnames(dgs), dgs);
            dlgs <- dlgs[, -c(2, 3)];
            col.links2 <- sapply(1:ncol(dlgs), function(j, n) rep(NA, n), ncol(dlgs), simplify=F);
            col.links2[[1]] <- c(NA, paste(dgs[, 1], ".html", sep=''));
            col.links2[[1]][!col.links2[[1]]%in%list.files(gseadir)] <- NA;
            names(col.links2) <- colnames(dlgs) <- paste("v", 1:ncol(dlgs), sep='');        

            p <- openPage(paste(fdir, "Gene_sets.html", sep=''));
            hwrite('<br><h2>SELECTED GENE SETS FOR LEADING-EDGE ANALYSIS</h2><br><br>', p);
            hwrite(dlgs, page=p, center=F, row.names=F, col.names=F, col.width='500px',
                   col.style=c("text-align:left", rep("text-align:center", ncol(dlgs)-1)),
                   row.style="text-align:center",
                   col.link=col.links2);
            closePage(p);


   ### HTML index
            
            p <- openPage(paste(fdir, "leading_edge_index.html", sep=''));
            hwrite(paste("<br><h1>LEADING EDGE ANALYSIS", "</h1>", sep=''), p);
            hwrite(paste("<h2>Analysis label: ", label, "<br>", "</h2><br><br>", sep=''), p);
            hwrite(paste("<h2>- Selected genes for leading-edge analysis: ",
                         "<A HREF='GSEA_results.txt' ",
                         "download='", label, "_Leading_edge_genes_genesets.txt'>plain text</A><br>",
                         "</h2>", sep=''), p);
            hwrite(dlgs, page=p, center=F, row.names=F, col.names=F, col.width='500px',
                   col.style=c("text-align:left", rep("text-align:center", ncol(dlgs)-1)),
                   row.style="text-align:center",
                   col.link=col.links2);
            
            hwrite("<br><br>", p);
            hwrite(paste("<h2>- Number of gene sets for each gene: <A HREF='gene_in_genesets.pdf'>view</A>",
                         "<A> | </A><A HREF='gene_in_genesets.txt' ",
                         "download='", label, "_gene_in_genesets.txt'>",
                         "plain text</A>","<br>",
                         "</h2>", sep=''), p);
            hwrite(paste("<h2>- Jacquard coefficients between gene sets: ",
                         "<A HREF='set_to_set_plot.pdf'>heatmap</A>",
                         "<A> | </A><A HREF='jacquard_freq.pdf'>barplot</A>",
                         "<A> | </A><A HREF='Jacquard_indexes.txt' ",
                         "download='", label, "_Jacquard_indexes.txt'>",
                         "plain text</A><br>",                     
                         "</h2>", sep=''), p);
            hwrite(paste("<h2>- Genes in each geneset (plot):  ",
                         "<A HREF='heatmap_v1.pdf'>genes vs gene sets</A><A> | </A>",
                         "<A HREF='heatmap_v2.pdf'> gene sets vs genes</A><br>",
                         "</h2>", sep=''), p);                        
            hwrite(paste("<h2>- Genes in each geneset (table):  ",
                         "<A HREF='Leading_edge_genes_genesets.html'>view</A>",
                         "<A> | </A><A HREF='Leading_edge_genes_genesets.txt' ",
                         "download='", label, "_Leading_edge_genes_genesets.txt'>plain text</A>",
                         "<A> | </A><A HREF='Leading_edge_genes_genesets.html' ",
                         "download='", label, "_Leading_edge_genes_genesets.xls'>xls like</A>",
                         "</h2>", sep=''), p);                        
            closePage(p);


        }
    }
}

#################################################################################################################
   # Simple like heatmap image (without clustering, only translation of number to colours)

heatmap.gsea <- function(mat, fn=paste(.wd, "heatmap.gsea", sep=''),
                         bk=NULL, colp=NULL, xlab='', ylab='', main='', #vals=NULL,
                         mar=c(35, 10, 5, 5), 
                         width=1000, height=1000, type='square')
{

    mat <- t(apply(mat, 1, function(o) rev(o)));
#    if (!is.null(vals)) vals <- t(apply(vals, 1, function(o) rev(o)));

    x <- 1:nrow(mat);
    y <- 1:ncol(mat);

    if (is.null(colp))
    {
        colp <- colorpanel(n=length(bk)-1, low="gray90", mid="yellow", high='red3');
    }

    pdf(paste(fn, ".pdf", sep=''), width=width/100, height=height/100);
    
    if (type == 'portrait')
    {
        layout(rbind(c(1, 2), c(1, 3)), widths = c(5, 1), heights = c(1, 5));
        
    }else if (type == 'landscape')
    {
        layout(rbind(c(1, 2), c(1, 3)), widths = c(25, 1), heights = c(2, 1));
        
    }else if (type =='square')
    {
        layout(rbind(c(1, 2), c(1, 3)), widths = c(8, 1), heights = c(1, 3));
    }
    
    par(mar=mar);
    image(x=x, y=y, z=mat, las=2, axes=F,
          breaks=bk, col=colp, xlab=xlab, ylab=ylab, main=main,
          cex.main=1.2);
    cex.row <- 1;
    cex.col <- 1;
    mtext(rownames(mat), 1, .1, at=x, las=2, cex=cex.col);
    mtext(colnames(mat), 2, .2, at=y, las=2, cex=cex.row);
    segments(c(x[1]-1, x)+.5, min(y)-.5,
             c(x[1]-1, x)+.5, min(y)+c(max(y), max(y))-.5, xpd=TRUE)
    segments(min(x)-.5, c(y[1]-1, y)+.5,  min(x)+c(max(x), max(x))-.5, c(y[1]-1, y)+.5, xpd=TRUE)
#    if (!is.null(vals))
#    {
#        vals[is.na(vals) | regexpr("NA$", vals) > 0] <- '';
#        text(rep(x, ncol(mat)), rep(y, each=nrow(mat)), vals, cex=1/max(c(ncol(mat)/nrow(mat))));
#    }

    z <-  t(bk[-length(bk)]+diff(bk)/2)
    xl <- 1:nrow(z);
    yl <- 1:ncol(z);
    par(mar=c(1,3,1,1));    
    image(x=xl, y=yl, z=z, las=2, axes=F, breaks=bk, col=colp, xlab=xlab, ylab=ylab, main='');
    segments(0.6, c(yl[1]-1, yl)+.5,
                 1.4, c(yl[1]-1, yl)+.5, xpd=TRUE)
    segments(c(xl[1]-0.8, xl)+.4,  min(yl)-.5, c(xl[1]-0.8, xl)+.4, min(yl)+c(max(yl), max(yl))-.5, xpd=TRUE)
    mtext(formatC(bk, format='f', dig=2), 2, .2, at=c(0, yl)+0.5, cex=max(c(cex.col, cex.row))*0.8, las=2);

    plot(1, pch='', xlab='', ylab='', axes=F);
    
    dev.off();
}



#################################################################################################################
#################################################################################################################
#################################################################################################################
