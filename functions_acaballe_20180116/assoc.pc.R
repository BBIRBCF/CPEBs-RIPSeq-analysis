#################################################################################################################
#################################################################################################################
### assoc.pc.R
#################################################################################################################
#################################################################################################################


assoc.vars <- function(d,  pd, nv, fn, vadd=NULL)
{

### Association of components with variables

    lr <- sapply(nv, function(v, pd)
             {
                 apply(d, 2, function(o, v, pd)
                   {
                     m <- lm(o ~ pd[, v])
                     pv <- drop1(m, test='F')[2, "Pr(>F)"]
                     r2 <- summary(m)$adj.r.squared
                     c(r2, pv)
                   }, v, pd)
             }, pd, simplify=F)

    pvs <- sapply(lr, function(o) o[2, ])
    r2 <- sapply(lr, function(o) o[1, ])

    dr <- NULL
    for (j in 1:ncol(pvs))
    {
      dr <- cbind(dr, paste(formatC(r2[, j], format='f', dig=3), "<br>(", base::format.pval(pvs[, j]), ")", sep=''))
    }
    colnames(dr) <- colnames(pvs)
    rownames(dr) <- rownames(pvs)    

    dr <- cbind(pc=rownames(dr), dr)
    dr <- rbind(pc=colnames(dr), dr)    
    for (j in 1:ncol(dr)) dr[, j] <- as.character(dr[, j]);
    bgc <- pvs
    rgb2 <- function(col)
    {
        cl <- col2rgb(col);
        rgb(cl[1], cl[2], cl[3], maxColorValue=255);
    }
    bgc[, ] <- rgb2("white");
    bgc[pvs < 0.05] <- rgb2("red");
    bgc <- rbind(rep(rgb2("white"), ncol(bgc)), bgc);
    bgc <- cbind(rep(rgb2("white"), nrow(bgc)), bgc);

    if (!is.null(vadd))
    {
      dr <- cbind(dr, c("pcvar", vadd))
      bgc <- cbind(bgc, rep(rgb2("white"), nrow(bgc)))
    }
    
    p <- openPage(paste(fn, '.html', sep=''));
    hwrite('<br><br><br>R2 (F pval)', p);
    hwrite(dr, page=p, center=F, row.names=F, col.names=F, 
           col.width=rep('120px', ncol(dr)),  col.style=rep(c("text-align:center"), ncol(dr)), bgcolor=bgc);
    closePage(p);

    invisible(NULL)
}

assoc.pc <- function(pc, pd, nv, rdir, nlab=NULL, show.labs=F, vadd=NULL, npc=9)
{
  require("hwriter");
    
  ### Create directories in case they don t exist

  rdir <- paste(rdir, "/", sep='')
  dir.create(rdir, F)

  
  ### Show labs-values or points

  if (show.labs)
    {
      pch <- ''
    }else
  {
    pch <- 16
  }

  ### Compute proporcions of variance

  pcvar <- pc$sdev^2/sum(pc$sdev^2)*100
  png(paste(rdir, "pc_barplot.png", sep=''), width=800, height=800)
  barplot(pcvar)
  dev.off()

  sink(paste(rdir, "summary_pc.txt", sep=''));
  print(summary(pc));
  sink();


  ### Association with PCAs

  pcvar <- formatC(pcvar, format='f', dig=1);
  assoc.vars(d=pc$x[, 1:npc],  pd=pd, nv=nv, fn=paste(rdir, "Assoc_pca", sep=''), vadd=paste(pcvar[1:npc], "%", sep=''))
  

### Plots for interest variables

    rgs  <- apply(pc$x, 2, range);
    rgs <- rbind(rgs,  apply(rgs, 2, diff));
    rgs <- apply(rgs, 2, function(o) c(o[1]-o[3]*0.05, o[2]+o[3]*0.05))

    for (v in nv)
    {
        x <- pd[, v];
        if (is.character(x)) x <- factor(x)
        sm <- table(x);
        if (is.factor(x))
        {
            ncols <- c("red3", "blue", "green3", "violet", "orange", "black",
                       "cyan", "brown", "yellow2", "lightblue", "darkolivegreen", "deepskyblue4",
                       "deeppink", "bisque", "yellowgreen")[1:(nlevels(x))];
            cols <- ncols[as.numeric(x)];
            lv <- levels(x);
        }else{
            ncols <- cols <- 'black';
        }
        labs <- pd[, nlab];        
        for (j in 2:npc)
        {
            if (is.factor(x))
            {
                if ((j%%npc) == 2)
                {
                    png(paste(rdir, v, "_pc1_pc", j, "_", j + 7, ".png", sep=''),
                        width=800, height=800);
                    par(mfrow=c(3, 3));
                    plot(1:10, axes=F, pch='', xlab='', ylab='');
                    legend(x='top', fill=ncols, legend=paste(levels(x), " (", sm, ")", sep=''),
                           horiz=F, xpd=T, cex=1.5);
                }
                plot(pc$x[, 1], pc$x[, j], pch=pch, cex.main=1.5, cex.lab=1.2, cex.axis=1.2,
                     main=paste('Principal Components 1 and ', j),
                     xlab=paste("PC1 (", pcvar[1], "%)", sep=''),
                     ylab=paste("PC", j, " (", pcvar[j], "%)", sep=''),
                     xlim=rgs[1:2, 1], ylim=rgs[1:2, j], col=cols);
                if (show.labs)
                {
                    text(pc$x[, 1], pc$x[, j], label=labs, col=cols);
                }
                abline(h=0, v=0, lty=3);
                if ((j%%8) == 1 | j == npc)
                {
                    dev.off();
                }
            }
            else
            
            {
                if ((j%%npc) == 2)
                {
                    png(paste(rdir, v, "_pc", j-1, "_", j + 7, ".png", sep=''),
                        width=800, height=800);
                    par(mfrow=c(3, 3));
                }
                rgv <- range(x, na.rm=T);
                amp <- diff(rgv);
                rgv[1] <- rgv[1] - amp*0.05;
                rgv[2] <- rgv[2] + amp*0.05;
                plot(pd[, v], pc$x[, j-1],
                     xlab=v, ylab=paste("PC", j-1, " (", pcvar[j-1], "%)", sep=''),
                     main=paste("Principal Component", j-1, "and", v),
                     xlim=rgv, ylim=rgs[1:2, j-1],
                     pch=pch, cex.main=1.5, cex.lab=1.2, cex.axis=1.2, col=cols);
                if (show.labs)
                {
                    text(x, pc$x[, j-1], label=labs, col=cols);
                }                

                if ((j%%8) == 1 | j == npc)
                {
                    dev.off();
                }
            }
        }
    }
    NULL;
  
}

#################################################################################################################

varcomp <- function(d,  pd, nv, fn, pcvar=NULL)
{


  ### Components of variance

  vc.pc <- t(apply(d, 2, function(o, pd, nv)
              {
                pd$pcaj <- o
                fm <- paste("pcaj  ~ ", paste("(1|", nv, ")", sep='', collapse=' + '), sep='')
                m <- lmer(as.formula(fm), REML=T, data=pd)
                r <- c(unlist(VarCorr(m))[nv], resid=summary(m)$sigma^2)
                r/sum(r)*100
              }, pd, nv))

  if (!is.null(pcvar))
  {
    vc.pc.w <- vc.pc*pcvar/100
    vc.pc.t <- apply(vc.pc.w, 2, sum)
  }else
  {
    vc.pc.w <- NULL
    vc.pc.t <- NULL
  }

 
  tab <- NULL
  
  for (j in 1:ncol(vc.pc))
  {
    if (!is.null(vc.pc.w))
    {  
      tab <- cbind(tab, paste(paste(formatC(vc.pc[, j], format='f', dig=1), "%", sep=''),
                              paste(formatC(vc.pc.w[, j], format='f', dig=1), "%", sep=''), sep='<br>'))
    }else
    {
      tab <- cbind(tab, paste(formatC(vc.pc[, j], format='f', dig=1), "%", sep=''))
    }
    
  }
  rownames(tab) <- rownames(vc.pc)
  colnames(tab) <- colnames(vc.pc)
  if (!is.null(vc.pc.t))
    tab <- rbind(tab, TOT=paste(formatC(vc.pc.t, format='f', dig=1), "%", sep=''))
  if (!is.null(pcvar))
    tab <- cbind(tab, pcvar=paste(formatC(c(pcvar, sum(pcvar)), format='f', dig=1), "%", sep=''))
  tab <- rbind(colnames(tab), tab)
  tab <- cbind(rownames(tab), tab)
  
  p <- openPage(fn);
  hwrite('<br><br><br>Variance proportion within PC / relative to  total variance<br><br>', p);
  hwrite(tab, page=p, center=F, row.names=F, col.names=F, 
         col.width=rep('120px', ncol(tab)),  col.style=rep(c("text-align:center"), ncol(tab)))
  closePage(p);

  vc.pc


}

#################################################################################################################
#################################################################################################################
#################################################################################################################
