##########################################################################
## PROJECT: BDURAN RIPSEQ
## DESCRIPTION OF FILE: GENES DIFFERENTIAL BINDING ANALYSIS VIA DESEQ2
## USE SAMPLES FROM SECOND BATCH
## DATE: NOV 2019
## NOTES: RUN IN OSX CATALINA
###########################################################################

library(GenomicRanges)
library(Rsamtools)
library(DESeq2)
library(GenomicAlignments)
library(DESeq2)
library(Rsubread)
library(openxlsx)
library(R2HTML)
library(vioplot)
library(pheatmap)
library(pvclust)

# Compute counts or RPKM over samples
getCounts <- function(peaks,x,rpkm=FALSE,dolog=FALSE)
    {
        mnames <- intersect(names(peaks),names(x))
        peaks <- peaks[mnames]
        ans <- countOverlaps(query=peaks,subject=x[names(peaks)])
        if (rpkm)
            {
                nreads <- nrow(x)
                ans <- 10^9*unlist(ans)/(nreads*(end(peaks)-start(peaks)+1))
                if (dolog) ans <- log2(ans+1)
                return(ans)
            }
        else return(unlist(ans))
    }

###########################################################################

datadir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/'
datadir2 <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/'
datadir1 <- '/Volumes/biostats/consulting/raul_mendez/bduran_201810_ripseq/data/'
reportsdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports'
tablesdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables'
tdfdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables/files_4_igv'
figsdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/figs'

###########################################################################
# 1. DESEQ2, RIP vs INPUT. RAW READS. GENE LEVEL (MULTIOVERLAP=TRUE)
###########################################################################

## # Extract GENESs
## xl92 <- read.delim('/Volumes/biostats/databases/RefGenome/xlaevis/xlaevis_92_ucsc/XENLA_9.2_Xenbase.gff3',comment.char='#',header=FALSE,as.is=TRUE)

## # Build annotation
## # Use makeTxDbFromGFF function from GenomicFeatures package
## library(GenomicFeatures)
## xl92 <- makeTxDbFromGFF('/Volumes/biostats/databases/RefGenome/xlaevis/xlaevis_92_ucsc/XENLA_9.2_Xenbase.gff3',organism='Xenopus laevis')
## xl92.genes <- genes(xl92)
## xl92.exons <- exons(xl92)
## xl92.cds <- cds(xl92)
## xl92.transcripts <- transcripts(xl92)
## xl92.5utr <- fiveUTRsByTranscript(xl92)
## xl92.genes <- threeUTRsByTranscript(xl92)
## xl92.introns <- intronsByTranscript(xl92)
## xl92.txdb <- list(Genes=xl92.genes,Exons=xl92.exons,CDS=xl92.cds,Transcripts=xl92.transcripts,UTR5=xl92.5utr,GENES=xl92.genes,Introns=xl92.introns)
## save(xl92.txdb,file=file.path(datadir,'xl92.txdb.GRanges.annot.RData'))

# Use precomputed gene information
load(file=file.path(datadir,'xl92.txdb.GRanges.annot.RData'))
genes <- as.data.frame(xl92.txdb$Genes)
genes <- genes[,c('gene_id','seqnames','start','end','strand')]
colnames(genes) <- c('GeneID','Chr','Start','End','Strand')

# Load precomputed Count reads over GENES, both batches
## bamfiles1 <- list.files(file.path(datadir1,'sequences/bam'),pattern='*.fastq.bam$',full.names=TRUE)
## names(bamfiles1) <- substr(basename(bamfiles1),1,17)
## bamfiles2 <- list.files(file.path(datadir2,'sequences/bam'),pattern='*.fastq.bam$',full.names=TRUE)
## names(bamfiles2) <- substr(basename(bamfiles2),1,17)
## bamfiles <- c(bamfiles1,bamfiles2)
## fcounts.mapq1 <- featureCounts(bamfiles,annot.ext=genes,isGTFAnnotationFile=FALSE,allowMultiOverlap=TRUE,countMultiMappingReads=FALSE,minMQS=1,nthreads=16)
## save(fcounts.mapq1,file=file.path(datadir,'featureCounts_bowtie2_genes_allbatches_multioverlap_minq1_notmultimap.RData'))
     
# Read sampleinfo, prepare for compare IP vs Input, join both batches
load(file=file.path(datadir,'featureCounts_bowtie2_genes_allbatches_multioverlap_minq1_notmultimap.RData'))
load(file=file.path(datadir1,'sampleinfo.RData'))
sampleinfo1 <- sampleinfo
load(file=file.path(datadir,'sampleinfo_forDESeq2_may19.RData'))
sampleinfo2 <- sampleinfo

# Load premade sampleinfo for Join datasets
## table(sampleinfo2$Group,sampleinfo2$Replicate)
## sampleinfo1 <- sampleinfo1[-c(grep('+',sampleinfo1$Sample.Name,fixed=TRUE)),]
## sampleinfo1$Sample.Name <- gsub(' ','_',sampleinfo1$Sample.Name)
## sampleinfo1$Type <- 'IP'; sampleinfo1$Type[grep('Input',sampleinfo1$Sample.Name)] <- 'Input'
## sampleinfo1$Factor <- c('NI','CPEB1','CPEB3','NI','CPEB1','CPEB3')
## sampleinfo1$Replicate <- 0
## sampleinfo1$Sample.Name <- gsub('Elu','IP',sampleinfo1$Sample.Name)
## sampleinfo1$Name <- sampleinfo1$Group <- unlist(lapply(strsplit(sampleinfo1$Sample.Name,'_'),function(x) paste(x[2],x[1],sep='_')))
## sampleinfo1$Batch <- 1
## rownames(sampleinfo1) <- sampleinfo1$Name <- paste(sampleinfo1$Group,sampleinfo1$Replicate,sep='_')
## sampleinfo2$Batch <- 2
## sampleinfo <- rbind(sampleinfo1,sampleinfo2)
## sampleinfo$Group <- factor(sampleinfo$Group)
## sampleinfo$Replicate <- factor(sampleinfo$Replicate)
## sampleinfo$Batch <- factor(sampleinfo$Batch)
load(file=file.path(datadir,'sampleinfo_allbatches_forDESeq2_may19.RData'))

# Subset for CPEBs from Batch 2
sampleinfo.bak <- sampleinfo
sampleinfo <- sampleinfo[sampleinfo$Factor!='NI' & sampleinfo$Type!='Input' & sampleinfo$Replicate %in% c(3,4),]

# Build DESeq dataset, paired by replicate, use counts mapq1
de.counts <- fcounts.mapq1$counts
colnames(de.counts) <- gsub('X.Volumes.biostats.consulting.raul_mendez.bduran_201810_ripseq.data..sequences.bam.','',gsub('.fastq.bam','',colnames(de.counts)))
colnames(de.counts) <- gsub('X.Volumes.biostats.consulting.raul_mendez.bduran_201903_ripseq.data..sequences.bam.','',gsub('.fastq.bam','',colnames(de.counts)))
de.counts <- de.counts[,sampleinfo$ID]
colnames(de.counts) <- sampleinfo$Name

# Add CPEB1 vs CPEB234
sampleinfo$CPEBGroup <- ifelse(sampleinfo$Factor=='CPEB1','CPEB1','CPEB234')
dds <- DESeqDataSetFromMatrix(countData=de.counts,colData=sampleinfo,design= ~ Replicate + CPEBGroup)
dds <- DESeq(dds)
summary(results(dds))
save(dds,file=file.path(datadir,'bowtie2_genes_allbatches_dds_DESeq2_diffbind_fcounts_minq1_nov19.RData'))

# Plot PCA
rld <- rlog(dds)

pdf(file.path(figsdir,'bowtie2_genes_allbatches_pca_DESeq2_diffbind_fcounts_minq1_nov19.pdf'))
plotPCA(rld,intgroup='Type')
plotPCA(rld,intgroup='Factor')
plotPCA(rld,intgroup='Replicate')
plotPCA(rld,intgroup='Group')
plotPCA(rld,intgroup='CPEBGroup')
plotPCA(rld,intgroup='Batch')
plotPCA(rld,intgroup=c('Batch','Group'))
dev.off()

# Extract results for contrasts of interest
conts <- list(CPEB1_vs_CPEB234=c('CPEBGroup','CPEB1','CPEB234'))

# Output, shrink and non-shrink
reslist <- lapply(conts,function(x) results(dds,pAdjustMethod='BH',cooksCutoff=FALSE,contrast=x))
reslist.adj <- mclapply(conts,function(x) lfcShrink(dds,contrast=x,pAdjustMethod='BH',cooksCutoff=FALSE),mc.cores=length(conts))

# Add rej column
reslist.adj <- lapply(reslist.adj,function(x) { x$rej <- sign(x$log2FoldChange) * (x$padj<0.05) * (abs(x$log2FoldChange)>log2(2)); x })
library(gtools)

# Import computed foldchanges for IP vs Input and IP vs IgG
# Files
xout <- list.files(tablesdir,'bowtie2_genes_allbatches',full.names=TRUE)
xout <- xout[grep('_DESeq2_adj_fcounts_minq1_may19.csv',basename(xout))]
names(xout) <- gsub('bowtie2_genes_allbatches_','',gsub('_DESeq2_adj_fcounts_minq1_may19.csv','',basename(xout)))
# Read
xout <- lapply(xout,read.csv)
xout <- lapply(xout,function(x) { rownames(x) <- x$GeneID; x })
allgenes.rejany <- sort(unique(unlist(lapply(xout,function(x) { sel <- !is.na(x$rej) & x$rej==1; rownames(x)[sel] }))))
allgenes.rej1 <- TRUE
# xout normcounts and log2FC
xout.rlog <- xout[[1]][allgenes.rej1,28:33]
xout.log2FC <- as.data.frame(do.call(cbind,lapply(xout,function(x) x[allgenes.rej1,'log2FoldChange'])))
xout.padj <- as.data.frame(do.call(cbind,lapply(xout,function(x) x[allgenes.rej1,'padj'])))
xout.target <- as.data.frame(do.call(cbind,lapply(xout,function(x) x[allgenes.rej1,'rej'])))
xout.target[is.na(xout.target)] <- 0
xout.target <- xout.target[,grep('CPEB',colnames(xout.target))]
rownames(xout.target) <- rownames(xout.rlog) <- rownames(xout.log2FC) <- rownames(xout[[1]])

# Some sanity checks
xx <- xout.log2FC
sel.up <- rownames(xx)[!is.na(reslist.adj[[1]]$rej) & reslist.adj[[1]]$rej==1]
sel.dn <- rownames(xx)[!is.na(reslist.adj[[1]]$rej) & reslist.adj[[1]]$rej==-1]
idx.up <- intersect(sel.up,rownames(xout.log2FC))
idx.dn <- intersect(sel.dn,rownames(xout.log2FC))
idx.up.target <- intersect(sel.up,allgenes.rejany)
idx.dn.target <- intersect(sel.dn,allgenes.rejany)

# Target sanity checks
kk <- xout.target
kk[kk==-1] <- 0
colSums(kk)

# Add exhaustive CPEB1 and CPEB234 targets
idx.up.strict.input <- intersect(sel.up,rownames(xout.target)[xout.target$CPEB1_vs_Input==1])
idx.up.strict.ni <- intersect(sel.up,rownames(xout.target)[xout.target$CPEB1_vs_NI==1])
idx.dn.strict.input <- intersect(sel.dn,rownames(xout.target)[xout.target$CPEB2_vs_Input==1 & xout.target$CPEB3_vs_Input==1 & xout.target$CPEB4_vs_Input==1])
idx.dn.strict.ni <- intersect(sel.dn,rownames(xout.target)[xout.target$CPEB2_vs_NI==1 & xout.target$CPEB3_vs_NI==1 & xout.target$CPEB4_vs_NI==1])

# Plots and Violin plots
pdf(file.path(figsdir,'DESEq2_Diff_CPEB1_vs_CPEB234_Violin_v2.pdf'),width=10,height=10)
# CPEB1 IP vs Input
stripchart(xx[idx.up,c(1,3,5,7)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB1\n(n=%d genes)',length(idx.up)),col='grey',ylab='log2FC IP/Input')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.up,1],xx[idx.up,3],xx[idx.up,5],xx[idx.up,7],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
# CPEB234 IP vs Input
stripchart(xx[idx.dn,c(1,3,5,7)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB234\n(n=%d genes)',length(idx.dn)),col='grey',ylab='log2FC IP/Input')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.dn,1],xx[idx.dn,3],xx[idx.dn,5],xx[idx.dn,7],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
# CPEB1 IP vs IP_NI
stripchart(xx[idx.up,c(2,4,6,8)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB1\n(n=%d genes)',length(idx.up)),col='grey',ylab='log2FC IP/IP_NI')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.up,2],xx[idx.up,4],xx[idx.up,6],xx[idx.up,8],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
# CPEB234 IP vs IP_NI
stripchart(xx[idx.dn,c(2,4,6,8)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB234\n(n=%d genes)',length(idx.dn)),col='grey',ylab='log2FC IP/IP_NI')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.dn,2],xx[idx.dn,4],xx[idx.dn,6],xx[idx.dn,8],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
dev.off()

# Plots and Violin plots (only genes with rej=1 in some CPEB)
pdf(file.path(figsdir,'DESEq2_Diff_CPEB1_vs_CPEB234_Targets_Violin_v2.pdf'),width=10,height=10)
# CPEB1 IP vs Input
stripchart(xx[idx.up.target,c(1,3,5,7)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB1\n(n=%d genes)',length(idx.up.target)),col='grey',ylab='log2FC IP/Input')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.up.target,1],xx[idx.up.target,3],xx[idx.up.target,5],xx[idx.up.target,7],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
# CPEB234 IP vs Input
stripchart(xx[idx.dn.target,c(1,3,5,7)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB234\n(n=%d genes)',length(idx.dn.target)),col='grey',ylab='log2FC IP/Input')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.dn.target,1],xx[idx.dn.target,3],xx[idx.dn.target,5],xx[idx.dn.target,7],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
# CPEB1 IP vs IP_NI
stripchart(xx[idx.up.target,c(2,4,6,8)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB1\n(n=%d genes)',length(idx.up.target)),col='grey',ylab='log2FC IP/IP_NI')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.up.target,2],xx[idx.up.target,4],xx[idx.up.target,6],xx[idx.up.target,8],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
# CPEB234 IP vs IP_NI
stripchart(xx[idx.dn.target,c(2,4,6,8)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB234\n(n=%d genes)',length(idx.dn.target)),col='grey',ylab='log2FC IP/IP_NI')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.dn.target,2],xx[idx.dn.target,4],xx[idx.dn.target,6],xx[idx.dn.target,8],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
dev.off()

# Plots and Violin plots (only genes with rej=1 in some CPEB)
pdf(file.path(figsdir,'DESEq2_Diff_CPEB1_vs_CPEB234_StrictTargets_Violin_v2.pdf'),width=10,height=10)
# CPEB1 IP vs Input
stripchart(xx[idx.up.strict.input,c(1,3,5,7)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB1\n(n=%d genes)',length(idx.up.strict.input)),col='grey',ylab='log2FC IP/Input')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.up.strict.input,1],xx[idx.up.strict.input,3],xx[idx.up.strict.input,5],xx[idx.up.strict.input,7],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
# CPEB234 IP vs Input
stripchart(xx[idx.dn.strict.input,c(1,3,5,7)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB234\n(n=%d genes)',length(idx.dn.strict.input)),col='grey',ylab='log2FC IP/Input')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.dn.strict.input,1],xx[idx.dn.strict.input,3],xx[idx.dn.strict.input,5],xx[idx.dn.strict.input,7],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
# CPEB1 IP vs IP_NI
stripchart(xx[idx.up.strict.ni,c(2,4,6,8)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB1\n(n=%d genes)',length(idx.up.strict.ni)),col='grey',ylab='log2FC IP/IP_NI')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.up.strict.ni,2],xx[idx.up.strict.ni,4],xx[idx.up.strict.ni,6],xx[idx.up.strict.ni,8],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
# CPEB234 IP vs IP_NI
stripchart(xx[idx.dn.strict.ni,c(2,4,6,8)],method='jitter',vertical=TRUE,pch=19,xlim=c(0.5,4.5),ylim=c(-10,10),main=sprintf('CPEB234\n(n=%d genes)',length(idx.dn.strict.ni)),col='grey',ylab='log2FC IP/IP_NI')
grid()
abline(h=0,lty=2)
vioplot(xx[idx.dn.strict.ni,2],xx[idx.dn.strict.ni,4],xx[idx.dn.strict.ni,6],xx[idx.dn.strict.ni,8],names=c('TSS','CDS'),col='#00000020',lwd=1.5,add=TRUE)
dev.off()

# Heatmap sel.up and down
pdf(file.path(figsdir,'DESEq2_Diff_CPEB1_vs_CPEB234_Heatmap.pdf'),width=5,height=10)
pheatmap(xx[c(idx.up,idx.dn),c(1,3,5,7)],scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,clustering_distance_rows='correlation',main='DESeq2 log2FC IP/Input \n(lfcShrink) ')
pheatmap(xx[c(idx.up,idx.dn),c(2,4,6,8)],scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,clustering_distance_rows='correlation',main='DESeq2 log2FC IP/IP_NI \n(lfcShrink) ')
pheatmap(xx[c(idx.up,idx.dn),c(1,3,5,7,2,4,6,8)],scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,main='DESeq2 log2FC IP/Input | IP/IP_NI \n(lfcShrink) ')
dev.off()

# Heatmap sel.up and down AND target
pdf(file.path(figsdir,'DESEq2_Diff_CPEB1_vs_CPEB234_Targets_Heatmap.pdf'),width=5,height=10)
pheatmap(xx[c(idx.up.target,idx.dn.target),c(1,3,5,7)],scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,clustering_distance_rows='correlation',main='DESeq2 log2FC IP/Input \n(lfcShrink) ')
pheatmap(xx[c(idx.up.target,idx.dn.target),c(2,4,6,8)],scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,clustering_distance_rows='correlation',main='DESeq2 log2FC IP/IP_NI \n(lfcShrink) ')
pheatmap(xx[c(idx.up.target,idx.dn.target),c(1,3,5,7,2,4,6,8)],scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,main='DESeq2 log2FC IP/Input | IP/IP_NI \n(lfcShrink) ')
dev.off()

# Heatmap sel.up and down AND strict target
pdf(file.path(figsdir,'DESEq2_Diff_CPEB1_vs_CPEB234_StrictTargets_Heatmap.pdf'),width=5,height=10)
pheatmap(xx[c(idx.up.target.input,idx.dn.target.input),c(1,3,5,7)],scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,clustering_distance_rows='correlation',main='DESeq2 log2FC IP/Input \n(lfcShrink) ')
pheatmap(xx[c(idx.up.target.ni,idx.dn.target.ni),c(2,4,6,8)],scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,clustering_distance_rows='correlation',main='DESeq2 log2FC IP/IP_NI \n(lfcShrink) ')
pheatmap(xx[c(unique(c(idx.up.strict.input,idx.up.strict.ni)),unique(c(idx.dn.strict.input,idx.dn.strict.ni))),c(1,3,5,7,2,4,6,8)],
         scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,main='DESeq2 log2FC IP/Input OR IP/IP_NI \n(lfcShrink) ')
pheatmap(xx[c(intersect(idx.up.strict.input,idx.up.strict.ni),intersect(idx.dn.strict.input,idx.dn.strict.ni)),c(1,3,5,7,2,4,6,8)],
         scale='row',method='ward',cluster_cols=FALSE,show_rownames=FALSE,main='DESeq2 log2FC IP/Input AND IP/IP_NI \n(lfcShrink) ')
dev.off()

# Heatmap of log2FC correlations
sel <- apply(xx,1,function(x) sum(is.na(x))==0)

# Pvclust for support
xx.clus <- xx[sel,]

pvclus.1 <- pvclust(xx[,c(1,3,5,7)],method.hclust='ward',method.dist='correlation',nboot=1000)
pvclus.2 <- pvclust(xx[,c(2,4,6,8)],method.hclust='ward',method.dist='correlation',nboot=1000)
set.seed(149)
samplegenes <- allgenes.rejany[sample(1:length(allgenes.rejany),1000)]

pdf(file.path(figsdir,'DESEq2_Diff_CPEB1_vs_CPEB234_Targets_Clustering.pdf'))
pheatmap(cor(xx[sel,c(1:8)]),main='log2FC IP/Control correlation')
plot(pvclus.1,main='logFC CPEB any target (IP vs Input)\nPvclust cor, Ward.D, B=1000')
pvrect(pvclus.1,alpha=0.95)
plot(pvclus.2,main='log2FC CPEB any target (IP vs IP_NI)\nPvclust cor, Ward.D, B=1000')
pvrect(pvclus.2,alpha=0.95)
dev.off()

# Make final tables
means <- xout[[1]][,29:34]
reslist.df <- reslist.adj[[1]]
colnames(reslist.df) <- paste(colnames(reslist.df),'CPEB1_vs_CPEB234',sep='.')
unique(rownames(means)==rownames(xout.target))
unique(rownames(means)==rownames(xout.log2FC))
unique(rownames(means)==rownames(reslist.df))
#unique(rownames(means)==rownames(xout.padj))

xout.df <- data.frame(xout[[1]][,2:6],
                      NormCounts=means,
                      log2FC=xout.log2FC,
                      padj=xout.padj,
                      rej=xout.target,
                      reslist.df[,c(2,4:7)])

library(openxlsx)
write.xlsx(xout.df,file=file.path(tablesdir,'bowtie2_genes_allbatches_CPEB1_vs_CPEB234_DESeq2_Diff_fcounts_minq1_nov19_v2.xlsx'),keepNA=TRUE)

# Save BED files
for (i in names(xoutlist))
    {
        xx <- as.data.frame(xoutlist[[i]])
        xx <- xx[!is.na(xx$rej),]
        write.table(xx[xx$rej==1,c('Chr','Start','End','log2FoldChange')],file=file.path(tdfdir,sprintf('bowtie2_genes_allbatches_%s_DESeq2_fcounts_minq1_may19.bed',i)),
                    sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE,dec='.')
    }

### THE END
