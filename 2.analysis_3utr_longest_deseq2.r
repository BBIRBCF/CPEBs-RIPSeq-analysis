##########################################################################
## PROJECT: BDURAN RIPSEQ
## DESCRIPTION OF FILE: 3UTR DIFFERENTIAL ANALYSIS VIA DESEQ2
## DATE: MAY 2019
###########################################################################

library(GenomicRanges)
library(Rsamtools)
library(parallel)
library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(biomaRt)
library(DESeq2)
library(GenomicAlignments)
library(MDA)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DESeq2)
library(Rsubread)

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
reportsdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports'
tablesdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables'
tdfdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables/files_4_igv'
figsdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/figs'

###########################################################################
# 1. DESEQ2, RIP vs INPUT. RAW READS. ALL UTRS (MULTIOVERLAP=TRUE)
###########################################################################

## # Extract 3UTRs
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
## xl92.3utr <- threeUTRsByTranscript(xl92)
## xl92.introns <- intronsByTranscript(xl92)
## xl92.txdb <- list(Genes=xl92.genes,Exons=xl92.exons,CDS=xl92.cds,Transcripts=xl92.transcripts,UTR5=xl92.5utr,UTR3=xl92.3utr,Introns=xl92.introns)
## save(xl92.txdb,file=file.path(datadir,'xl92.txdb.GRanges.annot.RData'))

# Map 3UTR exons to genes
load(file=file.path(datadir,'xl92.txdb.GRanges.annot.RData'))
idx <- findOverlaps(xl92.txdb$UTR3,xl92.txdb$Genes)
exon2gene <- data.frame(Exon=as.data.frame(xl92.txdb$UTR3)[queryHits(idx),],Gene=as.data.frame(xl92.txdb$Genes)[subjectHits(idx),])
sum(is.na(exon2gene$Gene.gene_id))
exon2gene <- exon2gene[exon2gene$Exon.strand==exon2gene$Gene.strand,]
#exon2gene <- exon2gene[order(exon2gene$Gene.gene_id,exon2gene$Exon.width,decreasing=TRUE),]
#head(exon2gene,20)
#exon2gene <- exon2gene[!duplicated(exon2gene$Gene.gene_id),]
rownames(exon2gene) <- paste0(exon2gene$Gene.gene_id,'.',exon2gene$Exon.exon_name)

# make utr3
utr3 <- exon2gene[,c('Exon.exon_name','Exon.seqnames','Exon.start','Exon.end','Exon.strand')]
colnames(utr3) <- c('GeneID','Chr','Start','End','Strand')
utr3$GeneID <- rownames(utr3)

# Count reads over 3UTR
bamfiles <- list.files(file.path(datadir,'sequences/bam'),pattern='*.fastq.bam$',full.names=TRUE)
names(bamfiles) <- substr(basename(bamfiles),1,17)
fcounts.mapq1 <- featureCounts(bamfiles,annot.ext=utr3,isGTFAnnotationFile=FALSE,allowMultiOverlap=TRUE,countMultiMappingReads=FALSE,minMQS=1,nthreads=16)
save(fcounts.mapq1,file=file.path(datadir,'featureCounts_bowtie2_multioverlap_minq1_notmultimap.RData'))
     
# Read sampleinfo, prepare for compare IP vs Input, join both lanes
load(file.path(datadir,'sampleinfo.RData'))
sampleinfo$Group <- gsub('_[1-4]','',sampleinfo$Name)
colnames(sampleinfo)[8:10] <- c('Type','Factor','Replicate')
for (i in colnames(sampleinfo)) sampleinfo[,i] <- as.factor(as.character(sampleinfo[,i]))
sampleinfo$Group <- relevel(sampleinfo$Group,'Input_NI')
rownames(sampleinfo) <- sampleinfo$Name
save(sampleinfo,file=file.path(datadir,'sampleinfo_forDESeq2_may19.RData'))
load(file=file.path(datadir,'sampleinfo_forDESeq2_may19.RData'))

# Build DESeq dataset, paired by replicate, use counts mapq1
de.counts <- fcounts.mapq1$counts
colnames(de.counts) <- gsub('X.Volumes.biostats.consulting.raul_mendez.bduran_201903_ripseq.data..sequences.bam.','',gsub('.fastq.bam','',colnames(de.counts)))
de.counts <- de.counts[,sampleinfo$ID]
colnames(de.counts) <- sampleinfo$Name
dds <- DESeqDataSetFromMatrix(countData=de.counts,colData=sampleinfo,design= ~ Replicate + Group)
dds <- DESeq(dds)
summary(results(dds))
save(dds,file=file.path(datadir,'bowtie2_3utr_dds_DESeq2_fcounts_minq1_may19.RData'))

# Plot PCA
pdf(file.path(figsdir,'bowtie2_3utr_pca_DESeq2_fcounts_minq1_may19.pdf'))
rld <- rlog(dds)
plotPCA(rld,intgroup='Type')
plotPCA(rld,intgroup='Factor')
plotPCA(rld,intgroup='Replicate')
plotPCA(rld,intgroup='Group')
dev.off()

# Extract results for contrasts of interest
conts <- list(CPEB1_vs_Input=c('Group','IP_CPEB1','Input_NI'),
              CPEB2_vs_Input=c('Group','IP_CPEB2','Input_NI'),
              CPEB3_vs_Input=c('Group','IP_CPEB3','Input_NI'),
              CPEB4_vs_Input=c('Group','IP_CPEB4','Input_NI'),
              IP_NI_vs_Input=c('Group','IP_NI','Input_NI'))
reslist <- lapply(conts,function(x) results(dds,pAdjustMethod='BH',cooksCutoff=FALSE,contrast=x))

# Add rej column
reslist <- lapply(reslist,function(x) { x$rej <- sign(x$log2FoldChange) * (x$padj<0.05) * (abs(x$log2FoldChange)>log2(2)); x })
library(gtools)
do.call(smartbind,lapply(reslist,function(x) table(x$rej)))
sink(file.path(tablesdir,'bowtie2_3utr_DESeq2_fcounts_minq1_may19.txt'))
do.call(smartbind,lapply(reslist,function(x) table(x$rej)))
sink()

# Make final tables
means <- sapply(unique(unlist(conts))[-1], function(g) rowMeans(log2(counts(dds,normalized=TRUE)[,dds$Group == g]+1)))
xoutlist <- lapply(reslist,function(x) cbind(exon2gene[rownames(means),-1:-2],de.counts,means,x))
save(reslist,xoutlist,file=file.path(datadir,'bowtie2_3utr_res_DESeq2_fcounts_minq1_may19.RData'))
for (i in names(xoutlist)) write.csv(as.data.frame(xoutlist[[i]]),file=file.path(tablesdir,sprintf('bowtie2_3utr_%s_DESeq2_fcounts_minq1_may19.csv',i)))
library(openxlsx)
for (i in names(xoutlist)) write.xlsx(as.data.frame(xoutlist[[i]]),file=file.path(tablesdir,sprintf('bowtie2_3utr_%s_DESeq2_fcounts_minq1_may19.xlsx',i)))

# Save BED files
for (i in names(xoutlist))
    {
        xx <- as.data.frame(xoutlist[[i]])
        xx <- xx[!is.na(xx$rej),]
        write.table(xx[xx$rej==1,c('Exon.seqnames','Exon.start','Exon.end','log2FoldChange')],file=file.path(tdfdir,sprintf('bowtie2_3utr_%s_DESeq2_fcounts_minq1_may19.bed',i)),
                    sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE,dec='.')
    }

###########################################################################
# 2. DESEQ2, RIP vs INPUT. LONGEST 3UTR BY GENE
###########################################################################

# Map 3UTR exons to genes, keep longest
load(file=file.path(datadir,'xl92.txdb.GRanges.annot.RData'))
idx <- findOverlaps(xl92.txdb$UTR3,xl92.txdb$Genes)
exon2gene <- data.frame(Exon=as.data.frame(xl92.txdb$UTR3)[queryHits(idx),],Gene=as.data.frame(xl92.txdb$Genes)[subjectHits(idx),])
sum(is.na(exon2gene$Gene.gene_id))
exon2gene <- exon2gene[exon2gene$Exon.strand==exon2gene$Gene.strand,]
exon2gene <- exon2gene[order(exon2gene$Gene.gene_id,exon2gene$Exon.width,decreasing=TRUE),]
head(exon2gene,20)
exon2gene <- exon2gene[!duplicated(exon2gene$Gene.gene_id),]
rownames(exon2gene) <- paste0(exon2gene$Gene.gene_id,'.',exon2gene$Exon.exon_name)

# make utr3
utr3 <- exon2gene[,c('Exon.exon_name','Exon.seqnames','Exon.start','Exon.end','Exon.strand')]
colnames(utr3) <- c('GeneID','Chr','Start','End','Strand')
utr3$GeneID <- rownames(utr3)

# Annotate 3'UTR
bamfiles <- list.files(file.path(datadir,'sequences/bam'),pattern='*.fastq.filterdup.bam$',full.names=TRUE)
names(bamfiles) <- substr(basename(bamfiles),1,17)
fcounts.mapq1 <- featureCounts(bamfiles,annot.ext=utr3,isGTFAnnotationFile=FALSE,allowMultiOverlap=TRUE,countMultiMappingReads=FALSE,minMQS=1,nthreads=16)
save(fcounts.mapq1,file=file.path(datadir,'featureCounts_bowtie2_longest3utr_filterdup_singleoverlap_mapq1_notmultimap.RData'))
     
# Read sampleinfo
load(file=file.path(datadir,'sampleinfo_forDESeq2_may19.RData'))

# Build DESeq dataset, paired by replicate, use counts mapq1
de.counts <- fcounts.mapq1$counts
colnames(de.counts) <- gsub('X.Volumes.biostats.consulting.raul_mendez.bduran_201903_ripseq.data..sequences.bam.','',gsub('.fastq.filterdup.bam','',colnames(de.counts)))
de.counts <- de.counts[,sampleinfo$ID]
colnames(de.counts) <- sampleinfo$Name
dds <- DESeqDataSetFromMatrix(countData=de.counts,colData=sampleinfo,design= ~ Replicate + Group)
dds <- DESeq(dds)
summary(results(dds))
save(dds,file=file.path(datadir,'bowtie2_longest_3utr_dds_DESeq2_fcounts_mapq1_may19.RData'))

# Plot PCA
pdf(file.path(figsdir,'bowtie2_longest_3utr_pca_DESeq2_fcounts_mapq1_may19.pdf'))
rld <- rlog(dds)
plotPCA(rld,intgroup='Type')
plotPCA(rld,intgroup='Factor')
plotPCA(rld,intgroup='Replicate')
plotPCA(rld,intgroup='Group')
dev.off()

# Extract results for contrasts of interest
conts <- list(CPEB1_vs_Input=c('Group','IP_CPEB1','Input_NI'),
              CPEB2_vs_Input=c('Group','IP_CPEB2','Input_NI'),
              CPEB3_vs_Input=c('Group','IP_CPEB3','Input_NI'),
              CPEB4_vs_Input=c('Group','IP_CPEB4','Input_NI'),
              IP_NI_vs_Input=c('Group','IP_NI','Input_NI'))
reslist <- lapply(conts,function(x) results(dds,pAdjustMethod='BH',cooksCutoff=FALSE,contrast=x))

# Add rej column
reslist <- lapply(reslist,function(x) { x$rej <- sign(x$log2FoldChange) * (x$padj<0.05) * (abs(x$log2FoldChange)>log2(2)); x })
library(gtools)
do.call(smartbind,lapply(reslist,function(x) table(x$rej)))
sink(file.path(tablesdir,'DESeq2_fcounts_longest3utr_mapq1_may19.txt'))
do.call(smartbind,lapply(reslist,function(x) table(x$rej)))
sink()

# Make final tables
means <- sapply(unique(unlist(conts))[-1], function(g) rowMeans(log2(counts(dds,normalized=TRUE)[,dds$Group == g]+1)))
xoutlist <- lapply(reslist,function(x) cbind(exon2gene[rownames(means),-1:-2],de.counts,means,x))
save(reslist,xoutlist,file=file.path(datadir,'bowtie2_longest3utr_res_DESeq2_fcounts_mapq1_may19.RData'))
for (i in names(xoutlist)) write.csv(as.data.frame(xoutlist[[i]]),file=file.path(tablesdir,sprintf('bowtie2_longest3utr_%s_DESeq2_fcounts_mapq1_may19.csv',i)))
library(openxlsx)
for (i in names(xoutlist)) write.xlsx(as.data.frame(xoutlist[[i]]),file=file.path(tablesdir,sprintf('bowtie2_longest3utr_%s_DESeq2_fcounts_mapq1_may19.xlsx',i)))

for (i in names(xoutlist))
    {
        xx <- as.data.frame(xoutlist[[i]])
        xx <- xx[!is.na(xx$rej),]
        write.table(xx[xx$rej==1,c('Exon.seqnames','Exon.start','Exon.end','log2FoldChange')],file=file.path(tdfdir,sprintf('bowtie2_longest3utr_%s_DESeq2_fcounts_mapq1_may19.bed',i)),
                    sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE,dec='.')
    }

# Save main effect tables for GSEA
rnks <- xoutlist
rnks <- lapply(rnks,function(x) aggregate(x$log2FoldChange,by=list(symbol=x$symbol),FUN=function(y) mean(y,na.rm=TRUE)))
rnks <- lapply(rnks,function(x) x[!is.na(x$x),])
for (i in names(rnks)) write.table(rnks[[i]],file=file.path(tablesdir,sprintf('bowtie2_3utr_%s_DESeq2_fcounts_mapq1_190618.rnk',i)),sep='\t',quote=FALSE,dec='.',row.names=FALSE,col.names=FALSE)

# Check some candidate genes
mygenes <- c('ccnd1','rspo1','tnfsf11','cpeb2','cpeb4')
lapply(mygenes,function(g) lapply(xoutlist,function(x) x[grep(g,tolower(x$symbol)),]))

#####################################################################################################
# 3. DESEQ2, RIP.WT VS INPUT.WT / RIP.KO VS INPUT.KO (INTERACTION ANALYSIS). LONGEST 3UTR BY GENE
#####################################################################################################

# Refactor Genotype so that actually in this case KO is our control for WT
sampleinfo$Genotype <- relevel(sampleinfo$Genotype,'KO')

dds <- DESeqDataSetFromMatrix(countData=de.counts,colData=sampleinfo,design= ~ Replicate + Type + Genotype + Type:Genotype)
dds <- DESeq(dds)
resultsNames(dds)
save(dds,file=file.path(datadir,'bowtie2_3utr_dds_Interaction_DESeq2_fcounts_mapq1_190618.RData'))

# Make final table
xout.int <- results(dds,pAdjustMethod='BH',cooksCutoff=FALSE,name='TypeIP.GenotypeWT')
unique(rownames(annot)==rownames(xout.int))
means <- sapply(c('WT_IP','WT_Input','KO_IP','KO_Input'), function(g) rowMeans(log2(counts(dds,normalized=TRUE)[,dds$Group == g]+1)))
xout.int <- cbind(as.data.frame(annot),de.counts,means,xout.int)
xout.int$rej <- sign(xout.int$log2FoldChange) * (xout.int$pvalue<0.10) * (abs(xout.int$log2FoldChange)>log2(1.5))
table(xout.int$rej)

write.csv(as.data.frame(xout.int),file=file.path(tablesdir,'bowtie2_3utr_Interaction_DESeq2_fcounts_mapq1_190618.csv'))
library(openxlsx)
write.xlsx(as.data.frame(xout.int),file=file.path(tablesdir,'bowtie2_3utr_Interaction_DESeq2_fcounts_mapq1_190618.xlsx'))

###########################################################################
# 4. GENE SET ENRICHMENT ANALYSIS MAIN EFFECTS
###########################################################################

# List rnks
rnks <- list.files(tablesdir,'.rnk$',full.names=TRUE)
file.exists(rnks)
names(rnks) <- gsub('bowtie2_3utr_','',gsub('_DESeq2_fcounts_mapq1_190618.rnk','',basename(rnks)))

# Genesets, all V3
gsets <- paste(file.path('/Volumes/biostats/databases/BroadGSEA/genesets/mmusculus',c('GOBP','GOCC','GOMF','Broad_Hallmarks','GOSLIM','KEGG')),'v3.gmt',sep='_')
names(gsets) <- gsub('.gmt','',basename(gsets))
file.exists(gsets)

# Comps to make
comps <- data.frame(rep(rnks,each=length(gsets)),gsets)
comps <- apply(comps,2,as.character)
rownames(comps) <- apply(data.frame(rep(names(rnks),each=length(gsets)),names(gsets)),1,paste,collapse='.')

# Do it
source("/Volumes/biostats/oreina/scripts/R/GSEA_mod3.R")
system(sprintf('mkdir %s',file.path(tablesdir,'GSEA')))
ans <- mclapply(rownames(comps),function(i) runGSEApreRanked.mod(comps[i,1],customgeneset=comps[i,2],score="weighted",minSize=10,label=i,numplots=50,outdir=file.path(tablesdir,'GSEA'),local=TRUE),mc.cores=12)

# Rename output dirs
dir.in <- list.files(file.path(tablesdir,'GSEA'),pattern='v3',full.names=TRUE)
dir.out <- gsub('GseaPreranked.[0-9]+','GseaPreranked',dir.in)
for (i in 1:length(dir.in)) system(sprintf('mv %s %s',dir.in[i],dir.out[i]))

# Make GSEA heatmaps
library(gplots)
hclustfun <- function(x) hclust(dist(x),method='ward')
setwd(file.path(tablesdir,'GSEA'))
patterns <- c('Broad_Hallmarks_v3','GOSLIM_v3','KEGG_v3')

# For every geneset collection considered
for (pattern in patterns)
    {
        gseadirs <- system(sprintf('ls | grep %s',pattern),intern=TRUE)
        names(gseadirs) <- gsub('.GseaPreranked','',gseadirs)
        nes <- nes.bak <- getGseaNES(gseadirs)
        pdf(file.path(figsdir,sprintf('NES_GSEA_%s.pdf',pattern)),height=25,width=15)
        fdr <- fdrFormat(nes[,grep('FDR',colnames(nes))])
        colnames(nes) <- gsub(pattern,'',colnames(nes))
        heatmap.2(as.matrix(nes[,grep('NES',colnames(nes))]),hclustfun=hclustfun,Colv=TRUE,Rowv=TRUE,col=bluered(100),trace='none',margins=c(25,35),cellnote=fdr,notecol='black',notecex=1.2,dendrogram='none',keysize=.7,cexCol=1.5)
        title(sprintf('GSEA NES Heatmap %s\nFDR + / * / ** / *** for 0.25 / 0.10 / 0.05 / 0.01',pattern))
        dev.off()
    }

# Make summary table
tt <- matrix(nrow=length(gsets),ncol=length(rnks),'')
rownames(tt) <- names(gsets)
colnames(tt) <- names(rnks)
for (i in rownames(tt)) for (j in colnames(tt)) tt[i,j] <- sprintf('<a href="GSEA/%s.%s.GseaPreranked/index.html">[html]</a>',j,i)
tt <- as.data.frame(tt)
library(R2HTML)
HTML(tt,file=file.path(tablesdir,'GSEA_290618.html'),append=FALSE)

### THE END

