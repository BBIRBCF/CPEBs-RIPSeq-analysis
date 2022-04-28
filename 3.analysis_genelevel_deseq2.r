##########################################################################
## PROJECT: BDURAN RIPSEQ
## DESCRIPTION OF FILE: GENES DIFFERENTIAL ANALYSIS VIA DESEQ2
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

# Use gene information
load(file=file.path(datadir,'xl92.txdb.GRanges.annot.RData'))
genes <- as.data.frame(xl92.txdb$Genes)
genes <- genes[,c('gene_id','seqnames','start','end','strand')]
colnames(genes) <- c('GeneID','Chr','Start','End','Strand')

# Count reads over GENES
bamfiles <- list.files(file.path(datadir,'sequences/bam'),pattern='*.fastq.bam$',full.names=TRUE)
names(bamfiles) <- substr(basename(bamfiles),1,17)
fcounts.mapq1 <- featureCounts(bamfiles,annot.ext=genes,isGTFAnnotationFile=FALSE,allowMultiOverlap=TRUE,countMultiMappingReads=FALSE,minMQS=1,nthreads=16)
save(fcounts.mapq1,file=file.path(datadir,'featureCounts_bowtie2_genes_multioverlap_minq1_notmultimap.RData'))
     
# Read sampleinfo, prepare for compare IP vs Input, join both lanes
load(file=file.path(datadir,'featureCounts_bowtie2_genes_multioverlap_minq1_notmultimap.RData'))
load(file=file.path(datadir,'sampleinfo_forDESeq2_may19.RData'))

# Build DESeq dataset, paired by replicate, use counts mapq1
de.counts <- fcounts.mapq1$counts
colnames(de.counts) <- gsub('X.Volumes.biostats.consulting.raul_mendez.bduran_201903_ripseq.data..sequences.bam.','',gsub('.fastq.bam','',colnames(de.counts)))
de.counts <- de.counts[,sampleinfo$ID]
colnames(de.counts) <- sampleinfo$Name
dds <- DESeqDataSetFromMatrix(countData=de.counts,colData=sampleinfo,design= ~ Replicate + Group)
dds <- DESeq(dds)
summary(results(dds))
save(dds,file=file.path(datadir,'bowtie2_genes_dds_DESeq2_fcounts_minq1_may19.RData'))

# Plot PCA
pdf(file.path(figsdir,'bowtie2_genes_pca_DESeq2_fcounts_minq1_may19.pdf'))
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
sink(file.path(tablesdir,'bowtie2_genes_DESeq2_fcounts_minq1_may19.txt'))
do.call(smartbind,lapply(reslist,function(x) table(x$rej)))
sink()

# Make final tables
means <- sapply(unique(unlist(conts))[-1], function(g) rowMeans(log2(counts(dds,normalized=TRUE)[,dds$Group == g]+1)))
xoutlist <- lapply(reslist,function(x) cbind(genes[rownames(means),],de.counts,means,x))
save(reslist,xoutlist,file=file.path(datadir,'bowtie2_genes_res_DESeq2_fcounts_minq1_may19.RData'))
for (i in names(xoutlist)) write.csv(as.data.frame(xoutlist[[i]]),file=file.path(tablesdir,sprintf('bowtie2_genes_%s_DESeq2_fcounts_minq1_may19.csv',i)))
library(openxlsx)
for (i in names(xoutlist)) write.xlsx(as.data.frame(xoutlist[[i]]),file=file.path(tablesdir,sprintf('bowtie2_genes_%s_DESeq2_fcounts_minq1_may19.xlsx',i)))

# Save BED files
for (i in names(xoutlist))
    {
        xx <- as.data.frame(xoutlist[[i]])
        xx <- xx[!is.na(xx$rej),]
        write.table(xx[xx$rej==1,c('Chr','Start','End','log2FoldChange')],file=file.path(tdfdir,sprintf('bowtie2_genes_%s_DESeq2_fcounts_minq1_may19.bed',i)),
                    sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE,dec='.')
    }

### THE END
