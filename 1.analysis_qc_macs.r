###########################################################################
# BDURAN RIPSEQ 201903
###########################################################################

library(parallel)
library(Rsamtools)
library(htSeqTools)
library(org.Mm.eg.db)
library(R2HTML)
library(htSeqTools)
library(ChIPpeakAnno)
library(biomaRt)
library(openxlsx)

routinesdir <- "/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/routines"
datadir <- "/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data"
fastqdir1 <- "/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/sequences/fastq"
fastqdir2 <- "/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/sequences/bam_rRNA"
bamdatadir <- "/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/sequences/bam_rRNA_clean"
tablesdir <- "/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables"
tdfdir1 <- "/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables/files_4_igv"
figsdir <- "/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/figs"

###########################################################################
# 1. READ SAMPLE INFO FILE, ALIGNED READS AND STORE
###########################################################################

# Read sampleinfo
sampleinfo <- read.xlsx(file.path(datadir,'sequences/BduranMar19_RIPSeq Barcodes ID.xlsx'))
sampleinfo$Sample.Name <- gsub(' ','_',sampleinfo$Sample.Name)
sampleinfo$Sample.Name <- gsub('Input_','Input',sampleinfo$Sample.Name)
sampleinfo$Type <- 'IP'
sampleinfo$Type[grep('INPUT',toupper(sampleinfo$Sample.Name))] <- 'Input'
sampleinfo$Type[grep('IGG',toupper(sampleinfo$Sample.Name))] <- 'IgG'
sampleinfo$Genotype <- unlist(lapply(strsplit(sampleinfo[,2],'_'),function(x) x[2]))
sampleinfo$Prefix <- gsub('X','',make.names(unlist(lapply(strsplit(sampleinfo[,2],'_'),function(x) x[1]))))
sampleinfo$ID <- apply(sampleinfo[,c(3,5,7)],1,paste,collapse='_')
sampleinfo$Name <- apply(sampleinfo[,c(10,9,8)],1,paste,collapse='_')
rownames(sampleinfo) <- sampleinfo$ID

## read and store raw aligned reads
## BAM parameters
bams <- list.files(file.path(bamdatadir),'*fastq.un.m1.fq.bam$',full.names=TRUE)
names(bams) <- substr(basename(bams),1,17)
what=c("flag","rname", "strand", "pos", "qwidth","mapq")
param=ScanBamParam(what = what)
bams

# scan Files
bams <- mclapply(bams,function(x) scanBam(x,param=param)[[1]],mc.cores=8)
x <- lapply(bams,function(x) { sel <- !is.na(x$rname);as(GRanges(seqnames=as.character(Rle(x$rname[sel])),ranges=IRanges(x$pos[sel],width=x$qwidth[sel]),strand=x$strand[sel],mapq=x$mapq[sel]),'RangedData') })
save(x,file=file.path(datadir,'seqs_bowtie2.RData'))

###########################################################################
# 2. GENERATE MDS PLOT
###########################################################################

# prepare plotinfo
plotinfo <- sampleinfo
plotinfo <- plotinfo[names(x),]

# Plot MDS
mymds <- cmds(as.list(x),mc.cores=length(x)/2)

pdf(file.path(figsdir,'cmds.pdf'))
# Full
xy <- mymds@points
plot(xy[,1],xy[,2],pch=as.numeric(as.factor(plotinfo$Type)),col=factor(plotinfo$Genotype),xlim=2*range(xy),ylim=2*range(xy),cex=2,xlab='',ylab='',main='MDS (Bowtie2)')
text(plotinfo$Name,x=xy[,1],y=xy[,2],pos=4,cex=.75)
legend('bottomleft',pch=c(1:3),legend=levels(factor(plotinfo$Type)))
legend('topleft',col=unique(sort(as.numeric(factor(unique(plotinfo$Genotype))))),pch=21,legend=levels(factor(plotinfo$Genotype)))
legend('bottomright',legend=sprintf('R2=%.3f',mymds@R.square))
dev.off()

save(sampleinfo,plotinfo,file=file.path(datadir,'sampleinfo.RData'))

###########################################################################
# 3. RUN MACS AND PROCESS PEAKS
###########################################################################

# Read sampleinfo
load(file=file.path(datadir,'sampleinfo.RData'))
sampleinfo$Filename <- paste0(rownames(sampleinfo),'.fastq.un.m1.fq.bam')
file.exists(file.path(bamdatadir,sampleinfo$Filename))
colnames(sampleinfo)[8:10] <- c('Type','Factor','Replicate')

# Prepare contrasts
inputs <- sampleinfo[sampleinfo$Type=='Input',]
conts <- data.frame(Sample=sampleinfo[sampleinfo$Type=='IP',],Control=inputs)
rownames(conts) <- paste0(conts$Sample.Name,'_vs_Input')
conts$Sample.Factor==conts$Control.Factor
conts$Sample.Replicate==conts$Control.Replicate
write.csv(conts,file=file.path(datadir,'conts_bowtie2.csv'))

# Run MACS
# Call regions with MACS 1.4.2
setwd(bamdatadir)
file.exists(conts$Sample.Filename)
file.exists(conts$Control.Filename)
outdir <- file.path(tablesdir,'macs_bowtie2')
dir.create(outdir)

macs.call <- sprintf('/software/MACS/MACS-1.4.2/bin/macs14 -t %s -c %s -f BAM --keep-dup=all -g mm -n %s/%s',conts$Sample.Filename,conts$Control.Filename,outdir,gsub(' ','_',rownames(conts)))
macs.out <- mclapply(macs.call,function(command) system(command,intern=TRUE,wait=FALSE),mc.cores=nrow(conts))

# Load MACS peak tables and store as rdata
peaks <- list.files(file.path(tablesdir,'macs_bowtie2'),'peaks.xls$',full.names=TRUE)
peaks <- peaks[-c(grep('negative',peaks))]
names(peaks) <- gsub('_peaks.xls','',basename(peaks))
peaks <- lapply(peaks,read.delim,header=TRUE,as.is=TRUE,comment.char='#',skip=1)
do.call(rbind,lapply(peaks,function(x) summary(x[,9])))
peaks <- lapply(peaks,function(x) { colnames(x) <- c('space','start','end','width','summit','nreads','log10.pval','fold','fdr'); RangedData(x) })
lapply(peaks,nrow)
save(peaks,file=file.path(datadir,'peaks_bowtie2_macs_rangeddata.RData'))

###########################################################################
# 4. ANNOTATE PEAKS WITH CHIPSEEKER
###########################################################################

load(file=file.path(datadir,'peaks_bowtie2_macs_rangeddata.RData'))

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Annotate
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peaksanno <- lapply(peaks,function(x) annotatePeak(as(x,'GRanges'),TxDb=txdb,tssRegion=c(-1000,1000)))
save(peaksanno,file=file.path(datadir,'peaks_bowtie2_macs_chipseeker.RData'))
load(file=file.path(datadir,'peaks_bowtie2_macs_chipseeker.RData'))

# Save Excel tables for evaluation
for (i in names(peaksanno)) write.xlsx(as.data.frame(peaksanno[[i]]),file=file.path(tablesdir,sprintf('%s_peaks_macs_chipseeker.xlsx',i)))

# Plot peak distribution per sample
# upsetplot by Sample
for (i in names(peaksanno))
    {
        pdf(file.path(figsdir,sprintf('%s_PeakDistribution_mm10_ChIPseeker.pdf',i)),width=10,height=5)
        upsetplot(peaksanno[[i]],vennpie=TRUE); text(x=-3.5,y=-.25,i,font=2)
        dev.off()
    }

###########################################################################
# 5. GENERATE LOG2RPKM GENE AND UTR LEVEL TABLES
###########################################################################

# Load annotation
load(file='/Volumes/biostats/databases/RefGenome/xlaevis/xlaevis_92_ucsc/XENLA_9.2_Xenbase.RData')

# Load sequences
load(file=file.path(datadir,'seqs_bowtie2.RData'))

# getRPKM function
getRPKM <- function(peaks,x)
    {
        mnames <- intersect(names(peaks),names(x))
        peaks <- peaks[mnames]
        rpm <- countOverlaps(query=peaks,subject=x[names(peaks)])
        counts <- unlist(rpm)
        nreads <- nrow(x)
        rpm <- 10^9*unlist(rpm)/(nreads*(end(peaks)-start(peaks)+1))
        return(log2(rpm+1))
    }
     
# getCounts function
getCounts <- function(peaks,x)
    {
        mnames <- intersect(names(peaks),names(x))
        peaks <- peaks[mnames]
        rpm <- countOverlaps(query=peaks,subject=x[names(peaks)])
        counts <- unlist(rpm)
        return(counts)
    }

# Compute RPKM per sample
rpkm <- as.data.frame(do.call(cbind,mclapply(x,function(y) getRPKM(as(xl92.genes,'RangedData'),y),mc.cores=length(x))))
counts <- as.data.frame(do.call(cbind,mclapply(x,function(y) getCounts(as(xl92.genes,'RangedData'),y),mc.cores=length(x))))
unique(rownames(counts)==rownames(rpkm))
rownames(rpkm) <- rownames(counts) <- unlist(lapply(strsplit(rownames(rpkm),'.',fixed=TRUE),function(x) paste(unlist(x[2:length(x)]),collapse='.')))

# Join with gene info
genes.info <- as.data.frame(xl92.genes)
table(rownames(rpkm) %in% rownames(genes.info))
counts$ID <- rpkm$ID <- rownames(counts)
counts <- data.frame(genes.info[rownames(counts),],counts)
rpkm <- data.frame(genes.info[rownames(rpkm),],rpkm)
# ensure no partial matching is done
unique(counts$gene_id==counts$ID)
unique(rpkm$gene_id==rpkm$ID)

# Save csv and XLSX
library(openxlsx)
write.csv(counts,file=file.path(tablesdir,'samplecounts_xl92_genes.csv'),row.names=TRUE)
write.xlsx(counts,file=file.path(tablesdir,'samplecounts_xl92_genes.xlsx'),row.names=TRUE)
write.csv(rpkm,file=file.path(tablesdir,'samplerpkm_xl92_genes.csv'),row.names=TRUE)
write.xlsx(rpkm,file=file.path(tablesdir,'samplerpkm_xl92_genes.xlsx'),row.names=TRUE)

###########################################################################
# 6. RUN DIFFBIND
###########################################################################

library(DiffBind)

# Read contrast info
conts <- read.csv(file=file.path(datadir,'conts_m1_filterdup.csv'),header=TRUE,as.is=TRUE)
conts$Factor <- unlist(lapply(strsplit(conts$Sample.Sample.Name,' '),function(x) x[length(x)]))
conts$Treatment <- unlist(lapply(strsplit(conts$Sample.Sample.Name,' '),function(x) x[2]))
conts$Condition <- paste(conts$Factor,conts$Treatment,sep='.')
beds <- list.files(file.path(tablesdir,'macs_m1'),pattern='_peaks.bed',full.names=TRUE)
names(beds) <- gsub('_peaks.bed','',basename(beds))
conts$Peaks <- beds[gsub(' ','_',paste(conts$Control.Condition,'IP',conts$Factor))]
conts$PeakCaller='macs'
conts$Replicate <- c(rep('R1',4),rep('R2',4))
conts$Tissue <- 'S2'
conts <- conts[,c('Sample.Sample.Name','Tissue','Factor','Condition','Treatment','Replicate','Sample.Filename','Control.Sample.Name','Control.Filename','Peaks','PeakCaller')]
conts$Sample.Filename <- file.path(bamdatadir,conts$Sample.Filename)
conts$Control.Filename <- file.path(bamdatadir,conts$Control.Filename)
colnames(conts) <- c('SampleID','Tissue','Factor','Condition','Treatment','Replicate','bamReads','ControlID','bamControl','Peaks','PeakCaller')
for (i in 1) conts[,i] <- make.names(conts[,i])
write.csv(conts,file=file.path(datadir,'conts_m1_filterdup_diffbind.csv'))

# Perform diffbind
setwd(datadir)
samples <- dba(sampleSheet='conts_m1_filterdup_diffbind.csv',peakCaller='macs',peakFormat='macs')
samples <- dba.count(samples,minOverlap=1,score=DBA_SCORE_RPKM_FOLD,bLog=TRUE)
samples <- dba.contrast(samples,categories=DBA_CONDITION,minMembers=2)
samples
samples.bak <- samples

samples <- dba.analyze(samples,method=c(DBA_EDGER,DBA_DESEQ2),bSubControl=FALSE,bTagwise=FALSE)

sink(file.path(tablesdir,'diffbind_results.txt'))
dba.show(samples,bContrast=TRUE)
sink()

conts <- dba.show(samples,bContrast=TRUE)
conts

reslist.edger <- lapply(1:nrow(conts),function(i) dba.report(samples,method=DBA_EDGER,contrast=i))
reslist.deseq2 <- lapply(1:nrow(conts),function(i) dba.report(samples,method=DBA_DESEQ2,contrast=i))
names(reslist.edger) <- names(reslist.deseq2) <- paste(conts$Group1,conts$Group2,sep='_')
for (i in names(reslist.edger)) write.table(as.data.frame(reslist.edger[[i]])[,c('seqnames','start','end','Fold')],file=file.path(tablesdir,sprintf('diffbind_%s_edger.bedgraph',i)),sep='\t',quote=FALSE,dec='.',row.names=F,col.names=F)

### THE END
