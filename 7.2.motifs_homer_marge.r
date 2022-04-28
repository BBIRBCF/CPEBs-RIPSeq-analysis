##########################################################################
## PROJECT: BDURAN RIPSEQ
## DESCRIPTION OF FILE: MOTIF SEARCHING VIA HOMER+MARGE
## DATE: JAN 2020
###########################################################################

library(marge)
library(openxlsx)
library(BSgenome)

.datadir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/'
.datadir2 <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/'
.datadir1 <- '/Volumes/biostats/consulting/raul_mendez/bduran_201810_ripseq/data/'
.reportsdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports'
.tablesdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables'
.tdfdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables/files_4_igv'
.figsdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/figs'

###########################################################################
## chunk 0.1 QUICK MARGE TEST
###########################################################################

## Extract UTRs
## library(GenomicFeatures)
## gff <- makeTxDbFromGFF(file.path(datadir,'genomes/XENLA_9.2_Xenbase.gff3'))
## genes <- transcriptsBy(gff)
## genes <- mclapply(1:length(genes),function(i) genes[[i]]$GeneID <<- i,mc.cores=10)
## geneids <- rep(names(genes),unlist(lapply(genes,length)))
## genes <- do.call(rbind,mclapply(genes,as.data.frame,mc.cores=10))
## genes$GeneID <- geneids
## utr3 <- threeUTRsByTranscript(gff)
## utr3 <- do.call(rbind,mclapply(utr3,as.data.frame,mc.cores=10))
## utr3$tx_id <- rownames(utr3)
## genes.utr <- merge(genes,utr3,by.x='tx_id',by.y='tx_id',all.x=TRUE,all.y=TRUE)
## save(genes,utr3,genes.utr,file=file.path(datadir,'genes_utr3_xenlae92.RData'))
load(file=file.path(.datadir,'genes_utr3_xenlae92.RData'))

## Load output
xout <- read.xlsx(file.path(.tablesdir,'bowtie2_genes_allbatches_CPEB1_vs_CPEB234_DESeq2_Diff_fcounts_minq1_nov19_v2.xlsx'))
table(xout$GeneID %in% genes.utr$GeneID)

## Dummy set, 1kb downstream of CPEB1 specific genes
cpeb1 <- xout[xout$rej.CPEB1_vs_Input==1 & xout$rej.CPEB1_vs_NI==1 & xout$rej.CPEB1_vs_CPEB234==1,]
cpeb1 <- cpeb1[!is.na(cpeb1$rej.CPEB1_vs_Input) & !is.na(cpeb1$rej.CPEB1_vs_NI) & !is.na(cpeb1$rej.CPEB1_vs_CPEB234),]

## Extract cpeb1 genes 3utrs, remove genes without known utr, keep longest utr per gene
genes.utr <- genes.utr[!is.na(genes.utr$width.y),]
genes.utr <- genes.utr[order(genes.utr$GeneID,genes.utr$width.y,decreasing=TRUE),]
genes.utr <- genes.utr[!duplicated(genes.utr$GeneID),]
cpeb1.utr <- genes.utr[genes.utr$GeneID %in% cpeb1$GeneID,]
xx <- cpeb1.utr[,c(2,3,4,8)]

## Check homer
check_homer()

## Test marge, full genes
dir.create(.reportsdir,'motifs_homer_cpeb1_140120')
out <- find_motifs_genome(xx[1,],path=file.path(.reportsdir,'motifs_homer_cpeb1_140120'),genome='xenLae92',scan_size=100,cores=10,only_denovo=TRUE,overwrite=TRUE)

## It works !

###########################################################################
## chunk 1. RUN HOMER FOR THE DIFFERENT TARGET LISTS
###########################################################################

rm(list=ls())

## Load output
xout <- read.xlsx(file.path(.tablesdir,'bowtie2_genes_allbatches_CPEB1_vs_CPEB234_DESeq2_Diff_fcounts_minq1_nov19_v2.xlsx'))

# Define individual targets for all 4 CPEBs according to Berta's criteria
targets <- as.data.frame(do.call(cbind,mclapply(1:4,function(i)
                            {
                                log2fc.1 <- xout[,sprintf('log2FC.CPEB%s_vs_Input',i)]
                                log2fc.2 <- xout[,sprintf('log2FC.CPEB%s_vs_NI',i)]
                                padj.1 <- xout[,sprintf('padj.CPEB%s_vs_Input',i)]
                                padj.2 <- xout[,sprintf('padj.CPEB%s_vs_NI',i)]
                                ans <- log2fc.1 >= 2 & padj.1 <= 0.05 & log2fc.2 >= 1 & padj.2 <= 0.05
                                ans[is.na(ans)] <- FALSE
                                as.numeric(ans)
                            },mc.cores=4)))
colnames(targets) <- paste0('CPEB',1:4)
rownames(targets) <- xout$GeneID
head(targets)
colSums(targets)
length(rownames(targets)[rowSums(targets[,1:4])>0])

# Add CPEB234
targets$CPEB1234 <- as.numeric(rowSums(targets[,c(1,2,3,4)])>0)
targets$CPEB234 <- as.numeric(rowSums(targets[,c(2,3,4)])>0)
colSums(targets)
targets$rej.CPEB1_vs_CPEB234 <- sign(xout$log2FoldChange.CPEB1_vs_CPEB234) * as.numeric(abs(xout$log2FoldChange.CPEB1_vs_CPEB234) >= 1) * as.numeric(xout$padj.CPEB1_vs_CPEB234 <= 0.05)
targets$rej.CPEB1_vs_CPEB234[is.na(targets$rej.CPEB1_vs_CPEB234)] <- 0
colSums(abs(targets)) # Total Diff enriched
targets$rej.CPEB1_vs_CPEB234 <- sign(targets$rej.CPEB1_vs_CPEB234) * as.numeric(targets$rej.CPEB1_vs_CPEB234 & rowSums(targets[,1:4])>0)
colSums(abs(targets)) # 649 CPEB targets (of at least one CPEB) are differentially enriched in 1 vs 234
table(targets$rej.CPEB1_vs_CPEB234)

# Build binary matrix of targets for sanity checks, add Gene Info and longest 3'UTR info
load(file=file.path(.datadir,'genes_utr3_xenlae92.RData'))
head(genes.utr)
table(rownames(targets) %in% genes.utr$GeneID)
sum(is.na(genes.utr$width.y))
targets$hasUTR3 <- as.numeric(rownames(targets) %in% genes.utr$GeneID[!is.na(genes.utr$width.y)])
colSums(abs(targets))
sel.homer <- targets$rej.CPEB1_vs_CPEB234!=0 & targets$hasUTR3==1
table(sel.homer)

# Filter for longest UTR and sanity checks
genes.utr <- genes.utr[!is.na(genes.utr$GeneID) & !is.na(genes.utr$width.y),]
genes.utr <- genes.utr[order(genes.utr$GeneID,genes.utr$width.y,decreasing=TRUE),]
genes.utr <- genes.utr[!duplicated(genes.utr$GeneID),]
rownames(genes.utr) <- genes.utr$GeneID
table(rownames(targets) %in% genes.utr$GeneID)
table(rownames(targets) %in% genes$GeneID)

# Save
save(genes.utr,file=file.path(.datadir,'genes_3utr_xl92_jan20.RData'))

# Export genes.utr fasta, use Homer
## homerTools extract <peak/BED file> <FASTA directory or file location> [-fa]
xl92 <- '/Volumes/biostats/databases/RefGenome/xlaevis/xlaevis_92_ucsc/xenLae2_061118.fa'
write.table(utr3.all,file=file.path(.tablesdir,'xlaevis_ucsc_92_3utr_longest.bed'),row.names=F,col.names=F,quote=F,sep='\t',dec='.')
utr3.bed <- file.path(.tablesdir,'xlaevis_ucsc_92_3utr_longest.bed')
utr3.fa <- file.path(.tablesdir,'xlaevis_ucsc_92_3utr_longest.fa')
extract <- '/software/homer/bin/homerTools'
cmd <- sprintf('%s extract %s %s -fa > %s',extract,utr3.bed,xl92,utr3.fa)
system(cmd)

# Prepare targets table for interim use
targets.xout <- targets
targets.xout$GeneID <- rownames(targets)
targets.xout <- merge(targets.xout,genes.utr,by.x='GeneID',by.y='GeneID',all.x=TRUE,all.y=FALSE)
colnames(targets.xout) <- gsub('.x','.tx_of_3utr',colnames(targets.xout),fixed=TRUE)
colnames(targets.xout) <- gsub('.y','.3utr',colnames(targets.xout),fixed=TRUE)
rownames(xout) <- xout$GeneID

targets.xout <- cbind(xout[targets.xout$GeneID,1:5],targets.xout)
unique(targets.xout[,1]==targets.xout[,6])
write.xlsx(targets.xout[,-6],file=file.path(.tablesdir,'xl92_genes_longest3UTR_cpebtargets_jan20.xlsx'))

# Split for homer UP (CPEB1) or DOWN (CPEB234)
targets.bak <- targets
targets <- targets[targets$hasUTR3==1,]
sel.homer.cpeb1 <- targets$rej.CPEB1_vs_CPEB234==1
sel.homer.cpeb234 <- targets$rej.CPEB1_vs_CPEB234==-1
sel.homer.cpeb1234 <- targets$CPEB1234==1

# Generate dfs
utr3.all <- genes.utr[,c('seqnames.y','start.y','end.y','GeneID')]
utr3.cpeb1 <- genes.utr[rownames(targets)[sel.homer.cpeb1],c('seqnames.y','start.y','end.y','GeneID')]
utr3.cpeb234 <- genes.utr[rownames(targets)[sel.homer.cpeb234],c('seqnames.y','start.y','end.y','GeneID')]
utr3.cpeb1234 <- genes.utr[rownames(targets)[sel.homer.cpeb1234],c('seqnames.y','start.y','end.y','GeneID')]

## Check homer
check_homer()

## Launch Homer, CPEB1234 vs 3UTR BG (targets of at least one of them), updated 210620 to use BG of all longest 3UTR
dir.create(.reportsdir,'motifs_homer_cpeb1234_vs_bg_210620')
out <- find_motifs_genome(utr3.cpeb1234,background=utr3.all,path=file.path(.reportsdir,'motifs_homer_cpeb1234_vs_bg_210620'),genome='xenLae92',scan_size='given',cores=10,only_known=FALSE,only_denovo=FALSE,overwrite=TRUE)

## Launch Homer, CPEB1 vs BG
dir.create(.reportsdir,'motifs_homer_cpeb1_vs_bg_210620')
out <- find_motifs_genome(utr3.cpeb1,background=utr3.all,path=file.path(.reportsdir,'motifs_homer_cpeb1_vs_bg_210620'),genome='xenLae92',scan_size='given',cores=10,only_known=FALSE,only_denovo=FALSE,overwrite=TRUE)

## Launch Homer, CPEB234 vs BG
dir.create(.reportsdir,'motifs_homer_cpeb234_vs_bg_210620')
out <- find_motifs_genome(utr3.cpeb234,background=utr3.all,path=file.path(.reportsdir,'motifs_homer_cpeb234_vs_bg_210620'),genome='xenLae92',scan_size='given',cores=10,only_known=FALSE,only_denovo=FALSE,overwrite=TRUE)

## Launch Homer, CPEB1 vs CPEB234
dir.create(.reportsdir,'motifs_homer_cpeb1_vs_cpeb234_170120')
out <- find_motifs_genome(utr3.cpeb1,background=utr3.cpeb234,path=file.path(.reportsdir,'motifs_homer_cpeb1_vs_cpeb234_170120'),genome='xenLae92',scan_size='given',cores=10,only_known=FALSE,only_denovo=FALSE,overwrite=TRUE)

## Launch Homer, CPEB234 vs CPEB1
dir.create(.reportsdir,'motifs_homer_cpeb234_vs_cpeb1_170120')
out <- find_motifs_genome(utr3.cpeb234,background=utr3.cpeb1,path=file.path(.reportsdir,'motifs_homer_cpeb234_vs_cpeb1_170120'),genome='xenLae92',scan_size='given',cores=10,only_known=FALSE,only_denovo=FALSE,overwrite=TRUE)

## Launch Homer, CPEB1 vs CPEB1234
dir.create(.reportsdir,'motifs_homer_cpeb1_vs_cpeb1234_170120')
out <- find_motifs_genome(utr3.cpeb1,background=utr3.cpeb1234,path=file.path(.reportsdir,'motifs_homer_cpeb1_vs_cpeb1234_170120'),genome='xenLae92',scan_size='given',cores=10,only_known=FALSE,only_denovo=FALSE,overwrite=TRUE)

## Launch Homer, CPEB234 vs CPEB1
dir.create(.reportsdir,'motifs_homer_cpeb234_vs_cpeb1234_170120')
out <- find_motifs_genome(utr3.cpeb234,background=utr3.cpeb1234,path=file.path(.reportsdir,'motifs_homer_cpeb234_vs_cpeb1234_170120'),genome='xenLae92',scan_size='given',cores=10,only_known=FALSE,only_denovo=FALSE,overwrite=TRUE)

# Save dfs for homer
save(utr3.cpeb1,utr3.cpeb234,utr3.cpeb1234,file=file.path(.datadir,'dfs_utr_cpebs_4Homer_jan20.RData'))

###########################################################################
## chunk 2. LOAD RESULTS AND GENERATE MOTIF ENRICHMENT PLOTS AGAINST WG
###########################################################################

rm(list=ls())

# Load known results
dirs <- list.files(.reportsdir,'motifs_homer',full.names=TRUE)
names(dirs) <- gsub('motifs_homer_','',gsub('_170120','',basename(dirs)))
dirs
motifs <- lapply(dirs,function(d) list(known=read_known_results(d),denovo=read_denovo_results(d)))

known <- lapply(dirs,function(d) read.table(file.path(d,'knownResults.txt'),sep='\t',comment.char='',header=TRUE))

names(known)

# Plot cpeb 1 vs 4
x <- motifs[['cpeb1_vs_bg']]$denovo
y <- motifs[['cpeb234_vs_bg']]$denovo
z <- motifs[['cpeb1_vs_cpeb234']]$denovo
w <- motifs[['cpeb234_vs_cpeb1']]$denovo

sel <- intersect(x[,1], y[,1])
df <- data.frame(cpeb1=x[match(sel, x[,1]), 5], cpeb234=y[match(sel, y[,1]), 5], motif=sub("/.*", "", x[match(sel, x[,1]), 1]), orig=x[match(sel, x[,1]), 1])

df$type <- rep("both", nrow(df))
sel <- (df$cpeb1 > 1 & df$cpeb234 < 0.1 )
df$type[sel] <- "cpeb234"
sel <- (df$cpeb1 < 0.1 & df$cpeb234 > 1 )
df$type[sel] <- "cpeb1"

library(ggplot2)
library(ggrepel)

pdf(file.path(.figsdir, "motifs_homer_cpeb1_vs_234_log10.pdf"), width=12, height=8)
q <- ggplot(aes(x=-log10(cpeb1), y=-log10(cpeb234), label=motif, colour=type), data=df)
sel <- df$type != "both"
# plot
q + geom_point() + geom_vline(xintercept = -log10(0.05), colour="#77777f", linetype = "longdash") + geom_hline(yintercept = -log10(0.05), colour="#77777f", linetype = "longdash") + geom_text_repel(data=df[sel,], aes(x=-log10(cpeb1), y=-log10(cpeb234), label=motif, colour=type), show.legend = FALSE, size = 4, box.padding = unit(0.25, 'lines'),point.padding = unit(1, 'lines'), segment.color = '#77777f', segment.size = 0.5, force = 2) + scale_colour_manual(breaks=c("both","cpeb1", "cpeb234"), values=c("black", "red", "green"))
dev.off()

###########################################################################
## chunk 3. SCAN TARGET UTRS FOR GIVEN CPE MOTIFS
###########################################################################

# Prepare motif files
motifs <- list(
    PAS=c('AATAAA','ATTAAA','AAGAAA'),
    CPE=c('TTTTAT','TTTTAAT','TTTTACT','TTTTAAAT','TTTTAAGT','TTTTCAT'))

dir.create(file.path(.tablesdir,'motifs_CPE_PAS'))
motifs.pwm <- lapply(names(motifs),function(i)
                     {
                         lapply(1:length(motifs[[i]]),function(j)
                                {
                                    motif <- motifs[[i]][[j]]
                                    print(sprintf('%s %s',i,motif))
                                    pwm <- as.data.frame(do.call(cbind,lapply(c('A','C','G','T'),function(nt) as.numeric(strsplit(motif,'')[[1]]==nt))))
                                    colnames(pwm) <- c('A','C','G','T')
                                    rownames(pwm) <- 1:nrow(pwm)
                                    as_tibble(pwm)
                                    write_homer_motif(as_tibble(pwm),motifs[[i]][[j]],log_odds_detection=0,file=file.path(.tablesdir,sprintf('motifs_CPE_PAS/motifs_%s.pwm',i)),append='TRUE')
                                })
                     })

# Run, will actually run for whole genome...
load(file=file.path(.datadir,'genes_3utr_xl92_jan20.RData'))
utr3.all <- genes.utr[,c('seqnames.y','start.y','end.y','GeneID')]

# PAS
resdir <- file.path(.reportsdir,'motifs_homer_PAS_290120')
dir.create(resdir)
find_motifs_instances(utr3.all,genome='xenLae92',path=file.path(.reportsdir,'motifs_homer_PAS_290120/utr3_xl92_PAS_match.txt'),motif_file=file.path(.tablesdir,'motifs_CPE_PAS/motifs_PAS.pwm'),scan_size='given',cores=6)

# CPE
resdir <- file.path(.reportsdir,'motifs_homer_CPE_290120')
dir.create(resdir)
find_motifs_instances(utr3.all,genome='xenLae92',path=file.path(.reportsdir,'motifs_homer_CPE_290120/utr3_xl92_CPE_match.txt'),motif_file=file.path(.tablesdir,'motifs_CPE_PAS/motifs_CPE.pwm'),scan_size='given',cores=6)

# Load results and keep only exact matches, PAS
match.pas <- read.delim(file.path(.reportsdir,'motifs_homer_PAS_290120/utr3_xl92_PAS_match.txt'),header=TRUE,as.is=TRUE)
match.pas2 <- match.pas[match.pas$Sequence==match.pas$Motif.Name,]
write.table(match.pas2,file.path(.reportsdir,'motifs_homer_PAS_290120/utr3_xl92_PAS_perfectmatch.txt'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE,dec='.')

# Load results and keep only exact matches, CPE
match.cpe <- read.delim(file.path(.reportsdir,'motifs_homer_CPE_290120/utr3_xl92_CPE_match.txt'),header=TRUE,as.is=TRUE)
match.cpe2 <- match.cpe[match.cpe$Sequence==match.cpe$Motif.Name,]
write.table(match.cpe2,file.path(.reportsdir,'motifs_homer_CPE_290120/utr3_xl92_CPE_perfectmatch.txt'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE,dec='.')

# Collapse and transpose
match.pas2.col <- aggregate(match.pas2[,c('Offset','Strand')],by=list(GeneID=match.pas2$PositionID,Motif=match.pas2$Motif.Name),FUN=function(x) paste(x,collapse=' / '))
rownames(match.pas2.col) <- paste(match.pas2.col$GeneID,match.pas2.col$Motif,sep='_')

# Integrate with our ref table
idx <- intersect(xout$GeneID,match.pas2.col$GeneID)
match.df <- do.call(cbind,mclapply(motifs$PAS,function(m) match.pas2.col[paste(idx,m,sep='_'),],mc.cores=6)) # probably faster with dplyr or smth like that...


### THE END

