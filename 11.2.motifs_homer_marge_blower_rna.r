##########################################################################
## PROJECT: BDURAN RIPSEQ
## DESCRIPTION OF FILE: MOTIF SEARCHING VIA HOMER+MARGE FOR BLOWER UTRs
## DATE: JUN 2020
###########################################################################

library(marge)
library(openxlsx)
.datadir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/'
.datadir2 <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/'
.datadir1 <- '/Volumes/biostats/consulting/raul_mendez/bduran_201810_ripseq/data/'
.reportsdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports'
.tablesdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables'
.tdfdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables/files_4_igv'
.figsdir <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/figs'

###########################################################################
## chunk 1. TEST TO RUN HOMER DIRECTLY
###########################################################################

## taken from find_motifs_genome function
##system(paste("mkdir -p", path))
##homer_base <- get <- homer <- bin()
##cmd <- paste(paste0(homer <- base, "findMotifsGenome.pl"), target <- bed,
##             genome, path, "-len", paste0(motif <- length, collapse = ","),
##             "-size", scan <- size, "-S", optimize <- count, "-p", cores,
##             "-cache", cache, "-fdr", fdr <- num)

## our parameters
path <- file.path(.reportsdir,'homer_output_blower_rna')
system(paste("mkdir -p", path))
homer_base <- get_homer_bin()
fasta <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/all_fasta.fa'
genome <- '/Volumes/biostats/databases/RefGenome/xlaevis/xlaevis_92_ucsc/xenLae2_061118.fa'
motif_length <- c(8,10,12)
scan_size <- 'given' ## equals option -chopify
optimize_count <- 8
fdr_num <- 0
cores <- 10

run_homer_fasta <- function(homer_base, fasta, path, genome, motif_length=c(8,10,12), optimize_count=9, fdr_num=0, cores=10, run=TRUE, ...)
    {
        cmd <- paste(paste0(homer_base, "findMotifs.pl"), fasta, "fasta",
                     path, "-fasta", genome, "-len", paste0(motif_length, collapse = ","),
                     "-chopify", "-S", optimize_count, "-p", cores,
                     "-fdr", fdr_num)
        if (run) system(cmd,...)
        cmd
    }

##out <- system(cmd,intern=FALSE,wait=FALSE,show.output.on.console=FALSE)

## It works !

###########################################################################
## chunk 1. RUN HOMER FOR THE DIFFERENT TARGET LISTS
###########################################################################

rm(list=ls())
library(parallel)

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
length(rownames(targets)[rowSums(targets[,1:4])>0]) # 1799

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

# Prepare targets table for interim use
targets.xout <- targets
targets.xout$GeneID <- rownames(targets)
targets.xout <- merge(targets.xout,genes.utr,by.x='GeneID',by.y='GeneID',all.x=TRUE,all.y=FALSE)
colnames(targets.xout) <- gsub('.x','.tx_of_3utr',colnames(targets.xout),fixed=TRUE)
colnames(targets.xout) <- gsub('.y','.3utr',colnames(targets.xout),fixed=TRUE)
rownames(xout) <- xout$GeneID

targets.xout <- cbind(xout[targets.xout$GeneID,1:5],targets.xout)
unique(targets.xout[,1]==targets.xout[,6])
##write.xlsx(targets.xout[,-6],file=file.path(.tablesdir,'xl92_genes_longest3UTR_cpebtargets_jan20.xlsx'))

# Split for homer UP (CPEB1) or DOWN (CPEB234)
targets.bak <- targets
targets <- targets[targets$hasUTR3==1,]
sel.homer.cpeb1 <- targets$rej.CPEB1_vs_CPEB234==1
sel.homer.cpeb234 <- targets$rej.CPEB1_vs_CPEB234==-1
sel.homer.cpeb1234 <- targets$CPEB1234==1

# Generate dfs
utr3.cpeb1 <- genes.utr[rownames(targets)[sel.homer.cpeb1],c('GeneID')]
utr3.cpeb234 <- genes.utr[rownames(targets)[sel.homer.cpeb234],c('GeneID')]
utr3.cpeb1234 <- genes.utr[rownames(targets)[sel.homer.cpeb1234],c('GeneID')]

## Check homer
check_homer()

## Load Blower 3'UTRs fasta, re-export for target selection
library(Biostrings)
fasta <- readDNAStringSet(file.path(.datadir,'blower_pas_seq/most_reads_fasta.fa'))

## Subset fastas
fasta.cpeb1 <- fasta[names(fasta) %in% utr3.cpeb1]
fasta.cpeb234 <- fasta[names(fasta) %in% utr3.cpeb234]
fasta.cpeb1234 <- fasta[names(fasta) %in% utr3.cpeb1234]

## Save
writeXStringSet(fasta.cpeb1,filepath=file.path(.datadir,'blower_pas_seq/most_reads_fasta_cpeb1.fa'))
writeXStringSet(fasta.cpeb234,filepath=file.path(.datadir,'blower_pas_seq/most_reads_fasta_cpeb234.fa'))
writeXStringSet(fasta.cpeb1234,filepath=file.path(.datadir,'blower_pas_seq/most_reads_fasta_cpeb1234.fa'))

## Redefine homer function to call options with fasta
## our global parameters
path <- file.path(.reportsdir,'homer_output_blower_rna')
homer_base <- get_homer_bin()
motif_length <- c(8,10,12)
scan_size <- 'given' ## equals option -chopify
optimize_count <- 8
fdr_num <- 0
cores <- 10

run_homer_fasta <- function(homer_base, fasta, path, background, motif_length=c(8,10,12), optimize_count=9, fdr_num=0, cores=20, run=TRUE, ...)
    {
        system(paste("mkdir -p", path))
        cmd <- paste(paste0(homer_base, "findMotifs.pl"), fasta, "fasta",
                     path, "-fasta", background, "-len", paste0(motif_length, collapse = ","),
                     "-chopify", "-rna", "-S", optimize_count, "-p", cores,
                     "-fdr", fdr_num)
        if (run) system(cmd,...)
        cmd
    }

## Launch Homer, CPEB1234 vs whole BG (targets of at least one of them)
fasta <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb1234.fa'
background <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta.fa'
path <- file.path(.reportsdir,'motifs_homer_blower_cpeb1234_vs_bg_120620_rna')
out <- run_homer_fasta(homer_base, fasta, path, background, run=TRUE)

## Launch Homer, CPEB1 vs All UTRs
fasta <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb1.fa'
background <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta.fa'
path <- file.path(.reportsdir,'motifs_homer_blower_cpeb1_vs_bg_120620_rna')
out <- run_homer_fasta(homer_base, fasta, path, background, run=TRUE)

## Launch Homer, CPEB234 vs All UTRs
fasta <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb234.fa'
background <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta.fa'
path <- file.path(.reportsdir,'motifs_homer_blower_cpeb234_vs_bg_120620_rna')
out <- run_homer_fasta(homer_base, fasta, path, background, run=TRUE)

## Launch Homer, CPEB1 vs CPEB234
fasta <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb1.fa'
background <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb234.fa'
path <- file.path(.reportsdir,'motifs_homer_blower_cpeb1_vs_cpeb234_120620_rna')
out <- run_homer_fasta(homer_base, fasta, path, background, run=TRUE)

## Launch Homer, CPEB234 vs CPEB1
fasta <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb234.fa'
background <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb1.fa'
path <- file.path(.reportsdir,'motifs_homer_blower_cpeb234_vs_cpeb1_120620_rna')
out <- run_homer_fasta(homer_base, fasta, path, background, run=TRUE)

## Launch Homer, CPEB1 vs CPEB1234
fasta <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb1.fa'
background <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb1234.fa'
path <- file.path(.reportsdir,'motifs_homer_blower_cpeb1_vs_cpeb1234_120620_rna')
out <- run_homer_fasta(homer_base, fasta, path, background, run=TRUE)

## Launch Homer, CPEB1234 vs CPEB1
fasta <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb1234.fa'
background <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta_cpeb1.fa'
path <- file.path(.reportsdir,'motifs_homer_blower_cpeb1234_vs_cpeb1_120620_rna')
out <- run_homer_fasta(homer_base, fasta, path, background, run=TRUE)

###########################################################################
## chunk 2. SCAN TARGET UTRS FOR GIVEN CPE MOTIFS
###########################################################################

# Prepare motif files, this is done in a previous script, not needed
##motifs <- list(
##    PAS=c('AATAAA','ATTAAA','AAGAAA'),
##    CPE=c('TTTTAT','TTTTAAT','TTTTACT','TTTTAAAT','TTTTAAGT','TTTTCAT'))

##dir.create(file.path(.tablesdir,'motifs_CPE_PAS'))
##motifs.pwm <- lapply(names(motifs),function(i)
##                     {
##                         lapply(1:length(motifs[[i]]),function(j)
##                                {
##                                    motif <- motifs[[i]][[j]]
##                                    print(sprintf('%s %s',i,motif))
##                                    pwm <- as.data.frame(do.call(cbind,lapply(c('A','C','G','T'),function(nt) as.numeric(strsplit(motif,'')[[1]]==nt))))
##                                    colnames(pwm) <- c('A','C','G','T')
##                                    rownames(pwm) <- 1:nrow(pwm)
##                                    as_tibble(pwm)
##                                    write_homer_motif(as_tibble(pwm),motifs[[i]][[j]],log_odds_detection=0,file=file.path(.tablesdir,sprintf('motifs_CPE_PAS/motifs_%s.pwm',i)),append='TRUE')
##                                })
##                     })

## Usage
##findMotifs.pl targets.fa fasta motifResults/ -find <motif file> (This will cause the program to only scan for motifs)
## findMotifs.pl <targetSequences.fa> fasta <output directory> -fasta <background.fa> [options] -find motif1.motif > outputfile.txt

## taken from find_motifs_instances
run_homer_fasta_motifs <- function(homer_base, fasta, path, background, motif_file, cores=10, run=TRUE, ...)
    {
        dir.create(path)
        cmd <- paste(paste0(homer_base, "findMotifs.pl"), fasta, "fasta", path, "-p", cores, "-find", motif_file, ">", outfile)
        if (run) system(cmd,...)
        cmd
    }
## Launch Homer, PAS motifs
fasta <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta.fa'
background <- '/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/blower_pas_seq/most_reads_fasta.fa'
motif_file <- file.path(.tablesdir,'motifs_CPE_PAS/motifs_PAS.pwm')
path <- file.path(.reportsdir,'motifs_homer_blower_PAS')
outfile <- file.path(.reportsdir,'motifs_homer_blower_PAS.txt')

###############
## pending...
###############

out <- run_homer_fasta_motifs(homer_base, fasta, path, background, motif_file, run=TRUE)

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


### THE END...

