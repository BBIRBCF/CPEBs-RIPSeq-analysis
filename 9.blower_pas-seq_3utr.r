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
# 1. READ PASSEQ AND EXPORT 3'UTR BED FILES FOR EVALUATION
###########################################################################

list.files(file.path(datadir,'yang_pas-seq'))
utr <- read.delim(file.path(datadir,'yang_pas-seq/2020_blower-et-al-PAS-Seq-calls_Supplemental_Table_S11.txt'),as.is=TRUE)

# check with our genome
load(file.path(datadir,'genes_3utr_xl92_jan20.RData')) # genes.utr, x is for gene, y is for UTR
table(utr$GeneID %in% genes.utr$GeneID) # 1351 not present in our genome...
table(genes.utr$GeneID %in% utr$GeneID)

# Export all new UTRs
utr.new <- data.frame(seqnames=utr$Chr,
                      start=as.numeric(unlist(lapply(strsplit(utr[,3],':'),function(x) x[2]))),
                      end=as.numeric(unlist(lapply(strsplit(utr[,3],':'),function(x) x[3]))),
                      Name=paste(utr$GeneID,utr$Number.of.reads,utr$Gene.strand,sep='_'),
                      GeneID=utr$GeneID,
                      Nreads=utr$Number.of.reads,utr$Gene.strand,sep='_')
utr.new$width <- abs(utr.new$start - utr.new$end)

# All
write.table(utr.new[,1:4],file=file.path(tablesdir,'blower_xl92_passeq_utr_all.bed'),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)

# Longest
utr.new <- utr.new[order(utr.new$GeneID,utr.new$width,decreasing=TRUE),]
write.table(utr.new[!duplicated(utr.new$GeneID),1:4],file=file.path(tablesdir,'blower_xl92_passeq_utr_longest.bed'),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)

# Most reads
utr.new <- utr.new[order(utr.new$GeneID,utr.new$Nreads,decreasing=TRUE),]
write.table(utr.new[!duplicated(utr.new$GeneID),1:4],file=file.path(tablesdir,'blower_xl92_passeq_utr_mostreads.bed'),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)

### THE END




