library(openxlsx)

## Read pblat results
out <- read.delim(file.path('../data/blower_pas_seq/output_blat_all_fasta_noheader.psl'),as.is=TRUE,header=FALSE) ## no header file in case read.psl worked, bad luck.

## colnames
cnames <- readLines(file.path('../data/blower_pas_seq/output.psl'),n=5)[3:4]
cnames1 <- unlist(strsplit(cnames[1],'\t'))
cnames2 <- unlist(strsplit(cnames[2],'\t'))
cnames <- gsub(' ','',paste(cnames1,cnames2))
cnames[20:21] <- gsub('match','',cnames[20:21])
colnames(out) <- make.names(cnames)

## filter by score, Score is calculated by matches-misMatches-qBaseInsert-tBaseInsert.
out$Score <- out$match - out$mis.match - out$Qgapbases - out$Tgapbases
out <- out[order(out$Qname,out$Score,decreasing=TRUE),]
out.fil <- out[!duplicated(out$Qname),]

## now filter to keep highest score per gene
out.fil$GeneID <- unlist(lapply(strsplit(out.fil$Qname,'_'),function(x) x[1]))
out.fil2 <- out.fil[order(out.fil$GeneID,out.fil$Score,decreasing=TRUE),]
out.fil2 <- out.fil2[!duplicated(out.fil2$GeneID),]

## Save Excels and BEDs
write.xlsx(out.fil,file.path('../data/blower_pas_seq/output_blat_all_fasta_best_score_per_UTR.xlsx'))
write.xlsx(out.fil2,file.path('../data/blower_pas_seq/output_blat_all_fasta_best_score_per_GeneID.xlsx'))

write.xlsx(out.fil,file.path('../reports/tables/blower_pas_seq_output_blat_all_fasta_best_score_per_UTR.xlsx'))
write.xlsx(out.fil2,file.path('../reports/tables/blower_pas_seq_output_blat_all_fasta_best_score_per_GeneID.xlsx'))

write.table(out.fil[,c('Tname','Tstart','Tend','Qname')],file.path('../reports/tables/blower_pas_seq_output_blat_all_fasta_best_score_per_UTR.bed'),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)
write.table(out.fil2[,c('Tname','Tstart','Tend','Qname')],file.path('../reports/tables/blower_pas_seq_output_blat_all_fasta_best_score_per_GeneID.bed'),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)

### THE END
