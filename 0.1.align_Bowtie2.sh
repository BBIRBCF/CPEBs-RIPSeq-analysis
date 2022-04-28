
# Binary paths...
#####################
fastqc='/software/FastQC_011/fastqc'
sambamba='/software/sambamba'
bowtie2='/software/bowtie2-2.2.2/bowtie2'
igvtools='/Volumes/biostats/soft/IGVTools2/igvtools'

# Genomes
#####################
genomebowtie2='/Volumes/biostats/databases/bowtie2_indexes/xenLae2_061118'
genomeigv='xenLae2_061118'

# Procs
#####################
nprocfastqc=8
nprocbowtie=16
nprocsambamba=8

# Data paths
#####################
fastqdir='/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/sequences/fastq'
bamdir='/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/data/sequences/bam'
tdfdir='/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/tables/files_4_igv'
fastqcdir='/Volumes/biostats/consulting/raul_mendez/bduran_201903_ripseq/reports/fastqc'

#################################################################################################

# Main thing...

mkdir $bamdir

cd $fastqdir

for file in *.fastq
do {

   # Bowtie2 alignment
   $bowtie2 -q -N 1 -p $nprocbowtie --phred33 -x $genomebowtie2 -U ./$file --no-unal -S $bamdir/$file.sam --un $bamdir/$file.un.fq

   cd $bamdir

   # sambamba SAM to BAM and sort
   $sambamba view  -S -f bam -t $nprocsambamba -p $file.sam -o $file.bam
   $sambamba sort -t $nprocsambamba -p $file.bam
   # rename sorted files and keep sorted only
   rm $file.sam $file.bam
   mv $file.sorted.bam $file.bam
   mv $file.sorted.bam.bai $file.bam.bai

   # Generation of TDF tracks for IGV of aligned bam
   $igvtools count -z 7 -w 25 -e 0 $file.bam $file.tdf $genomeigv

   # [Optional] Picard based filtering and removal of duplicates using sambamba
   $sambamba markdup -p -r -t $nprocsambamba $file.bam $file.filterdup.bam 

   # [Optional] Generation of TDF tracks for IGV of filterdup aligned bam
   $igvtools count -z 7 -w 25 -e 0 $file.filterdup.bam $file.filterdup.tdf $genomeigv

   # Move tdfs
   mv *.tdf $tdfdir

   cd $fastqdir

} done
