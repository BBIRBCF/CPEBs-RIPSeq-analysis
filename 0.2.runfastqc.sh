/software/FastQC_011/fastqc ../data/sequences/fastq/*.fastq -f fastq -t 16 --noextract -o ../reports/fastqc/
/software/FastQC_011/fastqc ../data/sequences/bam/*.fq -f fastq -t 16 --noextract -o ../reports/fastqc/
/software/FastQC_011/fastqc ../data/sequences/bam/*.bam -f bam_mapped -t 16 --noextract -o ../reports/fastqc/
