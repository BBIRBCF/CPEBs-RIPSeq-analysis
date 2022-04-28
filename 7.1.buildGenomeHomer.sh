/software/cufflinks-2.2.1.Linux_x86_64/gffread XENLA_9.2_Xenbase.gff3 -T -o XENLA_9.2_Xenbase.gtf
perl /software/homer/bin/loadGenome.pl -name xenLae92 -org laevis -fasta ./XL9_2.fa -gtf ./XENLA_9.2_Xenbase.gtf -force

