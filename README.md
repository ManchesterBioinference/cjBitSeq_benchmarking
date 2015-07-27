# cjBitSeq_benchmarking
benchmarking experiments for cjBitSeq


The following software is required:

Spanki 0.5.0

cufflinks-2.1.1.Linux_x86_64

bowtie2-2.2.5

tophat-2.1.0.Linux_x86_64

BitSeq 0.7.5

rsem-1.2.15

cjBitSeq

Also boost libraries and samtools should be available on your system. 

1. CD to the drosophila/Annotation directory
2. Run ./makeAnnotation.sh which will extract and build the necessary reference annotations

The run the test example:

1. CD to the drosophila/testExample directory and run `R CMD BATCH groundTruth.R` which will simulate ground truth files.
2. Run `runMe.sh` in order to:

    2.a simulate reads with spanki
    
    2.b align reads with bowtie2
    
    3.c run cjBitSeq
    
    3.d run BitSeq stage1/stage2
    
    3.e run Cufflinks/Cuffdif
    
    3.f run RSEM/EBSeq
    
    3.g produce graphs in *.pdf format.
