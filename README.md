# cjBitSeq_benchmarking
benchmarking experiments for cjBitSeq


The following software is required:

* Spanki 0.5.0
* cufflinks-2.1.1.Linux_x86_64
* bowtie2-2.2.5
* tophat-2.1.0.Linux_x86_64
* BitSeq 0.7.5
* rsem-1.2.15
* cjBitSeq
* R

Also boost libraries and samtools should be available on your system. Required R libraries: 

* Matrix
* foreach
* doMC
* EBSeq
* ROCR
* pracma
* fields

1. CD to the drosophila/Annotation directory
2. Run ./makeAnnotation.sh which will extract and build the necessary reference annotations

To run the test example:

1. CD to the drosophila/testExample directory and run `R CMD BATCH groundTruth.R` which will simulate ground truth files.
2. Run `runMe.sh` in order to:
* simulate reads with spanki
*  align reads with bowtie2
*  run cjBitSeq
*  run BitSeq stage1/stage2
*  run Cufflinks/Cuffdif
*  run RSEM/EBSeq
*  produce graphs in *.pdf format.
