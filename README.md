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

We provide drosophila and human examples. The following steps apply to the drosophila annotation but they are essentially the same for the human annotation. 

1. CD to the drosophila/Annotation directory
2. Run ./makeAnnotation.sh which will extract and build the necessary reference annotations (it needs to be done only once).

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

The graphs include the following:

* SAR curve
* ROC curve
* precision-recall curve
* power-to-achieved FDR plots
* venn graph of DE transcripts at the 0.05 level
* scatterplot of the true relative transcript expression (in log-scale) for both conditions coloured according to the DE evidence of each method.


SAR-curve tends to be the most informative plot: it combines three different performance measures, that is, (a) Accuracy, (b) Area under Curve and (c) Root mean square Error, as a function of the cut-off threshold (which is controlled by the user). The highlighted area of the SAR-curve corresponds to the range of the cut-off values used in practice (that is, \alpha < 0.1).


