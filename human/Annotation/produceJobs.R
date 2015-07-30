
conOut1 <- file("simulateReads.sh",open = "w")
for (i in c(paste("A",1:n1,sep=""),paste("B",1:n2,sep=""))){
	myCommand <- paste("spankisim_transcripts -o ",i," -g ../Annotation/genesNew.gtf -f ../Annotation/WholeGenomeFasta/genome.fa -bp 100 -t tr_File_",i,".tr -ends 2",sep="")
	cat(myCommand,"\n",file = conOut1,append=TRUE)
}
close(conOut1)


conOut1 <- file("alignReads.sh",open = "w")
for (i in c(paste("A",1:n1,sep=""),paste("B",1:n2,sep=""))){
	myCommand1 <- paste("cd ",i,sep="")
	myCommand2 <- "bowtie2 -q -k 100 --threads 4 --no-discordant -x ../../Annotation/transcriptome_data/known -1 sim_1.fastq -2 sim_2.fastq -S data.sam"
	myCommand3 <- "parseAlignment data.sam -o data.prob --trSeqFile ../../Annotation/transcriptome_data/known.fa --trInfoFile data.tr --uniform"
	myCommand4 <- "R CMD BATCH ../../Annotation/getTranscriptNames.R"
	myCommand5 <- "cd .."
	cat(myCommand1,"\n",file = conOut1,append=TRUE)
	cat(myCommand2,"\n",file = conOut1,append=TRUE)
	cat(myCommand3,"\n",file = conOut1,append=TRUE)
	cat(myCommand4,"\n",file = conOut1,append=TRUE)
	cat(myCommand5,"\n\n",file = conOut1,append=TRUE)
}
close(conOut1)

conOut1 <- file("cjBitSeq.sh",open = "w")
l1 <- paste("A",1:n1,"/data.prob",sep="",collapse = " ")
l2 <- paste("B",1:n1,"/data.prob",sep="",collapse = " ")
commandLine <- paste(l1,l2,sep = " C ")
myCommand <- paste("cjBitSeq cjBitSeqOutput ",commandLine,sep="")
cat(myCommand,"\n",file = conOut1,append=TRUE)
close(conOut1)

conOut1 <- file("BitSeq.sh",open = "w")
for (i in c(paste("A",1:n1,sep=""),paste("B",1:n2,sep=""))){
	myCommand1 <- paste("cd ",i,sep="")
	myCommand2 <- "estimateExpression data.prob -o data --outType RPKM -p data.txt -t data.tr -P 4"
	cat(myCommand1,"\n",file = conOut1,append=TRUE)
	cat(myCommand2,"\n",file = conOut1,append=TRUE)
	myCommand1 <- "cd ../"
	cat(myCommand1,"\n",file = conOut1,append=TRUE)
}
cat("mkdir bitseq","\n",file = conOut1,append=TRUE)
cat("cd bitseq","\n",file = conOut1,append=TRUE)
l1 <- paste("../A",1:n1,"/data.rpkm",sep="",collapse = " ")
l2 <- paste("../B",1:n1,"/data.rpkm",sep="",collapse = " ")
commandLine <- paste(l1,l2,sep = " ")
myCommand <- paste("getVariance --log -o data.Lmean ",commandLine,sep="")
cat(myCommand,"\n",file = conOut1,append=TRUE)
commandLine <- paste(l1,l2,sep = " C ")
myCommand <- paste("estimateHyperPar --lowess-f=0.2 --meanFile data.Lmean -o data.param ",commandLine,sep="") #0.2 is the default
cat(myCommand,"\n",file = conOut1,append=TRUE)
myCommand <- paste("estimateDE -o data -p data.param ",commandLine,sep="") #0.2 is the default
cat(myCommand,"\n",file = conOut1,append=TRUE)
close(conOut1)

conOut1 <- file("cufflinks.sh",open = "w")
cat("mkdir cufflinks","\n",file = conOut1,append=TRUE)
cat("cd cufflinks","\n",file = conOut1,append=TRUE)	
for (i in c(paste("A",1:n1,sep=""),paste("B",1:n2,sep=""))){
	myCommand1 <- paste("tophat2 -o ",i,".tophat.out --no-novel-juncs -g 100 -i 20 -m 2 --microexon-search --min-anchor-length 3 --read-mismatches 10 --read-gap-length 10 --read-edit-dist 10 -x 100 -T -p 4 --transcriptome-index=../../Annotation/transcriptome_data/known ../../Annotation/WholeGenomeFasta/genome ../",i,"/sim_1.fastq ../",i,"/sim_2.fastq",sep="")
	cat(myCommand1,"\n",file = conOut1,append=TRUE)
	myCommand1 <- paste("cufflinks -q -o ",i,".cufflinks.out -u -p 4 -G ../../Annotation/genes.gtf ",i,".tophat.out/accepted_hits.bam",sep="")
	cat(myCommand1,"\n",file = conOut1,append=TRUE)
}
#cat("cp ../../Annotation/deleteDuplicates.sh deleteDuplicates.sh","\n\n",file = conOut1,append=TRUE)
#cat("./deleteDuplicates.sh","\n\n",file = conOut1,append=TRUE)
l1 <- paste("A",1:n1,".cufflinks.out/transcripts.gtf",sep="",collapse = " ")
l2 <- paste("B",1:n1,".cufflinks.out/transcripts.gtf",sep="",collapse = " ")
commandLine <- paste(l1,l2,sep = " ")
myCommand <- paste("cuffcompare -o cuff_compare ",commandLine,sep="")
cat(myCommand,"\n",file = conOut1,append=TRUE)
l1 <- paste("A",1:n1,".tophat.out/accepted_hits.bam",sep="",collapse = ",")
l2 <- paste("B",1:n1,".tophat.out/accepted_hits.bam",sep="",collapse = ",")
commandLine <- paste(l1,l2,sep = " ")
myCommand <- paste("cuffdiff -o cuffdif_out -p 4 cuff_compare.combined.gtf ",commandLine,sep="")
cat(myCommand,"\n",file = conOut1,append=TRUE)
cat("R CMD BATCH ../../Annotation/matchNames.R","\n",file = conOut1,append=TRUE)
close(conOut1)

conOut1 <- file("rsem.sh",open = "w")
cat("mkdir rsem","\n",file = conOut1,append=TRUE)
cat("cd rsem","\n",file = conOut1,append=TRUE)
for (i in c(paste("A",1:n1,sep=""),paste("B",1:n2,sep=""))){
	myCommand1 <- paste("rsem-calculate-expression -p 4 --bowtie2 --paired-end --no-bam-output ../",i,"/sim_1.fastq ../",i,"/sim_2.fastq ../../Annotation/rsemReference/ref expression",i,sep="")
	cat(myCommand1,"\n",file = conOut1,append=TRUE)
}
cat("R CMD BATCH ../ebSeq.R","\n",file = conOut1,append=TRUE)
close(conOut1)

conOut1 <- file("ebSeq.R",open = "w")
myCommand1 <- "library('EBSeq')"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
myCommand1 <- "rawData <- read.table('expressionA1.isoforms.results',header = TRUE)[,5]"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
if(n1 > 1){
	for (i in 2:n1){
		myFile <- paste("read.table('expressionA",i,".isoforms.results',header = TRUE)[,5]",sep="")
		myCommand1 <- paste("rawData <- cbind(rawData,",myFile,")",sep="")
		cat(myCommand1,"\n",file = conOut1,append=TRUE)
	}
}
for (i in 1:n2){
	myFile <- paste("read.table('expressionB",i,".isoforms.results',header = TRUE)[,5]",sep="")
	myCommand1 <- paste("rawData <- cbind(rawData,",myFile,")",sep="")
	cat(myCommand1,"\n",file = conOut1,append=TRUE)
}
myCommand1 <- "K <- 48009"
myCommand2 <- "myFile <- file('../../Annotation/FastaNames.txt',open='r')"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
cat(myCommand2,"\n",file = conOut1,append=TRUE)
myCommand1 <- "sailNames <- array(data = NA, dim = c(K,2))"
myCommand2 <- "i <- 0"
myCommand3 <- "while (length(oneLine <- readLines(myFile, n = 1, warn = FALSE)) > 0){i <- i + 1"
myCommand4 <- "sailNames[i,1] <- as.numeric(strsplit(strsplit(oneLine,split = ' ')[[1]][1],split = '>')[[1]][2])"
myCommand5 <- "sailNames[i,2] <- strsplit(oneLine,split = ' ')[[1]][2]}"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
cat(myCommand2,"\n",file = conOut1,append=TRUE)
cat(myCommand3,"\n",file = conOut1,append=TRUE)
cat(myCommand4,"\n",file = conOut1,append=TRUE)
cat(myCommand5,"\n",file = conOut1,append=TRUE)
myCommand1 <- "close(myFile)"
myCommand2 <- "A1 <- read.table('expressionA1.isoforms.results',header=TRUE)"
myCommand3 <- "groundTruth <- read.table('../A1/transcriptNames.txt')"
myCommand4 <- "rawData2 <- rawData"
myCommand5 <- "genes <- read.table('../../Annotation/geneTranscriptNames.txt',fill=TRUE)"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
cat(myCommand2,"\n",file = conOut1,append=TRUE)
cat(myCommand3,"\n",file = conOut1,append=TRUE)
cat(myCommand4,"\n",file = conOut1,append=TRUE)
cat(myCommand5,"\n",file = conOut1,append=TRUE)
myCommand1 <- "for (i in 1:dim(groundTruth)[1]){j <- which(sailNames[,2] == as.character(genes[i,2]));if(length(j)>0){k <- which(A1[,1] == sailNames[j,1]);rawData2[i,] <- rawData[k,]}else{"
myCommand2 <- "rawData2[i,] <- rep(0,dim(rawData2)[2])}"
myCommand3 <- "if(i%%1000 == 0){print(i)}}"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
cat(myCommand2,"\n",file = conOut1,append=TRUE)
cat(myCommand3,"\n",file = conOut1,append=TRUE)
myCommand1 <- "IsoMat <- rawData2"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
myCommand1 <- "rownames(IsoMat) <- genes[,2]"
myCommand2 <- "IsoNames <- genes[,2]"
myCommand3 <- "IsosGeneNames <- genes[,1]"
myCommand4 <- "IsoSizes=MedianNorm(IsoMat)"
myCommand5 <- "NgList=GetNg(IsoNames, IsosGeneNames)"
myCommand6 <- "IsoNgTrun=NgList$IsoformNgTrun"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
cat(myCommand2,"\n",file = conOut1,append=TRUE)
cat(myCommand3,"\n",file = conOut1,append=TRUE)
cat(myCommand4,"\n",file = conOut1,append=TRUE)
cat(myCommand5,"\n",file = conOut1,append=TRUE)
cat(myCommand6,"\n",file = conOut1,append=TRUE)
myCommand1 <- paste("IsoEBOut=EBTest(Data=IsoMat, NgVector=IsoNgTrun, Conditions=as.factor(rep(c('C1','C2'),c(",n1,",",n2,"))),sizeFactors=IsoSizes, maxround=10)",sep="")
myCommand2 <- "IsoPP=GetPPMat(IsoEBOut)"
myCommand3 <- "ebSeqDe <- as.numeric(IsoPP[,'PPDE'])"
myCommand4 <- "ebSeqDe2 <- numeric(K)"
myCommand5 <- "tested <- which(is.na(match(IsoNames,rownames(IsoPP)))==FALSE)"
myCommand6 <- "ebSeqDe2[tested] <- ebSeqDe"
myCommand7 <- "write.table(ebSeqDe2,file = 'ebSeqDe.txt')"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
cat(myCommand2,"\n",file = conOut1,append=TRUE)
cat(myCommand3,"\n",file = conOut1,append=TRUE)
cat(myCommand4,"\n",file = conOut1,append=TRUE)
cat(myCommand5,"\n",file = conOut1,append=TRUE)
cat(myCommand6,"\n",file = conOut1,append=TRUE)
cat(myCommand7,"\n",file = conOut1,append=TRUE)


close(conOut1)



conOut1 <- file("runMe.sh",open = "w")
myCommand1 <- "./simulateReads.sh"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
myCommand1 <- "./alignReads.sh"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
myCommand1 <- "./cjBitSeq.sh"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
myCommand1 <- "./BitSeq.sh"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
myCommand1 <- "./cufflinks.sh"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
myCommand1 <- "./rsem.sh"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
myCommand1 <- "R CMD BATCH ../Annotation/graphics.R"
cat(myCommand1,"\n",file = conOut1,append=TRUE)
close(conOut1)

system("chmod u+x runMe.sh")
system("chmod u+x simulateReads.sh")
system("chmod u+x alignReads.sh")
system("chmod u+x cjBitSeq.sh")
system("chmod u+x BitSeq.sh")
system("chmod u+x cufflinks.sh")
system("chmod u+x rsem.sh")







