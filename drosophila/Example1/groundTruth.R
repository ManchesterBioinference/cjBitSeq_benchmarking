# load names
txid <- read.table("../Annotation/trNames.tr")[,1]
K <- length(txid)

#conditionA1
mus <- rep(65,K)
mus[1001:2000] <- rep(20,1000)
mus[10001:11000] <- rep(100,1000)
rpk <- rpois(K,mus)
tr_File_1 <- data.frame(txid,rpk)
write.table(tr_File_1,file = "tr_File_A1.tr",quote = FALSE,row.names = FALSE,sep = "\t")
#conditionA2
mus <- rep(65,K)
mus[1001:2000] <- rep(20,1000)
mus[10001:11000] <- rep(100,1000)
rpk <- rpois(K,mus)
tr_File_1 <- data.frame(txid,rpk)
write.table(tr_File_1,file = "tr_File_A2.tr",quote = FALSE,row.names = FALSE,sep = "\t")

#conditionB1
mus <- rep(65,K)
mus[1001:2000] <- rep(100,1000)
mus[10001:11000] <- rep(20,1000)
rpk <- rpois(K,mus)
tr_File_2 <- data.frame(txid,rpk)
write.table(tr_File_2,file = "tr_File_B1.tr",quote = FALSE,row.names = FALSE,sep = "\t")
#conditionB2
mus <- rep(65,K)
mus[1001:2000] <- rep(100,1000)
mus[10001:11000] <- rep(20,1000)
rpk <- rpois(K,mus)
tr_File_2 <- data.frame(txid,rpk)
write.table(tr_File_2,file = "tr_File_B2.tr",quote = FALSE,row.names = FALSE,sep = "\t")

nDE <- 2000
colInd <- rep(1,K)
perm <- c(1001:2000,10001:11000)
colInd[perm] <- rep(2,nDE)
write.table(colInd - 1,file = "realDE.txt")



#produce jobScripts:
n1 <- 2 #replicates for 1st condition
n2 <- 2 #replicates for 2nd condition
conOut1 <- file("simulateReads.sh",open = "w")
for (i in c(paste("A",1:n1,sep=""),paste("B",1:n2,sep=""))){
	myCommand <- paste("spankisim_transcripts -o ",i," -g ../Annotation/genesNew.gtf -f ../Annotation/WholeGenomeFasta/genome.fa -bp 76 -t tr_File_",i,".tr -ends 2",sep="")
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




