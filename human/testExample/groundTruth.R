# load names
txid <- read.table("../Annotation/trNames.tr")[,1]
trLengths <- read.table("../Annotation/data.tr")
smallTranscripts <- which(trLengths[,4] < 200)
K <- length(txid)
mus0 <- rep(0,K)
activeTranscripts <- 1:K
activeTranscripts <- activeTranscripts[-smallTranscripts]
Kactive <- length(activeTranscripts)
set.seed(1)
expressedIndex <- activeTranscripts[sample(Kactive,10000,replace = FALSE)]
nDE <- 1000
deIndex <- activeTranscripts[sample(Kactive,nDE,replace = FALSE)]
mus0[expressedIndex] <- rep(65,length(expressedIndex))

#conditionA1
mus <- mus0
mus[deIndex[1:(nDE/2)]] <- rep(20,nDE/2)
mus[deIndex[(1+nDE/2):nDE]] <- rep(100,nDE/2)
rpk <- rpois(K,mus)
tr_File_1 <- data.frame(txid,rpk)
write.table(tr_File_1,file = "tr_File_A1.tr",quote = FALSE,row.names = FALSE,sep = "\t")
#conditionA2
mus <- mus0
mus[deIndex[1:(nDE/2)]] <- rep(20,nDE/2)
mus[deIndex[(1+nDE/2):nDE]] <- rep(100,nDE/2)
rpk <- rpois(K,mus)
tr_File_1 <- data.frame(txid,rpk)
write.table(tr_File_1,file = "tr_File_A2.tr",quote = FALSE,row.names = FALSE,sep = "\t")

#conditionB1
mus <- mus0
mus[deIndex[1:(nDE/2)]] <- rep(100,nDE/2)
mus[deIndex[(1+nDE/2):nDE]] <- rep(20,nDE/2)
rpk <- rpois(K,mus)
tr_File_2 <- data.frame(txid,rpk)
write.table(tr_File_2,file = "tr_File_B1.tr",quote = FALSE,row.names = FALSE,sep = "\t")
#conditionB2
mus <- mus0
mus[deIndex[1:(nDE/2)]] <- rep(100,nDE/2)
mus[deIndex[(1+nDE/2):nDE]] <- rep(20,nDE/2)
rpk <- rpois(K,mus)
tr_File_2 <- data.frame(txid,rpk)
write.table(tr_File_2,file = "tr_File_B2.tr",quote = FALSE,row.names = FALSE,sep = "\t")


colInd <- rep(1,K)
perm <- deIndex
colInd[perm] <- rep(2,nDE)
write.table(colInd - 1,file = "realDE.txt")



#produce jobScripts:
n1 <- 2 #replicates for 1st condition
n2 <- 2 #replicates for 2nd condition
source(file = "../Annotation/produceJobs.R")

