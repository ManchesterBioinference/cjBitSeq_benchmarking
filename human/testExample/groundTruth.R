# load names
txid <- read.table("../Annotation/trNames.tr")[,1]
trLengths <- read.table("../Annotation/data.tr")
smallTranscripts <- which(trLengths[,3] < 200)
K <- length(txid)

smalls <- which(l<200)
probs = rep(1/(K - length(smalls)),K)
probs[smalls] <- rep(0,length(smalls))
mus0 <- rep(0,K)
set.seed(1)
expressedIndex <- sample(K,10000,replace = FALSE,prob = probs)
nDE <- 2000
deIndex <- sample(K,nDE,replace = FALSE,prob = probs)
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

