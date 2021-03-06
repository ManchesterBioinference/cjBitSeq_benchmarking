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
source(file = "../Annotation/produceJobs.R")


