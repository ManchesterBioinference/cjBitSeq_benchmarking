myNames <- read.table("../../Annotation/trNames.tr")[,1]
cuffTranscripts <- read.table("cuff_compare.combined.gtf")[,19]
K <- length(myNames)

perm <- match(myNames,cuffTranscripts)


# myNames[i] = cuffTranscripts[perm[i]] if transcript i is contained in cufftranscripts.
length(which(is.na(perm)==TRUE))

indexNA <- which(is.na(perm)==TRUE) #for this we have no test.

cuffTCONS <- read.table("cuff_compare.combined.gtf")[,c(13,19)]

cufflinksDE <- numeric(K)
indexNA <- which(is.na(perm)==TRUE) #for this we have no test.
cufflinksDE[indexNA] <- rep(0,length(indexNA)) #call them as non-de (?)


# Note that: K = length(unique((cuffTCONS[perm]))) + length(which(is.na(perm)==TRUE)) - 1
cuffResults <- read.table("cuffdif_out/isoform_exp.diff",header = TRUE)[,c(1,13)]

for(i in 1:K){
	if(is.na(perm[i]) == TRUE){
		cufflinksDE[i] <- 0
	}else{

	#	print(cuffTCONS[match(cuffTranscripts[perm[i]],cuffTCONS[,2]),])
	#	cat("\n")
	#	print(myNames[i])
	#	cat("\n")
	#	cuffResults[match(cuffTCONS[match(cuffTranscripts[perm[i]],cuffTCONS[,2]),][,1],cuffResults[,1]),]
		cufflinksDE[i] <- 1 - cuffResults[match(cuffTCONS[match(cuffTranscripts[perm[i]],cuffTCONS[,2]),][,1],cuffResults[,1]),][,2]
	}


}




if(is.na(perm[i]) == TRUE){
	cufflinksDE[i] <- 0
}else{

#	print(cuffTCONS[match(cuffTranscripts[perm[i]],cuffTCONS[,2]),])
#	cat("\n")
#	print(myNames[i])
#	cat("\n")
#	cuffResults[match(cuffTCONS[match(cuffTranscripts[perm[i]],cuffTCONS[,2]),][,1],cuffResults[,1]),]
	cufflinksDE[i] <- 1 - cuffResults[match(cuffTCONS[match(cuffTranscripts[perm[i]],cuffTCONS[,2]),][,1],cuffResults[,1]),][,2]
}
write.table(cufflinksDE,file = "cufflinksDE.txt")











