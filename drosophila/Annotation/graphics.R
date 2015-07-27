pplr <- read.table(file = "bitseq/data.pplr")[,1]
cufflinksDE <- read.table("cufflinks/cufflinksDE.txt")[,1]
collapsed <- read.table("cjBitSeqOutput/estimates.txt",header=TRUE)
bitseqPredictions <- 2*abs(0.5 - pplr)
cufflinksPredictions <- cufflinksDE
ebSeqDE2 <- as.numeric(read.table("rsem/ebSeqDe.txt")[,1])
ebseqPredictions2 <- ebSeqDE2

n1 <- length(list.files(pattern = "tr_File_A"))
n2 <- length(list.files(pattern = "tr_File_B"))

theta <- read.table("A1/transcriptNames.txt")[,2]
if(n1 > 1)
for (k in 2:n1){
theta <- cbind(theta,read.table(paste("A",k,"/transcriptNames.txt",sep=""))[,2])
}
w <- read.table("B1/transcriptNames.txt")[,2]
for (k in 2:n2){
w <- cbind(w,read.table(paste("B",k,"/transcriptNames.txt",sep=""))[,2])
}

for(k in 1:n1){
theta[,k] <- theta[,k]/sum(theta[,k])
}
for(k in 1:n2){
w[,k] <- w[,k]/sum(w[,k])
}

K <- dim(theta)[1]
trueDE <- as.numeric(read.table("realDE.txt")[,1])

##################################################################################################

# rocr package
library(ROCR)
library(pracma)
library(fields)

cufflinksPredictions <- cufflinksDE
bitseqPredictions <- 2*abs(0.5 - pplr)
collapsed2 <- collapsed[,3]
ebseqPredictions <- ebSeqDE2
trueClass <- trueDE

pred <- prediction(bitseqPredictions, trueClass)
perf <- performance(pred,"sar")
myMin <- min(mean(unlist(perf@y.values)))
myMax <- max(mean(unlist(perf@y.values)))
perfMean <- mean(unlist(perf@y.values))
perfMeanUP <- mean(unlist(perf@y.values)[which(unlist(perf@x.values)>0.9)])
predCol <- prediction(collapsed2, trueClass)
perfCol <- performance(predCol,"sar")
myMin <- min(myMin,min(mean(unlist(perfCol@y.values))))
myMax <- max(myMax,max(mean(unlist(perfCol@y.values))))
perfColMean <- mean(unlist(perfCol@y.values))
perfColMeanUP <- mean(unlist(perfCol@y.values)[which(unlist(perfCol@x.values)>0.9)])
predC <- prediction(cufflinksPredictions, trueClass)
perfC <- performance(predC,"sar")
myMin <- min(myMin,min(mean(unlist(perfC@y.values))))
myMax <- max(myMax,max(mean(unlist(perfC@y.values))))
perfCMean <- mean(unlist(perfC@y.values))
perfCMeanUP <- mean(unlist(perfC@y.values)[which(unlist(perfC@x.values)>0.9)])
predE <- prediction(ebseqPredictions, trueClass)
perfE <- performance(predE,"sar")
myMin <- min(myMin,min(mean(unlist(perfE@y.values))))
myMax <- max(myMax,max(mean(unlist(perfE@y.values))))
perfEMean <- mean(unlist(perfE@y.values))
perfEMeanUP <- mean(unlist(perfE@y.values)[which(unlist(perfE@x.values)>0.9)])
pdf(file = "SAR-curves.pdf",width = 9,height = 6,pointsize=12)
par(mfrow=c(1,1),mar = c(4,4,2,2));
plot(c(0,1),c(0,1),ylim = c(myMin-0.1,myMax+0.02),type = "n",ylab = "SAR-measure", xlab = "cut-off",main = "SAR", 
     panel.first = {
         usr <- par('usr')
         rect(0.9,0,1,1, col='gray', border=NA)
     } )
legendMeans <- paste(c(round(c(perfMean,perfColMean,perfCMean,perfEMean),3)),c(round(c(perfMeanUP,perfColMeanUP,perfCMeanUP,perfEMeanUP),3)),sep=", ")
legend("bottomright",c("mean: OVERALL, CONSTRAINED",paste(c("BitSeq","cjBitSeq","CuffDif","EbSeq"),legendMeans,sep=": ")),col = c("white","red","darkorange","green","purple"),lty=1.5,bty="n")
plot(perfCol,avg = "vertical",col = "darkorange",add = TRUE,lwd = 1.5)
plot(perfC,avg = "vertical",col = "green",add = TRUE,lwd = 1.5)
plot(perf,avg = "vertical",col = "red",add = TRUE,lwd = 1.5)
plot(perfE,avg = "vertical",col = "purple",add = TRUE,lwd = 1.5)
dev.off()


trueClass <- trueDE
#trueClass <- newDE
bitseqPredictions <- 2*abs(0.5 - pplr)
cufflinksPredictions <- cufflinksDE
pred <- prediction(bitseqPredictions, trueClass)
perf <- performance(pred,"tpr","fpr")
pred <- prediction(cufflinksPredictions, trueClass)
perf2 <- performance(pred,"tpr","fpr")
pred <- prediction(ebseqPredictions2, trueClass)
perf00 <- performance(pred,"tpr","fpr")
strike <- perf@x.values[[1]]
volatility <- perf@y.values[[1]]
bitseqAUC = round(trapz(strike,volatility),3)
strike <- perf2@x.values[[1]]
volatility <- perf2@y.values[[1]]
cuffdifAUC = round(trapz(strike,volatility),3)
strike <- perf00@x.values[[1]]
volatility <- perf00@y.values[[1]]
ebseqAUC2 = round(trapz(strike,volatility),3)

predCol <- prediction(collapsed[,3], trueClass)
perfCol <- performance(predCol,"tpr","fpr")
strike <- perfCol@x.values[[1]]
volatility <- perfCol@y.values[[1]]
ColAUC = round(trapz(strike,volatility),3)

pdf(file = "ROC-curves.pdf",width = 9,height = 6,pointsize=12)
par(mfrow=c(1,1),mar = c(4,4,2,2));
plot(perf@x.values[[1]],perf@y.values[[1]],type= "l",xlab = "False Positive Rate",ylab = "True Positive Rate",col = 2,xlim = c(0,0.4),lwd =1.5)
points(perf2@x.values[[1]],perf2@y.values[[1]],type= "l",col = 3,lwd =1.5)
points(perf00@x.values[[1]],perf00@y.values[[1]],type= "l",col = 6,lwd =1.5)
points(perfCol@x.values[[1]],perfCol@y.values[[1]],type= "l",col = "darkorange",lwd =1.5)
legend("bottomright",c("area under curve",paste("bitseq: ",bitseqAUC),paste("cuffdif: ",cuffdifAUC),paste("collapsed: ",ColAUC),paste("ebseq: ",ebseqAUC2)), col = c("white",2,3,"darkorange",6),lty = 1,bty="n",lwd =1.5)
dev.off()


#

pdf(file = "prec-rec-curves.pdf",width = 9,height = 6,pointsize = 12)
par(mfrow=c(1,1),mar = c(4,4,2,2));
Collapsed.perf <- performance(predCol, "prec", "rec");
plot(Collapsed.perf@x.values[[1]],Collapsed.perf@y.values[[1]],ylim = c(0,1),xlim = c(0,1),type = "l",col = "white",ylab = "precision",xlab = "recall",lwd =1.5,lty = 2)
points(Collapsed.perf@x.values[[1]],Collapsed.perf@y.values[[1]],type = "l",col = "darkorange",lwd =1.5)
pred <- prediction(2*abs(0.5-pplr), trueClass)
RP.bs <- performance(pred, "prec", "rec");
points(RP.bs@x.values[[1]],RP.bs@y.values[[1]],type = "l",col = "red",lwd =1.5)
pred <- prediction(ebseqPredictions2, trueClass)
RP.bs <- performance(pred, "prec", "rec");
points(RP.bs@x.values[[1]],RP.bs@y.values[[1]],type = "l",col = "purple",lwd =1.5)
pred <- prediction(cufflinksPredictions, trueClass)
RP.bs <- performance(pred, "prec", "rec");
points(RP.bs@x.values[[1]],RP.bs@y.values[[1]],type = "l",col = "green",lwd =1.5)
legend("bottomleft",c("bitseq","cuffdif","collapsed","ebseq"),col = c("red","green","darkorange","purple"),lty = 1,bty="n",lwd =1.5)
dev.off()




pdf(file = "trueRanking2.pdf",width = 9,height = 6,pointsize = 12)
par(mfrow=c(2,2),mar = c(4,4,2,2));
d <- cbind(log(rowMeans(theta)),log(rowMeans(w))) 
myCEX = 0.5
cs <- 2*abs(0.5-pplr)
perm <- order(cs) #,decreasing = TRUE)
plot(d[perm,],col = color.scale(cs[perm]),main = "bitseq",xlab = "log(theta)",ylab = "log(w)",pch = 16,cex=myCEX);
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
points(d[perm,],col = color.scale(cs[perm]),main = "bitseq",xlab = " ",ylab = " ",pch = 16,cex=myCEX);
colorbar.plot( x=-11.5, y=-7, strip =cs[perm] ,col = color.scale(seq(0,1,length=100)))
cs <- cufflinksPredictions
perm <- order(cs)#,decreasing = TRUE)
plot(d[perm,],col = color.scale(cs[perm]),main = "cuffdiff",xlab = "log(theta)",ylab = "log(w)",pch = 16,cex=myCEX);
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
points(d[perm,],col = color.scale(cs[perm]),main = "cuffdiff",xlab = " ",ylab = " ",pch = 16,cex=myCEX);
colorbar.plot( x=-11.5, y=-7, strip =cs[perm] ,col = color.scale(seq(0,1,length=100)))
cs <- collapsed[,3]
perm <- order(cs)#,decreasing = TRUE)
plot(d[perm,],col = color.scale(cs[perm]),main = "cjBitSeq",xlab = "log(theta)",ylab = "log(w)",pch = 16,cex=myCEX);
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
points(d[perm,],col = color.scale(cs[perm]),main = "rjmcmc",xlab = " ",ylab = " ",pch = 16,cex=myCEX);
colorbar.plot( x=-11.5, y=-7, strip =cs[perm] ,col = color.scale(seq(0,1,length=100)))
cs <- ebseqPredictions2
perm <- order(cs)#,decreasing = TRUE)
plot(d[perm,],col = color.scale(cs[perm]),main = "EBSeq",xlab = "log(theta)",ylab = "log(w)",pch = 16,cex=myCEX);
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
points(d[perm,],col = color.scale(cs[perm]),main = "EBSeq",xlab = " ",ylab = " ",pch = 16,cex=myCEX);
colorbar.plot( x=-11.5, y=-7, strip =cs[perm] ,col = color.scale(seq(0,1,length=100)))
dev.off()










myCut <- c(0.01,0.025,0.05,0.1,0.2,0.4)
l <- length(myCut)

realDE <- array(trueDE,dim = c(length(trueDE),1)) #read.table("realDE.txt")
p <- read.table("cjBitSeqOutput/estimates.txt",header=TRUE)
perm <- order(p[,3],decreasing = TRUE)
orderedP <- p[perm,3]
nDE <- length(which(realDE==1))


aoua <- array(data = NA, dim =c(l,3))
iter <- 0
for (alpha in myCut){
	iter <- iter + 1 

	K <- dim(p)[1]
	myList <-  1 - orderedP[1]
	k <- 1 
	criterion <- myList
	while (criterion < alpha){
		k <- k + 1 
		myList <- myList + 1 - orderedP[k]
		criterion <- myList/k
	}

	ind <- perm[1:(k-1)]
	if (dim(table(realDE[ind,1])) > 1){
		point1 <- as.numeric(table(realDE[ind,1])[1]/length(ind))  #achieved fdr
		point2 <- as.numeric(table(realDE[ind,1])[2]/nDE)  #achieved tpr
	}else{
		point1 <- 0
		point2 <- as.numeric(length(ind)/nDE)  #achieved tpr
	}
	aoua[iter,] <- c(point1,point2,alpha)

}
pdf(file = "fdr.pdf",width = 9,height = 6,pointsize=12)
par(mfrow=c(1,1),mar = c(4,4,2,2));
plot(aoua[,c(1,2)],xlab = "achieved FDR",ylab = "power",xlim = c(0,0.6),ylim = c(0,1),type = "l",col = "darkorange",lwd = 2)
points(aoua[,c(1,2)],pch = 1:l,col =  "darkorange",lwd = 2)


#ebseq

p <- read.table("rsem/ebSeqDe.txt")
perm <- order(p[,1],decreasing = TRUE)
#plot(p[perm,1])
orderedP <- p[perm,1]

aoua <- array(data = NA, dim =c(l,3))
iter <- 0
for (alpha in myCut){
	iter <- iter + 1 

	K <- dim(p)[1]
	myList <-  1 - orderedP[1]
	k <- 1 
	criterion <- myList
	while (criterion < alpha){
		k <- k + 1 
		myList <- myList + 1 - orderedP[k]
		criterion <- myList/k
	#	cat(paste("k = ",k,". Criterion = ",criterion,sep=""),"\n")
	}

	ind <- perm[1:(k-1)]
	ind <- which(p[,1] > 1-alpha)
	if (dim(table(realDE[ind,1])) > 1){
		point1 <- as.numeric(table(realDE[ind,1])[1]/length(ind))  #achieved fdr
		point2 <- as.numeric(table(realDE[ind,1])[2]/nDE)  #achieved tpr
	}else{
		point1 <- 0
		point2 <- as.numeric(length(ind)/nDE)  #achieved tpr
	}
	aoua[iter,] <- c(point1,point2,alpha)

}
points(aoua[,1:2],col = "purple",type="l",lwd = 2)
points(aoua[,1:2],col = "purple",type="p",pch = 1:l,lwd = 2)


#cufflinks

p <- read.table("cufflinks/cufflinksDE.txt")
perm <- order(p[,1],decreasing = TRUE)
#plot(p[perm,1])
orderedP <- p[perm,1]

aoua <- array(data = NA, dim =c(l,3))
iter <- 0
for (alpha in myCut){
	iter <- iter + 1 

	K <- dim(p)[1]
	myList <-  1 - orderedP[1]
	k <- 1 
	criterion <- myList
	while (criterion < alpha){
		k <- k + 1 
		myList <- myList + 1 - orderedP[k]
		criterion <- myList/k
	#	cat(paste("k = ",k,". Criterion = ",criterion,sep=""),"\n")
	}

	ind <- perm[1:(k-1)]
	ind <- which(p[,1] > 1-alpha)
	if (dim(table(realDE[ind,1])) > 1){
		point1 <- as.numeric(table(realDE[ind,1])[1]/length(ind))  #achieved fdr
		point2 <- as.numeric(table(realDE[ind,1])[2]/nDE)  #achieved tpr
	}else{
		point1 <- 0
		point2 <- as.numeric(length(ind)/nDE)  #achieved tpr
	}
	aoua[iter,] <- c(point1,point2,alpha)

}
points(aoua[,1:2],col = "green",type="l",lwd = 2)
points(aoua[,1:2],col = "green",type="p",pch = 1:l,lwd = 2)

#bitseq
p <- 2*abs(0.5 - pplr)
perm <- order(p,decreasing = TRUE)
orderedP <- p[perm]
nDE <- length(which(realDE==1))


aoua <- array(data = NA, dim =c(l,3))
iter <- 0
for (alpha in myCut){
	iter <- iter + 1 

	K <- dim(p)[1]
	myList <-  1 - orderedP[1]
	k <- 1 
	criterion <- myList
	while (criterion < alpha){
		k <- k + 1 
		myList <- myList + 1 - orderedP[k]
		criterion <- myList/k
	}

	ind <- perm[1:(k-1)]
	if (dim(table(realDE[ind,1])) > 1){
		point1 <- as.numeric(table(realDE[ind,1])[1]/length(ind))  #achieved fdr
		point2 <- as.numeric(table(realDE[ind,1])[2]/nDE)  #achieved tpr
	}else{
		point1 <- 0
		point2 <- as.numeric(length(ind)/nDE)  #achieved tpr
	}
	aoua[iter,] <- c(point1,point2,alpha)

}
points(aoua[,c(1,2)],type = "l",col = "red",lwd = 2)
points(aoua[,c(1,2)],pch = 1:l,col =  "red",lwd = 2)




#legend("bottomleft",c("rjmcmc (jeffreys prior)","rjmcmc (0.5 prior)","ebseq","cuffdif"),col = c("blue","blue","purple","green"),lty = c(1,2,1,1),lwd=2,bty="n")

#legend("bottomleft",c("rjmcmc","ebseq","cuffdif"),col = c("blue","purple","green"),lty = 1,lwd=2,bty="n")
legend("bottomright",paste("eFDR = ",myCut,sep=""),pch = 1:l,bty="n")

dev.off()





library('VennDiagram')
alpha <- 0.05
ind_col <- which(collapsed[,5] == "DE")
ind_rj <- which(trueClass == 1) #true
ind_b <- c(which(pplr < alpha/2),which(pplr > 1 - alpha/2))
ind_eb <- which(ebSeqDE2 > 1-alpha)
ind_cuf <- which(cufflinksDE > 1-alpha)

#1 : rj
#2 : bitseq
#3 : ebseq
#4 : cuffdif
#5 : collapsed

n12 <- intersect(ind_rj,ind_b)
n13 <- intersect(ind_rj,ind_eb)
n14 <- intersect(ind_rj,ind_cuf)
n15 <- intersect(ind_rj,ind_col)
n23 <- intersect(ind_b,ind_eb)
n24 <- intersect(ind_b,ind_cuf)
n25 <- intersect(ind_b,ind_col)
n34 <- intersect(ind_eb,ind_cuf)
n35 <- intersect(ind_eb,ind_col)
n45 <- intersect(ind_cuf,ind_col)
n123 <- intersect(n12,ind_eb)
n124 <- intersect(n12,ind_cuf)
n125 <- intersect(n12,ind_col)
n134 <- intersect(n13,ind_cuf)
n135 <- intersect(n13,ind_col)
n145 <- intersect(n14,ind_col)
n234 <- intersect(n23,ind_cuf)
n235 <- intersect(n23,ind_col)
n245 <- intersect(n24,ind_col)
n345 <- intersect(n34,ind_col)
n1234 <- intersect(n123,ind_cuf)
n1235 <- intersect(n123,ind_col)
n1245 <- intersect(n124,ind_col)
n1345 <- intersect(n134,ind_col)
n2345 <- intersect(n234,ind_col)
n12345 <- intersect(n1234,ind_col)
catNames <- c("ground-truth","BitSeq","EbSeq","Cuffdiff","cjBitSeq")

pdf(file = "vennDiagram.pdf",width = 6,height=6,pointsize=9)
#plot(c(0,1),c(0,1),type = "n",main = paste("HiSeq"),xaxt="n",yaxt="n",xlab = "",ylab = "",bty="n",cex.main=2)
venn.plot <- draw.quintuple.venn(
area1 = length(ind_rj),
area2 = length(ind_b),
area3 = length(ind_eb),
area4 = length(ind_cuf),
area5 = length(ind_col),
n12 = length(n12),
n13 = length(n13),
n14 = length(n14),
n15 = length(n15),
n23 = length(n23),
n24 = length(n24),
n25 = length(n25),
n34 = length(n34),
n35 = length(n35),
n45 = length(n45),
n123 = length(n123),
n124 = length(n124),
n125 = length(n125),
n134 = length(n134),
n135 = length(n135),
n145 = length(n145),
n234 = length(n234),
n235 = length(n235),
n245 = length(n245),
n345 = length(n345),
n1234 = length(n1234),
n1235 = length(n1235),
n1245 = length(n1245),
n1345 = length(n1345),
n2345 = length(n2345),
n12345 = length(n12345),
category = catNames,
fill = c("dodgerblue", "red", "purple", "green", "darkorange"),
cat.col = c("dodgerblue", "red", "purple", "green", "darkorange"),
cat.cex = 2,
margin = 0.1,
cex = 1.5*c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
ind = TRUE
);
dev.off()







