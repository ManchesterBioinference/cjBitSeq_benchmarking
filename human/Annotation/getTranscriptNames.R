con <- file("../../Annotation/transcriptome_data/known.fa") 
conOut <- file("transcriptNames.txt",open = "w")
spanki <- read.table("transcript_sims.txt",header = TRUE)
open(con);
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
	if(length(grep(">", line))>0){
		trName <- strsplit(line,split = " ")[[1]][2]
		nReads <- as.numeric(spanki[match(trName,spanki[,2]),][5])
		cat(trName," ",nReads, "\n",file = conOut,append=TRUE)
		print(trName)
	}
} 
close(con)
close(conOut)


