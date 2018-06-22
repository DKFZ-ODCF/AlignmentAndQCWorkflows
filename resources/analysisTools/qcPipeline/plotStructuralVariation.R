#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

# example script call: R -f structVarPlot.R --no-save --no-restore --args "/ngs/s_1_paired_mappedToDiffChr.txt" "/plots/" "52167542"

cmdArgs = commandArgs(TRUE)
print(cmdArgs)
if (length(cmdArgs) != 3) print(paste("Incorrect number of arguments (3 expected): ",length(cmdArgs)))
sample = cmdArgs[1]
outputPath = cmdArgs[2]
totalReads = cmdArgs[3]

# import libraries
library(gtools)

# read data
sampleRawTable=read.table(sample,header=F,sep="\t", comment.char = "%")

chrFreq <- table(as.vector(sampleRawTable$V2))

matrixNames=matrix(names(chrFreq))
matrixValues=matrix(chrFreq)
matrixChrFreq=matrix(c(matrixNames,matrixValues),length(matrixNames))
matrixChrFreq=matrixChrFreq[mixedorder(matrixChrFreq[,1]),]

values = as.integer(as.vector(matrixChrFreq[,2]))
names(values) = as.vector(matrixChrFreq[,1])


# prepare names for output files and plot labels
sampleNameTemp = unlist(strsplit(sample,paste('.txt',sep="")))[1] # also contains the file path
sampleName = tail(unlist(strsplit(sampleNameTemp,"/")), 1) # filename without path

# plot barplot into a .png file
filename=paste(outputPath,sampleName,".png",sep="")
png(filename, width = 680, height = 680)

captionPercentage= round(sum(chrFreq) / as.integer(totalReads) *100,2)

if (captionPercentage < 1){
	par(las=2)
	barplot(values,col=rainbow(length(values)),main= paste(sampleName,"\n(",captionPercentage,"% of total reads mapped have mate mapped to a different chr)",sep=" "),font.main=1,xlab="chromosome",ylab="#reads",ylim=c(0,max(matrixValues)))
}else{
	par(las=2)
	barplot(values,col=rainbow(length(values)),main=paste(sampleName,"\n(",captionPercentage,"% of total reads mapped have mate mapped to a different chr)",sep=" "),font.main=2,col.main="red",xlab="chromosome",ylab="#reads",ylim=c(0,max(matrixValues)))
}

dev.off()

# write the statistic values to a simple text file
sink(paste(filename,"_qcValues.txt",sep=""), append=FALSE, split=FALSE)
cat(paste(captionPercentage, "\n"))
sink()