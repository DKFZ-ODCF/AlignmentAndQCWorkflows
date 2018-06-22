#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

# example call:      R -f /home/hutter/workspace_ngs/Roddy/analysisTools/qcPipeline/coveragePlot.R --no-save --no-restore --args /icgc/dkfzlsdf/analysis/hipo/hipo_021/results_per_pid/P021-PTV8/coverage/buffy_coat_P021-PTV8_readCoverage_10kb_windows.txt /icgc/dkfzlsdf/analysis/hipo/hipo_021/results_per_pid/P021-PTV8/coverage/tumor_P021-PTV8_readCoverage_10kb_windows.txt read 10 /icgc/dkfzlsdf/analysis/hipo/hipo_021/results_per_pid/P021-PTV8/coverage/readDepth.png



# OBTAIN CONFIGURATION PARAMETERS FOR PLOTTING THE READ/BASE COVERAGE (WITH R)

cmdArgs = commandArgs(TRUE)
print(cmdArgs)
if (length(cmdArgs) < 4) print(paste("Incorrect number of arguments (at least 4 expected): ",length(cmdArgs)))
sampleOne = cmdArgs[1]
twoSamples = FALSE
if (length(cmdArgs) == 4) {
  countType = cmdArgs[2]  # specifies whether reads or bases have been counted for coverage within the windowSize
  windowSize = cmdArgs[3]
  outputFile = cmdArgs[4]
} else {
  sampleTwo = cmdArgs[2]
  countType = cmdArgs[3]  # specifies whether reads or bases have been counted for coverage within the windowSize
  windowSize = cmdArgs[4]
  outputFile = cmdArgs[5] 
  twoSamples = TRUE
}

# convert windowSize given in kb to bp
bp_scale=(1000/as.numeric(windowSize))


#input coverage data
sampleOneRawUnordered=read.table(sampleOne,sep="\t",as.is=TRUE)
if (twoSamples) sampleTwoRawUnordered=read.table(sampleTwo,sep="\t",as.is=TRUE)

# strip off "chr" prefix and ".fa" suffix if present
a = as.character(sampleOneRawUnordered[[1]])
a = gsub("^chr", "", a, ignore.case = TRUE)
a = gsub("\\.fa$", "", a, ignore.case = TRUE)
a = gsub("^(\\d)$", "0\\1", a)
sampleOneRawUnordered[[1]] = a

if (twoSamples) {
    a = as.character(sampleTwoRawUnordered[[1]])
    a = gsub("^chr", "", a, ignore.case = TRUE)
    a = gsub("\\.fa$", "", a, ignore.case = TRUE)
    a = gsub("^(\\d)$", "0\\1", a)
    sampleTwoRawUnordered[[1]] = a
}

# get numerically sorted chromosomes
chr.index = sort(unique(sampleOneRawUnordered[[1]]))
sampleOneRaw <- c()
for( i in chr.index){
  sampleOneRaw <- rbind(sampleOneRaw, sampleOneRawUnordered[which(sampleOneRawUnordered$V1 == i), ])
  
}

if (twoSamples) {
  sampleTwoRaw <- c()
  for( i in sort(unique(sampleOneRawUnordered[[1]]))){
    sampleTwoRaw <- rbind(sampleTwoRaw, sampleTwoRawUnordered[which(sampleTwoRawUnordered$V1 == i), ])
    
  }
  
}


#conversion function for coverage data
con=function(x){
	x=data.frame(ref=x[,1], start=x[,2], reads=x[,3])
	x=split(x,x$ref)
}

sampleOneCon=con(sampleOneRaw)
if (twoSamples) sampleTwoCon=con(sampleTwoRaw)

#calculate total number of reads for normalization
sumSampleOne=sum(as.numeric(sapply(sampleOneCon,function(x)sum(x$reads))))
if (twoSamples) sumSampleTwo=sum(as.numeric(sapply(sampleTwoCon,function(x)sum(x$reads))))


# individual chromosomes
#get names for plot labels
sampleOneNameTemp = unlist(strsplit(sampleOne,paste('_',countType,'Coverage_', windowSize,'kb_windows.txt',sep="")))[1] # also contains the file path
sampleOneName = tail(unlist(strsplit(sampleOneNameTemp,"/")), 1) # filename without path
if (twoSamples) {
  sampleTwoNameTemp = unlist(strsplit(sampleTwo,paste('_',countType,'Coverage_', windowSize,'kb_windows.txt',sep="")))[1]
  sampleTwoName = tail(unlist(strsplit(sampleTwoNameTemp,"/")), 1)
  analysisName = paste(sampleOneName,"_vs_",sampleTwoName,sep="")
} else {
  analysisName = sampleOneName
}

#plot coverage per chromosome 	
for( i in names(sampleOneCon)){
    
    chr.name = i
    chr.name = gsub("^0*", "", i)
    chr.name = paste("chr", chr.name, sep = "")
    
    fname = outputFile
    fname = gsub("\\.png", paste("_", chr.name, ".png", sep = ""), fname)
    
	png(fname, width=1100,height=1000)
	par(mfrow=c(3,1),cex=1.5,mar=c(2,4,2,0.1))#
	#sampleOne#
	plot(log(sampleOneCon[[i]]$reads+1,2),pch=".",main=paste(sampleOneName,"_", chr.name,sep=""),ylab=paste("log2 #",countType," per ",windowSize,"kb",sep=""),axes=FALSE)
	box()#
	axis(2)#
	axis(1,at=pretty(c(0,nrow(sampleOneCon[[i]]))),labels=paste(pretty(c(0,nrow(sampleOneCon[[i]]))/bp_scale),"Mb",sep=""))#
	#sampleTwo#
	if (twoSamples) {
	  plot(log(sampleTwoCon[[i]]$reads+1,2),pch=".",main=paste(sampleTwoName,"_", chr.name,sep=""),ylab=paste("log2 #",countType," per ",windowSize,"kb",sep=""),axes=FALSE)
	  box()#
	  axis(2)#for( i in names(sampleOneCon)){
	  axis(1,at=pretty(c(0,nrow(sampleTwoCon[[i]]))),labels=paste(pretty(c(0,nrow(sampleTwoCon[[i]]))/bp_scale),"Mb",sep=""))#
	  #
	  #Ratio between the 2 samples normalized to #reads all chromosomes
	  plot(log(((sampleTwoCon[[i]]$reads)/sumSampleTwo)/((sampleOneCon[[i]]$reads)/sumSampleOne),2),pch=".",main=paste("log2 ratio ",sampleTwoName," vs. ", sampleOneName, " ", chr.name,sep=""),ylab=paste("#",countType," per ", windowSize, "kb normalized to total reads",sep=""),axes=FALSE)#
	  box()#
	  axis(2)#
	  axis(1,at=pretty(c(0,nrow(sampleOneCon[[i]]))),labels=paste(pretty(c(0,nrow(sampleOneCon[[i]]))/bp_scale),"Mb",sep=""))#
	  abline(h=0,col=2)#
	  }
	dev.off()#
	}

#plot read depth for complete genome into one plot
png(outputFile ,width=2100,height=1000)
par(mfrow=c(3,1),mar=c(2,4,2,2))

#sampleOne (= control sample)
plot(log(sampleOneRaw[,3]+1,2),pch=".",main=paste(sampleOneName,sep="  "),cex.main = 1.5,xaxs="i",ylab=paste("log2 #",countType," per ",windowSize,"kb",sep=""),cex.lab=1.5,axes=FALSE)
box()
axis(2, cex.axis = 1.5)
abline(h=log(mean(sampleOneRaw[,3])+1,2),col=2)  # red line at the mean coverage level
abline(v=0,col="grey3")

chrLengthSum=0  # to separate each chromosome by a grey line


#for( i in names(sampleOneCon)){
for( i in chr.index){    
    chr.name = i
    chr.name = gsub("^0*", "", i)
    chr.name = paste("chr", chr.name, sep = "")
	text(x = chrLengthSum+round((nrow(sampleOneCon[[i]])/2)), y = 13, labels = chr.name, pos = 3, cex = 1.4)

	chrLengthSum=chrLengthSum+nrow(sampleOneCon[[i]])
	abline(v=chrLengthSum,col="grey3")
}
axis(1,at=pretty(c(0,chrLengthSum)),labels=paste(pretty(c(0,chrLengthSum/bp_scale)),"Mb",sep=""), cex.axis = 1.5)



#sampleTwo (= tumor sample)
if (twoSamples) {
plot(log(sampleTwoRaw[,3]+1,2),pch=".",main=paste(sampleTwoName,sep="  "),cex.main = 1.5,xaxs="i",ylab=paste("log2 #",countType," per ",windowSize,"kb",sep=""),cex.lab=1.5,axes=FALSE)
box()
axis(2, cex.axis = 1.5)
abline(h=log(mean(sampleTwoRaw[,3])+1,2),col=2)  # red line at the mean coverage level
abline(v=0,col="grey3")

chrLengthSum=0  # to separate each chromosome by a grey line

#for( i in names(sampleTwoCon)){
for( i in chr.index){
    
    chr.name = i
    chr.name = gsub("^0*", "", i)
    chr.name = paste("chr", chr.name, sep = "")
	text(x = chrLengthSum+round((nrow(sampleTwoCon[[i]])/2)), y = 13, labels = chr.name, pos = 3, cex = 1.4)

	chrLengthSum=chrLengthSum+nrow(sampleTwoCon[[i]])
	abline(v=chrLengthSum,col="grey3")
}
axis(1,at=pretty(c(0,chrLengthSum)),labels=paste(pretty(c(0,chrLengthSum/bp_scale)),"Mb",sep=""), cex.axis = 1.5)


#Ratio between the 2 samples normalized to #reads all chromosomes
plot(log(((sampleTwoRaw[,3])/sumSampleTwo)/((sampleOneRaw[,3])/sumSampleOne),2),pch=".",main=paste("log2 ratio ",sampleTwoName," vs. ",sampleOneName," genomewide",sep=""),cex.main = 1.5,xaxs="i",ylab=paste("#",countType," per ", windowSize, "kb normalized to total reads",sep=""),cex.lab=1.5,axes=FALSE)
box()
axis(2, cex.axis = 1.5)

abline(v=0,col="grey3")

chrLengthSum=0  # to separate each chromosome by a grey line

#for( i in names(sampleTwoCon)){
for( i in chr.index){
    chr.name = i
    chr.name = gsub("^0*", "", i)
    chr.name = paste("chr", chr.name, sep = "")
	text(x = chrLengthSum+round((nrow(sampleTwoCon[[i]])/2)), y = 2, labels = chr.name, pos = 3, cex = 1.4)

	chrLengthSum=chrLengthSum+nrow(sampleTwoCon[[i]])
	abline(v=chrLengthSum,col="grey3")
}


axis(1,at=pretty(c(0,chrLengthSum)),labels=paste(pretty(c(0,chrLengthSum/bp_scale)),"Mb",sep=""), cex.axis = 1.5)
abline(h=0,col=2)

}

dev.off()
