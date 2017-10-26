# OBTAIN CONFIGURATION PARAMETERS FOR PLOTTING THE READ/BASE COVERAGE (WITH R)
cmdArgs = commandArgs(TRUE)
print(cmdArgs)
if (length(cmdArgs) < 4) print(paste("Incorrect number of arguments (at least 4 expected): ",length(cmdArgs)))
sampleOne = cmdArgs[1]
twoSamples = FALSE
if (length(cmdArgs) == 4) {
  countType = cmdArgs[2]  # specifies whether reads or bases have been counted for coverage within the windowSize
  windowSize = cmdArgs[3]
  outputPath = cmdArgs[4]
} else {
  sampleTwo = cmdArgs[2]
  countType = cmdArgs[3]  # specifies whether reads or bases have been counted for coverage within the windowSize
  windowSize = cmdArgs[4]
  outputPath = cmdArgs[5] 
  twoSamples = TRUE
}

#input coverage data
sampleOneRaw=read.table(sampleOne,sep="\t",as.is=TRUE)
if (twoSamples) sampleTwoRaw=read.table(sampleTwo,sep="\t",as.is=TRUE)

#conversion function for coverage data
con=function(x){
	x=data.frame(ref=x[,1], start=x[,2], reads=x[,3])
	x=split(x,x$ref)
}

sampleOneCon=con(sampleOneRaw)
if (twoSamples) sampleTwoCon=con(sampleTwoRaw)

#calculate total number of reads for normalization
sumSampleOne=sum(sapply(sampleOneCon,function(x)sum(x$reads)))
if (twoSamples) sumSampleTwo=sum(sapply(sampleTwoCon,function(x)sum(x$reads)))

#get/prepare names for output files and plot labels
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
	png(paste(outputPath,analysisName,"_coverage_",i,".png",sep=""),width=1100,height=1000)
	par(mfrow=c(3,1),cex=1.5,mar=c(2,4,2,0.1))#
	#sampleOne#
	plot(log(sampleOneCon[[i]]$reads+1,2),pch=".",main=paste(sampleOneName,i,sep="  "),ylab=paste("log2 #",countType," per ",windowSize,"kb",sep=""),axes=FALSE)
	box()#
	axis(2)#
	axis(1,at=pretty(c(0,nrow(sampleOneCon[[i]]))),labels=paste(pretty(c(0,nrow(sampleOneCon[[i]]))/1000),"Mb",sep=""))#
	#sampleTwo#
	if (twoSamples) {
	  plot(log(sampleTwoCon[[i]]$reads+1,2),pch=".",main=paste(sampleTwoName,i,sep="  "),ylab=paste("log2 #",countType," per ",windowSize,"kb",sep=""),axes=FALSE)
	  box()#
	  axis(2)#for( i in names(sampleOneCon)){
	  axis(1,at=pretty(c(0,nrow(sampleTwoCon[[i]]))),labels=paste(pretty(c(0,nrow(sampleTwoCon[[i]]))/1000),"Mb",sep=""))#
	  #
	  #Ratio between the 2 samples normalized to #reads all chromosomes
	  plot(log(((sampleTwoCon[[i]]$reads)/sumSampleTwo)/((sampleOneCon[[i]]$reads)/sumSampleOne),2),pch=".",main=paste("log2 ratio ","tumor"," vs. ","control  ",i,sep=""),ylab=paste("#",countType," per 1kb normalized to total reads",sep=""),axes=FALSE)#
	  box()#
	  axis(2)#
	  axis(1,at=pretty(c(0,nrow(sampleOneCon[[i]]))),labels=paste(pretty(c(0,nrow(sampleOneCon[[i]]))/1000),"Mb",sep=""))#
	  abline(h=0,col=2)#
	  dev.off()#
	  }
	}


