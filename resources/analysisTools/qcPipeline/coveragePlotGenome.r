#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

# example call:      R -f /home/jaegern/rWorkspace/coveragePlotGenome.r --no-save --no-restore --args "/icgc/lsdf/mb/analysis/QC_pipeline/results_per_pid/179743/coverage/control_179743_merged.bam.rmdup_readCoverage_1kb_windows_sorted.txt" "/icgc/lsdf/mb/analysis/QC_pipeline/results_per_pid/179743/coverage/tumor_179743_merged.bam.rmdup_readCoverage_1kb_windows_sorted.txt" "read" "1" "/icgc/lsdf/mb/analysis/QC_pipeline/results_per_pid/179743/coverage/readDepth.png"
#                    R -f /home/jaegern/rWorkspace/coveragePlotGenome.r --no-save --no-restore --args "/icgc/lsdf/mb/analysis/QC_pipeline/results_per_pid/198151/coverage/control_198151_merged.bam.rmdup_readCoverage_1kb_windows_sorted.txt" "/icgc/lsdf/mb/analysis/QC_pipeline/results_per_pid/198151/coverage/tumor_198151_merged.bam.rmdup_readCoverage_1kb_windows_sorted.txt" "read" "1" "/icgc/lsdf/mb/analysis/QC_pipeline/results_per_pid/198151/coverage/readDepth.png"
        


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

# for testing
#setwd("/ibios/temp2/coveragePlots")
#sampleOne = "tumor_AvD_59208_readCoverage_1kb_windows.txt"
#sampleTwo = "control_AvD_59208_readCoverage_1kb_windows.txt"
#sampleOne = "tumor_mouse.fa_VU_361_readDepth_1kb_windows.txt"
#sampleTwo = "control_mouse.fa_VU_361_readDepth_1kb_windows.txt"
#countType = "read" 
#windowSize = "1"
#outputFile = "readDepth.png"
#twoSamples = TRUE



#input coverage data
sampleOneRawUnordered=read.table(sampleOne,sep="\t",as.is=TRUE)
if (twoSamples) sampleTwoRawUnordered=read.table(sampleTwo,sep="\t",as.is=TRUE)

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

# please add check if "chr" added in chrom identifier, then strip off "chr" to have one script working for all human genome versions, also strip off ".fa" as suffix, just in case it is EMBL data

# mouse genome would also possible to include by just making the chrom number list shorter, could also be added as option and vector set above
# reorder to standard chromosome order  # the following 3 lines are the new part to deal with hg19_1000genomes version of the genome

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



#plot read depth for complete genome into one plot
png(outputFile ,width=2100,height=1000)
par(mfrow=c(3,1),mar=c(2,4,2,2))

#sampleOne (= control sample)
plot(log(sampleOneRaw[,3]+1,2),pch=".",main=paste(sampleOne,sep="  "),cex.main = 1.5,xaxs="i",ylab=paste("log2 #",countType," per ",windowSize,"kb",sep=""),cex.lab=1.5,axes=FALSE)
box()
axis(2, cex.axis = 1.5)
abline(h=log(mean(sampleOneRaw[,3])+1,2),col=2)  # red line at the mean coverage configurationLevel
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
axis(1,at=pretty(c(0,chrLengthSum)),labels=paste(pretty(c(0,chrLengthSum/1000)),"Mb",sep=""), cex.axis = 1.5)



#sampleTwo (= tumor sample)
if (twoSamples) {
plot(log(sampleTwoRaw[,3]+1,2),pch=".",main=paste(sampleTwo,sep="  "),cex.main = 1.5,xaxs="i",ylab=paste("log2 #",countType," per ",windowSize,"kb",sep=""),cex.lab=1.5,axes=FALSE)
box()
axis(2, cex.axis = 1.5)
abline(h=log(mean(sampleTwoRaw[,3])+1,2),col=2)  # red line at the mean coverage configurationLevel
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
axis(1,at=pretty(c(0,chrLengthSum)),labels=paste(pretty(c(0,chrLengthSum/1000)),"Mb",sep=""), cex.axis = 1.5)



# ONLY THIS PART NEC FOR PLOTTING WITH COMBINED BAF PLOT
#Ratio between the 2 samples normalized to #reads all chromosomes
plot(log(((sampleTwoRaw[,3])/sumSampleTwo)/((sampleOneRaw[,3])/sumSampleOne),2),pch=".",main=paste("log2 ratio ","tumor"," vs. ","control genomewide",sep=""),cex.main = 1.5,xaxs="i",ylab=paste("#",countType," per 1kb normalized to total reads",sep=""),cex.lab=1.5,axes=FALSE)
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


axis(1,at=pretty(c(0,chrLengthSum)),labels=paste(pretty(c(0,chrLengthSum/1000)),"Mb",sep=""), cex.axis = 1.5)
abline(h=0,col=2)

}

dev.off()

