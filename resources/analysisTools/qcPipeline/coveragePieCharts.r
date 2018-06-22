#!/usr/bin/env Rscript
#
# This file is part of the AlignmentAndQCWorkflow plugin.
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this script.  If not, see <https://www.gnu.org/licenses/>.
#
#
# example script call: R -f coveragePieCharts.r --no-save --no-restore --args "/ngs/s_1_bwa.sorted_readCoverage_1kb_windows.txt" "read" "1" "/ngs/plots/"
#
# OBTAIN CONFIGURATION PARAMETERS FOR PLOTTING THE READ/BASE COVERAGE AS PIE CHARTS (WITH R)

cmdArgs = commandArgs(TRUE)
print(cmdArgs)
if (length(cmdArgs) != 4) print(paste("Incorrect number of arguments (4 expected): ",length(cmdArgs)))
sample = cmdArgs[1]   
countType = cmdArgs[2]  # specifies whether reads or bases have been counted for coverage within the given windowSize
windowSize = cmdArgs[3]
outputPath = cmdArgs[4]

# import libraries
library(RColorBrewer)
library(grDevices)


# read coverage data
sampleRawTable=read.table(sample,header=F,sep="\t",as.is=TRUE)

# prepare names for output files and plot labels
sampleNameTemp = unlist(strsplit(sample,paste('_',countType,'Coverage_', windowSize,'kb_windows.txt',sep="")))[1] # also contains the file path
sampleName = tail(unlist(strsplit(sampleNameTemp,"/")), 1) # filename without path

# plot pie chart into a .jpeg file
jpeg(filename=paste(outputPath,sampleName,"_CoveragePieChart",".jpg",sep=""), width = 680, height = 680)


sortedValues = sort(sampleRawTable$V3)
maxPos = length(sortedValues)

a=signif(quantile(sampleRawTable$V3)[[1]],digits=0) #a=min
b=signif(sortedValues[round(maxPos/8)],digits=0)
c=signif(quantile(sampleRawTable$V3)[[2]],digits=0) #c=q1
d=signif(sortedValues[round(maxPos/8*3)],digits=0)
e=signif(quantile(sampleRawTable$V3)[[3]],digits=0) #e=median
f=signif(sortedValues[round(maxPos/8*5)],digits=0)
g=signif(quantile(sampleRawTable$V3)[[4]],digits=0) #g=q3
h=signif(sortedValues[round(maxPos/8*7)],digits=0)
i=signif(quantile(sampleRawTable$V3)[[5]],digits=0) #i=max

# if elements <10 (except 0 and 1) should not be included:
#constant=c(0,1,Inf)
#oldBreaksDefault = c(a,b,c,d,e,f,g,h,i)
#oldBreaksDefaultUnique=unique(oldBreaksDefault)
#which(oldBreaksDefaultUnique<10)
#oldBreaksDefaultUnique=oldBreaksDefaultUnique[oldBreaksDefaultUnique>=10]
#breaksDefault= sort(cbind(t(constant),t(oldBreaksDefaultUnique)))

breaksDefault=unique(sort(c(0,a,1,b,c,d,e,f,g,h,i,Inf)))

coverageFreqs = table(cut(sampleRawTable$V3, breaks = breaksDefault, right=F))

# check if ranges for coverage frequencies contain enough reads (>1% of total reads), if not resize ranges
lowPercentage = c()
for (i in 1:length(coverageFreqs)) {# aktuell 8 Elemente, also 8 x i
	percent = round(coverageFreqs[[i]]/ sum(coverageFreqs) *100,2) 
	if (percent <= 1) (lowPercentage = c(lowPercentage, i))
}
if (length(lowPercentage) != 0) {
	coverageFreqs = table(cut(sampleRawTable$V3, breaks = unique(breaksDefault)[-lowPercentage], right=F))
	}


newNames = c()
for (i in 1:length(coverageFreqs)) {# aktuell 8 Elemente, also 8 x i
	t=strsplit(names(coverageFreqs[i]),"\\) ")[[1]]
	s=strsplit(t[1],"\\,")[[1]]
	r=strsplit(s[1],"\\[")[[1]]
	u=strsplit(s[2],"\\)")[[1]]

	percent = round(coverageFreqs[[i]]/ sum(coverageFreqs) *100,2) 
	if(r[2]==0){ # umbenennen in 'no coverage'
      		newNames = c(newNames, paste("No coverage\n","(",percent,"%)"))
	}else{
		if(u=="Inf"){
      			newNames = c(newNames, paste(as.numeric(r[2]),"- Inf","\n (",percent,"%)"))
		}else{
		newNames = c(newNames, paste(paste(as.numeric(r[2]),"-",as.numeric(u)-1),"\n (",percent,"%)"))
		}
	}
}

names(coverageFreqs) = newNames
pie(coverageFreqs, clockwise=T,cex=0.7, col=brewer.pal(length(coverageFreqs), "YlOrRd"), main=paste(sampleName,"-",countType,"coverage (in",windowSize,"KB windows)",sep=" "))

dev.off() # close .jpeg file
