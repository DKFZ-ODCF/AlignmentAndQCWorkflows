# example script call: R -f insertSizeDistribution.r --no-save --no-restore --args "/ibios/temp1/2010-06-07QCGATC/insertsize_distribution/tumor/s_1_insertsizes.txt" "/ibios/temp1/2010-06-07QCGATC/insertsize_distribution/tumor/s_1_isize_testplot.png" "isize tumor s1 gatc"
# OBTAIN CONFIGURATION PARAMETERS FOR PLOTTING THE INSERT SIZE DISTRIBUTION (WITH R)

cmdArgs = commandArgs(TRUE)
#print(cmdArgs)
if (length(cmdArgs) != 3) print(paste("Incorrect number of arguments (3 expected): ",length(cmdArgs)))

flatfilepath = cmdArgs[1]   
outputfilePathAndName = cmdArgs[2] 
xlab1 = cmdArgs[3]

#sample <- read.table(flatfilepath)
# read from the pipe instead!
sample <- read.table("stdin")
MEDIAN=median(abs(sample[,1]))
SD=round(sd(abs(sample[,1])))
#SDpercent=round(sd(abs(sample[,1]))/median(abs(sample[,1]))*100)
# argh, we have the values already, why calculate them again?!
SDpercent=round(SD/MEDIAN*100)
png(outputfilePathAndName)
hist(sample[,1], breaks=seq(min(sample), max(sample)+10, by=10),main=paste(outputfilePathAndName,"\nmedian:",MEDIAN,"sd:",SD,"(", SDpercent ,"%)"), xlab=xlab1, cex.main=0.7)
dev.off()

# write the statistic values to a simple text file
sink(paste(outputfilePathAndName,"_qcValues.txt",sep=""), append=FALSE, split=FALSE)
cat(paste(MEDIAN, "\n"))
cat(paste(SDpercent, "\n"))
sink()