#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

# Rscript ${TOOL_INSERT_SIZE_PLOT_SCRIPT} ${FILENAMED} ${FILENAMEP}_qcValues.txt ${FILENAMEP} "PE insertsize of ${FILEINFO} (rmdup)"
# Rscript ${TOOL_INSERT_SIZE_PLOT_SCRIPT} tumor_ICGC_GBM21_insertsizes.txt tumor_ICGC_GBM21_insertsize_plot.png_qcValues.txt tumor_ICGC_GBM21_insertsize_plot.png "PE insertsize of tumor_ICGC_GBM21 (rmdup)"

cmdArgs = commandArgs(TRUE)

if (length(cmdArgs) < 4)
{
	print(paste("Incorrect number of arguments (4 expected): 1.insert size bin file - 2.file with median and standard deviation values - 3.output plot file (png format) - 4.label"))
	quit()
}
print(cmdArgs)
binfile = cmdArgs[1]
msfile = cmdArgs[2]
outfile = cmdArgs[3]
xlabel = cmdArgs[4]

isizebins=read.table(binfile,sep="\t",as.is=TRUE)

#msfile contains 3 lines: median, standard deviation/median*100 ("stdpercent"), standard deviation
isizevalues=read.table(msfile,sep="\t",as.is=TRUE)
MEDIAN=isizevalues[1,1]
SDpercent=isizevalues[2,1]
SD=isizevalues[3,1]

png(outfile)
plot(isizebins[,1],isizebins[,2],main=paste(binfile,"\nmedian:",MEDIAN,"sd:",SD,"(", SDpercent ,"%)"),xlab=xlabel,ylab="Frequency")
dev.off()
