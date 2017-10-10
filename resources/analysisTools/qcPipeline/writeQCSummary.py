#! /usr/bin/python
#
# Use this script to
# - get all the different (statistic) values from the QC analysis pipeline
# - write tab-separated output text file (QC summary file) containing all values
#
# example call: python writeQcSummary.py --qcsummaryFile=$ANALYSIS_DIR/QC_summary.txt --pid=${pid} --sample=${sample} --runID=${run} --lane=${lane} --qcDir=$BASE_DIR --filename=$FILENAME
# 


import os
import sys


basesAll = 0

#dpid2samp={}
#for samp,pid in dsamp2pid.items():
#		dpid2samp[pid]=samp
		
		
def listToTabsep(listItems, sep='\t'):
	return sep.join(listItems)

def readFile(filename, interpretationMode):
	text = ""
	try:
		file = open(filename, 'r')
		lineList = file.readlines()
		file.close()
		if interpretationMode == "flagstat":
			text = lineList[0].split(" ")[0] +'\t' + lineList[2].split(" ")[-1].split("%")[0].split("(")[1] +'\t' + lineList[6].split(" ")[-1].split("%")[0].split("(")[1] +'\t' + lineList[8].split(" ")[-1].split("%")[0].split("(")[1] +'\t'
			text += str(int(int(lineList[0].split(" ")[0]) * (float(lineList[2].split(" ")[-1].split("%")[0].split("(")[1]) / 100))) + "\t"
		#elif interpretationMode == "flagstat_unique":
		elif interpretationMode == "metrics":
			info_line = ""
			counter = 0
			for elem in lineList:
				if(elem.split("\t")[0] == "LIBRARY"):
					info_line = lineList[counter+1]
					break
				counter += 1
			splitted_info_line = info_line.rstrip().split("\t")
			text += str(float(splitted_info_line[-2])*100.)+"\t"+splitted_info_line[-1]+"\t"
#		elif interpretationMode == "flagstat_rmdup":
#			readsRmdup = int(int(lineList[0].split(" ")[0]) * (float(lineList[2].split(" ")[-1].split("%")[0].split("(")[1]) / 100))
#			text = str( (unique - readsRmdup) / float(unique) * 100 ) +'\t' + str(unique - readsRmdup)+ '/'+ str(unique) +'\t'
		elif interpretationMode == "structural_variation":
			text = lineList[0].rstrip()+'\t'
		elif interpretationMode == "insertsize_distribution":
			text = lineList[1][0:-1] +'\t'+lineList[0][0:-1] +'\t'
		elif interpretationMode == "coverage":
			coverageQC = lineList[-1][4:-1].split("\t")
			text += coverageQC[14] +'\t' + coverageQC[15] +'\t'+ listToTabsep(coverageQC[0:14]) + '\t'
			for line in lineList:
				if line[0:4] == 'chrX' or line[0] == 'X':
					coverageChrX = line[5:-1].split('\t')[14] +'\t'
					text += coverageChrX
			for line in lineList:
				if line[0:4] == 'chrY' or line[0] == 'Y':
					coverageChrY = line[5:-1].split('\t')[14] +'\n'
					text += coverageChrY
		elif interpretationMode == "coverage_base":
			coverageQC = lineList[-1][4:-1].split("\t")
			text += coverageQC[14] +'\t' + coverageQC[15] +'\t'
			basesAll = int(coverageQC[15].split("/")[0])
		elif interpretationMode == "coverage_target":
			coverageQC = lineList[-1][4:-1].split("\t")
			basesOnTarget = float(coverageQC[15].split("/")[0])
			basesAll = int(coverageQC[15].split("/")[1])
			onTargetRatio=basesOnTarget/basesAll
			text += str(onTargetRatio) +'\t'

			text += coverageQC[14] +'\t' + coverageQC[15] +'\t'+ listToTabsep(coverageQC[2:14]) + '\t'
			for line in lineList:
				if line[0:4] == 'chrX' or line[0] == 'X':
					coverageChrX = line[5:-1].split('\t')[14] +'\t'
					text += coverageChrX
			for line in lineList:
				if line[0:4] == 'chrY' or line[0] == 'Y':
					coverageChrY = line[5:-1].split('\t')[14] +'\n'
					text += coverageChrY
	except IndexError, ei:
		print "Warning: File", filename, " is empty"
	except IOError, e:
		print "Warning: File", filename, " with mode '", interpretationMode, "' could not be found: ", e
	return text

# MAIN ANALYSIS PROCEDURE
def performAnalysis(options):
	# Write Header
	outfile = open(options.qcsummaryFile, 'w')
	outfile.write('PID\tPID_ICGC\tSAMPLE_TYPE\tRUN_ID\tLANE\tTOTAL_READ_COUNT (flagstat)\t%TOTAL_READ_MAPPED_BWA (flagstat)\t%properly_paired (flagstat)\t%singletons (flagstat)\tALIGNED_READ_COUNT\t%DUPLICATES (Picard metrics file)\tESTIMATED_LIBRARY_SIZE (Picard metrics file)\t%PE_reads_on_diff_chromosomes (mapq>0)\t%sd_PE_insertsize (mapq>0)\tPE_insertsize (mapq>0)\tcoverage QC bases w/o N\tQC bases/ total bases w/o N\tcoverage QC bases\tQC bases/ total bases\tmapq=0 read1\tmapq=0 read2\tmapq>0,readlength<minlength read1\tmapq>0,readlength<minlength read2\tmapq>0,BaseQualityMedian<basequalCutoff read1\tmapq>0,BaseQualityMedian<basequalCutoff read2\tmapq>0,BaseQualityMedian>=basequalCutoff read1\tmapq>0,BaseQualityMedian>=basequalCutoff read2\t%incorrect PE orientation\t#incorrect proper pair\t#duplicates read1 (excluded from coverage analysis\t#duplicates read2 (excluded from coverage analysis)\tChrX coverage QC bases\tChrY coverage QC bases\n')

	qcSummaryString = options.pid +'\t'+ options.pid +'\t'+ options.sample+'\t'+ options.runID +'\t'+ options.lane +'\t' 
	
	print
	
	qcSummaryString += readFile(options.filenameFlagstat, "flagstat")
	#qcSummaryString += readFile(options.filenameFlagstatUnique, "flagstat_unique")
	if(options.filenameMetrics == ""):
		qcSummaryString += "NA"+"\t"+"NA"+"\t"
	else:
		qcSummaryString += readFile(options.filenameMetrics, "metrics")
	qcSummaryString += readFile(options.filenameDiffChrom, "structural_variation")
	qcSummaryString += readFile(options.filenameInsertSize, "insertsize_distribution")

	if options.targetCaptureFilename == "":
		
		qcSummaryString += readFile(options.filenameDepthOfC, "coverage")
		
	else:
		# get coverage from the complete BAM file (including reads mapped to non-target regions)
		qcSummaryString += readFile(options.filenameDepthOfC, "coverage_base")
		qcSummaryString += readFile(options.targetCaptureFilename, "coverage_target")

	qcSummaryString += "\n"
	outfile.write(qcSummaryString)
	outfile.close() 


if __name__ == '__main__':
	print "Starting program...\n" 
	print "Arguments: ", len(sys.argv), sys.argv
	
	# constructing command line parser
	import optparse
	parser = optparse.OptionParser()
	parser.add_option('--qcsummaryFile',action='store',type='string',dest='qcsummaryFile',help='Specify the name of the QC summary file to which the results will be written to',default='')
	parser.add_option('--pid',action='store',type='string',dest='pid',help='Specify the pid',default='')
	parser.add_option('--sample',action='store',type='string',dest='sample',help='Specify the sample type (tumor or control)',default='')
	parser.add_option('--runID',action='store',type='string',dest='runID',help='Specify the runID (of the sequencing machine)',default='merged')
	parser.add_option('--lane',action='store',type='string',dest='lane',help='Specify the lane of the sequencing run',default='merged')
	#parser.add_option('--qcDir',action='store',type='string',dest='qcDir',help='Specify the directory containing the QC results per PID',default='')
	#parser.add_option('--filenameBam', action='store', type='string', dest='filenameBam', help='Input bam file', default='')
	parser.add_option('--filenameFlagstat', action='store', type='string', dest='filenameFlagstat', help='Specify flagstat file', default='')
	#parser.add_option('--filenameFlagstatUnique', action='store', type='string', dest='filenameFlagstatUnique', help='Specify unique flagstat file', default='')
	#parser.add_option('--filenameFlagstatRmDup', action='store', type='string', dest='filenameFlagstatRmDup', help='Specify flagstat file with removed duplicates', default='')
	parser.add_option('--filenameDiffChrom', action='store', type='string', dest='filenameDiffChrom', help='Specify file with diffed chroms', default='')
	parser.add_option('--filenameInsertSize', action='store', type='string', dest='filenameInsertSize', help='Specify file with insert sizes', default='')
	parser.add_option('--filenameDepthOfC', action='store', type='string', dest='filenameDepthOfC', help='Specify coverage file', default='')
	parser.add_option('--filenameMetrics', action='store', type='string', dest='filenameMetrics', help='Specify metrics file', default='')
	#parser.add_option('--filename',action='store',type='string',dest='filename',help='Specify the filename',default='')
	parser.add_option('--targetCaptureFilename',action='store',type='string',dest='targetCaptureFilename',help='Specify the targetCapture regions onyl - BAM filename',default='')
		
	
	(options,args) = parser.parse_args()
	
	if options.qcsummaryFile == '':
		print "Mandatory parameters missing or wrong. Program will terminate now."
		print "\nYour parameter settings:"
		print options	   
		raise SystemExit
	performAnalysis(options)	 
	print "\nProgram successfully terminating...."  
