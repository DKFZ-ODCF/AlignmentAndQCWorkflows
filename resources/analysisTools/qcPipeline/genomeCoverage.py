#! /usr/bin/python
#
# Use this script to
# - determine genome wide read or base coverage by using user-defined window sizes (default window size is 1KB)
# - reads/bases are counted when aligned and not a duplicate
# - coverage analysis for 2 input alignment files possible (typically tumor&normal), in that case plotting logRatio between tumor and normal is possible (in addition to the standard coverage plot per sample) 
# - R script for plotting must be in the current working directory
# - multithreading possible, fastest when using as many CPUs as chromosomes in the given genome
# 
# example call: 
# python genomeCoverage.py --alignmentFile=sample_paired.sorted.bam --outputDir=GenomeCoverage --processors=24 --countReads --windowSize=1 --plotCoverage --ignore_chrRandom_chrM_hap_chrUn
#



import os
import sys
import pysam
import multiprocessing

oneKb = 1000


def sequenceDataProcessing(chrom, queue, alignmentFile, outfileTempName):
    samfile = pysam.Samfile( alignmentFile, "rb" )
    outfileTemp = open(outfileTempName, 'w')
    # determine genomic unit borders for the windows for which read coverage will be determined
    chromUnits = (chromSizeLookup[chrom] / int(options.windowSize * oneKb)) + 1  # plus one more unit for the rest of the chromosome not covered by units of MBs (the modulo rest)
    totalCount = 0
    # analyze a single chromosome in units of multiples (or fractions) of kilobases
    for window in range(chromUnits):
        unitSize = options.windowSize * oneKb
        unitStart = window * unitSize
        if unitStart > chromSizeLookup[chrom]: unitStart = chromSizeLookup[chrom]
        unitEnd = unitStart + unitSize
        if unitEnd > chromSizeLookup[chrom]: unitEnd = chromSizeLookup[chrom]
        outfileTemp.write(chrom.split('.fa')[0] +'\t'+ str(unitStart) + '\t')
        
        readCount = 0
        if options.countReads:
            # start counting alignedReads in the current genome window
            for alignedRead in samfile.fetch(chrom, unitStart, unitEnd):
                if alignedRead.is_duplicate: continue
                if alignedRead.mapq < 1: continue
                if alignedRead.pos < unitStart: continue # do not count read since it was counted in previous window
                readCount += 1
            
        if options.countBases:
            # start counting aligned bases in the current genome window (1 base = 1 pileupcolumn)
            for pileupcolumn in samfile.pileup(chrom, unitStart, unitEnd):
                if pileupcolumn.pos <= unitStart: continue
                if pileupcolumn.pos > unitEnd: continue
                readCount += pileupcolumn.n

        outfileTemp.write(str(readCount) + '\n')
        totalCount += readCount

    if options.baseCoverage: 
        for pileupcolumn in samfile.pileup(chrom):
            outfileTemp.write(chrom.split('.fa')[0] +'\t'+ str(pileupcolumn.pos) + '\t' + str(pileupcolumn.n) + '\n')
            totalCount += pileupcolumn.n
            
    outfileTemp.close()
    samfile.close()
    queue.put(totalCount)
    

#def oufileName(inputAlignmentFile):
#    # prepare output directory and file names
#    if options.countReads: type = 'read'
#    else: type = 'base'
#    windowSize = str(options.windowSize)
#    
#    #sample = os.path.split(inputAlignmentFile)[0]  #TO DO
#    if options.outputDir !='':
#        if not os.path.exists(options.outputDir): os.mkdir(options.outputDir)
#        filenameOnly = os.path.split(inputAlignmentFile)[1]   # returns filename only (without path)
#        sampleName = os.path.splitext(filenameOnly)[0]
#        outfileName = options.outputDir + os.sep +sampleName + '_'+type+'Coverage_'+ windowSize +'kb_windows.txt'
#    else:
#        sampleName = os.path.splitext(inputAlignmentFile)[0]
#        outfileName = sampleName + '_'+type+'Coverage_'+ windowSize +'kb_windows.txt'
#    return outfileName

        
# MAIN ANALYSIS PROCEDURE
def performAnalysis(options):
    global chromSizeLookup
    chromSizeLookup = {}
    chromSortLookup = [] # sorting chromosome names from the SAM header is necessary since some alignment programs sort them differently (for instance, ELAND starts with chr10, while bwa starts with chr1; thus the coverage output files will have a different chr order --> R plotting will fail)
    excludedChromosomeList = options.excludedChromosomes.split(",")
    filesRplot = ''
    totalBases = 0
    
    for alignmentFile in options.alignmentFile:
	io_pair = alignmentFile.rstrip().split(",")
	alignmentFile = io_pair[0]
	output_file = io_pair[1]
        print '\nprocessing alignment file: ' + alignmentFile
        # index the BAM alignment file to enable fast random access (file has to be sorted before indexing)
        if not os.path.exists(alignmentFile+'.bai'):
            sys.stdout.flush()
            print 'indexing BAM file... (BAM file has to be sorted before indexing, otherwise indexing will fail. If you read this message and rest of this program failed, sort BAM file before restarting!)'
            pysam.index(alignmentFile)
        # use the pysam module to read a file in BAM format
        samfile = pysam.Samfile(alignmentFile, "rb")
        for item in samfile.lengths:
            totalBases += item
#        outfile = oufileName(alignmentFile)
	outfile = io_pair[1]
        filesRplot += ' "'+ outfile +'"'
        
        queue = multiprocessing.JoinableQueue()
        filesToCat = ''
        processAliveChecklist = []
        
        # retrieve all chromosomes to be analyzed and their respective sizes
        referenceCounter = 0
        for reference in samfile.references:
            if reference not in chromSizeLookup.keys():
                chrSize = samfile.lengths[referenceCounter]
                chromSizeLookup[reference] = chrSize
                chromSortLookup.append(reference)
            referenceCounter += 1
        chromSortLookup.sort()

        for chrom in chromSortLookup:
            if options.excludedChromosomes != "":
                if chrom in excludedChromosomeList:
                    continue
            if options.ignore_chrRandom_chrM_hap_chrUn:    
                if chrom == "chrM" or chrom == "chrM.fa" or chrom.find("random") >= 0 or chrom.find("hap") >=0 or chrom.find("chrUn") >=0 or chrom.find("GL") >= 0 or chrom.find("JH") >= 0 or chrom.find("KI") >= 0 or chrom == "MT" or chrom=="NC_007605" or chrom == "hs37d5":
                    continue
    
            outfileTempName = outfile + '_'+ chrom
            filesToCat += outfileTempName + ' '
    
            p = multiprocessing.Process(target = sequenceDataProcessing, args=(chrom, queue, alignmentFile, outfileTempName))
            if len(processAliveChecklist) < options.processors:
                p.start()
            processAliveChecklist.append(p)
        samfile.close()

        totalReads = 0
        aliveCounter = 0
        runningProcs = options.processors
        while (aliveCounter != len(processAliveChecklist)):
            for procID in range(len(processAliveChecklist)):
                proc = processAliveChecklist[procID]
                if proc == 'finishedProcess': continue
                elif proc.is_alive(): continue
                elif (proc.pid == None) and (runningProcs < options.processors): # a process that has not been started yet has the process ID 'None'
                    proc.start()
                    runningProcs += 1
                    continue
                elif (proc.pid != None): # a finished process has a process ID
                    readCount = queue.get()
                    queue.task_done()    # received all data from current process via the queue
                    totalReads += readCount
                    proc.join()
                    aliveCounter += 1
                    processAliveChecklist[procID] = 'finishedProcess'
                    runningProcs -= 1
        
        if options.countReads: countType = 'read'
        else: countType = 'base'  
        sys.stdout.flush()
        print 'total number of '+countType+'s counted for coverage analysis:'
        print totalReads

        catCmd = 'cat ' + filesToCat + '> ' + outfile
        sys.stdout.flush()
        print 'Concatenate temp files... ' + catCmd
        os.system(catCmd)
        # delete temp chromosome-specific coverage files 
        sys.stdout.flush()
        os.system('rm ' + outfile + '_*')
        print 'removed temp files...\n'
        
#    if options.plotCoverage:
#        sys.stdout.flush()
#        print '\nplotting coverage using R...\n'
#        filesRplotSplit = filesRplot.split(' "')
#        outputPath = os.path.split(filesRplotSplit[1])[0] + os.sep
#        sys.stdout.flush()
#        os.system('R -f '+ options.Rpath +'/coveragePlotsPerChrom.r --no-save --no-restore --args' + filesRplot +' "'+ countType +'" "'+ str(options.windowSize) +'" "'+ outputPath +'"')
#        
#        filesPieChart = filesRplotSplit[1:len(filesRplotSplit)]
#        for file in filesPieChart:
#            sys.stdout.flush()
#            os.system('R -f '+ options.Rpath +'/coveragePieCharts.r --no-save --no-restore --args "' + file +' "'+ countType +'" "'+ str(options.windowSize) +'" "'+ outputPath +'"')


if __name__ == '__main__':
    print "Starting program..." 
    # constructing command line parser
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--alignmentFile',action='append',type='string',dest='alignmentFile',help='Specify the name of the input and output file(s) (comma separated) containing short read alignments',default=[])
#    parser.add_option('--outputDir',action='store',type='string',dest='outputDir',help='Specify the name of the output directory',default='')
    parser.add_option('--processors',action='store',type='int',dest='processors',help='Specify the number of processors to use',default=1)
    parser.add_option('--windowSize',action='store',type='int',dest='windowSize',help='Specify the size of the window (in multiples of 1KB) for which the read coverage should be determined (default window size: 1KB)',default=1)
    parser.add_option('--ignore_chrRandom_chrM_hap_chrUn',action='store_true',dest='ignore_chrRandom_chrM_hap_chrUn',help='Specify whether reads with extension random, hap or chrM or chrUn should be ignored', default=False)
    parser.add_option('--countReads',action='store_true',dest='countReads',help='Specify whether to count reads within the windowSize (mutually exclusive from counting bases in the window)', default=False)
    parser.add_option('--countBases',action='store_true',dest='countBases',help='Specify whether to count bases within the windowSize (mutually exclusive from counting reads in the window)', default=False)
    parser.add_option('--baseCoverage',action='store_true',dest='baseCoverage',help='Specify whether to determine the genome-wide single base coverage (mutually exclusive from counting reads or bases per window)', default=False)
#    parser.add_option('--plotCoverage',action='store_true',dest='plotCoverage',help='Specify whether to generate a plot using R to show read (or base) coverage per chromosome (R script for plotting must be in the current working directory)', default=False)
#    parser.add_option('--Rpath',action='store',type='string',dest='Rpath',help='Specify the path to the R script for plotting',default='')
    parser.add_option('--excludedChromosomes',action='store',type='string',dest='excludedChromosomes',help='Specify in a comma-separated string which chromosomes should be excluded (it is useful to exclude chrY when aligning to the female genome)', default="")

    (options,args) = parser.parse_args()
    if len(options.alignmentFile) < 1 or len(options.alignmentFile) > 2 or (options.countBases==False and options.countReads==False) or (options.countBases and options.countReads):
        print "Mandatory parameters missing or wrong. Program will terminate now."
        print "\nYour parameter settings:"
        print options        
        raise SystemExit
    performAnalysis(options)     
    print "Program successfully terminating...."  
