#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#
#
# Use this script to
# - determine the genome wide base coverage (will also report base coverage per chromosome)
# - bases are counted when uniquely aligned, and not marked as duplicates
# - multi-threading possible, fastest when using as many CPUs as chromosomes in the given genome
# 
# example call: 
# python coverageQc.py --alignmentFile=${FILENAME_PREFIX}_paired.bam.sorted.unique.dupmark.bam --outputDir=coverage --processors=24 --basequalCutoff=30 --qualityScore=illumina --ignore_chrRandom_chrM_hap_chrUn
#


import os
import sys
import numpy
import pysam
import multiprocessing


manualGenomeCorrection = {"1":"chr1", "2":"chr2", "3":"chr3", "4":"chr4", "5":"chr5", "6":"chr6", "7":"chr7", "8":"chr8", "9":"chr9", "10":"chr10", "11":"chr11", "12":"chr12", "13":"chr13", "14":"chr14", "15":"chr15", "16":"chr16", "17":"chr17", "18":"chr18", "19":"chr19", "20":"chr20", "21":"chr21", "22":"chr22", "X":"chrX", "Y":"chrY", "chrMT":"chrM", "gi|149288852|ref|NC_000067.5|NC_000067 Mus musculus chromosome 1, reference assembly (C57BL/6J)":"chr1","gi|149338249|ref|NC_000068.6|NC_000068 Mus musculus chromosome 2, reference assembly (C57BL/6J)":"chr2","gi|149352351|ref|NC_000069.5|NC_000069 Mus musculus chromosome 3, reference assembly (C57BL/6J)":"chr3","gi|149354223|ref|NC_000070.5|NC_000070 Mus musculus chromosome 4, reference assembly (C57BL/6J)":"chr4","gi|149354224|ref|NC_000071.5|NC_000071 Mus musculus chromosome 5, reference assembly (C57BL/6J)":"chr5","gi|149361431|ref|NC_000072.5|NC_000072 Mus musculus chromosome 6, reference assembly (C57BL/6J)":"chr6","gi|149361432|ref|NC_000073.5|NC_000073 Mus musculus chromosome 7, reference assembly (C57BL/6J)":"chr7","gi|149361523|ref|NC_000074.5|NC_000074 Mus musculus chromosome 8, reference assembly (C57BL/6J)":"chr8","gi|149361524|ref|NC_000075.5|NC_000075 Mus musculus chromosome 9, reference assembly (C57BL/6J)":"chr9","gi|149288869|ref|NC_000076.5|NC_000076 Mus musculus chromosome 10, reference assembly (C57BL/6J)":"chr10","gi|149288871|ref|NC_000077.5|NC_000077 Mus musculus chromosome 11, reference assembly (C57BL/6J)":"chr11","gi|149292731|ref|NC_000078.5|NC_000078 Mus musculus chromosome 12, reference assembly (C57BL/6J)":"chr12","gi|149292733|ref|NC_000079.5|NC_000079 Mus musculus chromosome 13, reference assembly (C57BL/6J)":"chr13","gi|149292735|ref|NC_000080.5|NC_000080 Mus musculus chromosome 14, reference assembly (C57BL/6J)":"chr14","gi|149301884|ref|NC_000081.5|NC_000081 Mus musculus chromosome 15, reference assembly (C57BL/6J)":"chr15","gi|149304713|ref|NC_000082.5|NC_000082 Mus musculus chromosome 16, reference assembly (C57BL/6J)":"chr16","gi|149313536|ref|NC_000083.5|NC_000083 Mus musculus chromosome 17, reference assembly (C57BL/6J)":"chr17","gi|149321426|ref|NC_000084.5|NC_000084 Mus musculus chromosome 18, reference assembly (C57BL/6J)":"chr18","gi|149323268|ref|NC_000085.5|NC_000085 Mus musculus chromosome 19, reference assembly (C57BL/6J)":"chr19","gi|149361525|ref|NC_000086.6|NC_000086 Mus musculus chromosome X, reference assembly (C57BL/6J)":"chrX","gi|149361526|ref|NC_000087.6|NC_000087 Mus musculus chromosome Y, reference assembly (C57BL/6J)":"chrY" }
# from http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/index.shtml: Ungapped lengths are calculated by summing the lenth of the sequenced bases only. No 'Ns' are included in the count.
chromSizesUngapped = {"chr1":225280621, "chr2":238207373, "chr3":194797140, "chr4":187661676, "chr5":177695260, "chr6":167395067, "chr7":155353663, "chr8":142888922, "chr9":120143431, "chr10":131314747, "chr11":131129516, "chr12":130481895, "chr13":95589878, "chr14":88289540, "chr15":81694769, "chr16":78884753, "chr17":77795210, "chr18":74657248, "chr19":55808983, "chr20":59505520, "chr21":35108702, "chr22":34894566, "chrX":151100560, "chrY":25653566, "chrM":0}
# hg19 ungapped total size:  2861327131 without unplaced scaffolds (used here per chrom)
# hg19 ungapped total size:  2863768222 including unplaced scaffolds


# strips off extra whitespace
def striplist(l):
    return([x.strip() for x in l])  


def qualFromASCII(ch):
    return(ord(ch) - qualScoreOffset)


def transformQualStr(s):
    return map(qualFromASCII,s)


def sequenceDataProcessing(chrom, queue, alignmentFile):
    qcResults = []
    baseQualHash = { "mapq=0":{"1":0, "2":0}, "mapq>0,readlength<minlength":{"1":0, "2":0}, "mapq>0,BaseQualityMedian<basequalCutoff":{"1":0, "2":0}, "mapq>0,BaseQualityMedian>=basequalCutoff":{"1":0, "2":0} }
    badOrientationCounter = 0
    totalOrientationCounter = 0
    QCbases = 0
    incorrectProperPairs = 0
    duplicateCount = [0, 0] # first list entry to count duplicates for read1, second one to count read2
    samfile = pysam.Samfile( alignmentFile, "rb" )
    
    for alignedRead in samfile.fetch(chrom):
        if alignedRead.is_duplicate: 
            if alignedRead.is_read1: duplicateCount[0] += 1
            elif alignedRead.is_read2: duplicateCount[1] += 1
            else: duplicateCount[0] += 1  # for single-end reads
            continue

#  rname
#  DEPRECATED from pysam-0.4 - use tid in the future. The rname field caused a lot of confusion as it returns the target ID instead of the reference sequence name.    
        if alignedRead.is_proper_pair:
            if alignedRead.rname == alignedRead.mrnm:
                if alignedRead.is_reverse == alignedRead.mate_is_reverse:
                    # should return NO reads at all, since alignedRead.is_proper_pair does not allow mapping to different chromosomes or to the same strand
                    #print 'Error: read ' +  alignedRead.qname +' is flagged as a proper pair, but is not on opposite strands.'
                    incorrectProperPairs +=1
            
        # check if criterion 4 of 'Leistungsbeschreibung' is fulfilled:
        if not alignedRead.is_unmapped:
            if alignedRead.rname == alignedRead.mrnm:
                if not alignedRead.mate_is_unmapped:
                    totalOrientationCounter += 1
                    if alignedRead.is_reverse == alignedRead.mate_is_reverse:
                        badOrientationCounter += 1
                      
        if alignedRead.is_read1: peRead="1"
        elif alignedRead.is_read2: peRead="2"
        else: peRead="1"  # for single-end reads
        
        if (alignedRead.mapq > 0) and (not alignedRead.is_unmapped):
            if alignedRead.qlen >= options.minReadLength:
                    if round(numpy.mean(transformQualStr(alignedRead.qqual))) >= options.basequalCutoff:
                        baseQualHash["mapq>0,BaseQualityMedian>=basequalCutoff"][peRead]+=1
                        QCbases += alignedRead.qlen
                    else:
                        baseQualHash["mapq>0,BaseQualityMedian<basequalCutoff"][peRead]+=1                    
            else:
                baseQualHash["mapq>0,readlength<minlength"][peRead]+=1
        else:
            baseQualHash["mapq=0"][peRead]+=1
            
    samfile.close()
    qcResults.append(QCbases)


    if totalOrientationCounter != 0: orientationPE = str(round(badOrientationCounter/float(totalOrientationCounter)*100, 3))
    else: orientationPE = 'NA'
    baseQualHashStr = str(baseQualHash["mapq=0"]["1"]) +'\t'+ str(baseQualHash["mapq=0"]["2"]) +'\t'+str(baseQualHash["mapq>0,readlength<minlength"]["1"]) +'\t'+ str(baseQualHash["mapq>0,readlength<minlength"]["2"]) +'\t'+str(baseQualHash["mapq>0,BaseQualityMedian<basequalCutoff"]["1"]) +'\t'+ str(baseQualHash["mapq>0,BaseQualityMedian<basequalCutoff"]["2"]) +'\t'+str(baseQualHash["mapq>0,BaseQualityMedian>=basequalCutoff"]["1"]) +'\t'+ str(baseQualHash["mapq>0,BaseQualityMedian>=basequalCutoff"]["2"]) +'\t'
    resultString = chrom +'\t'+ str(round(float(QCbases)/chromSizeLookup[chrom], 2))+'x'+'\t'+ str(QCbases)+'/'+str(chromSizeLookup[chrom])+'\t'+ baseQualHashStr + orientationPE +'\t' + str(incorrectProperPairs)+ '\t' + str(duplicateCount[0]) +'\t' + str(duplicateCount[1])+'\t'+ str(round(float(QCbases)/chromSizesUngapped[chrom], 2))+'x' + '\t'+ str(QCbases)+'/'+str(chromSizesUngapped[chrom]) +'\n'
    
    qcResults.append(badOrientationCounter)
    qcResults.append(totalOrientationCounter)
    qcResults.append(incorrectProperPairs) # should always be 0
    qcResults.append(resultString)
    qcResults.append(baseQualHash)
    qcResults.append(duplicateCount)
    queue.put(qcResults)
    

        
# MAIN ANALYSIS PROCEDURE
def performAnalysis(options):
    global qualScoreOffset
    if options.qualityScore == 'illumina': qualScoreOffset = 64
    # default would be phred, but set to "null" in Roddy
    # elif options.qualityScore == 'phred': qualScoreOffset = 33
    else: qualScoreOffset = 33
    
    excludedChromosomeList = options.excludedChromosomes.split(",")
    
    
    global chromSizesUngapped
    
    if options.chrLengthExcludingNsFile != '':
        chromLengthFile = open(options.chrLengthExcludingNsFile, "r")
        chrLengthExcludingNsSum=0
        chromSizesUngapped={}
        for line in chromLengthFile:
            if line[0] == '#': continue  # skip header
            lineSplit=line.split('\t')
            lineSplitPlain=striplist(lineSplit)
            chromosome=lineSplitPlain[0]

            # check if to exclude chromosome/ contig / scaffold
            if chromosome in excludedChromosomeList: continue

            chromPlain=chromosome.split('chr')[-1]              # strips off the "chr" if it is the chromosome prefix
            if chromPlain not in options.chromsUsedForCoverageCalc:
                    continue

            chromSizesUngapped[chromosome] = int(lineSplitPlain[1])
            chrLengthExcludingNsSum += int(lineSplitPlain[1])
            
        chromLengthFile.close()
        
    else:
        chromSizesUngapped = {"chr1":225280621, "chr2":238207373, "chr3":194797140, "chr4":187661676, "chr5":177695260, "chr6":167395067, "chr7":155353663, "chr8":142888922, "chr9":120143431, "chr10":131314747, "chr11":131129516, "chr12":130481895, "chr13":95589878, "chr14":88289540, "chr15":81694769, "chr16":78884753, "chr17":77795210, "chr18":74657248, "chr19":55808983, "chr20":59505520, "chr21":35108702, "chr22":34894566, "chrX":151100560, "chrY":25653566, "chrM":16568}
    
    
    
    
    global chromSizeLookup
    chromSizeLookup = {}
    chromSortLookup = [] # sorting chromosome names from the SAM header is necessary since some alignment programs sort them differently (for instance, ELAND starts with chr10, while bwa starts with chr1; thus the coverage output files will have a different chr order --> R plotting will fail)
    
    
    baseQualHashAll = { "mapq=0":{"1":0, "2":0}, "mapq>0,readlength<minlength":{"1":0, "2":0}, "mapq>0,BaseQualityMedian<basequalCutoff":{"1":0, "2":0}, "mapq>0,BaseQualityMedian>=basequalCutoff":{"1":0, "2":0} }
    badOrientationCounterAll = 0
    totalOrientationCounterAll = 0
    incorrectProperPairsAll = 0
    duplicateCount = [0, 0] # first list entry to count duplicates for read1, second one to count read2
    basesPassedQc = 0
    referenceBases = 0    # 3095693983 #source NCBI 37 Chr1-22,X,Y,M counted with fasta_length.pl 

    for alignmentFile in options.alignmentFile:
        # index the BAM alignment file to enable fast random access (file has to be sorted before indexing)
        if not os.path.exists(alignmentFile+'.bai'):
            sys.stdout.flush()
            print 'indexing BAM file... (BAM file has to be sorted before indexing, otherwise indexing will fail. If you read this message and rest of this program failed, sort BAM file before restarting!)'
            pysam.index(alignmentFile)
        # use the pysam module to read a file in BAM format
        samfile = pysam.Samfile(alignmentFile, "rb")
        outfileName = options.outputDir + os.sep + os.path.split(alignmentFile)[1] + '.DepthOfCoverage.txt_temp'

        outfileHeader = open(outfileName + '.header', "w")
        fileHeader = 'interval\tcoverage QC bases\t#QC bases/#total bases\tmapq=0 read1\tmapq=0 read2\tmapq>0,readlength<minlength read1\tmapq>0,readlength<minlength read2\tmapq>0,BaseQualityMedian<basequalCutoff read1\tmapq>0,BaseQualityMedian<basequalCutoff read2\tmapq>0,BaseQualityMedian>=basequalCutoff read1\tmapq>0,BaseQualityMedian>=basequalCutoff read2\t%incorrect PE orientation\t#incorrect proper pair\t#duplicates read1 (excluded from coverage analysis)\t#duplicates read2 (excluded from coverage analysis)\tgenome_w/o_N coverage QC bases\t#QC bases/#total not_N bases\n'
        outfileHeader.write(fileHeader)
        outfileHeader.close()
        
        queue = multiprocessing.JoinableQueue()
        processAliveChecklist = []
        
        # retrieve all chromosomes to be analyzed and their respective sizes in bp
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
    
            chromPlain=chrom.split('chr')[-1]              # strips off the "chr" if it is the chromosome prefix
            if chromPlain not in options.chromsUsedForCoverageCalc:
                    continue

            if options.ignore_chrRandom_chrM_hap_chrUn:    
   		if chrom == "chrM" or chrom == "chrM.fa" or chrom.find("random") >= 0 or chrom.find("hap") >=0 or chrom.find("chrUn") >=0 or chrom.find("GL") >= 0 or chrom.find("JH") >= 0 or chrom == "MT" or chrom=="NC_007605" or chrom == "hs37d5":
                    continue

            referenceBases += chromSizeLookup[chrom]
    
            p = multiprocessing.Process(target = sequenceDataProcessing, args=(chrom, queue, alignmentFile))
            if len(processAliveChecklist) < options.processors:
                p.start()
            processAliveChecklist.append(p)
        samfile.close()
        
        outfile = open(outfileName, "w")

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
                    qcResultsPerChrom = queue.get()
                    queue.task_done()    # received all data from current process via the queue
                    # write QC results from chromosome to output file
                    basesPassedQc += qcResultsPerChrom[0]
                    badOrientationCounterAll += qcResultsPerChrom[1]
                    totalOrientationCounterAll += qcResultsPerChrom[2]
                    incorrectProperPairsAll += qcResultsPerChrom[3]
                    outfile.write(qcResultsPerChrom[4])
                    
                    baseQualHashChrom = qcResultsPerChrom[5]
                    for criterion in baseQualHashChrom.keys():
                        for peRead in baseQualHashChrom[criterion].keys():
                            baseQualHashAll[criterion][peRead] += baseQualHashChrom[criterion][peRead]
                    
                    duplicateCount[0] += qcResultsPerChrom[6][0]
                    duplicateCount[1] += qcResultsPerChrom[6][1]

                    proc.join()
                    aliveCounter += 1
                    processAliveChecklist[procID] = 'finishedProcess'
                    runningProcs -= 1
            
        outfile.close()
        sys.stdout.flush()
        os.system('sort '+outfileName+' > '+outfileName+'.sorted_temp')
        os.system('cat '+outfileName + '.header '+outfileName+'.sorted_temp > '+outfileName)
        os.system('rm '+outfileName + '.header '+outfileName+'.sorted_temp')
        
        if totalOrientationCounterAll != 0: percentBadOrientation = round(badOrientationCounterAll/float(totalOrientationCounterAll)*100, 3)
        else: percentBadOrientation = 'NA'
        badOrientationCount = str(percentBadOrientation)
        #if percentBadOrientation > 0.2:
            #badOrientationCount = str(percentBadOrientation) + ' - criterion 4 not fulfilled!'
        
        baseQualHashAllStr = str(baseQualHashAll["mapq=0"]["1"]) +'\t'+ str(baseQualHashAll["mapq=0"]["2"]) +'\t'+str(baseQualHashAll["mapq>0,readlength<minlength"]["1"]) +'\t'+ str(baseQualHashAll["mapq>0,readlength<minlength"]["2"]) +'\t'+str(baseQualHashAll["mapq>0,BaseQualityMedian<basequalCutoff"]["1"]) +'\t'+ str(baseQualHashAll["mapq>0,BaseQualityMedian<basequalCutoff"]["2"]) +'\t'+str(baseQualHashAll["mapq>0,BaseQualityMedian>=basequalCutoff"]["1"]) +'\t'+ str(baseQualHashAll["mapq>0,BaseQualityMedian>=basequalCutoff"]["2"]) +'\t'
        
        outfile = open(outfileName, "a")
        
        basesForCoverageCalculation=2861327131
        if options.chrLengthExcludingNsFile != '': basesForCoverageCalculation = chrLengthExcludingNsSum
        if options.targetEnrichmentSize > 0: basesForCoverageCalculation = options.targetEnrichmentSize
        if options.TruSeqExomeEnrichment: basesForCoverageCalculation=62000000      # see http://www.illumina.com/products/truseq_exome_enrichment_kit.ilmn
        if options.costomMBEnrichment: basesForCoverageCalculation=13712280
        
        outfile.write('all\t' +str(round(float(basesPassedQc)/referenceBases, 2))+'x'+'\t'+str(basesPassedQc)+' / '+str(referenceBases)+'\t'+ baseQualHashAllStr + badOrientationCount + '\t'+str(incorrectProperPairsAll) + '\t' + str(duplicateCount[0]) +'\t' + str(duplicateCount[1]) +'\t'+ str(round(float(basesPassedQc)/basesForCoverageCalculation, 2))+'x' + '\t'+ str(basesPassedQc)+'/'+str(basesForCoverageCalculation)+'\n')   # write results for whole genome
        outfile.close()


if __name__ == '__main__':
    print "Starting program..." 
    # constructing command line parser
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--alignmentFile',action='append',type='string',dest='alignmentFile',help='Specify the name of the input file containing short read alignments (duplicates should either be marked or removed, non-unique reads should either be removed)',default=[])
    parser.add_option('--outputDir',action='store',type='string',dest='outputDir',help='Specify the name of the output directory',default='')
    parser.add_option('--processors',action='store',type='int',dest='processors',help='Specify the number of processors to use',default=1)
    parser.add_option('--qualityScore',action='store',type='string',dest='qualityScore',help='Specify whether the per base  quality score is given in phred or illumina format (default is Illumina score: ASCII offset of 64, while PHRED scores have an ASCII offset of 33)', default='phred')
    parser.add_option('--basequalCutoff',action='store',type='int',dest='basequalCutoff',help='Specify the minimum mean of quality scores of the complete read length used as a quality cutoff (default: 30 (in phred score))',default=30)
    parser.add_option('--minReadLength',action='store',type='int',dest='minReadLength',help='Specify the minimum read length (default: 36)',default=36)
    parser.add_option('--chromsUsedForCoverageCalc',action='append',type='string',dest='chromsUsedForCoverageCalc',help='Specify which chromosomes to include for coverage calculation, other chromosomes/ contigs / scaffolds will be ignored for coverage calculation, like chrM or hs37d5 for example. Default is chr1-22,X,Y. For mouse the same default can be used since non-existing chroms will just be skipped. If BAM file contains "chr" or not does not matter, always use chrom IDs without the "chr" prefix.', default=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"])
    parser.add_option('--ignore_chrRandom_chrM_hap_chrUn',action='store_true',dest='ignore_chrRandom_chrM_hap_chrUn',help='Specify whether reads with extension random, hap or chrM or chrUn should be ignored', default=False)
    parser.add_option('--excludedChromosomes',action='store',type='string',dest='excludedChromosomes',help='Specify in a comma-separated string which chromosomes should be excluded (it is useful to exclude chrY when aligning to the female genome)', default='')
    parser.add_option('--TruSeqExomeEnrichment',action='store_true',dest='TruSeqExomeEnrichment',help='Set true when Illumia TruSeq Exome Enrichment Kit was used', default=False)
    parser.add_option('--costomMBEnrichment',action='store_true',dest='costomMBEnrichment',help='Set true when Illumia the 2700 genes costomMBEnrichment was used', default=False)
    parser.add_option('--targetEnrichmentSize',action='store',type='int',dest='targetEnrichmentSize',help='Specify the number of bases used for enrichment, value will be used for coverage calculation here', default=0)
    parser.add_option('--chrLengthExcludingNsFile',action='store',type='string',dest='chrLengthExcludingNsFile',help='Specify the name including full path of the file containing chromosome lengths for each chromosome excluding Ns; tab-separated text file with two columns, first column chrom name, second column chrom length, if header in file, then use # as start of header)',default='')
    #parser.add_option('--chrLengthIncludingNsFile',action='store',type='string',dest='chrLengthIncludingNsFile',help='Specify the name including full path of the file containing chromosome lengths for each chromosome including Ns as bases; tab-separated text file with two columns, first column chrom name, second column chrom length, if header in file, then use # as start of header)',default='')

    (options,args) = parser.parse_args()
    if (len(options.alignmentFile) < 1) or (options.outputDir == ''):
        print "Mandatory parameters missing or wrong. Program will terminate now."
        print "\nYour parameter settings:"
        print options        
        raise SystemExit
    performAnalysis(options)     
    print "Program successfully terminating...."  
