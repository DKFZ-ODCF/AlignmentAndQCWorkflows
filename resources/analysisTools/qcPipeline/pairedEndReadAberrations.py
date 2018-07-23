#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

#
# Use this script to
# - find reads of a PE read pair that have been mapped to different chromosomes; writes tab-separated output text file
#
# example call: 
# python2.6 pairedEndReadAberrations.py --alignmentFile=s_1_paired.sorted.bam --outputFile=/mapcomp/peReadsOnDiffChrms.txt --printSummaryToStdOut --plotDistribution
#
# contact n.jaeger@dkfz.de for further questions


import os
import sys
import pysam

# MAIN ANALYSIS PROCEDURE
def performAnalysis(options):
    excludedChromosomeList = options.excludedChromosomes.split(",")    
    alignmentFile = options.alignmentFile[0]
    outfile = open(options.outputFile, 'w')
    # index the BAM alignment file to enable fast random access
    if not os.path.exists(alignmentFile+'.bai'):
        pysam.index(alignmentFile)
    # use the pysam module to read a file in BAM format
    samfile = pysam.Samfile(alignmentFile, "rb")
    # retrieve all chromosomes to be analyzed
    totalReadCount = 0
    readCount = 0
    for chrom in samfile.references:
        if options.excludedChromosomes != "":
            if chrom in excludedChromosomeList:
                continue
        if options.ignore_chrRandom_chrM_hap_chrUn:    
            if chrom == "chrM" or chrom.find("random") >= 0 or chrom.find("hap") >=0 or chrom.find("chrUn") >=0:
                continue
        
        for alignedRead in samfile.fetch(chrom):
            if not alignedRead.is_unmapped: totalReadCount += 1
            if options.mapqFiltered:
                if alignedRead.mapq == 0: continue
            
            if not alignedRead.is_proper_pair: # reads is not part of a proper PE read pair
                if not alignedRead.is_unmapped: # but read is mapped
                    if alignedRead.mrnm != (-1): # and mate is mapped, too
                        if alignedRead.rname != alignedRead.mrnm: # but both reads have been mapped to different chromosomes (= reference names)               
                            readCount += 1
                            outfile.write("%s\t%s\t%i\t%s\t%s\t%i\n" % (alignedRead.qname, samfile.getrname(alignedRead.rname),alignedRead.pos, 'mate:', samfile.getrname(alignedRead.mrnm), alignedRead.mpos))
    samfile.close()
    outfile.close() 
    
    if options.printSummaryToStdOut:
        sys.stdout.flush()
        print options.alignmentFile[0]
        print '#reads of a PE pair mapped to different chromosomes:'
        print readCount
        print '# total reads in file:'
        print totalReadCount
        
    if options.plotDistribution:
        sys.stdout.flush()
        if os.path.split(options.outputFile)[0] != '': outputPath = os.path.split(options.outputFile)[0]+ os.sep
        else: outputPath = ''
        os.system('R -f '+ options.Rpath +'/plotStructuralVariation.R --no-save --no-restore --args "'+options.outputFile+'" "'+outputPath+'" "'+str(totalReadCount)+'"')


if __name__ == '__main__':
    print "Starting program...\n" 
    # constructing command line parser
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--alignmentFile',action='append',type='string',dest='alignmentFile',help='Specify the name of the input file containing short read alignments',default=[])
    parser.add_option('--outputFile',action='store',type='string',dest='outputFile',help='Specify the name of the output file',default='')
    parser.add_option('--plotDistribution',action='store_true',dest='plotDistribution',help='Specify whether to plot the distribution of structural variation as a histogram using R', default=False)
    parser.add_option('--Rpath',action='store',type='string',dest='Rpath',help='Specify the path to the R script for plotting',default='')
    parser.add_option('--mapqFiltered',action='store_true',dest='mapqFiltered',help='Specify whether to filter out reads with mapping quality of 0 (which means non-unique alignment for bwa)', default=False)
    parser.add_option('--printSummaryToStdOut',action='store_true',dest='printSummaryToStdOut',help='Specify whether to print the total insertion and deletion count summary to standard out (i.e. the .o file for cluster jobs)', default=False)
    parser.add_option('--excludedChromosomes',action='store',type='string',dest='excludedChromosomes',help='Specify in a comma-separated string which chromosomes should be excluded (it is useful to exclude chrY when aligning a female genome)', default="")
    parser.add_option('--ignore_chrRandom_chrM_hap_chrUn',action='store_true',dest='ignore_chrRandom_chrM_hap_chrUn',help='Specify whether reads with extension random, hap or chrM or chrUn should be ignored', default=False)
    
    (options,args) = parser.parse_args()
    if len(options.alignmentFile) < 1 or options.outputFile == '':
        print "Mandatory parameters missing or wrong. Program will terminate now."
        print "\nYour parameter settings:"
        print options        
        raise SystemExit
    performAnalysis(options)     
    print "\nProgram successfully terminating...."  
