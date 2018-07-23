#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

# Natalie Jaeger (n.jaeger@dkfz.de)
# version 1.0
# 16 Feb 2012



def snpComp(sysargv):
    usage = "insilicoGenotyper.py <mpileupOf2samples.vcf> <outfile.txt> <MINCOV> <MAXCOV>"
    import sys
    print(sysargv)

    
    if len(sysargv) < 4:                       
        sys.exit("usage: %s" % usage)

    try:                                        
        
        if sysargv[-4] == "-": fin = sys.stdin                                  # mpileup input
        else: fin = open(sysargv[-4], "r")
        
        if sysargv[-3] == "-": fout = sys.stdout                                # SNP discordance output
        else: fout = open(sysargv[-3], "w")   
        
        minCov=int(sysargv[-2])
        maxCov=int(sysargv[-1])

    except IndexError:
        sys.exit("usage: %s" % usage)
    except IOError as (strerror):
        sys.exit("I/O error: %s" % strerror)


    if "-Q" in sysargv: bq = int(sysargv[sysargv.index("-Q")+1])                # skip bases with baseQ smaller than INT (Phred+33), default [13]
    else: bq = 30
    
    
    
    
    def mp2bc(mps, mpq, mr, mq):                                                # count bases in mpileup sequence string
        mps = mps.upper()
        mr = mr.upper() # reference base
        if mr == "A": mrn = 0
        elif mr == "C": mrn = 1
        elif mr == "G": mrn = 2
        else: mrn = 3
        bc = [0, 0, 0, 0]  # [A, C, G, T]
                                                                                                                                                                    
        i1 = 0                                                                    # index for seq
        i2 = 0                                                                    # index for qual
        
        try:
            while i1 < len(mps):                                                # parse mpileup ..
                if mps[i1] in [".", ","]:   # equals reference base
                    if ord(mpq[i2]) >= mq+33:                                    # only consider bases with baseQ >= mq
                        bc[mrn] += 1
                    i1 += 1
                    i2 += 1
                elif mps[i1] == "A":
                    if ord(mpq[i2]) >= mq+33:
                        bc[0] += 1
                    i1 += 1
                    i2 += 1
                elif mps[i1] == "T":
                    if ord(mpq[i2]) >= mq+33:
                        bc[3] += 1
                    i1 += 1
                    i2 += 1
                elif mps[i1] == "C":
                    if ord(mpq[i2]) >= mq+33:
                        bc[1] += 1
                    i1 += 1
                    i2 += 1
                elif mps[i1] == "G":
                    if ord(mpq[i2]) >= mq+33:
                        bc[2] += 1
                    i1 += 1
                    i2 += 1
                elif mps[i1] == "^":                                            # start of read
                    i1 += 2                                                        # .. jump over alignment quality
                elif mps[i1] in ["$", "*"]:                                        # end of read
                    i1 += 1
                elif mps[i1] in ["+", "-"]:                                        # insertion/deletion
                    ni = 1
                    nis = ""
                    while mps[i1+ni] in [`n` for n in range(10)]:                # .. number of inserted/deleted bases (might be >9)
                        nis += mps[i1+ni]
                        ni += 1
                    i1 += int(nis)+1+len(nis)                                    # .. baseQ of inserted bases are not in mpileup 
                elif mps[i1] == "N":
                    i1 += 1
                    i2 += 1
                else:                                                            # should not happen
                    sys.exit(mps)
        except IndexError:
            print mps
            print mpq
        return bc                         # return base count



    def ratioCheck(x):
        if x < 0.3: return 1
        

        
    #######################################
    # MAIN


    totalPos = 0
    matchPos = 0
    
    l = fin.readline()
    while l:
        print(l)
        l = l.rstrip().split("\t")
        
        b1 = mp2bc(l[4], l[5], l[2], bq)        # Control                            # calculate base count
        b2 = mp2bc(l[7], l[8], l[2], bq)        # Tumor
        print(b1, b2)        
        if (not (minCov < sum(b1) < maxCov)) or (not (minCov < sum(b2) < maxCov)): 
            l = fin.readline()
            continue
    
        b1ratio = [elem/float(sum(b1)) for elem in b1]
        b2ratio = [elem/float(sum(b2)) for elem in b2]
        ratioCompzip = zip(b1ratio, b2ratio)
        ratioComp = [abs(elem[1]-elem[0]) for elem in ratioCompzip]

        totalPos += 1
        rCompMap = map(ratioCheck, ratioComp)
        print(rCompMap)
        if rCompMap == [1,1,1,1]:
            matchPos += 1
                    
        l = fin.readline()

    chr = sysargv[-3].split('_')[-1].split('.')[0]
    if not 0 == totalPos:
        genoMatch = (matchPos/float(totalPos)*100)
        fout.write("%s\t%.2f\n" % (chr+" genotype match: "+ str(matchPos)+"/"+ str(totalPos) + " = ", genoMatch))
    
        if genoMatch <= 98:
            foutWarn  = "_".join(sysargv[-3].split('_')[0:-1]) + "_WARNING"
            foutWarning = open(foutWarn, "a") 
            foutWarning.write("%s%.2f%s\n" % ("Genotypes possibly not matching for " + chr + " --> only ", genoMatch, "% genotype match"))
            foutWarning.close()
    
    else:
        foutWarn  = "_".join(sysargv[-3].split('_')[0:-1]) + "_WARNING"
        foutWarning = open(foutWarn, "a") 
        foutWarning.write("%s\n" % ("Genotypes possibly not matching for " + chr + ", position could not be found (totalPos=0)"))
        foutWarning.close()

    fin.close()
    fout.close()

if __name__ == "__main__":
    import sys
    snpComp(sys.argv)

