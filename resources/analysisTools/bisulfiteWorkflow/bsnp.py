#!/usr/bin/env python

#######################################
# methylCtools bsnp
# v0.9.2
# 30 march 2012
#
# volker hovestadt
# german cancer research center
# v.hovestadt@dkfz.de
#
#
# tool for simple genotyping of a predefined set of 40 SNPs.
# SNPs are either A/G or C/T with minor allele frequency ~0.5. 
# returns base counts for A&T and C&G.
#
# based on illumina 450k methylation array typeII rs probes.
# only works for whole-genome (bisulfite) sequencing alignments to hg19.


def mod_bsnp(sysargv):
	import sys
	import pysam
	import argparse

			
	#######################################
	# arguments, filehandles

	parser = argparse.ArgumentParser(prog="methylCtools bsnp", version="0.9.2", description="simple tool for genotyping of a predefined set of 40 SNPs (only for hg19)")

	parser.add_argument("-q", "--mapQ", metavar="INT", dest="mapQ", type=int, default=1, action="store", help="skip alignments with mapQ smaller than INT [1]")
	parser.add_argument("-Q", "--baseQ", metavar="INT", dest="baseQ", type=int, default=20, action="store", help="skip bases with baseQ smaller than INT (Phred+33) [20]")
	
	parser.add_argument("-t", "--trimPE", dest="trimPE", action="store_true", help="only consider read1 in overlapping PE reads")
	parser.add_argument("-b", "--bisulfite", dest="bisulfite", action="store_true", help="bisulfite mode, only one strand is considered")
	
	groupinput = parser.add_argument_group("input files, required")
	groupinput.add_argument("inBAM", metavar="aln.bam", action="store", default=False, help="alignment file (sorted, indexed)")
	
	args = parser.parse_args(sysargv)

	try:										
		samfileIN = pysam.Samfile(args.inBAM, "rb")
		out = sys.stdout
	except IOError as strerror:
		sys.exit("methylCtools bsnp: error: %s" % strerror)


	######################################
	# snp list (typeII rs probes from illumina 450k methylation array)
	
	snps = {}
	snps["rs877309"] = ["chr1", 11489678, "+"]
	snps["rs11249206"] = ["chr1", 25277982, "-"]
	snps["rs3818562"] = ["chr1", 110300441, "+"]
	snps["rs2804694"] = ["chr1", 181331833, "+"]
	snps["rs6426327"] = ["chr1", 246089811, "+"]
	snps["rs1945975"] = ["chr11", 105299724, "-"]
	snps["rs2468330"] = ["chr12", 43198722, "+"]
	snps["rs1495031"] = ["chr12", 63498136, "-"]
	snps["rs2959823"] = ["chr15", 76413529, "+"]
	snps["rs1510189"] = ["chr16", 56100524, "-"]
	snps["rs1941955"] = ["chr18", 35162838, "-"]
	snps["rs966367"] = ["chr2", 12148220, "-"]
	snps["rs4331560"] = ["chr2", 49084673, "+"]
	snps["rs1510480"] = ["chr2", 60825441, "+"]
	snps["rs6546473"] = ["chr2", 69260357, "+"]
	snps["rs264581"] = ["chr2", 159971363, "+"]
	snps["rs2235751"] = ["chr20", 1969934, "+"]
	snps["rs845016"] = ["chr21", 33998284, "-"]
	snps["rs2032088"] = ["chr21", 38477330, "+"]
	snps["rs1467387"] = ["chr22", 25931372, "-"]
	snps["rs133860"] = ["chr22", 26144760, "-"]
	snps["rs739259"] = ["chr22", 27356579, "+"]
	snps["rs2857639"] = ["chr22", 30055674, "+"]
	snps["rs939290"] = ["chr3", 14658866, "-"]
	snps["rs9839873"] = ["chr3", 86662155, "-"]
	snps["rs2125573"] = ["chr4", 128607832, "-"]
	snps["rs7660805"] = ["chr4", 131044050, "+"]
	snps["rs348937"] = ["chr5", 112825677, "-"]
	snps["rs1019916"] = ["chr5", 146648359, "+"]
	snps["rs7746156"] = ["chr6", 47130819, "-"]
	snps["rs9363764"] = ["chr6", 68232042, "+"]
	snps["rs10457834"] = ["chr6", 149042800, "-"]
	snps["rs6982811"] = ["chr8", 36547313, "-"]
	snps["rs1484127"] = ["chr8", 51725654, "+"]
	snps["rs472920"] = ["chr8", 96848744, "+"]
	snps["rs4742386"] = ["chr9", 7711758, "+"]
	snps["rs2521373"] = ["chrX", 9476990, "+"]
	snps["rs798149"] = ["chrX", 15885431, "+"]
	snps["rs5931272"] = ["chrX", 137077857, "+"]
	snps["rs1416770"] = ["chrX", 145219165, "-"]


	#######################################
	# main

	for snp in sorted(snps):
		#print snp, snps[snp]
		p = snps[snp][1]

		for pileupcolumn in samfileIN.pileup(snps[snp][0], start=snps[snp][1]-5, end=snps[snp][1]+3):
			pos = pileupcolumn.pos
			if pos == snps[snp][1]-1:

				if args.trimPE:												# trim overlapping PE reads: check if qname occurs more than once
					once, twice = [], []
					for pileupread in pileupcolumn.pileups:
						if pileupread.alignment.qname in once: twice += [pileupread.alignment.qname]
						else: once += [pileupread.alignment.qname]
						
				c = [0, 0]													# ["red", "green"]
				for pileupread in pileupcolumn.pileups:
					if args.trimPE:											# trim overlapping PE reads
						if pileupread.alignment.qname in twice and pileupread.alignment.is_read2: continue
				
					if args.bisulfite:										# only use one strand in BS mode
						if snps[snp][2] == "+":
							if pileupread.alignment.is_read1 == pileupread.alignment.is_reverse: continue
						else:
							if pileupread.alignment.is_read1 != pileupread.alignment.is_reverse: continue
		
					if pileupread.alignment.qual[pileupread.qpos] < (args.baseQ+33): continue
					if pileupread.alignment.mapq < args.mapQ: continue		# check for minimum mapQ and baseQ
													
					
					if pileupread.alignment.seq[pileupread.qpos] in ["A", "T"]: c[0] += 1
					elif pileupread.alignment.seq[pileupread.qpos] in ["C", "G"]: c[1] += 1
				
				if sum(c) > 0:
					out.write("%s\t%i\t%s\t%s\t%i\t%i\t%.3f\n" % (snps[snp][0], snps[snp][1], snps[snp][2], snp, c[0], c[1], c[1]/float(c[0]+c[1])))


if __name__ == "__main__":
    import sys
    mod_bsnp(sys.argv[1:])

