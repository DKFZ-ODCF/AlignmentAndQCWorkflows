#!/usr/bin/env python
#
# MIT License
#
# Copyright (c) 2018 Volker Hovestadt
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
# methylCtools bcall
# v0.9.3
# 18 july 2012
#
# Current version at https://github.com/hovestadt/methylCtools
#
# calls methylation ratio for inidividual cytosines from BS sequencing.
# input is sorted/indexed bam file.
# also requires cytosine position index as input (created by methylCtools fapos).
#
# writes one row per cytosine position in defined region.
# output columns are chromosome name, coordinate, strand, context/genomic sequence, snp rate, unconverted cytosines, converted cytosines.
# use "-" to write to stdout. can be piped to bgzip/tabix.
#
# methylCtools examines the available sequencing information to evaluate if a cytosine position differs from the reference (argument -x).
# for example, if in a sample the true genotype of a reference CpG would be TpG, it would likely appear as an unmethylated CpG since they cannot be distinguished
# from bisulfite sequencing data. however, on the opposite strand there would be clear evidence for this SNV. the algorithm examines all positions within the
# context (in this case CG, or similarly C/CH/CHG/CHH) on both strands and reports on the maximum rate of evidence from any of these positions.
# positive evidence for a SNV would be:
#   for C: sense A/G,   antisense A/T/G
#   for G: sense A/C/T, antisense C/T
#   for H: sense G,     antisense G


def mod_bcall(sysargv):
	import sys
	import argparse
	from argparse import RawTextHelpFormatter
	import pysam
	import datetime
	def nicetime(): return datetime.datetime.now().strftime("[bcall %Y-%m-%d %H:%M:%S]")
	
			
	#######################################
	# arguments, filehandles

	parser = argparse.ArgumentParser(prog="methylCtools bcall", version="0.9.2", description="calls single-basepair methylation ratios from aligned BS sequencing data", formatter_class=RawTextHelpFormatter)
	parser.add_argument("-s", "--silent", dest="qf", action="store_false", help="do not show status messages")

	parser.add_argument("-t", "--trimPE", dest="to", action="store_true", help="only consider read1 in overlapping PE reads")
	parser.add_argument("-q", "--mapQ", metavar="INT", dest="aq", type=int, default=1, action="store", help="skip alignments with mapQ smaller than INT [1]")
	parser.add_argument("-Q", "--baseQ", metavar="INT", dest="bq", type=int, default=20, action="store", help="skip bases with baseQ smaller than INT (Phred+33) [20]")
	parser.add_argument("-e", "--skipend", metavar="INT", dest="se", type=int, default=0, action="store", help="skip bases within INT bp from the 5' or 3' end of the read [0]")
	
	parser.add_argument("-x", "--snv", dest="mmr", action="store_true", help="perform SNV detection")
	parser.add_argument("-y", "--snvcov", metavar="INT", dest="mmrc", type=int, default=5, action="store", help="require INT fold coverage for SNV detection [5]")
	parser.add_argument("-f", "--snvfilter", metavar="FLOAT", dest="mmrf", type=float, default=False, action="store", help="do not write positions with SNV rate >= FLOAT [off]")

	parser.add_argument("-g", "--genomic", dest="ge", action="store_true", help="write 3nt genomic sequence instead of context")
	parser.add_argument("-z", "--zero", dest="pz", action="store_true", help="write positions without coverage")
	parser.add_argument("-r", "--region", metavar="STR", dest="reg", action="store", default=False, help="process positions in specified region [all sites]")
		
	parser.add_argument("-m", "--metrics", metavar="metrics.txt", dest="me", action="store", default=False, help="write metrics output")

	groupinput = parser.add_argument_group("input files, required")
	groupinput.add_argument("inPOS", metavar="ref.pos.gz", action="store", default=False, help="cytosine position index (tabix indexed)")
	groupinput.add_argument("inBAM", metavar="sample.bam", action="store", default=False, help="alignment file (sorted, indexed)")
	
	groupoutput = parser.add_argument_group("output files, will be created")
	groupoutput.add_argument("outMETH", metavar="sample.call", action="store", default=False, help="methylation calls (tab format, \"-\" for stdout)")

	args = parser.parse_args(sysargv)

	try:								
		tabixPositions = pysam.Tabixfile(args.inPOS, "r")						# positions
		samfileIN = pysam.Samfile(args.inBAM, "rb")
		if args.outMETH == "-": out = sys.stdout								# output
		else: out = open(args.outMETH, "w")
		if args.me: meout = open(args.me, "w")									# metrics output

	except IOError as strerror:
		sys.exit("methylCtools bcall: error: %s" % strerror)

	mmrff = False
	if type(args.mmrf) is not type(False):
		mmrff = True
		if args.mmrf > 1 or args.mmrf < 0:
			sys.exit("methylCtools bcall: error: -f must be between 0 and 1.")	# mismatch rate filter

	ge = args.ge + 1 # if it is 1 takes the next column for context dp[p][1] is CH/CG, 2 is e.g. CAG

	reg = args.reg
	reg_ref = False
	reg_range = False
	if reg:
		reg = "".join(reg.split(","))											# remove 1000 separator (,)
		if ":" not in reg and "-" not in reg:
			reg_ref = reg
		elif ":" in reg and "-" in reg:
			reg_ref, reg2 = reg.split(":")
			reg_range = [int(r) for r in reg2.split("-")]
		else:
			sys.exit("methylCtools bcall: error: " + reg + " is not a valid region.")

	if not reg_ref: reg_ref_list = sorted(tabixPositions.contigs)				# use all contigs if -r is not defined
	else:
		if not reg_ref in tabixPositions.contigs:
			sys.exit("methylCtools bcall: error: " + reg_ref + " not present in position index.")
		reg_ref_list = [reg_ref]
	
	if args.qf: sys.stderr.write("%s command: %s" % (nicetime(), " ".join(sys.argv)))
	

	#######################################
	# functions
	
	def qpos2bpos(read, qp):
		if not read.is_reverse: return qp
		else: return read.rlen-qp -1

	def dometrics(m,qf): 															# m==0: unconverted, m==1: converted
		metqual[int(pr.alignment.is_read2)][dp[p][1]][bqb-33][m] += 1
		if qf == True: # not filtered by quality
			metqpos[int(pr.alignment.is_read2)][dp[p][1]][qpos2bpos(pr.alignment, pr.qpos)][m] += 1

	#######################################
	# main

	window = 1000000															# window size
	
	# for C: sense A/G,   antisense A/T/G
	# for G: sense A/C/T, antisense C/T
	# for H: sense G,     antisense G
	illWs = {"C": ["A", "G"], "G": ["A", "C", "T"], "H": ["G"]}					# context
	illWa = {"C": ["A", "T", "G"], "G": ["C", "T"], "H": ["G"]}
	illW = {"C": [[illWs["C"]], [illWa["C"]]],
		"CG": [[illWs["C"], illWs["G"]], [illWa["C"], illWa["G"]]],
		"CH": [[illWs["C"], illWs["H"]], [illWa["C"], illWa["H"]]],
		"CHG": [[illWs["C"], illWs["H"], illWs["G"]], [illWa["C"], illWa["H"], illWa["G"]]],
		"CHH": [[illWs["C"], illWs["H"], illWs["H"]], [illWa["C"], illWa["H"], illWa["H"]]]}
	illCs = {"C": ["T", "C"], "G": ["T", "G", "A"], "H": ["C"]}
	illCa = {"C": ["T", "A", "C"], "G": ["G", "A"], "H": ["C"]}
	illC = {"C": [[illCs["C"]], [illCa["C"]]],
		"CG": [[illCs["C"], illCs["G"]], [illCa["C"], illCa["G"]]],
		"CH": [[illCs["C"], illCs["H"]], [illCa["C"], illCa["H"]]],
		"CHG": [[illCs["C"], illCs["H"], illCs["G"]], [illCa["C"], illCa["H"], illCa["G"]]],
		"CHH": [[illCs["C"], illCs["H"], illCs["H"]], [illCa["C"], illCa["H"], illCa["H"]]]}

	### loop over chromosomes/contigs
	for c in reg_ref_list:
		
		if args.me:																# empty structures for metrics
			metqpos = [dict(zip(["C", "CG", "CH", "CHG", "CHH"], [[[0, 0] for x in range(151)], [[0, 0] for x in range(151)], [[0, 0] for x in range(151)], [[0, 0] for x in range(151)], [[0, 0] for x in range(151)]])),
				   dict(zip(["C", "CG", "CH", "CHG", "CHH"], [[[0, 0] for x in range(151)], [[0, 0] for x in range(151)], [[0, 0] for x in range(151)], [[0, 0] for x in range(151)], [[0, 0] for x in range(151)]]))]
			metqual = [dict(zip(["C", "CG", "CH", "CHG", "CHH"], [[[0, 0] for x in range(43)], [[0, 0] for x in range(43)], [[0, 0] for x in range(43)], [[0, 0] for x in range(43)], [[0, 0] for x in range(43)]])),
				   dict(zip(["C", "CG", "CH", "CHG", "CHH"], [[[0, 0] for x in range(43)], [[0, 0] for x in range(43)], [[0, 0] for x in range(43)], [[0, 0] for x in range(43)], [[0, 0] for x in range(43)]]))]
			metcov = dict(zip(["C", "CG", "CH", "CHG", "CHH"], [[0]*201, [0]*201, [0]*201, [0]*201, [0]*201]))

		if args.qf: sys.stderr.write("\n%s start: calling methylation status" % nicetime())
		if args.qf: sys.stderr.write("\n%s status: processing %s" % (nicetime(), c))
		clen = samfileIN.lengths[samfileIN.gettid(c)]
		if args.me:																# metrics per contig
			met = {"C": [0, 0], "CG": [0, 0], "CH": [0, 0], "CHG": [0, 0], "CHH": [0, 0]}
						
		if not reg_range:
			reg_range_split = range(1, clen, window)
			if clen not in reg_range_split: reg_range_split += [clen]
		else:
			reg_range_split = range(reg_range[0], reg_range[1], window)
			if reg_range[1] not in reg_range_split: reg_range_split += [reg_range[1]]
	
		### loop over windows
		for i in range(len(reg_range_split)-1):
			reg_range_window = reg_range_split[i:(i+2)]
			if args.qf: sys.stderr.write("\n%s status: loading positions %s:%i-%i" % (nicetime(), c, reg_range_window[0], reg_range_window[1]-1))	# 1-indexed
			dp = dict([(int(t[1]), [t[2], t[3], t[4], 0, 0, 0]) for t in tabixPositions.fetch(reference=c, start=reg_range_window[0]-1, end=reg_range_window[1]-1, parser=pysam.asTuple())])	# 0-indexed
			#if args.qf: sys.stderr.write("\n%s status: %i positions loaded" % (nicetime(), len(dp)))
			
			if dp:
				### loop over positions
				for pileupcolumn in samfileIN.pileup(reference=c, start=min(dp), end=max(dp)):
					p = pileupcolumn.pos										# dp 1-indexed, pileup 0-indexed
					### check if in position list
					if p in dp:
						
						### check for overlapping PE reads
						if args.to:
							# may make the search more efficient if it is dictionary instead of list structure
							once, twice = {}, {}
							for pr in pileupcolumn.pileups:
								if pr.alignment.mapq < args.aq: continue
								if once.has_key(pr.alignment.qname): twice[pr.alignment.qname]=''
								else: once[pr.alignment.qname]=''
						
						### check context
						cx = dp[p][1]	# it is 3nt context I think or 2nt?												# context
						if not args.mmr: cx = "C"										# if no mmr caluclation (SNP calculation) test only C
						cxlen = len(cx)													# context length
						illcount = [[[0]*cxlen, [0]*cxlen], [[0]*cxlen, [0]*cxlen]]		# context illegal base count [ok][strand][q]
						ctcount = [0, 0]												# methylation [C, T]
						if dp[p][0] == "+":												# position strand
							### Watson
							ill = illW[cx]												# [strand][q]
							### loop over reads
							for pr in pileupcolumn.pileups:
								if pr.alignment.mapq < args.aq: continue
								if pr.alignment.is_paired:								# ok if single end read
									if not pr.alignment.is_proper_pair:					# ok if proper pair
										if not pr.alignment.mate_is_unmapped:			# ok if mate unmapped
											if not (pr.alignment.tid == pr.alignment.mrnm and pr.alignment.pos == pr.alignment.mpos):	# ok if sharing same pos (bwa does not flag these as proper pairs)
												continue								# otherwise continue
								if args.to:												# check for overlapping PE reads
									if  twice.has_key(pr.alignment.qname) and pr.alignment.is_read2: continue
								SKIP=False	
								if args.se:												# skip bases close to end of read
									if pr.qpos < args.se or pr.alignment.rlen-pr.qpos <= args.se: SKIP=True									
								strand = int(pr.alignment.is_reverse != pr.alignment.is_read2)	# 0: sense, 1: antisense
								### loop over context within read
								for q in range(cxlen):									# q: position in context
									if pr.qpos+q >= pr.alignment.qend: continue			# check for end of read (also softclip)
									bqb = ord(pr.alignment.qual[pr.qpos+q])				# check for base quality
									if bqb < (args.bq+33):
										qf = False
										if not args.me: continue						# continue if no metrics need to be generated
									else: qf = True
									if pr.alignment.seq[pr.qpos+q] in ill[strand][q]:	# check for illegal bases
										if qf: illcount[1][strand][q] += 1				# 1: illegal, 0: ok
									else:
										if qf: illcount[0][strand][q] += 1
										if q == 0 and not strand:
											if pr.alignment.seq[pr.qpos] == "C":
												if qf and not SKIP : ctcount[0] += 1
												if args.me: dometrics(0,qf)
											else:										# already checked for illegal bases, only C/T is possible
												if qf and not SKIP : ctcount[1] += 1
												if args.me: dometrics(1,qf)
												
						else:
							### Crick
							ill = illC[cx]
							for pr in pileupcolumn.pileups:
								if pr.alignment.mapq < args.aq: continue
								if pr.alignment.is_paired:
									if not pr.alignment.is_proper_pair:
										if not pr.alignment.mate_is_unmapped:
											if not (pr.alignment.tid == pr.alignment.mrnm and pr.alignment.pos == pr.alignment.mpos):
												continue
								if args.to:
									if twice.has_key(pr.alignment.qname) and pr.alignment.is_read2: continue
								SKIP=False	
								if args.se:
									if pr.qpos < args.se or pr.alignment.rlen-pr.qpos <= args.se: SKIP=True		
								strand = int(pr.alignment.is_reverse == pr.alignment.is_read2)
								for q in range(cxlen):
									if pr.qpos-q < pr.alignment.qstart: continue		# check for start of read (also softclip)
									bqb = ord(pr.alignment.qual[pr.qpos-q])
									if bqb < (args.bq+33):
										qf = False
										if not args.me: continue
									else: qf = True
									if pr.alignment.seq[pr.qpos-q] in ill[strand][q]:
										if qf: illcount[1][strand][q] += 1
									else:
										if qf: illcount[0][strand][q] += 1
										if q == 0 and not strand:
											if pr.alignment.seq[pr.qpos] == "G":
												if qf and not SKIP: ctcount[0] += 1
												if args.me : dometrics(0,qf)
											else:
												if qf and not SKIP: ctcount[1] += 1
												if args.me: dometrics(1,qf)				
						
						### calcluate max mismatch rate
						maxmmr = 0
						for q in range(cxlen):									# context illegal base count [ok][strand][q]
							for strand in [0, 1]:
								illcountsum = illcount[0][strand][q] + illcount[1][strand][q]
								if illcountsum >= max(1, args.mmrc):			# mmr coverage check
									illcountrate = illcount[1][strand][q]/float(illcountsum)
									if illcountrate > maxmmr:
										maxmmr = illcountrate
						
						dp[p][3:] = [maxmmr, ctcount[0], ctcount[1]]			# append mmr and C/T count to dp

				### write for all positions in window at once
				for p in sorted(dp):											# sort by position
					if args.pz or dp[p][4]+dp[p][5]:
						if args.mmr:											# mismatch rate calcluated
							if mmrff:
								if dp[p][3] < args.mmrf:						# pass mismatch rate filter (mmr not shown)
									if args.me:
										met[dp[p][1]][0] += dp[p][4]
										met[dp[p][1]][1] += dp[p][5]
										metcov[dp[p][1]][min(dp[p][4]+dp[p][5], 200)] += 1
									out.write("%s\t%i\t%s\t%s\t%i\t%i\n" % (c, p, dp[p][0], dp[p][ge], dp[p][4], dp[p][5]))
							else:												# no filter (mmr shown)
								if args.me:
									met[dp[p][1]][0] += dp[p][4]
									met[dp[p][1]][1] += dp[p][5]
									metcov[dp[p][1]][min(dp[p][4]+dp[p][5], 200)] += 1
								out.write("%s\t%i\t%s\t%s\t%.2f\t%i\t%i\n" % (c, p, dp[p][0], dp[p][ge], dp[p][3], dp[p][4], dp[p][5]))
						else:													# no mmr calculated; option turned off
							if args.me:
								met[dp[p][1]][0] += dp[p][4]
								met[dp[p][1]][1] += dp[p][5]
								metcov[dp[p][1]][min(dp[p][4]+dp[p][5], 200)] += 1						
							out.write("%s\t%i\t%s\t%s\t%i\t%i\n" % (c, p, dp[p][0], dp[p][ge], dp[p][4], dp[p][5]))
								
		### write metrics out
		if args.me:
			if sum(met['C'] + met['CG'] + met['CH'] + met['CHG'] + met['CHH']) > 0:
				meout.write("### global methylation\n")
				meout.write("region\tcontext\tmC\tC\tratio\n")
				cxs = []
				for cx in ["C", "CG", "CH", "CHG", "CHH"]:
					if sum(met[cx]):
						cxs += [cx]
						if not reg_range:
							meout.write("%s\t%s\t%i\t%i\t%.6f\n" % (c, cx, met[cx][0], met[cx][1], met[cx][0]/float(sum(met[cx]))))
						else:
							meout.write("%s:%i-%i\t%s\t%i\t%i\t%.6f\n" % (c, reg_range[0], reg_range[1], cx, met[cx][0], met[cx][1], met[cx][0]/float(sum(met[cx]))))
						
				meout.write("\n### methylation vs. position\n")
				meout.write("mate\tpos")
				for cx in cxs:
					meout.write("\t%s.mC\t%s.C" % (cx, cx))
				meout.write("\n")
				for m in [0, 1]:
					for q in range(151):											# break if true read length < 150
						if sum([sum(x) for x in metqpos[m][cxs[0]][q:]]) == 0: break
						meout.write("%i\t%i" % (m+1, q+1))
						for cx in cxs:
							meout.write("\t%i\t%i" % (metqpos[m][cx][q][0], metqpos[m][cx][q][1]))
						meout.write("\n")
					
				meout.write("\n### methylation vs. baseQ\n")			
				meout.write("mate\tbaseQ")
				for cx in cxs:
					meout.write("\t%s.mC\t%s.C" % (cx, cx))
				meout.write("\n")
				for m in [0, 1]:
					for q in range(42):
						meout.write("%i\t%i" % (m+1, q))
						for cx in cxs:
							meout.write("\t%i\t%i" % (metqual[m][cx][q][0], metqual[m][cx][q][1]))
						meout.write("\n")
					
				meout.write("\n### coverage\n")			
				meout.write("cov")
				for cx in cxs:
					meout.write("\t%s" % cx)
				meout.write("\n")
				for cov in range(201):
					meout.write("%i" % cov)
					for cx in cxs:
						meout.write("\t%i" % metcov[cx][cov])
					meout.write("\n")
				meout.write("\n")
															
	#tabixPositions.close()														# no tabix.close in pysam ..
	samfileIN.close()
	out.close()
	if args.me: meout.close()

	if args.qf: sys.stderr.write("\n%s end: methylation status called\n" % nicetime())


if __name__ == "__main__":
    import sys
    mod_bcall(sys.argv[1:])

