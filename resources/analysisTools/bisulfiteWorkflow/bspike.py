#!/usr/bin/env python

#######################################
# methylCtools bspike
# v0.9.2
# 30 march 2012
#
# volker hovestadt
# dkfz (german cancer research center)
# v.hovestadt@dkfz.de


def mod_bspike(sysargv):
	import sys
	import argparse
	import pysam
	import datetime
	def nicetime(): return datetime.datetime.now().strftime("[bspike %Y-%m-%d %H:%M:%S]")
	

	#######################################
	# arguments, filehandles
	
	parser = argparse.ArgumentParser(prog="methylCtools bspike", version="0.9.2", description="extract base information from spiked-in DNA (e.g. lambda phage genome)")
	parser.add_argument("-s", "--silent", dest="qf", action="store_false", help="do not show status messages")

	parser.add_argument("-q", "--mapQ", metavar="INT", dest="aq", type=int, default=1, action="store", help="skip alignments with mapQ smaller than INT [1]")
	parser.add_argument("-Q", "--baseQ", metavar="INT", dest="bq", type=int, default=20, action="store", help="skip bases with baseQ smaller than INT (Phred+33) [20]")
	parser.add_argument("-r", "--region", metavar="STR", dest="spikechr", action="store", default="chrL", help="spike fasta id [chrL]")
	
	groupinput = parser.add_argument_group("input files, required")
	groupinput.add_argument("fafile", metavar="ref.fa", action="store", default=False, help="reference sequence (indexed)")
	groupinput.add_argument("samfileIN", metavar="aln.bam", action="store", default=False, help="alignment file (sorted, indexed)")
	
	groupoutput = parser.add_argument_group("output files, will be created")
	groupoutput.add_argument("out", metavar="aln.spike", action="store", default=False, help="spike base information")	
	
	args = parser.parse_args(sysargv)
	bq = args.bq + 33
	
	try:
		fafile = pysam.Fastafile(args.fafile)
	
		samfileIN = pysam.Samfile(args.samfileIN, "rb")
		
		if args.out == "-":
			out = sys.stdout
		else:
			out = open(args.out, "w")
						
	except IOError as strerror:
		sys.exit("methylCtools bspike: error: %s" % strerror)
	
	if args.qf: sys.stderr.write("%s command: %s\n" % (nicetime(), " ".join(sys.argv)))				
	

	#######################################
	# functions
	
	def mapq(t):
		if ord(t[1]) < bq: return "N"
		else: return t[0]
	
	complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
	def comp(x):
		return complement.get(x, "N")
	

	#######################################
	# loop file
	
	out.write("chr\tpos\tref\tA.m1W\tA.m1C\tA.m2W\tA.m2C\tC.m1W\tC.m1C\tC.m2W\tC.m2C\tG.m1W\tG.m1C\tG.m2W\tG.m2C\tT.m1W\tT.m1C\tT.m2W\tT.m2C\n")
	f = fafile.fetch(reference=args.spikechr).upper()
	for pileupcolumn in samfileIN.pileup(reference=args.spikechr):
		out.write("%s\t%i\t%s" % (args.spikechr, pileupcolumn.pos , f[pileupcolumn.pos]))
		mp = dict(zip(["A", "C", "G", "T"], [[0] *4, [0] *4, [0] *4, [0] *4, [0] *4]))
		
		for pileupread in pileupcolumn.pileups:
			if pileupread.alignment.mapq < args.aq: continue
			if ord(pileupread.alignment.qual[pileupread.qpos]) < bq: continue

			if pileupread.alignment.is_read1:
				if not pileupread.alignment.is_reverse:
					mp[pileupread.alignment.seq[pileupread.qpos]][0] += 1			# 99
				else:
					mp[comp(pileupread.alignment.seq[pileupread.qpos])][1] += 1		# 83
			elif pileupread.alignment.is_read2:
				if pileupread.alignment.is_reverse:
					mp[comp(pileupread.alignment.seq[pileupread.qpos])][2] += 1		# 147
				else:
					mp[pileupread.alignment.seq[pileupread.qpos]][3] += 1			# 163
		
		for b in ["A", "C", "G", "T"]:
			out.write("\t%i\t%i\t%i\t%i" % (mp[b][0], mp[b][1], mp[b][2], mp[b][3]))
		out.write("\n")
	
	out.close()	
	samfileIN.close()


if __name__ == "__main__":
    import sys
    mod_bspike(sys.argv)

