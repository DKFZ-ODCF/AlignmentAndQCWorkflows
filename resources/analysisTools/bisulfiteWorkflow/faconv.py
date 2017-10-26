#!/usr/bin/env python

#######################################
# methylCtools faconv
# v0.9.2
# 30 march 2012
#
# volker hovestadt
# german cancer research center
# v.hovestadt@dkfz.de


def mod_faconv(sysargv):
	import sys
	import argparse
	import datetime
	def nicetime(): return datetime.datetime.now().strftime("[faconv %Y-%m-%d %H:%M:%S]")
	

	#######################################
	# arguments, filehandles
	
	parser = argparse.ArgumentParser(prog="methylCtools faconv", version="0.9.2", description="converts reference sequence C to T (Watson strand) and G to A (Crick strand)")
	parser.add_argument("-s", "--silent", dest="qf", action="store_false", help="do not show status messages")
	
	groupinput = parser.add_argument_group("input files, required")
	groupinput.add_argument("inREF", metavar="ref.fa", action="store", default=False, help="reference sequence (fasta format)")
	
	groupoutput = parser.add_argument_group("output files, will be created")
	groupoutput.add_argument("outCONV", metavar="ref.conv.fa", action="store", default=False, help="converted reference sequence")	
	
	args = parser.parse_args(sysargv)
	
	try:
		inFA = open(args.inREF, "r")
		outFA = open(args.outCONV, "w")
		
	except IOError as strerror:
		sys.exit("methylCtools faconv: error: %s" % strerror)

	if args.qf: sys.stderr.write("%s command: %s\n" % (nicetime(), " ".join(sys.argv)))


	#######################################
	# main

	if args.qf: sys.stderr.write("%s start: converting reference genome\n" % nicetime())

	cname = False
	for line in inFA:
		if line[0] == ">":
			if cname:
				outFA.write(">%s_W\n" % cname)
				outFA.write(seqW)
				outFA.write(">%s_C\n" % cname)
				outFA.write(seqC)
			cname = line.rstrip()[1:]
			seqW = ""
			seqC = ""
			if args.qf: sys.stderr.write("%s status: processing %s\n" % (nicetime(), cname))
		else:
			seqW += line.replace("C", "T").replace("c", "t")
			seqC += line.replace("G", "A").replace("g", "a")

	outFA.write(">%s_W\n" % cname)
	outFA.write(seqW)
	outFA.write(">%s_C\n" % cname)
	outFA.write(seqC)
	

	#######################################
	# end

	inFA.close()
	outFA.close()
	if args.qf: sys.stderr.write("%s end: reference genome converted\n" % nicetime())


if __name__ == "__main__":
    import sys
    mod_faconv(sys.argv[1:])
	
