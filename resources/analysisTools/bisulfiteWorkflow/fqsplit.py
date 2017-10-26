#!/usr/bin/env python

#######################################
# methylCtools bfq
# v0.9.2
# 30 march 2012
#
# volker hovestadt
# german cancer research center
# v.hovestadt@dkfz.de


def mod_fqsplit(sysargv):
	import sys
	import argparse
	import gzip
	import datetime
	def nicetime(): return datetime.datetime.now().strftime("[fqsplit %Y-%m-%d %H:%M:%S]")
	
	
	#######################################
	# arguments, filehandles
	
	parser = argparse.ArgumentParser(prog="methylCtools fqsplit", version="0.9.2", description="splits fastq file into multiple fastq files of defined size")
	parser.add_argument("-s", "--silent", dest="qf", action="store_false", help="do not show status messages")
	parser.add_argument("-l", "--reads", metavar="INT", dest="reads", type=int, default=10000000, action="store", help="split into INT reads per output file [10000000]")
	
	groupinput = parser.add_argument_group("input files, required")
	groupinput.add_argument("inFQ", metavar="reads.fq", action="store", default=False, help="input file (fastq format, \"-\" for stdin, supports gzip-compression)")
	
	groupoutput = parser.add_argument_group("output files, will be created")
	groupoutput.add_argument("outFQ", metavar="reads.split.fq", action="store", default=False, help="output file (creates multiple files, supports gzip-compression)")
	
	args = parser.parse_args(sysargv)

	gzflag = False
	if args.outFQ.split(".")[-1].lower() in ["gz", "gzip"]: gzflag = True
	nlc = 0
	
	try:
		if args.inFQ.split(".")[-1].lower() in ["gz", "gzip"]:
			in1 = gzip.open(args.inFQ, "rb")
		else:
			if args.inFQ == "-": in1 = sys.stdin
			else: in1 = open(args.inFQ, "r")

	except IOError as (errno, strerror):
		sys.exit("I/O error: %s" % strerror)

	if args.qf: sys.stderr.write("%s command: %s\n" % (nicetime(), " ".join(sys.argv)))


	#######################################
	# main

	if args.qf: sys.stderr.write("%s start: splitting into files containing %i reads each\n" % (nicetime(), args.reads))
	n = 0	
	for l1 in in1:
		n += 1
		if n%(args.reads*4) == 1:
			if nlc > 0: out1.close()
			nlc += 1
			try:
				if gzflag: out1 = gzip.open(".".join(args.outFQ.split(".")[:-2] + [format(nlc, "04d")] + [args.outFQ.split(".")[-2]] + [args.outFQ.split(".")[-1]]), "wb", 1)
				else: out1 = open(".".join(args.outFQ.split(".")[:-1] + [format(nlc, "04d")] + [args.outFQ.split(".")[-1]]), "w")
			except IOError as (errno, strerror):
				sys.exit("I/O error: %s" % strerror)
			if args.qf: sys.stderr.write("%s status: writing file %04i\n" % (nicetime(), nlc))
			
		out1.write(l1)
		

	#######################################
	# end

	in1.close()
	out1.close()

	if args.qf: sys.stderr.write("%s end: %i reads processed, %i files created\n" % (nicetime(), n/4, nlc))


if __name__ == "__main__":
    import sys
    mod_fqsplit(sys.argv)
	
	