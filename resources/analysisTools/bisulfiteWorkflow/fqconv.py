#!/usr/bin/env python

#######################################
# methylCtools fqconv
# v0.9.2
# 30 march 2012
#
# volker hovestadt
# german cancer research center
# v.hovestadt@dkfz.de
#
#
# converts illumina reads from BS sequencing to 3-letter alphabet.
# reads must be in fastq format. supports gzip-compressed input files.
# to read an uncompressed file from stdin, use "-".
# to write to stdout, use "-".
#
# for SE: read is C to T converted.
# for PE: 1st read is C to T converted, 2nd read is G to A converted.
# read ids are trimmed to the first occurrence of # or space character.
# conversion positions are appended to read id (MD-tag style).
# if the resulting read id becomes longer than 250 chars (violates SAM
# format specifications, happens if there are many converted positions,
# sequencing artifact), the read is not converted.


def mod_fqconv(sysargv):
	import sys
	import argparse
	import gzip
	import datetime
	def nicetime(): return datetime.datetime.now().strftime("[fqconv %Y-%m-%d %H:%M:%S]")
	
	
	#######################################
	# arguments, filehandles

	parser = argparse.ArgumentParser(prog="methylCtools fqconv", version="0.9.2", description="converts bisulfite sequencing reads to fully converted state")
	parser.add_argument("-s", "--silent", dest="qf", action="store_false", help="do not show status messages")

	groupinput = parser.add_argument_group("input files, required")
	groupinput.add_argument("inFQ", metavar="reads.fq", action="store", default=False, help="sequencing reads (fastq format, \"-\" for stdin, supports gzip compressed files)")

	groupoutput = parser.add_argument_group("output files, will be created")
	groupoutput.add_argument("outFQ", metavar="reads.conv.fq", action="store", default=False, help="converted sequencing reads, \"-\" for stdout")

	groupcontext = parser.add_mutually_exclusive_group()
	groupcontext.add_argument("-1", dest="m1", action="store_true", help="single end reads or first read in pair (default)")
	groupcontext.add_argument("-2", dest="m2", action="store_true", help="second read in pair")

	args = parser.parse_args(sysargv)
			
	if args.m1: cfrom, cto = "C", "T"
	elif args.m2: cfrom, cto = "G", "A"
	else: cfrom, cto = "C", "T"
				
	try:
		if args.inFQ.split(".")[-1].lower() in ["gz", "gzip"]:
			in1 = gzip.open(args.inFQ, "rb")
		else:
			if args.inFQ == "-": in1 = sys.stdin
			else: in1 = open(args.inFQ, "r")
			
		if args.outFQ == "-": out1 = sys.stdout
		else: out1 = open(args.outFQ, "w")
		
	except IOError as strerror:
		sys.exit("methylCtools fqconv: error: %s" % strerror)

	if args.qf: sys.stderr.write("%s command: %s\n" % (nicetime(), " ".join(sys.argv)))


	#######################################
	# main

	n = 0
	id1 = ""
	seq1 = ""

	if args.qf: sys.stderr.write("%s start: converting reads %s to %s\n" % (nicetime(), cfrom, cto))
	c = [0, 0]																# [reads, bases converted]
	wc = 0																	# warning count
	
	o1 = ""
	for l1 in in1:
		l1 = l1.rstrip()
		n += 1
				
		if n == 1:															# line 1: name
			c[0] += 1
			id1 = l1.split()[0].split("#")[0]								# BWA only uses id up to first space, clip # (not picard-compatible?)
			
		elif n == 2:														# line 2: sequence
			seq1 = l1
																			# (line 3: quality name, replace to "+")		
		elif n == 4:														# line 4: quality
			qual1 = l1

			id1 += "."														# convert sequence and append positions to id (MD-tag style)
			scount = 0
			for s in seq1:
				if s == cfrom:
					if scount > 0: id1 += "%i" % scount
					scount = 0
					id1 += "%s" % cfrom
					c[1] += 1
				else: scount += 1
			if scount > 0: id1 += "%i" % scount
			
			if len(id1) >= 250:												# if id too long
				id1 = id1.split(".")[0] + "."
				id1 += "%i" % len(seq1)
				if wc < 100: sys.stderr.write("%s warning: %s is not converted (%s)\n" % (nicetime(), id1, seq1))
				if wc == 99: sys.stderr.write("%s warning: only showing 100 warnings\n" % nicetime())
				wc += 1
			else:
				seq1 = seq1.replace(cfrom, cto)								# convert
			
			o1 += ("%s\n%s\n+\n%s\n" % (id1, seq1, qual1))
			n = 0															# reset counter

			if c[0]%10**6 == 0:
				out1.write(o1)												# write
				o1 = ""
				if args.qf: sys.stderr.write("%s status: %i reads processed\n" % (nicetime(), c[0]))

	out1.write(o1)
		

	#######################################
	# end

	in1.close()
	out1.close()

	if args.qf: sys.stderr.write("%s end: %i reads processed, %i bases converted\n" % (nicetime(), c[0], c[1]))


if __name__ == "__main__":
    import sys
    mod_fqconv(sys.argv[1:])

