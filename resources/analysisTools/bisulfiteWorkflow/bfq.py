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
# methylCtools bfq
# v0.9.2
# 30 march 2012
#
# Current version at https://github.com/hovestadt/methylCtools
#
# extracts fastq file from bam.
# paired end reads are written to the same file.


def mod_bfq(sysargv):
	import sys
	import argparse
	import pysam
	import datetime
	def nicetime(): return datetime.datetime.now().strftime("[bfq %Y-%m-%d %H:%M:%S]")


	#######################################
	# arguments, filehandles

	parser = argparse.ArgumentParser(prog="methylCtools bfq", version="0.9.2", description="bam to fastq converter")
	parser.add_argument("-s", "--silent", dest="qf", action="store_false", help="do not show status messages")

	groupinput = parser.add_argument_group("input files, required")
	groupinput.add_argument("inBAM", metavar="reads.bam", action="store", default=False, help="bam alignment, \"-\" for stdin")
	
	groupoutput = parser.add_argument_group("output files, will be created")
	groupoutput.add_argument("outFQ", metavar="reads.fq", action="store", default=False, help="fastq output, \"-\" for stdout")

	args = parser.parse_args(sysargv)

	try:
		samfileIN = pysam.Samfile(args.inBAM, "rb")
		if args.outFQ == "-": fastqOUT = sys.stdout
		else: fastqOUT = open(args.outFQ, "w")

	except IOError as strerror:
		sys.exit("methylCtools bfq: error: %s" % strerror)
	
	if args.qf: sys.stderr.write("%s command: %s\n" % (nicetime(), " ".join(sys.argv)))
	

	#######################################
	# define functions	
	
	complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
	def rcomp(st):
		return "".join([complement.get(nt, "N") for nt in st[::-1]])		# default: "N", converts all other characters to N


	#######################################
	# main

	if args.qf: sys.stderr.write("%s start: extracting fastq from alignments\n" % nicetime())
	
	while 1:
		try: read = samfileIN.next()
		except StopIteration: break
			
		if read.is_reverse:
			fastqOUT.write("@%s\n%s\n+\n%s\n" % (read.qname, rcomp(read.seq), read.qual[::-1]))
		else:
			fastqOUT.write("@%s\n%s\n+\n%s\n" % (read.qname, read.seq, read.qual))


	#######################################
	# end
			
	samfileIN.close()
	fastqOUT.close()

	if args.qf: sys.stderr.write("%s end: fastq extracted\n" % nicetime())


if __name__ == "__main__":
    import sys
    mod_bfq(sys.argv[1:])
	
