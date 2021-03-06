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
# methylCtools fapos
# v0.9.2
# 30 march 2012
#
# Current version at https://github.com/hovestadt/methylCtools


def mod_fapos(sysargv):
	import sys
	import argparse
	import re
	import datetime
	def nicetime(): return datetime.datetime.now().strftime("[fapos %Y-%m-%d %H:%M:%S]")


	#######################################
	# arguments, filehandles
	
	parser = argparse.ArgumentParser(prog="methylCtools fapos", version="0.9.2", description="creates cytosine position index for defined context")
	parser.add_argument("-s", "--silent", dest="qf", action="store_false", help="do not show status messages")
	
	groupinput = parser.add_argument_group("input files, required")
	groupinput.add_argument("inREF", metavar="ref.fa", action="store", default=False, help="reference sequence (fasta format)")
	
	groupoutput = parser.add_argument_group("output files, will be created")
	groupoutput.add_argument("outPOS", metavar="ref.pos", action="store", default=False, help="cytosine position index (tab format, \"-\" for stdout)")
	
	groupcontext = parser.add_mutually_exclusive_group()
	groupcontext.add_argument("-1", dest="c1", action="store_true", help="output positions for C context")
	groupcontext.add_argument("-2", dest="c2", action="store_true", help=".. CG context")
	groupcontext.add_argument("-3", dest="c3", action="store_true", help=".. CG/CH context (default)")
	groupcontext.add_argument("-4", dest="c4", action="store_true", help=".. CG/CHG/CHH context")
	
	args = parser.parse_args(sysargv)
	
	try:
		inFA = open(args.inREF, "r")
		if args.outPOS == "-": outPOS = sys.stdout
		else: outPOS = open(args.outPOS, "w")
		
	except IOError as strerror:
		sys.exit("methylCtools fapos: error: %s" % strerror)

	t = [args.c1, args.c2, args.c3, args.c4]
	if not sum(t): t[2] = True

	if args.qf: sys.stderr.write("%s command: %s\n" % (nicetime(), " ".join(sys.argv)))


	#######################################
	# functions

	complement = {"A":"T", "C":"G", "G":"C", "T":"A", "a":"t", "c":"g", "g":"c", "t":"a"}
	def rcomp(st):
		return "".join([complement.get(nt, "N") for nt in st[::-1]])		# default: "N", converts all other characters to N
			
	def pp(i, l):
		if seq[i].upper() == "C":
			if i+2 < l:
				if t[0]:
					outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "+", "C", seq[i:(i+3)]))
				elif seq[i+1].upper() == "G":
					outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "+", "CG", seq[i:(i+3)]))
				elif t[2]:
					outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "+", "CH", seq[i:(i+3)]))
				elif t[3]:
					if seq[i+2].upper() == "G":
						outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "+", "CHG", seq[i:(i+3)]))
					else:
						outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "+", "CHH", seq[i:(i+3)]))	
		else:
			if i-2 >= 0:
				if t[0]:
					outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "-", "C", rcomp(seq[(i-2):(i+1)])))
				elif seq[i-1].upper() == "C":
					outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "-", "CG", rcomp(seq[(i-2):(i+1)])))
				elif t[2]:
					outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "-", "CH", rcomp(seq[(i-2):(i+1)])))
				elif t[3]:
					if seq[i-2].upper() == "C":
						outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "-", "CHG", rcomp(seq[(i-2):(i+1)])))
					else:
						outPOS.write("%s\t%i\t%s\t%s\t%s\n" % (cname, i, "-", "CHH", rcomp(seq[(i-2):(i+1)])))		
	

	#######################################
	# main
							
	if t[0]: tc = "C"
	elif t[1]: tc = "CG"
	elif t[2]: tc = "CG/CH"
	else: tc = "CG/CHG/CHH"
	if args.qf: sys.stderr.write("%s start: generating position index for %s context\n" % (nicetime(), tc))
	
	cname = False
	p = re.compile("[C,G]", re.IGNORECASE)
	for line in inFA:
		if line[0] == ">":
			if cname:
				seqlen = len(seq)
				for m in p.finditer(seq):
					pp(m.start(), seqlen)
			cname = line.rstrip()[1:]
			seq = ""
			if args.qf: sys.stderr.write("%s status: processing %s\n" % (nicetime(), cname))
		else:
			seq += line.rstrip()
	
	seqlen = len(seq)
	for m in p.finditer(seq):
		pp(m.start(), seqlen)


	#######################################
	# end

	inFA.close()
	outPOS.close()

	if args.qf: sys.stderr.write("%s end: position index generated\n" % nicetime())


if __name__ == "__main__":
    import sys
    mod_fapos(sys.argv[1:])

