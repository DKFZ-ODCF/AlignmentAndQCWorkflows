#!/usr/bin/python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#


# Convert simple tab separated file with header to JSON.
import argparse
import sys
import json
from gccorrection_python_modules import Tabfile

parser=argparse.ArgumentParser(description='Convert tab separated CNV calls to vcf file' )
parser.add_argument('--file',	  	'-f', type=file, help='Tabfile with headerline')
parser.add_argument('--id', 	  	'-i', type=str, help='Patient identifier printed in json')
parser.add_argument('--out', 	  	'-o', type=str, help='Output file')
parser.add_argument('--key', 	  	'-k', type=str, help='key for 1st level')

parser.add_argument
args = parser.parse_args()



if __name__=='__main__':

	if not args.id or not args.file:
		sys.stderr.write("ERROR: Please specify patient id and input file. For more information, use -h.\n\n\n")
		sys.exit(2)

	if not args.out:
		out=sys.stdout
	else:
		try:
			out=open(args.out, 'w')
		except IOError as (errno, strerr ):
			sys.stderr.write("WARNING: Specified outputfile cannot be written. Please check given path.\n")
			sys.exit("IOError %i:%s\n" % (errno, strerr))


	try:
		tabfile = Tabfile.Input(args.file)
	except IOError as (errno, strerr ):
		sys.exit("IOError %i:%s\n" % (errno, strerr))


	for line in tabfile:
		if "pid" in line:
			del line['pid']
#		jsonID = { args.id : line }
		jsonMain = { args.key : line }
		out.write( json.dumps(jsonMain, indent=2, separators=(",",":")) )

