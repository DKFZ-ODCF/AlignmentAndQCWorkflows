#!/usr/bin/python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#
"""
Code to reconstruct correct R1-R2 reads relation from base ratio:

R1R2_layout = (R1_T_C > R2_T_C) and (R2_A_G > R1_A_G)
R2R1_layout = (R1_T_C < R2_T_C) and (R2_A_G < R1_A_G)

read pairs which don't fall into these categories are ignored

Charles Imbusch (c.imbusch@dkfz.de), G200

"""

import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description='Program trying to reconstruct correct R1-R2 reads relation from base ratio in WGBS data.')
parser.add_argument('--R1_in', help='R1 input file', required=True)
parser.add_argument('--R2_in', help='R2 input file', required=True)
parser.add_argument('--gzip_input', help='Set if input fastq files are compressed', action='store_true', required=False, default=False)
parser.add_argument('--R1_out', help='R1 output file', required=True)
parser.add_argument('--R2_out', help='R2 output file', required=True)
parser.add_argument('--R1_unassigned', help='Filename R1 for unassigned read pairs', required=True)
parser.add_argument('--R2_unassigned', help='Filename R2 for unassigned read pairs', required=True)
parser.add_argument('--gzip_output', help='Set if output fastq files should be compressed', action='store_true', required=False, default=False)
parser.add_argument('--log', help='Write basic statistics to this file', required=False, default=False)
parser.add_argument('--debug', help='Debug and write more output to console', action='store_true', required=False, default=False)

args = parser.parse_args()
args_dict = vars(args)

# buffer size when writing out uncompressed fastq files
bufsize = 0

# init some variables
f = None
f2 = None
f_R1_out = None
f_R2_out = None
R1_unassigned = None
R2_unassigned = None

# file handlers for input files
if args_dict['gzip_input'] == True:
    f = open(args_dict['R1_in'], 'r')
    f2 = open(args_dict['R2_in'], 'r')
else:
    f = gzip.open(args_dict['R1_in'], 'r')
    f2 = gzip.open(args_dict['R2_in'], 'r')

# file handlers for output files
if args_dict['gzip_output'] == True:
    f_R1_out = gzip.open(args_dict['R1_out'], 'wb')
    f_R2_out = gzip.open(args_dict['R2_out'], 'wb')
    R1_unassigned = gzip.open(args_dict['R1_unassigned'], 'wb')
    R2_unassigned = gzip.open(args_dict['R2_unassigned'], 'wb')
else:
    f_R1_out = open(args_dict['R1_out'], 'wb', bufsize)
    f_R2_out = open(args_dict['R2_out'], 'wb', bufsize)
    R1_unassigned = open(args_dict['R1_unassigned'], 'wb', bufsize)
    R2_unassigned = open(args_dict['R2_unassigned'], 'wb', bufsize)

if args_dict['log'] == True:
    f_log = open(args_dict['log'], 'wb')

# for some statistics
num_R1R2 = 0
num_R2R1 = 0
num_unsure = 0

### count A/T/C/G content in a given sequence
def getContent(DNAsequence):
    A = 0
    T = 0
    C = 0
    G = 0
    N = 0
    for c in DNAsequence:
        if c == 'A':
            A+=1
        if c == 'T':
            T+=1
        if c == 'C':
            C+=1
        if c == 'G':
            G+=1
        if c == 'N':
            N+=1
    return (A, T, C, G, N)

### count base content masking CG/TG content
def getContentR1(DNAsequence):
    DNAsequence = DNAsequence.replace("CG", "")
    DNAsequence = DNAsequence.replace("TG", "")
    return getContent(DNAsequence)


### count base content masking GC/GT content
def getContentR2(DNAsequence):
        DNAsequence = DNAsequence.replace("GC", "")
        DNAsequence = DNAsequence.replace("GT", "")
        return getContent(DNAsequence)


if args_dict['debug'] == True:
    print "R1_A\tR1_T\tR1_C\tR1_G\tR1_N\tR2_A\tR2_T\tR2_C\tR2_G\tR2_N\tR1_header"

print "Processing reads..."

# do it forever...
while True:
    header1 = f.readline().strip()
    # ...ah only if the first line is not empty we go on
    if not header1:
        break
    sequence1 = f.readline().strip()
    spacer1 = f.readline().strip()
    qual_scores1 = f.readline().strip()
    header2 = f2.readline().strip()
    sequence2 = f2.readline().strip()
    spacer2 = f2.readline().strip()
    qual_scores2 = f2.readline().strip()
    assert header1.split(" ")[0] == header2.split(" ")[0], "headers are not the same: %r" % header1 + "\t" + header2
    R1_content = getContentR1(sequence1)
    R2_content = getContentR2(sequence2)
    ### expert knowledge approach: R1: T:C,A:G versus R2: T:C,A:G;
    R1_T_C = (R1_content[1]+1)/(R1_content[2]+1)
    R1_A_G = (R1_content[0]+1)/(R1_content[3]+1)
    R2_T_C = (R2_content[1]+1)/(R2_content[2]+1)
    R2_A_G = (R2_content[0]+1)/(R2_content[3]+1)
    # for debugging
    if args_dict['debug'] == True:
        print "=="
        print "R1:R2\t" + str(R1_content) + "\t" + str(R2_content) + "\tR1:T:C|R1:A:G\t" + str(R1_T_C) + "\t" +str(R1_A_G)
        print "R2:R1\t" + str(R2_content) + "\t" + str(R1_content) + "\tR2:T:C|R2:A:G\t" + str(R2_T_C) + "\t" +str(R2_A_G)
    # third idea
    #R1R2_layout = ( (R1_T_C > R2_T_C) or (R1_T_C > R1_A_G) ) and ( (R2_A_G > R1_A_G) or (R2_A_G > R2_T_C) )
    #R2R1_layout = ( (R1_T_C < R2_T_C) or (R1_T_C < R1_A_G) ) and ( (R2_A_G < R1_A_G) or (R2_A_G < R2_T_C) )
    # second idea
    R1R2_layout = (R1_T_C > R2_T_C) and (R2_A_G > R1_A_G)
    R2R1_layout = (R1_T_C < R2_T_C) and (R2_A_G < R1_A_G)
    # first idea
    #R1R2_layout = (R1_T_C > R1_A_G) and (R2_A_G > R2_T_C)
    #R2R1_layout = (R1_T_C < R1_A_G) and (R2_A_G < R2_T_C)
    if R1R2_layout:
        R1_output = header1 + "\n" + sequence1 + "\n" + spacer1 + "\n" + qual_scores1 + "\n"
        R2_output = header2 + "\n" + sequence2 + "\n" + spacer2 + "\n" + qual_scores2 + "\n"
        num_R1R2 += 1
        f_R1_out.write(R1_output)
        f_R2_out.write(R2_output)
    if R2R1_layout:
        # swapping sequences but keeping the header the same
        R1_output = header1 + "\n" + sequence2 + "\n" + spacer2 + "\n" + qual_scores2 + "\n"
        R2_output = header2 + "\n" + sequence1 + "\n" + spacer1 + "\n" + qual_scores1 + "\n"
        num_R2R1 += 1
        f_R1_out.write(R1_output)
        f_R2_out.write(R2_output)
    if R1R2_layout==False and R2R1_layout==False:
        # write unassigned read pairs in the same order as they came in
        R1_output = header1 + "\n" + sequence1 + "\n" + spacer1 + "\n" + qual_scores1 + "\n"
        R2_output = header2 + "\n" + sequence2 + "\n" + spacer2 + "\n" + qual_scores2 + "\n"
        R1_unassigned.write(R1_output)
        R2_unassigned.write(R2_output)
        num_unsure += 1

# write basic statistics
reads_total = num_R1R2 + num_R2R1 + num_unsure * 1.0
reads_written = num_R1R2 + num_R2R1
statistics_summary = "number of R1R2 read-pairs:\t" + str(num_R1R2) + ", " + str(num_R1R2/reads_total * 100) + "%\n" + "number of R2R1 read-pairs:\t" + str(num_R2R1) + ", " + str(num_R2R1/reads_total * 100) + "%\n" + "number of unsure read-pairs:\t" + str(num_unsure) + ", " + str(num_unsure/reads_total * 100) + "%\n" + "number of reads-pairs written:\t" + str(reads_written) + ", " + str(reads_written/reads_total * 100) + "%\n" + "number of total reads-pairs:\t" + str(reads_total) + "\n"
if args_dict['log'] == True:
    f_log.write(statistics_summary)

print statistics_summary

# close file handlers
f.close()
f2.close()
f_R1_out.close()
f_R2_out.close()
if args_dict['log'] == True:
    f_log.close()
R1_unassigned.close()
R2_unassigned.close()
