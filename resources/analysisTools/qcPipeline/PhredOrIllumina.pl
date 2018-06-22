#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#


use strict;
use warnings;

# pipe in a fastq file to check for "diagnostic" base quality score to determine whether it's Phred or Illumina

my $count = 0;
while (<>)
{
	$count++;
	if ($count == 4)	# each 4th line contains base quality scores
	{
		$count = 0;	# reset => more efficient than if ($count % 4 == 3)
		#print;
		if (/[\x22-\x3f]/)	# ASCII 34 - 63 can only be phred
		{
			print "phred\n";
			last;
		}
		elsif (/[\x5e-\x7e]/)	# ASCII 94 - 126 can only be illumina
		{
			print "illumina\n";
			last;
		}
	}
}
# gunzip or awk always complain about broken pipe, so maybe that helps?
# "close" Closes the file or pipe associated with the file handle.
# close (STDIN);
exit;
