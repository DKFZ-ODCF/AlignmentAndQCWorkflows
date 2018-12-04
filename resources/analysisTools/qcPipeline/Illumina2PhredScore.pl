#!/usr/bin/perl 
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#
# Barbara Hutter 2011, see http://seqanswers.com/forums/showthread.php?t=8164

use strict;
use warnings;


if (@ARGV < 1)
{
	die "USAGE: $0 Illumina _sequence.txt file to convert into PHRED score.\n";
}
my $file = shift;

open (R, $file) or die "Could not open $file\n";
# count all lines and reads
my $ctr = 0;
my $reads = 0;
# every 4th line is a quality string
my $c4 = 0;

#my @help = ();
#my $qual = "";
#my $line = "";

while (<R>)
{
	$c4++;
	$ctr++;
	#if ($reads >= 3){last;}
	if ($c4 == 4)	# quality score line, transform it
	{
		chomp;	# otherwise it would transform the newline into "!"
		#$line = $_;
		#print STDERR "Original: $_\n";
		# suggestion of nicolallias at http://seqanswers.com/forums/showthread.php?t=8164
		# $_ =~ tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/;
		$_ =~ tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/;
		print $_, "\n";

		# Which could be written in Python as
		# q_line = "".join([chr(ord(i)-31) for i in q_line])
		# or in Perl:
		#@help = split('', $line);
		#$qual = "";
		#foreach (@help)
		#{
		#	$qual.= chr(ord($_)-31);
		#}
		#print STDERR "$qual\n";
		# since above $_ overwrote the $_ that contained the original line string:
		#$line =~ tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/;
		#print STDERR "shortcut: ", $line, "\n";
		
		#print "Reconverted:";
		#@help = split('', $_);
		#$qual = "";
		#foreach (@help)
		#{
		#	print STDERR (ord($_)-33);
		#	print STDERR " ";
		#}
		#print STDERR "\n";

		# reset the flag
		$c4 = 0;
		# count reads
		$reads++;
	}
	else	# print the current non-quality line
	{
		print;
	}
}
close R;
print STDERR "$ctr lines, $reads reads\n";
exit;
