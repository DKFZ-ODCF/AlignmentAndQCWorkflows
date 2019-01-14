#!/usr/env perl

use strict;
use warnings;

# pipe in a fastq file to discard "Illumina chastity filter failed reads"


# this is a normal read, with "N" (N stands for did not fail)
#@39V34V1:78:C2UEDACXX:1:1101:1274:2170 1:N:0:CAGATC
#TCTGATGCCCTCTTCTGGTGTGTCTGAAGACAGCTACCGGGTACTTACAT
#+
#BBBFFFFFFFFFFIIIIIFFIIIIIIFIIIIIIIIFIIIIIBFFFIIIII

# this is a failed read (Y for yes it did fail)
#@39V34V1:78:C2UEDACXX:1:1101:1337:2204 1:Y:0:CAGATC
#TCCGCTCTAGTAGTGTTCGAGGCGTTTCCCTCTCGGGGAGCTGGGGGCGG
#+
#BBBFFFFFFFFFFIIIIFIFFIIIIIIIIIIIIIIIIIFIIIFFFFFBFF

# original grep commands:
# cat ~/DEEP/tmp.fastq | grep -A 3 '^@.* [^:]*:N:[^:]*:'  | grep -v '^--$'
# first grep gets the read ID and following 3 lines of not-failed reads, second omits the -- lines with which grep separates these matchesOn

my $counter = 0;
my $failed = 0;
while (<>)
{
	$counter++;
	if ($_ =~ /^\@.*[^:]*:N:[^:]*:/)	# not failed
	{
		print $_;	# read ID
		$_ = <>;
		print;		# sequence
		$_ = <>;
		print;		# + line
		$_ = <>;
		print;		# quality
	}
	else
	{
		$failed++;
		$_ = <>;$_ = <>;$_ = <>;	# get rid of the 3 following lines
	}
}
# gunzip or awk always complain about broken pipe, so maybe that helps?
# "close" Closes the file or pipe associated with the file handle.
close (STDIN);
print STDERR "$counter reads, $failed failed\n";
exit;
