#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

use strict;
use warnings;

if (@ARGV < 1)
{
	die "Usage: $0 output of coverageBed -d\n";
}

my $file = shift;

open (FH, $file) or die "Cannot open $file: $!\n";

# coverageBed -d gives depth of coverage for each position in a range
# Each position and depth follow the complete B feature
# the coverage is the average: sum of all reads / length of region

# e.g.
# chr1    10244   10444   1       1	# first
# chr1    10244   10444   2       1
# chr1    10244   10444   3       1
# chr1    10244   10444   4       1
# ...
# chr1    10244   10444   199     0
# chr1    10244   10444   200     0
# chr1    28500   28565   1       0	# next
# chr1    28500   28565   2       0
# chr1    28500   28565   3       0
# chr1    28500   28565   4       0
# ...
# chr1    28500   28565   64      1
# chr1    28500   28565   65      6	# end


# sum reads over positions; length is the last entry before a new one; calculate average coverage
# the last 2 columns will be changed in the output; before that, 3+ columns have to be kept
# original coverageBed entry appends for complete interval:
# number of reads, length, covered length, coverage ratio.
# read numbers: if coverage increases by X at a pos (compared to previous pos.), reads+=X

my $sum = 0;	# of reads overlapping positions
my $length = 0;	# of the interval
my $average = 0;	# coverage
my $basescov = 0;	# pos. with non-zero coverage
my $ratio = 0;	# covered bases/length

# start with first line to initialize
my $line1 = <FH>;
chomp $line1;
my @help = split ("\t", $line1);

my $cov = pop @help;
my $pos = pop @help;

$sum+=$cov;
if ($cov)	# > 0
{
	$basescov++;
}
$length = $pos;	# will store the last entry!
my $interval = join ("\t", @help);	# keep current interval
my $reads = $cov;	# reads in interval
my $oldnum = $cov;	# keep read number of old position. if the difference between this
			# and the current pos.'s covearge is > 0, increase read counter
while (<FH>)
{
	chomp;
	@help = split ("\t", $_);
	$cov = pop @help;
	$pos = pop @help;
	if ($pos == 1)	# a new entry starts, print out the old entry
	{
		$ratio = $basescov/$length;
		$average = $sum/$length;
		print "$interval\t$reads\t$basescov\t$length\t$ratio\t$average\n";
		#print STDERR "$sum / $length\n";
		# reset the evil global variables
		$sum = 0;
		$basescov = 0;
		$oldnum = $cov;
		$reads = $cov;
	}
	if ($cov > $oldnum)	# 1 or > 1 reads have started at that pos.
	{
		$reads+= ($cov-$oldnum);
	}
	if ($cov)	# > 0
	{
		$basescov++;
	}
	$sum+=$cov;
	$length = $pos;
	$interval = join ("\t", @help);	# keep old data
	$oldnum = $cov;
}
close FH;

# print out the last interval
$ratio = $basescov/$length;
$average = $sum/$length;
print "$interval\t$reads\t$basescov\t$length\t$ratio\t$average\n";
#print STDERR "$sum / $length\n";

exit;
