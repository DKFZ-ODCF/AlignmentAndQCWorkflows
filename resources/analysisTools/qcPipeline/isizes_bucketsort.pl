#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

# bucketsort insert sizes, calculate average, median, standard deviation and "SDpercent" as standard deviation/median*100
# print out the insert size frequencies on stdout for R:
# insertsize_bin\tnumber_reads
# optional, print out median and "SDpercent" to a file (mimicking the text output of insertSizeDistribution.r)

use strict;
use warnings;

if (@ARGV < 1)
{
	die "Anwendung: $0 1.file with insert sizes - optional 2.: file to write the median and standard deviation/median*100 out to (..._plot.png_qcValues.txt) - optional 3.:maximal insert size (default 1000)\n";
}

my $file = shift;

open (IN, $file) or die "Cannot open $file: $!\n";
my $outfile = shift;
if (defined $outfile)
{
	open (OUT, ">$outfile") or die "could not open $outfile for writing: $!\n";
}

my $maxisize = shift;
if (defined $maxisize)
{
	if ($maxisize !~ /^\d+$/)
	{
		print STDERR "maximal insert size $maxisize is not a number, setting to default 1000\n";
		$maxisize = 1000;
	}
}
else
{
	$maxisize = 1000;
}

my @isizes = ();
my $entry = 0;	# absolute value of insert size
my $min = 10000000000000000000;
my $max = 0;
my $ctr = 0;	# number of entries
my $sum = 0;	# sum of insert sizes
while (<IN>)
{
	chomp;
	$entry = abs($_);
	if ($entry < $min)
	{
		$min = $entry;
	}
	if ($entry > $maxisize)	# bwa mem seems to give flag "proper pair" to each read => insert size can go over 200 million bp, makes no sense at all
	{
		next;
	}
	if ($entry > $max)
	{
		$max = $entry;
	}
	$isizes[$entry]++;
	$ctr++;
	$sum+=$entry;
}
close IN;
my $average = $sum/$ctr;
print STDERR "$ctr values. minimum: $min; maximum: $max; average: $average\n";

# the R script calculates:
#MEDIAN=median(abs(sample[,1]))
#SD=round(sd(abs(sample[,1])))
#SDpercent=round(sd(abs(sample[,1]))/median(abs(sample[,1]))*100)

# The median is in a sorted list the value in the middle. The list is actually sorted but contains the same values several times.
# I know how many entries are there. So with some smart calculation, I can "jump" to it. In case of a even number, the median is
# the average of the 2 middle values. But for isizes that will be the same value as middle-1 anyways -
# here for simplicity, just take the value when read number is > 1/2 of all
# Also get the standard deviation by the sqrt of the deviations from the average by going through the list efficiently

my $median = 0;
my $middle = int($ctr/2);
my $howmany = 0;	# number of entries for each isize bin
my $ctr2 = 0;
my $docollect = 1;	# a flag for the median determination
my $devsqsum = 0;		# the sum of the squared deviations from the average
# the modes = peaks; should be only 1 (ideally, the median), but might be 2 in case of adapter dimers, then it's bad!
# a real peak must have several 1000 reads to really stand out
my $mode1 = 0;
my $mode = 0;
my $modnr = 0;
my $oldnr = 0;	# number of reads in previous bin
my $newmode = 1;	# a flag for setting the first mode

# print on stdout for R:
# insertsize_bin\tnumber_reads
for (my $i = $min; $i <= $max; $i++)
{
	$howmany = $isizes[$i] || 0;
	print "$i\t$howmany\n";
	# find highest peak
	if ($howmany > 1000 && $howmany > $modnr)
	{
		$modnr = $howmany;
		$mode1 = $i;
	}
	# drop after a peak => keep this peak as the first mode
	# this is not optimal because after a small dip, there could be the real first peak
	if ($howmany > 1000 && $howmany < $oldnr && $newmode > 0)
	{
		$mode = $mode1;
		$mode1 = 0;
		$modnr = 0;
		$newmode = 0;
	}
	$ctr2+=$howmany;
	$oldnr = $howmany;
	if ($ctr2 >= $middle && $docollect > 0)	# we have hit the middle, fix the median!
	{
		$median = $i;
		print STDERR "hit middle ($middle) at $ctr2 entries: median is $i\n";
		$docollect = 0;
	}
	# quared deviations from the average
	$devsqsum += $howmany * (($i - $average)**2);
}

my $stddev = sqrt($devsqsum/($ctr-1));
my $sdrd=int($stddev+0.5);
print STDERR "stddev: $stddev; rounded: $sdrd\n";
my $sdperc = $stddev/$median*100;
my $rounded = int($sdperc+0.5);
print STDERR "median: $median; SDpercent: $sdperc; rounded: $rounded\n";

# mode(s)
if ($mode1 > 0)
{
	print STDERR "2 peaks: $mode and $mode1";
	# if the difference between the 2 peaks is small, it was just a random fluctuation,
	# i.e. a small dip before the real peak
	if (($mode1 - $mode) > 10)
	{
		print STDERR " - caution, bimodal distribution, might be adapter dimers!\n";
	}
	else
	{
		print STDERR " - small difference does not hint at critical binomial distribution.\n";
	}
}
else
{
	print STDERR "1 peak: $mode\n";
}

if (defined $outfile)
{
	print OUT "$median\n$rounded\n$sdrd\n";
}

exit;
