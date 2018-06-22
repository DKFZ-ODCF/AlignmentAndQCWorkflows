#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

use strict;
use warnings;
use Getopt::Std;

my %opts = (p=>1000,n=>1);
getopts("c:d:f:i:l:m:p:r:s:t:w:h", \%opts);

my $coverage = $opts{c};
my $diffchroms = $opts{d};
my $flagstats = $opts{f};
my $isizes = $opts{i};
my $lane = $opts{l};
my $metrics = $opts{m};
my $pid = $opts{p};
my $runID = $opts{r};
my $sample = $opts{s};
my $target = $opts{t};
my $warnfile = $opts{w};


if (defined $opts{h} || ! defined $pid || ! defined $sample  || ! defined $runID || ! defined $lane || ! defined $flagstats || ! defined $diffchroms || ! defined $isizes || ! defined $coverage)
{
	die "USAGE: $0 [options]
	-c FILE	coverage file (coverage/*DepthOfCoverage_Genome.txt)
	-d FILE	diffchrom stats file (structural_variation/*_DiffChroms.png_qcValues.txt)
	-f FILE flagstat file (flagstats/*_flagstats.txt )
	-i FILE	insert size stats file (insertsize_distribution/*_insertsize_plot.png_qcValues.txt)
	-l STRING lane identifier (in case of merged-mardupped BAM: 'genome' for WGS to check for 30x coverage; 'exome' for WES)
	-m FILE	Picard metrics file (alignment/*.dupmark_metrics.txt; optional)
	-p STRING PID
	-r STRING run-ID ('all_merged' for merged-markdupped BAM)
	-s STRING sample-ID
	-t FILE (coverage/*DepthOfCoverage_Target.txt; optional)
	-w FILE output file for warnings, if any occur (warnings also got to STDERR)
	-h help\n";
}

open (FL, $flagstats) or die "could not open flagstat file $flagstats: $!\n";
open (IS, $isizes) or die "could not open insert size stats file $isizes: $!\n";
open (DC, $diffchroms) or die "could not open aberrant reads = diffchroms stats file $diffchroms: $!\n";
open (CO, $coverage) or die "could not open genome coverage file $coverage: $!\n";

my ($totalreads, $mappedreads, $proper, $singletons, $aligned, $dups, $libsize, $aberrant, $iszieSD, $isizemed) = ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA");

# from output of coverageQC
my ($covNoN, $basesNoN, $covX, $covY, $coverageQCstring, $basesall, $ontarget, $targetcov, $qcbases_targetsize) = ("NA", "NA", "NA", "NA", "NA", 0, 0, "NA", "NA");

# collect warnings to print to STDERR and if wished into a warnings file
# if the insert stats file has an entry about 2 peaks; from isize stats
# if the number of singletons is > W%; from flagstats
my $singlimit = 2;
# if the duplication rate is > X%; from Picard metrics
my $duplimit = 20;
# if the fraction of improper pairs is > Y%; from PEaberrations output
my $improplimit = 3;	# some GBMs have that and are really strange
# proper pairs from flagstats correlate with this
my $prpmin = 70;
# if the mapping rate is < Z%; taken from flagstats
my $mapmin = 85;
# coverage of a "all_merged" for genome not 30x
my $mingenomecov = 30;
# and for exomes below 80x
my $minexoncov = 80;
# on target rate is usually > 0.7 but never > 0.75
my $mintargetrate = 0.7;
my $covvalue = 0;	# without the x at the end
my $warnstring = "";

# flagstat => TOTAL_READ_COUNT %TOTAL_READ_MAPPED_BWA %properly_paired %singletons ALIGNED_READ_COUNT


# usual content (samtools flagstat, sambamba_v0.4.6):
# 45762344 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 duplicates
# 45577551 + 0 mapped (99.60%:-nan%)
# 45762344 + 0 paired in sequencing
# 22881172 + 0 read1
# 22881172 + 0 read2
# 45368164 + 0 properly paired (99.14%:-nan%)
# 45489202 + 0 with itself and mate mapped
# 88349 + 0 singletons (0.19%:-nan%)
# 68095 + 0 with mate mapped to a different chr
# 59832 + 0 with mate mapped to a different chr (mapQ>=5)

# but from sambamba v0.5.9 on bwa mem mapped additional lines for secondary and supplementary
# 46022245 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 secondary
# 259901 + 0 supplementary
# 0 + 0 duplicates
# 45954051 + 0 mapped (99.85%:N/A)
# 45762344 + 0 paired in sequencing
# 22881172 + 0 read1
# 22881172 + 0 read2
# 45288038 + 0 properly paired (98.96%:N/A)
# 45631658 + 0 with itself and mate mapped
# 62492 + 0 singletons (0.14%:N/A)
# 194236 + 0 with mate mapped to a different chr
# 127606 + 0 with mate mapped to a different chr (mapQ>=5)

my $line = <FL>;	# in total	
($totalreads) = $line =~ (/^(\d+)/);# total number of reads
$line = <FL>;	# duplicates - //or// secondary, followed by supplementary: in this case, skip over 2 lines
if ($line !~ /duplicates/)
{
	$line = <FL>;
	$line = <FL>;
}
$line = <FL>;	# mapped
($aligned) = $line =~ /^(\d+)/;	# number of mapped reads
($mappedreads) = $line =~ /\((.+)\%:/;	# mapped reads %
if ($mappedreads < $mapmin)
{
	$warnstring .= "MAPPING RATE ${mappedreads}%\n";
}

$line = <FL>;	# paired in sqcing
$line = <FL>;	# read1
$line = <FL>;	# read2
$line = <FL>;	# properly paired
($proper) = $line =~ /\((.+)\%:/;	# properly paired %
if ($proper < $prpmin)
{
	$warnstring .= "PROPERLY PAIRED READS: ${proper}%\n";
}


$line = <FL>;	# with itself and mate mapped 
$line = <FL>;	# singletons
($singletons) = $line =~ /\((.+)\%:/;	# singletons %
if ($singletons > $singlimit)
{
	$warnstring .= "SINGLETON READS ${singletons}%\n";
}

# last 2 lines are not wanted
#with mate mapped to a different chr
#with mate mapped to a different chr (mapQ>=5)
close FL;

my @sline = ();
if (defined $metrics)
{
	open (ME, $metrics) or die "could not open Picard metrics file $metrics: $!\n";
	# Picard metrics => %DUPLICATES  ESTIMATED_LIBRARY_SIZE
	while (<ME>)
	{
		if ($_ =~ /^LIBRARY/)
		{
			$line = <ME>;
			chomp $line;
			# PERCENT_DUPLICATION is in fact a fraction, has to be converted to %
			@sline = split ("\t", $line);
			$dups = sprintf ("%2.2f", ($sline[7]*100));
			if ($dups > $duplimit)
			{
				$warnstring .= "DUPLICATION RATE ${dups}%\n";
			}
			if (scalar @sline >= 8) {
				## Picard may not produce a library size estimate and leave the field free.
				$libsize = $sline[8];    # fake metrics from postmerge QC or sambamba does not have a number
			}
			last;
		}
	}
	close ME;
}

# insert size stats => %sd_PE_insertsize PE_insertsize
my @isstatscontent = <IS>;
close IS;
$isizemed = $isstatscontent[0];
chomp $isizemed;
$iszieSD = $isstatscontent[1];
chomp $iszieSD;
# 3rd line is standarddeviation
# if there is a 4th line, it's the bimodality warning
if (@isstatscontent > 3)
{
	$warnstring .= $isstatscontent[3];	# INSERT SIZE DISTRIBUTION POSSIBLY BIMODAL: (first) 2 peaks at $mode and $mode1\n"
}

# diffchroms stats  => %PE_reads_on_diff_chromosomes
$aberrant = <DC>;
chomp $aberrant;
if ($aberrant ne "NA" && $aberrant > $improplimit)
{
	$warnstring .= "PAIRED END READS MAPPING ON DIFFERENT CHROMOSOMES: ${aberrant}%\n";
}

# coverageQC is the most complicated:
#"coverage" is the genome-wide coverage for genome
#"coverage_base" is the genome-wide coverage but with different interpretation mode for target
#"coverage_target" is the target coverage file, only used for target

while (<CO>)
{
	chomp;
	@sline = split ("\t", $_);
	# X and Y are of special interest - % coverage without N
	if ($sline[0] =~ /X/)
	{
	# the second last element
		$covX = $sline[-2];
	}
	if ($sline[0] =~ /Y/)
	{
		$covY = $sline[-2];
	}
	if ($sline[0] =~ /^all/)
	# need almost the whole line except the outcrossed things:
	#--interval--        coverage QC bases       #QC bases/#total bases  mapq=0 read1    mapq=0 read2    mapq>0,readlength<minlength read1       mapq>0,readlength<minlength read2       mapq>0,BaseQualityMedian<basequalCutoff read1 mapq>0,BaseQualityMedian<basequalCutoff read2    mapq>0,BaseQualityMedian>=basequalCutoff read1  mapq>0,BaseQualityMedian>=basequalCutoff read2  %incorrect PE orientation       #incorrect proper pair  #duplicates read1 (excluded from coverage analysis)    #duplicates read2 (excluded from coverage analysis)     --genome_w/o_N coverage QC bases--  --#QC bases/#total not_N bases--
	{
		shift @sline;
		$basesNoN = pop @sline;
		($basesall) = $basesNoN	=~ /^(\d+)/;	# QC bases
		#print STDERR "all bases from genome QC summary: $basesall\n";
		$covNoN = pop @sline;
		$coverageQCstring = join ("\t", @sline);
		if ($lane eq "genome" && $runID eq "all_merged")
		{
			# this would always raise a warning for the genome-wide coverage file
			# from an exome capture, because in there it is genome-wide coverage
			($covvalue = $covNoN) =~ s/x$//;
			if ($covvalue < $mingenomecov)
			{
				$warnstring .= "GENOME COVERAGE $covNoN\n";
			}
		}		
	}
}
close CO;


# the average target coverage is in "all", the number of QC bases here has to be
# divided by that of the whole genome to get the on target ratio
if (defined $target)
{
	open (T, $target) or warn "could not open target coverage file $target: $!\n";
	while (<T>)
	{
		chomp;
		@sline = split ("\t", $_);
		if ($sline[0] =~ /^all/)
		# need almost the whole line except the outcrossed things:
		#--interval--        coverage QC bases       #QC bases/#total bases  mapq=0 read1    mapq=0 read2    mapq>0,readlength<minlength read1       mapq>0,readlength<minlength read2       mapq>0,BaseQualityMedian<basequalCutoff read1 mapq>0,BaseQualityMedian<basequalCutoff read2    mapq>0,BaseQualityMedian>=basequalCutoff read1  mapq>0,BaseQualityMedian>=basequalCutoff read2  %incorrect PE orientation       #incorrect proper pair  #duplicates read1 (excluded from coverage analysis)    #duplicates read2 (excluded from coverage analysis)     --genome_w/o_N coverage QC bases--  --#QC bases/#total not_N bases--
		{
			shift @sline;
			$qcbases_targetsize = pop @sline;
			($ontarget) = $qcbases_targetsize =~ /^(\d+)/;	# QC bases
			#print STDERR "bases on target: $ontarget\n";
			$targetcov = pop @sline;
			$coverageQCstring = join ("\t", @sline);
			($covvalue = $targetcov) =~ s/x$//;
			# print STDERR "target coverage: $targetcov\n";
			if ($covvalue < $minexoncov)
			{
				$warnstring .= "AVERAGE EXOME COVERAGE $targetcov\n";
			}
		}
	}
}

# finally write the output:
# hardcoded header ...

print "PID\tPID\tSAMPLE_TYPE\tRUN_ID\tLANE\tTOTAL_READ_COUNT (flagstat)\t%TOTAL_READ_MAPPED_BWA (flagstat)\t%properly_paired (flagstat)\t%singletons (flagstat)\tALIGNED_READ_COUNT\t%DUPLICATES (Picard metrics file)\tESTIMATED_LIBRARY_SIZE (Picard metrics file)\t%PE_reads_on_diff_chromosomes (mapq>0)\t%sd_PE_insertsize (mapq>0)\tPE_insertsize (mapq>0)\t";
if (defined $target)
{
	print "coverage QC bases whole genome w/o N\tQC bases/ total bases whole genome w/o N\tOn Target ratio\tcoverage QC bases On Target\tQC bases/ total bases On Target\t";
}
else
{
	print "coverage QC bases w/o N\tQC bases/ total bases w/o N\t";
}
print "coverage QC bases\tQC bases/ total bases\tmapq=0 read1\tmapq=0 read2\tmapq>0,readlength<minlength read1\tmapq>0,readlength<minlength read2\tmapq>0,BaseQualityMedian<basequalCutoff read1\tmapq>0,BaseQualityMedian<basequalCutoff read2\tmapq>0,BaseQualityMedian>=basequalCutoff read1\tmapq>0,BaseQualityMedian>=basequalCutoff read2\t%incorrect PE orientation\t#incorrect proper pair\t#duplicates read1 (excluded from coverage analysis\t#duplicates read2 (excluded from coverage analysis)\tChrX coverage QC bases\tChrY coverage QC bases\n";

print "$pid\t$pid\t$sample\t$runID\t$lane\t$totalreads\t$mappedreads\t$proper\t$singletons\t$aligned\t$dups\t$libsize\t$aberrant\t$iszieSD\t$isizemed\t$covNoN\t$basesNoN\t";

if (defined $target)
{
	# on target ratio
	my $targetratio = sprintf ("%2.2f", ($ontarget/$basesall));
	if ($targetratio < $mintargetrate)
	{
		$warnstring .= "ON TARGET RATE $targetratio\n";
	}
	# also print on target average coverage and QC bases divided by size of target region
	print "$targetratio\t$targetcov\t$qcbases_targetsize\t";
}
print $coverageQCstring, "\t$covX\t$covY\n";

if ($warnstring ne "")
{
	print STDERR $warnstring;
	if (defined $warnfile)
	{
		open  (WF, ">$warnfile") or warn "could not open output file for writing the warnings to, $warnfile: $!\n";
		print WF $warnstring;
		close WF;
	}
}
exit;
