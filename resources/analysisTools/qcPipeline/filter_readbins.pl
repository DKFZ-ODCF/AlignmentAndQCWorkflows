#!/usr/bin/env perl

use strict;
use warnings;

if (scalar(@ARGV) < 2) {
	die "1. read bins file (use - for reading from STDIN) - 2. file with chromosomes to keep in first column  - 3. optional: ignore chr prefix and .fa suffix\n";
}

my $binsFile= shift;
my $chromfile = shift;
my $ignore = shift;

open (my $chromFh, $chromfile) or die "could not open chromlength file $chromfile: $!\n";
open (my $binsFh, $binsFile) or die "could not open read bins file $binsFile: $!\n";

my @fields = ();
my %chroms = ();
my $chr;
while (! eof($chromFh)) {
	my $line = readline($chromFh)
		|| die "Couldn't read chromosomes file '$chromfile': $!";
	if ($line =~ /^#/)	{
		next;
	}
	chomp $line;
	@fields = split ("\t", $line);
	# Get rid of possible chr and .fa if wanted
	if ($ignore) {
		($chr = $fields[0]) =~ s/^chr//;
		$chr =~ s/^\.fa$//;
		$chroms{$chr} = 1;
	} else {
		$chroms{$fields[0]} = 1;
	}
}
close $chromFh;

print STDERR (scalar keys %chroms) . " chromosomes to keep (from '$chromfile'): " . join(" ", sort keys %chroms) . "\n";

my $all = 0;
my $kept = 0;
my %foundChromosomesBins = ();
# Read bins: chromosome is in first column
while (! eof($binsFh)) {
	my $line = readline($binsFh)
		|| "Couldn't read bins file '$binsFile': $!";
	++$all;
	@fields = split ("\t", $line);
	my $chr = $fields[0];
	if ($ignore) {
		$chr =~ s/^chr//;
		$chr =~ s/^\.fa$//;
	}
	++$foundChromosomesBins{$chr};
	if (exists $chroms{$chr}) {
		print $line;
		++$kept;
	}
}
close $binsFh;


if ($all == 0) {
	die "Empty read-bins file '$binsFile'!\n";
} elsif ($kept < 1) {
	die "Error in filtering read bins (filter_readbins.pl): No lines were kept from read bins input. Probably the prefixes of the chromosome file and the BAM file are not compatible. Read bins file '$binsFile' contained: " . join(" ", sort keys %foundChromosomesBins) . "\n";
} else {
	print STDERR "From $all lines, kept $kept with selected chromosomes in '$binsFile'.\n";
}


