#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

use strict;
use warnings;
use Carp;

our $USAGE = <<ENDE;
    1. read bins file (use - for reading from STDIN)
    2. file with chromosomes to keep in first column
    3. optional: ignore 'chr' prefixes and '.fa' suffixes in both read bins file and chromosomes files;
ENDE


sub maybeRemovePrefixSuffix($$) {
    my ($chr, $doRemove) = @_;
    if ($doRemove) {
        $chr =~ s/^chr//;
        $chr =~ s/^\.fa$//;
    }
    return $chr;
}

if (scalar(@ARGV) < 2) {
    print STDERR $USAGE;
    exit 1;
}

my ($binsfile,
    $chromfile,
    $ignore) = @ARGV;

my $chromFh = IO::File->new($chromfile)
    || croak "Could not open chromosome length file '$chromfile': $!";
my $binsFh = IO::File->new($binsfile)
    || croak "Could not open read bins file '$binsfile': $!";

my %chroms = ();
while (! eof($chromFh)) {
    my $line = readline($chromFh)
        || croak "Couldn't read chromosomes file '$chromfile': $!";
    if ($line !~ /^#/) {
        chomp $line;
        my @fields = split("\t", $line);
        my $chr = maybeRemovePrefixSuffix($fields[0], $ignore);
        $chroms{$chr} = 1;
    }
}
$chromFh->close;

print STDERR "Chromosomes to keep: " . join(" ", sort keys %chroms) . "\n";

my %readBinsChr = ();
my $all = 0;
my $kept = 0;

my $timestamp = localtime();
print STDERR "($timestamp) filter_readbins.pl is waiting for input...\n";

my $line = "";
while (! eof($binsFh)) {
    $line = readline($binsFh)
            || croak "Couldn't read bins file '$binsfile': $!";
    $all++;
    my @fields = split("\t", $line);
    # read bins: chromosome is in the first column
    my $chr = maybeRemovePrefixSuffix($fields[0], $ignore);
    $readBinsChr{$chr} = 1;
    if (exists $chroms{$chr}) {
        print $line;
        $kept++;
    }
}
$binsFh->close;

$timestamp = localtime();
print STDERR "($timestamp) filter_readbins.pl finished reading bins file\n";
print STDERR "Content of last line was ---'$line'---\n";
print STDERR "Chromosomes in read-bins file: " . join(" ", sort keys %readBinsChr) . "\n";
print STDERR "From $all lines, kept $kept with selected chromosomes\n";

if ($kept < 1) {
    croak "Error in filtering read bins (filter_readbins.pl): no lines were extracted from input. Probably the prefixes of the chromosome file and the BAM file are not compatible";
}

