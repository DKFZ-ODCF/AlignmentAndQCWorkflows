#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

use strict;
use warnings FATAL => 'all';
use IO::File;
use Carp;
use Scalar::Util qw(looks_like_number);
use Pod::Usage;
use Test::Simple;
use Test::More;
use Data::Dumper;

my $chromosomeColumn = "interval";
my $coverageAllColumn = "#QC bases/#total bases";
my $coverageNotNColumn = "#QC bases/#total not_N bases";

our $USAGE = "
Read in the output table of genomeCoverage.d and sum up the values in specific columns for the rows (chromosomes) matching or not matching
the '\$prefix\$chr\$suffix' pattern with '\$chr' from a list of chromosomes. This way it is possible to get separate values for human and mouse
chromosomes in xenograft samples. Whether the matching or non-matching chromosome-set are mouse or human depends on the assembly. The only
assumption is that in xenograft, of host and xenograft genomes one has chromosome names matching, the other has names not matching the pattern.

Run

$0 test

to run the tests, or

$0 coverageTable prefix chromosomeIndexList [suffix]

coverageTable         The output file of coverageQc[.d]. Usually suffixed by '.DepthOfCoverage_Genome.txt'
prefix                e.g. '' or 'chr'
chromosomeIndexList   e.g. '1,2,3', chromosome identifiers without any prefix
suffix                optional suffix. Defaults to ''

";


sub runTests() {
    require "Test/More.pm";
    Test::More->import("no_plan");

    my @header = (
        "interval", "coverage QC bases", "#QC bases/#total bases",
        "mapq=0 read1", "mapq=0 read2",
        "mapq>0,readlength<minlength read1", "mapq>0,readlength<minlength read2",
        "mapq>0,BaseQualityMedian<basequalCutoff read1", "mapq>0,BaseQualityMedian<basequalCutoff read2",
        "mapq>0,BaseQualityMedian>=basequalCutoff read1", "mapq>0,BaseQualityMedian>=basequalCutoff read2",
        "%incorrect PE orientation", "#incorrect proper pair",
        "#duplicates read1 (excluded from coverage analysis)", "#duplicates read2 (excluded from coverage analysis)",
        "genome_w/o_N coverage QC bases",
        "#QC bases/#total not_N bases");

    my @line1 = qw(1 0.00x 20000/40000 9       9       1       1       0       0       79      80      0.000   0       0       0       0.00x   30000/60000);
    my %row = (
        "interval"                                            => "1",
        "coverage QC bases"                                   => "0.00x",
        "#QC bases/#total bases"                              => "20000/40000",
        "mapq=0 read1"                                        => "9",
        "mapq=0 read2"                                        => "9",
        "mapq>0,readlength<minlength read1"                   => "1",
        "mapq>0,readlength<minlength read2"                   => "1",
        "mapq>0,BaseQualityMedian<basequalCutoff read1"       => "0",
        "mapq>0,BaseQualityMedian<basequalCutoff read2"       => "0",
        "mapq>0,BaseQualityMedian>=basequalCutoff read1"      => "79",
        "mapq>0,BaseQualityMedian>=basequalCutoff read2"      => "80",
        "%incorrect PE orientation"                           => "0.000",
        "#incorrect proper pair"                              => "0",
        "#duplicates read1 (excluded from coverage analysis)" => "0",
        "#duplicates read2 (excluded from coverage analysis)" => "0",
        "genome_w/o_N coverage QC bases"                      => "0.00x",
        "#QC bases/#total not_N bases"                        => "30000/60000"
    );
    is_deeply(getRow(\@header, \@line1), \%row, "header/row zipping");

    my @line2 = qw(chr1 0.00x 30000/30000 9 9 1 1 0 0 79 80 0.000 0 0 0 0.00x 50000/30000);
    my @line3 = qw(chr2 0.00x 10000/60000 9 9 1 1 0 0 79 80 0.000 0 0 0 0.00x 30000/70000);
    my @line4 = qw(1 0.00x 30000/70000 9 9 1 1 0 0 79 80 0.000 0 0 0 0.00x 30000/70000);
    my %table = (
        "1"    => getRow(\@header, \@line1),
        "chr1" => getRow(\@header, \@line2),
        "chr2" => getRow(\@header, \@line3),
        "2"    => getRow(\@header, \@line4)
    );

    subtest "correctly count bases for all-base coverage and non-N-base coverage" => sub {
        my %chromosomeSets = ("theShorts" => [ "1", "2" ], "theLongs" => [ "chr1", "chr2" ]);
        my $result = aggregateByChromoseSet(\%table, \%chromosomeSets);
        is_deeply($result->{"theShorts"},
            {   -coverageAll  => { -coverageBases => 50000, -totalBases => 110000 },
                -coverageNotN => { -coverageBases => 60000, -totalBases => 130000 }
            },
            "short sums");
        is_deeply($result->{"theLongs"},
            {   -coverageAll  => { -coverageBases => 40000, -totalBases => 90000 },
                -coverageNotN => { -coverageBases => 80000, -totalBases => 100000 }
            },
            "long sums");
    };

    subtest "special cases" => sub {
        # Count chromosomes in two sets to both sets.
        # Produce 0-results for sets that are not found.
        my %chromosomeSetsWithEmptyAndOverlappingSets = ("theShorts" => [ "1", "2" ], "theLongs" => [ "chr1", "chr2", "1", "blub" ], "theNonExistent" => [ "bla", "blub" ]);
        my $resultsWithEmptySet = aggregateByChromoseSet(\%table, \%chromosomeSetsWithEmptyAndOverlappingSets);

        subtest "count bases in two sets to each" => sub {
            is_deeply($resultsWithEmptySet->{"theShorts"},
                {   -coverageAll  => { -coverageBases => 50000, -totalBases => 110000 },
                    -coverageNotN => { -coverageBases => 60000, -totalBases => 130000 }
                },
                "short sums");
            is_deeply($resultsWithEmptySet->{"theLongs"},
                {   -coverageAll  => { -coverageBases => 60000, -totalBases => 130000 },
                    -coverageNotN => { -coverageBases => 110000, -totalBases => 160000 }
                },
                "long sums");
        };

        is_deeply($resultsWithEmptySet->{"theNonExistent"},
            {   -coverageAll  => { -coverageBases => 0, -totalBases => 0 },
                -coverageNotN => { -coverageBases => 0, -totalBases => 0 }
            },
            "count non-existent sets as 0");
    };

    is_deeply([parseChromosomeIndexList("1,2,3,MT")], ["1", "2", "3", "MT"], "parsing chromosome index list");
    is_deeply([parseChromosomeIndexList("1,2,3,2,MT,3")], ["1", "2", "3", "MT"], "remove duplicates from parsed chromosome index list");

    return 0;
}

sub unique (@) {
    my @list = @_;
    my %index = map { $_ => 1 } @list;
    return grep { exists $index{$_}; $index{$_}-- > 0; } @list;
}

sub parseChromosomeIndexList($) {
    my ($listString) = @_;
    my @result = split(",", $listString);
    if (!scalar(@result)) {
        carp "Chromosome index list is empty";
    }
    return unique @result;
}


sub getRow($$) {
    my ($headerRef, $lineRef) = @_;
    my %row = map { $headerRef->[$_] => $lineRef->[$_] } @{[0 .. (scalar(@$headerRef) - 1)]};
    return \%row;
}

sub parseLine($) {
    my ($line) = @_;
    chomp $line;
    return [split("\t", $line)];
}

sub getValue($$) {
    my ($row, $columnName) = @_;
    my $value = $row->{$columnName};
    if (!defined $value) {
        croak "Undefined '$columnName' value in row: " . join(", ", %$row);
    }
    return $value;
}

sub chromosome($) {
    return getValue(shift(@_), $chromosomeColumn);
}


sub newQuotient() {
    return { -coverageBases => 0, -totalBases => 0 };
}

sub parseCoverageQuotient($$) {
    my ($row, $field) = @_;
    my $value = getValue($row, $field);
    my ($coveredBases, $totalBases) = $value =~ /(\d+)\/(\d+)/;
    if (!defined $coveredBases) {
        croak "Could not parse coverage from '$field': " . join(", ", @$row);
    }
    return {
        -coverageBases => $coveredBases,
        -totalBases    => $totalBases
    };
}

sub coverageQuotientString($) {
    my ($hash) = @_;
    return $hash->{-coverageBases} . "/" . $hash->{-totalBases};
}

sub sumCoverageQuotient($$) {
    my ($first, $second) = @_;
    return {
        -coverageBases => $first->{-coverageBases} + $second->{-coverageBases},
        -totalBases    => $first->{-totalBases} + $second->{-totalBases}
    };
}

sub coverageNotNQuotient($) {
    return parseCoverageQuotient(shift(@_), $coverageNotNColumn);
}


sub coverageAllQuotient($) {
    return parseCoverageQuotient(shift(@_), $coverageAllColumn);
}

sub assertTableColumns($$) {
    my ($header, $fileName) = @_;
    my %header = map { $_ => 1 } @$header;
    if (! exists($header{$chromosomeColumn})) {
        croak "File does not contain '$chromosomeColumn' column: '$fileName'"
    } elsif (! exists($header{$coverageNotNColumn})) {
        croak "File does not contain '$coverageNotNColumn' column: '$fileName'"
    } elsif (!exists($header{$coverageAllColumn})) {
        croak "File does not contain '$$coverageAllColumn' column: '$fileName'"
    }
}

=head1

Create a hash representing a table. The table keys are the chromosome names, the values are hashes with the column names as keys.

=cut
sub parseTsvFh($) {
    my ($fh) = @_;
    my $header = parseLine(scalar <$fh>);
    assertTableColumns($header, "filehandle");
    my %table = ();
    while (my $line = <$fh>) {
        my $row = getRow($header, parseLine($line));
        $table{chromosome($row)} = $row;
    }
    return %table;
}

sub readTsvFile($) {
    my ($file) = @_;
    my $fh = IO::File->new($file, "r")
        || croak "Could not open '$file'";
    my %table = parseTsvFh($fh);
    close($fh);
    return %table;
}

sub getLongChromosomeName($$$) {
    my ($prefix, $name, $suffix) = @_;
    return "$prefix$name$suffix";
}

sub getLongChromosomeNames($$$) {
    my ($prefix, $names, $suffix) = @_;
    return map { getLongChromosomeName($prefix, $_, $suffix) } @$names;
}

=head1

The aggregator data structure is just a hash-reference. Could have made this a class, but that would have been a bit more complicated
and would require more boilerplate code in Perl, which I did not want to introduce here for a two-method data type. Note that above
the sumCoverage() function deals with a similar structure (another Monoid-like type with different field names).

=cut
sub newAggregator(;$$) {
    my ($coverageNotNQuotient, $coverageAllQuotient) = @_;
    defined $coverageNotNQuotient or $coverageNotNQuotient = newQuotient();
    defined $coverageAllQuotient or $coverageAllQuotient = newQuotient();
    ref($coverageNotNQuotient) eq "HASH" || confess "Expect HASH. Got: '$coverageNotNQuotient'";
    ref($coverageAllQuotient) eq "HASH" || confess "Expect HASH. Got: '$coverageAllQuotient'";
    return {
        -coverageNotN => $coverageNotNQuotient,
        -coverageAll => $coverageAllQuotient
    };
}

sub aggregate($@) {
    my ($aggregator, %other) = @_;
    return newAggregator(
        sumCoverageQuotient($aggregator->{-coverageNotN}, $other{-coverageNotN}),
        sumCoverageQuotient($aggregator->{-coverageAll}, $other{-coverageAll}));
}

sub aggregateRow($$) {
    my ($aggregator, $row) = @_;
    return aggregate($aggregator,
        -coverageAll => coverageAllQuotient($row),
        -coverageNotN => coverageNotNQuotient($row));
}

sub chromosomeSetToIndexHash($) {
    my $chromSet = shift @_;
    return map { $_ => 1; } @$chromSet;
}

=head1

Take a table (see parseTsvFh()) and a hash mapping from set names to lists of chromosome names.

Return a hash with the set names as keys and the aggregated counts for each set according to aggregate($$$) as values.

=cut
sub aggregateByChromoseSet($$) {
    my ($table, $chromosomeSets) = @_;
    my %result = ();
    foreach my $chrSetName (keys %$chromosomeSets) {
        $result{$chrSetName} = newAggregator();
        my %chromosomeIndex = chromosomeSetToIndexHash($chromosomeSets->{$chrSetName});
        foreach my $chr (keys %$table) {
            if (exists $chromosomeIndex{$chr}) {
                $result{$chrSetName} = aggregateRow($result{$chrSetName}, $table->{$chr});
            }
        }
    }
    return \%result;
}

# Produce a reduced version of the input table with the aggregated values for matching and non-matching rows.
sub printAggregatedValues($$) {
    my ($long, $short) = @_;
    my $delimiter = "\t";
    print join($delimiter, ($chromosomeColumn,                       $coverageAllColumn,                             $coverageNotNColumn))  . "\n";
    print join($delimiter, ("long",  coverageQuotientString( $long->{-coverageAll}), coverageQuotientString( $long->{-coverageNotN}))) . "\n";
    print join($delimiter, ("short", coverageQuotientString($short->{-coverageAll}), coverageQuotientString($short->{-coverageNotN}))) . "\n";
}





my ($coverageTableFile,
    $prefix,
    $shortChromosomeNames,
    $suffix) = @ARGV;

if (defined($coverageTableFile) && lc($coverageTableFile) eq "test") {
    exit runTests();
} elsif (!defined $shortChromosomeNames) {
    print STDERR "Not enough arguments:\n$USAGE";
    exit 1;
}
if (!defined $suffix) {
    $suffix = "";
}

my %chromosomeSets = (
  "long" =>  [getLongChromosomeNames($prefix, [parseChromosomeIndexList($shortChromosomeNames)], $suffix)],
  "short" => [parseChromosomeIndexList($shortChromosomeNames)]
);

my %table = readTsvFile($coverageTableFile);

my $aggregationResult = aggregateByChromoseSet(\%table, \%chromosomeSets);
printAggregatedValues($aggregationResult->{"long"}, $aggregationResult->{"short"});


1;