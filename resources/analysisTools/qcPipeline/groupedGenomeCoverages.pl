#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
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
use IO::String;

my $chromosomeColumn = "interval";
my $coverageAllColumn = "coverage QC bases";
my $coverageAllQuotientColumn = "#QC bases/#total bases";
my $coverageNotNColumn = "genome_w/o_N coverage QC bases";
my $coverageNotNQuotientColumn = "#QC bases/#total not_N bases";

our $USAGE = "
Read in the output table of genomeCoverage.d and sum up the values in specified groups of chromosomes. This way it is possible to get separate values
for human and mouse chromosomes in xenograft samples, as one of the groups usually (in our use cases) has a 'chr' or 'chrMmu' prefix and the other not.
Which groups has which chromosome name format is defined by the calling script.

Run

$0 test

to run the tests, or

$0 coverageTable groupSpec groupSpec*

coverageTable         The output file of coverageQc[.d]. Usually suffixed by '.DepthOfCoverage_Genome.txt'
groupSpec             Of the format groupName=chromosomeName{,chromosomeName}*

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

    subtest "parseGroupSpec" => sub {
        my @specs = ("short=1,2,3", "long=chr1,chr2,chrX");
        is_deeply(parseGroupSpec($specs[0]), { "short" => ["1", "2", "3"] }, "parse short group spec");
        is_deeply(parseGroupSpec($specs[1]), { "long" => ["chr1", "chr2", "chrX"] }, "parse short group spec");
    };

    is_deeply([parseChromosomeIndexList("1,2,3,MT")], ["1", "2", "3", "MT"], "parsing chromosome index list");
    is_deeply([parseChromosomeIndexList("1,2,3,2,MT,3")], ["1", "2", "3", "MT"], "remove duplicates from parsed chromosome index list");

    my $testFile =<<THE_END;
interval	coverage QC bases	#QC bases/#total bases	mapq=0 read1	mapq=0 read2	mapq>0,readlength<minlength read1	mapq>0,readlength<minlength read2	mapq>0,BaseQualityMedian<basequalCutoff read1	mapq>0,BaseQualityMedian<basequalCutoff read2	mapq>0,BaseQualityMedian>=basequalCutoff read1	mapq>0,BaseQualityMedian>=basequalCutoff read2	%incorrect PE orientation	#incorrect proper pair	#duplicates read1 (excluded from coverage analysis)	#duplicates read2 (excluded from coverage analysis)	genome_w/o_N coverage QC bases	#QC bases/#total not_N bases
1	44.01x	10970144211/249250621	1763667	3631605	214324	617373	0	0	36741487	36527982	0.191	51408	4881994	4881927	48.70x	10970144211/225280621
2	48.50x	11795106623/243199373	1856166	3886692	225491	641494	0	0	39426883	39450876	0.252	58129	5378280	5381861	49.52x	11795106623/238204518
3	54.78x	10847876685/198022430	833647	2358591	175894	518065	0	0	36284661	36094829	0.119	38744	4883260	4883856	55.69x	10847876685/194797135
4	48.60x	9290925293/191154276	894880	2355481	165530	489927	0	0	31067768	30912423	0.119	31009	4238801	4238781	49.51x	9290925293/187661676
5	47.43x	8580279977/180915260	1131651	2504345	156913	460994	0	0	28698776	28550559	0.133	27109	3888998	3889078	48.29x	8580279977/177695260
6	48.39x	8280057056/171115067	819130	2149758	152573	447623	0	0	27704730	27558047	0.106	26216	3801299	3801572	49.46x	8280057056/167395066
7	78.72x	12527978168/159138663	1318041	2612053	156293	434179	0	0	41954282	41705948	0.166	42368	5629882	5629805	80.64x	12527978168/155353663
8	64.84x	9490437369/146364022	1105515	2261133	134194	394431	0	0	31745743	31578594	0.102	25677	4362416	4362794	66.42x	9490437369/142888922
9	43.93x	6203563145/141213431	1534576	2458104	107160	309493	0	0	20766912	20650535	0.191	15057	2812124	2812581	51.63x	6203563145/120143431
10	46.59x	6315208372/135534747	861382	1942630	124306	358351	0	0	21150752	21032334	0.110	17495	2802068	2802486	48.09x	6315208372/131314738
11	46.98x	6342946221/135006516	634880	1674593	125235	359130	0	0	21232661	21115957	0.101	16963	2813713	2813720	48.37x	6342946221/131129516
12	62.46x	8361020771/133851895	604571	1667839	128709	367741	0	0	27999889	27836957	0.098	23705	3721624	3721626	64.08x	8361020771/130481393
13	58.51x	6738224563/115169878	429609	1193258	92092	265169	0	0	22541667	22421767	0.072	13316	3074781	3075117	70.49x	6738224563/95589878
14	40.55x	4353239298/107349540	514548	1249206	82783	241113	0	0	14573043	14491923	0.078	8184	1945149	1946128	49.31x	4353239298/88289540
15	51.39x	5269560837/102531392	928992	1568771	79671	219271	0	0	17647140	17542090	0.128	9868	2347185	2346961	64.50x	5269560837/81694766
16	57.22x	5169931526/90354753	1475003	2172086	84678	233594	0	0	17365019	17223066	0.513	13597	2698740	2698643	65.54x	5169931526/78884753
17	45.49x	3693396776/81195210	583268	1262624	84665	233985	0	0	12408789	12318170	0.132	10293	1589959	1590015	47.48x	3693396776/77795210
18	44.29x	3457945360/78077248	337253	946952	69286	203303	0	0	11569402	11509343	0.053	5311	1556587	1556561	46.32x	3457945360/74657229
19	43.84x	2592481529/59128983	470021	1012679	73500	196517	0	0	8745770	8664426	0.150	9534	1086896	1086875	46.45x	2592481529/55808983
20	90.76x	5720067103/63025520	345108	860611	69312	188841	0	0	19172464	19050107	0.059	9446	2486981	2486785	96.13x	5720067103/59505520
21	37.51x	1805194387/48129895	306244	620914	36967	101417	0	0	6049744	6013128	0.036	1666	813787	813593	51.42x	1805194387/35106642
22	29.65x	1521257924/51304566	384196	689847	38411	106513	0	0	5113990	5075024	0.124	2380	646656	646584	43.60x	1521257924/34894545
X	44.76x	6950264113/155270560	1527410	2648956	125687	365157	0	0	23251925	23127592	0.127	18897	3263077	3263393	46.00x	6950264113/151100560
Y	5.42x	321786276/59373566	1516296	1680345	15394	37983	0	0	1083749	1078861	0.319	397	271719	271689	14.00x	321786276/22984529
chrMmu1	20.77x	4060606717/195471971	1493134	2948906	189825	509368	0	0	13647035	13561644	0.251	22448	1790176	1790882	21.16x	4060606717/191909192
chrMmu10	20.63x	2696333452/130694993	1036773	2040257	133083	351271	0	0	9058293	9009380	0.178	11045	1184860	1185360	21.22x	2696333452/127067662
chrMmu11	20.49x	2501183854/122082543	895243	1882308	147742	363221	0	0	8424725	8367836	0.181	11621	1076693	1076986	21.06x	2501183854/118745945
chrMmu12	20.43x	2454469006/120129022	1443695	2000977	411012	297934	0	0	8276098	8179023	0.413	10959	2211005	2086632	20.99x	2454469006/116922420
chrMmu13	20.57x	2476957202/120421639	940306	1851951	117999	311252	0	0	8335451	8276292	0.287	8784	1094334	1095769	21.15x	2476957202/117121193
chrMmu14	19.81x	2473690427/124902244	1380473	2317450	116982	312897	0	0	8327050	8280624	0.316	9199	1378733	1355338	20.37x	2473690427/121442110
chrMmu15	20.59x	2142151231/104043685	664530	1472049	105592	281603	0	0	7196088	7153130	0.152	6734	941532	942746	21.28x	2142151231/100653315
chrMmu16	20.69x	2031532860/98207768	668752	1397766	93618	253325	0	0	6823094	6780110	0.162	5442	899242	899212	21.38x	2031532860/95019758
chrMmu17	20.37x	1935055828/94987271	727142	1402179	100840	255835	0	0	6517847	6472004	0.200	6443	843917	844454	21.10x	1935055828/91707462
chrMmu18	20.60x	1868567701/90702639	581693	1256727	91759	239696	0	0	6280475	6239107	0.141	5170	820047	820097	21.37x	1868567701/87452634
chrMmu19	19.97x	1226996482/61431566	373529	840861	65161	167150	0	0	4127043	4099378	0.085	2489	533083	533065	21.08x	1226996482/58205856
chrMmu2	22.13x	4030174160/182113224	1990803	3425124	189063	500150	0	0	13762659	13403444	1.394	76055	13249105	13305353	22.60x	4030174160/178326651
chrMmu3	20.81x	3331032107/160039680	1041330	2218884	148579	404476	0	0	11190240	11119267	0.203	14058	1478756	1478892	21.30x	3331032107/156398855
chrMmu4	20.25x	3168872286/156508116	1283306	2455729	160605	415051	0	0	10660389	10590887	0.301	15254	1393565	1393846	20.84x	3168872286/152055611
chrMmu5	20.17x	3062415798/151834684	1127493	2285852	165037	419605	0	0	10301993	10236947	0.248	16348	1336018	1335999	20.70x	3062415798/147919674
chrMmu6	20.72x	3102379845/149736546	1036570	2142294	145596	387755	0	0	10438120	10365856	0.222	15022	1903453	1896286	21.20x	3102379845/146336543
chrMmu7	19.23x	2796538027/145441459	1260320	2287415	136959	359482	0	0	9407746	9345395	0.273	12012	1225333	1225744	19.71x	2796538027/141855407
chrMmu8	20.34x	2631568180/129401213	963006	1960280	132048	344194	0	0	8852250	8794189	0.230	10832	1154591	1154284	20.95x	2631568180/125611432
chrMmu9	22.01x	2741828158/124595110	3255474	4409972	131975	344567	0	0	9285506	9141785	1.121	15354	7438656	7420342	22.63x	2741828158/121157018
chrMmuX	19.41x	3319651260/171031299	2037615	3105321	119047	338707	0	0	11122472	11065103	0.295	10486	1525691	1530184	20.31x	3319651260/163487995
chrMmuY	0.27x	25107036/91744698	219888	326874	6967	14475	0	0	91531	91661	0.642	87	13991	14047	0.28x	25107036/88124698
all	36.19x	210676005200/5821198782	46597129	89438249	5628557	14663678	0	0	706423351	702103600	0.206	762611	114488757	114387949	38.26x	210676005200/5506179525
THE_END
;;
    my %table2 = parseTsvFh(IO::String->new($testFile));
    my $aggregationResult = aggregateByChromoseSet(\%table2, {
        "long" =>  [qw(chrMmu2 chrMmu1 chrMmu3 chrMmu4 chrMmu5 chrMmu6 chrMmu7 chrMmuX chrMmu8 chrMmu10 chrMmu11 chrMmu12 chrMmu9 chrMmu13 chrMmu14 chrMmu15 chrMmu16 chrMmu17 chrMmu18 chrMmu19 chrMmu20 chrMmu21 chrMmu22 chrMmuY)],
        "short" => [qw(2 1 3 4 5 6 7 X 8 10 11 12 9 13 14 15 16 17 18 19 20 21 22 Y)],
    });
    is(coverageQuotientString($aggregationResult->{"long"}->{-coverageAll}), "54077111617/2725521370", "long chromosome names");
    is(coverageQuotientString($aggregationResult->{"short"}->{-coverageAll}), "156598893583/3095677412", "short chromosome names");

    is(coverageString($aggregationResult->{"short"}->{-coverageAll}), "50.59x", "coverage short chromosome names");

    is(coverageQuotientString($aggregationResult->{"long"}->{-coverageNotN}), "54077111617/2647521431", "long chromosome names (no N)");
    is(coverageQuotientString($aggregationResult->{"short"}->{-coverageNotN}), "156598893583/2858658094", "short chromosome names (no N)");

    is(coverageString($aggregationResult->{"short"}->{-coverageNotN}), "54.78x", "coverage short chromosome names (no N)");
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

sub coverageString($) {
    my ($hash) = @_;
    return sprintf("%1.2fx", $hash->{-coverageBases}  / $hash->{-totalBases});
}

sub sumCoverageQuotient($$) {
    my ($first, $second) = @_;
    return {
        -coverageBases => $first->{-coverageBases} + $second->{-coverageBases},
        -totalBases    => $first->{-totalBases} + $second->{-totalBases}
    };
}

sub coverageNotNQuotient($) {
    return parseCoverageQuotient(shift(@_), $coverageNotNQuotientColumn);
}


sub coverageAllQuotient($) {
    return parseCoverageQuotient(shift(@_), $coverageAllQuotientColumn);
}

sub assertTableColumns($$) {
    my ($header, $fileName) = @_;
    my %header = map { $_ => 1 } @$header;
    if (! exists($header{$chromosomeColumn})) {
        croak "File does not contain '$chromosomeColumn' column: '$fileName'"
    } elsif (! exists($header{$coverageNotNQuotientColumn})) {
        croak "File does not contain '$coverageNotNQuotientColumn' column: '$fileName'"
    } elsif (!exists($header{$coverageAllQuotientColumn})) {
        croak "File does not contain '$$coverageAllQuotientColumn' column: '$fileName'"
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

# Produce a reduced version of the input table with the aggregated values for each group.
sub printAggregatedValues($) {
    my ($aggregatedValues) = @_;

    sub headerLine() {
        join("\t", (
             $chromosomeColumn,
             $coverageAllColumn,
             $coverageAllQuotientColumn,
             $coverageNotNColumn,
             $coverageNotNQuotientColumn)) . "\n";
    }

    sub groupLine($$) {
        my ($name, $values) = @_;
        join("\t", ($name,
             coverageString($values->{-coverageAll}),
             coverageQuotientString($values->{-coverageAll}),
             coverageString($values->{-coverageNotN}),
             coverageQuotientString($values->{-coverageNotN}))) . "\n";
    }

    print headerLine();
    print join("", map { groupLine($_, $aggregatedValues->{$_}) } sort keys %$aggregatedValues );
}

sub parseGroupSpec($) {
    my ($specString) = @_;
    if ($specString =~ /([^=]+)=(.+)/) {
        return { $1 => [parseChromosomeIndexList($2)] };
    } else {
        croak "Could not parse group name from '$specString'";
    }
}

my ($coverageTableFile,
    @groupSpecs) = @ARGV;

if (defined($coverageTableFile) && lc($coverageTableFile) eq "test") {
    exit runTests();
} elsif (!scalar(@groupSpecs)) {
    print STDERR "Not enough arguments:\n$USAGE";
    exit 1;
}

my %chromosomeSets = map { %{ parseGroupSpec($_) } } @groupSpecs;

my %table = readTsvFile($coverageTableFile);

my $aggregationResult = aggregateByChromoseSet(\%table, \%chromosomeSets);
printAggregatedValues($aggregationResult);


1;