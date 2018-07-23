#!/usr/bin/perl
# Read in a number of statistics files and compile some of the content into a JSON file.
# The JSON file is required by OTP.
# Philip R. Kensche <p.kensche@dkfz.de>


use strict;
use warnings;
use Carp;
use IO::File;
use Scalar::Util qw(looks_like_number);
use Pod::Usage;
use Test::Simple;
use Test::More;
use constant { INDENT_INCREMENT => 2 };
use Data::Dumper;


sub runTests () {
    require "Test/More.pm";
    Test::More->import('no_plan');

    my @testGenomeCoverage = (
		"interval	coverage QC bases	#QC bases/#total bases	mapq=0 read1	mapq=0 read2	mapq>0,readlength<minlength read1	mapq>0,readlength<minlength read2	mapq>0,BaseQualityMedian<basequalCutoff read1	mapq>0,BaseQualityMedian<basequalCutoff read2	mapq>0,BaseQualityMedian>=basequalCutoff read1	mapq>0,BaseQualityMedian>=basequalCutoff read2	%incorrect PE orientation	#incorrect proper pair	#duplicates read1 (excluded from coverage analysis)	#duplicates read2 (excluded from coverage analysis)	genome_w/o_N coverage QC bases	#QC bases/#total not_N bases",
		"1	0.01x	2885357/249250621	1373	1509	168	177	0	0	14400	14305	0.301	8	25	27	0.01x	2885357/225280621",
		"Y	0.00x	114475/59373566	809	820	2	7	0	0	571	571	0.273	0	0	0	0.02x	114475/22984529",
		"all	0.01x	35857865/3095677412	16470	20185	2003	2149	0	0	178934	178019	0.293	119	384	397	46.03x	35857865/2858658094"
	);


    my @readLines = readFile($0);
    is($readLines[0], "#!/usr/bin/perl", "readLines");

    subtest 'parseGenomeCoverage' => sub {
		my $docRes = parseGenomeCoverage(@testGenomeCoverage);
		is(ref($docRes), "HASH", "return hash");
		ok(eq_set([keys %$docRes], [qw(1 all Y)]), "contigs");
		ok(exists($docRes->{"all"}->{"genome_w/o_N coverage QC bases"}), "genome coverage w/o N");
		is($docRes->{"all"}->{"genome_w/o_N coverage QC bases"}, "46.03x", "genome coverage w/o N");
    };

    subtest 'parseInsertSize' => sub {
	eval {
	    parseInsertSize(qw(1 2));
	};
	like($@, qr/^Parse error/, "parse error");
	my $stats = parseInsertSize(qw(399 23 93));
	is($stats->{-median}, 399, "median");
	is($stats->{-stddev}, 93, "standard deviation");
	is($stats->{-varcoef}, 23, "coefficient of variation");
    };


    my @testFlagstats = (
		"421309 + 0 in total (QC-passed reads + QC-failed reads)",
		"805 + 0 duplicates",
		"420369 + 0 mapped (99.78%:-nan%)",
		"421309 + 0 paired in sequencing",
		"209146 + 0 read1",
		"212163 + 0 read2",
		"384766 + 0 properly paired (91.33%:-nan%)",
		"419289 + 0 with itself and mate mapped",
		"1080 + 0 singletons (0.26%:-nan%)",
		"33635 + 0 with mate mapped to a different chr",
		"6161 + 0 with mate mapped to a different chr (mapQ>=5)"
	);

    subtest 'parseFlagstats' => sub {
	eval {
	    parseFlagstats("421XXXX309 + 0 in total (QC-passed reads + QC-failed reads)",
			   "805 + 0 duplicates",
			   "420369 + 0 mapped (99.78%:-nan%)",
			   "421309 + 0 paired in sequencing",
			   "209146 + 0 read1",
			   "212163 + 0 read2",
			   "384766 + 0 properly paired (91.33%:-nan%)",
			   "419289 + 0 with itself and mate mapped",
			   "1080 + 0 singletons (0.26%:-nan%)",
			   "33635 + 0 with mate mapped to a different chr",
			   "6161 + 0 with mate mapped to a different chr (mapQ>=5)"
		);
	};
	like($@, qr/Parse error/, "parse error");
	my $stats = parseFlagstats(@testFlagstats);
	is($stats->{-total}, 421309, "total");
	is($stats->{-duplicates}, 805, "duplicates");
	is($stats->{-mapped}, 420369, "mapped");
	is($stats->{-percentage_mapped}, 99.78, "percentage mapped");
	is($stats->{-paired_in_sequencing}, 421309, "paired in sequencing");
	is($stats->{-read1}, 209146, "read1");
	is($stats->{-read2}, 212163, "read2");
	is($stats->{-properly_paired}, 384766, "properly paired");
	is($stats->{-percentage_properly_paired}, 91.33, "percentage properly paired");
	is($stats->{-with_itself_and_mate_mapped}, 419289, "with itself and mate mapped");
	is($stats->{-singletons}, 1080, "singletons");
	is($stats->{-percentage_singletons}, 0.26, "percentage singletons");
	is($stats->{-with_mate_mapped_to_different_chr}, 33635, "with mate mapped to a different chr");
	is($stats->{-with_mate_mapped_to_different_chr_mapQge5}, 6161, "with mate mapped to a different chr (mapQ>=5)");
    };


    subtest 'parseDiffChrom' => sub {
		my $stats = parseDiffChrom(1.55);
		is($stats->{-percentage_mates_on_different_chr}, 1.55, "diffChrom");
    };

    subtest 'toJSON' => sub {
		my $res1 = toJSON({});
		$res1 =~ s/\s//g;
		like($res1, qr/^{}*$/, "{}");

		my $res2 = toJSON({ "a" => 1 });
		$res2 =~ s/\s//g;
		like($res2, qr/^{"a":1}$/, "{ \"a\": 1 }");

		my $res3 = toJSON([ "a", 1 ]);
		$res3 =~ s/\s//g;
		like($res3, qr/^\["a",1\]$/, "[ \"a\", 1 ]");

		my $res4 = toJSON({ "a" => [ "b", "c" ] });
		$res4 =~ s/\s//g;
		like($res4, qr/^{"a":\["b","c"\]}$/, '{ "a": [ "b", "c" ] }');


		my $res5 = toJSON([ "a", { "b" => "c" } ]);
		unlike($res5, qr/"b"\s:/, "No space before ':' key-value separator (for Groovy JSON parser (OTP))");
		$res5 =~ s/\s//g;
		like($res5 , qr/^\["a",{"b":"c"}\]$/, '[ "a", { "b": "c" } ]');

		my $res6 = toJSON({ "1" => 5, "100" => 3 , "10" => 4, "100000" => 1, "1000" => "b"});
		$res6 =~ s/\s//g;
		like($res6, qr/{"1":5,"10":4,"100":3,"1000":"b","100000":1}/, "alphabetic key order")
    };

    subtest 'mapQcValues' => sub {
	my $res1 = mapQcValues(-genomeCoverage => parseGenomeCoverage(@testGenomeCoverage),
					       -insertSize => parseInsertSize(qw(399 23 93)),
			       		   -flagstats => parseFlagstats(@testFlagstats),
			       		   -diffChrom => parseDiffChrom(1.55));
	my $exp1 = {
          '1' => {
                   'genomeWithoutNCoverageQcBases' => '0.01',
                   'chromosome' => '1',
                   'referenceLength' => '249250621',
                   'qcBasesMapped' => '2885357'
                 },
          'Y' => {
                   'chromosome' => 'Y',
                   'referenceLength' => '59373566',
                   'genomeWithoutNCoverageQcBases' => '0.02',
                   'qcBasesMapped' => '114475'
                 },
          'all' => {
                     'withMateMappedToDifferentChr' => '33635',
                     'singletons' => '1080',
                     'insertSizeCV' => '23',
                     'withItselfAndMateMapped' => '419289',
                     'insertSizeMedian' => '399',
                     'properlyPaired' => '384766',
                     'pairedRead1' => '209146',
                     'chromosome' => 'all',
                     'qcFailedReads' => '0',
                     'pairedInSequencing' => '421309',
                     'totalReadCounter' => '421309',
                     'duplicates' => '805',
                     'insertSizeSD' => '93',
                     'totalMappedReadCounter' => '420369',
                     'qcBasesMapped' => '35857865',
                     'percentageMatesOnDifferentChr' => '1.55',
                     'referenceLength' => '3095677412',
                     'genomeWithoutNCoverageQcBases' => '46.03',
                     'withMateMappedToDifferentChrMaq' => '6161',
                     'pairedRead2' => '212163'
                   }
        };
		is_deeply($res1, $exp1, "mapQcValues");
    };

    exit 0;
}


sub readFile ($) {
    my ($infile) = @_;
    my $ifh = IO::File->new("<" . $infile)
		|| confess("Cannot open file '$infile'");
    my @lines = <$ifh>;
    close($ifh);
    chomp @lines;
    return @lines;
}


sub parseGenomeCoverage (@) {
    my @lines = @_;
    if (scalar(@lines) <= 1) {
		croak "Empty genome coverage file";
    }
    my $colIndex = 0;
    my %colNames = map { $_ => $colIndex++ } split(/\t/, shift(@lines));
    my $interval = $colNames{"interval"};
    my %chrDat;
    foreach my $line (@lines) {
		my @line = split(/\t/, $line);
		if (scalar(@line) != scalar(keys %colNames)) {
		    carp "Unexpected number of columns in depth of coverage file '$line'";
		}
		$chrDat{$line[$interval]} = {
		    map {
				$_ => $line[$colNames{$_}]
		    } keys %colNames
		};
    }
    return \%chrDat;
}

sub ensureDef (@) {
    my ($val) = @_;
    if (!defined $val) {
		confess("Parse error");
    } else {
		return $val;
    }
}


sub extract ($$) {
    my ($str, $pattern) = @_;
    my @res = $str =~ $pattern;
    if (scalar(@res) == 0) {
		confess("Parse error in '$str' with pattern /$pattern/");
    }
    return @res;
}

=head1

    Take two list refs and return a list ref with the elements of the input lists interleaved, starting with\
    the first element of the first list. The two input lists need to have identical lengths.

=cut
sub zip ($$) {
    my ($alist, $blist) = @_;
    my @res = ();
    if (scalar(@$alist) != scalar(@$blist)) {
		confess "Different numbers of elements in lists";
    }
    for (my $idx = 0; $idx < scalar(@$alist); ++$idx) {
		push(@res, ($alist->[$idx], $blist->[$idx]));
    }
    return \@res;
}


sub ensureInt ($) {
    my ($val) = @_;
    if ($val !~ /^\d$/) {
		confess("Expected integer, got '$val'!");
    }
    return $val;
}


sub ensureFloat ($) {
    my ($val) = @_;
    if ($val !~ /^(?:-?\d+(?:\.\d+)?|\d(?:\.\d+)?(?:[eE]-?\d+)?)$/) {
		confess("Expected float, got '$val'!");
    }
    return $val;
}


sub parseInsertSize (@) {
    my @lines = @_;
    return {
		-median => ensureFloat(ensureDef($lines[0])),
		-varcoef => ensureFloat(ensureDef($lines[1])),
		-stddev => ensureFloat(ensureDef($lines[2]))
    };
}


sub parseDipStatistics (@) {
    return {};
}


sub parseFlagstats (@) {
    my @lines = @_;
    return {
		@{zip([-total, -total_qcfailed],
		      [extract($lines[0], qr/^(\d+) \+ (\d+) in total \(QC-passed reads \+ QC-failed reads\)$/)])
		},
		@{zip([-duplicates, -duplicates_qcfailed],
		      [extract($lines[1], qr/^(\d+) \+ (\d+) duplicates$/)])
		},
		@{zip([-mapped, -mapped_qcfailed, -percentage_mapped],
		      [extract($lines[2], qr/^(\d+) \+ (\d+) mapped \((\d+\.\d+)%:-nan%\)$/)])
		},
		@{zip([-paired_in_sequencing, -paired_in_sequencing_qcfailed],
		      [extract($lines[3], qr/^(\d+) \+ (\d+) paired in sequencing$/)])
		},
		@{zip([-read1, -read1_qcfailed],
		      [extract($lines[4], qr/^(\d+) \+ (\d+) read1$/)])
		},
		@{zip([-read2, -read2_qcfailed],
		      [extract($lines[5], qr/^(\d+) \+ (\d+) read2$/)])
		},
        @{zip([-properly_paired, -properly_paired_qcfailed, -percentage_properly_paired],
		    [extract($lines[6], qr/^(\d+) \+ (\d+) properly paired \((\d+\.\d+)%:-nan%\)$/)])
		},
		@{zip([-with_itself_and_mate_mapped, -with_itself_and_mate_mapped_qcfailed],
		      [extract($lines[7], qr/^(\d+) \+ (\d+) with itself and mate mapped$/)])
		},
        @{zip([-singletons, -singletons_qcfailed, -percentage_singletons],
		    [extract($lines[8], qr/^(\d+) \+ (\d+) singletons \((\d+\.\d+)%:-nan%\)$/)])
		},
        @{zip([-with_mate_mapped_to_different_chr, -with_mate_mapped_to_different_chr_qcfailed],
		      [extract($lines[9], qr/^(\d+) \+ (\d+) with mate mapped to a different chr$/)])
		},
        @{zip([-with_mate_mapped_to_different_chr_mapQge5, -with_mate_mapped_to_different_chr_mapQge5_qcfailed],
		      [extract($lines[10], qr/^(\d+) \+ (\d+) with mate mapped to a different chr \(mapQ>=5\)$/)])
		}
    };
}



sub parseDiffChrom (@) {
    my @lines = @_;
    return {
		-percentage_mates_on_different_chr => ensureFloat($lines[0])
    };
}

sub mergeAllAndGroupedGenomeCoverages ($$) {
	my ($all, $grouped) = @_;
	map {
		if (exists $all->{$_}) {
			croak "Chromosome group '$_' exist in among the chromosomes";
		}
	} keys %$grouped;
	return { %$all, %$grouped };
}

=head

    Here, the parsed data is mapped into a structure of hashes and arrays that is going to
    to be written out as JSON.

=cut
sub mapQcValues (%) {
    my (%parseResults) = @_;
    my %results = (
		map {
		    my $chrData = $parseResults{-genomeCoverage}->{$_};

			my ($qcBasesNotN, $referenceLengthNotN) = $chrData->{"#QC bases/#total not_N bases"} =~ /^(\d+)\/(\d+)/;
			my $coverageNotN = $chrData->{"genome_w/o_N coverage QC bases"};
			$coverageNotN =~ s/x$//;

			my ($qcBases, $referenceLength) = $chrData->{"#QC bases/#total bases"} =~ /^(\d+)\/(\d+)/;
			my $coverage = $chrData->{"coverage QC bases"};
			$coverage =~ s/x$//;

	    	$_ => {
				"chromosome"                    => $chrData->{"interval"},

				"qcBasesMapped"                 => ensureDef($qcBases),
				"referenceLength"               => ensureDef($referenceLength),
				"coverageQcBases"               => ensureDef($coverage),

				"genomeWithoutNQcBasesMapped"   => ensureDef($qcBasesNotN),
				"genomeWithoutNReferenceLength" => ensureDef($referenceLengthNotN),
				"genomeWithoutNCoverageQcBases" => ensureFloat($coverageNotN),
			}
		} keys %{ $parseResults{-genomeCoverage} }
    );

    $results{"all"} = {
		%{$results{"all"}},
		"totalReadCounter" => $parseResults{-flagstats}->{-total},
		"qcFailedReads" => $parseResults{-flagstats}->{-total_qcfailed},
		"duplicates" => $parseResults{-flagstats}->{-duplicates},
		"totalMappedReadCounter" => $parseResults{-flagstats}->{-mapped},
		"pairedInSequencing" => $parseResults{-flagstats}->{-paired_in_sequencing},
		"pairedRead2" => $parseResults{-flagstats}->{-read2},
		"pairedRead1" => $parseResults{-flagstats}->{-read1},
		"properlyPaired" => $parseResults{-flagstats}->{-properly_paired},
		"withItselfAndMateMapped" => $parseResults{-flagstats}->{-with_itself_and_mate_mapped},
		"withMateMappedToDifferentChrMaq" => $parseResults{-flagstats}->{-with_mate_mapped_to_different_chr_mapQge5},
		"withMateMappedToDifferentChr" =>    $parseResults{-flagstats}->{-with_mate_mapped_to_different_chr},
		"singletons" => $parseResults{-flagstats}->{-singletons},
		"insertSizeMedian" => $parseResults{-insertSize}->{-median},
		"insertSizeSD" => $parseResults{-insertSize}->{-stddev},
		"insertSizeCV" => $parseResults{-insertSize}->{-varcoef},        ## New value. Kind of coefficient of variation but using median.
		"percentageMatesOnDifferentChr" => $parseResults{-diffChrom}->{-percentage_mates_on_different_chr} ## New value.
    };
    return \%results;
}


=head1

    Recursive JSON building.

=cut
sub jsonIndent ($) {
    my ($indent) = @_;
    return join("", (" ") x $indent);
}


sub jsonArray ($$) {
    my ($value, $indent) = @_;
    return join("\n",
		"[",
		join(",\n",
		     map {
				 json($_, $indent + INDENT_INCREMENT)
		     } @$value),
			 jsonIndent($indent) . "]",
			 "");
}


sub jsonHash ($$) {
    my ($value, $indent) = @_;
    return join("\n",
				"{",
				join(",\n",
		     		 map {
			 			jsonKeyValue($_, $value->{$_}, $indent + INDENT_INCREMENT)
		     		 }
		     		 sort { $a cmp $b }
		     		 keys %$value),
				jsonIndent($indent) . "}");
}


sub jsonKeyValue ($$) {
    my ($key, $value, $indent) = @_;
    return join("",
				(" ") x $indent,
				"\"$key\": ", json($value, $indent));
}


sub json ($$) {
    my ($value, $indent) = @_;
    if (ref($value) eq "ARRAY") {
		return jsonArray($value, $indent);
    } elsif (ref($value) eq "HASH") {
		return jsonHash($value, $indent);
    } elsif (looks_like_number($value)) {
		## Numbers are returned as is, i.e. as the originally parsed string and thus with unchanged precision.
		return $value;
    } else {
		if (!defined $value) {
		    confess "Undefined value!";
		}
		return '"' . $value . '"';
    }
}


sub toJSON ($) {
    my ($qcdata) = @_;
    return json($qcdata, 0) . "\n";
}

=head

    MAIN

=cut
my ($genomeCoverageFile,
	$groupedGenomeCoverageFile,
    $insertSizeFile,
    $flagstatsFile,
    $diffChromFile) = @ARGV;

if (defined $genomeCoverageFile && $genomeCoverageFile eq "test") {
    runTests();
} elsif (!defined $diffChromFile) {
    pod2usage(1);
}

my %parse_results = (
	-genomeCoverage => mergeAllAndGroupedGenomeCoverages(
		parseGenomeCoverage(readFile($genomeCoverageFile)),
		parseGenomeCoverage(readFile($groupedGenomeCoverageFile))
	),
    -insertSize => parseInsertSize(readFile($insertSizeFile)),
    -flagstats => parseFlagstats(readFile($flagstatsFile)),
    -diffChrom => parseDiffChrom(readFile($diffChromFile))
    );

print STDOUT toJSON(mapQcValues(%parse_results));


__END__

=head1 NAME

qcJson.pl     Read in a number of QC and statistics files and compile (some of) the data into a quality-control JSON.

=head1 SYNOPSIS

qcJson.pl test

qcJson.pl genomeCoverageFile groupedGenomeCoverage insertSizeFile flagstatsFile diffChromFile > outFile

genomeCoverageFile          output of genomeCoverage.d
groupedGenomeCoverage       derived from output of genomeCovarege.d, with coverage values for chromosomes matching and non-matching $prefix$chr$suffix
                            aggregated
flagstatsFile               samtools flagstat output
diffChromFile               _DiffChroms.png_qcValues.txt file

=head1 BACKGROUND

From discussion of Ivo, Kortine, Philip 13-07-2015:

=item coverage* (decided per project) $readgroup/coverage/control_ICGC_MB99.DepthOfCoverage_Genome.txt
=item insertsize
           - median*		$readgroup/insertsize_distribution/control_ICGC_MB99_insertsize_plot.png_qcValues.txt (Median, SD, ?)
           - 10 and 90% quartile ( not done yet ) ==> add to Joachims skript
           - (p-value Joachim Weischenfedt)* ==> ToDo $readgroup/insertsize_distribution/${pid}_HartigansDip.txt
=item % reads mapped*             $readgroup/flagstats/control_ICGC_MB99_merged.mdup.bam_flagstats.txt (mapped)
=item proper pairs*	        $readgroup/flagstats/control_ICGC_MB99_merged.mdup.bam_flagstats.txt (properly paired, convert to % of paired)
=item singletons*	                $readgroup/flagstats/control_ICGC_MB99_merged.mdup.bam_flagstats.txt (singletons, convert to % of total mapped)
=item % diffChrom*	        $readgroup/structural_variation/tumor_ICGC_MB99_merged.mdup.bam_DiffChroms.png_qcValues.txt
=item total number of reads	$readgroup/flagstats/control_ICGC_MB99_merged.mdup.bam_flagstats.txt ()
=item duplication rate	        $readgroup/flagstats/control_ICGC_MB99_merged.mdup.bam_flagstats.txt (convert to % of total mapped)

=head1 AUTHOR

Philip R. Kensche <p.kensche@dkfz.de>

