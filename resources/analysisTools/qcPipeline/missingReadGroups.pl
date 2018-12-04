#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

use strict;
use warnings;

if (@ARGV < 4)
{
        die "need 3 arguments: 1) String with : separated BAM paths, 2) string with the RG tags including the names, 3) paired bam suffix 4) sample type\n";
}

my $bampaths = shift;
my $readGroups = shift;
my $pairedBamSuffix = shift;
my $sampletype = shift;

#print STDERR "BAMS: $bampaths\nRGs: $rgs\n";

# BAM paths:
#/icgc/dkfzlsdf/analysis/hipo_021/results_per_pid/P021-W00D/STAR_NEW/run140113_SN7001397_0115_BH871TADXX_P021-W00D-T1-R1_ACTGAT_L002.bam:/icgc/dkfzlsdf/analysis/hipo_021/results_per_pid/P021-W00D/STAR_NEW/run140113_SN7001397_0115_BH871TADXX_P021-W00D-T1-R1_ACTGAT_L001.bam:

# RG tags:
#@RG     ID:run140113_SN7001397_0115_BH871TADXX_P021-W00D-T1-R1_ACTGAT_L001      PL:ILLUMINA     PU:run140113_SN7001397_0115_BH871TADXX_P021-W00D-T1-R1_ACTGAT_L001      LB:P021-W00D_tumor  SM:P021-W00D_tumor @RG     ID:run140113_SN7001397_0115_BH871TADXX_P021-W00D-T1-R1_ACTGAT_L002      PL:ILLUMINA     PU:run140113_SN7001397_0115_BH871TADXX_P021-W00D-T1-R1_ACTGAT_L002      LB:P021-W00D_tumor  SM:P021-W00D_tumor

# already present lanes from RG:
my %present = ();
my $lanes = 0;
while ($readGroups =~ /ID:([\w+-]+)\t/g)
{
        #print STDERR "$1\n";
        $present{$1} = 1;
        $lanes++;
}

print STDERR "$lanes lanes are already present in the merged BAM\n";
my @bams = split (/:/, $bampaths);

# get basenames for them and remove .bam, look up if it's already merged
my $missing = "";
for (my $lb = 0; $lb < @bams; $lb++)
{
        # convert BAM name. e.g.
        #/icgc/dkfzlsdf/analysis/hipo_021/results_per_pid/P021-SBD8/alignment/blood_run140212_SN7001393_0144_AH8KAFADXX_P021-SBD8-B1-D1_ACTTGA_L001_paired.bam.sorted.bam
        # into ID, e.g.
        # run140212_SN7001393_0144_AH8KAFADXX_P021-SBD8-B1-D1_ACTTGA_L001
        my $name =`basename $bams[$lb]`;
        chomp $name;
        # the ID is composed of RUN and LANE (which often includes the PID)
        # but does not have the sampletype (tumor, blood, ...) as a prefix like the BAM
        # $name =~ s/^\w+_run/run/;	# evil hardcoding: "run" prefix is only always present for data management's view-by-pid
	$name =~ s/^${sampletype}_//;
        # and it also does not have the suffix _paired.bam.sorted.bam
        # the suffix is probably somewhere in the Roddy setup but not in runtimeConfig.sh
        # so it is hardcoded here!
        $name =~ s/_$pairedBamSuffix$//;
        print STDERR "ID: $name";
        if (exists $present{$name})
        {
                print STDERR " is already merged";
        }
        else
        {
               print STDERR " is missing";
               $missing.=$bams[$lb].":";
        }
        print STDERR "\n";
}
print $missing;
exit;
