// Copyright (C) 2014      Artem Tarasov (lomereiter@gmail.com)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
// associated documentation files (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial
// portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
// LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH
// THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

// Command-line for compilation:
// rdmd --compiler=ldmd2 --build-only -I<path/to/BioD> -O -release -inline coverageQc.d

import bio.bam.reader;
import std.range, std.algorithm, std.array, std.string, std.typetuple, std.math;
import std.getopt, std.stdio, std.parallelism, std.path;

struct Statistics {
  string reference_name;
  ulong reference_length;
  ulong ungapped_ref_size;
  ulong[2] duplicate_count;
  ulong incorrect_proper_pairs;
  ulong total_orientation_counter;
  ulong bad_orientation_counter;
  ulong qc_bases;

  private ulong[2][4] _data;
  alias BaseQualStats = TypeTuple!("mapq=0",
                                   "mapq>0,readlength<minlength",
                                   "mapq>0,BaseQualityMedian<basequalCutoff",
                                   "mapq>0,BaseQualityMedian>=basequalCutoff");

  ref inout(ulong[2]) base_qual(string description)() @property inout {
    enum index = staticIndexOf!(description, BaseQualStats);
    static assert(index != -1, "invalid description");
    return _data[index];
  }

  static string[] columns =
    ["interval",
     "coverage QC bases",
     "#QC bases/#total bases",
     "mapq=0 read1",
     "mapq=0 read2",
     "mapq>0,readlength<minlength read1",
     "mapq>0,readlength<minlength read2",
     "mapq>0,BaseQualityMedian<basequalCutoff read1",
     "mapq>0,BaseQualityMedian<basequalCutoff read2",
     "mapq>0,BaseQualityMedian>=basequalCutoff read1",
     "mapq>0,BaseQualityMedian>=basequalCutoff read2",
     "%incorrect PE orientation",
     "#incorrect proper pair",
     "#duplicates read1 (excluded from coverage analysis)",
     "#duplicates read2 (excluded from coverage analysis)",
     "genome_w/o_N coverage QC bases",
     "#QC bases/#total not_N bases"];

  string toString() const {
    auto orientationPE = total_orientation_counter == 0 ? "NA" :
      format("%.3f", bad_orientation_counter*100.0/total_orientation_counter);

    return std.string.format("%s\t%.2fx\t%s/%s\t%s%s\t%s\t%s\t%s\t%.2fx\t%s/%s",
                             reference_name,
                             qc_bases.to!float / reference_length,
                             qc_bases, reference_length,
                             baseQualString(),
                             orientationPE,
                             incorrect_proper_pairs,
                             duplicate_count[0], duplicate_count[1],
                             qc_bases.to!float / ungapped_ref_size,
                             qc_bases, ungapped_ref_size);
  }

  private string baseQualString() const {
    string result;
    foreach (stat; BaseQualStats)
      foreach (index; 0 .. 2)
        result ~= base_qual!stat[index].to!string() ~ "\t";
    return result;
  }
}

struct Options { int min_read_length = 36; uint base_qual_cutoff = 25; }

void addRead(ref Statistics stat, BamRead read, Options options) {
  with (stat) {
    if (read.is_duplicate) {
      if (read.is_first_of_pair)
        duplicate_count[0] += 1;
      else if (read.is_second_of_pair)
        duplicate_count[1] += 1;
      else
        duplicate_count[0] += 1;
      return;
    }

    if (read.proper_pair) 
      if (read.ref_id == read.mate_ref_id)
        if (read.is_reverse_strand == read.mate_is_reverse_strand)
          incorrect_proper_pairs += 1;

    if (!read.is_unmapped)
      if (read.ref_id == read.mate_ref_id)
        if (!read.mate_is_unmapped) {
          total_orientation_counter += 1;
          if (read.is_reverse_strand == read.mate_is_reverse_strand)
	  // the case that both map to the forward strand is not considered in the original script either
	  // this should anyways return NO reads at all, since alignedRead.is_proper_pair
	  // does not allow mapping to different chromosomes or to the same strand
            bad_orientation_counter += 1;
        }

    size_t peRead = read.is_second_of_pair ? 1 : 0;

    if (read.mapping_quality > 0 && !read.is_unmapped) {
      auto qualities = read.base_qualities;
      int sum_qual, non_softclipped_length;
      foreach (op; read.cigar) {
        if (!op.is_query_consuming) continue; // skip D/N
        if (op.type != 'S') {
        foreach (i; 0 .. op.length) sum_qual += qualities[i];
          non_softclipped_length += op.length;
        }  
        qualities = qualities[op.length .. $];
      }
      if (non_softclipped_length >= options.min_read_length) {
        auto avg_qual = round(sum_qual.to!double / non_softclipped_length);
        if (avg_qual >= options.base_qual_cutoff) {
          base_qual!"mapq>0,BaseQualityMedian>=basequalCutoff"[peRead] += 1;
          //qc_bases += read.sequence_length;
          	
		 qc_bases += non_softclipped_length;	
        } else {
          base_qual!"mapq>0,BaseQualityMedian<basequalCutoff"[peRead] += 1;
        }
      } else {
        base_qual!"mapq>0,readlength<minlength"[peRead] += 1;
      }
    } else {
      base_qual!"mapq=0"[peRead] += 1;
    }
  }
}

void printUsage() {
  stderr.writeln("usage: coverageQc OPTIONS\n");
  stderr.writeln("OPTIONS:  --alignmentFile=BAMFILE");
  stderr.writeln("             (Required) Specify the name of the input file containing short read alignments, use /dev/stdin to read from stdin");
  stderr.writeln("          --outputFile=OUTPUTFILE");
  stderr.writeln("             (Required) Specify the name of the output file, use /dev/stdin to write to stout");
  stderr.writeln("          --processors=N");
  stderr.writeln("             Specify how many threads to use (default: 1)");
  stderr.writeln("          --ungappedSizes=UNGAPPED_SIZE_LIST");
  stderr.writeln("             Specify the file with ungapped sizes of the chromosomes (default: ungapped_sizes_default.lst)");
  stderr.writeln("          --basequalCutoff=CUTOFF");
  stderr.writeln("             Specify the minimum mean of quality scores of the complete read length used as a quality cutoff (default: 25 (in phred score))");
  stderr.writeln("          --minReadLength=K");
  stderr.writeln("             Specify the minimum read length (default: 36)");
  stderr.writeln("          --targetsize=T");
  stderr.writeln("             Specify the target region size (default: 0, will use sum of ungapped chromosome sizes)");
  stderr.writeln("          --excludedChromosomes=CHROMOSOMES");
  stderr.writeln("             Specify in a comma-separated string which chromosomes should be excluded (it is useful to exclude chrY when aligning to the female genome)");
}

int main(string[] args) {

  immutable exe_dir = dirName(args[0]);
    
  string[] bam_filenames;
  string ungapped_sizes_fn = buildPath(exe_dir, "ungapped_sizes_default.lst");
  string output_file;
  ushort processors = 0;
  ulong targetsize = 0;
  string skiplist;
  bool[string] chromosomes_to_skip;

  Options options;

  if (args.length == 1) {
    printUsage();
    return -1;
  }

  getopt(args,
         std.getopt.config.caseSensitive,
         "alignmentFile", &bam_filenames,
         "ungappedSizes", &ungapped_sizes_fn,
         "outputFile", &output_file,
         "processors", &processors,
         "excludedChromosomes", &skiplist,
         "basequalCutoff", &options.base_qual_cutoff,
         "minReadLength", &options.min_read_length,
	 "targetsize", &targetsize);

  foreach (chrom; skiplist.split(","))
    chromosomes_to_skip[chrom] = true;

  if (bam_filenames.length == 0) {
      stderr.writeln("specify BAM filenames with --alignmentFile option, use /dev/stdin to read from stdin");
      return -1;
  }

  if (output_file == "") {
      stderr.writeln("specify output file with --outputFile option, use /dev/stdout to write to stdout");
      return -1;
  }

  ulong[string] ungapped_sizes;
  auto ungapped_f = File(ungapped_sizes_fn);
  foreach (line; ungapped_f.byLine()) {
    auto fields = std.array.split(line, "\t");
    if (fields.length < 2) continue;
    ungapped_sizes[fields[0].idup] = fields[1].to!ulong;
  }

  defaultPoolThreads = processors == 0 ? 0 : processors - 1;

  foreach (bam_fn; bam_filenames) {
    auto bam = new BamReader(bam_fn);
    bam.assumeSequentialProcessing();
    auto per_ref_stats = new Statistics[](bam.reference_sequences.length);
    foreach (i; 0 .. per_ref_stats.length) {
      per_ref_stats[i].reference_name = bam.reference_sequences[i].name;
      per_ref_stats[i].reference_length = bam.reference_sequences[i].length;
      auto name = per_ref_stats[i].reference_name;
      if (name in ungapped_sizes)
        per_ref_stats[i].ungapped_ref_size = ungapped_sizes[name];
    }

    foreach (read; bam.reads) {
      if (read.ref_id < 0 || read.ref_id >= per_ref_stats.length) continue;
      if (per_ref_stats[read.ref_id].ungapped_ref_size == 0) continue;
      auto chrom = per_ref_stats[read.ref_id].reference_name;
      if (chrom in chromosomes_to_skip) continue;
      per_ref_stats[read.ref_id].addRead(read, options);
    }

    auto f = File(buildPath(output_file), "w+");
    f.writeln(Statistics.columns.join("\t"));
    foreach (stats; per_ref_stats)
      if (stats.ungapped_ref_size > 0)
        f.writeln(stats);

    Statistics overall_stats;
    with (overall_stats) {
      foreach (stats; per_ref_stats) {
        if (stats.reference_name in chromosomes_to_skip) continue;
        if (stats.ungapped_ref_size == 0) continue;
        reference_length += stats.reference_length;
        ungapped_ref_size += stats.ungapped_ref_size;
        qc_bases += stats.qc_bases;
        bad_orientation_counter += stats.bad_orientation_counter;
        total_orientation_counter += stats.total_orientation_counter;
        incorrect_proper_pairs += stats.incorrect_proper_pairs;
        foreach (stat; Statistics.BaseQualStats)
          base_qual!stat[] += stats.base_qual!stat[];
        duplicate_count[] += stats.duplicate_count[];


      }
      // hack: for exome, divide by the target size to get average target coverage over all chrosomomes
      if (targetsize != 0){
        ungapped_ref_size = targetsize;
      }
      reference_name = "all";
    }
    f.writeln(overall_stats);
    
    f.close();
  }

  return 0;
}
