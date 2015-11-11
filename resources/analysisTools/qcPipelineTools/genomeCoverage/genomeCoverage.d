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
// rdmd --compiler=ldmd2 --build-only -I<path/to/BioD> -O -release -inline genomeCoverage.d
import bio.bam.reader;
import bio.bam.pileup;

import std.stdio, std.range, std.algorithm, std.getopt, std.conv, 
       std.parallelism, std.path, std.file, std.process, std.traits;

enum Mode { countReads, countBases }

struct Options {
    string output_file;
    
    size_t window_size_kb = 1;
    size_t window_size() @property const { return window_size_kb * 1000; }

    Mode mode = Mode.countReads;

    string countType() @property const { return mode == Mode.countReads ? "read" : "base"; }
}

size_t writeCounts(alias posGetter, alias countGetter, alias refIdGetter, R)
    (R items, BamReader reader, string out_fn, Options options)
{
    auto f = File(out_fn, "w+");
    f.setvbuf(1_024_576);

    long prev_window = long.min; string prev_chrom;
    size_t curr_count; long curr_window; string curr_chrom;

    size_t total;

    long lastWindow(string chrom) {
        return reader[chrom].length / options.window_size;
    }

    void writeEmptyWindows(string chrom, long from, lazy long to) {
        if (chrom == "") return; // handles first chromosome
        foreach (window; from .. to + 1)
            f.writeln(chrom, "\t", window * options.window_size, "\t0");
    }

    void dump() {
        if (prev_window != long.min)
            f.writeln(prev_chrom, "\t", prev_window * options.window_size, "\t", curr_count);
        total += curr_count;
        curr_count = 0;
        if (prev_chrom == curr_chrom) {
            writeEmptyWindows(prev_chrom, prev_window + 1, curr_window - 1);
        } else {
            writeEmptyWindows(prev_chrom, prev_window + 1, lastWindow(prev_chrom));
            writeEmptyWindows(curr_chrom, 0, curr_window - 1);
        }
        prev_window = curr_window;
        prev_chrom = curr_chrom;
    }
    
    foreach (item; items) {
        auto pos = posGetter(item).to!long;
        curr_window = pos / options.window_size;
        auto curr_ref_id = refIdGetter(item);
        if (curr_ref_id < 0)
            break;
        curr_chrom = reader.reference_sequences[curr_ref_id].name;
        if (curr_window > prev_window || curr_chrom != prev_chrom)
            dump();

        curr_count += countGetter(item);
    }

    dump();
    writeEmptyWindows(prev_chrom, prev_window + 1, lastWindow(prev_chrom));

    return total;
}

void process(string bam_fn, Options options) {
  writeln("processing ", bam_fn);
  auto bam = new BamReader(bam_fn);
  bam.setBufferSize(16_000_000);

  auto out_fn = buildPath(options.output_file);

  auto reads = bam.reads().filter!(r => !r.is_duplicate);
  size_t total;

  final switch (options.mode) {
  case Mode.countReads:
    total = writeCounts!(r => r.position, _ => 1, r => r.ref_id)
                        (reads.filter!(r => r.mapping_quality >= 1), bam, out_fn, options);
    break;
  case Mode.countBases:
    auto columns = reads.filter!(
           r => !r.is_secondary_alignment && !r.failed_quality_control && !r.is_unmapped
         ).pileupColumns();
    total = writeCounts!(c => c.position, c => c.coverage, c => c.ref_id)(columns, bam, out_fn, options);
    break;
  }
  writeln("total number of ", options.countType, "s counted for coverage analysis:");
  writeln(total);
}

void printUsage() {
    stderr.writeln("   usage: genomeCoverage OPTIONS\n"
                   "\n"
                   "      OPTIONS:\n"
                   "        --alignmentFile=FILENAME\n"
                   "          BAM file(s) to process\n"
                   "        --outputFile=OUTPUT_FILE\n"
                   "          Output File\n"
                   "        --processors=N\n"
                   "          How many threads to use (1 by default)\n"
                   "        --mode=countReads|countBases\n"
                   "          Operation mode (default is countReads)\n"
                   "        --windowSize=SIZE_IN_KBASES\n"
                   "          Window width (1kb by default)\n");
}

int main(string[] args) {

  string[] bam_filenames;
  ushort processors = 0;

  Options options;
 

  if (args.length == 1) {
    printUsage();
    return -1;
  }

  getopt(args,
         std.getopt.config.caseSensitive,
         "alignmentFile", &bam_filenames,
         "outputFile", &options.output_file,
         "mode", &options.mode,
         "windowSize", &options.window_size_kb,
         "processors", &processors);

  

  if (bam_filenames.length == 0) {
      stderr.writeln("specify BAM filenames with --alignmentFile option");
      return -1;
  }

  if (options.output_file.length == 0) {
      stderr.writeln("specify Output filename with --outputFile option");
      return -1;
  }

  defaultPoolThreads = processors == 0 ? 0 : processors - 1;

  foreach (bam_fn; bam_filenames)
    process(std.path.expandTilde(bam_fn), options);

  return 0;
}
