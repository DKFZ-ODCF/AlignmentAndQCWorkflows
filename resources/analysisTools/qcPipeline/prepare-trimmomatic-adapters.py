#!/usr/bin/env python3
#
# Philip R. Kensche (2021)
#
# Requirements:
#
# * Python 3.7
#
# Run the test with:
#
#  pytest prepare-trimmomatic-adapters.py
#
# Run:
#
# ./prepare-trimmomatic-adapters.py -h
#
from __future__ import annotations  # Require Python >= 3.7

import bz2
import gzip
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from concurrent import futures
from concurrent.futures import Future
from contextlib import contextmanager
from copy import deepcopy
from enum import Enum
from io import TextIOBase
from pathlib import Path
from textwrap import dedent
from typing import List, Optional, Dict, Union

PREFIX = "Prefix_"

# See https://www.bioinformatics.org/sms/iupac.html.
dna_complementary_base = {
    'a': 't', 'A': 'T',
    't': 'a', 'T': 'A',
    'c': 'g', 'C': 'G',
    'g': 'c', 'G': 'C',
    's': 's', 'S': 'S',  # strong
    'w': 'w', 'W': 'W',  # weak
    'n': 'n', 'N': 'N',
    'r': 'y', 'R': 'Y',  # purine, pyrimidine
    'y': 'r', 'Y': 'R',
    'k': 'm', 'K': 'M',  # keto, amino
    'm': 'k', 'M': 'K'
}


def revcom(seq: str) -> str:
    try:
        return "".join(list(map(lambda base: dna_complementary_base[base], reversed(seq))))
    except KeyError as e:
        raise ValueError("Not a DNA sequence", e)


def test_revcom():
    assert revcom("ATCGSWNRYKM") == "KMRYNWSCGAT"
    assert revcom("atcgswnrykm") == "kmrynwscgat"
    assert revcom("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGAAA") == "TTTCTCTTTCCCTACACGACGCTCTTCCGATCT"
    import pytest
    with pytest.raises(ValueError) as ex:
        revcom("x")


class ConstructEnd(Enum):
    both = 0
    one = 1
    two = 2


class Adapter:

    def __init__(self, name: str, sequence: str, end: ConstructEnd):
        if any(map(lambda base: base not in dna_complementary_base.keys(),
                   sequence)):
            raise ValueError(f"Non-DNA symbol in sequence: {name}={sequence}")
        self._name = name
        self._sequence = sequence
        self._end = end

    @property
    def name(self) -> str:
        return self._name

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def end(self) -> ConstructEnd:
        return self._end

    def to_fasta(self) -> str:
        fasta_id = None
        if self.end == ConstructEnd.one:
            fasta_id = self.name + "/1"
        elif self.end == ConstructEnd.two:
            fasta_id = self.name + "/2"
        elif self.end == ConstructEnd.both:
            fasta_id = self.name
        return f">{fasta_id}\n{self.sequence}\n"

    def reverse_complement(self) -> Adapter:
        return Adapter(name=self.name, sequence=revcom(self.sequence), end=self.end)

    def copy(self, name: Optional[str] = None,
             sequence: Optional[str] = None,
             end: Optional[ConstructEnd] = None):
        return Adapter(name if name is not None else self.name,
                       sequence if sequence is not None else self.sequence,
                       end if end is not None else self.end)


def test_adapter():
    import pytest
    with pytest.raises(ValueError):
        Adapter(name="name", sequence="!!!!", end=ConstructEnd.both)

    adapter = Adapter(name="name",
                      sequence="NSYTCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                      end=ConstructEnd.one)
    rev_adapter = adapter.reverse_complement()
    assert rev_adapter.name == adapter.name
    assert rev_adapter.end == adapter.end
    assert adapter.sequence == "NSYTCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    assert rev_adapter.sequence == "ACACTCTTTCCCTACACGACGCTCTTCCGARSN"

    assert adapter.to_fasta() == ">name/1\nNSYTCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n"
    assert rev_adapter.to_fasta() == ">name/1\nACACTCTTTCCCTACACGACGCTCTTCCGARSN\n"

    adapter2 = Adapter(name="name",
                       sequence="NNNTCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                       end=ConstructEnd.two)
    assert adapter2.end == ConstructEnd.two
    assert adapter2.to_fasta() == ">name/2\nNNNTCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n"

    assert adapter.copy(name="newName").name == "newName"
    assert adapter.copy(sequence="ACG").sequence == "ACG"
    assert adapter.copy(end=ConstructEnd.both).end == ConstructEnd.both


class AdapterPair:

    def __init__(self, read1: Adapter, read2: Adapter):
        if read1.name != read2.name:
            raise ValueError(f"read1 and read2 must have same name: {read1.name} vs. {read2.name}")
        if read1.end != ConstructEnd.one:
            raise ValueError(f"read1 is not a Read.one: {read1.end}")
        if read2.end != ConstructEnd.two:
            raise ValueError(f"read1 is not a Read.two: {read2.end}")
        self._read1 = read1
        self._read2 = read2

    @property
    def read1(self) -> Adapter:
        return deepcopy(self._read1)

    @property
    def palindrome_read1(self) -> Adapter:
        return Adapter(name=f"{PREFIX}{self.read2.name}",     # Add the "Prefix"
                       sequence=revcom(self.read2.sequence),  # Reverse complement
                       end=ConstructEnd.one)                  # Switch read 1 and read 2

    @property
    def read2(self) -> Adapter:
        return deepcopy(self._read2)

    @property
    def palindrome_read2(self) -> Adapter:
        return Adapter(name=f"{PREFIX}{self.read1.name}",     # Add the "Prefix"
                       sequence=revcom(self.read1.sequence),  # Reverse complement
                       end=ConstructEnd.two)                  # Switch read 1 and read 2

    @property
    def name(self) -> str:
        return self._read1.name

    def to_palindrome_fasta(self) -> str:
        """
        The palindrome sequence has the /1 and /2 labels switched, and both sequences have to
        be provided as reverse-complement of what is found in the FASTQ files. Finally, both
        read identifiers must be prefixed with "Prefix".
        """
        return self.palindrome_read1.to_fasta() + self.palindrome_read2.to_fasta()


def test_adapter_pair_to_fasta():
    import pytest
    adapter1 = Adapter(name="TheName", sequence="ACTG", end=ConstructEnd.one)
    adapter2 = Adapter(name="TheName", sequence="RYTW", end=ConstructEnd.two)
    pair = AdapterPair(adapter1, adapter2)
    assert pair.to_palindrome_fasta() == \
        f">{PREFIX}TheName/1\nWARY\n" +\
        f">{PREFIX}TheName/2\nCAGT\n"

    with pytest.raises(ValueError) as ex:
        AdapterPair(adapter2, adapter1)

    with pytest.raises(ValueError):
        AdapterPair(adapter1, Adapter(name="other", sequence="", end=ConstructEnd.two))


def parse_parameters(args: List[str]):
    parser = ArgumentParser(
        description=dedent("""
        Create an adapter file for Trimmomatic. Supports output for paired-end adapters, 
        read 1 only, or read 2 only. 

        For paired-end adapters, both simple-mode and palindrome-mode adapter sequences can be
        created. Note that the script does *not* prevent you from creating a Trimmomatic adapter 
        file with both, simple-mode and palindrome-mode adapter sequences, if you request them.
        Usually you only want the palindrome-mode sequences. 

        IMPORTANT: Sequences must be provided in exactly the way they will be found at the end
                   of reads.
        """),
                           formatter_class=RawTextHelpFormatter)
    parser.add_argument('-a', '--adapter', dest="adapters", action="append", default=[],
                        help="Simple mode adapters to be trimmed from *both* read's ends. " +
                        "Format: $name=$sequence")
    parser.add_argument('-a1', '--adapter1', dest="adapters1", action="append", default=[],
                        help="Simple mode adapters at end of read 1. " +
                        "Format: $name=$sequence")
    parser.add_argument('-a2', '--adapter2', dest="adapters2", action="append", default=[],
                        help="Simple mode adapters at end of read 2. " +
                        "Format: $name=$sequence")
    parser.add_argument('-p', '--paired-adapters', dest="paired_adapters",
                        action="append", default=[],
                        help="Palindrome mode adapter sequences. " +
                        "Format: $name=$read1,$read2")
    parser.add_argument('-f1', '--fastq1', dest="fastq1", default=None,
                        help="Optional R1 FASTQ. Input adapter 1 will be checked against this.")
    parser.add_argument('-f2', '--fastq2', dest="fastq2", default=None,
                        help="Optional R2 FASTQ. Input adapter 2 will be checked against this.")
    parser.add_argument('-T', '--truncate', dest="truncate", type=int, default=None,
                        help="Truncate adapters to this length. Useful for long inserts or " +
                        "problems with adapters that contain variable index sequences.")
    parser.add_argument('-N', '--first-n', dest='first_n', type=int, default=100000,
                        help="Read the first N reads of the input FASTQs. Default 100000.")
    if len(args) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args(args[1:])


def parse_simple_adapter(adapter_opt: str, read: ConstructEnd) -> Adapter:
    split_res = adapter_opt.split("=")
    if len(split_res) != 2:
        raise ValueError(f"Could not parse: '{adapter_opt}")
    if len(split_res[1].split(",")) != 1:
        raise ValueError(f"Simple mode adapter are not in pairs. Found comma in {adapter_opt}")
    return Adapter(name=split_res[0], sequence=split_res[1], end=read)


def test_parse_simple_adapter():
    import pytest
    adapter = parse_simple_adapter("name=ACGT", ConstructEnd.both)
    assert adapter.name == "name"
    assert adapter.sequence == "ACGT"
    assert adapter.end == ConstructEnd.both

    adapter1 = parse_simple_adapter("namo=TGC", ConstructEnd.one)
    assert adapter1.name == "namo"
    assert adapter1.sequence == "TGC"
    assert adapter1.end == ConstructEnd.one

    with pytest.raises(ValueError):
        parse_simple_adapter("ACGT", ConstructEnd.one)

    with pytest.raises(ValueError):
        parse_simple_adapter("namo=ACGT,CGT", ConstructEnd.both)


def parse_palindrome_adapters(adapter_opt: str) -> AdapterPair:
    name_split = adapter_opt.split("=")
    if len(name_split) != 2:
        raise ValueError(f"Could not parse name from: '{adapter_opt}")
    name = name_split[0]
    seq_split = name_split[1].split(",")
    if len(seq_split) != 2:
        raise ValueError(f"Could not parse read 1 and 2 from: '{adapter_opt}")
    seq1, seq2 = seq_split
    return AdapterPair(Adapter(name=name, sequence=seq1, end=ConstructEnd.one),
                       Adapter(name=name, sequence=seq2, end=ConstructEnd.two))


def test_parse_palindrome_adapters():
    adapters = parse_palindrome_adapters("Name=ACGT,GTCA")
    assert adapters.name == "Name"
    assert adapters.to_palindrome_fasta() == \
        f">{PREFIX}{adapters.name}/1\nTGAC\n" +\
        f">{PREFIX}{adapters.name}/2\nACGT\n"

    adapters2 = parse_palindrome_adapters(
        "WH=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC,CTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT")
    assert adapters2.to_palindrome_fasta() == \
        f">{PREFIX}{adapters2.name}/1\n" + \
        "AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG\n" + \
        f">{PREFIX}{adapters2.name}/2\n" + \
        "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG\n"


@contextmanager
def open_decompressed(file: Union[str, Path]) -> TextIOBase:
    _file = Path(file) if isinstance(file, str) else file
    if _file.suffix == ".gz":
        yield gzip.open(_file, "rt", newline="\n")
    elif _file.suffix == ".bz2" or _file.suffix == ".bzip2":
        yield bz2.open(_file, "rt", newline="\n")
    else:
        yield open(_file, "r")


def next_read_sequence(file_stream: TextIOBase) -> Optional[str]:
    """
    Returns None if the stream has ended. E.g. relevant for truncated input files.
    """
    id = file_stream.readline()
    seq = file_stream.readline()
    plus = file_stream.readline()
    quals = file_stream.readline()
    if seq == "":
        # readline() returns an empty string, then EOF has been reached.
        if id != "":
            print("Possibly truncated file. Could read read-ID but not sequence.", file=sys.stderr)
        else:
            return None
    else:
        return seq


class AdapterCounter:

    def __init__(self, adapters: List[Adapter], truncate: Optional[int] = None):
        if len(adapters) == 0:
            raise ValueError("No adapters")
        self._ensure_unique_adapter_name(adapters)
        self._ensure_compatible_read_ends(adapters)
        self._truncate = truncate
        if truncate is not None:
            self._match_seqs = \
                list(map(lambda a: (a.name, a.sequence[0:min(len(a.sequence), truncate)]),
                         adapters))
        else:
            self._match_seqs = \
                list(map(lambda a: (a.name, a.sequence),
                         adapters))
        self._result = \
            dict(map(lambda a: (a.name, 0),
                     adapters))

    def _ensure_compatible_read_ends(self, adapters: List[Adapter]) -> None:
        relevant_adapters = list(filter(lambda a: a.end != ConstructEnd.both, adapters))
        if len(relevant_adapters) > 0:
            if any(map(lambda a: a.end != relevant_adapters[0].end, relevant_adapters)):
                raise RuntimeError("Adapters with incompatible read-association")

    def _ensure_unique_adapter_name(self, adapters: List[Adapter]) -> None:
        # We check that all adapters have different identifiers. Otherwise, the counts would
        # be pooled.
        adapter_name_seen: Dict[str, Adapter] = {}
        for adapter in adapters:
            if adapter.name in adapter_name_seen.keys():
                raise ValueError(f"Adapter name '{adapter.name}' used for multiple adapters: " +
                                 adapter.sequence + " and " +
                                 adapter_name_seen[adapter.name].sequence)
            else:
                adapter_name_seen[adapter.name] = adapter

    def match_sequences(self) -> List[(str, str)]:
        return deepcopy(self._match_seqs)

    def count_matches(self, query: str) -> None:
        for (name, seq) in self._match_seqs:
            self._result[name] += 0 if query.find(seq) == -1 else 1

    def result(self) -> Dict[str, int]:
        return self._result


def test_adapter_counter():
    import pytest
    cnt1 = AdapterCounter([Adapter(name="a", sequence="ACTG", end=ConstructEnd.one),
                           Adapter(name="b", sequence="AAAA", end=ConstructEnd.one)])
    assert cnt1.match_sequences() == [("a", "ACTG"), ("b", "AAAA")]
    assert cnt1.result() == {"a": 0, "b": 0}

    cnt1.count_matches("ACTGAAAA")
    assert cnt1.result() == {"a": 1, "b": 1}

    cnt1.count_matches("NNNACTGAA")
    assert cnt1.result() == {"a": 2, "b": 1}

    cnt1.count_matches("NNAAAA")
    assert cnt1.result() == {"a": 2, "b": 2}

    cnt2 = AdapterCounter([Adapter(name="a", sequence="ACTGCCCC", end=ConstructEnd.one),
                           Adapter(name="b", sequence="AAAACCC", end=ConstructEnd.one)],
                          truncate=4)
    assert cnt2.match_sequences() == [("a", "ACTG"), ("b", "AAAA")]

    with pytest.raises(ValueError):
        AdapterCounter([Adapter(name="a", sequence="ACTGCCCC", end=ConstructEnd.one),
                        Adapter(name="a", sequence="AAAACCC", end=ConstructEnd.one)])

    with pytest.raises(RuntimeError):
        AdapterCounter([Adapter(name="a", sequence="ACTGCCCC", end=ConstructEnd.one),
                        Adapter(name="b", sequence="AAAACCC", end=ConstructEnd.two)])

    AdapterCounter([Adapter(name="a", sequence="ACTGCCCC", end=ConstructEnd.both),
                    Adapter(name="b", sequence="AAAACCC", end=ConstructEnd.both)])

    AdapterCounter([Adapter(name="a", sequence="ACTGCCCC", end=ConstructEnd.two),
                    Adapter(name="b", sequence="AAAACCC", end=ConstructEnd.both)])

    AdapterCounter([Adapter(name="a", sequence="ACTGCCCC", end=ConstructEnd.both),
                    Adapter(name="b", sequence="AAAACCC", end=ConstructEnd.one)])


def count_adapters_in_fastq_stream(file_stream: TextIOBase, adapters: List[Adapter],
                                   first_n: Optional[int] = 100000,
                                   truncate: Optional[int] = None) -> Dict[str, int]:
    adapter_counter = AdapterCounter(adapters, truncate)
    n_reads = 0
    sequence = next_read_sequence(file_stream)
    while n_reads < first_n and sequence is not None:
        adapter_counter.count_matches(sequence)
        n_reads += 1
        sequence = next_read_sequence(file_stream)
    return adapter_counter.result()


def count_adapters_in_fastq_file(file: Union[str, Path], *args, **kwargs) -> Dict[str, int]:
    with open_decompressed(file) as f:
        return count_adapters_in_fastq_stream(f, *args, **kwargs)


def report_adapter_counts(message: str, counts: Dict[str, int]) -> str:
    return "\n".join([message,
                     *list(map(lambda cnt: f"\t{cnt[0]} = {cnt[1]}",
                               counts.items()))])


if __name__ == "__main__":
    opts = parse_parameters(sys.argv)

    r12_adapters = list(map(lambda opt: parse_simple_adapter(opt, ConstructEnd.both),
                            opts.adapters))
    r1_adapters = list(map(lambda opt: parse_simple_adapter(opt, ConstructEnd.one),
                           opts.adapters1))
    r2_adapters = list(map(lambda opt: parse_simple_adapter(opt, ConstructEnd.two),
                           opts.adapters2))
    p_adapters = list(map(lambda opt: parse_palindrome_adapters(opt),
                          opts.paired_adapters))

    with futures.ThreadPoolExecutor(max_workers=2) as executor:
        matches1F: Optional[Future[Dict[str, int]]] = None
        matches2F: Optional[Future[Dict[str, int]]] = None
        if opts.fastq1 is not None:
            print(f"Searching adapters in first {opts.first_n} reads of {opts.fastq1} ...", file=sys.stderr)
            matches1F = \
                executor.submit(count_adapters_in_fastq_file,
                                opts.fastq1,
                                r12_adapters + r1_adapters +
                                # For searching, we use the raw read, but for reporting the name
                                # with the prefix added.
                                list(map(lambda p: p.read1.copy(name=f"{PREFIX}{p.read1.name}"),
                                         p_adapters)),
                                first_n=opts.first_n,
                                truncate=opts.truncate)
        if opts.fastq2 is not None:
            print(f"Searching adapters in first {opts.first_n} reads of {opts.fastq2} ...", file=sys.stderr)
            matches2F = \
                executor.submit(count_adapters_in_fastq_file,
                                opts.fastq2,
                                r12_adapters + r2_adapters +
                                # For searching, we use the raw read, but for reporting the name
                                # with the prefix added.
                                list(map(lambda p: p.read2.copy(name=f"{PREFIX}{p.read2.name}"),
                                         p_adapters)),
                                first_n=opts.first_n,
                                truncate=opts.truncate)

        if matches1F is not None:
            print(report_adapter_counts(f"Matching read 1 adapters in {opts.fastq1}", matches1F.result()),
                  file=sys.stderr)
        if matches2F is not None:
            print(report_adapter_counts(f"Matching read 2 adapters in {opts.fastq2}", matches2F.result()),
                  file=sys.stderr)

    try:
        print("".join(list(map(lambda a: a.to_fasta(), r12_adapters))),
              end="")
        print("".join(list(map(lambda a: a.to_fasta(), r1_adapters))),
              end="")
        print("".join(list(map(lambda a: a.to_fasta(), r2_adapters))),
              end="")
        print("".join(list(map(lambda p: p.to_palindrome_fasta(), p_adapters))),
              end="")

    except ValueError as ex:
        print(str(ex), file=sys.stderr)
        sys.exit(1)
