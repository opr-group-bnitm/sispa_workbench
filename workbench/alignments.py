"""Where each read of a single-organism FASTQ aligns on its reference.

For every FASTQ it makes, add_from_ref.py writes an alignment file to
data/alignments/, in the same layout as the read folders
(data/alignments/virus_reads/COVID.fastq.gz.alignments.tsv.gz for
data/virus_reads/COVID.fastq.gz):

    #sequence  <name>  <length>      once per reference sequence
    read  sequence  start  end       column names
    0     LASV_S    0      1532      one row per alignment, primary and supplementary
    ...
    #reads  <n>                      number of reads in the FASTQ

``read`` is the 0-based position of the read in the FASTQ; start and end are
0-based, end exclusive. create_sispa_run.py uses it to pick reads for a
minimum depth without mapping again.
"""

from __future__ import annotations

import gzip
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from . import paths
from .fastq import open_write

SUFFIX = ".alignments.tsv.gz"


class AlignmentsError(RuntimeError):
    pass


def alignments_path(data_dir: Path, filename: str) -> Path:
    """The alignment file of a FASTQ, given as in organisms.tsv (virus_reads/COVID.fastq.gz)."""
    return paths.alignments_dir(data_dir) / f"{filename}{SUFFIX}"


class AlignmentWriter:
    def __init__(self, path: Path, lengths: Dict[str, int]) -> None:
        self.out = open_write(path)
        for name, length in lengths.items():
            self.out.write(f"#sequence\t{name}\t{length}\n".encode())
        self.out.write(b"read\tsequence\tstart\tend\n")

    def add(self, read: int, alignments: List[Tuple[str, int, int]]) -> None:
        for sequence, start, end in alignments:
            self.out.write(f"{read}\t{sequence}\t{start}\t{end}\n".encode())

    def close(self, n_reads: int) -> None:
        self.out.write(f"#reads\t{n_reads}\n".encode())
        self.out.close()


@dataclass
class ReadAlignments:
    lengths: Dict[str, int]  # reference sequences
    n_reads: int  # reads in the FASTQ
    read: array  # per alignment: the read's position in the FASTQ
    sequence: array  # per alignment: index into sequences
    start: array
    end: array

    @property
    def sequences(self) -> List[str]:
        return list(self.lengths)


def read_alignments(path: Path) -> ReadAlignments:
    lengths: Dict[str, int] = {}
    index: Dict[str, int] = {}
    read, sequence, start, end = array("q"), array("q"), array("q"), array("q")
    n_reads = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "#sequence":
                index[fields[1]] = len(lengths)
                lengths[fields[1]] = int(fields[2])
            elif fields[0] == "#reads":
                n_reads = int(fields[1])
            elif fields[0] != "read":
                if fields[1] not in index:
                    raise AlignmentsError(f"{path}: sequence {fields[1]!r} is not declared")
                read.append(int(fields[0]))
                sequence.append(index[fields[1]])
                start.append(int(fields[2]))
                end.append(int(fields[3]))
    if n_reads is None:
        raise AlignmentsError(f"{path} is incomplete")
    return ReadAlignments(lengths, n_reads, read, sequence, start, end)
