"""Minimal, dependency free FASTQ reading, writing and statistics."""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterator, Tuple

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

# (header without '@', sequence, quality), all without line endings
Record = Tuple[bytes, bytes, bytes]


class FastqError(ValueError):
    pass


def is_fastq(path: Path | str) -> bool:
    return str(path).lower().endswith(FASTQ_SUFFIXES)


def strip_fastq_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in FASTQ_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _is_gzip(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_read(path: Path | str) -> BinaryIO:
    path = Path(path)
    if _is_gzip(path):
        return io.BufferedReader(gzip.GzipFile(path, "rb"), buffer_size=1 << 20)
    return open(path, "rb", buffering=1 << 20)


def open_write(path: Path | str) -> BinaryIO:
    """Open for writing, gzip compressed when the name ends in .gz."""
    path = Path(path)
    if path.name.lower().endswith(".gz"):
        return io.BufferedWriter(gzip.GzipFile(path, "wb", compresslevel=4), buffer_size=1 << 20)
    return open(path, "wb", buffering=1 << 20)


def read_fastq(path: Path | str) -> Iterator[Record]:
    """Yield the records of a four-line FASTQ file, plain or gzipped."""
    with open_read(path) as fh:
        line_no = 0
        while True:
            header = fh.readline()
            line_no += 1
            if not header:
                return
            if not header.strip():
                continue  # tolerate blank lines between or after records
            seq = fh.readline().rstrip(b"\r\n")
            plus = fh.readline()
            qual = fh.readline().rstrip(b"\r\n")
            if not header.startswith(b"@") or not plus.startswith(b"+") or len(seq) != len(qual):
                raise FastqError(
                    f"{path}: malformed FASTQ record at line {line_no} "
                    "(expected four-line records: @header, sequence, +, quality)"
                )
            line_no += 3
            yield header[1:].rstrip(b"\r\n"), seq, qual


def write_record(out: BinaryIO, record: Record) -> None:
    header, seq, qual = record
    out.write(b"@" + header + b"\n" + seq + b"\n+\n" + qual + b"\n")


def read_name(header: bytes) -> str:
    """The read id as aligners report it: the header up to the first whitespace."""
    return header.split(None, 1)[0].decode()


@dataclass
class FastqStats:
    n_reads: int = 0
    min_read_length: int = 0
    max_read_length: int = 0
    total_bases: int = 0

    @property
    def avg_read_length(self) -> float:
        return self.total_bases / self.n_reads if self.n_reads else 0.0

    def add(self, length: int) -> None:
        if self.n_reads == 0 or length < self.min_read_length:
            self.min_read_length = length
        if length > self.max_read_length:
            self.max_read_length = length
        self.n_reads += 1
        self.total_bases += length


def fastq_stats(path: Path | str) -> FastqStats:
    stats = FastqStats()
    for _, seq, _ in read_fastq(path):
        stats.add(len(seq))
    return stats
