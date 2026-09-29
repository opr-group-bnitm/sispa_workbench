"""Mapping reads with minimap2 the way the vimop pipeline does.

The settings follow vimop (github.com/opr-group-bnitm/vimop, lib/processes.nf
and nextflow.config), and reads are selected from minimap2's SAM output.
"""

from __future__ import annotations

import re
import shlex
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .coverage import DepthCounter


class MappingError(RuntimeError):
    pass


# minimap2 settings of the vimop pipeline (github.com/opr-group-bnitm/vimop,
# lib/processes.nf and nextflow.config)
VIMOP_ARGS = ["--secondary=no"]
# map_to_ref: reads against a virus reference (map_to_target_minimap_kmer_size / _window_size)
VIMOP_VIRUS_ARGS = ["-k", "11", "-w", "5"]
# filter_virus_target: references above target_filter_minimap_index_split_thresh
# get a split index of at most 2G bases per part
VIMOP_SPLIT_THRESHOLD = 500 * 1024 * 1024
VIMOP_SPLIT_ARGS = ["-I", "2G"]

def minimap2_args(reference: Path, category: str, split_prefix: Path) -> List[str]:
    """vimop's settings for this reference: -k 11 -w 5 for viruses, like its
    map_to_ref, and plain map-ont for hosts and bacteria, like its
    filter_contaminants; large FASTA references are indexed in 2G parts."""
    args = list(VIMOP_ARGS)
    if category == "virus":
        args += VIMOP_VIRUS_ARGS
    if not reference.name.endswith(".mmi") and reference.stat().st_size > VIMOP_SPLIT_THRESHOLD:
        args += VIMOP_SPLIT_ARGS + ["--split-prefix", str(split_prefix)]
    return args


def covered_length(intervals: List[Tuple[int, int]]) -> int:
    covered, end = 0, -1
    for start, stop in sorted(intervals):
        if stop <= end:
            continue
        covered += stop - max(start, end)
        end = stop
    return covered


REF_OPS = re.compile(rb"(\d+)[MDN=X]")      # CIGAR operations that consume the reference
QUERY_OPS = re.compile(rb"(\d+)[MIS=XH]")   # ... that consume the read, clips included
LEAD_CLIP = re.compile(rb"^(\d+)[SH]")
TRAIL_CLIP = re.compile(rb"(\d+)[SH]$")


@dataclass
class Alignment:
    sequence: str
    start: int  # 0-based, on the reference
    end: int
    query_start: int  # on the read as sequenced
    query_end: int
    query_length: int


def parse_alignment(fields: List[bytes]) -> Alignment:
    """An alignment from the first six SAM fields of a mapped record."""
    _, flag, rname, pos, _, cigar = fields
    start = int(pos) - 1
    end = start + sum(map(int, REF_OPS.findall(cigar)))
    query_length = sum(map(int, QUERY_OPS.findall(cigar)))
    lead = LEAD_CLIP.match(cigar)
    trail = TRAIL_CLIP.search(cigar)
    lead_clip = int(lead.group(1)) if lead else 0
    trail_clip = int(trail.group(1)) if trail else 0
    if int(flag) & 0x10:  # reverse strand: the CIGAR runs against the read
        lead_clip, trail_clip = trail_clip, lead_clip
    return Alignment(rname.decode(), start, end, lead_clip, query_length - trail_clip, query_length)


class Mapping:
    """Maps FASTQs to one reference and collects the reads that pass the filters."""

    def __init__(self, reference: Path, preset: str, threads: int, min_mapq: int,
                 min_fraction: float, extra_args: List[str]) -> None:
        self.reference = reference
        self.preset = preset
        self.threads = threads
        self.min_mapq = min_mapq
        self.min_fraction = min_fraction
        self.extra_args = extra_args
        self.lengths: Dict[str, int] = {}  # reference sequences, from the SAM header
        self.depth = DepthCounter()

    def keep(self, alignments: List[Alignment]) -> bool:
        if self.min_fraction <= 0:
            return True
        length = alignments[0].query_length
        spans = [(a.query_start, a.query_end) for a in alignments]
        return length == 0 or covered_length(spans) / length >= self.min_fraction

    def read_names(self, fastq: Path) -> Dict[str, List[Tuple[str, int, int]]]:
        """The reads in fastq to keep, by name, with their alignments as
        (sequence, start, end); these count towards the depth."""
        cmd = ["minimap2", "-a", "-x", self.preset, "--sam-hit-only", "-t", str(self.threads),
               *self.extra_args, str(self.reference), str(fastq)]
        print(f"[map]     {' '.join(shlex.quote(c) for c in cmd)}")
        names: Dict[str, List[Tuple[str, int, int]]] = {}

        def flush(qname: Optional[bytes], alignments: List[Alignment]) -> None:
            if alignments and self.keep(alignments):
                names[qname.decode()] = [(a.sequence, a.start, a.end) for a in alignments]
                for a in alignments:
                    self.depth.add(a.sequence, a.start, a.end)

        with tempfile.TemporaryFile() as err:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, bufsize=1 << 20)
            current: Optional[bytes] = None
            group: List[Alignment] = []
            for line in proc.stdout:
                if line.startswith(b"@"):
                    if line.startswith(b"@SQ"):
                        tags = dict(f.split(b":", 1) for f in line.rstrip(b"\r\n").split(b"\t")[1:] if b":" in f)
                        self.lengths[tags[b"SN"].decode()] = int(tags[b"LN"])
                    continue
                fields = line.split(b"\t", 6)[:6]
                # unmapped or secondary alignments, or below the mapping quality cut-off
                if int(fields[1]) & 0x104 or int(fields[4]) < self.min_mapq:
                    continue
                # minimap2 writes all alignments of a read one after the other
                if fields[0] != current:
                    flush(current, group)
                    current, group = fields[0], []
                group.append(parse_alignment(fields))
            flush(current, group)
            if proc.wait() != 0:
                err.seek(0)
                tail = err.read().decode(errors="replace").strip().splitlines()[-15:]
                raise MappingError(f"minimap2 failed on {fastq} (exit code {proc.returncode}):\n  " + "\n  ".join(tail))
        return names
