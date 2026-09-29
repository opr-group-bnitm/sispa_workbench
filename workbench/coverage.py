"""Breadth of coverage and depth range of alignments on a reference.

Depth follows vimop's ``samtools depth -aa -J``: the depth at a reference
position is the number of alignments (primary and supplementary) spanning it,
deletions included, and every position of every reference sequence counts,
covered or not. Only alignment start and end points are kept, so memory grows
with the number of alignments, not with the reference length.
"""

from __future__ import annotations

import random
from array import array
from bisect import bisect_left
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from .alignments import ReadAlignments


@dataclass
class DepthSummary:
    length: int
    covered: int  # positions with depth >= 1
    min_depth: int
    max_depth: int
    depth_sum: int = 0  # the depths of all positions added up

    @property
    def breadth(self) -> float:
        """Percentage of positions with depth >= 1."""
        return 100 * self.covered / self.length if self.length else 0.0

    @property
    def mean_depth(self) -> float:
        """Mean depth over all positions, covered or not."""
        return self.depth_sum / self.length if self.length else 0.0


def summarize(length: int, starts: Sequence[int], ends: Sequence[int]) -> DepthSummary:
    """Depth summary over [0, length) for the half-open intervals [starts[i], ends[i])."""
    starts, ends = sorted(starts), sorted(ends)
    n = len(starts)
    covered = max_depth = depth = depth_sum = prev = 0
    min_depth = None
    i = j = 0
    while j < n:
        pos = min(starts[i], ends[j]) if i < n else ends[j]
        if pos > prev:  # [prev, pos) has constant depth
            if depth:
                covered += pos - prev
                depth_sum += depth * (pos - prev)
            min_depth = depth if min_depth is None else min(min_depth, depth)
            max_depth = max(max_depth, depth)
            prev = pos
        while i < n and starts[i] == pos:
            depth += 1
            i += 1
        while j < n and ends[j] == pos:
            depth -= 1
            j += 1
    if prev < length:  # uncovered tail, or the whole sequence without alignments
        min_depth = 0
    return DepthSummary(length, covered, min_depth or 0, max_depth, depth_sum)


class DepthCounter:
    """Collects alignment spans per reference sequence."""

    def __init__(self) -> None:
        self.starts: Dict[str, array] = {}
        self.ends: Dict[str, array] = {}

    def add(self, sequence: str, start: int, end: int) -> None:
        if sequence not in self.starts:
            self.starts[sequence] = array("q")
            self.ends[sequence] = array("q")
        self.starts[sequence].append(start)
        self.ends[sequence].append(end)

    def summarize(self, lengths: Dict[str, int]) -> Dict[str, DepthSummary]:
        """One summary per reference sequence, in the order of lengths."""
        return {
            name: summarize(length, self.starts.get(name, ()), self.ends.get(name, ()))
            for name, length in lengths.items()
        }


def combine(summaries: List[DepthSummary]) -> DepthSummary:
    """The summary over all reference sequences together."""
    if not summaries:
        return DepthSummary(0, 0, 0, 0)
    return DepthSummary(
        sum(s.length for s in summaries),
        sum(s.covered for s in summaries),
        min(s.min_depth for s in summaries),
        max(s.max_depth for s in summaries),
        sum(s.depth_sum for s in summaries),
    )


def depth_of_reads(alignments: ReadAlignments, reads: Optional[Sequence[int]] = None) -> DepthSummary:
    """Depth over the whole reference of the given reads (positions in the
    FASTQ), or of all reads."""
    keep = None if reads is None else set(reads)
    names = alignments.sequences
    counter = DepthCounter()
    for read, k, start, end in zip(alignments.read, alignments.sequence, alignments.start, alignments.end):
        if keep is None or read in keep:
            counter.add(names[k], start, end)
    return combine(list(counter.summarize(alignments.lengths).values()))


@dataclass
class DepthSelection:
    reads: List[int]  # positions in the FASTQ of the chosen reads, ascending
    min_depth: int
    full_bases: int  # reference positions where min_depth is reached
    short_bases: int  # covered positions where all reads together stay below min_depth
    total_bases: int  # reference length

    @property
    def uncovered_bases(self) -> int:
        return self.total_bases - self.full_bases - self.short_bases


def select_for_depth(alignments: ReadAlignments, min_depth: int, rng: random.Random) -> DepthSelection:
    """Reads, taken in random order, that bring every reference position to
    min_depth, or to the depth all reads together reach there if that is lower.

    A read is taken only if it covers a position still below its target, so
    reads over positions that are already deep enough are skipped, and taking
    stops once every position is served. Depth counts alignments as in
    summarize(). The reference is cut into the stretches between alignment
    ends, which have a constant depth, so memory and time grow with the number
    of alignments, not with the reference length."""
    al = alignments
    # the stretches: between consecutive alignment ends of each reference sequence
    points = [set() for _ in al.lengths]
    for k, start, end in zip(al.sequence, al.start, al.end):
        points[k].update((start, end))
    breaks = [sorted(p) for p in points]
    offsets, n_stretches = [], 0
    for b in breaks:
        offsets.append(n_stretches)
        n_stretches += max(len(b) - 1, 0)

    # every alignment as the range [first, last) of stretches it covers
    first, last = array("q"), array("q")
    change = [0] * (n_stretches + 1)
    for k, start, end in zip(al.sequence, al.start, al.end):
        b, offset = breaks[k], offsets[k]
        first.append(offset + bisect_left(b, start))
        last.append(offset + bisect_left(b, end))
        change[first[-1]] += 1
        change[last[-1]] -= 1

    # need: how many more alignments each stretch wants
    need = [0] * n_stretches
    full = short = depth = i = 0
    for b in breaks:
        for left, right in zip(b, b[1:]):
            depth += change[i]
            need[i] = min(min_depth, depth)
            if depth >= min_depth:
                full += right - left
            elif depth:
                short += right - left
            i += 1

    # next_open[i] leads to the first stretch at or after i that still needs
    # alignments (n_stretches if none); served stretches are skipped for good
    next_open = list(range(n_stretches + 1))
    open_stretches = 0
    for i in range(n_stretches):
        if need[i]:
            open_stretches += 1
        else:
            next_open[i] = i + 1

    def find(i: int) -> int:
        root = i
        while next_open[root] != root:
            root = next_open[root]
        while i != root:
            next_open[i], i = root, next_open[i]
        return root

    rows_by_read: Dict[int, List[int]] = {}
    for row, read in enumerate(al.read):
        rows_by_read.setdefault(read, []).append(row)
    reads = list(rows_by_read.items())
    rng.shuffle(reads)

    chosen = []
    for read, rows in reads:
        if not open_stretches:
            break
        if not any(find(first[row]) < last[row] for row in rows):
            continue  # every stretch it covers is served already
        chosen.append(read)
        for row in rows:
            i = find(first[row])
            while i < last[row]:
                need[i] -= 1
                if not need[i]:
                    next_open[i] = i + 1
                    open_stretches -= 1
                i = find(i + 1)
    return DepthSelection(sorted(chosen), min_depth, full, short, sum(al.lengths.values()))
