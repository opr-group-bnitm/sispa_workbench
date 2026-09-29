#!/usr/bin/env python3
"""Mix single-organism FASTQs into one artificial sequencing run.

The composition says, for each organism (the ids in data/organisms.tsv, made
with add_from_ref.py), which of its reads go into the run:

    n_reads        that many reads, drawn at random without replacement, or "all"
    minimum_depth  just enough reads, drawn at random, to cover every position
                   of the organism's reference at least this deep; where all its
                   reads together are not that deep, every read there is used.
                   "all" takes all reads, as for n_reads

It is either a CSV/TSV file with the column organism_id and the columns
n_reads and/or minimum_depth (a row with both asks for the depth), or a list
such as "COVID 5000, human all, LASV 20x", where 20x asks for a minimum depth
of 20. data/input_compositions/template.csv is a template for the file.

Reads of the different organisms are interleaved in random order, as they
would be in a real run. The result goes to output/fastqs/<name>.fastq.gz, and
what it contains to output/compositions/<name>.csv: per organism, the number
of reads it got and the min_depth, max_depth and mean_depth these reach on its
reference, counted like in organisms.tsv. The same composition and the
printed --seed recreate the run exactly.

The name defaults to the CSV's file name, or for a list to e.g.
COVID-5000_human-7000_DENV2-80000_LASV-all.

Examples:
    ./create_sispa_run.py COVID 5000, human 7000, DENV2 80000, LASV all
    ./create_sispa_run.py LASV 20x, human 50000 --seed 42
    ./create_sispa_run.py "COVID 5000, human all" --name covid_in_human
    ./create_sispa_run.py my_mix.csv
"""

from __future__ import annotations

import argparse
import csv
import difflib
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Sequence

from workbench.alignments import AlignmentsError, alignments_path, read_alignments
from workbench.coverage import DepthSelection, DepthSummary, depth_of_reads, select_for_depth
from workbench.fastq import FastqError, Record, open_write, read_fastq, strip_fastq_suffix, write_record
from workbench.organisms import Organism, RegistryError, load_organisms, unregistered_fastqs
from workbench.paths import default_data_dir, default_output_dir


# what output/compositions/<name>.csv records: what the run contains
COMPOSITION_COLUMNS = ["organism_id", "n_reads", "min_depth", "max_depth", "mean_depth"]


class RunError(RuntimeError):
    pass


@dataclass
class Request:
    organism_id: str
    n_reads: Optional[int] = None  # with no min_depth either: all reads
    min_depth: Optional[int] = None

    @property
    def label(self) -> str:
        if self.min_depth is not None:
            return f"{self.min_depth}x"
        return "all" if self.n_reads is None else str(self.n_reads)


DEPTH = re.compile(r"^(\d+)[xX]?$")


def positive(token: str) -> int:
    try:
        return max(int(token), 0)
    except ValueError:
        return 0


def parse_amount(token: str, organism_id: str) -> Request:
    token = token.strip()
    if token.lower() == "all":
        return Request(organism_id)
    if token[-1:] in ("x", "X") and positive(token[:-1]):
        return Request(organism_id, min_depth=positive(token[:-1]))
    if positive(token):
        return Request(organism_id, n_reads=positive(token))
    raise RunError(f"n_reads for {organism_id} must be a positive number, 'all' or a minimum depth "
                   f"such as 20x, got {token!r}")


def parse_list(items: Sequence[str]) -> List[Request]:
    tokens = " ".join(items).replace(",", " ").split()
    if not tokens or len(tokens) % 2:
        raise RunError("expected pairs of 'organism_id amount', e.g.: COVID 5000, human all, LASV 20x")
    return [parse_amount(tokens[i + 1], tokens[i]) for i in range(0, len(tokens), 2)]


def parse_table(path: Path) -> List[Request]:
    delimiter = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    with open(path, newline="") as fh:
        lines = [line for line in fh if line.strip() and not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines, delimiter=delimiter)
    columns = [c.strip() for c in reader.fieldnames or []]
    if "organism_id" not in columns or not {"n_reads", "minimum_depth"} & set(columns):
        raise RunError(f"{path} needs the column organism_id and the column n_reads and/or minimum_depth "
                       f"(found: {', '.join(columns) or 'none'})")
    requests = []
    for row in reader:
        row = {(k or "").strip(): (v or "").strip() for k, v in row.items() if k is not None}
        organism_id = row["organism_id"]
        depth = row.get("minimum_depth", "")
        if depth.lower() == "all":
            requests.append(Request(organism_id))
        elif depth:
            match = DEPTH.match(depth)
            if not match or not positive(match.group(1)):
                raise RunError(f"minimum_depth for {organism_id} must be a positive number or 'all', got {depth!r}")
            requests.append(Request(organism_id, min_depth=positive(match.group(1))))
        elif row.get("n_reads"):
            requests.append(parse_amount(row["n_reads"], organism_id))
        else:
            raise RunError(f"{path}: {organism_id} has neither n_reads nor minimum_depth")
    if not requests:
        raise RunError(f"{path} lists no organisms")
    return requests


def default_name(requests: List[Request]) -> str:
    return "_".join(f"{r.organism_id}-{r.label}" for r in requests)


@dataclass
class Source:
    organism: Organism
    request: Request
    path: Path
    alignments: Path
    n_reads: int = 0  # set by plan()
    wanted: Optional[List[int]] = None  # positions of the reads to take, ascending; None for all
    depth: Optional[DepthSelection] = None  # for a minimum depth: how the reads were picked
    coverage: Optional[DepthSummary] = None  # depth of the reads taken


def resolve(requests: List[Request], organisms: List[Organism], data_dir: Path) -> List[Source]:
    by_id = {o.organism_id: o for o in organisms}
    seen = set()
    sources = []
    for request in requests:
        organism_id = request.organism_id
        if organism_id in seen:
            raise RunError(f"{organism_id} is listed more than once")
        seen.add(organism_id)
        organism = by_id.get(organism_id)
        if organism is None:
            close = difflib.get_close_matches(organism_id, by_id, n=3)
            hint = f" -- did you mean {', '.join(close)}?" if close else ""
            known = ", ".join(sorted(by_id)) or "none yet; extract reads with add_from_ref.py"
            stray = [f for f in unregistered_fastqs(data_dir, organisms)
                     if strip_fastq_suffix(Path(f).name) == organism_id]
            if stray:
                hint += (f"\n  {stray[0]} was not created by add_from_ref.py: put its reads in "
                         "data/raw_data/<dataset_id>/ and extract them with add_from_ref.py")
            raise RunError(f"unknown organism_id {organism_id!r}{hint}\n  known ids: {known}")
        path = organism.path(data_dir)
        if not path.is_file():
            raise RunError(f"{organism_id}: {organism.filename} is missing; extract it again with add_from_ref.py")
        if organism.n_reads == 0:
            raise RunError(f"{organism_id}: {organism.filename} holds no reads")
        alignments = alignments_path(data_dir, organism.filename)
        if not alignments.is_file():
            raise RunError(f"{organism_id}: {organism.filename} has no alignment file ({alignments}), which "
                           "is needed for its depth; extract it again with add_from_ref.py --force")
        if request.n_reads is not None and request.n_reads > organism.n_reads:
            raise RunError(
                f"{organism_id}: {request.n_reads} reads requested but {organism.filename} "
                f"holds only {organism.n_reads}"
            )
        sources.append(Source(organism, request, path, alignments))
    return sources


def stale(source: Source) -> RunError:
    return RunError(f"{source.path} no longer holds the {source.organism.n_reads} reads organisms.tsv "
                    "lists for it; extract it again with add_from_ref.py --force")


def plan(sources: List[Source], rng: random.Random) -> None:
    """Decide which reads of each source go into the run, and how deep they cover its reference."""
    for source in sources:
        total = source.organism.n_reads
        alignments = read_alignments(source.alignments)
        if alignments.n_reads != total or (alignments.read and max(alignments.read) >= total):
            raise stale(source)
        if source.request.min_depth is not None:
            source.depth = select_for_depth(alignments, source.request.min_depth, rng)
            source.n_reads = len(source.depth.reads)
            source.wanted = None if source.n_reads == total else source.depth.reads
        else:
            source.n_reads = total if source.request.n_reads is None else source.request.n_reads
            # a uniform random subset of read positions, or None to take every read
            source.wanted = None if source.n_reads == total else sorted(rng.sample(range(total), source.n_reads))
        source.coverage = depth_of_reads(alignments, source.wanted)


def draw(source: Source) -> Iterator[Record]:
    """Yield the records at the positions in source.wanted, or every record."""
    if source.wanted is None:
        count = 0
        for record in read_fastq(source.path):
            count += 1
            if count > source.n_reads:
                raise stale(source)
            yield record
        if count < source.n_reads:
            raise stale(source)
        return

    wanted, j = source.wanted, 0
    for i, record in enumerate(read_fastq(source.path)):
        if i == wanted[j]:
            yield record
            j += 1
            if j == len(wanted):
                return
    raise stale(source)


def mix(sources: List[Source], out_path: Path, rng: random.Random) -> None:
    """Interleave the sources in random order: each next read comes from a
    source with probability proportional to the reads it still has to give."""
    streams = [draw(source) for source in sources]
    remaining = [source.n_reads for source in sources]
    total = sum(remaining)
    tmp_path = out_path.with_name(f".partial.{out_path.name}")  # keeps the .gz suffix
    try:
        with open_write(tmp_path) as out:
            while total:
                pick = rng.randrange(total)
                for i, left in enumerate(remaining):
                    if pick < left:
                        break
                    pick -= left
                write_record(out, next(streams[i]))
                remaining[i] -= 1
                total -= 1
            for stream in streams:  # lets the "all" streams check the file did not grow
                for _ in stream:
                    pass
        tmp_path.replace(out_path)
    finally:
        tmp_path.unlink(missing_ok=True)


def depth_report(source: Source) -> str:
    d = source.depth

    def pct(bases: int) -> str:
        return f"{100 * bases / d.total_bases if d.total_bases else 0:.2f}%"

    return (f"[depth] {source.organism.organism_id}: {source.n_reads} of {source.organism.n_reads} reads for "
            f"{d.min_depth}x; {pct(d.full_bases)} of the reference at {d.min_depth}x or more, "
            f"{pct(d.short_bases)} below with every read there, {pct(d.uncovered_bases)} not covered")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("composition", nargs="+",
                        help="a CSV/TSV file with organism_id and n_reads and/or minimum_depth, "
                             "or a list like: COVID 5000, human all, LASV 20x")
    parser.add_argument("--name", help="run name (default: CSV file name, or built from the list)")
    parser.add_argument("--seed", type=int, help="random seed, for a reproducible run (default: random, printed)")
    parser.add_argument("--fastq-dir", type=Path, default=default_output_dir() / "fastqs",
                        help="where the run FASTQ is written (default: %(default)s)")
    parser.add_argument("--composition-dir", type=Path, default=default_output_dir() / "compositions",
                        help="where the composition CSV is written (default: %(default)s)")
    parser.add_argument("--force", action="store_true", help="overwrite an existing run of the same name")
    parser.add_argument("--data-dir", type=Path, default=default_data_dir(),
                        help="workbench data directory (default: %(default)s)")
    args = parser.parse_args(argv)

    try:
        return run(args)
    except (RunError, RegistryError, FastqError, AlignmentsError) as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1


def run(args) -> int:
    table: Optional[Path] = None
    first = Path(args.composition[0])
    if len(args.composition) == 1 and (first.is_file() or first.suffix.lower() in (".csv", ".tsv", ".tab")):
        if not first.is_file():
            raise RunError(f"composition file not found: {first}")
        table = first
        requests = parse_table(table)
    else:
        requests = parse_list(args.composition)

    name = args.name or (table.stem if table else default_name(requests))
    if not name or "/" in name or name.startswith("."):
        raise RunError(f"invalid run name {name!r}")
    fastq_path = args.fastq_dir / f"{name}.fastq.gz"
    composition_path = args.composition_dir / f"{name}.csv"
    for path in (fastq_path, composition_path):
        if path.exists() and not args.force:
            raise RunError(f"{path} already exists (use --force to overwrite, or pick another --name)")

    sources = resolve(requests, load_organisms(args.data_dir), args.data_dir)
    seed = args.seed if args.seed is not None else random.SystemRandom().randrange(2**32)
    rng = random.Random(seed)
    plan(sources, rng)

    total = sum(s.n_reads for s in sources)
    print(f"[mix]  {name}: {total} reads, seed {seed}")
    for s in sources:
        c = s.coverage
        print(f"       {s.organism.organism_id:<24} {s.n_reads:>10} reads  {100 * s.n_reads / total:6.2f}%"
              f"  depth {c.min_depth}-{c.max_depth}, mean {c.mean_depth:.2f}  from {s.organism.filename}")
    for s in sources:
        if s.depth:
            print(depth_report(s))

    args.fastq_dir.mkdir(parents=True, exist_ok=True)
    mix(sources, fastq_path, rng)
    print(f"[done] {fastq_path}")

    args.composition_dir.mkdir(parents=True, exist_ok=True)
    with open(composition_path, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(COMPOSITION_COLUMNS)
        for s in sources:
            writer.writerow([s.organism.organism_id, s.n_reads,
                             s.coverage.min_depth, s.coverage.max_depth, f"{s.coverage.mean_depth:.2f}"])
    print(f"[done] {composition_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
