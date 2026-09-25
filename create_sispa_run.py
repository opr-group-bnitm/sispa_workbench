#!/usr/bin/env python3
"""Mix single-organism FASTQs into one artificial sequencing run.

The composition is either a CSV/TSV file with the columns organism_id,n_reads
or a list of "organism_id n_reads" pairs. n_reads is a number of reads drawn
at random (without replacement) from that organism's FASTQ, or "all". The
organism ids are the ones in data/organisms.tsv (see sync_reads.py).

Reads of the different organisms are interleaved in random order, as they
would be in a real run. The result goes to output/fastqs/<name>.fastq.gz. When
the composition was given as a list it is also saved as
output/compositions/<name>.csv, with "all" replaced by the actual read count,
so the run can be recreated with the CSV and the printed --seed.

The name defaults to the CSV's file name, or for a list to e.g.
COVID-5000_human-7000_DENV2-80000_LASV-all.

Examples:
    ./create_sispa_run.py COVID 5000, human 7000, DENV2 80000, LASV all
    ./create_sispa_run.py "COVID 5000, human all" --name covid_in_human --seed 42
    ./create_sispa_run.py my_mix.csv
"""

from __future__ import annotations

import argparse
import csv
import difflib
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Tuple

from workbench.fastq import FastqError, Record, open_write, read_fastq, write_record
from workbench.organisms import Organism, SyncError, sync_reads
from workbench.paths import default_data_dir, default_output_dir

Composition = List[Tuple[str, Optional[int]]]  # (organism_id, n_reads or None for all)


class RunError(RuntimeError):
    pass


def parse_amount(token: str, organism_id: str) -> Optional[int]:
    if token.strip().lower() == "all":
        return None
    try:
        n = int(token)
    except ValueError:
        n = -1
    if n <= 0:
        raise RunError(f"n_reads for {organism_id} must be a positive number or 'all', got {token!r}")
    return n


def parse_list(items: Sequence[str]) -> Composition:
    tokens = " ".join(items).replace(",", " ").split()
    if not tokens or len(tokens) % 2:
        raise RunError("expected pairs of 'organism_id n_reads', e.g.: COVID 5000, human 7000, LASV all")
    return [(tokens[i], parse_amount(tokens[i + 1], tokens[i])) for i in range(0, len(tokens), 2)]


def parse_table(path: Path) -> Composition:
    delimiter = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    with open(path, newline="") as fh:
        lines = [line for line in fh if line.strip() and not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines, delimiter=delimiter)
    columns = [c.strip() for c in reader.fieldnames or []]
    if "organism_id" not in columns or "n_reads" not in columns:
        raise RunError(f"{path} needs the columns organism_id and n_reads (found: {', '.join(columns) or 'none'})")
    composition = []
    for row in reader:
        row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
        composition.append((row["organism_id"], parse_amount(row["n_reads"], row["organism_id"])))
    if not composition:
        raise RunError(f"{path} lists no organisms")
    return composition


def default_name(composition: Composition) -> str:
    return "_".join(f"{oid}-{'all' if n is None else n}" for oid, n in composition)


@dataclass
class Source:
    organism: Organism
    n_reads: int
    path: Path


def resolve(composition: Composition, organisms: List[Organism], data_dir: Path) -> List[Source]:
    by_id = {o.organism_id: o for o in organisms}
    seen = set()
    sources = []
    for organism_id, n in composition:
        if organism_id in seen:
            raise RunError(f"{organism_id} is listed more than once")
        seen.add(organism_id)
        organism = by_id.get(organism_id)
        if organism is None:
            close = difflib.get_close_matches(organism_id, by_id, n=3)
            hint = f" -- did you mean {', '.join(close)}?" if close else ""
            known = ", ".join(sorted(by_id)) or "none yet; add FASTQs and run sync_reads.py"
            raise RunError(f"unknown organism_id {organism_id!r}{hint}\n  known ids: {known}")
        if n is None:
            n = organism.n_reads
        elif n > organism.n_reads:
            raise RunError(
                f"{organism_id}: {n} reads requested but {organism.filename} holds only {organism.n_reads}"
            )
        if n == 0:
            raise RunError(f"{organism_id}: {organism.filename} holds no reads")
        sources.append(Source(organism, n, organism.path(data_dir)))
    return sources


def draw(source: Source, wanted: Optional[List[int]]) -> Iterator[Record]:
    """Yield the records at the sorted indices in wanted, or every record if wanted is None."""
    stale = RunError(f"{source.path} changed while it was being read; run again")
    if wanted is None:
        count = 0
        for record in read_fastq(source.path):
            count += 1
            if count > source.n_reads:
                raise stale
            yield record
        if count < source.n_reads:
            raise stale
        return

    j = 0
    for i, record in enumerate(read_fastq(source.path)):
        if i == wanted[j]:
            yield record
            j += 1
            if j == len(wanted):
                return
    raise stale


def mix(sources: List[Source], out_path: Path, rng: random.Random) -> None:
    """Interleave the sources in random order: each next read comes from a
    source with probability proportional to the reads it still has to give."""
    streams = []
    for source in sources:
        total_in_file = source.organism.n_reads
        # a uniform random subset of read indices, or None to take every read
        wanted = None if source.n_reads == total_in_file else sorted(rng.sample(range(total_in_file), source.n_reads))
        streams.append(draw(source, wanted))
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


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("composition", nargs="+",
                        help="a CSV/TSV file with organism_id,n_reads, or a list like: COVID 5000, human all")
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
    except (RunError, SyncError, FastqError) as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1


def run(args) -> int:
    table: Optional[Path] = None
    first = Path(args.composition[0])
    if len(args.composition) == 1 and (first.is_file() or first.suffix.lower() in (".csv", ".tsv", ".tab")):
        if not first.is_file():
            raise RunError(f"composition file not found: {first}")
        table = first
        composition = parse_table(table)
    else:
        composition = parse_list(args.composition)

    name = args.name or (table.stem if table else default_name(composition))
    if not name or "/" in name or name.startswith("."):
        raise RunError(f"invalid run name {name!r}")
    fastq_path = args.fastq_dir / f"{name}.fastq.gz"
    composition_path = None if table else args.composition_dir / f"{name}.csv"
    for path in (fastq_path, composition_path):
        if path and path.exists() and not args.force:
            raise RunError(f"{path} already exists (use --force to overwrite, or pick another --name)")

    # make sure the read counts in organisms.tsv describe the files as they are now
    organisms = sync_reads(args.data_dir, log=lambda msg: print(f"[sync] {msg}"))
    sources = resolve(composition, organisms, args.data_dir)

    seed = args.seed if args.seed is not None else random.SystemRandom().randrange(2**32)
    total = sum(s.n_reads for s in sources)
    print(f"[mix]  {name}: {total} reads, seed {seed}")
    for s in sources:
        print(f"       {s.organism.organism_id:<24} {s.n_reads:>10} reads  {100 * s.n_reads / total:6.2f}%"
              f"  from {s.organism.filename}")

    args.fastq_dir.mkdir(parents=True, exist_ok=True)
    mix(sources, fastq_path, random.Random(seed))
    print(f"[done] {fastq_path}")

    if composition_path:
        args.composition_dir.mkdir(parents=True, exist_ok=True)
        with open(composition_path, "w", newline="") as fh:
            writer = csv.writer(fh, lineterminator="\n")
            writer.writerow(["organism_id", "n_reads"])
            writer.writerows([s.organism.organism_id, s.n_reads] for s in sources)
        print(f"[done] {composition_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
