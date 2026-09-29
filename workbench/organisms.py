"""The organisms.tsv registry of single-organism FASTQ files.

organisms.tsv lists every FASTQ in data/{virus,bacteria,host}_reads/ with its
read statistics. ``filename`` is relative to the data directory, e.g.
``virus_reads/COVID.fastq.gz``. A file's organism_id defaults to its name
without the FASTQ suffix; an id given explicitly (by add_from_ref.py, or by
editing the table) is kept on later syncs.

``reference`` and ``source_fastq`` record where add_from_ref.py took the reads
from: the reference they were mapped to and the FASTQ(s) they were extracted
from (';'-separated), relative to the data directory when inside it. They are
empty for files added by hand, and kept on later syncs like explicit ids.

Read statistics are cached in data/.organisms_cache.tsv together with each
file's size and modification time, so a sync only reads files that are new or
have changed since the last one.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Tuple

from .fastq import FastqStats, fastq_stats, is_fastq, strip_fastq_suffix
from .paths import ORGANISMS_TSV, READ_CATEGORIES, SYNC_CACHE, reads_dir

COLUMNS = [
    "organism_id", "filename", "avg_read_length", "max_read_length", "min_read_length", "n_reads",
    "reference", "source_fastq",
]
CACHE_COLUMNS = ["filename", "size", "mtime_ns", "n_reads", "min_read_length", "max_read_length", "total_bases"]

# ids are used unquoted in composition lists ("COVID 5000, human all") and in
# run names, so keep them free of whitespace, commas and path separators
ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class SyncError(RuntimeError):
    pass


@dataclass
class Organism:
    organism_id: str
    filename: str
    avg_read_length: float
    max_read_length: int
    min_read_length: int
    n_reads: int
    reference: str = ""
    source_fastq: str = ""

    def path(self, data_dir: Path) -> Path:
        return data_dir / self.filename


def invalid_id_reason(organism_id: str) -> Optional[str]:
    if ID_PATTERN.match(organism_id):
        return None
    return (
        f"organism_id {organism_id!r} may only contain letters, digits, '.', '_' and '-' "
        "and must start with a letter or digit"
    )


def load_organisms(data_dir: Path) -> List[Organism]:
    path = data_dir / ORGANISMS_TSV
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return [
            Organism(
                organism_id=row["organism_id"],
                filename=row["filename"],
                avg_read_length=float(row["avg_read_length"]),
                max_read_length=int(row["max_read_length"]),
                min_read_length=int(row["min_read_length"]),
                n_reads=int(row["n_reads"]),
                reference=row.get("reference") or "",
                source_fastq=row.get("source_fastq") or "",
            )
            for row in csv.DictReader(fh, delimiter="\t")
        ]


def _write_tsv(path: Path, columns: List[str], rows: List[list]) -> None:
    tmp = path.with_name(f".{path.name}.tmp")
    with open(tmp, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
    os.replace(tmp, path)


def write_organisms(data_dir: Path, organisms: List[Organism]) -> None:
    _write_tsv(
        data_dir / ORGANISMS_TSV,
        COLUMNS,
        [
            [
                o.organism_id, o.filename, f"{o.avg_read_length:.2f}", o.max_read_length,
                o.min_read_length, o.n_reads, o.reference, o.source_fastq,
            ]
            for o in organisms
        ],
    )


Fingerprint = Tuple[int, int]


def _load_cache(data_dir: Path) -> Dict[str, Tuple[Fingerprint, FastqStats]]:
    path = data_dir / SYNC_CACHE
    if not path.exists():
        return {}
    cache = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cache[row["filename"]] = (
                (int(row["size"]), int(row["mtime_ns"])),
                FastqStats(
                    int(row["n_reads"]), int(row["min_read_length"]),
                    int(row["max_read_length"]), int(row["total_bases"]),
                ),
            )
    return cache


def _write_cache(data_dir: Path, cache: Dict[str, Tuple[Fingerprint, FastqStats]]) -> None:
    _write_tsv(
        data_dir / SYNC_CACHE,
        CACHE_COLUMNS,
        [
            [name, fp[0], fp[1], s.n_reads, s.min_read_length, s.max_read_length, s.total_bases]
            for name, (fp, s) in sorted(cache.items())
        ],
    )


def scan_read_files(data_dir: Path) -> List[str]:
    """FASTQ files in the single-organism read folders, relative to data_dir."""
    found = []
    for category in READ_CATEGORIES:
        folder = reads_dir(data_dir, category)
        if not folder.is_dir():
            continue
        for path in sorted(folder.iterdir()):
            if path.is_file() and not path.name.startswith(".") and is_fastq(path.name):
                found.append(path.relative_to(data_dir).as_posix())
    return found


def sync_reads(
    data_dir: Path,
    explicit_ids: Optional[Mapping[str, str]] = None,
    sources: Optional[Mapping[str, Tuple[str, str]]] = None,
    rescan: bool = False,
    log: Callable[[str], None] = print,
) -> List[Organism]:
    """Bring organisms.tsv in line with the FASTQs on disk and return its rows.

    explicit_ids maps a filename (relative to data_dir) to the organism_id it
    should get instead of the default one, and sources maps it to the
    (reference, source_fastq) its reads were extracted with.
    """
    explicit_ids = dict(explicit_ids or {})
    sources = dict(sources or {})
    for category in READ_CATEGORIES:
        reads_dir(data_dir, category).mkdir(parents=True, exist_ok=True)

    existing = {o.filename: o for o in load_organisms(data_dir)}
    cache = _load_cache(data_dir)
    new_cache: Dict[str, Tuple[Fingerprint, FastqStats]] = {}
    rows: List[Organism] = []
    added, updated, skipped = [], [], []

    for filename in scan_read_files(data_dir):
        if filename in explicit_ids:
            organism_id = explicit_ids[filename]
        elif filename in existing:
            organism_id = existing[filename].organism_id
        else:
            organism_id = strip_fastq_suffix(Path(filename).name)
        reason = invalid_id_reason(organism_id)
        if reason:
            skipped.append(f"{filename}: {reason}; rename the file")
            continue

        stat = (data_dir / filename).stat()
        fingerprint = (stat.st_size, stat.st_mtime_ns)
        cached = cache.get(filename)
        if not rescan and cached and cached[0] == fingerprint:
            stats = cached[1]
        else:
            log(f"[stats]   reading {filename}")
            stats = fastq_stats(data_dir / filename)
        new_cache[filename] = (fingerprint, stats)

        if filename in sources:
            reference, source_fastq = sources[filename]
        elif filename in existing:
            reference, source_fastq = existing[filename].reference, existing[filename].source_fastq
        else:
            reference, source_fastq = "", ""
        row = Organism(
            organism_id, filename, round(stats.avg_read_length, 2),
            stats.max_read_length, stats.min_read_length, stats.n_reads,
            reference, source_fastq,
        )
        if filename not in existing:
            added.append(row)
        elif existing[filename] != row:
            updated.append(row)
        rows.append(row)

    by_id: Dict[str, List[str]] = {}
    for row in rows:
        by_id.setdefault(row.organism_id, []).append(row.filename)
    clashes = {oid: names for oid, names in by_id.items() if len(names) > 1}
    if clashes:
        detail = "\n".join(f"  {oid}: {', '.join(names)}" for oid, names in sorted(clashes.items()))
        raise SyncError(
            "organism_id must be unique, but several files map to the same id:\n"
            f"{detail}\nrename one of the files (or edit its id in {ORGANISMS_TSV}) and sync again"
        )

    removed = sorted(set(existing) - {row.filename for row in rows})
    write_organisms(data_dir, rows)
    _write_cache(data_dir, new_cache)

    for row in added:
        log(f"[added]   {row.organism_id:<24} {row.filename} ({row.n_reads} reads)")
    for row in updated:
        log(f"[updated] {row.organism_id:<24} {row.filename} ({row.n_reads} reads)")
    for filename in removed:
        log(f"[removed] {existing[filename].organism_id:<24} {filename} (file no longer present)")
    for message in skipped:
        log(f"[skipped] {message}")
    return rows
