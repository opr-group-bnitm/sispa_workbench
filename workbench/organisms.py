"""The organisms.tsv registry of single-organism FASTQ files.

Every row is written by add_from_ref.py, which extracts the reads of one
organism from datasets in data/raw_data/ into data/{virus,bacteria,host}_reads/.
FASTQs put into those folders any other way are not registered and cannot be
used.

Columns:
    organism_id       id used in compositions ("COVID 5000, human all")
    filename          the FASTQ, relative to the data directory
    avg_read_length, max_read_length, min_read_length, n_reads
    breadth_coverage  % of reference positions covered by at least one read
    min_depth, max_depth
                      depth range over all reference positions (see coverage.py)
    reference         the reference the reads were mapped to, relative to the
                      data directory when inside it
    source_dataset    the dataset(s) the reads were extracted from, ';'-separated
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import astuple, dataclass
from pathlib import Path
from typing import Callable, List, Optional

from .alignments import alignments_path
from .fastq import is_fastq
from .paths import ORGANISMS_TSV, READ_CATEGORIES, reads_dir

COLUMNS = [
    "organism_id", "filename", "avg_read_length", "max_read_length", "min_read_length", "n_reads",
    "breadth_coverage", "min_depth", "max_depth", "reference", "source_dataset",
]

# ids are used unquoted in composition lists ("COVID 5000, human all") and in
# run names, so keep them free of whitespace, commas and path separators
ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class RegistryError(RuntimeError):
    pass


@dataclass
class Organism:
    organism_id: str
    filename: str
    avg_read_length: float
    max_read_length: int
    min_read_length: int
    n_reads: int
    breadth_coverage: float
    min_depth: int
    max_depth: int
    reference: str
    source_dataset: str

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
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise RegistryError(
                f"{path} lacks the column(s) {', '.join(missing)}: it was written by an older "
                "version. Delete it and extract the organisms again with add_from_ref.py"
            )
        return [
            Organism(
                organism_id=row["organism_id"],
                filename=row["filename"],
                avg_read_length=float(row["avg_read_length"]),
                max_read_length=int(row["max_read_length"]),
                min_read_length=int(row["min_read_length"]),
                n_reads=int(row["n_reads"]),
                breadth_coverage=float(row["breadth_coverage"]),
                min_depth=int(row["min_depth"]),
                max_depth=int(row["max_depth"]),
                reference=row["reference"],
                source_dataset=row["source_dataset"],
            )
            for row in reader
        ]


def write_organisms(data_dir: Path, organisms: List[Organism]) -> None:
    path = data_dir / ORGANISMS_TSV
    tmp = path.with_name(f".{path.name}.tmp")
    with open(tmp, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for o in organisms:
            row = list(astuple(o))
            row[2] = f"{o.avg_read_length:.2f}"
            row[6] = f"{o.breadth_coverage:.2f}"
            writer.writerow(row)
    os.replace(tmp, path)


def prune_missing(data_dir: Path, log: Callable[[str], None] = print) -> List[Organism]:
    """Drop the rows whose FASTQ was deleted, and their alignment files, and
    return the remaining rows."""
    organisms = load_organisms(data_dir)
    present = [o for o in organisms if o.path(data_dir).is_file()]
    for o in organisms:
        if o not in present:
            alignments_path(data_dir, o.filename).unlink(missing_ok=True)
            log(f"[removed] {o.organism_id}: {o.filename} no longer exists")
    if len(present) != len(organisms):
        write_organisms(data_dir, present)
    return present


def register(data_dir: Path, organism: Organism, log: Callable[[str], None] = print) -> None:
    """Add a row to organisms.tsv, replacing the row of the same file."""
    organisms = prune_missing(data_dir, log)
    for other in organisms:
        if other.organism_id == organism.organism_id and other.filename != organism.filename:
            raise RegistryError(f"organism_id {organism.organism_id!r} is already used by {other.filename}")
    replaced = [i for i, o in enumerate(organisms) if o.filename == organism.filename]
    if replaced:
        organisms[replaced[0]] = organism
    else:
        organisms.append(organism)
    write_organisms(data_dir, organisms)
    log(f"[{'updated' if replaced else 'added'}] {organism.organism_id}: {organism.filename} "
        f"({organism.n_reads} reads) in {data_dir / ORGANISMS_TSV}")


def unregistered_fastqs(data_dir: Path, organisms: List[Organism]) -> List[str]:
    """FASTQs in the read folders that add_from_ref.py did not create."""
    known = {o.filename for o in organisms}
    found = []
    for category in READ_CATEGORIES:
        folder = reads_dir(data_dir, category)
        if folder.is_dir():
            for path in sorted(folder.iterdir()):
                name = path.relative_to(data_dir).as_posix()
                if path.is_file() and not path.name.startswith(".") and is_fastq(path.name) and name not in known:
                    found.append(name)
    return found
