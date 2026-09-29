"""The two dataset tables, and where each dataset's FASTQs are.

publicly_available_datasets.tsv (in git) lists public datasets, which
download_datasets.sh fetches into data/raw_data/<dataset_id>/.

data/own_datasets.tsv lists your own datasets; it is local to your checkout,
not in git, and created when first needed. An own dataset's FASTQs are in
data/raw_data/<dataset_id>/ unless its ``path`` column points elsewhere on the
machine, to a folder or a single FASTQ (absolute, or relative to the data
directory).

Dataset ids are unique across both tables. Symlinks are followed wherever
FASTQs are looked for.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from . import paths
from .fastq import is_fastq

OWN_COLUMNS = ["dataset_id", "sample_type", "organism_target", "notes", "path"]

# the same rule download_datasets.sh applies: ids are folder names
ID_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")


class DatasetError(RuntimeError):
    pass


@dataclass
class Dataset:
    dataset_id: str
    table: Path  # the table that lists it
    path: str = ""  # where an own dataset lives, if not in data/raw_data/

    @property
    def is_public(self) -> bool:
        return self.table.name != paths.OWN_DATASETS_TSV

    def location(self, data_dir: Path) -> Path:
        if self.path:
            return data_dir / os.path.expanduser(self.path)  # an absolute path replaces data_dir
        return paths.raw_data_dir(data_dir) / self.dataset_id


def own_datasets_tsv(data_dir: Path) -> Path:
    return data_dir / paths.OWN_DATASETS_TSV


def _header_index(lines: List[str]) -> Optional[int]:
    for i, line in enumerate(lines):
        if line.strip() and not line.startswith("#"):
            return i
    return None


def _columns(line: str) -> List[str]:
    return [c.strip() for c in line.rstrip("\r\n").split("\t")]


def ensure_own_table(data_dir: Path) -> Path:
    """Create data/own_datasets.tsv with just its header if it does not exist
    yet, and add columns an older table lacks."""
    path = own_datasets_tsv(data_dir)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\t".join(OWN_COLUMNS) + "\n")
        return path
    lines = path.read_text().splitlines(keepends=True)
    i = _header_index(lines)
    if i is None:
        path.write_text("".join(lines) + "\t".join(OWN_COLUMNS) + "\n")
        return path
    missing = [c for c in OWN_COLUMNS if c not in _columns(lines[i])]
    if missing:  # rows written before a column existed read as empty there
        lines[i] = lines[i].rstrip("\r\n") + "".join("\t" + c for c in missing) + "\n"
        path.write_text("".join(lines))
    return path


def _rows(path: Path) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        lines = [line for line in fh if line.strip() and not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    if "dataset_id" not in [(c or "").strip() for c in reader.fieldnames or []]:
        raise DatasetError(f"{path} has no dataset_id column")
    rows = []
    for row in reader:
        row = {(k or "").strip(): (v or "").strip() for k, v in row.items() if k is not None}
        if row.get("dataset_id"):
            rows.append(row)
    return rows


def load_datasets(data_dir: Path) -> Dict[str, Dataset]:
    """Every dataset of both tables by id; creates the own table if needed."""
    datasets: Dict[str, Dataset] = {}
    own = ensure_own_table(data_dir)
    for table in (paths.PUBLIC_DATASETS_TSV, own):
        if not table.exists():
            continue
        for row in _rows(table):
            dataset_id = row["dataset_id"]
            if not ID_PATTERN.match(dataset_id):
                raise DatasetError(f"dataset_id {dataset_id!r} in {table} may only contain "
                                   "letters, digits, '.', '_' and '-'")
            if dataset_id in datasets:
                first = datasets[dataset_id].table
                where = f"twice in {table}" if first == table else f"in both {first} and {table}"
                raise DatasetError(f"dataset_id {dataset_id!r} is listed {where}; ids must be unique")
            datasets[dataset_id] = Dataset(dataset_id, table, row.get("path", "") if table == own else "")
    return datasets


def add_own_dataset(data_dir: Path, dataset_id: str, path: str = "", notes: str = "") -> Path:
    """Append a row for dataset_id, described only by notes, to own_datasets.tsv."""
    if not ID_PATTERN.match(dataset_id):
        raise DatasetError(f"dataset_id {dataset_id!r} may only contain letters, digits, '.', '_' and '-'")
    table = ensure_own_table(data_dir)
    lines = table.read_text().splitlines(keepends=True)
    values = {"dataset_id": dataset_id, "path": path, "notes": notes.replace("\t", " ")}
    row = "\t".join(values.get(column, "") for column in _columns(lines[_header_index(lines)]))
    with open(table, "a") as fh:
        if not lines[-1].endswith("\n"):
            fh.write("\n")
        fh.write(row + "\n")
    return table


def _broken_link(path: Path) -> DatasetError:
    return DatasetError(f"{path} links to {os.readlink(path)}, which does not exist (is the drive mounted?)")


def fastqs_at(location: Path) -> List[Path]:
    """The FASTQs at location: the file itself, or every .fastq/.fq(.gz) below
    the folder, skipping hidden files and folders (download bookkeeping,
    partial files).

    Symlinks to files and folders are followed at any depth. A broken link is an
    error rather than a silently smaller dataset, and a file reached through
    several links is used once."""
    if location.is_symlink() and not location.exists():
        raise _broken_link(location)
    if location.is_file():
        return [location]
    if not location.is_dir():
        return []
    found: Dict[str, Path] = {}  # real path -> first path it was reached by
    seen_dirs = set()
    for root, dirs, files in os.walk(location, followlinks=True):
        real_root = os.path.realpath(root)
        if real_root in seen_dirs:  # a link loop, or a folder linked twice
            dirs[:] = []
            continue
        seen_dirs.add(real_root)
        dirs[:] = sorted(d for d in dirs if not d.startswith("."))
        for name in sorted(files):
            if name.startswith("."):
                continue
            path = Path(root) / name
            if path.is_symlink() and not path.exists():  # os.walk lists broken links as files
                raise _broken_link(path)
            if is_fastq(name):
                found.setdefault(os.path.realpath(path), path)
    return sorted(found.values())
