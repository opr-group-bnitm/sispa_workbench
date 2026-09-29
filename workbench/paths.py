"""Default locations, all relative to the repository root.

``SISPA_DATA_DIR`` and ``SISPA_OUTPUT_DIR`` override the defaults, and every
tool also takes ``--data-dir`` / output flags.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parent.parent

# curated public datasets, tracked in git
PUBLIC_DATASETS_TSV = REPO_ROOT / "publicly_available_datasets.tsv"

# the three folders whose FASTQs hold reads of exactly one organism; only
# add_from_ref.py writes to them
READ_CATEGORIES = ("virus", "bacteria", "host")

# local tables inside the data directory, not tracked in git
ORGANISMS_TSV = "organisms.tsv"
OWN_DATASETS_TSV = "own_datasets.tsv"


def default_data_dir() -> Path:
    return Path(os.environ.get("SISPA_DATA_DIR", REPO_ROOT / "data"))


def default_output_dir() -> Path:
    return Path(os.environ.get("SISPA_OUTPUT_DIR", REPO_ROOT / "output"))


def reads_dir(data_dir: Path, category: str) -> Path:
    return data_dir / f"{category}_reads"


def raw_data_dir(data_dir: Path) -> Path:
    return data_dir / "raw_data"


def alignments_dir(data_dir: Path) -> Path:
    return data_dir / "alignments"


def inside(path: Path, root: Path) -> Optional[str]:
    """path relative to root ("" for root itself) if it lies inside root, else None."""
    for p, r in ((os.path.abspath(path), os.path.abspath(root)),
                 (os.path.realpath(path), os.path.realpath(root))):
        if p == r:
            return ""
        if p.startswith(r.rstrip(os.sep) + os.sep):
            return Path(os.path.relpath(p, r)).as_posix()
    return None
