"""Default locations, all relative to the repository root.

``SISPA_DATA_DIR`` and ``SISPA_OUTPUT_DIR`` override the defaults, and every
tool also takes ``--data-dir`` / output flags.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# the three folders whose FASTQs hold reads of exactly one organism
READ_CATEGORIES = ("virus", "bacteria", "host")

ORGANISMS_TSV = "organisms.tsv"
SYNC_CACHE = ".organisms_cache.tsv"


def default_data_dir() -> Path:
    return Path(os.environ.get("SISPA_DATA_DIR", REPO_ROOT / "data"))


def default_output_dir() -> Path:
    return Path(os.environ.get("SISPA_OUTPUT_DIR", REPO_ROOT / "output"))


def reads_dir(data_dir: Path, category: str) -> Path:
    return data_dir / f"{category}_reads"
