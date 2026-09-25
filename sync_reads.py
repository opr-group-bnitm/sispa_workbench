#!/usr/bin/env python3
"""Register the single-organism FASTQs in data/{virus,bacteria,host}_reads/.

Every FASTQ (.fastq, .fq, optionally .gz) in those folders becomes a row of
data/organisms.tsv:

    organism_id  filename  avg_read_length  max_read_length  min_read_length  n_reads

A file you drop in yourself gets its name without the suffix as organism_id
(virus_reads/COVID.fastq.gz -> COVID). Files that disappeared are dropped from
the table, and only new or changed files are read, so re-running is cheap.

Examples:
    ./sync_reads.py
    ./sync_reads.py --rescan          # recompute statistics for every file
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from workbench.fastq import FastqError
from workbench.organisms import SyncError, sync_reads
from workbench.paths import ORGANISMS_TSV, default_data_dir


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, default=default_data_dir(),
                        help="workbench data directory (default: %(default)s)")
    parser.add_argument("--rescan", action="store_true",
                        help="recompute read statistics even for files that did not change")
    args = parser.parse_args(argv)

    try:
        rows = sync_reads(args.data_dir, rescan=args.rescan)
    except (SyncError, FastqError) as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1
    print(f"{len(rows)} organism(s) in {args.data_dir / ORGANISMS_TSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
