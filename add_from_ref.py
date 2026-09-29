#!/usr/bin/env python3
"""Extract the reads of one organism by mapping a FASTQ against its reference.

Every read with a primary minimap2 alignment to REFERENCE is written, unchanged,
to data/<category>_reads/FILENAME, and the file is registered in
data/organisms.tsv under ORGANISM_ID (see sync_reads.py), together with the
reference and the input FASTQ(s) it came from.

minimap2 runs with the settings of the vimop pipeline
(https://github.com/opr-group-bnitm/vimop): -x map-ont --secondary=no, plus
-k 11 -w 5 for --category virus (vimop's map_to_ref); references over 500 MB
are indexed in 2G parts (vimop's filter_virus_target). Like vimop, every read
with a primary alignment is kept, whatever its mapping quality.

REFERENCE is a FASTA (optionally gzipped) or a prebuilt minimap2 .mmi index; a
bare name is also looked up in data/references/. Several FASTQs can be given,
e.g. all runs of a dataset in data/raw_data/<dataset_id>/.

Examples:
    ./add_from_ref.py NC_045512.2.fasta COVID COVID.fastq.gz data/raw_data/SARS2-BC/*.fastq.gz
    ./add_from_ref.py --category host human.mmi human human.fastq.gz data/mixed_reads/sample1.fastq.gz
    ./add_from_ref.py --min-mapq 20 --min-aligned-fraction 0.8 lasv.fasta LASV LASV.fastq.gz run1.fastq.gz run2.fastq.gz
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Set, Tuple

from workbench.fastq import FastqError, is_fastq, open_write, read_fastq, read_name, write_record
from workbench.organisms import SyncError, invalid_id_reason, load_organisms, sync_reads
from workbench.paths import READ_CATEGORIES, default_data_dir, reads_dir


class AddError(RuntimeError):
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


def mapped_read_names(
    reference: Path, fastq: Path, preset: str, threads: int,
    min_mapq: int, min_fraction: float, extra_args: List[str],
) -> Set[str]:
    """Names of the reads in fastq whose primary alignments pass the filters."""
    # -c: base-level alignment, so a read counts as mapped exactly when it
    # would be a mapped record in minimap2's SAM output (what vimop uses)
    cmd = ["minimap2", "-c", "-x", preset, "-t", str(threads), *extra_args, str(reference), str(fastq)]
    print(f"[map]     {' '.join(shlex.quote(c) for c in cmd)}")

    # per read: query length and the query intervals of its passing alignments
    hits: Dict[str, Tuple[int, List[Tuple[int, int]]]] = {}
    with tempfile.TemporaryFile() as err:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True, bufsize=1 << 20)
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or "tp:A:S" in fields[12:] or int(fields[11]) < min_mapq:
                continue
            qname, qlen = fields[0], int(fields[1])
            hits.setdefault(qname, (qlen, []))[1].append((int(fields[2]), int(fields[3])))
        if proc.wait() != 0:
            err.seek(0)
            tail = err.read().decode(errors="replace").strip().splitlines()[-15:]
            raise AddError(f"minimap2 failed on {fastq} (exit code {proc.returncode}):\n  " + "\n  ".join(tail))

    return {
        name for name, (qlen, intervals) in hits.items()
        if qlen == 0 or covered_length(intervals) / qlen >= min_fraction
    }


def data_relative(path: Path, data_dir: Path) -> str:
    """path relative to data_dir when it lies inside it, else absolute."""
    path = path.resolve()
    try:
        return path.relative_to(data_dir.resolve()).as_posix()
    except ValueError:
        return str(path)


def resolve_reference(reference: str, data_dir: Path) -> Path:
    path = Path(reference)
    if path.exists():
        return path
    candidate = data_dir / "references" / reference
    if candidate.exists():
        return candidate
    raise AddError(f"reference not found: {reference} (also looked in {candidate.parent})")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference", help="reference FASTA or minimap2 index of the organism")
    parser.add_argument("organism_id", help="id to register the reads under in organisms.tsv")
    parser.add_argument("filename", help="output file name, e.g. COVID.fastq.gz (.fastq.gz is added if missing)")
    parser.add_argument("fastq", nargs="+", type=Path, help="FASTQ file(s) to extract reads from")
    parser.add_argument("--category", choices=READ_CATEGORIES, default="virus",
                        help="folder the reads go to: data/<category>_reads/ (default: %(default)s)")
    parser.add_argument("--preset", default="map-ont",
                        help="minimap2 preset: map-ont, lr:hq, map-pb, map-hifi, sr, ... (default: %(default)s)")
    parser.add_argument("--threads", type=int, default=min(4, os.cpu_count() or 1),
                        help="minimap2 threads (default: %(default)s)")
    parser.add_argument("--min-mapq", type=int, default=0,
                        help="ignore alignments below this mapping quality (default: %(default)s)")
    parser.add_argument("--min-aligned-fraction", type=float, default=0.0, metavar="FRACTION",
                        help="keep a read only if at least this fraction of it aligns, e.g. 0.8 to drop "
                             "chimeras with long unaligned stretches (default: %(default)s = any alignment)")
    parser.add_argument("--minimap2-args", default="", metavar="ARGS",
                        help="extra minimap2 arguments; they override the vimop defaults, "
                             "e.g. \"-k 15 -w 10\" for map-ont's own k-mer and window size")
    parser.add_argument("--force", action="store_true", help="overwrite an existing output file")
    parser.add_argument("--data-dir", type=Path, default=default_data_dir(),
                        help="workbench data directory (default: %(default)s)")
    args = parser.parse_args(argv)

    try:
        return run(args)
    except (AddError, SyncError, FastqError) as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1


def run(args) -> int:
    data_dir: Path = args.data_dir
    if shutil.which("minimap2") is None:
        raise AddError("minimap2 not found on PATH (conda env create -f environment.yml)")
    reason = invalid_id_reason(args.organism_id)
    if reason:
        raise AddError(reason)
    if Path(args.filename).name != args.filename:
        raise AddError(f"filename must be a plain file name, not a path: {args.filename} "
                       "(use --category to pick the folder)")
    if not 0.0 <= args.min_aligned_fraction <= 1.0:
        raise AddError("--min-aligned-fraction must be between 0 and 1")
    reference = resolve_reference(args.reference, data_dir)
    for fastq in args.fastq:
        if not fastq.is_file():
            raise AddError(f"FASTQ not found: {fastq}")

    filename = args.filename if is_fastq(args.filename) else f"{args.filename}.fastq.gz"
    out_dir = reads_dir(data_dir, args.category)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / filename
    rel_name = out_path.relative_to(data_dir).as_posix()
    if out_path.exists() and not args.force:
        raise AddError(f"{out_path} already exists (use --force to overwrite)")
    for organism in load_organisms(data_dir):
        if organism.organism_id == args.organism_id and organism.filename != rel_name:
            raise AddError(f"organism_id {args.organism_id!r} is already used by {organism.filename}")

    split_dir = tempfile.TemporaryDirectory(prefix="add_from_ref_")
    extra = minimap2_args(reference, args.category, Path(split_dir.name) / "map_index")
    extra += shlex.split(args.minimap2_args)
    tmp_path = out_dir / f".partial.{filename}"  # hidden from sync, keeps the .gz suffix
    n_in = n_out = 0
    try:
        with open_write(tmp_path) as out:
            for fastq in args.fastq:
                names = mapped_read_names(reference, fastq, args.preset, args.threads,
                                          args.min_mapq, args.min_aligned_fraction, extra)
                file_in = file_out = 0
                for record in read_fastq(fastq):
                    file_in += 1
                    name = read_name(record[0])
                    # minimap2 drops a /1 or /2 mate suffix from read names
                    if name in names or (name[-2:] in ("/1", "/2") and name[:-2] in names):
                        write_record(out, record)
                        file_out += 1
                print(f"[extract] {fastq}: {file_out} of {file_in} reads map to {reference.name}")
                n_in += file_in
                n_out += file_out
        if n_out == 0:
            raise AddError(f"no reads mapped to {reference}; nothing written")
        os.replace(tmp_path, out_path)
    finally:
        tmp_path.unlink(missing_ok=True)
        split_dir.cleanup()

    print(f"[written] {out_path} ({n_out} of {n_in} reads, {100 * n_out / n_in:.2f}%)")
    source = (
        data_relative(reference, data_dir),
        ";".join(data_relative(fastq, data_dir) for fastq in args.fastq),
    )
    sync_reads(data_dir, explicit_ids={rel_name: args.organism_id}, sources={rel_name: source})
    return 0


if __name__ == "__main__":
    sys.exit(main())
