#!/usr/bin/env python3
"""Extract the reads of one organism from datasets by mapping them to its reference.

This is the only way to add an organism: every read of the datasets' FASTQs
with a primary minimap2 alignment to REFERENCE is written, unchanged, to
data/<category>_reads/FILENAME, and the file is registered in data/organisms.tsv
under ORGANISM_ID with its read statistics, the breadth of coverage and the
depth range on the reference, the reference, and the source dataset(s). Where
each read aligns is saved in data/alignments/<category>_reads/FILENAME.alignments.tsv.gz,
which create_sispa_run.py uses for depths.

Each DATASET is one of:
  - a dataset id: a public dataset from publicly_available_datasets.tsv (fetch
    it with download_datasets.sh), or one of yours from data/own_datasets.tsv;
  - a folder in data/raw_data/ that is in neither table yet;
  - a FASTQ file or a folder of FASTQs anywhere on the machine.
Data that is in no table yet is added to data/own_datasets.tsv (created when
needed) under an id made from its name, with its path if it is outside
data/raw_data/; using the same path again reuses that row. A path inside a
known dataset, e.g. one run of a public dataset, uses just that part of it.
organisms.tsv records the source as <dataset_id> or <dataset_id>/<part>.

minimap2 runs as in the vimop pipeline (https://github.com/opr-group-bnitm/vimop):
-ax map-ont --secondary=no, plus -k 11 -w 5 for --category virus (vimop's
map_to_ref); references over 500 MB are indexed in 2G parts (vimop's
filter_virus_target). Like vimop, every read with a primary alignment is kept,
whatever its mapping quality. Depth is counted like vimop's
samtools depth -aa -J, over all reference sequences together.

REFERENCE is a FASTA (optionally gzipped) or a prebuilt minimap2 .mmi index; a
bare name is also looked up in data/references/.

Examples:
    ./add_from_ref.py NC_045512.2.fasta COVID COVID.fastq.gz SARS2-BC
    ./add_from_ref.py NC_045512.2.fasta COVID-1 COVID-1.fastq.gz data/raw_data/SARS2-BC/SRR15356294_1.fastq.gz
    ./add_from_ref.py --category host ARS-UCD2.0.fna cattle cattle.fastq.gz BOV-6760 BOV-6763
    ./add_from_ref.py lasv.fasta LASV LASV.fastq.gz /Volumes/runs/2025-06-01/fastq_pass/barcode05
"""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from workbench.alignments import AlignmentWriter, alignments_path
from workbench.coverage import combine
from workbench.datasets import (
    Dataset, DatasetError, add_own_dataset, fastqs_at, load_datasets, own_datasets_tsv,
)
from workbench.fastq import (
    FastqError, FastqStats, is_fastq, open_write, read_fastq, read_name, strip_fastq_suffix, write_record,
)
from workbench.organisms import (
    Organism, RegistryError, invalid_id_reason, load_organisms, prune_missing, register, unregistered_fastqs,
)
from workbench.mapping import Mapping, MappingError, minimap2_args
from workbench.paths import READ_CATEGORIES, default_data_dir, inside, raw_data_dir, reads_dir


class AddError(RuntimeError):
    pass


# print coverage per reference sequence only for references with few of them
MAX_SEQUENCES_SHOWN = 10


def data_relative(path: Path, data_dir: Path) -> str:
    """path relative to data_dir when it lies inside it, else absolute."""
    for inner, outer in ((Path(os.path.abspath(path)), Path(os.path.abspath(data_dir))),
                         (path.resolve(), data_dir.resolve())):
        try:
            return inner.relative_to(outer).as_posix()
        except ValueError:
            pass
    return os.path.abspath(path)


def resolve_reference(reference: str, data_dir: Path) -> Path:
    path = Path(reference)
    if path.exists():
        return path
    candidate = data_dir / "references" / reference
    if candidate.exists():
        return candidate
    raise AddError(f"reference not found: {reference} (also looked in {candidate.parent})")


@dataclass
class Source:
    label: str  # what organisms.tsv records: <dataset_id> or <dataset_id>/<part>
    files: List[Path]


def containing_dataset(path: Path, datasets: Dict[str, Dataset], data_dir: Path) -> Tuple[Optional[str], str]:
    """The dataset a path lies in and the part of it the path points to."""
    raw = raw_data_dir(data_dir)
    rel = inside(path, raw)
    if rel == "":
        raise AddError(f"give a dataset in {raw}/, not the whole folder")
    if rel is not None:
        dataset_id, _, part = rel.partition("/")
        return dataset_id, part
    best: Tuple[Optional[str], str] = (None, "")
    for dataset in datasets.values():
        if dataset.path:
            rel = inside(path, dataset.location(data_dir))
            if rel is not None and (best[0] is None or len(rel) < len(best[1])):
                best = (dataset.dataset_id, rel)
    return best


def new_dataset_id(path: Path, datasets: Dict[str, Dataset], data_dir: Path) -> str:
    """An unused dataset id made from the name of path."""
    name = strip_fastq_suffix(Path(os.path.abspath(path)).name)
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-") or "dataset"
    raw = raw_data_dir(data_dir)
    taken = set(datasets) | ({p.name for p in raw.iterdir()} if raw.is_dir() else set())
    dataset_id, n = base, 1
    while dataset_id in taken:
        n += 1
        dataset_id = f"{base}-{n}"
    return dataset_id


def resolve_sources(items: List[str], data_dir: Path) -> List[Source]:
    """The FASTQs behind each DATASET argument; data in no table yet is added
    to own_datasets.tsv."""
    datasets = load_datasets(data_dir)
    raw = raw_data_dir(data_dir)
    sources: List[Source] = []
    used: Dict[str, str] = {}  # real path of each FASTQ -> label it came with
    for item in items:
        if "/" not in item and (item in datasets or os.path.lexists(raw / item)):
            dataset_id, part = item, ""
        else:
            path = Path(os.path.expanduser(item))
            if not os.path.lexists(path):
                raise AddError(
                    f"{item!r} is neither a dataset id (see publicly_available_datasets.tsv, "
                    f"{own_datasets_tsv(data_dir)} and {raw}/) nor an existing file or folder"
                )
            dataset_id, part = containing_dataset(path, datasets, data_dir)
            if dataset_id is None:  # FASTQs outside every known dataset
                dataset_id = new_dataset_id(path, datasets, data_dir)
                location = os.path.abspath(path)
                table = add_own_dataset(data_dir, dataset_id, location)
                datasets[dataset_id] = Dataset(dataset_id, table, location)
                print(f"[dataset] {dataset_id}: added to {table} for {location} -- describe it there")
        if dataset_id not in datasets:  # a folder in data/raw_data/ that is in neither table
            table = add_own_dataset(data_dir, dataset_id)
            datasets[dataset_id] = Dataset(dataset_id, table)
            print(f"[dataset] {dataset_id}: added to {table} as your own dataset -- describe it there")

        dataset = datasets[dataset_id]
        location = dataset.location(data_dir)
        label = f"{dataset_id}/{part}" if part else dataset_id
        if dataset.path and not os.path.lexists(location):
            raise AddError(f"dataset {dataset_id} is at {location}, which does not exist (is the drive mounted?)")
        files = fastqs_at(location / part if part else location)
        if not files:
            if dataset.is_public and not part:
                raise AddError(f"dataset {dataset_id} has no FASTQ files in {location}/; "
                               f"download it with ./download_datasets.sh {dataset_id}")
            raise AddError(f"{label}: no FASTQ files (.fastq/.fq, optionally .gz) in {location / part}")
        for fastq in files:
            real = os.path.realpath(fastq)
            if real in used:
                raise AddError(f"{fastq} is included twice, via {used[real]} and {label}")
            used[real] = label
        kind = "public" if dataset.is_public else f"own, {location}" if dataset.path else "own"
        print(f"[dataset] {label} ({kind}): {len(files)} FASTQ file(s)")
        sources.append(Source(label, files))
    return sources


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference", help="reference FASTA or minimap2 index of the organism")
    parser.add_argument("organism_id", help="id to register the reads under in organisms.tsv")
    parser.add_argument("filename", help="output file name, e.g. COVID.fastq.gz (.fastq.gz is added if missing)")
    parser.add_argument("dataset", nargs="+",
                        help="where to extract the reads from: dataset ids, or FASTQ files / folders anywhere")
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
    except (AddError, MappingError, RegistryError, DatasetError, FastqError) as err:
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

    filename = args.filename if is_fastq(args.filename) else f"{args.filename}.fastq.gz"
    out_dir = reads_dir(data_dir, args.category)
    out_path = out_dir / filename
    rel_name = out_path.relative_to(data_dir).as_posix()
    if out_path.exists() and not args.force:
        raise AddError(f"{out_path} already exists (use --force to overwrite)")
    for organism in prune_missing(data_dir):
        if organism.organism_id == args.organism_id and organism.filename != rel_name:
            raise AddError(f"organism_id {args.organism_id!r} is already used by {organism.filename}")
    sources = resolve_sources(args.dataset, data_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    split_dir = tempfile.TemporaryDirectory(prefix="add_from_ref_")
    extra = minimap2_args(reference, args.category, Path(split_dir.name) / "map_index")
    extra += shlex.split(args.minimap2_args)
    mapping = Mapping(reference, args.preset, args.threads, args.min_mapq, args.min_aligned_fraction, extra)
    tmp_path = out_dir / f".partial.{filename}"  # hidden, keeps the .gz suffix
    final_alignments = alignments_path(data_dir, rel_name)
    final_alignments.parent.mkdir(parents=True, exist_ok=True)
    tmp_alignments = final_alignments.with_name(f".partial.{final_alignments.name}")
    alignments: Optional[AlignmentWriter] = None
    stats = FastqStats()
    n_in = 0
    try:
        with open_write(tmp_path) as out:
            for source in sources:
                for fastq in source.files:
                    kept = mapping.read_names(fastq)
                    if alignments is None:  # the reference sequences are known after the first mapping
                        alignments = AlignmentWriter(tmp_alignments, mapping.lengths)
                    file_in = file_out = 0
                    for record in read_fastq(fastq):
                        file_in += 1
                        name = read_name(record[0])
                        if name not in kept and name[-2:] in ("/1", "/2"):
                            name = name[:-2]  # minimap2 drops a /1 or /2 mate suffix from read names
                        if name in kept:
                            alignments.add(stats.n_reads, kept[name])
                            write_record(out, record)
                            stats.add(len(record[1]))
                            file_out += 1
                    print(f"[extract] {data_relative(fastq, data_dir)}: {file_out} of {file_in} reads "
                          f"map to {reference.name}")
                    n_in += file_in
        if stats.n_reads == 0:
            raise AddError(f"no reads mapped to {reference}; nothing written")
        alignments.close(stats.n_reads)
        alignments = None
        os.replace(tmp_alignments, final_alignments)
        os.replace(tmp_path, out_path)
    finally:
        if alignments is not None:
            alignments.out.close()
        tmp_path.unlink(missing_ok=True)
        tmp_alignments.unlink(missing_ok=True)
        split_dir.cleanup()
    print(f"[written] {out_path} ({stats.n_reads} of {n_in} reads, {100 * stats.n_reads / n_in:.2f}%)")

    per_sequence = mapping.depth.summarize(mapping.lengths)
    if len(per_sequence) <= MAX_SEQUENCES_SHOWN:
        for name, s in per_sequence.items():
            print(f"[coverage] {name} ({s.length} bp): breadth {s.breadth:.2f}%, "
                  f"depth {s.min_depth}-{s.max_depth}")
    total = combine(list(per_sequence.values()))
    print(f"[coverage] {reference.name}, {len(per_sequence)} sequence(s), {total.length} bp: "
          f"breadth {total.breadth:.2f}%, depth {total.min_depth}-{total.max_depth}")

    register(data_dir, Organism(
        organism_id=args.organism_id,
        filename=rel_name,
        avg_read_length=round(stats.avg_read_length, 2),
        max_read_length=stats.max_read_length,
        min_read_length=stats.min_read_length,
        n_reads=stats.n_reads,
        breadth_coverage=round(total.breadth, 2),
        min_depth=total.min_depth,
        max_depth=total.max_depth,
        reference=data_relative(reference, data_dir),
        source_dataset=";".join(source.label for source in sources),
    ))
    for name in unregistered_fastqs(data_dir, load_organisms(data_dir)):
        print(f"[ignored] {name} was not created by add_from_ref.py and cannot be used; "
              f"put its reads in {raw_data_dir(data_dir)}/<dataset_id>/ and extract them from there")
    return 0


if __name__ == "__main__":
    sys.exit(main())
