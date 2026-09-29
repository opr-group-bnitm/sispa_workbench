#!/usr/bin/env python3
"""Remove the viral reads from a sequencing run, as found by vimop.

Runs the vimop pipeline (https://github.com/opr-group-bnitm/vimop) on the run
with its default settings,

    nextflow run opr-group-bnitm/vimop --fastq <run> --out_dir <vimop output> -resume

and removes the reads vimop mapped to the viruses it found: for every reference
in <sample>/tables/consensus.tsv whose consensus reached --min-recovery percent
(the table's Coverage column: positions called, not N), every read in
<sample>/consensus/<reference>.reads.bam. All other reads are written
unchanged. --min-recovery 0 removes the reads of every virus vimop reported.

FASTQ is one or more FASTQ files, folders of them (e.g. fastq_pass/barcode05)
or dataset ids (folders in data/raw_data/), all treated as one sample. The
cleaned run goes to data/output/background_fastqs/<name>_no_viral.fastq.gz,
<name> being the first input's dataset or file name, unless -o says otherwise
(written into a new folder in data/raw_data/, it is added to
data/own_datasets.tsv). The reads removed per virus are listed next to the
output in <output name>.removed_reads.tsv.

vimop's input, output, nextflow work folder and log go to
output/vimop/<output name>/ (--vimop-dir) and are kept, for vimop's report and
to reuse them; --discard-vimop deletes them once the cleaned run is written
(a failed run is always kept). --vimop-output uses the output of an earlier
vimop run of the same reads instead of running vimop; every sample in it
counts, and it is never deleted.

Examples:
    ./remove_viral_reads.py BOV-6760
    ./remove_viral_reads.py BOV-6760 --discard-vimop
    ./remove_viral_reads.py fastq_pass/barcode05 -o clean.fastq.gz --vimop-args "--targets LASV"
    ./remove_viral_reads.py run.fastq.gz --min-recovery 0 --vimop-output output/vimop/run_no_viral/output
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys
from collections import Counter
from contextlib import nullcontext
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Set, Tuple

from workbench.datasets import DatasetError, add_own_dataset, fastqs_at, load_datasets
from workbench.fastq import FastqError, open_write, read_fastq, read_name, strip_fastq_suffix, write_record
from workbench.paths import background_fastqs_dir, default_data_dir, default_output_dir, inside, raw_data_dir

VIMOP = "opr-group-bnitm/vimop"
SUFFIX = "_no_viral"


class RemoveError(RuntimeError):
    pass


def input_paths(items: List[str], data_dir: Path) -> List[Path]:
    """The inputs as paths: a dataset id stands for its folder in data/raw_data/."""
    paths = []
    for item in items:
        path = Path(os.path.expanduser(item))
        if not os.path.lexists(path) and "/" not in item and os.path.lexists(raw_data_dir(data_dir) / item):
            path = raw_data_dir(data_dir) / item
        if not os.path.lexists(path):
            raise RemoveError(f"not found: {item} (neither a path nor a dataset in {raw_data_dir(data_dir)}/)")
        paths.append(path)
    return paths


def input_fastqs(paths: List[Path]) -> List[Path]:
    files: Dict[str, Path] = {}  # a file given twice (directly and in a folder) counts once
    for path in paths:
        found = fastqs_at(path)
        if not found:
            raise RemoveError(f"no FASTQ files (.fastq/.fq, optionally .gz) in {path}")
        for fastq in found:
            files.setdefault(os.path.realpath(fastq), fastq)
    return list(files.values())


def default_out(first: Path, data_dir: Path) -> Path:
    """data/output/background_fastqs/<name>_no_viral.fastq.gz, <name> being the
    dataset the first input belongs to, or else its file name."""
    rel = inside(first, raw_data_dir(data_dir))
    name = rel.split("/")[0] if rel else strip_fastq_suffix(Path(os.path.abspath(first)).name)
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-") or "run"
    return background_fastqs_dir(data_dir) / f"{name}{SUFFIX}.fastq.gz"


def check_output(out: Path, files: List[Path], data_dir: Path) -> None:
    if os.path.realpath(out) in {os.path.realpath(f) for f in files}:
        raise RemoveError(f"{out} is one of the input files")
    raw = raw_data_dir(data_dir)
    rel = inside(out, raw)
    if rel:
        dataset = rel.split("/")[0]
        if any(inside(f, raw / dataset) is not None for f in files):
            raise RemoveError(f"{out} would be in dataset {dataset} next to its own input, which would then hold "
                              f"these reads twice; write it to a new folder such as {raw / (dataset + SUFFIX)}/")


def stage_input(files: List[Path], folder: Path) -> None:
    """vimop reads a folder of FASTQs: hard-link the run's files into one (or copy
    them where the file system cannot link). The folder name becomes vimop's
    sample name."""
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    for i, fastq in enumerate(files):
        staged = folder / (fastq.name if len(files) == 1 else f"{i + 1:03d}_{fastq.name}")
        try:
            os.link(os.path.realpath(fastq), staged)
        except OSError:
            print(f"[vimop]   cannot link {fastq} into {folder}, copying it")
            shutil.copyfile(fastq, staged)


def run_vimop(args, files: List[Path], vimop_dir: Path, name: str, output: Path) -> None:
    nextflow = shutil.which(args.nextflow)
    if nextflow is None:
        raise RemoveError(f"{args.nextflow} not found; install nextflow or give its path with --nextflow")
    if args.config and not args.config.is_file():
        raise RemoveError(f"nextflow config not found: {args.config}")
    fastq_dir = vimop_dir / "input" / name
    stage_input(files, fastq_dir)
    cmd = [nextflow, "run", args.pipeline, "--fastq", str(fastq_dir.resolve()), "--out_dir", str(output.resolve())]
    if args.config:
        cmd += ["-c", str(args.config.resolve())]
    cmd += ["-resume", *shlex.split(args.nextflow_args), *shlex.split(args.vimop_args)]
    print(f"[vimop]   {' '.join(shlex.quote(c) for c in cmd)}\n[vimop]   in {vimop_dir}")
    result = subprocess.run(cmd, cwd=vimop_dir)
    if result.returncode != 0:
        raise RemoveError(f"vimop failed (exit code {result.returncode}); see {vimop_dir / '.nextflow.log'}")


# what run_vimop and nextflow create in the vimop folder
VIMOP_RUN_FILES = ("input", "output", "work", ".nextflow", ".nextflow.log*")


def discard_vimop(vimop_dir: Path) -> None:
    """Delete a vimop run made by run_vimop, and the folder if nothing else is in it."""
    for pattern in VIMOP_RUN_FILES:
        for path in vimop_dir.glob(pattern):
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            else:
                path.unlink()
    if not any(vimop_dir.iterdir()):
        vimop_dir.rmdir()
    print(f"[vimop]   discarded the vimop run in {vimop_dir}")


def bam_read_names(path: Path) -> Set[str]:
    """Names of the mapped reads in a BAM file (BGZF blocks are gzip members,
    so the standard library can read it)."""
    names = set()
    with io.BufferedReader(gzip.GzipFile(path, "rb"), buffer_size=1 << 20) as fh:
        if fh.read(4) != b"BAM\1":
            raise RemoveError(f"{path} is not a BAM file")
        (l_text,) = struct.unpack("<i", fh.read(4))
        fh.read(l_text)
        (n_ref,) = struct.unpack("<i", fh.read(4))
        for _ in range(n_ref):
            (l_name,) = struct.unpack("<i", fh.read(4))
            fh.read(l_name + 4)  # name and length
        while True:
            head = fh.read(4)
            if len(head) < 4:
                return names
            (block_size,) = struct.unpack("<i", head)
            record = fh.read(block_size)
            # refID, pos, l_read_name, mapq, bin, n_cigar_op, flag, l_seq, next refID/pos, tlen: 32 bytes
            l_read_name, flag = record[8], struct.unpack_from("<H", record, 14)[0]
            if not flag & 0x4:
                names.add(record[32:32 + l_read_name - 1].decode())


@dataclass
class Virus:
    name: str  # <sample>/<reference>
    description: str
    recovery: float  # % of the consensus positions called
    reads: Set[str] = field(default_factory=set)  # the reads vimop mapped to it, if it counts


def recovery_of(row: Dict[str, str]) -> float:
    if row.get("Coverage"):
        return float(row["Coverage"])
    length = float(row.get("ConsensusLength") or 0)  # tables of older vimop versions
    return 100 * (length - float(row.get("Ambiguous positions") or 0)) / length if length else 0.0


def found_viruses(output: Path, min_recovery: float) -> List[Virus]:
    """Every virus in the samples' tables/consensus.tsv, with the reads of its
    BAM file if its consensus reached min_recovery."""
    viruses = []
    for table in sorted(output.glob("*/tables/consensus.tsv")):
        sample = table.parent.parent
        with open(table, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                reference = row["Reference"]
                virus = Virus(f"{sample.name}/{reference}", row.get("Description", ""), recovery_of(row))
                viruses.append(virus)
                if virus.recovery < min_recovery:
                    continue
                # the files are named after the reference without its version (MZ092003.1 -> MZ092003)
                candidates = [sample / "consensus" / f"{ref}.reads.bam" for ref in (reference, re.sub(r"\.\d+$", "", reference))]
                bam = next((c for c in candidates if c.is_file()), None)
                if bam is None:
                    raise RemoveError(f"no {candidates[-1].name} in {sample / 'consensus'} for {reference}")
                virus.reads = bam_read_names(bam)
    return viruses


def remove(files: List[Path], viral: Dict[str, List[Virus]], out: Path, viral_out) -> Tuple[int, int, Counter]:
    """Write the reads vimop did not map to a virus to out (and the others to viral_out)."""
    removed: Counter = Counter()  # per virus
    n_in = n_removed = 0
    tmp = out.with_name(f".partial.{out.name}")  # keeps the .gz suffix
    tmp_viral = viral_out.with_name(f".partial.{viral_out.name}") if viral_out else None
    try:
        with open_write(tmp) as keep, (open_write(tmp_viral) if tmp_viral else nullcontext()) as viral_fh:
            for fastq in files:
                for record in read_fastq(fastq):
                    n_in += 1
                    name = read_name(record[0])
                    if name not in viral and name[-2:] in ("/1", "/2"):
                        name = name[:-2]
                    if name in viral:
                        n_removed += 1
                        for virus in viral[name]:
                            removed[virus.name] += 1
                        if viral_fh:
                            write_record(viral_fh, record)
                    else:
                        write_record(keep, record)
        os.replace(tmp, out)
        if tmp_viral:
            os.replace(tmp_viral, viral_out)
    finally:
        tmp.unlink(missing_ok=True)
        if tmp_viral:
            tmp_viral.unlink(missing_ok=True)
    return n_in, n_removed, removed


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fastq", nargs="+",
                        help="FASTQ files, folders of them, or dataset ids (folders in data/raw_data/) of one run")
    parser.add_argument("-o", "--out", type=Path,
                        help=f"the run without its viral reads (default: data/output/background_fastqs/<name>{SUFFIX}.fastq.gz)")
    parser.add_argument("--min-recovery", type=float, default=50.0, metavar="PERCENT",
                        help="only remove the reads of viruses whose consensus reached this recovery "
                             "(default: %(default)g; 0 = every virus vimop reported)")
    parser.add_argument("-c", "--config", type=Path,
                        help="optional nextflow config for vimop, passed as -c (default: none, vimop's defaults)")
    parser.add_argument("--vimop-dir", type=Path,
                        help="where vimop's input, output, work folder and log go (default: output/vimop/<output name>)")
    parser.add_argument("--discard-vimop", action="store_true",
                        help="delete the vimop run (input, output, work folder, logs) once the cleaned run is "
                             "written (default: keep it)")
    parser.add_argument("--vimop-output", type=Path, metavar="DIR",
                        help="use this output of an earlier vimop run of the same reads instead of running vimop")
    parser.add_argument("--vimop-args", default="", metavar="ARGS", help='extra vimop parameters, e.g. "--targets LASV"')
    parser.add_argument("--nextflow-args", default="", metavar="ARGS",
                        help='extra nextflow options, e.g. "-profile docker -r v1.1.0"')
    parser.add_argument("--pipeline", default=VIMOP, help="the vimop pipeline to run (default: %(default)s)")
    parser.add_argument("--nextflow", default="nextflow", help="the nextflow executable (default: %(default)s)")
    parser.add_argument("--viral-out", type=Path, metavar="FASTQ", help="also write the removed reads here")
    parser.add_argument("--force", action="store_true", help="overwrite existing output files")
    parser.add_argument("--data-dir", type=Path, default=default_data_dir(),
                        help="workbench data directory (default: %(default)s)")
    args = parser.parse_args(argv)

    try:
        return run(args)
    except (RemoveError, DatasetError, FastqError) as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1


def register_output(out: Path, source: Path, min_recovery: float, data_dir: Path) -> None:
    """Add a new dataset folder in data/raw_data/ to own_datasets.tsv."""
    raw = raw_data_dir(data_dir)
    rel = inside(out, raw)
    if not rel:
        return
    dataset = rel.split("/")[0]
    if dataset in load_datasets(data_dir):
        return
    from_dataset = inside(source, raw)
    origin = from_dataset.split("/")[0] if from_dataset else os.path.abspath(source)
    notes = f"{origin} without the reads of viruses vimop found (min recovery {min_recovery:g}%)"
    table = add_own_dataset(data_dir, dataset, notes=notes)
    print(f"[dataset] {dataset}: added to {table}")


def run(args) -> int:
    if not 0 <= args.min_recovery <= 100:
        raise RemoveError("--min-recovery must be between 0 and 100")
    sources = input_paths(args.fastq, args.data_dir)
    files = input_fastqs(sources)
    out = args.out or default_out(sources[0], args.data_dir)
    check_output(out, files, args.data_dir)
    name = strip_fastq_suffix(out.name)
    report = out.with_name(f"{name}.removed_reads.tsv")
    for path in (out, args.viral_out, report):
        if path and path.exists() and not args.force:
            raise RemoveError(f"{path} already exists (use --force to overwrite)")

    vimop_dir = args.vimop_dir or default_output_dir() / "vimop" / name
    if args.vimop_output:
        output = args.vimop_output
        if not output.is_dir():
            raise RemoveError(f"vimop output not found: {output}")
    else:
        vimop_dir.mkdir(parents=True, exist_ok=True)
        output = vimop_dir / "output"
        run_vimop(args, files, vimop_dir, name, output)

    viruses = found_viruses(output, args.min_recovery)
    used = [v for v in viruses if v.recovery >= args.min_recovery]
    print(f"[viral]   vimop reports {len(viruses)} virus(es) in {output}, {len(used)} with at least "
          f"{args.min_recovery:g}% recovery")
    viral: Dict[str, List[Virus]] = {}
    for virus in used:
        for read in virus.reads:
            viral.setdefault(read, []).append(virus)

    out.parent.mkdir(parents=True, exist_ok=True)
    if args.viral_out:
        args.viral_out.parent.mkdir(parents=True, exist_ok=True)
    n_in, n_removed, removed = remove(files, viral, out, args.viral_out)

    with open(report, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["virus", "description", "recovery", "removed_reads"])
        for v in viruses:
            writer.writerow([v.name, v.description, f"{v.recovery:.2f}", removed[v.name]])
    for v in viruses:
        if v.recovery < args.min_recovery:
            note = f"kept (below {args.min_recovery:g}%)"
        else:
            note = f"{removed[v.name]} of the {len(v.reads)} reads vimop mapped to it removed"
            if v.reads and not removed[v.name]:
                note += " -- none of them is in the input: is the vimop output from these reads?"
        print(f"[virus]   {v.name:<32} recovery {v.recovery:6.2f}%  {note}  {v.description[:50]}")
    print(f"[done]    {out}: {n_in - n_removed} of {n_in} reads kept, {n_removed} viral reads removed")
    print(f"[done]    {report}")
    if args.viral_out:
        print(f"[done]    {args.viral_out}: the {n_removed} removed reads")
    register_output(out, sources[0], args.min_recovery, args.data_dir)
    if args.discard_vimop:
        if args.vimop_output:
            print(f"[vimop]   kept {args.vimop_output}: --discard-vimop only deletes vimop runs made here")
        else:
            discard_vimop(vimop_dir)
    else:
        print(f"[vimop]   the vimop run is kept in {vimop_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
