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
cleaned run goes to output/background_fastqs/<name>_no_viral.fastq.gz,
<name> being the first input's dataset or file name, unless -o says otherwise
(written into a new folder in data/raw_data/, it is added to
data/own_datasets.tsv). The reads removed per virus are listed next to the
output in <output name>.removed_reads.tsv.

--keep-viral turns the removed reads into virus organisms. They are written to
output/vimop_viral_fastqs/<name>_viral.fastq.gz (--viral-out), which is added
to data/own_datasets.tsv, and every virus above --min-recovery is extracted
from it with add_from_ref.py, mapped to the reference genome vimop used, which
is copied to data/references/<accession>.fasta. A curated virus becomes one
organism per segment, from the reference vimop marks as best. Organisms are
named <handle>_<accession>_<name>: the handle is vimop's label for a curated
virus, with its segment if it has several (COVID_OX637002_BOV-6760,
LASV_L_MG812631_BOV-6760), and else the first word of vimop's organism name
(Paenibacillus_MZ092003_BOV-6760). The sample's consensus is kept as
data/references/<organism>.consensus.fasta.

vimop's input, output, nextflow work folder and log go to
output/vimop/<output name>/ (--vimop-dir) and are kept, for vimop's report and
to reuse them; --discard-vimop deletes them once the cleaned run is written
(a failed run is always kept). --vimop-output uses the output of an earlier
vimop run of the same reads instead of running vimop; every sample in it
counts, and it is never deleted.

Examples:
    ./remove_viral_reads.py BOV-6760
    ./remove_viral_reads.py BOV-6760 --discard-vimop
    ./remove_viral_reads.py BOV-6760 --keep-viral
    ./remove_viral_reads.py fastq_pass/barcode05 -o clean.fastq.gz --vimop-args "--targets LASV"
    ./remove_viral_reads.py run.fastq.gz --min-recovery 0 --vimop-output output/vimop/run_no_viral/output
"""

from __future__ import annotations

import argparse
import csv
import filecmp
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
from typing import Dict, List, Optional, Set, Tuple

import add_from_ref
from workbench.datasets import DatasetError, add_own_dataset, fastqs_at, load_datasets
from workbench.fastq import FastqError, open_write, read_fastq, read_name, strip_fastq_suffix, write_record
from workbench.organisms import RegistryError, load_organisms
from workbench.paths import (
    background_fastqs_dir, default_data_dir, default_output_dir, inside, raw_data_dir, reads_dir,
    vimop_viral_fastqs_dir,
)

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
    """output/background_fastqs/<name>_no_viral.fastq.gz, <name> being the
    dataset the first input belongs to, or else its file name."""
    rel = inside(first, raw_data_dir(data_dir))
    name = rel.split("/")[0] if rel else strip_fastq_suffix(Path(os.path.abspath(first)).name)
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-") or "run"
    return background_fastqs_dir(default_output_dir()) / f"{name}{SUFFIX}.fastq.gz"


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
    sample: Path  # vimop's folder for the sample
    curated: bool = False
    label: str = ""  # vimop's organism label for curated viruses, e.g. LASV
    segment: str = ""  # e.g. L, S, Unsegmented
    organism: str = ""  # vimop's organism name, e.g. Mammarenavirus lassaense
    best: bool = False  # vimop's pick for its label and segment
    target: str = ""  # the name of its files in consensus/, e.g. MG812631
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
                virus = Virus(f"{sample.name}/{reference}", row.get("Description", ""), recovery_of(row), sample,
                              curated=row.get("Curated") == "True", label=row.get("Organism Label", ""),
                              segment=row.get("Segment", ""), best=row.get("IsBest") == "True",
                              organism=row.get("Organism", ""))
                viruses.append(virus)
                if virus.recovery < min_recovery:
                    continue
                # the files are named after the reference without its version (MZ092003.1 -> MZ092003)
                candidates = [sample / "consensus" / f"{ref}.reads.bam" for ref in (reference, re.sub(r"\.\d+$", "", reference))]
                bam = next((c for c in candidates if c.is_file()), None)
                if bam is None:
                    raise RemoveError(f"no {candidates[-1].name} in {sample / 'consensus'} for {reference}")
                virus.target = bam.name[:-len(".reads.bam")]
                virus.reads = bam_read_names(bam)
    return viruses


@dataclass
class ViralOrganism:
    virus: Virus
    organism_id: str
    reference: Path  # data/references/<accession>.fasta
    consensus: Path  # data/references/<organism>.consensus.fasta


def organism_id(virus: Virus, name: str) -> Optional[str]:
    """<handle>_<accession>_<name>, the handle being vimop's label (with the
    segment, LASV_L) for a curated virus and else the first word of its organism
    name; None for a curated reference vimop did not pick as best."""
    if virus.curated:
        if not virus.best:
            return None
        handle = virus.label if virus.segment in ("", "Unknown", "Unsegmented") else f"{virus.label}_{virus.segment}"
    else:
        handle = (virus.organism.split() or [""])[0]
    joined = "_".join(part for part in (handle, virus.target, name) if part)
    return re.sub(r"[^A-Za-z0-9._-]+", "_", joined).strip("._-")


def plan_organisms(viruses: List[Virus], name: str, source: str, data_dir: Path, force: bool) -> List[ViralOrganism]:
    """The organisms --keep-viral makes; checked before anything is written."""
    registered = {o.organism_id: o for o in load_organisms(data_dir)}
    references = data_dir / "references"
    plan: List[ViralOrganism] = []
    for virus in viruses:
        oid = organism_id(virus, name)
        if oid is None:
            continue
        if oid in {o.organism_id for o in plan}:  # the same reference in another sample of the vimop output
            oid = f"{oid}_{virus.sample.name}"
        other = registered.get(oid)
        if other is not None and other.source_dataset != source:
            raise RemoveError(f"organism_id {oid} is already used by {other.filename} (from {other.source_dataset})")
        organism = ViralOrganism(virus, oid, references / f"{virus.target}.fasta", references / f"{oid}.consensus.fasta")
        for path in (reads_dir(data_dir, "virus") / f"{oid}.fastq.gz", organism.consensus):
            if path.exists() and not force:
                raise RemoveError(f"{path} already exists (use --force to overwrite)")
        plan.append(organism)
    return plan


def viral_dataset(viral_out: Path, data_dir: Path) -> Tuple[str, bool]:
    """The own dataset id of the viral reads file, and whether it is registered already."""
    datasets = load_datasets(data_dir)
    for dataset in datasets.values():
        if dataset.path and os.path.realpath(dataset.location(data_dir)) == os.path.realpath(viral_out):
            return dataset.dataset_id, True
    return add_from_ref.new_dataset_id(viral_out, datasets, data_dir), False


def install_references(organism: ViralOrganism) -> None:
    """Copy the reference genome vimop used, and the sample's consensus, to data/references/."""
    consensus_dir = organism.virus.sample / "consensus"
    reference = consensus_dir / f"{organism.virus.target}.reference.fasta"
    organism.reference.parent.mkdir(parents=True, exist_ok=True)
    if organism.reference.exists():
        if not filecmp.cmp(reference, organism.reference, shallow=False):
            raise RemoveError(f"{organism.reference} exists with other sequences than {reference}")
    else:
        shutil.copyfile(reference, organism.reference)
    consensus = consensus_dir / f"{organism.virus.target}.consensus.fasta"
    if consensus.is_file():
        with open(consensus) as fh:
            text = fh.read()
        # vimop names every consensus sequence "consensus": name it after the organism
        text = re.sub(r"^>consensus\b", f">{organism.organism_id}_consensus", text, flags=re.M)
        organism.consensus.write_text(text)


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
                        help=f"the run without its viral reads (default: output/background_fastqs/<name>{SUFFIX}.fastq.gz)")
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
    parser.add_argument("--keep-viral", action="store_true",
                        help="turn the removed reads into virus organisms with add_from_ref.py")
    parser.add_argument("--viral-out", type=Path, metavar="FASTQ",
                        help="also write the removed reads here (with --keep-viral default: "
                             "output/vimop_viral_fastqs/<name>_viral.fastq.gz)")
    parser.add_argument("--force", action="store_true", help="overwrite existing output files")
    parser.add_argument("--data-dir", type=Path, default=default_data_dir(),
                        help="workbench data directory (default: %(default)s)")
    args = parser.parse_args(argv)

    try:
        return run(args)
    except (RemoveError, DatasetError, FastqError, RegistryError) as err:
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


def keep_viral(plan: List[ViralOrganism], viral_out: Path, viral_id: str, registered: bool,
               source: Path, args) -> None:
    """Register the viral reads as an own dataset and extract each organism from them."""
    if not registered:
        from_dataset = inside(source, raw_data_dir(args.data_dir))
        origin = from_dataset.split("/")[0] if from_dataset else os.path.abspath(source)
        notes = f"reads of {origin} that vimop mapped to viruses with at least {args.min_recovery:g}% recovery"
        table = add_own_dataset(args.data_dir, viral_id, path=os.path.abspath(viral_out), notes=notes)
        print(f"[dataset] {viral_id}: added to {table}")
    if not plan:
        print(f"[viral]   no virus with at least {args.min_recovery:g}% recovery to keep as an organism")
    for organism in plan:
        install_references(organism)
        print(f"[viral]   {organism.virus.name} -> organism {organism.organism_id}")
        argv = ["--data-dir", str(args.data_dir), "--category", "virus", str(organism.reference),
                organism.organism_id, f"{organism.organism_id}.fastq.gz", viral_id] + (["--force"] if args.force else [])
        if add_from_ref.main(argv) != 0:
            raise RemoveError(f"add_from_ref.py could not extract {organism.organism_id}")


def run(args) -> int:
    if not 0 <= args.min_recovery <= 100:
        raise RemoveError("--min-recovery must be between 0 and 100")
    sources = input_paths(args.fastq, args.data_dir)
    files = input_fastqs(sources)
    out = args.out or default_out(sources[0], args.data_dir)
    check_output(out, files, args.data_dir)
    name = strip_fastq_suffix(out.name)
    base = name[:-len(SUFFIX)] if name.endswith(SUFFIX) and len(name) > len(SUFFIX) else name
    report = out.with_name(f"{name}.removed_reads.tsv")
    viral_out = args.viral_out
    if args.keep_viral:
        if shutil.which("minimap2") is None:
            raise RemoveError("minimap2 not found on PATH, needed by --keep-viral (conda env create -f environment.yml)")
        viral_out = viral_out or vimop_viral_fastqs_dir(default_output_dir()) / f"{base}_viral.fastq.gz"
        if inside(viral_out, raw_data_dir(args.data_dir)) is not None:
            raise RemoveError(f"with --keep-viral the viral reads cannot go into {raw_data_dir(args.data_dir)}/")
        viral_id, viral_registered = viral_dataset(viral_out, args.data_dir)
    for path in (out, viral_out, report):
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
    plan = plan_organisms(used, base, viral_id, args.data_dir, args.force) if args.keep_viral else []
    organism_of = {o.virus.name: o.organism_id for o in plan}

    out.parent.mkdir(parents=True, exist_ok=True)
    if viral_out:
        viral_out.parent.mkdir(parents=True, exist_ok=True)
    n_in, n_removed, removed = remove(files, viral, out, viral_out)

    with open(report, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["virus", "description", "recovery", "removed_reads", "organism"])
        for v in viruses:
            writer.writerow([v.name, v.description, f"{v.recovery:.2f}", removed[v.name], organism_of.get(v.name, "")])
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
    if viral_out:
        print(f"[done]    {viral_out}: the {n_removed} removed reads")
    register_output(out, sources[0], args.min_recovery, args.data_dir)
    if args.keep_viral:
        keep_viral(plan, viral_out, viral_id, viral_registered, sources[0], args)
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
