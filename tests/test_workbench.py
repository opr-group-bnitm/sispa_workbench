import csv
import gzip
import os
import random
import shutil
import sys
import textwrap
from array import array
from pathlib import Path

import pytest

import add_from_ref
import create_sispa_run
import remove_viral_reads
from workbench import mapping, paths
from workbench.alignments import AlignmentWriter, ReadAlignments, alignments_path, read_alignments
from workbench.coverage import DepthCounter, combine, select_for_depth, summarize
from workbench.datasets import DatasetError, fastqs_at, load_datasets
from workbench.fastq import FastqError, fastq_stats, read_fastq, strip_fastq_suffix
from workbench.organisms import COLUMNS, Organism, RegistryError, load_organisms, prune_missing, register


def write_fastq(path: Path, reads, gz=None):
    """reads: list of (name, sequence)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(f"@{name} some comment\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in reads)
    if gz if gz is not None else path.name.endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def names_in(path: Path):
    return [header.split()[0].decode() for header, _, _ in read_fastq(path)]


def organisms(data_dir):
    return {o.organism_id: o for o in load_organisms(data_dir)}


def add_organism(data_dir, filename, reads, organism_id=None):
    """Register a FASTQ the way add_from_ref.py would, without mapping: every
    read aligns from the start of a 200 bp reference."""
    path = write_fastq(data_dir / filename, reads)
    alignments_path(data_dir, filename).parent.mkdir(parents=True, exist_ok=True)
    alignments = AlignmentWriter(alignments_path(data_dir, filename), {"ref": 200})
    for i, (_, seq) in enumerate(reads):
        alignments.add(i, [("ref", 0, min(len(seq), 200))])
    alignments.close(len(reads))
    s = fastq_stats(path)
    register(data_dir, Organism(
        organism_id or strip_fastq_suffix(path.name), filename, round(s.avg_read_length, 2),
        s.max_read_length, s.min_read_length, s.n_reads, 100.0, 1, 3, "references/ref.fasta", "PUB",
    ), log=lambda _: None)
    return path


@pytest.fixture(autouse=True)
def public_table(tmp_path, monkeypatch):
    """Tests use their own public dataset table, not the repository's."""
    table = tmp_path / "public.tsv"
    table.write_text(
        "dataset_id\tinclude\tstudy_accession\trepository\tdataset_accession\n"
        "PUB\tfalse\tPRJ1\tENA\t\n"
        "PUB-EMPTY\tfalse\tPRJ2\tENA\t\n"
    )
    monkeypatch.setattr(paths, "PUBLIC_DATASETS_TSV", table)
    return table


@pytest.fixture
def data_dir(tmp_path):
    d = tmp_path / "data"
    add_organism(d, "virus_reads/COVID.fastq.gz", [(f"covid_{i}", "ACGT" * (i + 1)) for i in range(50)])
    add_organism(d, "virus_reads/LASV.fq", [(f"lasv_{i}", "GGCC" * 3) for i in range(20)])
    add_organism(d, "host_reads/human.fastq", [(f"human_{i}", "T" * 100) for i in range(300)])
    return d


# ------------------------------------------------------------------ fastq ---

def test_stats(tmp_path):
    path = write_fastq(tmp_path / "x.fastq.gz", [("a", "ACG"), ("b", "ACGTACGTAC"), ("c", "AC")])
    stats = fastq_stats(path)
    assert (stats.n_reads, stats.min_read_length, stats.max_read_length) == (3, 2, 10)
    assert stats.avg_read_length == pytest.approx(5.0)


def test_gzip_detected_by_content_not_name(tmp_path):
    path = write_fastq(tmp_path / "x.fastq", [("a", "ACGT")], gz=True)
    assert fastq_stats(path).n_reads == 1


def test_malformed_fastq(tmp_path):
    path = tmp_path / "bad.fastq"
    path.write_text("@a\nACGT\n+\nII\n")
    with pytest.raises(FastqError):
        fastq_stats(path)


# --------------------------------------------------------------- coverage ---

def test_depth_summary():
    gap = summarize(20, [0, 5], [10, 15])
    assert (gap.covered, gap.min_depth, gap.max_depth, gap.breadth, gap.mean_depth) == (15, 0, 2, 75.0, 1.0)
    full = summarize(10, [0, 0, 4], [10, 6, 10])
    assert (full.covered, full.min_depth, full.max_depth) == (10, 2, 3)
    adjacent = summarize(10, [0, 5], [5, 10])
    assert (adjacent.covered, adjacent.min_depth, adjacent.max_depth) == (10, 1, 1)
    empty = summarize(10, [], [])
    assert (empty.covered, empty.min_depth, empty.max_depth, empty.breadth) == (0, 0, 0, 0.0)


def test_depth_summary_matches_per_base_counting():
    rng = random.Random(1)
    for _ in range(200):
        length = rng.randint(1, 60)
        spans = [sorted(rng.sample(range(length + 1), 2)) for _ in range(rng.randint(0, 8))]
        spans = [(s, e) for s, e in spans if e > s]
        depth = [sum(s <= p < e for s, e in spans) for p in range(length)]
        got = summarize(length, [s for s, _ in spans], [e for _, e in spans])
        assert (got.covered, got.min_depth, got.max_depth) == (sum(d > 0 for d in depth), min(depth), max(depth))
        assert got.mean_depth == pytest.approx(sum(depth) / length)


def test_depth_over_several_sequences():
    counter = DepthCounter()
    counter.add("S", 0, 50)
    per_sequence = counter.summarize({"S": 100, "L": 200})
    assert per_sequence["L"].covered == 0 and per_sequence["S"].breadth == 50.0
    total = combine(list(per_sequence.values()))
    assert (total.length, total.covered, total.min_depth, total.max_depth) == (300, 50, 0, 1)
    assert total.breadth == pytest.approx(100 / 6)


def make_alignments(lengths, reads):
    """reads: per read, its alignments as (sequence, start, end)."""
    names = list(lengths)
    columns = [array("q") for _ in range(4)]
    for read, alignments in enumerate(reads):
        for name, start, end in alignments:
            for column, value in zip(columns, (read, names.index(name), start, end)):
                column.append(value)
    return ReadAlignments(dict(lengths), len(reads), *columns)


def depth_of(reads, chosen, name, length):
    return [sum(s <= p < e for r in chosen for n, s, e in reads[r] if n == name) for p in range(length)]


def test_select_for_depth_meets_the_target_everywhere():
    rng = random.Random(5)
    for trial in range(300):
        lengths = {f"seg{k}": rng.randint(1, 50) for k in range(rng.randint(1, 3))}
        reads = []
        for _ in range(rng.randint(1, 25)):
            alignments = []
            for _ in range(rng.choice((1, 1, 1, 2))):  # some reads have a supplementary alignment
                name = rng.choice(list(lengths))
                start = rng.randrange(lengths[name])
                alignments.append((name, start, rng.randint(start + 1, lengths[name])))
            reads.append(alignments)
        min_depth = rng.randint(1, 6)
        selection = select_for_depth(make_alignments(lengths, reads), min_depth, random.Random(trial))
        assert selection.reads == sorted(set(selection.reads))
        full = short = 0
        for name, length in lengths.items():
            available = depth_of(reads, range(len(reads)), name, length)
            chosen = depth_of(reads, selection.reads, name, length)
            assert all(c >= min(min_depth, a) for c, a in zip(chosen, available))
            full += sum(a >= min_depth for a in available)
            short += sum(0 < a < min_depth for a in available)
        assert (selection.full_bases, selection.short_bases, selection.total_bases) == (
            full, short, sum(lengths.values()))


def test_select_for_depth_takes_no_more_than_needed():
    # 100 reads over the whole genome and 10 over its first half
    reads = [[("g", 0, 100)]] * 100 + [[("g", 0, 50)]] * 10
    for seed in range(5):
        selection = select_for_depth(make_alignments({"g": 100}, reads), 5, random.Random(seed))
        whole = [r for r in selection.reads if r < 100]
        assert len(whole) == 5  # the second half needs 5 whole reads, which also serve the first half
        assert selection.reads != select_for_depth(make_alignments({"g": 100}, reads), 5,
                                                   random.Random(seed + 10)).reads


def test_select_for_depth_uses_every_read_where_depth_is_short():
    reads = [[("g", 0, 60)]] * 3 + [[("g", 40, 100)]] * 20
    selection = select_for_depth(make_alignments({"g": 100}, reads), 10, random.Random(1))
    assert {0, 1, 2} <= set(selection.reads)  # 0-40 has only 3 reads
    assert (selection.full_bases, selection.short_bases) == (60, 40)


def test_parse_alignment_strands_and_clips():
    forward = mapping.parse_alignment([b"r", b"0", b"ref", b"11", b"60", b"5S30M10I2D3S"])
    assert (forward.start, forward.end, forward.query_length) == (10, 42, 48)
    assert (forward.query_start, forward.query_end) == (5, 45)
    reverse = mapping.parse_alignment([b"r", b"16", b"ref", b"11", b"60", b"5S30M10I2D3S"])
    assert (reverse.query_start, reverse.query_end) == (3, 43)
    supplementary = mapping.parse_alignment([b"r", b"2048", b"ref", b"61", b"60", b"20H20M"])
    assert (supplementary.start, supplementary.end) == (60, 80)
    assert (supplementary.query_start, supplementary.query_end, supplementary.query_length) == (20, 40, 40)


# --------------------------------------------------------------- registry ---

def test_registry_columns(data_dir):
    header = (data_dir / "organisms.tsv").read_text().splitlines()[0].split("\t")
    assert header == [
        "organism_id", "filename", "avg_read_length", "max_read_length", "min_read_length", "n_reads",
        "breadth_coverage", "min_depth", "max_depth", "reference", "source_dataset",
    ] == COLUMNS
    covid = organisms(data_dir)["COVID"]
    assert (covid.n_reads, covid.min_read_length, covid.max_read_length) == (50, 4, 200)
    assert covid.avg_read_length == pytest.approx(102.0)


def test_register_replaces_the_row_of_the_same_file(data_dir):
    add_organism(data_dir, "virus_reads/LASV.fq", [("x", "ACGT")], organism_id="Lassa")
    orgs = organisms(data_dir)
    assert "LASV" not in orgs and orgs["Lassa"].n_reads == 1
    with pytest.raises(RegistryError, match="already used by virus_reads/COVID.fastq.gz"):
        add_organism(data_dir, "virus_reads/other.fq", [("x", "ACGT")], organism_id="COVID")


def test_prune_missing(data_dir):
    (data_dir / "virus_reads" / "COVID.fastq.gz").unlink()
    removed = []
    assert {o.organism_id for o in prune_missing(data_dir, log=removed.append)} == {"LASV", "human"}
    assert "COVID" in removed[0] and "COVID" not in organisms(data_dir)


def test_old_organisms_table_is_rejected(tmp_path):
    (tmp_path / "organisms.tsv").write_text(
        "organism_id\tfilename\tavg_read_length\tmax_read_length\tmin_read_length\tn_reads\n"
    )
    with pytest.raises(RegistryError, match="older version"):
        load_organisms(tmp_path)


# --------------------------------------------------------------- datasets ---

OWN_HEADER = "dataset_id\tsample_type\torganism_target\tnotes\tpath\n"


def test_own_table_is_created_when_needed(tmp_path):
    data_dir = tmp_path / "data"
    assert set(load_datasets(data_dir)) == {"PUB", "PUB-EMPTY"}
    assert (data_dir / "own_datasets.tsv").read_text() == OWN_HEADER


def test_older_own_table_gets_the_path_column(tmp_path):
    table = tmp_path / "own_datasets.tsv"
    table.write_text("# my runs\ndataset_id\tsample_type\torganism_target\tnotes\nOLD\tswab\tDENV\tfirst run\n")
    datasets = load_datasets(tmp_path)
    assert datasets["OLD"].path == "" and not datasets["OLD"].is_public
    assert table.read_text() == (
        "# my runs\ndataset_id\tsample_type\torganism_target\tnotes\tpath\nOLD\tswab\tDENV\tfirst run\n"
    )


def test_dataset_ids_must_be_unique_across_tables(tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "own_datasets.tsv").write_text("dataset_id\tsample_type\torganism_target\tnotes\nPUB\t\t\t\n")
    with pytest.raises(DatasetError, match="in both"):
        load_datasets(data_dir)


def test_fastqs_skip_hidden_and_partial_files(tmp_path):
    folder = tmp_path / "raw_data" / "PUB"
    for name in ("a.fastq.gz", "barcode01/b.fq", ".done/c.fastq", ".tmp_run/d.fastq", "e.fastq.gz.part", "runs.tsv"):
        (folder / name).parent.mkdir(parents=True, exist_ok=True)
        (folder / name).write_text("")
    assert [p.relative_to(folder).as_posix() for p in fastqs_at(folder)] == ["a.fastq.gz", "barcode01/b.fq"]
    assert fastqs_at(folder / "a.fastq.gz") == [folder / "a.fastq.gz"]
    assert fastqs_at(tmp_path / "NOT-THERE") == []


def test_fastqs_follow_symlinks(tmp_path):
    elsewhere = tmp_path / "elsewhere"
    for name in ("a.fastq.gz", "barcode01/b.fastq.gz"):
        (elsewhere / name).parent.mkdir(parents=True, exist_ok=True)
        (elsewhere / name).write_text("")
    raw = tmp_path / "raw_data"
    raw.mkdir()
    (raw / "LINKED").symlink_to(elsewhere)  # a whole folder
    (raw / "PARTS").mkdir()
    (raw / "PARTS" / "a.fastq.gz").symlink_to(elsewhere / "a.fastq.gz")  # a file ...
    (raw / "PARTS" / "same.fastq.gz").symlink_to(elsewhere / "a.fastq.gz")  # ... linked twice
    (raw / "PARTS" / "barcode01").symlink_to(elsewhere / "barcode01")  # a sub-folder
    (raw / "PARTS" / "loop").symlink_to(raw / "PARTS")
    rel = lambda folder: [p.relative_to(folder).as_posix() for p in fastqs_at(folder)]
    assert rel(raw / "LINKED") == ["a.fastq.gz", "barcode01/b.fastq.gz"]
    assert rel(raw / "PARTS") == ["a.fastq.gz", "barcode01/b.fastq.gz"]
    (raw / "PARTS" / "barcode02").symlink_to(tmp_path / "unmounted" / "barcode02")
    with pytest.raises(DatasetError, match="is the drive mounted"):
        fastqs_at(raw / "PARTS")


# ------------------------------------------------------- create_sispa_run ---

def run_mix(data_dir, tmp_path, *args):
    return create_sispa_run.main([
        *args, "--data-dir", str(data_dir),
        "--fastq-dir", str(tmp_path / "out" / "fastqs"),
        "--composition-dir", str(tmp_path / "out" / "compositions"),
    ])


def test_mix_from_list(data_dir, tmp_path):
    assert run_mix(data_dir, tmp_path, "COVID", "10,", "human", "200,", "LASV", "all", "--seed", "1") == 0
    name = "COVID-10_human-200_LASV-all"
    reads = names_in(tmp_path / "out" / "fastqs" / f"{name}.fastq.gz")
    assert len(reads) == len(set(reads)) == 230
    assert sum(r.startswith("covid_") for r in reads) == 10
    assert sum(r.startswith("human_") for r in reads) == 200
    assert sorted(r for r in reads if r.startswith("lasv_")) == sorted(f"lasv_{i}" for i in range(20))

    with open(tmp_path / "out" / "compositions" / f"{name}.csv") as fh:
        rows = list(csv.DictReader(fh))
    # every read aligns from position 0 of a 200 bp reference: COVID read i is 4 * (i + 1) long
    covid = [4 * (int(r.split("_")[1]) + 1) for r in reads if r.startswith("covid_")]
    assert rows == [
        {"organism_id": "COVID", "n_reads": "10", "min_depth": str(int(max(covid) == 200)),
         "max_depth": "10", "mean_depth": f"{sum(covid) / 200:.2f}"},
        {"organism_id": "human", "n_reads": "200", "min_depth": "0", "max_depth": "200", "mean_depth": "100.00"},
        {"organism_id": "LASV", "n_reads": "20", "min_depth": "0", "max_depth": "20", "mean_depth": "1.20"},
    ]


def test_the_composition_template_is_valid():
    template = paths.REPO_ROOT / "data" / "input_compositions" / "template.csv"
    assert create_sispa_run.parse_table(template) == [
        create_sispa_run.Request("COVID", n_reads=5000),
        create_sispa_run.Request("DENV2", min_depth=20),
        create_sispa_run.Request("human", n_reads=70000),
        create_sispa_run.Request("LASV"),
    ]


def test_mix_single_quoted_list_and_name(data_dir, tmp_path):
    assert run_mix(data_dir, tmp_path, "COVID 5, LASV 5", "--name", "tiny", "--seed", "3") == 0
    assert len(names_in(tmp_path / "out" / "fastqs" / "tiny.fastq.gz")) == 10
    assert (tmp_path / "out" / "compositions" / "tiny.csv").exists()


def test_mix_from_csv_uses_its_name(data_dir, tmp_path):
    table = tmp_path / "spike_in.csv"
    table.write_text("organism_id,n_reads\nCOVID,all\nhuman,7\n")
    assert run_mix(data_dir, tmp_path, str(table)) == 0
    assert len(names_in(tmp_path / "out" / "fastqs" / "spike_in.fastq.gz")) == 57
    # all 50 COVID reads: the last position is covered by the 200 bp read alone
    assert (tmp_path / "out" / "compositions" / "spike_in.csv").read_text() == (
        "organism_id,n_reads,min_depth,max_depth,mean_depth\n"
        "COVID,50,1,50,25.50\n"
        "human,7,0,7,3.50\n"
    )


def test_mix_is_reproducible_and_interleaved(data_dir, tmp_path):
    run_mix(data_dir, tmp_path, "COVID 50, human 50", "--name", "a", "--seed", "7")
    run_mix(data_dir, tmp_path, "COVID 50, human 50", "--name", "b", "--seed", "7")
    run_mix(data_dir, tmp_path, "COVID 50, human 50", "--name", "c", "--seed", "8")
    fastqs = tmp_path / "out" / "fastqs"
    a, b, c = (names_in(fastqs / f"{n}.fastq.gz") for n in "abc")
    assert a == b and a != c
    labels = [r.split("_")[0] for r in a]
    switches = sum(x != y for x, y in zip(labels, labels[1:]))
    assert switches > 20  # organisms are interleaved, not concatenated blocks


def test_mix_subsample_is_uniform(tmp_path):
    data_dir = tmp_path / "data"
    add_organism(data_dir, "virus_reads/V.fastq", [(f"r{i}", "ACGT") for i in range(10)])
    hits = [0] * 10
    for seed in range(300):
        run_mix(data_dir, tmp_path, "V 3", "--name", f"s{seed}", "--seed", str(seed))
        for name in names_in(tmp_path / "out" / "fastqs" / f"s{seed}.fastq.gz"):
            hits[int(name[1:])] += 1
    assert min(hits) > 50 and max(hits) < 130  # expected 90 each


@pytest.mark.parametrize("composition, message", [
    (["COVIDD 5"], "did you mean COVID"),
    (["COVID 51"], "holds only 50"),
    (["COVID 5, COVID 6"], "more than once"),
    (["COVID five"], "positive number"),
    (["COVID 5, human"], "expected pairs"),
    (["missing.csv"], "not found"),
])
def test_mix_errors(data_dir, tmp_path, capsys, composition, message):
    assert run_mix(data_dir, tmp_path, *composition) == 1
    assert message in capsys.readouterr().err


def test_mix_only_uses_files_made_by_add_from_ref(data_dir, tmp_path, capsys):
    write_fastq(data_dir / "virus_reads" / "ROGUE.fastq", [("x", "ACGT")])
    assert run_mix(data_dir, tmp_path, "ROGUE 1") == 1
    err = capsys.readouterr().err
    assert "unknown organism_id 'ROGUE'" in err and "virus_reads/ROGUE.fastq was not created by add_from_ref.py" in err


def test_mix_checks_the_files_still_match_the_registry(data_dir, tmp_path, capsys):
    (data_dir / "virus_reads" / "COVID.fastq.gz").unlink()
    assert run_mix(data_dir, tmp_path, "COVID 1") == 1
    assert "virus_reads/COVID.fastq.gz is missing" in capsys.readouterr().err
    write_fastq(data_dir / "virus_reads" / "LASV.fq", [("lasv_0", "GGCC")])
    assert run_mix(data_dir, tmp_path, "LASV all") == 1
    assert "no longer holds the 20 reads" in capsys.readouterr().err
    assert not list((tmp_path / "out" / "fastqs").glob("*LASV*"))


def test_mix_refuses_to_overwrite(data_dir, tmp_path, capsys):
    assert run_mix(data_dir, tmp_path, "LASV 2", "--seed", "1") == 0
    assert run_mix(data_dir, tmp_path, "LASV 2", "--seed", "1") == 1
    assert "already exists" in capsys.readouterr().err
    assert run_mix(data_dir, tmp_path, "LASV 2", "--seed", "1", "--force") == 0


# ----------------------------------------------------------- add_from_ref ---

# Stand-in for `minimap2 -a --sam-hit-only` that "maps" reads by name, onto the
# first reference sequence:
#   v_*     end to end at position 1,                       mapq 60
#   half_*  first half at position 1, rest soft-clipped,    mapq 60
#   low_*   end to end at position 1,                       mapq 1
#   chim_*  first half at position 1, second half as a supplementary
#           alignment at position 61,                       mapq 60
#   at<N>_* end to end at position N,                       mapq 60
#   h_*     no alignment
FAKE_MINIMAP2 = """\
    #!{python}
    import gzip, sys
    ref, fastq = sys.argv[-2], sys.argv[-1]
    lengths, name = {}, None
    for line in open(ref):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            lengths[name] = 0
        elif name:
            lengths[name] += len(line)
    print("@HD\\tVN:1.6\\tSO:unsorted")
    for name, length in lengths.items():
        print(f"@SQ\\tSN:{name}\\tLN:{length}")
    target = next(iter(lengths))
    with (gzip.open if fastq.endswith(".gz") else open)(fastq, "rt") as fh:
        lines = fh.read().splitlines()
    for header, seq in zip(lines[0::4], lines[1::4]):
        qname, n = header[1:].split()[0], len(seq)
        qname = qname[:-2] if qname.endswith(("/1", "/2")) else qname  # like minimap2
        h = n // 2
        records = {
            "v": [(0, 1, 60, f"{n}M")],
            "half": [(0, 1, 60, f"{h}M{n - h}S")],
            "low": [(0, 1, 1, f"{n}M")],
            "chim": [(0, 1, 60, f"{h}M{n - h}S"), (2048, 61, 60, f"{h}H{n - h}M")],
        }.get(qname.split("_")[0], [])
        if qname.startswith("at"):
            records = [(0, int(qname[2:].split("_")[0]), 60, f"{n}M")]
        for flag, pos, mapq, cigar in records:
            print(f"{qname}\\t{flag}\\t{target}\\t{pos}\\t{mapq}\\t{cigar}\\t*\\t0\\t0\\t{seq}\\t*")
"""

MIXED_READS = [
    ("v_1", "ACGT" * 10), ("h_1", "TTTT"), ("half_1", "ACGT" * 10), ("low_1", "ACGT" * 10),
    ("v_2/1", "ACGT"), ("h_2", "TTTT"), ("chim_1", "ACGT" * 10),
]


@pytest.fixture
def fake_minimap2(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "minimap2"
    exe.write_text(textwrap.dedent(FAKE_MINIMAP2).replace("{python}", sys.executable))
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path(shutil.which('sh')).parent}")


@pytest.fixture
def raw(data_dir):
    """data/raw_data with the public dataset PUB and an unregistered folder MINE."""
    write_fastq(data_dir / "raw_data" / "PUB" / "run1.fastq", MIXED_READS)
    write_fastq(data_dir / "raw_data" / "MINE" / "barcode01" / "reads.fastq.gz", [("v_9", "ACGT"), ("h_9", "ACGT")])
    ref = data_dir / "references" / "ref.fasta"
    ref.parent.mkdir(parents=True, exist_ok=True)
    ref.write_text(">ref\n" + "A" * 60 + "\n" + "A" * 40 + "\n")  # 100 bp
    return data_dir / "raw_data"


def add(data_dir, *args):
    return add_from_ref.main(["--data-dir", str(data_dir), *args])


def test_add_from_ref_extracts_mapped_reads(data_dir, raw, fake_minimap2):
    assert add(data_dir, "ref.fasta", "DENV2", "dengue", "PUB") == 0
    out = data_dir / "virus_reads" / "dengue.fastq.gz"
    assert names_in(out) == ["v_1", "half_1", "low_1", "v_2/1", "chim_1"]
    # records are copied unchanged, comments included
    assert next(read_fastq(out))[0] == b"v_1 some comment"
    denv = organisms(data_dir)["DENV2"]
    assert (denv.filename, denv.n_reads) == ("virus_reads/dengue.fastq.gz", 5)
    assert (denv.reference, denv.source_dataset) == ("references/ref.fasta", "PUB")
    # depth 5 on 0-3, 4 on 4-19, 2 on 20-39, the supplementary alignment 1 on 60-79
    assert (denv.breadth_coverage, denv.min_depth, denv.max_depth) == (60.0, 0, 5)

    # where each read aligns is kept in data/alignments/, and only FASTQs in the read folder
    assert alignments_path(data_dir, denv.filename) == (
        data_dir / "alignments" / "virus_reads" / "dengue.fastq.gz.alignments.tsv.gz")
    assert not [p.name for p in (data_dir / "virus_reads").iterdir() if p.name.endswith(".alignments.tsv.gz")]
    with gzip.open(alignments_path(data_dir, denv.filename), "rt") as fh:
        assert fh.read().splitlines() == [
            "#sequence\tref\t100", "read\tsequence\tstart\tend",
            "0\tref\t0\t40", "1\tref\t0\t20", "2\tref\t0\t40", "3\tref\t0\t4",
            "4\tref\t0\t20", "4\tref\t60\t80", "#reads\t5",
        ]
    assert read_alignments(alignments_path(data_dir, denv.filename)).n_reads == 5


def test_add_from_ref_filters(data_dir, raw, fake_minimap2):
    assert add(data_dir, "--category", "host", "--min-mapq", "20", "--min-aligned-fraction", "0.8",
               "ref.fasta", "strict", "strict.fq", "PUB") == 0
    # half_1 aligns over 50%, low_1 only below mapq 20; chim_1 aligns fully with its two parts
    assert names_in(data_dir / "host_reads" / "strict.fq") == ["v_1", "v_2/1", "chim_1"]
    strict = organisms(data_dir)["strict"]
    assert (strict.breadth_coverage, strict.min_depth, strict.max_depth) == (60.0, 0, 3)


def test_add_from_ref_registers_own_datasets(data_dir, raw, fake_minimap2, capsys):
    own = data_dir / "own_datasets.tsv"
    assert not own.exists()
    assert add(data_dir, "ref.fasta", "X", "x", "PUB", "MINE") == 0
    assert "MINE: added to" in capsys.readouterr().out
    assert own.read_text() == OWN_HEADER + "MINE\t\t\t\t\n"
    assert names_in(data_dir / "virus_reads" / "x.fastq.gz")[-1] == "v_9"
    assert organisms(data_dir)["X"].source_dataset == "PUB;MINE"

    # a described own dataset is used as it is
    own.write_text(OWN_HEADER + "MINE\tswab\tDENV\tmy run\t\n")
    assert add(data_dir, "ref.fasta", "Y", "y", "MINE") == 0
    assert own.read_text().count("MINE") == 1


@pytest.mark.parametrize("args, message", [
    (["ref.fasta", "COVID", "other.fastq.gz", "PUB"], "already used by virus_reads/COVID.fastq.gz"),
    (["ref.fasta", "NEW", "COVID.fastq.gz", "PUB"], "already exists"),
    (["ref.fasta", "bad id", "x.fastq", "PUB"], "may only contain"),
    (["ref.fasta", "NEW", "sub/x.fastq", "PUB"], "plain file name"),
    (["nope.fasta", "NEW", "x.fastq", "PUB"], "reference not found"),
    (["ref.fasta", "NEW", "x.fastq", "PUB-EMPTY"], "./download_datasets.sh PUB-EMPTY"),
    (["ref.fasta", "NEW", "x.fastq", "NOWHERE"], "'NOWHERE' is neither a dataset id"),
    (["ref.fasta", "NEW", "x.fastq", "no/such/reads.fastq"], "nor an existing file or folder"),
    (["ref.fasta", "NEW", "x.fastq", "PUB", "PUB"], "included twice"),
])
def test_add_from_ref_errors(data_dir, raw, fake_minimap2, capsys, args, message):
    assert add(data_dir, *args) == 1
    assert message in capsys.readouterr().err


def test_add_from_ref_nothing_mapped(data_dir, raw, fake_minimap2, capsys):
    write_fastq(raw / "HOST-ONLY" / "r.fastq", [("h_1", "ACGT")])
    assert add(data_dir, "ref.fasta", "NEW", "new.fastq.gz", "HOST-ONLY") == 1
    assert "no reads mapped" in capsys.readouterr().err
    assert not list((data_dir / "virus_reads").glob("*new*"))
    assert "NEW" not in organisms(data_dir)


def test_add_from_ref_housekeeping(data_dir, raw, fake_minimap2, capsys):
    (data_dir / "virus_reads" / "LASV.fq").unlink()
    write_fastq(data_dir / "virus_reads" / "rogue.fastq", [("x", "ACGT")])
    assert add(data_dir, "ref.fasta", "Z", "z", "PUB") == 0
    out = capsys.readouterr().out
    assert "[removed] LASV" in out and "LASV" not in organisms(data_dir)
    assert "[ignored] virus_reads/rogue.fastq was not created by add_from_ref.py" in out
    assert "rogue" not in organisms(data_dir)


def test_minimap2_args_follow_vimop(tmp_path):
    small = tmp_path / "virus.fasta"
    small.write_text(">v\nACGT\n")
    split = tmp_path / "idx"
    assert mapping.minimap2_args(small, "virus", split) == ["--secondary=no", "-k", "11", "-w", "5"]
    assert mapping.minimap2_args(small, "host", split) == ["--secondary=no"]
    big = tmp_path / "host.fasta"
    with open(big, "wb") as fh:
        fh.truncate(mapping.VIMOP_SPLIT_THRESHOLD + 1)
    assert mapping.minimap2_args(big, "host", split) == ["--secondary=no", "-I", "2G", "--split-prefix", str(split)]


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not installed")
def test_add_from_ref_with_real_minimap2(tmp_path):
    rng = random.Random(0)
    genome = "".join(rng.choice("ACGT") for _ in range(20000))
    other = "".join(rng.choice("ACGT") for _ in range(20000))
    data_dir = tmp_path / "data"
    (data_dir / "references").mkdir(parents=True)
    (data_dir / "references" / "virus.fasta").write_text(f">virus\n{genome}\n")
    # 2 kb reads every 900 bp: up to 3 reads deep, the last 900 bp uncovered
    reads = [(f"virus_{i}", genome[i * 900:i * 900 + 2000]) for i in range(20)]
    reads += [(f"other_{i}", other[i * 900:i * 900 + 2000]) for i in range(20)]
    write_fastq(data_dir / "raw_data" / "REAL" / "mixed.fastq.gz", reads)
    assert add_from_ref.main(["--data-dir", str(data_dir), "virus.fasta", "V", "V.fastq.gz", "REAL"]) == 0
    assert sorted(names_in(data_dir / "virus_reads" / "V.fastq.gz")) == sorted(n for n, _ in reads[:20])
    v = organisms(data_dir)["V"]
    assert 95.0 <= v.breadth_coverage <= 95.5 and (v.min_depth, v.max_depth) == (0, 3)


def test_add_from_ref_with_fastqs_anywhere(data_dir, raw, fake_minimap2, tmp_path, capsys):
    run = tmp_path / "sequencer" / "run1"
    write_fastq(run / "barcode05" / "part0.fastq", MIXED_READS)
    write_fastq(run / "barcode06" / "reads.fastq.gz", [("v_7", "ACGT")])
    own = data_dir / "own_datasets.tsv"

    # a folder outside data/ is added to own_datasets.tsv with its path
    assert add(data_dir, "ref.fasta", "A", "a", str(run / "barcode05")) == 0
    assert f"barcode05: added to {own} for {run / 'barcode05'}" in capsys.readouterr().out
    assert own.read_text() == OWN_HEADER + f"barcode05\t\t\t\t{run / 'barcode05'}\n"
    assert organisms(data_dir)["A"].source_dataset == "barcode05"
    assert names_in(data_dir / "virus_reads" / "a.fastq.gz") == ["v_1", "half_1", "low_1", "v_2/1", "chim_1"]

    # the same path again, by path or by id, reuses the row; a file inside it is a part of it
    assert add(data_dir, "ref.fasta", "B", "b", str(run / "barcode05") + "/") == 0
    assert add(data_dir, "ref.fasta", "C", "c", "barcode05") == 0
    assert add(data_dir, "ref.fasta", "D", "d", str(run / "barcode05" / "part0.fastq")) == 0
    assert own.read_text().splitlines()[1:] == [f"barcode05\t\t\t\t{run / 'barcode05'}"]
    assert [organisms(data_dir)[o].source_dataset for o in "BCD"] == ["barcode05", "barcode05", "barcode05/part0.fastq"]

    # a single file gets an id from its name; one part of a dataset in data/raw_data/ is used alone
    assert add(data_dir, "ref.fasta", "E", "e", str(run / "barcode06" / "reads.fastq.gz"),
               str(raw / "PUB" / "run1.fastq")) == 0
    assert organisms(data_dir)["E"].source_dataset == "reads;PUB/run1.fastq"
    assert [line.split("\t")[0] for line in own.read_text().splitlines()[1:]] == ["barcode05", "reads"]


def test_add_from_ref_path_ids_stay_unique(data_dir, raw, fake_minimap2, tmp_path):
    write_fastq(tmp_path / "x" / "MINE" / "r.fastq", [("v_1", "ACGT")])  # same name as data/raw_data/MINE
    write_fastq(tmp_path / "y" / "MINE" / "r.fastq", [("v_1", "ACGT")])
    assert add(data_dir, "ref.fasta", "F", "f", str(tmp_path / "x" / "MINE"), str(tmp_path / "y" / "MINE")) == 0
    assert organisms(data_dir)["F"].source_dataset == "MINE-2;MINE-3"


def test_add_from_ref_rejects_overlapping_inputs(data_dir, raw, fake_minimap2, capsys):
    assert add(data_dir, "ref.fasta", "G", "g", "PUB", str(raw / "PUB" / "run1.fastq")) == 1
    assert "included twice, via PUB and PUB/run1.fastq" in capsys.readouterr().err


def test_add_from_ref_unmounted_own_dataset(data_dir, raw, fake_minimap2, tmp_path, capsys):
    (data_dir / "own_datasets.tsv").write_text(OWN_HEADER + f"EXT\t\t\t\t{tmp_path / 'unmounted' / 'run'}\n")
    assert add(data_dir, "ref.fasta", "H", "h", "EXT") == 1
    assert "which does not exist (is the drive mounted?)" in capsys.readouterr().err


# ------------------------------------------------ create_sispa_run by depth ---

@pytest.fixture
def tiled(data_dir, raw, fake_minimap2, tmp_path):
    """Organism TILE: 10 reads each on 0-40, 30-70 and 60-100 of a 100 bp reference."""
    reads = [(f"at{pos}_{i}", "ACGT" * 10) for pos in (1, 31, 61) for i in range(10)]
    write_fastq(raw / "TILE-DS" / "reads.fastq", reads)
    assert add(data_dir, "ref.fasta", "TILE", "tile", "TILE-DS") == 0
    return data_dir


def groups(names):
    return {pos: sum(n.startswith(f"at{pos}_") for n in names) for pos in (1, 31, 61)}


def test_mix_by_minimum_depth_from_a_table(tiled, tmp_path, capsys):
    table = tmp_path / "depth_mix.csv"
    table.write_text("organism_id,n_reads,minimum_depth\nTILE,,3\nhuman,7,\n")
    assert run_mix(tiled, tmp_path, str(table), "--seed", "4") == 0
    names = names_in(tmp_path / "out" / "fastqs" / "depth_mix.fastq.gz")
    # each third of the reference is only covered by its own reads: 3 of each is just enough
    assert groups(names) == {1: 3, 31: 3, 61: 3}
    assert sum(n.startswith("human_") for n in names) == 7
    assert "of the reference at 3x or more" in capsys.readouterr().out
    # what the run contains is recorded; TILE: 3x on each third, 6x where two thirds overlap
    composition = tmp_path / "out" / "compositions" / "depth_mix.csv"
    assert composition.read_text() == (
        "organism_id,n_reads,min_depth,max_depth,mean_depth\n"
        "TILE,9,3,6,3.60\n"
        "human,7,0,7,3.50\n"
    )
    # the same table and seed recreate the run
    assert run_mix(tiled, tmp_path, str(table), "--seed", "4", "--name", "again") == 0
    assert names_in(tmp_path / "out" / "fastqs" / "again.fastq.gz") == names


def test_mix_by_minimum_depth_from_a_list(tiled, tmp_path, capsys):
    assert run_mix(tiled, tmp_path, "TILE 15x", "--seed", "1") == 0
    names = names_in(tmp_path / "out" / "fastqs" / "TILE-15x.fastq.gz")
    assert groups(names) == {1: 10, 31: 10, 61: 10}  # 15x is more than any third has: every read
    out = capsys.readouterr().out
    assert "20.00% of the reference at 15x or more, 80.00% below with every read there" in out
    assert (tmp_path / "out" / "compositions" / "TILE-15x.csv").read_text() == (
        "organism_id,n_reads,min_depth,max_depth,mean_depth\nTILE,30,10,20,12.00\n")


def test_mix_all_reads_in_the_depth_column(tiled, tmp_path):
    table = tmp_path / "all_mix.csv"
    table.write_text("organism_id,minimum_depth\nTILE,all\nLASV,5\n")
    assert run_mix(tiled, tmp_path, str(table), "--seed", "2") == 0
    names = names_in(tmp_path / "out" / "fastqs" / "all_mix.fastq.gz")
    assert groups(names) == {1: 10, 31: 10, 61: 10}
    assert sum(n.startswith("lasv_") for n in names) == 5  # its 20 reads all cover 0-12: 5 give 5x
    assert (tmp_path / "out" / "compositions" / "all_mix.csv").read_text() == (
        "organism_id,n_reads,min_depth,max_depth,mean_depth\nTILE,30,10,20,12.00\nLASV,5,0,5,0.30\n")


def test_mix_needs_the_alignments(tiled, tmp_path, capsys):
    alignments_path(tiled, "virus_reads/COVID.fastq.gz").unlink()
    for composition in ("COVID 3x", "COVID 5"):
        assert run_mix(tiled, tmp_path, composition) == 1
        assert "has no alignment file" in capsys.readouterr().err
    for bad in ("TILE 0x", "TILE 3y"):
        assert run_mix(tiled, tmp_path, bad) == 1
        assert "minimum depth such as 20x" in capsys.readouterr().err
    table = tmp_path / "bad.csv"
    table.write_text("organism_id,minimum_depth\nTILE,deep\n")
    assert run_mix(tiled, tmp_path, str(table)) == 1
    assert "minimum_depth for TILE must be a positive number" in capsys.readouterr().err


def test_alignments_go_with_their_fastq(tiled, raw, fake_minimap2):
    tile = tiled / "virus_reads" / "tile.fastq.gz"
    assert alignments_path(tiled, "virus_reads/tile.fastq.gz").exists()
    tile.unlink()
    assert add(tiled, "ref.fasta", "Z", "z", "PUB") == 0
    assert not alignments_path(tiled, "virus_reads/tile.fastq.gz").exists() and "TILE" not in organisms(tiled)


# ------------------------------------------------------- remove_viral_reads ---

def write_bam(path, reference, names, unmapped=()):
    """A minimal BAM file (gzip is enough for the reader) holding the named reads."""
    import struct
    text = f"@SQ\tSN:{reference}\tLN:100\n".encode()
    data = b"BAM\1" + struct.pack("<i", len(text)) + text + struct.pack("<i", 1)
    data += struct.pack("<i", len(reference) + 1) + reference.encode() + b"\0" + struct.pack("<i", 100)
    for name, flag in [(n, 0) for n in names] + [(n, 4) for n in unmapped]:
        read = name.encode() + b"\0"
        record = struct.pack("<iiBBHHHiiii", 0, 0, len(read), 60, 0, 0, flag, 0, -1, -1, 0) + read
        data += struct.pack("<i", len(record)) + record
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as fh:
        fh.write(data)


CONSENSUS_TABLE = (
    "\tReference\tLength\tMapped reads\tAmbiguous positions\tConsensusLength\tAverage read coverage\t"
    "Description\tFamily\tOrganism\tSegment\tOrientation\tCurated\tOrganism Label\tPositions called\tCoverage\tIsBest\n"
    "0\tV1.1\t100\t3\t20\t100\t5.0\tV1.1 some virus\t\tsome virus\tUnknown\tUnknown\tFalse\tNon-Curated\t80\t80.0\tFalse\n"
    "1\tV2.1\t100\t2\t90\t100\t1.0\tV2.1 a phage\t\ta phage\tUnknown\tUnknown\tFalse\tNon-Curated\t10\t10.0\tFalse\n"
)


def vimop_result(sample_dir):
    """What vimop writes for a sample: V1 at 80% recovery with 3 reads, V2 at 10% with 2."""
    (sample_dir / "tables").mkdir(parents=True)
    (sample_dir / "tables" / "consensus.tsv").write_text(CONSENSUS_TABLE)
    write_bam(sample_dir / "consensus" / "V1.reads.bam", "V1.1", ["v_1", "half_1", "chim_1"], unmapped=["h_1"])
    write_bam(sample_dir / "consensus" / "V2.reads.bam", "V2.1", ["low_1", "v_2", "v_1"])
    return sample_dir


def test_bam_read_names(tmp_path):
    write_bam(tmp_path / "x.bam", "V1.1", ["a", "bb"], unmapped=["c"])
    assert remove_viral_reads.bam_read_names(tmp_path / "x.bam") == {"a", "bb"}


# Stand-in for `nextflow run ... --fastq DIR --out_dir OUT`: records how it was
# called and puts a prepared vimop result into OUT/<name of DIR>.
FAKE_NEXTFLOW = """\
    #!{python}
    import json, os, shutil, sys
    args = sys.argv[1:]
    fastq, out = args[args.index("--fastq") + 1], args[args.index("--out_dir") + 1]
    with open(os.environ["FAKE_NEXTFLOW_LOG"], "w") as fh:
        json.dump({"args": args, "cwd": os.getcwd(), "staged": sorted(os.listdir(fastq))}, fh)
    os.makedirs("work/ab/cdef")  # nextflow's work folder and log, in the launch folder
    open(".nextflow.log", "w").close()
    if os.environ.get("FAKE_NEXTFLOW_FAIL"):
        sys.exit(1)
    shutil.copytree(os.environ["FAKE_VIMOP_RESULT"], os.path.join(out, os.path.basename(fastq)))
"""


@pytest.fixture
def fake_nextflow(tmp_path, monkeypatch):
    exe = tmp_path / "bin" / "nextflow"
    exe.parent.mkdir(exist_ok=True)
    exe.write_text(textwrap.dedent(FAKE_NEXTFLOW).replace("{python}", sys.executable))
    exe.chmod(0o755)
    monkeypatch.setenv("FAKE_NEXTFLOW_LOG", str(tmp_path / "nextflow_call.json"))
    monkeypatch.setenv("FAKE_VIMOP_RESULT", str(vimop_result(tmp_path / "vimop_result")))
    return exe


def remove_viral(tmp_path, *args):
    return remove_viral_reads.main([*args, "--data-dir", str(tmp_path / "data")])


def test_remove_viral_reads_with_vimop(tmp_path, fake_nextflow, capsys):
    import json
    run = write_fastq(tmp_path / "run" / "run.fastq", MIXED_READS)
    config = tmp_path / "vimop.config"
    config.write_text("")
    out = tmp_path / "clean" / "clean.fastq.gz"
    vimop_dir = tmp_path / "vimop"
    assert remove_viral(tmp_path, str(run), "-o", str(out), "-c", str(config), "--vimop-dir", str(vimop_dir),
                        "--nextflow", str(fake_nextflow), "--viral-out", str(tmp_path / "viral.fastq"),
                        "--vimop-args", "--targets LASV") == 0

    # vimop ran on a folder named after the output, holding the run, with its output in the vimop folder
    call = json.loads((tmp_path / "nextflow_call.json").read_text())
    assert call["args"] == ["run", "opr-group-bnitm/vimop", "--fastq", str((vimop_dir / "input" / "clean").resolve()),
                            "--out_dir", str((vimop_dir / "output").resolve()), "-c", str(config.resolve()),
                            "-resume", "--targets", "LASV"]
    assert os.path.realpath(call["cwd"]) == os.path.realpath(vimop_dir) and call["staged"] == ["run.fastq"]
    assert os.path.samefile(vimop_dir / "input" / "clean" / "run.fastq", run)  # a hard link, not a copy

    # the reads vimop mapped to V1 (80% recovery) are removed; V2 (10%) keeps its reads
    assert names_in(out) == ["h_1", "low_1", "v_2/1", "h_2"]
    assert names_in(tmp_path / "viral.fastq") == ["v_1", "half_1", "chim_1"]
    assert (tmp_path / "clean" / "clean.removed_reads.tsv").read_text() == (
        "virus\tdescription\trecovery\tremoved_reads\n"
        "clean/V1.1\tV1.1 some virus\t80.00\t3\n"
        "clean/V2.1\tV2.1 a phage\t10.00\t0\n"
    )
    out_text = capsys.readouterr().out
    assert "4 of 7 reads kept, 3 viral reads removed" in out_text and "kept (below 50%)" in out_text


def test_remove_viral_reads_keeps_or_discards_the_vimop_run(tmp_path, fake_nextflow, monkeypatch, capsys):
    run = write_fastq(tmp_path / "run.fastq", MIXED_READS)
    common = ["--nextflow", str(fake_nextflow)]
    kept, discarded = tmp_path / "vimop_kept", tmp_path / "vimop_discarded"
    assert remove_viral(tmp_path, str(run), "-o", str(tmp_path / "a.fastq"), "--vimop-dir", str(kept), *common) == 0
    assert {p.name for p in kept.iterdir()} == {"input", "output", "work", ".nextflow.log"}
    assert remove_viral(tmp_path, str(run), "-o", str(tmp_path / "b.fastq"), "--vimop-dir", str(discarded),
                        "--discard-vimop", *common) == 0
    assert not discarded.exists() and names_in(tmp_path / "b.fastq") == names_in(tmp_path / "a.fastq")
    # an earlier vimop output is never deleted, and a failed run is kept to read its log
    assert remove_viral(tmp_path, str(run), "-o", str(tmp_path / "c.fastq"), "--vimop-output", str(kept / "output"),
                        "--discard-vimop") == 0
    assert (kept / "output").exists()
    monkeypatch.setenv("FAKE_NEXTFLOW_FAIL", "1")
    assert remove_viral(tmp_path, str(run), "-o", str(tmp_path / "d.fastq"), "--vimop-dir", str(discarded),
                        "--discard-vimop", *common) == 1
    assert (discarded / ".nextflow.log").exists()


def test_remove_viral_reads_min_recovery(tmp_path, fake_nextflow):
    run = write_fastq(tmp_path / "run.fastq", MIXED_READS)
    vimop_result(tmp_path / "earlier" / "run")
    common = ["--vimop-output", str(tmp_path / "earlier"), "--nextflow", "/does/not/exist"]  # vimop must not run
    assert remove_viral(tmp_path, str(run), "-o", str(tmp_path / "all.fastq"), "--min-recovery", "0", *common) == 0
    assert names_in(tmp_path / "all.fastq") == ["h_1", "h_2"]  # v_2/1 is v_2 in the BAM
    assert remove_viral(tmp_path, str(run), "-o", str(tmp_path / "none.fastq"), "--min-recovery", "90", *common) == 0
    assert names_in(tmp_path / "none.fastq") == [n for n, _ in MIXED_READS]


def test_remove_viral_reads_default_output(tmp_path, fake_nextflow, capsys):
    raw = tmp_path / "data" / "raw_data"
    write_fastq(raw / "RUN" / "run.fastq", MIXED_READS)
    assert remove_viral(tmp_path, "RUN", "--vimop-dir", str(tmp_path / "vimop"), "--nextflow", str(fake_nextflow)) == 0
    out = raw / "RUN_no_viral" / "RUN_no_viral.fastq.gz"
    assert names_in(out) == ["h_1", "low_1", "v_2/1", "h_2"]
    assert (raw / "RUN_no_viral" / "RUN_no_viral.removed_reads.tsv").exists()
    own = (tmp_path / "data" / "own_datasets.tsv").read_text().splitlines()
    assert own[1:] == ["RUN_no_viral\t\t\tRUN without the reads of viruses vimop found (min recovery 50%)\t"]
    # a run outside data/raw_data is named after its file
    write_fastq(tmp_path / "elsewhere" / "sample 7.fastq.gz", MIXED_READS)
    assert remove_viral(tmp_path, str(tmp_path / "elsewhere" / "sample 7.fastq.gz"),
                        "--vimop-output", str(tmp_path / "vimop" / "output")) == 0
    assert (raw / "sample_7_no_viral" / "sample_7_no_viral.fastq.gz").exists()


def test_remove_viral_reads_stages_several_files(tmp_path, fake_nextflow):
    import json
    first = write_fastq(tmp_path / "run" / "part1.fastq", MIXED_READS[:3])
    write_fastq(tmp_path / "run" / "barcode05" / "part2.fastq", MIXED_READS[3:])
    assert remove_viral(tmp_path, str(first), str(tmp_path / "run"), "-o", str(tmp_path / "out.fastq"),
                        "--vimop-dir", str(tmp_path / "vimop"), "--nextflow", str(fake_nextflow)) == 0
    staged = json.loads((tmp_path / "nextflow_call.json").read_text())["staged"]
    assert staged == ["001_part1.fastq", "002_part2.fastq"]  # part1 is given twice but used once
    assert names_in(tmp_path / "out.fastq") == ["h_1", "low_1", "v_2/1", "h_2"]


def test_remove_viral_reads_when_vimop_finds_nothing(tmp_path, capsys):
    run = write_fastq(tmp_path / "run.fastq", MIXED_READS)
    (tmp_path / "empty" / "run").mkdir(parents=True)
    assert remove_viral(tmp_path, str(run), "-o", str(tmp_path / "out.fastq"), "--vimop-output", str(tmp_path / "empty")) == 0
    assert names_in(tmp_path / "out.fastq") == [n for n, _ in MIXED_READS]
    assert "vimop reports 0 virus(es)" in capsys.readouterr().out


def test_remove_viral_reads_warns_about_foreign_vimop_output(tmp_path, capsys):
    run = write_fastq(tmp_path / "run.fastq", [("other_1", "ACGT")])
    vimop_result(tmp_path / "earlier" / "run")
    assert remove_viral(tmp_path, str(run), "-o", str(tmp_path / "out.fastq"), "--vimop-output", str(tmp_path / "earlier")) == 0
    assert "none of them is in the input" in capsys.readouterr().out


@pytest.mark.parametrize("setup, message", [
    ("fail", "vimop failed (exit code 1)"),
    ("same_dataset", "would be in dataset RUN next to its own input"),
    ("no_config", "nextflow config not found"),
    ("exists", "already exists"),
    ("no_bam", "no V1.reads.bam"),
])
def test_remove_viral_reads_errors(tmp_path, fake_nextflow, monkeypatch, capsys, setup, message):
    run = write_fastq(tmp_path / "data" / "raw_data" / "RUN" / "run.fastq", MIXED_READS)
    out = tmp_path / "data" / "raw_data" / "RUN-clean" / "clean.fastq"
    args = [str(run), "-o", str(out), "--vimop-dir", str(tmp_path / "vimop"), "--nextflow", str(fake_nextflow)]
    if setup == "fail":
        monkeypatch.setenv("FAKE_NEXTFLOW_FAIL", "1")
    elif setup == "same_dataset":
        args[2] = str(tmp_path / "data" / "raw_data" / "RUN" / "clean.fastq")
    elif setup == "no_config":
        args += ["-c", str(tmp_path / "missing.config")]
    elif setup == "exists":
        write_fastq(out, [("x", "ACGT")])
    elif setup == "no_bam":
        (tmp_path / "vimop_result" / "consensus" / "V1.reads.bam").unlink()
    assert remove_viral(tmp_path, *args) == 1
    assert message in capsys.readouterr().err
    if setup != "exists":
        assert not out.exists()
