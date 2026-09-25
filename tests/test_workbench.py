import csv
import gzip
import random
import shutil
import sys
import textwrap
from pathlib import Path

import pytest

import add_from_ref
import create_sispa_run
import sync_reads
from workbench.fastq import FastqError, fastq_stats, read_fastq
from workbench.organisms import SyncError, load_organisms


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


@pytest.fixture
def data_dir(tmp_path):
    d = tmp_path / "data"
    write_fastq(d / "virus_reads" / "COVID.fastq.gz", [(f"covid_{i}", "ACGT" * (i + 1)) for i in range(50)])
    write_fastq(d / "virus_reads" / "LASV.fq", [(f"lasv_{i}", "GGCC" * 3) for i in range(20)])
    write_fastq(d / "host_reads" / "human.fastq", [(f"human_{i}", "T" * 100) for i in range(300)])
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


# ------------------------------------------------------------- sync_reads ---

def test_sync_registers_files(data_dir):
    assert sync_reads.main(["--data-dir", str(data_dir)]) == 0
    orgs = organisms(data_dir)
    assert set(orgs) == {"COVID", "LASV", "human"}
    covid = orgs["COVID"]
    assert covid.filename == "virus_reads/COVID.fastq.gz"
    assert (covid.n_reads, covid.min_read_length, covid.max_read_length) == (50, 4, 200)
    assert covid.avg_read_length == pytest.approx(102.0)
    header = (data_dir / "organisms.tsv").read_text().splitlines()[0].split("\t")
    assert header == ["organism_id", "filename", "avg_read_length", "max_read_length", "min_read_length", "n_reads"]


def test_sync_only_reads_changed_files(data_dir, capsys):
    sync_reads.main(["--data-dir", str(data_dir)])
    capsys.readouterr()
    sync_reads.main(["--data-dir", str(data_dir)])
    assert "[stats]" not in capsys.readouterr().out

    write_fastq(data_dir / "virus_reads" / "LASV.fq", [("lasv_0", "GGCC")])
    sync_reads.main(["--data-dir", str(data_dir)])
    out = capsys.readouterr().out
    assert "reading virus_reads/LASV.fq" in out and "COVID" not in out
    assert organisms(data_dir)["LASV"].n_reads == 1


def test_sync_drops_removed_files_and_keeps_explicit_ids(data_dir):
    from workbench.organisms import sync_reads as sync
    sync(data_dir, explicit_ids={"virus_reads/LASV.fq": "Lassa"})
    (data_dir / "virus_reads" / "COVID.fastq.gz").unlink()
    sync(data_dir)
    assert set(organisms(data_dir)) == {"Lassa", "human"}


def test_sync_rejects_duplicate_ids(data_dir, capsys):
    write_fastq(data_dir / "bacteria_reads" / "COVID.fq", [("x", "ACGT")])
    assert sync_reads.main(["--data-dir", str(data_dir)]) == 1
    assert "COVID" in capsys.readouterr().err
    with pytest.raises(SyncError):
        from workbench.organisms import sync_reads as sync
        sync(data_dir)


def test_sync_skips_unusable_names(data_dir, capsys):
    write_fastq(data_dir / "virus_reads" / "my virus.fastq", [("x", "ACGT")])
    write_fastq(data_dir / "virus_reads" / ".partial.hidden.fastq", [("x", "ACGT")])
    (data_dir / "virus_reads" / "notes.txt").write_text("not a fastq")
    assert sync_reads.main(["--data-dir", str(data_dir)]) == 0
    assert "[skipped] virus_reads/my virus.fastq" in capsys.readouterr().out
    assert set(organisms(data_dir)) == {"COVID", "LASV", "human"}


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
    assert rows == [
        {"organism_id": "COVID", "n_reads": "10"},
        {"organism_id": "human", "n_reads": "200"},
        {"organism_id": "LASV", "n_reads": "20"},
    ]


def test_mix_single_quoted_list_and_name(data_dir, tmp_path):
    assert run_mix(data_dir, tmp_path, "COVID 5, LASV 5", "--name", "tiny", "--seed", "3") == 0
    assert len(names_in(tmp_path / "out" / "fastqs" / "tiny.fastq.gz")) == 10
    assert (tmp_path / "out" / "compositions" / "tiny.csv").exists()


def test_mix_from_csv_uses_its_name_and_writes_no_composition(data_dir, tmp_path):
    table = tmp_path / "spike_in.csv"
    table.write_text("organism_id,n_reads\nCOVID,all\nhuman,7\n")
    assert run_mix(data_dir, tmp_path, str(table)) == 0
    assert len(names_in(tmp_path / "out" / "fastqs" / "spike_in.fastq.gz")) == 57
    assert not (tmp_path / "out" / "compositions").exists()


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
    write_fastq(data_dir / "virus_reads" / "V.fastq", [(f"r{i}", "ACGT") for i in range(10)])
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


def test_mix_refuses_to_overwrite(data_dir, tmp_path, capsys):
    assert run_mix(data_dir, tmp_path, "LASV 2", "--seed", "1") == 0
    assert run_mix(data_dir, tmp_path, "LASV 2", "--seed", "1") == 1
    assert "already exists" in capsys.readouterr().err
    assert run_mix(data_dir, tmp_path, "LASV 2", "--seed", "1", "--force") == 0


# ----------------------------------------------------------- add_from_ref ---

FAKE_MINIMAP2 = """\
    #!{python}
    # Stand-in for minimap2 that "maps" reads by name:
    #   v_*    align end to end,            mapq 60
    #   half_* align over half their length, mapq 60
    #   low_*  align end to end,            mapq 1
    #   h_*    do not align
    import sys
    fastq = sys.argv[-1]
    with open(fastq) as fh:
        lines = fh.read().splitlines()
    for header, seq in zip(lines[0::4], lines[1::4]):
        name, n = header[1:].split()[0], len(seq)
        name = name[:-2] if name.endswith(("/1", "/2")) else name
        span, mapq = {{"v": (n, 60), "half": (n // 2, 60), "low": (n, 1)}}.get(name.split("_")[0], (0, 0))
        if span:
            print(f"{{name}}\\t{{n}}\\t0\\t{{span}}\\t+\\tref\\t1000\\t0\\t{{span}}\\t{{span}}\\t{{span}}\\t{{mapq}}\\ttp:A:P")
"""


@pytest.fixture
def fake_minimap2(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "minimap2"
    exe.write_text(textwrap.dedent(FAKE_MINIMAP2.format(python=sys.executable)))
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path(shutil.which('sh')).parent}")


@pytest.fixture
def mixed_fastq(tmp_path):
    return write_fastq(tmp_path / "mixed.fastq", [
        ("v_1", "ACGT" * 10), ("h_1", "TTTT"), ("half_1", "ACGT" * 10),
        ("low_1", "ACGT" * 10), ("v_2/1", "ACGT"), ("h_2", "TTTT"),
    ])


def add(data_dir, *args):
    ref = data_dir / "references" / "ref.fasta"
    ref.parent.mkdir(parents=True, exist_ok=True)
    ref.write_text(">ref\nACGT\n")
    return add_from_ref.main(["--data-dir", str(data_dir), *args])


def test_add_from_ref_extracts_mapped_reads(data_dir, fake_minimap2, mixed_fastq):
    assert add(data_dir, "ref.fasta", "DENV2", "dengue", str(mixed_fastq)) == 0
    out = data_dir / "virus_reads" / "dengue.fastq.gz"
    assert names_in(out) == ["v_1", "half_1", "low_1", "v_2/1"]
    # records are copied unchanged, comments included
    assert next(read_fastq(out))[0] == b"v_1 some comment"
    orgs = organisms(data_dir)
    assert orgs["DENV2"].filename == "virus_reads/dengue.fastq.gz" and orgs["DENV2"].n_reads == 4

    # the explicit id survives later syncs
    sync_reads.main(["--data-dir", str(data_dir)])
    assert "DENV2" in organisms(data_dir) and "dengue" not in organisms(data_dir)


def test_add_from_ref_filters(data_dir, fake_minimap2, mixed_fastq):
    assert add(data_dir, "--category", "host", "--min-mapq", "20", "--min-aligned-fraction", "0.8",
               "ref.fasta", "strict", "strict.fq", str(mixed_fastq)) == 0
    assert names_in(data_dir / "host_reads" / "strict.fq") == ["v_1", "v_2/1"]


def test_add_from_ref_several_inputs(data_dir, fake_minimap2, mixed_fastq, tmp_path):
    other = write_fastq(tmp_path / "other.fastq", [("v_9", "ACGT"), ("h_9", "ACGT")])
    assert add(data_dir, "ref.fasta", "X", "x", str(mixed_fastq), str(other)) == 0
    assert names_in(data_dir / "virus_reads" / "x.fastq.gz")[-1] == "v_9"


@pytest.mark.parametrize("args, message", [
    (["ref.fasta", "COVID", "other.fastq.gz"], "already used by virus_reads/COVID.fastq.gz"),
    (["ref.fasta", "NEW", "COVID.fastq.gz"], "already exists"),
    (["ref.fasta", "bad id", "x.fastq"], "may only contain"),
    (["ref.fasta", "NEW", "sub/x.fastq"], "plain file name"),
    (["nope.fasta", "NEW", "x.fastq"], "reference not found"),
])
def test_add_from_ref_errors(data_dir, fake_minimap2, mixed_fastq, capsys, args, message):
    sync_reads.main(["--data-dir", str(data_dir)])
    assert add(data_dir, *args, str(mixed_fastq)) == 1
    assert message in capsys.readouterr().err


def test_add_from_ref_nothing_mapped(data_dir, fake_minimap2, tmp_path, capsys):
    host_only = write_fastq(tmp_path / "host.fastq", [("h_1", "ACGT")])
    assert add(data_dir, "ref.fasta", "NEW", "new.fastq.gz", str(host_only)) == 1
    assert "no reads mapped" in capsys.readouterr().err
    assert not list((data_dir / "virus_reads").glob("*new*"))


def test_minimap2_args_follow_vimop(tmp_path):
    small = tmp_path / "virus.fasta"
    small.write_text(">v\nACGT\n")
    split = tmp_path / "idx"
    assert add_from_ref.minimap2_args(small, "virus", split) == ["--secondary=no", "-k", "11", "-w", "5"]
    assert add_from_ref.minimap2_args(small, "host", split) == ["--secondary=no"]
    big = tmp_path / "host.fasta"
    with open(big, "wb") as fh:
        fh.truncate(add_from_ref.VIMOP_SPLIT_THRESHOLD + 1)
    assert add_from_ref.minimap2_args(big, "host", split) == ["--secondary=no", "-I", "2G", "--split-prefix", str(split)]


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not installed")
def test_add_from_ref_with_real_minimap2(tmp_path):
    rng = random.Random(0)
    genome = "".join(rng.choice("ACGT") for _ in range(20000))
    other = "".join(rng.choice("ACGT") for _ in range(20000))
    data_dir = tmp_path / "data"
    (data_dir / "references").mkdir(parents=True)
    (data_dir / "references" / "virus.fasta").write_text(f">virus\n{genome}\n")
    reads = [(f"virus_{i}", genome[i * 900:i * 900 + 2000]) for i in range(20)]
    reads += [(f"other_{i}", other[i * 900:i * 900 + 2000]) for i in range(20)]
    fastq = write_fastq(tmp_path / "mixed.fastq.gz", reads)
    assert add_from_ref.main(["--data-dir", str(data_dir), "virus.fasta", "V", "V.fastq.gz", str(fastq)]) == 0
    assert sorted(names_in(data_dir / "virus_reads" / "V.fastq.gz")) == sorted(n for n, _ in reads[:20])
