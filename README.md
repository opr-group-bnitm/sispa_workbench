# sispa_workbench
A tool that can be used to digitally create viral samples from existing or simulated data

It keeps a library of FASTQ files that each hold the reads of a single
organism (a virus, a bacterium, a host) and mixes them in chosen proportions
into artificial sequencing runs with a known composition.

```
publicly_available_datasets.tsv ─ download_datasets.sh ─▶ data/raw_data/<dataset_id>/ ─┐
data/own_datasets.tsv: your FASTQs, in data/raw_data/ or anywhere on the machine ──────┤
data/references/*.fasta ───────────────────────────────────────────────────────────────┤
                                                                                       ▼
                       data/{virus,bacteria,host}_reads/ + data/organisms.tsv ◀─ add_from_ref.py
                                    │
                                    ▼ create_sispa_run.py
                 output/fastqs/<name>.fastq.gz + output/compositions/<name>.csv
```

The single-organism FASTQs are made by `add_from_ref.py` alone, from datasets
listed in the two dataset tables, so every organism can be traced back to a
dataset and a reference.

## Setup

```bash
conda env create -f environment.yml
conda activate sispa_workbench
```

The Python tools only need the standard library; `add_from_ref.py` needs
`minimap2`, and `download_datasets.sh` needs `curl`. `remove_viral_reads.py`
needs `nextflow` and a working [vimop](https://github.com/opr-group-bnitm/vimop)
setup (its database and containers).

## Folders

| folder | contents |
|---|---|
| `data/raw_data/<dataset_id>/` | one folder per dataset: public runs downloaded by `download_datasets.sh`, and your own FASTQs |
| `data/references/` | reference FASTAs (or minimap2 `.mmi` indexes) for `add_from_ref.py` |
| `data/virus_reads/` | reads of exactly one virus per FASTQ, made by `add_from_ref.py` |
| `data/bacteria_reads/` | reads of exactly one bacterium per FASTQ, made by `add_from_ref.py` |
| `data/host_reads/` | reads of exactly one host per FASTQ, made by `add_from_ref.py` |
| `data/output/background_fastqs/` | runs without their viral reads, made by `remove_viral_reads.py` |
| `data/input_compositions/` | your compositions for `create_sispa_run.py`; start from `template.csv` |
| `data/alignments/` | where each read of those FASTQs aligns on its reference, e.g. `virus_reads/COVID.fastq.gz.alignments.tsv.gz`, made by `add_from_ref.py`; `create_sispa_run.py` needs them for depths |
| `output/fastqs/` | artificial runs made by `create_sispa_run.py` |
| `output/compositions/` | what each artificial run is made of |
| `output/vimop/<name>/` | vimop runs of `remove_viral_reads.py`: input, output, nextflow work folder and log |

Sequencing data is git-ignored, and so are the two tables that describe your
local files, `data/own_datasets.tsv` and `data/organisms.tsv`. Composition CSVs
are kept.

## 1. Datasets

### Public datasets

`publicly_available_datasets.tsv` lists the candidate datasets:

| column | meaning |
|---|---|
| `dataset_id` | name, also the folder under `data/raw_data/` |
| `include` | `true` / `false`: download it by default |
| `study_accession` | BioProject / study |
| `repository` | `ENA`, `SRA` or `DDBJ` are downloaded (all via the ENA mirror); others such as `Zenodo` are skipped |
| `dataset_accession` | what to download: run(s), sample(s), experiment(s) or a study, separated by `;`. Empty means every run of `study_accession` |
| `sample_type`, `organism_target` | what was sequenced |
| `notes` | study type, platform, ground truth, caveats |

It was seeded from `../data/input/sispa_benchmark/sispa_benchmark_runs.tsv`,
keeping only SISPA datasets, and extended since; set `include` to choose what
`download_datasets.sh` fetches by default.

```bash
DRY_RUN=1 ./download_datasets.sh          # list runs and sizes, download nothing
./download_datasets.sh                    # every dataset with include=true
./download_datasets.sh BTV-REPS EVA71      # these, ignoring include
MAX_RUNS=2 PLATFORM=OXFORD_NANOPORE ./download_datasets.sh
```

Downloads are md5-checked, resumable and parallel (`JOBS=4`); rerunning skips
finished runs and retries failed ones. The full list of options is at the top
of the script.

### Your own datasets

Your own FASTQs (`.fastq`, `.fq`, optionally gzipped, also in subfolders such
as `fastq_pass/barcode05/`) can stay where they are, or go into
`data/raw_data/<dataset_id>/`. They are listed in `data/own_datasets.tsv`:

| column | meaning |
|---|---|
| `dataset_id` | letters, digits, `.`, `_`, `-`; unique across both tables |
| `sample_type`, `organism_target` | what was sequenced |
| `notes` | anything worth knowing about it |
| `path` | where the FASTQs are, a folder or a single file; empty means `data/raw_data/<dataset_id>/` |

This table is yours: it is not in git, and it is created when first needed. You
rarely need to add rows yourself: data you give `add_from_ref.py` that is in
neither table (a folder in `data/raw_data/`, or a FASTQ file or folder anywhere
on the machine) gets a row, with an id made from its name and its path. Fill in
its description afterwards. Symlinks are followed, and a link to a drive that
is not mounted is reported instead of being skipped.

### Removing viral reads from a run

`remove_viral_reads.py` runs vimop on a run and removes the reads of the
viruses it found, for example to get a virus-free background:

```bash
./remove_viral_reads.py BOV-6760
# -> data/output/background_fastqs/BOV-6760_no_viral.fastq.gz
```

It runs vimop with its default settings,
`nextflow run opr-group-bnitm/vimop --fastq ... --out_dir ... -resume`, keeping vimop's input, output, work folder and log in
`output/vimop/<name>/`. Then, for every virus in vimop's
`tables/consensus.tsv` whose consensus reached `--min-recovery` (default 50%,
the table's `Coverage`: positions called, not N), it removes every read in that
virus's `consensus/<reference>.reads.bam`; all other reads are written
unchanged. `--min-recovery 0` removes the reads of everything vimop reported,
which can include real background: a low-level hit like Epstein-Barr virus in
human cells, or phages matching bacteria. Viruses with too few reads for a
consensus stay in with the default.

The input can be dataset ids, FASTQ files or folders of them. Without `-o`, the
cleaned run goes to `data/output/background_fastqs/<name>_no_viral.fastq.gz`
(written with `-o` into a new folder in `data/raw_data/`, it becomes a dataset
in `data/own_datasets.tsv`). The reads removed per virus are listed in
`<name>_no_viral.removed_reads.tsv` next to it, and `--viral-out` also writes
the removed reads. `--vimop-output` reuses an earlier vimop output of the same
reads instead of running vimop (every sample in it counts); `--vimop-args` and
`--nextflow-args` pass options on, e.g. `--vimop-args "--targets LASV"` or
`--nextflow-args "-profile docker"`, and `-c` adds a nextflow config if your
setup needs one.

The vimop run is kept by default, for vimop's report and to reuse it.
`--discard-vimop` deletes it (input, output, nextflow's work folder and logs)
once the cleaned run is written; a failed run is always kept, and an output
given with `--vimop-output` is never deleted.

## 2. Build the single-organism library

Organisms are only ever added by extracting their reads from datasets:
`add_from_ref.py` maps every FASTQ of the given datasets against the organism's
reference and keeps the reads that align. A DATASET is a dataset id, or a path
to a FASTQ file or a folder of FASTQs, inside or outside the workbench; a path
inside a known dataset (say, one run of a public dataset) uses just that part.

```bash
# add_from_ref.py REFERENCE ORGANISM_ID FILENAME DATASET [DATASET ...]
./add_from_ref.py NC_045512.2.fasta COVID COVID.fastq.gz SARS2-BC
./add_from_ref.py NC_045512.2.fasta COVID-1 COVID-1.fastq.gz data/raw_data/SARS2-BC/SRR15356294_1.fastq.gz
./add_from_ref.py --category host ARS-UCD2.0.fna cattle cattle.fastq.gz BOV-6760 BOV-6763
./add_from_ref.py lasv.fasta LASV LASV.fastq.gz /Volumes/runs/2025-06-01/fastq_pass/barcode05
```

Reads with a primary minimap2 alignment are copied unchanged to
`data/<category>_reads/FILENAME` (`--category` is `virus` by default), and the
file is registered in `data/organisms.tsv` under ORGANISM_ID. Ids may contain
letters, digits, `.`, `_` and `-`, and must be unique. A bare reference name is
looked up in `data/references/`.

minimap2 runs with the settings of the
[vimop](https://github.com/opr-group-bnitm/vimop) pipeline:
`-ax map-ont --secondary=no`, plus `-k 11 -w 5` for `--category virus` (vimop's
`map_to_ref`); references over 500 MB are indexed in 2G parts (`-I 2G
--split-prefix`, vimop's `filter_virus_target`). As in vimop, every read with a
primary alignment is kept regardless of mapping quality. `--minimap2-args`
overrides these (e.g. `"-k 15 -w 10"`), `--preset` changes the preset, and
`--min-mapq` and `--min-aligned-fraction` (for example 0.8) make the extraction
stricter, for example to drop chimeric reads that are only partly from the
organism. For a large host genome, build the index once with
`minimap2 -x map-ont -d GRCh38.mmi GRCh38.fa` and pass the `.mmi`.

The library is listed in `data/organisms.tsv` (example values):

```
organism_id  filename                    avg_read_length  max_read_length  min_read_length  n_reads  breadth_coverage  min_depth  max_depth  reference                     source_dataset
COVID        virus_reads/COVID.fastq.gz  512.33           2911             87               48211    99.87             0          1532       references/NC_045512.2.fasta  SARS2-BC
cattle       host_reads/cattle.fastq.gz  1203.10          40122            52               1500000  0.42              0          311        references/ARS-UCD2.0.fna     BOV-6760;BOV-6763
```

- `breadth_coverage`: % of reference positions covered by at least one read.
- `min_depth`, `max_depth`: the depth range over all reference positions,
  counted like vimop's `samtools depth -aa -J` (every alignment spanning a
  position, deletions included). `min_depth` is 0 whenever breadth is below
  100%. A reference with several sequences (segments) is taken as a whole; the
  per-sequence numbers are printed while `add_from_ref.py` runs. With a
  reference holding several strains, breadth is diluted accordingly.
- `reference`: relative to `data/` (absolute if outside it).
- `source_dataset`: the dataset(s) the reads came from, separated by `;`, as
  `<dataset_id>` or, when only part of a dataset was used,
  `<dataset_id>/<file or folder>` (e.g. `SARS2-BC/SRR15356294_1.fastq.gz`).

FASTQs put into the read folders any other way are not registered and cannot be
used; `add_from_ref.py` lists them as ignored. To remove an organism, delete its
FASTQ: the next `add_from_ref.py` run drops its row.

## 3. Create an artificial run

```bash
./create_sispa_run.py COVID 5000, human 7000, DENV2 80000, LASV all
# -> output/fastqs/COVID-5000_human-7000_DENV2-80000_LASV-all.fastq.gz
#    output/compositions/COVID-5000_human-7000_DENV2-80000_LASV-all.csv

./create_sispa_run.py data/input_compositions/my_mix.csv   # -> output/fastqs/my_mix.fastq.gz
./create_sispa_run.py COVID 500, human all --name low_titre --seed 42
./create_sispa_run.py LASV 20x, human 50000         # LASV at a minimum depth of 20
```

The organism ids are those in `data/organisms.tsv`. Per organism, the
composition asks for either

- a number of reads (`5000`): drawn at random without replacement from its
  FASTQ, or all of them (`all`); or
- a minimum depth (`20x`): reads are drawn at random, and a read is only taken
  if it covers a position of the reference that is still below 20x. This stops
  once every position is at 20x or more, so only as many reads as needed are
  used. Where all of the organism's reads together are not 20x deep, every
  read there is taken. Depth is counted as in `organisms.tsv`, and the run
  reports how much of the reference reached the depth.

Reads from the different organisms are interleaved in random order, as in a
real run.

The composition file is comma separated (`.tsv` for tab separated), with the
column `organism_id` and the column `n_reads` and/or `minimum_depth`; each row
fills one of them, and lines starting with `#` are ignored. To write one, copy
`data/input_compositions/template.csv` under the name the run should get; it
explains the columns and has an example of each kind of row:

```
organism_id,n_reads,minimum_depth
LASV,,20
human,50000,
```

A row with a `minimum_depth` asks for that depth; an `n_reads` next to it is
ignored. `all` works in either column and takes every read of the organism, so
a table can also use `minimum_depth` alone:

```
organism_id,minimum_depth
LASV,20
human,all
```

Every run records what it contains in `output/compositions/<name>.csv`, one
row per organism:

```
organism_id,n_reads,min_depth,max_depth,mean_depth
LASV,563,20,92,31.42
human,50000,0,12,0.02
```

- `n_reads`: the reads it got.
- `min_depth`, `max_depth`, `mean_depth`: the depth these reads reach on its
  reference, counted like in `organisms.tsv` (over every reference position,
  so `min_depth` is 0 wherever the reads leave a gap).

The seed is printed; the same composition, seed and organisms reproduce a run
exactly.

Options: `--name`, `--seed`, `--fastq-dir`, `--composition-dir`, `--force` to
overwrite an existing run. All tools take `--data-dir`; the environment
variables `SISPA_DATA_DIR` and `SISPA_OUTPUT_DIR` change the defaults.

## Tests

```bash
python -m pytest tests
```
