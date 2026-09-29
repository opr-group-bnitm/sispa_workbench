# sispa_workbench
A tool that can be used to digitally create viral samples from existing or simulated data

It keeps a library of FASTQ files that each hold the reads of a single
organism (a virus, a bacterium, a host) and mixes them in chosen proportions
into artificial sequencing runs with a known composition.

```
publicly_available_datasets.tsv ──download_datasets.sh──▶ data/raw_data/
                                                               │
data/references/*.fasta ───────────────add_from_ref.py◀────────┤  (or data/mixed_reads/)
                                              │
your own single-organism FASTQs ─────▶ data/{virus,bacteria,host}_reads/
                                              │ sync_reads.py
                                              ▼
                                       data/organisms.tsv
                                              │ create_sispa_run.py
                                              ▼
                              output/fastqs/<name>.fastq.gz
                              output/compositions/<name>.csv
```

## Setup

```bash
conda env create -f environment.yml
conda activate sispa_workbench
```

The Python tools only need the standard library; `add_from_ref.py` needs
`minimap2`, and `download_datasets.sh` needs `curl`.

## Folders

| folder | contents |
|---|---|
| `data/raw_data/<dataset_id>/` | runs downloaded by `download_datasets.sh` |
| `data/references/` | reference FASTAs (or minimap2 `.mmi` indexes) for `add_from_ref.py` |
| `data/virus_reads/` | FASTQs with the reads of exactly one virus |
| `data/bacteria_reads/` | FASTQs with the reads of exactly one bacterium |
| `data/host_reads/` | FASTQs with the reads of exactly one host |
| `data/mixed_reads/` | your own FASTQs with reads of several organisms, to extract from with `add_from_ref.py` |
| `output/fastqs/` | artificial runs made by `create_sispa_run.py` |
| `output/compositions/` | what each artificial run is made of |

Sequencing data is git-ignored, and so is `data/organisms.tsv`, which describes
your local files. Composition CSVs are kept.

## 1. Download public datasets

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
keeping only SISPA datasets. Only `VEEV-TC83` (~1.6 GB, sequence known by
design) has `include=true`.

```bash
DRY_RUN=1 ./download_datasets.sh          # list runs and sizes, download nothing
./download_datasets.sh                    # every dataset with include=true
./download_datasets.sh BTV-REPS EVA71      # these, ignoring include
MAX_RUNS=2 PLATFORM=OXFORD_NANOPORE ./download_datasets.sh
```

Downloads are md5-checked, resumable and parallel (`JOBS=4`); rerunning skips
finished runs and retries failed ones. The full list of options is at the top
of the script.

## 2. Build the single-organism library

**Your own files:** copy a FASTQ (`.fastq`, `.fq`, optionally gzipped) into
`data/virus_reads/`, `data/bacteria_reads/` or `data/host_reads/`, then

```bash
./sync_reads.py
```

The file name without its suffix becomes the organism id
(`virus_reads/COVID.fastq.gz` → `COVID`). Ids may contain letters, digits,
`.`, `_` and `-`, and must be unique across the three folders.

**Extracted from a mixed run:** map a FASTQ against the organism's reference
and keep the reads that align:

```bash
# add_from_ref.py REFERENCE ORGANISM_ID FILENAME FASTQ [FASTQ ...]
./add_from_ref.py NC_045512.2.fasta COVID COVID.fastq.gz data/raw_data/SARS2-BC/*.fastq.gz
./add_from_ref.py --category host GRCh38.mmi human human.fastq.gz data/mixed_reads/sample1.fastq.gz
```

Reads with a primary minimap2 alignment are copied unchanged to
`data/<category>_reads/FILENAME` (`--category` is `virus` by default), and the
file is registered under ORGANISM_ID. A bare reference name is looked up in
`data/references/`.

minimap2 runs with the settings of the
[vimop](https://github.com/opr-group-bnitm/vimop) pipeline:
`-x map-ont --secondary=no`, plus `-k 11 -w 5` for `--category virus` (vimop's
`map_to_ref`); references over 500 MB are indexed in 2G parts (`-I 2G
--split-prefix`, vimop's `filter_virus_target`). As in vimop, every read with a
primary alignment is kept regardless of mapping quality. `--minimap2-args`
overrides these (e.g. `"-k 15 -w 10"`), `--preset` changes the preset, and
`--min-mapq` and `--min-aligned-fraction` (for example 0.8) make the extraction
stricter, for example to drop chimeric reads that are only partly from the
organism. For a large host genome, build the index once with
`minimap2 -x map-ont -d GRCh38.mmi GRCh38.fa` and pass the `.mmi`.

Either way the library is listed in `data/organisms.tsv`:

```
organism_id  filename                    avg_read_length  max_read_length  min_read_length  n_reads  reference                   source_fastq
COVID        virus_reads/COVID.fastq.gz  512.33           2911             87               48211    references/NC_045512.2.fasta raw_data/SARS2-BC/SRR15356294_1.fastq.gz
human        host_reads/human.fastq.gz   1203.10          40122            52               1500000
```

(example values). `reference` and `source_fastq` record what `add_from_ref.py`
mapped: paths relative to `data/` (absolute if outside it), several source
FASTQs separated by `;`. They stay empty for files you add yourself; you can
fill them in by hand and later syncs keep them. `sync_reads.py` adds new files, drops deleted ones and only
rereads files that changed; `add_from_ref.py` and `create_sispa_run.py` run it
for you.

## 3. Create an artificial run

```bash
./create_sispa_run.py COVID 5000, human 7000, DENV2 80000, LASV all
# -> output/fastqs/COVID-5000_human-7000_DENV2-80000_LASV-all.fastq.gz
#    output/compositions/COVID-5000_human-7000_DENV2-80000_LASV-all.csv

./create_sispa_run.py my_mix.csv                    # -> output/fastqs/my_mix.fastq.gz
./create_sispa_run.py COVID 500, human all --name low_titre --seed 42
```

For each organism the requested number of reads is drawn at random without
replacement from its FASTQ, or all of them for `all`. Reads from the different
organisms are interleaved in random order, as in a real run.

The composition file has the columns `organism_id,n_reads`, comma separated
(`.tsv` for tab separated). When the composition is given as a list,
`output/compositions/<name>.csv` records it with `all` resolved to the actual
count. When it comes from a CSV, that CSV already is the record and nothing is
written. The seed is printed; the same composition, seed and input files
reproduce a run exactly.

Options: `--name`, `--seed`, `--fastq-dir`, `--composition-dir`, `--force` to
overwrite an existing run. All tools take `--data-dir`; the environment
variables `SISPA_DATA_DIR` and `SISPA_OUTPUT_DIR` change the defaults.

## Tests

```bash
python -m pytest tests
```
