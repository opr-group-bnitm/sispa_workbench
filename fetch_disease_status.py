#!/usr/bin/env python3
"""Write the disease status of every run of a public BioProject as a TSV.

The run list (run -> BioSample) comes from the ENA Portal API, and each
BioSample's attributes from EBI BioSamples (which mirrors NCBI BioSample):

    ./fetch_disease_status.py PRJNA634678 -o PRJNA634678_disease_status.tsv

Columns: dataset_id (looked up by run accession in
publicly_available_datasets.tsv; empty if the run is not listed),
run_accession, sample_accession, sample_alias, biosample_title,
host_disease (the BioSample "host disease" attribute as submitted),
diseased (false if host_disease is "healthy", true for any other value,
empty if there is none), host_sex, collection_year.

For PRJNA634678 (Zhang et al. 2021, doi:10.1111/tbed.13873) host_disease is
BRD or Healthy, and refers to the outcome: swabs were taken on arrival at the
feedlot, and BRD animals developed BRD within the next 40 days.

Standard library only; about 300 requests, 8 at a time, retried on failure.
"""

import argparse
import csv
import json
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ENA_RUNS = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run"
            "&fields=run_accession,sample_accession,sample_alias&format=tsv&limit=0")
BIOSAMPLE = "https://www.ebi.ac.uk/biosamples/samples/{}"
TABLE = Path(__file__).resolve().parent / "publicly_available_datasets.tsv"
COLUMNS = ["dataset_id", "run_accession", "sample_accession", "sample_alias", "biosample_title",
           "host_disease", "diseased", "host_sex", "collection_year"]


def get(url, tries=4):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return response.read().decode()
        except Exception as error:
            if attempt == tries - 1:
                raise RuntimeError(f"{url}: {error}") from error
            time.sleep(2 ** attempt)


def runs(project):
    rows = list(csv.DictReader(get(ENA_RUNS.format(project)).splitlines(), delimiter="\t"))
    if not rows:
        sys.exit(f"ENA lists no runs for {project}")
    return sorted(rows, key=lambda row: row["run_accession"])


def attributes(sample):
    """The BioSample's attributes, e.g. {"host disease": "BRD", "host sex": "male", ...}."""
    characteristics = json.loads(get(BIOSAMPLE.format(sample)))["characteristics"]
    return {name: "|".join(value["text"] for value in values)
            for name, values in characteristics.items()}


def dataset_ids(table):
    """dataset_id of every accession listed in the table's dataset_accession column."""
    ids = {}
    if table.exists():
        with open(table, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                for accession in (row.get("dataset_accession") or "").split(";"):
                    if accession.strip():
                        ids[accession.strip()] = row["dataset_id"]
    return ids


def diseased(host_disease):
    if not host_disease:
        return ""
    return "false" if host_disease.strip().lower() == "healthy" else "true"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("project", help="BioProject or study accession, e.g. PRJNA634678")
    parser.add_argument("-o", "--out", help="output TSV (default: stdout)")
    parser.add_argument("--table", type=Path, default=TABLE,
                        help="where dataset_id is looked up (default: %(default)s)")
    args = parser.parse_args()

    run_rows = runs(args.project)
    samples = sorted({row["sample_accession"] for row in run_rows})
    with ThreadPoolExecutor(8) as pool:
        attrs = dict(zip(samples, pool.map(attributes, samples)))
    ids = dataset_ids(args.table)

    out = open(args.out, "w", newline="") if args.out else sys.stdout
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(COLUMNS)
    counts = {}
    for row in run_rows:
        sample = attrs[row["sample_accession"]]
        host_disease = sample.get("host disease", "")
        counts[host_disease or "(none)"] = counts.get(host_disease or "(none)", 0) + 1
        writer.writerow([ids.get(row["run_accession"], ""), row["run_accession"],
                        row["sample_accession"], row["sample_alias"], sample.get("title", ""),
                        host_disease, diseased(host_disease), sample.get("host sex", ""),
                        sample.get("collection date", "")[:4]])
    if args.out:
        out.close()
    summary = ", ".join(f"{value} {n}" for value, n in sorted(counts.items()))
    print(f"[done]    {len(run_rows)} runs of {args.project}: {summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
