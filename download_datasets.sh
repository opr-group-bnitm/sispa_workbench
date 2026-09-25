#!/usr/bin/env bash
#
# Download the FASTQ files of every dataset in publicly_available_datasets.tsv
# that has include=true, from the ENA mirror (ENA, SRA and DDBJ are all served
# there).
#
# Which runs belong to a dataset:
#   dataset_accession set   -> only that accession; it may be a run (SRR/ERR/DRR),
#                              experiment, sample (SAMN/SAMEA) or study, and
#                              several can be given separated by ';' or ','
#   dataset_accession empty -> every run of study_accession
#
# Layout produced:
#   data/raw_data/<dataset_id>/
#     ├── runs.tsv                   # run manifest fetched from ENA
#     ├── SRRxxxxxxx_1.fastq.gz      # files exactly as ENA names them
#     ├── ...
#     └── .done/SRRxxxxxxx           # written only after every md5 checks out
#
# Usage:
#   ./download_datasets.sh                         # every row with include=true
#   ./download_datasets.sh ID [ID ...]             # these rows, whatever include says
#   DRY_RUN=1 ./download_datasets.sh               # manifests + sizes only, no FASTQ
#   MAX_RUNS=3 ./download_datasets.sh              # first 3 runs of each dataset
#   PLATFORM=OXFORD_NANOPORE ./download_datasets.sh  # skip runs from other platforms
#   JOBS=8 ./download_datasets.sh                  # parallel downloads (default 4)
#   REFRESH=1 ./download_datasets.sh               # re-query ENA for the manifests
#   PROTO=ftp ./download_datasets.sh               # ENA over FTP instead of HTTPS
#   TSV=other.tsv RAW_DIR=/big/disk ./download_datasets.sh
#
# Safe to re-run and safe to interrupt: finished runs are skipped, interrupted
# transfers resume, and anything that fails its md5 is discarded. Runs that
# fail are listed in data/raw_data/failed_runs.txt; just run again.
#
# Runs that ENA holds no FASTQ for are fetched with fasterq-dump when sra-tools
# is installed. Rows from other repositories (Zenodo, ...) are skipped with a
# warning: put those files into data/raw_data/<dataset_id>/ by hand.
#
# Requires: curl, awk, gzip, md5sum (Linux) or md5 (macOS). Optional: fasterq-dump.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TSV="${TSV:-$REPO_ROOT/publicly_available_datasets.tsv}"
RAW_DIR="${RAW_DIR:-$REPO_ROOT/data/raw_data}"
JOBS="${JOBS:-4}"
DRY_RUN="${DRY_RUN:-0}"
MAX_RUNS="${MAX_RUNS:-0}"                 # 0 = no limit
PLATFORM="${PLATFORM:-ALL}"               # ENA instrument_platform value, or ALL
REFRESH="${REFRESH:-0}"
PROTO="${PROTO:-https}"                   # https (recommended), ftp
RETRIES="${RETRIES:-8}"
RETRY_DELAY="${RETRY_DELAY:-10}"

ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS="run_accession,study_accession,sample_accession,experiment_accession,sample_title,library_name,library_strategy,library_layout,instrument_platform,instrument_model,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp"

command -v curl >/dev/null || { echo "ERROR: curl not found" >&2; exit 1; }
[[ -f "$TSV" ]] || { echo "ERROR: dataset table not found: $TSV" >&2; exit 1; }
mkdir -p "$RAW_DIR"

# ---------------------------------------------------------- dataset table ---
# Emit one line per selected dataset: "dataset_id repository accessions",
# SPACE separated with "-" standing in for empty fields. Columns are resolved
# by header name so the table can be reordered or extended freely.
WANT="$(IFS=,; echo "$*")"

select_datasets() {
  awk -F'\t' -v want="$WANT" '
    function trim(s) { gsub(/^[ \r]+|[ \r]+$/, "", s); return s }
    function field(name) { return trim($(col[name])) }
    NR == 1 {
      for (i = 1; i <= NF; i++) col[trim($i)] = i
      n = split("dataset_id include study_accession repository dataset_accession", req, " ")
      for (i = 1; i <= n; i++) if (!(req[i] in col)) {
        print "ERROR: column \"" req[i] "\" missing from dataset table" > "/dev/stderr"; bad = 1
      }
      if (bad) exit 2
      if (want != "") { n = split(want, w, ","); for (i = 1; i <= n; i++) wanted[w[i]] = 1 }
      next
    }
    /^[ \t\r]*$/ || /^#/ { next }
    {
      id = field("dataset_id")
      if (id == "") next
      if (id ~ /[^A-Za-z0-9._-]/) {
        print "ERROR: dataset_id \"" id "\" may only contain letters, digits, . _ -" > "/dev/stderr"; exit 2
      }
      if (id in seen) { print "ERROR: duplicate dataset_id \"" id "\"" > "/dev/stderr"; exit 2 }
      seen[id] = 1

      if (want != "") { if (!(id in wanted)) next; found[id] = 1 }
      else if (tolower(field("include")) !~ /^(true|t|yes|y|1)$/) next

      acc = field("dataset_accession")
      if (acc == "" || acc == "-" || toupper(acc) == "NA") acc = field("study_accession")
      gsub(/[ \t]/, "", acc); gsub(/,/, ";", acc)
      repo = field("repository"); gsub(/[ \t]/, "_", repo)
      if (acc == "") acc = "-"; if (repo == "") repo = "-"
      print id, repo, acc
    }
    END {
      if (want != "") for (id in wanted) if (!(id in found)) {
        print "ERROR: dataset_id \"" id "\" not in dataset table" > "/dev/stderr"; exit 2
      }
    }' "$TSV"
}

DATASETS="$(select_datasets)"
if [[ -z "$DATASETS" ]]; then
  echo "nothing to do: no dataset in $(basename "$TSV") has include=true"
  exit 0
fi

# --------------------------------------------------------------- manifests ---
fetch_manifest() {   # fetch_manifest <accessions;...> <out.tsv>
  local accs="$1" out="$2" acc first=1 n
  : > "$out.tmp"
  IFS=';' read -r -a acc_list <<< "$accs"
  for acc in "${acc_list[@]}"; do
    [[ -n "$acc" ]] || continue
    if ! curl -fsSL --retry 5 --retry-delay 5 --get "$ENA_API" \
          --data-urlencode "accession=$acc" \
          --data-urlencode "result=read_run" \
          --data-urlencode "fields=$FIELDS" \
          --data-urlencode "format=tsv" \
          --data-urlencode "limit=0" \
          -o "$out.part"; then
      echo "  [WARN] ENA query failed for $acc" >&2
      rm -f "$out.part" "$out.tmp"
      return 1
    fi
    n=$(awk 'NR>1 && NF' "$out.part" | wc -l | tr -d ' ')
    if [[ "$n" -eq 0 ]]; then
      echo "  [WARN] ENA returned no runs for $acc" >&2
    elif [[ $first -eq 1 ]]; then
      cat "$out.part" >> "$out.tmp"; first=0
    else
      awk 'NR>1' "$out.part" >> "$out.tmp"
    fi
  done
  rm -f "$out.part"
  if [[ $first -eq 1 ]]; then rm -f "$out.tmp"; return 1; fi
  # the same run can be reached through several accessions: keep it once
  awk -F'\t' 'NR==1{print; next} !seen[$1]++' "$out.tmp" > "$out"
  rm -f "$out.tmp"
}

JOBLIST="$RAW_DIR/.jobs.txt"
: > "$JOBLIST"
SKIPPED=""
N_DATASETS=0

while read -r ID REPO ACCS; do
  N_DATASETS=$((N_DATASETS + 1))
  case "$(echo "$REPO" | tr '[:lower:]' '[:upper:]')" in
    ENA|SRA|DDBJ|INSDC|NCBI|EBI) ;;
    *) echo "[skip]     $ID: repository \"$REPO\" is not supported -- download it into $RAW_DIR/$ID/ by hand"
       SKIPPED="$SKIPPED $ID"; continue ;;
  esac
  if [[ "$ACCS" == "-" ]]; then
    echo "[skip]     $ID: neither dataset_accession nor study_accession is set"
    SKIPPED="$SKIPPED $ID"; continue
  fi

  DIR="$RAW_DIR/$ID"
  MANIFEST="$DIR/runs.tsv"
  mkdir -p "$DIR/.done"
  if [[ "$REFRESH" == "1" || ! -s "$MANIFEST" || "$(cat "$DIR/.accession" 2>/dev/null)" != "$ACCS" ]]; then
    echo "[manifest] $ID: querying ENA for $ACCS ..."
    if ! fetch_manifest "$ACCS" "$MANIFEST"; then
      echo "[skip]     $ID: no runs found on ENA for $ACCS" >&2
      SKIPPED="$SKIPPED $ID"; continue
    fi
    echo "$ACCS" > "$DIR/.accession"
  fi

  # resolve columns by header name -- ENA does not guarantee field order
  awk -F'\t' -v id="$ID" -v want="$PLATFORM" -v max="$MAX_RUNS" '
    NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    {
      plat = $(h["instrument_platform"]); nb = split($(h["fastq_bytes"]), b, ";")
      for (i = 1; i <= nb; i++) gb[plat] += b[i]
      runs[plat]++
      if (want != "ALL" && plat != want) next
      if (max > 0 && ++taken > max) next
      ftp = $(h["fastq_ftp"]); md5 = $(h["fastq_md5"]); byt = $(h["fastq_bytes"])
      gsub(/[ \t\r]/, "", ftp); gsub(/[ \t\r]/, "", md5); gsub(/[ \t\r]/, "", byt)
      if (ftp == "") ftp = "-"; if (md5 == "") md5 = "-"; if (byt == "") byt = "-"
      print id, $(h["run_accession"]), ftp, md5, byt >> jobs
    }
    END {
      for (p in runs) printf "           %-28s %-18s %5d runs %9.2f GB\n", id, p, runs[p], gb[p] / 1073741824
    }' jobs="$JOBLIST" "$MANIFEST"
done <<< "$DATASETS"

N_SEL=$(wc -l < "$JOBLIST" | tr -d ' ')
SEL_GB=$(awk '{n = split($5, a, ";"); for (i = 1; i <= n; i++) s += a[i]} END {printf "%.2f", s / 1073741824}' "$JOBLIST")
echo "[select]   PLATFORM=$PLATFORM MAX_RUNS=$MAX_RUNS -> $N_SEL runs, ~${SEL_GB} GB"
echo "[dest]     $RAW_DIR"
AVAIL_GB=$( { df -Pg "$RAW_DIR" 2>/dev/null || df -BG "$RAW_DIR"; } | awk 'NR==2{gsub("G","",$4); print $4}')
echo "[disk]     ~${AVAIL_GB} GB free on target volume"

if [[ "$N_SEL" -eq 0 ]]; then
  echo "ERROR: no runs selected (see breakdown above)" >&2
  exit 1
fi
if [[ "$DRY_RUN" == "1" ]]; then
  echo "[dry-run] stopping here; selected runs are listed in $JOBLIST"
  exit 0
fi

# ---------------------------------------------------------------- download ---
md5of() {
  if command -v md5sum >/dev/null; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi
}
sizeof() { wc -c < "$1" | tr -d ' '; }

# ENA has no FASTQ for some (mostly very recent) SRA runs; sra-tools can
# still convert them from the SRA archive.
fetch_with_sra_tools() {
  local run="$1" dir="$2" tmp="$2/.tmp_$1" f
  if ! command -v fasterq-dump >/dev/null; then
    echo "[WARN] $run: no FASTQ on ENA -- install sra-tools (fasterq-dump) or fetch it by hand" >&2
    return 1
  fi
  echo "[sra]  $run: no FASTQ on ENA, converting with fasterq-dump"
  rm -rf "$tmp"; mkdir -p "$tmp"
  if ! fasterq-dump --split-files --outdir "$tmp" --temp "$tmp" "$run" >/dev/null; then
    echo "[FAIL] $run: fasterq-dump error" >&2
    rm -rf "$tmp"; return 1
  fi
  for f in "$tmp"/*.fastq; do
    [[ -e "$f" ]] || continue
    gzip "$f"
    mv "$f.gz" "$dir/"
  done
  rm -rf "$tmp"
}

fetch_run() {
  local id="$1" run="$2" ftp_field="$3" md5_field="$4" bytes_field="$5"
  local dir="$RAW_DIR/$id"

  [[ -f "$dir/.done/$run" ]] && { echo "[skip] $id/$run"; return 0; }
  mkdir -p "$dir/.done"

  if [[ "$ftp_field" == "-" || -z "$ftp_field" ]]; then
    fetch_with_sra_tools "$run" "$dir" || return 1
    touch "$dir/.done/$run"
    echo "[done] $id/$run"
    return 0
  fi
  [[ "$md5_field"   == "-" ]] && md5_field=""
  [[ "$bytes_field" == "-" ]] && bytes_field=""

  local urls=() md5s=() sizes=()
  IFS=';' read -r -a urls  <<< "$ftp_field"
  [[ -n "$md5_field"   ]] && IFS=';' read -r -a md5s  <<< "$md5_field"
  [[ -n "$bytes_field" ]] && IFS=';' read -r -a sizes <<< "$bytes_field"

  local i url fname want want_size got part
  for i in "${!urls[@]}"; do
    # ENA reports paths as bare "ftp.sra.ebi.ac.uk/vol1/...". Fetch them over
    # HTTPS: the same paths are served there, and anonymous FTP often fails
    # with curl (9) "Server denied you to change to the given directory".
    url="${urls[$i]}"
    url="${url#ftp://}"; url="${url#https://}"; url="${url#http://}"
    url="$PROTO://$url"
    fname="$(basename "$url")"
    part="$dir/$fname.part"
    want="${md5s[$i]:-}"
    want_size="${sizes[$i]:-}"

    # keep an existing file only if it verifies; otherwise re-fetch it
    if [[ -f "$dir/$fname" ]]; then
      if [[ -z "$want" || "$(md5of "$dir/$fname")" == "$want" ]]; then
        echo "[ok]   $id/$fname (cached)"
        continue
      fi
      echo "[bad]  $id/$fname failed md5, re-downloading"
      rm -f "$dir/$fname"
    fi

    # a leftover .part is only resumable if it is a strict prefix of the target
    if [[ -f "$part" && -n "$want_size" ]] && (( $(sizeof "$part") >= want_size )); then
      rm -f "$part"
    fi

    echo "[get]  $id/$fname"
    # -sS: parallel progress meters interleave into noise; errors still print
    if ! curl -fsSL --retry "$RETRIES" --retry-delay "$RETRY_DELAY" --retry-all-errors \
              -C - --connect-timeout 30 -o "$part" "$url"; then
      echo "[FAIL] $id/$fname download error" >&2
      return 1
    fi

    if [[ -n "$want" ]]; then
      got="$(md5of "$part")"
      if [[ "$got" != "$want" ]]; then
        echo "[FAIL] $id/$fname md5 mismatch (got $got, want $want) -- discarded" >&2
        rm -f "$part"
        return 1
      fi
    fi
    mv "$part" "$dir/$fname"
  done

  touch "$dir/.done/$run"
  echo "[done] $id/$run"
}

FAILLOG="$RAW_DIR/failed_runs.txt"
: > "$FAILLOG"
export -f fetch_run fetch_with_sra_tools md5of sizeof
export RAW_DIR FAILLOG RETRIES RETRY_DELAY PROTO

# -n5 hands each worker exactly one run's five space-separated fields as
# $1..$5, which avoids -I{} (BSD xargs mangles tabs inside the replacement).
xargs -P "$JOBS" -n 5 "${BASH:-bash}" -c '
    fetch_run "$1" "$2" "$3" "$4" "$5" || echo "$1/$2" >> "$FAILLOG"
  ' _ < "$JOBLIST"

# ------------------------------------------------------------------ report ---
echo
echo "==== download summary ===="
awk '{n[$1]++} END {for (id in n) print id, n[id]}' "$JOBLIST" | sort | while read -r ID N; do
  DONE=0
  while read -r id run _; do
    [[ "$id" == "$ID" && -f "$RAW_DIR/$id/.done/$run" ]] && DONE=$((DONE + 1))
  done < "$JOBLIST"
  printf '%-32s %5d / %-5d runs complete\n' "$ID" "$DONE" "$N"
done
[[ -n "$SKIPPED" ]] && echo "skipped:$SKIPPED"
rm -f "$JOBLIST"
if [[ -s "$FAILLOG" ]]; then
  echo "failed:   $(wc -l < "$FAILLOG" | tr -d ' ') run(s), listed in $FAILLOG"
  echo "          re-run this script to retry them"
  exit 1
fi
rm -f "$FAILLOG"
echo "all selected runs downloaded and md5-verified -> $RAW_DIR"
