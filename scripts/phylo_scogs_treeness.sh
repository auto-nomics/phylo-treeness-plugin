#!/usr/bin/env bash
# phylo_scogs_treeness — per-gene alignment, trimming, tree inference, and
# treeness over a tarball of single-copy-ortholog protein FASTA files.
#
# Contract (see manifest doc): input 0 is a tar archive whose members are
# per-gene FASTA files named <gene_id>.faa (or .fasta). Outputs are a
# treeness TSV (one row per input gene) and a JSON summary with medians.
# Per-gene failures are recorded and never abort the node.

set -eu

fail() {
    echo "phylo_scogs_treeness: error: $*" >&2
    exit 1
}

# ---------------------------------------------------------------------------
# 1. Runtime-provided paths and parameters (validate early, fail fast).
# ---------------------------------------------------------------------------
input_tar="${AUTONOMICS_INPUT0:-}"
out_tsv="${AUTONOMICS_OUTPUT0:-}"
out_json="${AUTONOMICS_OUTPUT1:-}"
workdir="${AUTONOMICS_WORKDIR:-/work}"

[[ -n "$input_tar" ]] || fail "AUTONOMICS_INPUT0 is not set"
[[ -n "$out_tsv" ]] || fail "AUTONOMICS_OUTPUT0 is not set"
[[ -n "$out_json" ]] || fail "AUTONOMICS_OUTPUT1 is not set"
[[ -f "$input_tar" ]] || fail "input tar not found: $input_tar"

threads="${AUTONOMICS_THREADS:?AUTONOMICS_THREADS is not set}"
model="${AUTONOMICS_MODEL:?AUTONOMICS_MODEL is not set}"
min_taxa="${AUTONOMICS_MIN_TAXA:?AUTONOMICS_MIN_TAXA is not set}"
seed="${AUTONOMICS_SEED:?AUTONOMICS_SEED is not set}"

case "$threads" in ''|*[!0-9]*) fail "threads must be a positive integer: $threads" ;; esac
[[ "$threads" -ge 1 ]] || fail "threads must be >= 1"
case "$min_taxa" in ''|*[!0-9]*) fail "min_taxa must be an integer >= 2: $min_taxa" ;; esac
[[ "$min_taxa" -ge 2 ]] || fail "min_taxa must be >= 2"
case "$seed" in ''|*[!0-9]*) fail "seed must be a non-negative integer: $seed" ;; esac
[[ -n "$model" ]] || fail "model must be a non-empty string"

# ---------------------------------------------------------------------------
# 2. Workspace layout under the writable mount.
# ---------------------------------------------------------------------------
genes_dir="$workdir/genes"
aln_dir="$workdir/alignments"
tree_dir="$workdir/trees"
log_dir="$workdir/logs"
mkdir -p "$genes_dir" "$aln_dir" "$tree_dir" "$log_dir"

# ---------------------------------------------------------------------------
# 3. Extract and validate the gene tarball.
# ---------------------------------------------------------------------------
if ! tar -xf "$input_tar" -C "$genes_dir" 2>"$log_dir/untar.err"; then
    fail "input is not a readable tar archive: $(cat "$log_dir/untar.err" | tail -1)"
fi

mapfile -t gene_files < <(find "$genes_dir" -maxdepth 1 -type f \( -name '*.faa' -o -name '*.fasta' \) | sort)
[[ "${#gene_files[@]}" -ge 1 ]] || fail "tar archive contains no .faa/.fasta gene files"

# ---------------------------------------------------------------------------
# 4. Version banner (recorded into the summary for provenance).
# ---------------------------------------------------------------------------
mafft_version="$(mafft --version 2>&1 | grep -oE '[0-9]+\.[0-9]+' | head -1)"
clipkit_version="$(clipkit --version 2>/dev/null | tail -1 | awk '{print $NF}')"
iqtree_version="$(iqtree --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
# phykit's CLI has no --version; read the installed distribution metadata.
phykit_version="$(python3 -c 'from importlib.metadata import version; print(version("phykit"))' 2>/dev/null || echo unknown)"

# ---------------------------------------------------------------------------
# 5. Per-gene pipeline: MAFFT -> ClipKIT -> IQ-TREE -> PhyKIT treeness.
# ---------------------------------------------------------------------------
tsv_tmp="$workdir/treeness.rows.tsv"
: > "$tsv_tmp"

count_seqs() {
    # Number of FASTA records; a record needs a '>' header line.
    grep -c '^>' "$1" || true
}

for fa in "${gene_files[@]}"; do
    gene_id="$(basename "$fa")"
    gene_id="${gene_id%.faa}"
    gene_id="${gene_id%.fasta}"

    n_seqs="$(count_seqs "$fa")"
    if [[ "$n_seqs" -lt "$min_taxa" ]]; then
        printf '%s\t%s\t\ttoo_few_taxa\n' "$gene_id" "$n_seqs" >> "$tsv_tmp"
        continue
    fi

    aln="$aln_dir/$gene_id.mafft.fasta"
    trim="$aln_dir/$gene_id.mafft.clipkit.fasta"
    prefix="$tree_dir/$gene_id"

    if ! mafft --auto --thread "$threads" "$fa" > "$aln" 2> "$log_dir/$gene_id.mafft.err"; then
        printf '%s\t%s\t\tfailed_mafft\n' "$gene_id" "$n_seqs" >> "$tsv_tmp"
        continue
    fi
    if ! clipkit "$aln" -o "$trim" > "$log_dir/$gene_id.clipkit.out" 2>&1; then
        printf '%s\t%s\t\tfailed_clipkit\n' "$gene_id" "$n_seqs" >> "$tsv_tmp"
        continue
    fi
    if ! iqtree -s "$trim" -m "$model" -T "$threads" -seed "$seed" \
            --prefix "$prefix" -quiet > "$log_dir/$gene_id.iqtree.out" 2>&1; then
        printf '%s\t%s\t\tfailed_iqtree\n' "$gene_id" "$n_seqs" >> "$tsv_tmp"
        continue
    fi

    treefile="$prefix.treefile"
    if [[ ! -s "$treefile" ]]; then
        printf '%s\t%s\t\tfailed_no_treefile\n' "$gene_id" "$n_seqs" >> "$tsv_tmp"
        continue
    fi

    # phykit treeness prints the ratio on stdout; take the last float token.
    treeness="$(phykit treeness "$treefile" 2>"$log_dir/$gene_id.phykit.err" \
        | grep -oE '[-+]?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?' | tail -1)" || treeness=""
    if [[ -z "$treeness" ]]; then
        printf '%s\t%s\t\tfailed_treeness\n' "$gene_id" "$n_seqs" >> "$tsv_tmp"
        continue
    fi

    printf '%s\t%s\t%s\tok\n' "$gene_id" "$n_seqs" "$treeness" >> "$tsv_tmp"
done

# ---------------------------------------------------------------------------
# 6. Emit the TSV output (always, even when every gene failed).
# ---------------------------------------------------------------------------
{
    printf 'gene_id\tn_taxa\ttreeness\tstatus\n'
    sort -t "$(printf '\t')" -k1,1 "$tsv_tmp"
} > "$out_tsv"

# ---------------------------------------------------------------------------
# 7. Emit the JSON summary (always).
# ---------------------------------------------------------------------------
python3 - "$tsv_tmp" "$out_json" \
    "$threads" "$model" "$min_taxa" "$seed" \
    "$mafft_version" "$clipkit_version" "$iqtree_version" "$phykit_version" <<'PY'
import json
import statistics
import sys
from pathlib import Path

rows_path, out_path = Path(sys.argv[1]), sys.argv[2]
threads, model, min_taxa, seed = sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
versions = {
    "mafft": sys.argv[7],
    "clipkit": sys.argv[8],
    "iqtree": sys.argv[9],
    "phykit": sys.argv[10],
}

rows = []
for line in rows_path.read_text().splitlines():
    gene_id, n_taxa, treeness, status = (line.split("\t") + ["", "", "", ""])[:4]
    rows.append(
        {
            "gene_id": gene_id,
            "n_taxa": int(n_taxa) if n_taxa else None,
            "treeness": float(treeness) if treeness else None,
            "status": status,
        }
    )

ok = [row for row in rows if row["status"] == "ok" and row["treeness"] is not None]
values = [row["treeness"] for row in ok]
summary = {
    "input_genes": len(rows),
    "trees_built": len(ok),
    "treeness_computed": len(values),
    "skipped_too_few_taxa": sum(1 for r in rows if r["status"] == "too_few_taxa"),
    "failed": sum(1 for r in rows if r["status"].startswith("failed")),
    "median_treeness": statistics.median(values) if values else None,
    "mean_treeness": (sum(values) / len(values)) if values else None,
    "n_zero_treeness": sum(1 for v in values if v == 0.0),
    "trees_with_3_taxa": sum(1 for r in ok if r["n_taxa"] == 3),
    "trees_with_4_taxa": sum(1 for r in ok if r["n_taxa"] == 4),
    "params": {
        "threads": int(threads),
        "model": model,
        "min_taxa": int(min_taxa),
        "seed": int(seed),
    },
    "versions": versions,
}
out_path_text = sys.argv[2]
Path(out_path_text).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
PY
