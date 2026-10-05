#!/usr/bin/env bash
# Smoke test for the phylo_scogs_treeness node.
#
# Cases: happy path (4- and 3-taxon genes), a below-min_taxa gene, and a
# malformed tar. Runs the exact image by digest by default; PHYLO_RUNNER=local
# executes the script with the host toolchain instead.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
runner="${PHYLO_RUNNER:-container}"
image_tag="${PHYLO_IMAGE:-localhost/autonomics/phylo-treeness:0.1.0}"
image_ref="ghcr.io/auto-nomics/autonomics/phylo-treeness@sha256:7e938250ab0fe718e0fcfad120293f1a08814c67aecce3b0cb4b86c2f853a3cb"
script="$root/scripts/phylo_scogs_treeness.sh"

work="$(mktemp -d)"
chmod 0777 "$work"

# Container-written files are owned by uid 1001; fall back to a rootful
# container delete when the host user cannot remove them.
cleanup_all() {
  rm -rf "$work" 2>/dev/null || \
    podman run --rm --user 0 --volume "$work:/work" --entrypoint bash "$image_tag" \
      -c 'rm -rf /work/* /work/.[!.]* 2>/dev/null || true; rmdir /work 2>/dev/null || true'
}
trap cleanup_all EXIT

# Fixture tarball is materialized at test time from the committed FASTAs so
# the archive itself never needs to be reviewed as a binary artifact.
tar -cf "$work/scogs_genes.tar" -C "$root/fixtures/genes" .
cp "$root/fixtures/malformed.tar" "$work/malformed.tar"
chmod 0644 "$work/scogs_genes.tar" "$work/malformed.tar"

# run_node <case> <host-input> — outputs land at out_<case>.{tsv,json},
# each case gets its own scratch workdir so cases never interact.
run_node() {
  local case="$1" host_input="$2"
  local input="/work/$(basename "$host_input")"
  local tsv="/work/out_${case}.tsv"
  local json="/work/out_${case}.json"
  local wd="/work/wd_${case}"
  rm -rf "$work/wd_${case}" 2>/dev/null || true
  mkdir -p "$work/wd_${case}"
  chmod 0777 "$work/wd_${case}"
  if [[ "$runner" == local ]]; then
    AUTONOMICS_INPUT0="$host_input" \
      AUTONOMICS_OUTPUT0="$work/out_${case}.tsv" \
      AUTONOMICS_OUTPUT1="$work/out_${case}.json" \
      AUTONOMICS_WORKDIR="$work/wd_${case}" \
      AUTONOMICS_THREADS=2 \
      AUTONOMICS_MODEL=TEST \
      AUTONOMICS_MIN_TAXA=3 \
      AUTONOMICS_SEED=12345 \
      bash "$script"
  else
    podman run --rm \
      --network=none --read-only \
      --tmpfs /tmp:rw,size=512m \
      --volume "$work:/work" \
      --volume "$script:/work/.autonomics/script:ro" \
      --env AUTONOMICS_INPUT0="$input" \
      --env AUTONOMICS_OUTPUT0="$tsv" \
      --env AUTONOMICS_OUTPUT1="$json" \
      --env AUTONOMICS_WORKDIR="$wd" \
      --env AUTONOMICS_THREADS=2 \
      --env AUTONOMICS_MODEL=TEST \
      --env AUTONOMICS_MIN_TAXA=3 \
      --env AUTONOMICS_SEED=12345 \
      --entrypoint bash \
      "$image_ref" /work/.autonomics/script
  fi
}

echo "== happy path =="
run_node happy "$work/scogs_genes.tar"

[[ -s "$work/out_happy.tsv" ]]
[[ -s "$work/out_happy.json" ]]
head -1 "$work/out_happy.tsv" | grep -q '^gene_id	n_taxa	treeness	status$'

python3 - "$work/out_happy.tsv" "$work/out_happy.json" <<'PY'
import sys
from pathlib import Path

lines = Path(sys.argv[1]).read_text().splitlines()
header = lines[0].split("\t")
assert header == ["gene_id", "n_taxa", "treeness", "status"], header
by_id = {parts[0]: parts for parts in (line.split("\t") for line in lines[1:])}
assert set(by_id) == {"gene_001", "gene_002", "gene_003", "gene_004"}, set(by_id)
assert by_id["gene_001"][3] == "ok" and float(by_id["gene_001"][2]) > 0.0
assert by_id["gene_002"][3] == "ok" and float(by_id["gene_002"][2]) > 0.0
# 3-taxon unrooted tree: no internal branches, treeness exactly 0.
assert by_id["gene_003"][3] == "ok" and float(by_id["gene_003"][2]) == 0.0
assert by_id["gene_004"][3] == "too_few_taxa" and by_id["gene_004"][2] == ""
assert by_id["gene_004"][1] == "2"

import json
summary = json.loads(Path(sys.argv[2]).read_text())
assert summary["input_genes"] == 4
assert summary["trees_built"] == 3
assert summary["treeness_computed"] == 3
assert summary["skipped_too_few_taxa"] == 1
assert summary["failed"] == 0
assert summary["n_zero_treeness"] == 1
assert summary["trees_with_3_taxa"] == 1
assert summary["trees_with_4_taxa"] == 2
assert abs(summary["median_treeness"] - sorted([0.0, float(by_id["gene_001"][2]), float(by_id["gene_002"][2])])[1]) < 1e-12
assert summary["params"] == {"threads": 2, "model": "TEST", "min_taxa": 3, "seed": 12345}
for tool in ("mafft", "clipkit", "iqtree", "phykit"):
    assert summary["versions"][tool], f"missing {tool} version"
print("happy path OK: median_treeness =", summary["median_treeness"])
PY

echo "== malformed tar =="
if run_node bad "$work/malformed.tar" 2>"$work/bad.err"; then
  echo "malformed tar must fail" >&2
  exit 1
fi
grep -qi "not a readable tar" "$work/bad.err" || {
  echo "unexpected malformed-tar error:" >&2
  cat "$work/bad.err" >&2
  exit 1
}
[[ ! -e "$work/out_bad.tsv" ]]
echo "malformed tar OK"

echo "== digest-addressed image resolves locally =="
if [[ "$runner" == container ]]; then
  podman image exists "$image_ref"
  echo "digest resolution OK"
fi

echo "ALL SMOKE TESTS PASSED"
