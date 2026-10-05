# phylo-treeness

Autonomics container plugin for per-gene phylogenetic tree-quality analysis
over single-copy ortholog (SCOG) protein sets — the MAFFT → ClipKIT →
IQ-TREE → PhyKIT treeness pipeline used by BixBench-style comparative
genomics questions (e.g. "median treeness difference between fungi and
animals").

## Node

### `phylo_scogs_treeness`

**Input port 0 — `scogs_genes` (tar).** A tar archive whose members are
per-gene multi-FASTA protein files named `<gene_id>.faa` (or `.fasta`). Each
file holds one gene's ortholog sequences across the taxa of one clade; the
gene identifier is the file stem. Sequence headers are not interpreted —
treeness only depends on the inferred tree topology and branch lengths. A
one-clade-per-run convention: run the node once per clade (fungi, animals)
and combine the summaries downstream.

**Output 0 — `treeness.tsv`** (`phylo_scogs_treeness_tsv`), one row per
input gene, sorted by `gene_id`:

| column | meaning |
| --- | --- |
| `gene_id` | input file stem |
| `n_taxa` | sequence count in the input FASTA |
| `treeness` | internal-branch-length sum ÷ total tree length; empty on non-ok rows |
| `status` | `ok`, `too_few_taxa`, `failed_mafft`, `failed_clipkit`, `failed_iqtree`, `failed_no_treefile`, `failed_treeness` |

**Output 1 — `treeness_summary.json`** (`phylo_scogs_treeness_json`): counts
(`input_genes`, `trees_built`, `treeness_computed`,
`skipped_too_few_taxa`, `failed`, `n_zero_treeness`, `trees_with_3_taxa`,
`trees_with_4_taxa`), `median_treeness` / `mean_treeness` across
successfully treed genes (null when none), resolved `params`, and pinned
tool `versions`.

## Parameters

| param | default | semantics |
| --- | --- | --- |
| `threads` | 2 | Thread cap for MAFFT and IQ-TREE (`--thread` / `-T`). Keep low on shared hosts — alignment and tree-search peak memory scale with threads. |
| `model` | `TEST` | IQ-TREE substitution model. `TEST` runs ModelFinder per gene; pass a fixed model (e.g. `LG+G4`) for faster, fully deterministic runs. |
| `min_taxa` | 3 | Genes with fewer sequences are recorded `too_few_taxa` and never aligned. |
| `seed` | 12345 | Fixed IQ-TREE random seed for reproducible inference. |

## Scientific contract notes

- **Treeness** follows PhyKIT's definition: the sum of internal branch
  lengths divided by the total tree length of the inferred (unrooted) tree.
- A **3-taxon unrooted tree has no internal branches**, so its treeness is
  exactly `0.0`. Clade medians therefore depend on taxon completeness:
  comparing clades with different fractions of 3-taxon genes biases the
  comparison. Use the `n_taxa` column (or `trees_with_3_taxa` /
  `trees_with_4_taxa` in the summary) to restrict to complete-taxa trees
  when the scientific question requires it.
- Per-gene tool failures are recorded with `failed_*` status; the node
  still writes both outputs. An input tar with no FASTA members, or a
  non-tar input, fails the node.
- The container has **no network** (`isolated`) and a **read-only rootfs**;
  all scratch work happens under the runtime-provided workdir.

## Image provenance

Digest-pinned base and exact bioconda pins:

| component | version | source |
| --- | --- | --- |
| base | `condaforge/miniforge3:latest@sha256:be01a3e63d4af20d8a9937d1e12fa299d1c572edf601b961180a59c78bd7250d` | docker.io/condaforge |
| mafft | 7.526 | bioconda |
| clipkit | 2.14.0 | bioconda |
| iqtree | 3.1.4 | bioconda |
| phykit | 2.7.0 | bioconda |
| python | 3.12 | conda-forge |

Build: `podman build -t ghcr.io/auto-nomics/autonomics/phylo-treeness:0.1.0 .`

**Release status:** the `image.reference` digest is the **registry-returned
repository digest** of the GHCR push
(`sha256:7e938250ab0fe718e0fcfad120293f1a08814c67aecce3b0cb4b86c2f853a3cb`),
verified resolvable via
`podman pull …@sha256:<digest>`.

## Resource overrides

`cpus = 2.0`, `memory = "8Gi"`, `pids_limit = 512`, `shm_size = "1Gi"` are
explicit (matching the shared-host profile of other installed plugins);
`network` and `read_only_rootfs` keep the hardened defaults. The 4-hour
`timeout_secs` accommodates ~250-gene clades at 2 threads on one host.

## Tests

```bash
# Image + script smoke (happy path, too_few_taxa, malformed tar, digest pull)
./test_phylo_scogs_treeness.sh

# Manifest contract test (from the autonomics checkout)
NODE_PLUGINS_ROOT=/mnt/projects/node-plugins \
  cargo test -p container-plugin --test phylo_treeness_plugin
```

Fixture coverage: 4-taxon genes (treeness > 0), a 3-taxon gene (treeness
exactly 0), a 2-taxon gene (`too_few_taxa`), and a malformed tar (node
fails with a readable error and writes no outputs).

## Known differences from upstream CLIs

- MAFFT runs with `--auto` (algorithm auto-selection); no manual acurate
  variants are exposed.
- ClipKit runs in its default `gappy` mode.
- IQ-TREE writes per-gene artifact files (`.iqtree`, `.log`, `.ckp.gz`,
  `.model.gz`, `.mldist`, `.bionj`) into the scratch workdir; only the
  `.treefile` is consumed. Scratch is discarded after the run.
- PhyKIT is invoked per gene via its CLI (`phykit treeness`); the value is
  parsed as the last float on stdout.
