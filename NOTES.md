# Overnight run log

Running log for the unattended agent cycles (this instance: Adharv's machine). Honest,
no-overclaiming notes — see `CLAUDE.md` invariants and `PROGRESS.md` for the source of truth
on hypothesis status.

## Cycle 1 — 2026-07-08

- **Blocker resolved:** this machine had no Python 3.11/3.12 (Makefile requires it); installed
  Python 3.12 via Homebrew, created `.venv`, installed `requirements.txt`.
- **Sanity check:** `make test` equivalent (`pytest -m "not network"`) — 30/30 passed, environment
  healthy on a fresh clone.
- **Hygiene fix:** `data/visium/` (squidpy's raw download cache, ~47MB) wasn't covered by
  `.gitignore`'s data-exclusion patterns; added it so raw re-downloadable data doesn't get
  committed by accident.
- **Unblocked PROGRESS.md's stated next step:** kicked off
  `PYTHONPATH=. .venv/bin/python scripts/main.py --dataset breast_cancer` (the tumor sweep for
  H1/H2, previously "blocked on tumor sweep" per `PROGRESS.md`). 30 specs, 3798 spots x 36601
  genes. Running in background; each spec's Moran's I permutation step is the long pole, so this
  is expected to take on the order of hours, not minutes. No results yet — nothing to report on
  H1/H2 for tumor tissue this cycle. Will run `scripts/analyze_dataset.py --dataset breast_cancer`
  once the sweep aggregate exists, and report actual numbers (with permutation p-values and
  bootstrap CIs, per the repo's rigor guardrails) in the next entry — not before.

## Cycle 2 — 2026-07-09 (this instance: Sparsh's machine)

- **Root-caused the tumor sweep failures.** A breast_cancer sweep on this machine finished but
  6/30 specs *failed* (logged in the failure manifest, not silently): `cell_ranger` HVG hit
  duplicate quantile-bin edges, and `total`-count normalization overflowed the seurat flavor's
  dispersion binning. Cause: raw 10x tumor matrix has 36601 genes, ~17k all-zero. Mouse-brain
  demos ship pre-filtered so Year-1 never hit this. (Note: this was a *pipeline* bug on real
  tumor data, independent of the Cycle-1 venv/gitignore work.)
- **Fix (leakage-free, Year-1 unchanged):** minimal shared gene QC `filter_genes(min_cells=10)`
  for raw 10x datasets only (36601→19690), plus a clip fallback in HVG selection that fires only
  when the standard binning fails. Verified `norm_total` now runs end-to-end (5 clusters, 5000
  SVGs); re-ran the full 30-spec sweep — `hvg_cellranger` now succeeds (10 clusters).
- **Also landed this cycle:** Geary's C as a second SVG method (config axis, Task 2), the
  fragility-aware reporting demo (Task 8), and downstream calibration (Brier/ECE, Task 7).
- **Pending:** full tumor sweep re-running; will run `analyze_dataset.py --dataset breast_cancer`
  for tumor H1 (stakes gap) and H2 (does the dominant axis move off resolution?) on completion.
  No tumor H1/H2 numbers claimed yet.
