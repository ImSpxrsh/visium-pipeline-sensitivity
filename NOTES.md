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

## Cycle 3 — 2026-07-09 (this instance: Adharv's machine)

- **Converged with Sparsh's fix.** This machine independently root-caused the same 6 tumor-sweep
  failures in parallel (same conclusion: all-zero genes breaking scanpy's HVG binning) and had
  applied its own `filter_genes(min_counts=1)` inside `_normalize()`. `git pull --rebase` merged
  both changes cleanly (different code paths, no conflict). Verified this machine's fix directly:
  it resolved `hvg_cellranger` but *not* `norm_total` — traced that one further and confirmed the
  real mechanism: `seurat`/`cell_ranger` flavors call `expm1` internally assuming log-scale input;
  on this tissue `IGKC`/`IGLC1`/`IGLC2`/`IGHG3` (immunoglobulin/plasma-cell markers — real immune
  infiltration signal, not noise) reach ~2000x target-sum in single spots, so `expm1` overflows to
  infinity. The same combo passes on mouse-brain data by luck (no gene there is that extreme), not
  because it's actually valid — a genuine, tissue-dependent numerical-stability finding, not a bug
  to paper over.
- **Deduplicated:** Sparsh's `min_cells=10` gene QC (in `load_dataset`, gated to raw 10x datasets)
  strictly subsumes this machine's `min_counts=1` filter, so removed the redundant one rather than
  leave two overlapping QC mechanisms in the codebase. Offline tests (30/30) still pass after.
- **Correctness catch:** this machine had a *stale* local `results/breast_cancer/` from Cycle 1's
  pre-fix sweep (24/30 specs, computed on the un-filtered 36601-gene universe). The sweep script's
  incremental/content-hashed caching would have silently reused those 24 as "already done" and
  only computed the 6 new ones on the *now-filtered* 19690-gene universe — mixing two different
  gene universes inside one S³ aggregate, which would have been a real correctness bug, not just
  staleness. Moved the stale directory aside (not deleted) and started a fully fresh 30-spec sweep.
  Worth checking for the same trap if any other machine has a pre-fix local `results/breast_cancer/`
  lying around.
- **Kicked off:** fresh `scripts/main.py --dataset breast_cancer` sweep on the deduplicated,
  verified pipeline. Running in background; will report actual H1/H2 numbers (permutation
  p-values, bootstrap CIs) via `analyze_dataset.py` once it completes — none claimed yet.
