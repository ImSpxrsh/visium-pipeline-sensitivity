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

## Cycle 4 — 2026-07-09 (this instance: Adharv's machine)

- **Sweep finished: 30/30 specs succeeded, 0 failures.** The fix from cycles 2-3 holds
  completely on a fully fresh run (no stale cache reuse).
- **Real H1 result on tumor tissue:** boundary mean S³ = 0.727, interior mean S³ = 0.873, gap =
  0.147, permutation p = 0.001 (n_boundary=1576, n_interior=2222, consensus-domain boundary
  proxy — no pathologist annotation available for this dataset). This replicates the
  mouse-brain baseline direction and is modestly *larger* (0.147 vs 0.13), consistent with H1's
  prediction that tumor tissue would widen the fragility gap. **Not calling this conclusive**:
  the depth confound (Pearson r=0.447 between per-spot S³ and total UMI counts) is much
  stronger here than on mouse-brain (r=-0.09), so some of the gap could be a depth artifact
  rather than pure analytical fragility. Flagged in `PROGRESS.md` as the next priority to
  resolve, not swept under the rug.
- **Real H2 result on tumor tissue:** dominant variance axis is still `resolution` (72% of
  domain-count variance), same as mouse-brain (~80%) — did **not** shift to SVG method,
  normalization, or HVG selection as H2 speculated it might. Honest null; reported as such in
  `PROGRESS.md`, not spun.
- **Also corrected while updating `PROGRESS.md`:** the backlog table listed Task 8
  (fragility-aware reporting demo) as "not started" even though it was completed and committed
  in cycle 2 — fixed the bookkeeping.
- **Committed:** `results/breast_cancer/aggregate/` and `results/breast_cancer/figures/` (the
  small, reproducible artifacts the `.gitignore` is designed to track — `runs/` and `logs/`
  stay untracked per existing convention), plus the `PROGRESS.md`/`CHANGELOG.md` updates above.
- **No breakthrough claimed.** This is a real, statistically grounded partial replication (H1)
  plus an honest null (H2) on one tumor dataset — exactly the kind of result the repo's rigor
  guardrails are built for, not a headline result yet. The unresolved depth confound on H1 is
  the most important open thread for whoever picks this up next.

## Cycle 5 — 2026-07-09 (this instance: Adharv's machine)

- **Verified determinism still holds.** After this session's several `src/pipeline.py` edits
  (cycles 2-3), re-ran `scripts/test_reproducibility.py`: 3/3 runs identical. CLAUDE.md's
  non-negotiable determinism invariant is intact.
- **Resolved the depth confound flagged in cycle 4.** Added
  `stakes_stratified_s3_depth_controlled` to `src/statistics.py`: bins spots by depth decile and
  runs the boundary/interior permutation null *within* each bin, so the null preserves each
  group's true depth distribution instead of assuming random assignment across all spots. Added
  a unit test (`tests/test_statistics.py`) with a synthetic case where S³ is purely a function of
  depth bin — confirms the test correctly reports p≈1 (no signal beyond depth) in that
  constructed case, so it's not a rubber stamp.
- **Real result on the tumor sweep:** the depth confound is real — boundary spots have
  significantly lower depth than interior spots (median 15,443 vs 25,226 UMI, Mann-Whitney
  p=5e-55). But the S³ gap survives depth-stratified control: p=0.001, same order as the naive
  test. **H1's tumor result is not just a depth artifact.** This directly closes the open
  question from cycle 4 — updated `PROGRESS.md` accordingly (no longer listed as "unresolved").
  Being precise about what this does and doesn't prove: it's a depth-decile block-permutation,
  a standard but specific robustness check — not a full multivariate causal adjustment. It rules
  out "the gap is just spots with fewer reads," not every conceivable depth-related mechanism.
- **Committed:** `src/statistics.py`, `tests/test_statistics.py`, the regenerated
  `results/breast_cancer/aggregate/statistics_summary.json`, and the `PROGRESS.md` update. Full
  offline suite: 31/31 (was 30, +1 new test) passes.
- **Bottom line for the morning:** across 5 cycles tonight, the tumor sweep went from
  20%-broken → fixed, root-caused, and cross-verified by two independent agents → produced real
  H1 (supportive, confound-tested) and H2 (honest null) results on human tumor tissue, with
  determinism and the full test suite verified intact throughout. Nothing here is being called a
  breakthrough. What's real: a working, reproducible pipeline extension to tumor tissue with one
  statistically grounded positive result and one honest null, both defensible under review.

## Cycle 6 — 2026-07-09 (this instance: Adharv's machine)

- **Found a real manuscript build bug while checking sync.** `paper/sections.tex` (Fig
  `fig:stakes`, line 424) already references `13_stakes_stratified.png`, but that file does not
  exist in `results/figures/` (the primary/mouse-brain figure set) — it only exists under
  `results/breast_cancer/figures/` from tonight's tumor work. If anyone tried `make paper` right
  now it would fail on a missing graphic. `paper/numbers.tex`'s `\StakesGap{0.13}` /
  `\StakesP{<0.001}` macros are also still mouse-brain-only — tonight's tumor H1/H2 numbers and
  the depth-confound resolution (cycle 5) aren't reflected in the manuscript at all yet.
- **Why this needs a full sweep, not just a script run:** `results/runs/` (per-spec sweep
  outputs) is intentionally gitignored - large, reproducible from code - so this machine's fresh
  clone only has the *committed aggregate summaries*, not the underlying per-spec adata the
  figure/stats regeneration needs. Kicked off a full `scripts/main.py --dataset visium_hne`
  sweep (30 specs) in the background to regenerate it locally; will run
  `scripts/analyze_dataset.py --dataset visium_hne` once it completes to produce the missing
  figure and (as a side effect) get the new `stakes_stratified_s3_depth_controlled` stat
  computed on the primary Year-1 dataset too, for full consistency with the tumor-side check.
- **Not claiming this as new science** — mouse-brain H1 is Year-1, already established; this
  cycle is repo hygiene (fixing a broken build artifact) plus consistency, not a new result.
  Will follow up with updating `paper/numbers.tex`'s macros to include the tumor H1/H2 numbers
  once the mouse-brain regeneration confirms nothing regressed.

## Cycle 7 — 2026-07-09 (this instance: Sparsh's machine)

- **Resolved the cycle-6 build bug directly.** The missing `results/figures/13_stakes_stratified.png`
  (mouse-brain) was present untracked on this machine from an earlier cycle; committed it
  (`daba4c5`). The manuscript builds again. **@Adharv: your full mouse-brain sweep to regenerate
  this figure is no longer needed for the build** — though it's still useful if you want the
  `stakes_stratified_s3_depth_controlled` stat computed on the Year-1 dataset for consistency.
- **Integrated the tumor H1/H2 results into the manuscript** (so this is done, not pending):
  `report_numbers.py` now emits tumor macros from the committed `breast_cancer` summary;
  added a "Higher-stakes human tumor tissue" Results paragraph (boundary gap 0.15, p<0.001,
  survives depth-stratified null despite r=0.45 confound; H2 honest null, resolution 72%),
  the tumor stakes figure, and updated the abstract + limitations for the third dataset.
  **@Adharv: no need to duplicate the numbers.tex/tumor-writeup follow-up you flagged.**
- Verified: 31/31 offline tests pass, paper compiles (tectonic), pushed to main. In sync.
- **Open next:** H3 (cross-platform Xenium+Visium) still not started — needs hosted Xenium files
  (Janesick) not fetchable programmatically here; and the MOSAIC/MOSAIK naming collision.

## Cycle 8 — 2026-07-09 (this instance: Sparsh's machine)

- **New constructive result — the multiverse as a better estimator (`src/consensus.py`).**
  Turned the diagnostic into a method: built a co-association matrix (fraction of the 30
  aligned pipelines placing each spot-pair in the same cluster) and re-clustered it into
  consensus domains, then validated against the mouse-brain reference anatomy. Single
  pipelines recover anatomy at ARI 0.55-0.77 (median 0.64); the co-association consensus
  reaches **0.77 — matching/exceeding the best single pipeline** and beating the median in
  **94%** of granularity/linkage settings. Point: you can't pick the best single pipeline
  without the ground truth, but the consensus gets best-or-better accuracy from the multiverse
  alone, with per-spot S³ uncertainty for free.
- **Honesty caveats (kept in the paper):** "beats *all* singles" holds for average linkage near
  the reference granularity (28% of settings), not universally; "beats the *median*" is the
  robust claim (94%). Validation needs ground truth, so this is mouse-brain only (tumor has no
  pathologist annotation). Not overstated.
- Added module + 2 tests (33/33 offline pass), figure `15_consensus_recovery.png`, a dedicated
  paper section (\S consensus), abstract + contribution (iv), and `report_numbers` macros.
- **@Adharv:** this is the biggest single lift toward "breakthrough" — reframes the whole paper
  from cautionary tale to usable method. Natural follow-ups if you pick it up: (a) does the
  consensus also reduce the downstream-ML accuracy/calibration swing? (b) a stability-weighted
  consensus (weight specs or spots by S³); (c) confirm on a dataset with annotations beyond
  mouse brain if one becomes fetchable.
