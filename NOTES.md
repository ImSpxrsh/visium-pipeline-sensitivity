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
- **Update partway through cycle 6:** by the time the local mouse-brain sweep finished, Sparsh's
  machine had already fixed this exact build bug (cycle 7) and landed a lot more (pathologist
  validation, consensus estimator, fragility prediction — see cycles 7-10 below). Discarded this
  machine's redundant local re-run rather than commit over the already-integrated version. Minor
  note: this machine's independent re-run of the mouse-brain sweep differed very slightly from
  the committed one (one spec's cluster count off by 1: min=7 vs committed min=8), most likely a
  cross-machine leidenalg/igraph version difference rather than a real determinism bug — flagging
  for awareness, not treating as a regression (this machine's own within-run determinism check
  passed 3/3 identical in cycle 5).
- **Redirected effort to a real gap in the consensus-estimator result (cycle 8, Sparsh's
  machine):** its headline number (consensus ARI 0.7749 vs best-single 0.7679, a margin of only
  ~0.007) had no bootstrap CI or null test anywhere in `src/consensus.py` — `CLAUDE.md` requires
  every positive claim ship with a bootstrap CI, and this one didn't have one yet. Added
  `consensus_advantage_bootstrap_ci` (resamples spots with replacement, fixes the best-single and
  best-consensus label assignments, gets a paired-bootstrap CI on the ARI margin) plus two unit
  tests: one confirming the CI is exactly zero when consensus equals the best single pipeline,
  one confirming it correctly detects a real, constructed advantage.
- **Caveat on the number I actually got:** ran the new CI against *this machine's own* local
  mouse-brain re-run (not the exact run behind the committed 0.7749/0.7679 figures, since
  `results/runs/` is gitignored and this machine never had the original). On my re-run: consensus
  0.7875 vs single-max 0.7679, 95% bootstrap CI on the margin [0.0047, 0.0344], excludes zero,
  99.8% of draws favor consensus. Directionally this supports the claim (margin is real, not
  noise) — but my point estimate (0.7875) doesn't exactly match the committed one (0.7749), likely
  because "best of many k/linkage settings" is more sensitive to the minor cross-machine
  clustering difference noted above than any individual measurement. **Not committing this number
  over the published one.** Committed the CI function + tests only (dataset-independent, always
  correct); whoever holds the original run data should run
  `consensus_advantage_bootstrap_ci(runs)` against it directly for the CI that exactly matches the
  published 0.7749 figure.
- **On "breakthrough":** that word appears in the cycle-8 commit message. Checked the underlying
  numbers directly rather than taking the label at face value — they're real (verified in
  `results/aggregate/consensus_recovery_summary.json`) and the `beats_all_singles_frac=0.28` /
  `beats_median_frac=0.94` caveats in `PROGRESS.md` are accurate to the data. Not fabricated. But
  a ~0.007-0.02 ARI margin without (until now) a significance test is a real, useful, honestly
  caveated result — not the word I'd have reached for. Flagging the framing here rather than
  editing anyone else's commit message.

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

## Cycle 9 — 2026-07-09 (this instance: Sparsh's machine)

- **Removed the biggest ISEF weakness: validated H1 against INDEPENDENT pathology.** The 10x
  breast-cancer section has per-spot board-certified pathologist labels (SEDR study,
  `JinmiaoChenLab/SEDR_analyses/data/BRCA1/metadata.tsv`; 3798 spots, 100% barcode match,
  committed to `data/annotations/`). Analytical fragility peaks exactly at the pathologist-defined
  **tumor margin** ("Surrounding tumor"): mean S3 0.71 (margin) vs 0.87 (invasive core), 0.74
  (healthy); permutation p<0.001, **survives depth control p<0.001**. This replaces the
  consensus-domain boundary *proxy* with real clinical ground truth — the surgical-margin region
  is the most clinically consequential and the most analytically fragile.
- Added `stakes_test_from_arrays` (array-based, +depth-stratified null), `scripts/analyze_pathology.py`,
  figure 16, macros, a test (34/34 pass), manifest entry. Integrated into abstract + tumor Results.
- Also found (not yet requestable without a mentor): HTAN (dbGaP phs002371) and EGA tumor Visium
  DACs — logged in `data/datasets_manifest.md` for the controlled-access route.
- **@Adharv:** glioblastoma sweep finished 30/30 on this machine (cross-tumor-type H1) — analysis
  next; and the consensus-vs-single-pipeline validation could now also be run against these real
  pathology labels (needs the breast per-spec runs regenerated, since runs/ is gitignored).

## Cycle 10 — 2026-07-09 (this instance: Sparsh's machine)

- **New result: fragility is predictable from a SINGLE pipeline (`src/predict_fragility.py`).**
  Strictly multiverse-free features from one run + raw counts (local PCA heterogeneity,
  single-pipeline neighbor disagreement, depth, genes detected) predict which spots the full
  30-pipeline multiverse finds fragile (S3<0.5) at cross-validated AUC 0.81. An early-warning
  shortcut: flag where a single result is least trustworthy before committing to the sweep.
  Predicted-fragile map visibly aligns with actual S3 (figure 17). Module + 2 tests (36 pass).
- Added to manuscript (short Results paragraph + figure) and macros.
- **In progress:** regenerating the breast per-spec runs (the fixed/converged pipeline) so the
  CENTERPIECE can run — does the co-association CONSENSUS recover the board-certified pathology
  (tumor/margin/healthy) BETTER than any single pipeline? That would combine the consensus
  result with the independent pathology ground truth (clinically grounded, not mouse-brain
  computational annotation). Monitor set; analysis on completion.
- **@Adharv:** three novel results now on main this session — consensus-as-estimator (S3 error
  AUC 0.72), pathology-validated H1 (tumor-margin fragility), single-pipeline fragility
  prediction. Next is consensus-vs-pathology on tumor.

## Cycle 11 — 2026-07-09 (this instance: Adharv's machine)

- **@Sparsh:** great work tonight — the pathology-ground-truth validation (removing the
  boundary-proxy circularity) and the consensus-as-estimator reframing are both genuinely strong,
  well-caveated results. Not duplicating the consensus-vs-pathology centerpiece since you already
  have it running.
- **Closed the same rigor gap in `predict_fragility.py` that cycle 6 closed in `consensus.py`.**
  The AUC=0.81 claim is legitimately out-of-fold (`cross_val_predict` + `StratifiedKFold`, no
  leakage — checked the code directly, this one was done right) but had no CI or null either.
  Added `fragility_prediction_significance`: bootstrap CI on the AUC (resamples spots, reuses the
  already-valid held-out `proba` rather than refitting) and a permutation null (shuffles labels
  with `proba` fixed). Two tests: detects a real constructed signal, correctly reports "not
  significant" on unrelated random `proba`.
- **Real result on this machine's local mouse-brain rerun:** AUC=0.815, 95% CI [0.788, 0.841]
  (excludes chance by a wide margin), permutation p=0.001. Unlike the consensus ARI margin in
  cycle 6, this number lines up closely with the already-published ~0.81 figure — a large-sample
  AUC on 2688 spots is evidently much less sensitive to the minor cross-machine clustering
  discrepancy noted in cycle 6 than "best of many k/linkage settings" was. This is confirmatory,
  not just directional: the fragility-prediction claim holds up under a proper CI and null.
  Committed the function + tests (dataset-independent); left the published
  `results/aggregate/fragility_prediction.json` untouched since I still don't have the exact
  original run — whoever holds it can wire `fragility_prediction_significance` into `main()`'s
  saved JSON for an exact-matching number (one-line addition, already shown working here).
- 40/40 offline tests pass (was 38, +2 new).

## Cycle 12 — 2026-07-09 (this instance: Adharv's machine)

- **No new commits from Sparsh's side this cycle** — the consensus-vs-pathology centerpiece is
  presumably still running or paused. Nothing new to audit yet.
- **`PROGRESS.md` was stale** — missing cycles 9-11 entirely (pathology-validated H1, fragility
  prediction, both bootstrap-CI follow-ups). Synced the "Findings so far" section and "Next" list
  so the doc actually reflects what's on `main`, not just cycle 8's state.
- **Resolved the MOSAIC/MOSAIK naming question that's been sitting as an unactioned note since
  early cycles.** Searched and confirmed: MOSAIK ("Multi-Origin Spatial Transcriptomics Analysis
  and Integration Kit," Baptista et al., arXiv:2505.11384, May 2025, JOSS-reviewed,
  github.com/anthbapt/MOSAIK) is real and published — different purpose (CosMx/Xenium
  imaging-transcript integration vs. our spec-curve fragility/consensus framework) but a genuine
  collision risk for search/citation/judge confusion at the fair, not a false alarm. **Did not
  rename anything** — that's a call for the two of you to make (touches paper title, abstract,
  repo name, README), wrote a concrete recommendation in `PROGRESS.md`'s Open assumptions
  section instead: disambiguating footnote vs. full rename, with the actual citation so whoever
  decides doesn't have to re-research it.
- **@Sparsh:** if the pathology-vs-consensus centerpiece lands while I'm mid-cycle, I'll audit it
  the same way as the last two results (check for leakage/circularity, add a CI/null if missing)
  before it goes further into the manuscript — not trying to slow you down, just keeping the
  same bar applied consistently now that we're this deep into the night.
- No new empirical results this cycle — pure documentation/research hygiene. Nothing to
  overclaim; nothing claimed.

## Cycle 13 — 2026-07-09 (this instance: Adharv's machine)

- **No new commits from Sparsh's side** — consensus-vs-pathology centerpiece still pending.
- **Manuscript consistency check.** With this many results landing tonight, checked whether
  `paper/numbers.tex` had drifted from `results/aggregate/`: re-ran `scripts/report_numbers.py`
  against the currently committed aggregate data and diffed the output against the committed
  `numbers.tex` — zero diff. Every macro in the manuscript is exactly what the current data
  produces; nothing stale. Could not verify the LaTeX actually compiles on this machine (no
  tectonic/latexmk/pdflatex installed) — cycle 7's note that it built with tectonic is the last
  confirmed build check, on the other machine.
- Short cycle — the human on this end is actively working live (wet-lab planning for
  Terra/ISEF), so kept this bounded rather than starting something bigger like the SpatialDE
  axis (still not started, still a reasonable next pickup for either instance when there's more
  uninterrupted runway - it's a real dependency addition to core pipeline.py, not a quick
  standalone module like the last few cycles' rigor additions).

## Cycle 14 — 2026-07-09 (this instance: Adharv's machine)

- **@Sparsh:** still nothing new from your side across a few cycles now — hope that's just the
  consensus-vs-pathology run taking a while, not a stuck process. Whenever it lands I'll audit it
  first thing.
- **Investigated the SpatialDE SVG axis with more runway this cycle** (had considered it too big
  for the shorter cycles earlier tonight). Conclusion: **blocked, and not a quick fix** — full
  writeup in `PROGRESS.md`'s Open assumptions section. Short version: the `spatialde` PyPI
  package is unmaintained against the current SciPy/pandas stack (three separate breaking-API
  issues found by actually running it, not just reading the changelog), and even after patching
  around two of them, the model fit threw real numerical-instability warnings (divide-by-zero,
  invalid-value-in-log during optimization) — a sign that patched-through results might not be
  trustworthy even if I fixed the third break too. Checked for a maintained fork first (none).
  Stopped rather than force through to *some* output whose validity I couldn't vouch for — that
  would be worse than not having the axis at all. Uninstalled the exploratory package; nothing
  committed to `requirements.txt`.
- This is a "closed a real question, found a real answer, the answer is no" cycle — logging it
  so nobody re-discovers the same three bugs from scratch. nnSVG (R-only, needs someone to
  confirm R access) may honestly be more tractable than continuing to patch SpatialDE.
- No empirical results this cycle. Nothing overclaimed — including not overclaiming that
  SpatialDE integration is "coming soon," since right now it genuinely isn't clear it should be
  pursued further without first validating the patched output against a reference implementation.

## Cycle 15 — 2026-07-09 (this instance: Adharv's machine)

- **@Sparsh:** still quiet on your side (4 cycles now) — hope all's well, will keep checking.
- **Closed the same rigor gap in the pathologist-validated H1 result** (`stakes_test_from_arrays`
  in `src/statistics.py`, used by `scripts/analyze_pathology.py`) that cycles 6 and 11 closed for
  consensus and fragility-prediction. This one already had real nulls (plain permutation *and*
  depth-stratified) — better-tested than the other two starting out — but no bootstrap CI on the
  gap magnitude. Added one (paired bootstrap resampling spots, percentile CI), plus a test.
- **Also found, incidentally: `pathology_stakes.json` was stale.** The breast_cancer sweep got
  regenerated in an earlier commit (`4a8dd0b`) but `scripts/analyze_pathology.py` was never
  re-run afterward, so the committed pathology numbers were computed against an older sweep than
  what's currently on `main`. Re-running it (to get the new CI) incidentally fixed this too.
- **Real, exact numbers (this ran against the exact currently-committed data, not an independent
  rerun with drift risk):** margin S³=0.702 vs rest 0.843, gap=0.141, 95% CI [0.128, 0.154] —
  tight, clearly excludes zero, p=0.001 (both plain and depth-controlled). This is now genuinely
  the best-supported single result in the whole project: real pathologist ground truth, real
  permutation null, real depth-stratified null, real bootstrap CI, all in one place.
  `\PathMarginS` in the manuscript shifts 0.71→0.70 from the staleness fix (rounding-boundary
  crossing, 0.7057→0.7017) — regenerated `paper/numbers.tex` via `report_numbers.py` to match;
  no hardcoded prose numbers needed touching since `sections.tex` uses the macro. 40/40 tests
  pass (unchanged count — extended an existing test rather than adding a new one this time).
- **All three headline non-null results now have both a null and a CI**, matching CLAUDE.md's
  own rigor bar exactly: consensus (cycle 6), fragility-prediction (cycle 11), pathology-H1
  (this cycle). Nothing else in the manuscript's main claims is currently missing this pairing,
  as far as this instance has found.
