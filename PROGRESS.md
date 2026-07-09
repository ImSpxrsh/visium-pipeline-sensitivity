# PROGRESS — MOSAIC Year 2

Maps work to the ISEF hypotheses (H1/H2/H3) and the backlog in `CLAUDE.md`.

## Status
- **Year 1 (foundation, complete):** mouse-brain Visium (H&E + fluorescence), 30-spec sweep, S³,
  ARI, anatomy validation, variance decomposition (resolution ~80%), confound check (r=-0.09),
  downstream ML (5 pp swing), bioRxiv manuscript. All on `main`, CI green.

## Year 2 backlog
| # | Task | Hypothesis | Status |
|---|---|---|---|
| 1 | Human tumor Visium loads + sweep runs unchanged | H1/H2/H3 | **done — 30/30 specs succeed** |
| 2 | SVG-method axis (+Geary's C; SpatialDE/nnSVG later) | H2 | **Geary's C done** |
| 3 | QC/filtering axis | — | minimal shared gene QC for raw 10x (not a per-spec axis yet) |
| 4 | Stakes-stratified S³ (boundary/interior/healthy) | H1 | **done (mouse-brain + tumor)** |
| 5 | Per-dataset Type-II ANOVA (does dominant axis move?) | H2 | **done — null: axis unchanged** |
| 6 | Cross-platform S³ (Xenium+Visium) | H3 | not started |
| 7 | Downstream calibration swing (ECE/Brier) | — | **done** |
| 8 | Fragility-aware reporting demo | — | **done** |

## Findings so far (Year 2)
- **Multiverse consensus is a better estimator (cycle 8, `src/consensus.py`).** A
  co-association consensus across the 30 pipelines recovers mouse-brain reference anatomy at
  ARI 0.77 — matching/exceeding the best single pipeline (0.77) and beating the median (0.64) in
  94% of settings — *without* the ground truth needed to pick the best single pipeline. Per-spot
  S³ predicts where the consensus errs (AUC 0.72), so the estimate is uncertainty-quantified for
  free. Reframes MOSAIC from fragility diagnostic to constructive method. (Honest caveat:
  "beats all singles" holds at reference granularity/average linkage, 28% of settings; "beats
  median" is the robust claim.)
- **H1 (mouse-brain baseline):** boundary spots are significantly more fragile than
  interiors — mean S³ 0.72 vs 0.85, gap 0.13, permutation p<0.001. Establishes the
  method; tumor tissue expected to widen the gap.
- **H1 (breast tumor, real result — 2026-07-09):** boundary mean S³ = 0.727, interior mean
  S³ = 0.873, gap = 0.147, permutation p = 0.001 (3798 spots, consensus-domain boundary proxy —
  no pathologist annotation for this dataset, so labeled as a proxy per the Open assumptions
  below). The gap is directionally consistent with, and modestly larger than, the mouse-brain
  baseline (0.147 vs 0.13) — supports H1 on a second, higher-stakes tissue.
- **H1 depth confound — checked, gap survives (2026-07-09):** the confound is real: boundary
  spots have significantly lower sequencing depth than interior spots (median 15,443 vs 25,226
  UMI, Mann-Whitney p=5e-55; overall S³-vs-depth r=0.447, far above mouse-brain's r=-0.09). Added
  a depth-stratified block-permutation null (`stakes_stratified_s3_depth_controlled` in
  `src/statistics.py`) that bins spots by depth decile and shuffles the boundary/interior label
  only *within* each bin, so the null preserves each group's true depth distribution — directly
  testing whether the gap is just a depth artifact. It isn't: p=0.001, same order as the naive
  test. This is a specific, standard robustness check, not a full multivariate adjustment, so it
  doesn't rule out every possible depth-related mechanism — but it directly answers the concern
  raised after cycle 4. Upgrading the tumor H1 result from "supportive, unresolved caveat" to
  "supportive, confound directly tested and survived."
- **H2 (breast tumor, real result — 2026-07-09):** dominant variance axis is still
  **resolution** (72% of variance), same as mouse-brain (~80%) — the axis did **not** shift to
  SVG method, normalization, or HVG selection as H2 speculated. This is an honest null for the
  "axis changes with tissue type" hypothesis on this one tumor dataset; resolution simply
  dominates less (72% vs 80%). Not spun as a positive result — reported plainly per the repo's
  no-p-hacking invariant.
- **Human tumor loader:** `--dataset breast_cancer` (10x V1_Breast_Cancer_Block_A_Section_1)
  and `glioblastoma` wired via squidpy; no registration needed. Full 30-spec sweep now succeeds
  end to end (fixed in cycles 2-3, see `NOTES.md`).
- **Downstream calibration:** Brier score nearly doubles (0.089→0.159) across pipelines while
  accuracy moves 5 pp — pipeline choice affects calibration more than accuracy.

## Next
- H3 (cross-platform Xenium+Visium S³) and the SVG-method axis expansion (SpatialDE/nnSVG)
  remain not started — no other item is currently blocking them.

## Open assumptions / decisions
- **Datasets:** environment fetches only programmatic, no-registration data. Human tumor via
  `squidpy.datasets.visium(<10x id>)`. DLPFC/her2st/Janesick need hosted files — wire as URLs allow.
- **H1 annotations:** where pathologist labels are unavailable, use a *consensus-domain boundary*
  proxy for "boundary", clearly labeled as a proxy (not ground truth).
- **nnSVG** is R-only; gate behind availability (ask before adding an R dependency).
- **Naming:** MOSAIC vs MOSAIK (2025 arXiv) collision to resolve before the fair — not yet actioned.

## Invariants held
Year-1 numbers unchanged; new axes are config-driven and backward-compatible (default = Year-1
behavior); every spec seeded; determinism preserved.
