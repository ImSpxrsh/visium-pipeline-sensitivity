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
  median" is the robust claim.) **Rigor follow-up (cycle 6):** the ARI margin had no bootstrap CI
  — added `consensus_advantage_bootstrap_ci`; on an independent rerun the margin's 95% CI
  excludes zero (real effect, not noise), though the exact published point estimate needs
  re-running against the original run data for a matching number.
- **H1 validated against independent pathology, not just a proxy (cycle 9,
  `scripts/analyze_pathology.py`).** Board-certified pathologist labels (SEDR study,
  `data/annotations/`, 100% barcode match to the breast-cancer section) replace the
  consensus-domain boundary proxy used everywhere else. Fragility peaks exactly at the
  pathologist-defined **tumor margin** ("Surrounding tumor"): mean S³ 0.71 (margin) vs 0.87
  (invasive core) vs 0.74 (healthy), permutation p<0.001, survives depth control p<0.001. Removes
  the biggest remaining circularity concern in the tumor H1 result.
- **Fragility is predictable from a single pipeline, before running the multiverse (cycle 10,
  `src/predict_fragility.py`).** Multiverse-free features (local PCA heterogeneity,
  single-pipeline neighbor disagreement, depth, genes detected) predict which spots the full
  30-spec sweep will find fragile, cross-validated AUC 0.81 (genuine out-of-fold, no leakage).
  **Rigor follow-up (cycle 11):** had no CI/null either — added
  `fragility_prediction_significance`; independent rerun gives AUC 0.815, 95% CI [0.788, 0.841],
  permutation p=0.001 — closely matches the published figure, confirmatory not just directional.
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
- **In progress (cycle 10+):** the centerpiece combining cycles 8+9 — does the co-association
  consensus recover the board-certified pathology labels (tumor/margin/healthy) better than any
  single pipeline? Needs the breast per-spec runs regenerated on the converged pipeline first.
- Naming collision (MOSAIC vs MOSAIK) — verified real, recommendation above, needs a human
  decision on disambiguate-vs-rename.
- Exact-matching bootstrap CIs for the consensus ARI margin and (already-done) fragility AUC
  should ideally be re-run against whichever machine holds the original committed run data, since
  `results/runs/` isn't in git and independent re-runs show minor cross-machine numeric drift.

## Open assumptions / decisions
- **Datasets:** environment fetches only programmatic, no-registration data. Human tumor via
  `squidpy.datasets.visium(<10x id>)`. DLPFC/her2st/Janesick need hosted files — wire as URLs allow.
- **H1 annotations:** where pathologist labels are unavailable, use a *consensus-domain boundary*
  proxy for "boundary", clearly labeled as a proxy (not ground truth).
- **nnSVG** is R-only; gate behind availability (ask before adding an R dependency).
- **Naming — verified, not yet actioned (2026-07-09):** the collision is real. MOSAIK
  ("Multi-Origin Spatial Transcriptomics Analysis and Integration Kit," Baptista et al.,
  arXiv:2505.11384, May 2025, JOSS-reviewed, github.com/anthbapt/MOSAIK) is a published,
  actively-reviewed tool in the same subfield (spatial transcriptomics). Different purpose —
  MOSAIK integrates CosMx/Xenium imaging + transcript data into one object (re-segmentation, cell
  typing, domain ID, cell-cell communication); MOSAIC is a spec-curve fragility/consensus
  framework — but the name similarity is a real risk for literature search, citation, and judge
  confusion at the fair. **Recommendation (not actioned — needs a human call, not a unilateral
  rename):** either (a) keep MOSAIC but add an explicit disambiguating footnote/related-work
  citation to MOSAIK in the manuscript's intro or related-work section, or (b) rename before the
  fair if (a) feels insufficient. A full rename touches the paper title, abstract, repo name,
  README, and citation macros — bigger than a single cycle should decide unilaterally.

## Invariants held
Year-1 numbers unchanged; new axes are config-driven and backward-compatible (default = Year-1
behavior); every spec seeded; determinism preserved.
