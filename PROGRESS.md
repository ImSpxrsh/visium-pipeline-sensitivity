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
- **H1 (mouse-brain baseline):** boundary spots are significantly more fragile than
  interiors — mean S³ 0.72 vs 0.85, gap 0.13, permutation p<0.001. Establishes the
  method; tumor tissue expected to widen the gap.
- **H1 (breast tumor, real result — 2026-07-09):** boundary mean S³ = 0.727, interior mean
  S³ = 0.873, gap = 0.147, permutation p = 0.001 (3798 spots, consensus-domain boundary proxy —
  no pathologist annotation for this dataset, so labeled as a proxy per the Open assumptions
  below). The gap is directionally consistent with, and modestly larger than, the mouse-brain
  baseline (0.147 vs 0.13) — supports H1 on a second, higher-stakes tissue. **Caveat, not yet
  resolved:** the depth-confound check here is r=0.447 (S³ vs total UMI counts per spot), far
  higher than mouse-brain's r=-0.09. Some of the boundary/interior gap could be
  driven by depth differences between margin and interior spots rather than analytical
  fragility alone. Treat the tumor H1 result as supportive, not conclusive, until this is
  checked (e.g. depth-matched resampling or partialling out total_counts).
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
- **Resolve the depth confound** on the tumor H1 result before treating it as a clean
  replication (see caveat above) — this is the next priority, not a new dataset or axis.
- H3 (cross-platform Xenium+Visium S³) and the SVG-method axis expansion (SpatialDE/nnSVG)
  remain not started.

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
