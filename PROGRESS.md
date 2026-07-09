# PROGRESS — MOSAIC Year 2

Maps work to the ISEF hypotheses (H1/H2/H3) and the backlog in `CLAUDE.md`.

## Status
- **Year 1 (foundation, complete):** mouse-brain Visium (H&E + fluorescence), 30-spec sweep, S³,
  ARI, anatomy validation, variance decomposition (resolution ~80%), confound check (r=-0.09),
  downstream ML (5 pp swing), bioRxiv manuscript. All on `main`, CI green.

## Year 2 backlog
| # | Task | Hypothesis | Status |
|---|---|---|---|
| 1 | Human tumor Visium loads + sweep runs unchanged | H1/H2/H3 | **fixed 6 spec failures; full sweep re-running** |
| 2 | SVG-method axis (+Geary's C; SpatialDE/nnSVG later) | H2 | **Geary's C done** |
| 3 | QC/filtering axis | — | minimal shared gene QC for raw 10x (not a per-spec axis yet) |
| 4 | Stakes-stratified S³ (boundary/interior/healthy) | H1 | **done (mouse-brain baseline)** |

## Findings so far (Year 2)
- **H1 (mouse-brain baseline):** boundary spots are significantly more fragile than
  interiors — mean S³ 0.72 vs 0.85, gap 0.13, permutation p<0.001. Establishes the
  method; tumor tissue expected to widen the gap.
- **Human tumor loader:** `--dataset breast_cancer` (10x V1_Breast_Cancer_Block_A_Section_1)
  and `glioblastoma` wired via squidpy; no registration needed. Breast sweep (3798 spots,
  36601 genes) running to test H1/H2/H3 on tumor tissue.
- **Downstream calibration:** Brier score nearly doubles (0.089→0.159) across pipelines while
  accuracy moves 5 pp — pipeline choice affects calibration more than accuracy.

## Next (blocked on tumor sweep)
- Run `scripts/analyze_dataset.py --dataset breast_cancer` for tumor H1 (stakes gap) and
  **H2** (does the variance-dominating axis move off resolution on tumor tissue?).
| 5 | Per-dataset Type-II ANOVA (does dominant axis move?) | H2 | not started |
| 6 | Cross-platform S³ (Xenium+Visium) | H3 | not started |
| 7 | Downstream calibration swing (ECE/Brier) | — | **done** |
| 8 | Fragility-aware reporting demo | — | not started |

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
