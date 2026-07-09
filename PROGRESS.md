# PROGRESS — MOSAIC Year 2

Maps work to the ISEF hypotheses (H1/H2/H3) and the backlog in `CLAUDE.md`.

## Status
- **Year 1 (foundation, complete):** mouse-brain Visium (H&E + fluorescence), 30-spec sweep, S³,
  ARI, anatomy validation, variance decomposition (resolution ~80%), confound check (r=-0.09),
  downstream ML (5 pp swing), bioRxiv manuscript. All on `main`, CI green.

## Year 2 backlog
| # | Task | Hypothesis | Status |
|---|---|---|---|
| 1 | Human tumor Visium loads + sweep runs unchanged | H1/H2/H3 | in progress |
| 2 | SVG-method axis (+SpatialDE, +nnSVG) | H2 | not started |
| 3 | QC/filtering axis | — | not started |
| 4 | Stakes-stratified S³ (boundary/interior/healthy) | H1 | not started |
| 5 | Per-dataset Type-II ANOVA (does dominant axis move?) | H2 | not started |
| 6 | Cross-platform S³ (Xenium+Visium) | H3 | not started |
| 7 | Downstream calibration swing (ECE/Brier) | — | not started |
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
