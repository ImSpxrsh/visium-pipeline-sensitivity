# CLAUDE.md — MOSAIC (working guidance for the coding agent)

## What this repo is
MOSAIC is a specification-curve / multiverse framework for spatial transcriptomics. It runs a curated set of defensible analysis pipelines over fixed data and quantifies how much biological conclusions swing across analyst choices, via the **Spatial Stability Score (S³)** at spot / gene / headline levels, with cluster labels aligned across pipelines by the Hungarian algorithm.

## Current goal (ISEF Year 2)
Extend from clean mouse-brain Visium to **human tumor and disease tissue across platforms** and test three hypotheses:
- **H1:** fragility (low S³) is worse at tumor margins / disease boundaries than interiors or healthy tissue.
- **H2:** the variance-dominating analytical axis changes with tissue type (mouse brain = clustering resolution ~80%).
- **H3:** the S³ distribution's robust-core/fragile-tail shape conserves Visium→Xenium but absolute fragility rises in tumor.

See `ISEF_MOSAIC_YEAR2_PLAN.md` for the full plan. Treat that as the source of truth for priorities.

## Non-negotiable invariants — never break these
1. **Reproducibility.** Every specification is seeded; each pipeline re-derives its normalized matrix from stored raw counts. No information leaks between specs. `make repro` must always reproduce the primary result. Sweeps stay idempotent and content-hashed.
2. **Determinism check.** After any pipeline change, re-run one fixed spec twice and assert identical spot-level labels.
3. **No leakage / no p-hacking.** Do not tune choices to make a hypothesis land. Every positive claim ships with a null (permutation / label-shuffle) and a bootstrap CI.
4. **Honest nulls.** If a hypothesis fails, report it plainly in results and CHANGELOG. A null is a valid outcome, not a bug to fix.
5. **Config-driven.** New analytical axes and datasets are added via config, not by forking the pipeline. New datasets go behind a `--dataset` flag. Preserve the existing public API and function signatures unless explicitly asked.

## How to work
- **Read before writing.** Inspect the current module layout (data loader, `SpatialPipeline`, `StabilityScorer`, config, viz) and the existing tests before editing. Match existing style.
- **Small, tested increments.** One axis or one dataset per change. Add/extend unit tests with each change; run the test suite before declaring done.
- **Scientific tasks are not just code tasks.** When adding an axis (e.g. nnSVG, SpatialDE, QC thresholds), confirm it's comparable across specs (resolution-aware, seeded) the way Leiden/Louvain were made comparable.
- **Keep figures scripted.** Every figure regenerates from a script; no hand-edited outputs.
- **Track progress.** Maintain `CHANGELOG.md` (what changed, why) and a short `PROGRESS.md` mapping work to H1/H2/H3 and remaining plan items.

## Task backlog (pull from the top; update as you go)
1. Load DLPFC (spatialLIBD) + one human tumor Visium; confirm existing sweep runs unchanged on new data.
2. Add **SVG-method axis** (Moran's I → +nnSVG, +SpatialDE) with a shared significance convention.
3. Add **QC/filtering axis** (min counts, min genes/spot).
4. Implement **stakes-stratified S³**: ingest annotations, label boundary/interior/healthy, compare with permutation null (H1).
5. **Per-dataset Type-II ANOVA** on domain count; report if the dominant axis moves (H2).
6. **Cross-platform S³** on the Janesick Xenium+Visium serial sections; reuse 10x companion registration (H3).
7. Extend downstream classifier to a clinical-flavored label; report accuracy **and** calibration (ECE/Brier) swing.
8. **Fragility-aware reporting demo**: surface one headline (SVG or domain) that flips across defensible specs; show per-call S³ flags it.

## Guardrails on scope
- H1 alone is a complete result. Do not let H2/H3 or extra datasets block H1.
- Do not add heavy dependencies without asking. Prefer the existing stack (NumPy, SciPy, scikit-learn, SCANPY, Squidpy, leidenalg, igraph).
- Do not silently change defaults that affect published Year-1 numbers; if a change would, flag it and gate it behind a flag.

## When unsure
State the assumption in `PROGRESS.md`, make the smallest reversible change, and leave the reasoning in the commit message. Ask before anything that alters reproducibility guarantees or the S³ definition.
