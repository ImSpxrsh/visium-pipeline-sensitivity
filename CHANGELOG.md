# Changelog

## Year 2 (ISEF) — in progress
- Add Year-2 planning + tracking docs (`CLAUDE.md`, `ISEF_MOSAIC_YEAR2_PLAN.md`, `PROGRESS.md`).
  Goal: extend MOSAIC to human tumor/disease tissue across platforms and test H1/H2/H3.

## Year 1 (foundation) — complete
- Parametric spatial-transcriptomics pipeline (normalize → HVG → PCA → graph → Leiden/Louvain →
  Moran's I), fully seeded and deterministic.
- Spatial Stability Score (S³): per-spot (Hungarian label alignment), per-feature, headline;
  label-free pairwise ARI; ARI vs reference anatomy.
- 30-spec multiverse sweep with incremental (content-hashed) execution.
- Two mouse-brain datasets (H&E `visium_hne`, fluorescence `visium_fluo`); generalization confirmed.
- Statistics: permutation test, Type-II variance decomposition (resolution ~80%), bootstrap CIs,
  sequencing-depth confound check (r=-0.09), multiverse-size convergence.
- Downstream-ML impact: region-classification accuracy swings ~5 pp across pipelines.
- Real-histology figures (raw H&E, boundary insets, marker genes) + interactive explorer.
- bioRxiv-format manuscript with formal S³ math, complexity analysis, architecture diagram.
- 27 offline tests, CI on Python 3.11/3.12 + ruff.
