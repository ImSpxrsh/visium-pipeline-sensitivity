# Changelog

## Year 2 (ISEF) — in progress
- **Tumor sweep fixed**: raw 10x breast-cancer matrix broke 6/30 specs (all-zero genes breaking
  scanpy's HVG binning; expm1 overflow on immunoglobulin markers under non-log normalization).
  Fixed via upstream gene QC + a clip-fallback retry in HVG selection; 30/30 specs now succeed,
  Year-1 mouse-brain numbers unchanged.
- **H1 on tumor tissue**: boundary mean S³ 0.727 vs interior 0.873, gap 0.147, p=0.001 (3798
  spots) — replicates and modestly widens the mouse-brain gap (0.13), supporting H1. Depth
  confound (r=0.447) is notably higher than mouse-brain (r=-0.09) and not yet resolved — treat
  as supportive, not conclusive.
- **H2 on tumor tissue**: dominant variance axis stays resolution (72%, vs ~80% mouse-brain) —
  honest null for "the dominant axis shifts with tissue type."
- Add Year-2 planning + tracking docs (`CLAUDE.md`, `ISEF_MOSAIC_YEAR2_PLAN.md`, `PROGRESS.md`).
  Goal: extend MOSAIC to human tumor/disease tissue across platforms and test H1/H2/H3.
- **H1** stakes-stratified S³ (boundary vs interior, permutation null). Mouse-brain baseline:
  boundary 0.71 vs interior 0.85, gap 0.13, p<0.001. Added to manuscript (Fig `13_stakes`).
- **Human tumor loader**: `--dataset breast_cancer` / `glioblastoma` via public 10x Visium.
- `prepare_data` records `total_counts` (exact depth) for the confound check on new datasets.
- **Downstream calibration**: Brier/ECE per spec; Brier ~doubles across pipelines (0.089→0.159).
- Tooling: `scripts/analyze_dataset.py` runs the stats suite on any dataset's sweep.

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
