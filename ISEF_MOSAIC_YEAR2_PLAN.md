# MOSAIC — Year 2 (ISEF)
## Where the Stakes Are Highest: analytical fragility on human tumor & disease tissue across platforms

**One-line thesis.** Analytical fragility in spatial transcriptomics is *worst exactly where the
biological stakes are highest* — at tumor margins and disease boundaries — and the single most
consequential analytical choice *changes with the tissue*. Year 1 established the framework on clean
mouse brain; Year 2 stress-tests it where conclusions actually carry clinical weight.

This is the "this year" project (ISEF judges evaluate only the current year's work; file **Form 7**
and treat mouse-brain MOSAIC as the Year-1 foundation).

---

### Falsifiable hypotheses (each is a result either way)

- **H1 — Fragility tracks stakes.** Mean per-spot S³ is significantly *lower* (more fragile) at
  annotated tumor margins / tumor–immune interfaces and disease boundaries than in matched tissue
  interiors, and lower than in healthy tissue. *Refuting it is still a finding.*
- **H2 — The dominant axis shifts.** In mouse brain, clustering resolution explains ~80% of
  domain-count variance. In tumor tissue the variance-dominating axis is **different** (candidate:
  SVG method, normalization, or HVG selection). Test via per-dataset Type-II ANOVA.
- **H3 — Structure conserves, magnitude worsens.** The robust-core / fragile-tail *shape* of the S³
  distribution is conserved from Visium → Xenium, but *absolute* fragility is higher in tumor.

---

### Datasets (all public, most pre-annotated, all free)

| Role | Dataset | Why |
|---|---|---|
| Healthy human bridge | **DLPFC** (spatialLIBD, 12 sections, cortical-layer labels) | annotated *human* tissue; gold-standard benchmark |
| Primary tumor | **10x Visium human breast cancer** (IDC, fresh-frozen) | high-stakes tissue; community annotations |
| Second platform + margins | **Janesick 2023 Xenium + Visium serial breast cancer** | same tissue, two platforms → cross-platform S³ |
| Cross-patient tumor | **her2st** (HER2+ breast tumor ST, multiple patients) | tests H1 across patients |
| Optional stretch | Pathologist-annotated lung/kidney (Zenodo) or a Xenium cancer | breadth for generalization |

**Environment note.** This coding environment can fetch datasets only programmatically without
registration. 10x public Visium tumor datasets are reachable via `squidpy.datasets.visium(<id>)`.
DLPFC/her2st/Janesick require specific hosted files; wire them behind `--dataset` as URLs become
available, and never fabricate annotations.

---

### New analytical axes (fold in while you're there)

- **SVG method** as a first-class axis: Moran's I → add SpatialDE (Python) and nnSVG (R; gate behind
  availability). Closes the single-test limitation.
- **QC / filtering thresholds** (min counts, min genes/spot) as an axis — matters more on messy
  tumor data than on demo brain.
- Keep all existing axes (normalization, HVG method/n, PCA, graph k, algorithm, resolution).

---

### New analyses

1. **Stakes-stratified S³.** Label spots as *boundary* (tumor margin / interface / disease edge) vs
   *interior* vs *healthy*; compare mean S³ across strata with the region-size-preserving permutation
   null. → **H1**. (Where pathologist annotations are unavailable, a data-driven consensus-domain
   boundary proxy is used and labeled as such — never presented as pathologist ground truth.)
2. **Per-tissue variance decomposition.** Re-run Type-II ANOVA separately per dataset; report whether
   the dominant axis moves. → **H2**.
3. **Cross-platform S³.** Register Visium/Xenium serial sections; compare S³ distributions and, where
   registered, per-location agreement. → **H3**.
4. **Higher-stakes downstream swing.** Extend the classifier to a clinical label; report accuracy
   **and calibration** (ECE / Brier) swing across the multiverse.
5. **Fragility-aware reporting demo.** Find one headline that flips across defensible pipelines on
   tumor data and show how per-call S³ flags it.

---

### Rigor guardrails (keep Year-1 standards)
Every spec seeded; re-derive from raw counts (no leakage); determinism verified. Positive claims
paired with a null and bootstrap CIs. Depth-confound check per dataset. S³ convergence check.
Honest nulls are contributions.

### ISEF deliverables
Poster spine: one claim (fragility peaks where stakes peak) → per-spot S³ on tumor H&E → variance-shift
figure → cross-platform conservation figure. Category: **Computational Biology & Bioinformatics**
(Systems Software fallback). Resolve the **MOSAIC vs MOSAIK** naming collision before the fair.

### Priority
H1 alone is a complete, winnable project. H2/H3 are upside, not prerequisites.
