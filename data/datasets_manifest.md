# Datasets Manifest

MOSAIC is evaluated on independent public spatial-transcriptomics datasets. All are
free and download automatically via `squidpy.datasets`; no application is required.

## Datasets used

| ID | Accession / source | Platform | Spots/cells | Genes | Tissue | Reference annotation | Role |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `visium_hne` | 10x V1_Adult_Mouse_Brain (squidpy) | Visium (H&E) | 2,688 | 18,078 | Mouse brain, coronal | 15 anatomical regions (`obs['cluster']`) | Primary |
| `visium_fluo` | 10x Adult_Mouse_Brain (fluorescence, squidpy) | Visium (fluorescence) | 2,800 | 16,562 | Mouse brain, coronal | 15 anatomical regions (`obs['cluster']`) | Generalization |

Both ship the real H&E / fluorescence histology image and a curated anatomical
annotation, which MOSAIC uses as ground truth for the S³-vs-anatomy validation.

## Justification

- **`visium_hne` (primary):** canonical, well-characterized mouse-brain reference with
  clear anatomical structure (cortex layers, hippocampus, thalamus, striatum, fiber
  tracts), ideal for interpreting where instability concentrates.
- **`visium_fluo` (generalization):** an independently prepared section on a different
  imaging modality (fluorescence vs H&E), testing that the instability pattern is not
  specific to a single slide or staining protocol.

## Candidate additional datasets (future work)

These strengthen external validity but require manual download and/or format wrangling
and are **not** included in the current results (documented here for reproducibility and
to scope follow-up work, not claimed as evidence):

- **10x Genomics FFPE human Visium** (e.g. lymph node, lung adenocarcinoma) —
  https://www.10xgenomics.com/datasets — tests generalization beyond curated demo data.
- **10x Genomics Xenium** (breast/lung, single-cell resolution) — adds a cell
  segmentation axis and single-cell histology; highest engineering cost.
- **GEO Visium cohort** (e.g. colon or prostate adenocarcinoma) — an external, non-vendor
  cohort to guard against demo-data cherry-picking.

Notes: `squidpy.datasets.visium_hne_adata()` is the **mouse-brain** coronal section (its
gene symbols are mouse nomenclature), not a breast-cancer slide.
