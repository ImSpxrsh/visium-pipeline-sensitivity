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

## Independent annotations (ground truth)

- **Breast-cancer pathologist annotations** (`data/annotations/brca1_metadata.tsv`): per-spot
  board-certified pathologist labels for `V1_Breast_Cancer_Block_A_Section_1` (Invasive /
  Tumor in-situ / Surrounding-tumor margin / Healthy), from the SEDR study
  (`github.com/JinmiaoChenLab/SEDR_analyses`, `data/BRCA1/metadata.tsv`). 3,798 spots, 100%
  barcode match. Used as *independent* ground truth for H1 (fragility peaks at the tumor
  margin) — not a computational annotation, so it removes the circularity concern.

## Controlled-access datasets to pursue with a mentor (forms)

- **HTAN** (Human Tumor Atlas Network): Xenium/CosMx/Visium HD + pathology; controlled Level 1-2
  via dbGaP `phs002371` (open Level 3-4 via Synapse). Needs an institutional Signing Official + PI.
- **EGA** tumor Visium: HER2+ breast `EGAD00001011061`, melanoma `EGAD00001010921`, prostate
  `EGAD00001007921` — each via its Data Access Committee form.

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
