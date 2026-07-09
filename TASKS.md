# MOSAIC Task Tracker

## Sparsh (lead) + Claude Code

### Phase 1 — Foundation
- [x] Repo structure + dependencies
- [x] `data_loader.py`: load + validate Visium H&E demo
- [x] `utils.py`: config, logging, IO, grid expansion
- [ ] Smoke-test loader against real download (`python -m mosaic.data_loader`)

### Phase 2 — Core pipeline
- [x] `pipeline.py`: parametric pipeline (norm -> HVG -> PCA -> graph -> cluster -> SVG)
- [x] Determinism test (`scripts/test_reproducibility.py`)
- [ ] Eyeball 3–5 specs for sensible clusters/SVGs

### Phase 3 — Sweep + stability
- [x] `configs/multiverse_spec.yaml`: curated specs
- [x] `main.py`: sweep orchestration + failure manifest
- [x] `stability.py`: per-spot / per-feature S³, headline instability
- [ ] Run full sweep

### Phase 4 — Visualization
- [x] `visualization.py`: stability map, spec curve, instability barplot
- [x] `generate_figures.py`
- [ ] Pick 5–8 paper figures

### Phase 5 — Writing
- [ ] Draft manuscript

## Adharv (intern) — Specification Cartographer
- [ ] Literature review (12–15 papers), extract analytical choices
- [ ] Frequency table of choices
- [ ] Curate 30–40 defensible specs into `configs/multiverse_spec.yaml`
- [ ] One-page "Specification Space Definition" methods write-up

## Alex (intern) — Validation & Ground Truth
- [ ] Cross-check the reference anatomical annotation for the mouse-brain slide
- [ ] Build `data/reference/annotation_reference.csv`
- [ ] Hand-check S³ for 5 specs in a spreadsheet
- [ ] Define operational thresholds for "flipping"/instability bands
