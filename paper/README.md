# Paper: The Garden of Forking Tissues

NeurIPS-formatted manuscript for MOSAIC.

## Build

```bash
# 1. Generate quantitative macros + figures from a completed sweep
PYTHONPATH=. ../.venv/bin/python ../scripts/report_numbers.py
PYTHONPATH=. ../.venv/bin/python ../scripts/generate_figures.py

# 2. Compile
latexmk -pdf main.tex        # or: make -C .. paper
```

## Files
- `main.tex` — the manuscript (uses the NeurIPS style, `[preprint]` option).
- `numbers.tex` — auto-generated result macros (`scripts/report_numbers.py`).
- `references.bib` — bibliography.
- `neurips_2024.sty` — NeurIPS style file.

Figures are read from `../results/figures/`. Remove the `[preprint]` option in
`main.tex` for the anonymized submission version.
