.PHONY: help venv install test test-all lint data sweep figures stats downstream \
        tissue fluo cross-dataset numbers analyses repro paper clean

PY := .venv/bin/python
PIP := .venv/bin/pip

help:
	@echo "MOSAIC targets:"
	@echo "  make venv      Create the Python 3.12 virtualenv"
	@echo "  make install   Install dependencies into the venv"
	@echo "  make test      Run offline tests (no network)"
	@echo "  make test-all  Run the full test suite (downloads data)"
	@echo "  make lint      Run ruff"
	@echo "  make data      Download + cache the Visium demo"
	@echo "  make sweep     Run the full multiverse sweep"
	@echo "  make figures   Core publication figures from sweep outputs"
	@echo "  make analyses  figures + stats + downstream + tissue + numbers"
	@echo "  make fluo      Run the sweep on the second dataset (generalization)"
	@echo "  make repro     data -> sweep -> analyses (end to end)"
	@echo "  make paper     Build the bioRxiv PDF (needs LaTeX / tectonic)"

venv:
	python3.12 -m venv .venv
	$(PIP) install --upgrade pip

install: venv
	$(PIP) install -r requirements.txt

test:
	$(PY) -m pytest -m "not network" -q

test-all:
	$(PY) -m pytest -q

lint:
	.venv/bin/ruff check src scripts tests

data:
	PYTHONPATH=. $(PY) -m src.data_loader

sweep:
	PYTHONPATH=. $(PY) scripts/main.py

figures:
	PYTHONPATH=. $(PY) scripts/generate_figures.py

stats:
	PYTHONPATH=. $(PY) -m src.statistics

downstream:
	PYTHONPATH=. $(PY) -m src.downstream

tissue:
	PYTHONPATH=. $(PY) scripts/generate_tissue.py

fluo:
	PYTHONPATH=. $(PY) scripts/main.py --dataset visium_fluo

cross-dataset:
	PYTHONPATH=. $(PY) scripts/generate_cross_dataset.py

numbers:
	PYTHONPATH=. $(PY) scripts/report_numbers.py

# All analysis figures + populated paper macros from an existing sweep.
analyses: figures stats downstream tissue numbers

repro: data sweep analyses

paper:
	cd paper && (latexmk -pdf main.tex || tectonic main.tex)

clean:
	rm -rf results/runs/* results/logs/* __pycache__ src/__pycache__ tests/__pycache__ .pytest_cache
