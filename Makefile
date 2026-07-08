.PHONY: help venv install test test-all lint data sweep figures repro paper clean

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
	@echo "  make figures   Generate publication figures from sweep outputs"
	@echo "  make repro     data -> sweep -> figures (end to end)"
	@echo "  make paper     Build the NeurIPS PDF (needs LaTeX)"

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

repro: data sweep figures

paper:
	cd paper && latexmk -pdf main.tex

clean:
	rm -rf results/runs/* results/logs/* __pycache__ src/__pycache__ tests/__pycache__ .pytest_cache
