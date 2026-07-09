# Overnight run log

Running log for the unattended agent cycles (this instance: Adharv's machine). Honest,
no-overclaiming notes — see `CLAUDE.md` invariants and `PROGRESS.md` for the source of truth
on hypothesis status.

## Cycle 1 — 2026-07-08

- **Blocker resolved:** this machine had no Python 3.11/3.12 (Makefile requires it); installed
  Python 3.12 via Homebrew, created `.venv`, installed `requirements.txt`.
- **Sanity check:** `make test` equivalent (`pytest -m "not network"`) — 30/30 passed, environment
  healthy on a fresh clone.
- **Hygiene fix:** `data/visium/` (squidpy's raw download cache, ~47MB) wasn't covered by
  `.gitignore`'s data-exclusion patterns; added it so raw re-downloadable data doesn't get
  committed by accident.
- **Unblocked PROGRESS.md's stated next step:** kicked off
  `PYTHONPATH=. .venv/bin/python scripts/main.py --dataset breast_cancer` (the tumor sweep for
  H1/H2, previously "blocked on tumor sweep" per `PROGRESS.md`). 30 specs, 3798 spots x 36601
  genes. Running in background; each spec's Moran's I permutation step is the long pole, so this
  is expected to take on the order of hours, not minutes. No results yet — nothing to report on
  H1/H2 for tumor tissue this cycle. Will run `scripts/analyze_dataset.py --dataset breast_cancer`
  once the sweep aggregate exists, and report actual numbers (with permutation p-values and
  bootstrap CIs, per the repo's rigor guardrails) in the next entry — not before.
