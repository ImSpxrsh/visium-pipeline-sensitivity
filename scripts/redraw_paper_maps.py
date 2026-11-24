"""Redraw the two spatial-map figures in the paper (Fig. 1 montage and the
fragility-prediction map) from the curated runs, without recomputing stats.

Run: PYTHONPATH=. python scripts/redraw_paper_maps.py
"""
from pathlib import Path

import pandas as pd

from mosaic.predict_fragility import plot_predicted_vs_actual, predict_fragility
from mosaic.visualization import figure1_spec_montage, load_runs

PROJECT = Path(__file__).resolve().parent.parent

if __name__ == "__main__":
    runs = load_runs(PROJECT / "results")
    per_spot = pd.read_csv(PROJECT / "results" / "aggregate" / "per_spot_stability.csv")
    figs = PROJECT / "results" / "figures"
    figure1_spec_montage(runs, per_spot, figs / "fig1_spec_montage.png")
    result = predict_fragility(runs, per_spot)
    print("AUC", result.get("auc_combined"))
    plot_predicted_vs_actual(result, figs / "17_fragility_prediction.png")
