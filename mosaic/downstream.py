"""
Downstream machine-learning impact of multiverse instability.

A key consequence of pipeline fragility: the representation a model trains on is
itself pipeline-dependent. Here we quantify how much a standard downstream task --
predicting anatomical region from a pipeline's low-dimensional embedding -- varies
across the multiverse. For each specification we cross-validate a classifier on that
spec's ``X_pca`` and report accuracy/F1, so the spread across specs measures how much
the analytical pipeline alone alters model performance.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _brier_ece(y_true, proba, classes, n_bins: int = 10):
    """Multiclass Brier score and expected calibration error from CV probabilities."""
    idx = {c: k for k, c in enumerate(classes)}
    onehot = np.zeros_like(proba)
    for i, yt in enumerate(y_true):
        onehot[i, idx[yt]] = 1.0
    brier = float(np.mean(np.sum((proba - onehot) ** 2, axis=1)))

    conf = proba.max(axis=1)
    pred = classes[proba.argmax(axis=1)]
    correct = (pred == y_true).astype(float)
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (conf > lo) & (conf <= hi)
        if m.any():
            ece += (m.mean()) * abs(correct[m].mean() - conf[m].mean())
    return brier, float(ece)


def downstream_region_classification(
    runs: List[Dict[str, Any]],
    reference_key: str = "cluster",
    n_splits: int = 5,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Cross-validated region classification from each spec's PCA embedding.

    Returns a row per specification with mean CV accuracy, macro-F1, and calibration
    (Brier score and expected calibration error) computed from cross-validated
    predicted probabilities.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_validate
    from sklearn.preprocessing import StandardScaler

    rows = []
    for r in runs:
        a = r["adata"]
        if "X_pca" not in a.obsm or reference_key not in a.obs:
            continue
        x = StandardScaler().fit_transform(np.asarray(a.obsm["X_pca"]))
        # Force a plain numpy string array (pandas 3 arrow-backed arrays break
        # sklearn's fancy indexing during cross-validation).
        y = np.array(a.obs[reference_key].astype(str).tolist())
        clf = LogisticRegression(max_iter=1000)
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        scores = cross_validate(clf, x, y, cv=cv,
                                scoring=("accuracy", "f1_macro"), n_jobs=1)
        proba = cross_val_predict(clf, x, y, cv=cv, method="predict_proba", n_jobs=1)
        classes = np.unique(y)
        brier, ece = _brier_ece(y, proba, classes)
        rows.append({
            "spec_name": r["metadata"].get("spec_name", r.get("name", "")),
            "resolution": r["spec"].get("clustering", {}).get("resolution", np.nan),
            "n_pca": int(x.shape[1]),
            "n_domains": r["metadata"].get("n_clusters_found", np.nan),
            "accuracy": float(scores["test_accuracy"].mean()),
            "accuracy_std": float(scores["test_accuracy"].std()),
            "f1_macro": float(scores["test_f1_macro"].mean()),
            "brier": brier,
            "ece": ece,
        })
    df = pd.DataFrame(rows).sort_values("accuracy", ascending=False).reset_index(drop=True)
    if len(df):
        logger.info(
            "Downstream region classification accuracy across %d specs: "
            "%.3f-%.3f (range %.1f pp)",
            len(df), df["accuracy"].min(), df["accuracy"].max(),
            100 * (df["accuracy"].max() - df["accuracy"].min()),
        )
    return df


def summarize(df: pd.DataFrame) -> Dict[str, Any]:
    """Headline numbers for the downstream experiment."""
    if df.empty:
        return {}
    out = {
        "n_specs": int(len(df)),
        "accuracy_min": float(df["accuracy"].min()),
        "accuracy_max": float(df["accuracy"].max()),
        "accuracy_range_pp": float(100 * (df["accuracy"].max() - df["accuracy"].min())),
        "f1_min": float(df["f1_macro"].min()),
        "f1_max": float(df["f1_macro"].max()),
    }
    if "brier" in df.columns:
        out.update({
            "brier_min": float(df["brier"].min()),
            "brier_max": float(df["brier"].max()),
            "brier_range": float(df["brier"].max() - df["brier"].min()),
            "ece_min": float(df["ece"].min()),
            "ece_max": float(df["ece"].max()),
            "ece_range": float(df["ece"].max() - df["ece"].min()),
        })
    return out


def plot_downstream(df: pd.DataFrame, output_path: Optional[Path] = None) -> None:
    """Accuracy across specs (sorted) with the achievable range highlighted."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [1.3, 1]})

    d = df.sort_values("accuracy").reset_index(drop=True)
    axes[0].errorbar(range(len(d)), d["accuracy"], yerr=d["accuracy_std"],
                     marker="o", color="#377EB8", capsize=2, lw=1)
    axes[0].axhspan(d["accuracy"].min(), d["accuracy"].max(), color="#377EB8", alpha=0.1)
    axes[0].set_xlabel("Specification (sorted by accuracy)")
    axes[0].set_ylabel("Region-classification accuracy (5-fold CV)")
    rng = 100 * (d["accuracy"].max() - d["accuracy"].min())
    axes[0].set_title("A  Accuracy", fontsize=16, fontweight="bold")

    if "ece" in df.columns:
        de = df.sort_values("ece").reset_index(drop=True)
        axes[1].bar(range(len(de)), de["ece"], color="#E6550D")
        axes[1].set_xlabel("Specification (sorted by ECE)")
        axes[1].set_ylabel("Expected calibration error")
        rng_e = de["ece"].max() - de["ece"].min()
        axes[1].set_title("B  Calibration error", fontsize=16, fontweight="bold")
    else:
        axes[1].scatter(df["resolution"], df["accuracy"], c="#984EA3", s=50)
        axes[1].set_xlabel("Clustering resolution")
        axes[1].set_ylabel("Accuracy")
        axes[1].set_title("B  Accuracy vs upstream resolution", fontsize=16, fontweight="bold")

    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info("Saved figure to %s", output_path)
    plt.close(fig)


def main() -> None:
    from mosaic.utils import save_json
    from mosaic.visualization import load_runs

    logging.basicConfig(level=logging.INFO)
    project = Path(__file__).resolve().parent.parent
    runs = load_runs(project / "results")
    df = downstream_region_classification(runs)
    agg = project / "results" / "aggregate"
    df.to_csv(agg / "downstream_classification.csv", index=False)
    save_json(summarize(df), agg / "downstream_summary.json")
    plot_downstream(df, project / "results" / "figures" / "12_downstream_ml.png")


if __name__ == "__main__":
    main()
