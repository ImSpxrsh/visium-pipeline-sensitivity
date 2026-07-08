"""
Parametric spatial transcriptomics pipeline.

A single ``spec`` dict defines every analytical choice. Given the same seed and
spec, the pipeline is deterministic. Each run starts from raw counts stored in
``adata.layers["counts"]`` so specs never leak information into one another.

Pipeline stages:
    1. Normalization
    2. Feature selection (highly variable genes)
    3. Dimensionality reduction (PCA)
    4. Expression neighborhood graph
    5. Clustering (Leiden / Louvain)
    6. Spatially variable genes (Moran's I via squidpy)
    7. (Optional) cell-cell communication (placeholder)
"""

import logging
import time
from typing import Any, Dict

import anndata as ad
import numpy as np

logger = logging.getLogger(__name__)

# Map friendly HVG method names to valid scanpy flavors.
# NOTE: the master prompt used "scanpy" as a flavor, which scanpy does not accept;
# the default scanpy dispersion method is exposed under the "cell_ranger" flavor.
_HVG_FLAVOR = {
    "seurat": "seurat",
    "scanpy": "cell_ranger",
    "cell_ranger": "cell_ranger",
}


class SpatialPipeline:
    """Run one specification of the spatial transcriptomics pipeline."""

    def __init__(
        self,
        adata: ad.AnnData,
        spec: Dict[str, Any],
        random_state: int = 42,
        verbose: bool = True,
    ):
        """
        Args:
            adata: AnnData with raw counts (a ``counts`` layer is expected).
            spec: Specification dict (see ``configs/multiverse_spec.yaml``).
            random_state: Seed applied to every stochastic step.
            verbose: Print progress.
        """
        self.adata = adata.copy()
        self.spec = spec
        self.random_state = random_state
        self.verbose = verbose
        self.metadata: Dict[str, Any] = {
            "spec": spec,
            "random_state": random_state,
            "runtime_seconds": None,
        }
        np.random.seed(random_state)

    def run(self) -> Dict[str, Any]:
        """Execute the full pipeline and return a standardized result dict."""
        import scanpy as sc

        sc.settings.verbosity = 1 if self.verbose else 0
        start = time.time()

        try:
            self._normalize()
            self._feature_selection()
            self._pca()
            self._neighborhood_graph()
            self._clustering()
            self._svg_detection()
            if self.spec.get("cell_cell_communication", {}).get("enabled", False):
                self._cell_cell_communication()

            elapsed = time.time() - start
            self.metadata["runtime_seconds"] = elapsed
            if self.verbose:
                logger.info("Pipeline completed in %.1f s", elapsed)

            return {
                "adata": self.adata,
                "n_clusters": int(self.adata.obs["leiden"].nunique()),
                "svg_results": self.adata.uns.get("svg_results"),
                "metadata": self.metadata,
            }
        except Exception as e:  # record and re-raise so the sweep can log it
            logger.error("Pipeline failed: %s", e)
            self.metadata["error"] = str(e)
            raise

    # --- stages ---------------------------------------------------------------

    def _normalize(self) -> None:
        import scanpy as sc

        norm = self.spec.get("normalization", "log1p")
        if self.verbose:
            logger.info("Normalizing with: %s", norm)

        # Always start from raw counts so specs are independent.
        if "counts" in self.adata.layers:
            self.adata.X = self.adata.layers["counts"].copy()

        if norm == "log1p":
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            sc.pp.log1p(self.adata)
        elif norm == "sqrt":
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            self.adata.X = np.sqrt(
                self.adata.X.toarray() if hasattr(self.adata.X, "toarray") else self.adata.X
            )
        elif norm == "total":
            sc.pp.normalize_total(self.adata, target_sum=1e4)
        else:
            raise ValueError(f"Unknown normalization: {norm!r}")

        self.metadata["normalization_used"] = norm

    def _feature_selection(self) -> None:
        import scanpy as sc

        fs = self.spec.get("feature_selection", {})
        method = fs.get("method", "seurat")
        n_genes = fs.get("n_genes", 3000)
        flavor = _HVG_FLAVOR.get(method)
        if flavor is None:
            raise ValueError(f"Unknown HVG method: {method!r}")

        if self.verbose:
            logger.info("Feature selection: %s (flavor=%s), n=%d", method, flavor, n_genes)

        # "seurat"/"cell_ranger" flavors expect log-like data. Guard against the
        # degenerate case where n_top_genes exceeds the available gene count.
        n_top = min(n_genes, self.adata.n_vars)
        sc.pp.highly_variable_genes(self.adata, n_top_genes=n_top, flavor=flavor)

        self.adata = self.adata[:, self.adata.var["highly_variable"]].copy()
        if self.verbose:
            logger.info("Selected %d HVGs", self.adata.n_vars)
        self.metadata["n_hvgs"] = int(self.adata.n_vars)

    def _pca(self) -> None:
        import scanpy as sc

        n_comps = self.spec.get("dimensionality_reduction", {}).get("n_components", 30)
        # PCA components cannot exceed min(n_obs, n_vars) - 1.
        n_comps = int(min(n_comps, self.adata.n_vars - 1, self.adata.n_obs - 1))
        if self.verbose:
            logger.info("PCA with n_components=%d", n_comps)
        sc.tl.pca(self.adata, n_comps=n_comps, random_state=self.random_state)
        self.metadata["n_pca_components"] = n_comps

    def _neighborhood_graph(self) -> None:
        import scanpy as sc

        graph = self.spec.get("neighborhood_graph", {})
        graph_type = graph.get("type", "knn")
        param = graph.get("k_or_radius", 10)

        if self.verbose:
            logger.info("Neighborhood graph: %s (param=%s)", graph_type, param)

        # scanpy's expression graph is k-NN based; it has no radius argument.
        # For a "radius" spec we build a denser k-NN graph as a documented
        # approximation (a true radius graph belongs to the spatial domain, which
        # SVG detection handles separately via sq.gr.spatial_neighbors).
        if graph_type == "knn":
            n_neighbors = int(param)
        elif graph_type == "radius":
            n_neighbors = int(param)
            logger.warning(
                "radius graph approximated as k-NN with k=%d for the expression graph",
                n_neighbors,
            )
        else:
            raise ValueError(f"Unknown graph type: {graph_type!r}")

        sc.pp.neighbors(
            self.adata,
            n_neighbors=n_neighbors,
            use_rep="X_pca",
            random_state=self.random_state,
        )
        self.metadata["graph_type"] = graph_type
        self.metadata["graph_param"] = param

    def _clustering(self) -> None:
        import scanpy as sc

        clust = self.spec.get("clustering", {})
        method = clust.get("method", "leiden")
        resolution = clust.get("resolution", 0.5)
        if self.verbose:
            logger.info("Clustering: %s (resolution=%s)", method, resolution)

        # Both methods write to the "leiden" key so downstream code is uniform.
        if method == "leiden":
            # Pin the leidenalg flavor: it is resolution-aware and deterministic
            # given a seed. (scanpy's default flavor is changing to "igraph",
            # which we avoid to keep results stable across scanpy versions.)
            sc.tl.leiden(
                self.adata,
                resolution=resolution,
                random_state=self.random_state,
                flavor="leidenalg",
                key_added="leiden",
            )
        elif method == "louvain":
            # Drive Louvain through igraph's community_multilevel: it honors the
            # resolution parameter monotonically, unlike scanpy's igraph-flavor
            # louvain, and does not need the unbuildable `louvain` PyPI package.
            self._louvain_igraph(resolution)
        else:
            raise ValueError(f"Unknown clustering method: {method!r}")

        n_clusters = int(self.adata.obs["leiden"].nunique())
        if self.verbose:
            logger.info("Found %d clusters", n_clusters)
        self.metadata["clustering_method"] = method
        self.metadata["clustering_resolution"] = resolution
        self.metadata["n_clusters_found"] = n_clusters

    def _louvain_igraph(self, resolution: float) -> None:
        """Resolution-aware, deterministic Louvain via python-igraph."""
        import random

        import igraph
        import pandas as pd
        from scanpy._utils import get_igraph_from_adjacency

        # Seed igraph's RNG so the (order-dependent) Louvain pass is reproducible.
        # igraph expects a random.Random-style generator (needs .gauss), not a
        # numpy RandomState.
        igraph.set_random_number_generator(random.Random(self.random_state))
        g = get_igraph_from_adjacency(self.adata.obsp["connectivities"], directed=False)
        partition = g.community_multilevel(weights="weight", resolution=resolution)
        labels = np.asarray(partition.membership).astype(str)
        self.adata.obs["leiden"] = pd.Categorical(labels)

    def _svg_detection(self) -> None:
        import squidpy as sq

        method = self.spec.get("spatially_variable_genes", {}).get("method", "morans_i")
        if self.verbose:
            logger.info("SVG detection: %s", method)

        if method != "morans_i":
            raise ValueError(f"Unknown SVG method: {method!r}")

        # spatial_autocorr needs a spatial graph; the master prompt omitted this,
        # which would raise a KeyError on adata.obsp["spatial_connectivities"].
        sq.gr.spatial_neighbors(self.adata)
        sq.gr.spatial_autocorr(self.adata, mode="moran", n_perms=100, n_jobs=1, seed=self.random_state)

        # squidpy writes Moran's I results to adata.uns["moranI"] (not "morans_i").
        svg_df = self.adata.uns["moranI"].sort_values("I", ascending=False)
        self.adata.uns["svg_results"] = svg_df
        if self.verbose:
            logger.info("Computed Moran's I for %d genes", len(svg_df))
        self.metadata["n_svg_genes_tested"] = int(len(svg_df))

    def _cell_cell_communication(self) -> None:
        method = self.spec.get("cell_cell_communication", {}).get(
            "method", "distance_weighted_lr"
        )
        if self.verbose:
            logger.info("CCC detection: %s (placeholder)", method)
        self.metadata["ccc_method"] = method


def run_pipeline(
    adata: ad.AnnData,
    spec: Dict[str, Any],
    random_state: int = 42,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Convenience wrapper: build a SpatialPipeline and run it."""
    return SpatialPipeline(adata, spec, random_state=random_state, verbose=verbose).run()


def main() -> None:
    """Run one spec end-to-end as a smoke test."""
    from src.data_loader import load_dataset, prepare_data

    logging.basicConfig(level=logging.INFO)
    adata, _ = load_dataset("visium_hne", verbose=True)
    adata = prepare_data(adata, verbose=True)

    spec = {
        "normalization": "log1p",
        "feature_selection": {"method": "seurat", "n_genes": 3000},
        "dimensionality_reduction": {"n_components": 30},
        "neighborhood_graph": {"type": "knn", "k_or_radius": 10},
        "clustering": {"method": "leiden", "resolution": 0.5},
        "spatially_variable_genes": {"method": "morans_i"},
    }
    result = run_pipeline(adata, spec, verbose=True)
    logger.info("Found %d clusters", result["n_clusters"])
    logger.info("Pipeline test complete")


if __name__ == "__main__":
    main()
