from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from .preprocessing import prepare_subset_umap


def get_expr(
    adata: AnnData,
    gene: str,
) -> pd.Series:
    """指定遺伝子の発現量を1次元のSeriesで返す。"""
    if gene not in adata.var_names:
        return pd.Series(
            np.nan,
            index=adata.obs_names,
            dtype=float,
        )

    expression = adata[:, gene].X

    if hasattr(expression, "toarray"):
        expression = expression.toarray()

    return pd.Series(
        np.asarray(expression).ravel(),
        index=adata.obs_names,
        dtype=float,
    )


def show_zero_immune_samples(
    summary: pd.DataFrame,
    *,
    sample_col: str = "sample",
    group_col: str = "group",
    immune_col: str = "n_immune",
) -> pd.DataFrame:
    """免疫細胞数が0のサンプルを表示して返す。"""
    required_columns = {
        sample_col,
        group_col,
        immune_col,
    }

    missing_columns = required_columns - set(summary.columns)

    if missing_columns:
        raise ValueError(
            f"必須列がありません: {sorted(missing_columns)}"
        )

    zero_immune = summary.loc[
        summary[immune_col].fillna(0) == 0,
        [sample_col, group_col, immune_col],
    ].copy()

    if zero_immune.empty:
        print("免疫細胞数が0のサンプルはありません。")
    else:
        print("免疫細胞数が0のサンプル:")
        print(zero_immune.to_string(index=False))

    return zero_immune


def build_cdc2_subsets(
    adata_input: AnnData,
    *,
    immune_gene: str = "PTPRC",
    dc_genes: tuple[str, ...] = (
        "CLEC10A",
        "CD1C",
        "ITGAX",
        "FCER1A",
    ),
    main_gate_genes: tuple[str, ...] = (
        "PTPRC",
        "CLEC10A",
        "CD1C",
    ),
    n_top_genes: int = 2000,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    min_dist: float = 0.3,
    leiden_resolution: float = 0.5,
) -> tuple[AnnData, AnnData, AnnData]:
    """immune、DC-enriched、MAIN_GATE cDC2 subsetを作成する。"""
    adata_all = adata_input.copy()

    adata_all.var_names = (
        adata_all.var_names
        .astype(str)
        .str.upper()
    )
    adata_all.var_names_make_unique()

    if "group" in adata_all.obs.columns:
        adata_all.obs["group"] = (
            adata_all.obs["group"].astype(str)
        )

    # 1. Immune subset
    immune_expression = get_expr(
        adata_all,
        immune_gene,
    )

    immune_mask = (
        immune_expression.notna()
        & (immune_expression > 0)
    )

    adata_immune = adata_all[
        immune_mask.values
    ].copy()

    if adata_immune.n_obs == 0:
        raise RuntimeError(
            f"{immune_gene}陽性細胞がありません。"
        )

    print("immune subset:", adata_immune.shape)

    immune_embedding = prepare_subset_umap(
        adata_immune,
        n_top_genes=n_top_genes,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        min_dist=min_dist,
        leiden_resolution=leiden_resolution,
        leiden_key="leiden_immune",
    )

    adata_immune.obsm["X_umap_immune"] = (
        immune_embedding.obsm["X_umap"].copy()
    )

    adata_immune.obs["leiden_immune"] = (
        immune_embedding.obs["leiden_immune"]
        .astype(str)
        .values
    )

    # 2. DC-enriched subset
    dc_mask = pd.Series(
        False,
        index=adata_immune.obs_names,
    )

    available_dc_genes: list[str] = []

    for gene in dc_genes:
        if gene not in adata_immune.var_names:
            continue

        available_dc_genes.append(gene)
        dc_mask |= get_expr(adata_immune, gene) > 0

    if not available_dc_genes:
        raise RuntimeError(
            "DC抽出用遺伝子が1つも存在しません。"
        )

    adata_dc = adata_immune[
        dc_mask.values
    ].copy()

    if adata_dc.n_obs == 0:
        raise RuntimeError(
            "DC-enriched細胞がありません。"
        )

    print("DC-enriched subset:", adata_dc.shape)

    dc_embedding = prepare_subset_umap(
        adata_dc,
        n_top_genes=n_top_genes,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        min_dist=min_dist,
        leiden_resolution=leiden_resolution,
        leiden_key="leiden_dc",
    )

    adata_dc.obsm["X_umap_dc"] = (
        dc_embedding.obsm["X_umap"].copy()
    )

    adata_dc.obs["leiden_dc"] = (
        dc_embedding.obs["leiden_dc"]
        .astype(str)
        .values
    )

    # 3. MAIN_GATE
    missing_gate_genes = [
        gene
        for gene in main_gate_genes
        if gene not in adata_dc.var_names
    ]

    if missing_gate_genes:
        raise ValueError(
            "MAIN_GATEに必要な遺伝子がありません: "
            f"{missing_gate_genes}"
        )

    main_gate_mask = pd.Series(
        True,
        index=adata_dc.obs_names,
    )

    for gene in main_gate_genes:
        main_gate_mask &= get_expr(adata_dc, gene) > 0

    adata_dc.obs["MAIN_GATE"] = (
        main_gate_mask.values
    )

    print(
        "MAIN_GATE:",
        int(main_gate_mask.sum()),
    )

    adata_cdc2 = adata_dc[
        main_gate_mask.values
    ].copy()

    if adata_cdc2.n_obs == 0:
        raise RuntimeError(
            "MAIN_GATEに該当する細胞がありません。"
        )

    print("cDC2-like subset:", adata_cdc2.shape)

    cdc2_embedding = prepare_subset_umap(
        adata_cdc2,
        n_top_genes=n_top_genes,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        min_dist=min_dist,
        leiden_resolution=leiden_resolution,
        leiden_key="leiden_cdc2",
    )

    adata_cdc2.obsm["X_umap_cdc2"] = (
        cdc2_embedding.obsm["X_umap"].copy()
    )

    adata_cdc2.obs["leiden_cdc2"] = (
        cdc2_embedding.obs["leiden_cdc2"]
        .astype(str)
        .values
    )

    return (
        adata_immune,
        adata_dc,
        adata_cdc2,
    )