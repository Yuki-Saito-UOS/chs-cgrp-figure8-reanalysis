from __future__ import annotations

import numpy as np
import scanpy as sc
from anndata import AnnData


def run_qc_and_global_umap(
    adata_raw: AnnData,
    *,
    min_genes: int = 200,
    max_genes: int = 6000,
    max_pct_mt: float = 20.0,
    target_sum: float = 1e4,
    n_top_genes: int = 3000,
    hvg_batch_key: str | None = "sample",
    n_neighbors: int = 15,
    n_pcs: int = 30,
    pca_chunk_size: int = 1000,
    min_dist: float = 0.3,
    leiden_resolution: float = 0.5,
) -> AnnData:
    """QC、正規化、HVG、PCA、UMAP、Leidenを実行する。"""
    # 大規模データでは入力を複製するとメモリ使用量が倍増する。呼び出し元が
    # raw objectを後続で使わない全細胞解析では、入力をin-placeで前処理する。
    adata = adata_raw

    adata.var_names = adata.var_names.astype(str).str.upper()
    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.match(r"^HB[ABDGEMQZ]")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    qc_mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] < max_pct_mt)
    )

    # 全細胞データでの一時的な行列二重保持を避ける。
    adata._inplace_subset_obs(qc_mask.to_numpy())

    if adata.n_obs == 0:
        raise RuntimeError("QC後に細胞が残りませんでした。")

    print("QC後:", adata.shape)

    sc.pp.normalize_total(
        adata,
        target_sum=target_sum,
    )
    sc.pp.log1p(adata)

    # rawを複製すると全遺伝子行列がもう一つ必要になるため、全細胞の
    # exploratory UMAPでは正規化済みXをDotPlot/feature plotにそのまま使う。

    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor="cell_ranger",
            batch_key=hvg_batch_key,
        )
    except ValueError as error:
        # 少数の陽性細胞だけを再UMAPする場合、平均発現のbin境界が重複して
        # cell_ranger HVGが失敗する。分散上位遺伝子へ安全にフォールバックする。
        print(f"HVG binning failed; using variance fallback: {error}")
        matrix = adata.X
        if hasattr(matrix, "multiply"):
            mean = np.asarray(matrix.mean(axis=0)).ravel()
            mean_sq = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
        else:
            dense = np.asarray(matrix)
            mean = dense.mean(axis=0)
            mean_sq = np.square(dense).mean(axis=0)
        variance = np.maximum(mean_sq - np.square(mean), 0.0)
        n_select = min(int(n_top_genes), adata.n_vars)
        top = np.argpartition(variance, -n_select)[-n_select:]
        adata.var["highly_variable"] = False
        adata.var.iloc[top, adata.var.columns.get_loc("highly_variable")] = True

    if (
        "highly_variable" not in adata.var.columns
        or int(adata.var["highly_variable"].sum()) == 0
    ):
        raise RuntimeError("高変動遺伝子が選択されませんでした。")

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(
        adata_hvg,
        max_value=10,
        # 全細胞データを密行列化せず、外付けSSD上の大規模10x入力でも
        # メモリに収まるようにする。探索UMAPでは疎行列を維持する。
        zero_center=False,
    )

    usable_pcs = min(
        n_pcs,
        adata_hvg.n_obs - 1,
        adata_hvg.n_vars - 1,
    )

    sc.tl.pca(
        adata_hvg,
        n_comps=usable_pcs,
        zero_center=False,
        # IncrementalPCAにより、全細胞を含めたまま大規模行列の一括展開を避ける。
        chunked=True,
        chunk_size=pca_chunk_size,
    )

    sc.pp.neighbors(
        adata_hvg,
        n_neighbors=n_neighbors,
        n_pcs=usable_pcs,
    )

    sc.tl.umap(
        adata_hvg,
        min_dist=min_dist,
    )

    sc.tl.leiden(
        adata_hvg,
        resolution=leiden_resolution,
        key_added="leiden",
    )

    adata.obs["leiden"] = (
        adata_hvg.obs["leiden"]
        .astype(str)
        .values
    )

    adata.obsm["X_umap"] = (
        adata_hvg.obsm["X_umap"].copy()
    )

    adata.obsm["X_pca"] = (
        adata_hvg.obsm["X_pca"].copy()
    )

    adata.uns["neighbors"] = (
        adata_hvg.uns["neighbors"].copy()
    )

    adata.obsp["connectivities"] = (
        adata_hvg.obsp["connectivities"].copy()
    )

    adata.obsp["distances"] = (
        adata_hvg.obsp["distances"].copy()
    )

    return adata


def prepare_subset_umap(
    adata_input: AnnData,
    *,
    n_top_genes: int = 2000,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    min_dist: float = 0.3,
    leiden_resolution: float = 0.5,
    leiden_key: str = "leiden_subset",
) -> AnnData:
    """subsetにHVG、PCA、近傍グラフ、UMAP、Leidenを実行する。"""
    if adata_input.n_obs < 3:
        raise ValueError("UMAPには少なくとも3細胞必要です。")

    adata = adata_input.copy()

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=min(n_top_genes, adata.n_vars),
        flavor="seurat",
    )

    if (
        "highly_variable" in adata.var.columns
        and int(adata.var["highly_variable"].sum()) > 0
    ):
        adata = adata[:, adata.var["highly_variable"]].copy()
    else:
        print("HVG selection failed: all genes will be used.")

    sc.pp.scale(
        adata,
        max_value=10,
    )

    sc.tl.pca(
        adata,
        svd_solver="arpack",
    )

    usable_pcs = min(
        n_pcs,
        adata.obsm["X_pca"].shape[1],
    )

    usable_neighbors = min(
        n_neighbors,
        adata.n_obs - 1,
    )

    sc.pp.neighbors(
        adata,
        n_neighbors=usable_neighbors,
        n_pcs=usable_pcs,
    )

    sc.tl.umap(
        adata,
        min_dist=min_dist,
    )

    sc.tl.leiden(
        adata,
        resolution=leiden_resolution,
        key_added=leiden_key,
    )

    return adata
