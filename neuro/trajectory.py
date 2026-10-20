from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from .paths import get_project_dir
from .subsets import get_expr

def compute_cdc2_pseudotime(
    adata_cdc2: AnnData,
    *,
    selected_clusters: Sequence[str] | None = None,
    cluster_key: str = "leiden_cdc2",
    root_positive_gene: str = "CCR2",
    root_negative_genes: Sequence[str] = (
        "CCR7",
        "LAMP3",
    ),
    root_negative_weights: Sequence[float] = (
        3.0,
        4.0,
    ),
    n_top_genes: int = 2000,
    n_neighbors: int = 12,
    n_pcs: int = 30,
    min_dist: float = 0.3,
    pseudotime_key: str = "pt",
    umap_key: str = "X_umap_traj",
) -> AnnData:
    """cDC2 subsetについてDPT pseudotimeを計算する。

    rootは次のスコアが最大の細胞とする。

    root score =
        root_positive_gene
        - weighted root_negative_genes

    既定値では、
        CCR2 - 3 × CCR7 - 4 × LAMP3

    Parameters
    ----------
    adata_cdc2
        cDC2 subsetのAnnData。
    selected_clusters
        使用するcluster。Noneなら全clusterを使用。
    cluster_key
        cluster情報が保存されているobs列。
    root_positive_gene
        root側で高発現してほしい遺伝子。
    root_negative_genes
        root側で低発現してほしい遺伝子。
    root_negative_weights
        root_negative_genesに対応する重み。
    n_top_genes
        高変動遺伝子数。
    n_neighbors
        近傍グラフの近傍数。
    n_pcs
        近傍グラフに用いるPC数。
    min_dist
        UMAPのmin_dist。
    pseudotime_key
        pseudotimeの保存先。
    umap_key
        UMAP座標の保存先。

    Returns
    -------
    anndata.AnnData
        pseudotimeとUMAP座標を含むAnnData。
    """
    if cluster_key not in adata_cdc2.obs.columns:
        raise KeyError(
            f"adata_cdc2.obsに'{cluster_key}'がありません。"
        )

    if len(root_negative_genes) != len(root_negative_weights):
        raise ValueError(
            "root_negative_genesと"
            "root_negative_weightsの長さが一致しません。"
        )

    adata = adata_cdc2.copy()

    adata.var_names = (
        adata.var_names
        .astype(str)
        .str.upper()
    )
    adata.var_names_make_unique()

    adata.obs[cluster_key] = (
        adata.obs[cluster_key].astype(str)
    )

    if selected_clusters is not None:
        selected_clusters = tuple(
            str(cluster)
            for cluster in selected_clusters
        )

        mask = adata.obs[cluster_key].isin(
            selected_clusters
        )

        adata = adata[mask].copy()

    if adata.n_obs < 3:
        raise RuntimeError(
            "pseudotime計算に十分な細胞がありません。"
        )

    print("trajectory subset:", adata.shape)

    # root geneの存在確認
    required_root_genes = [
        root_positive_gene,
        *root_negative_genes,
    ]

    missing_root_genes = [
        gene
        for gene in required_root_genes
        if gene not in adata.var_names
    ]

    if missing_root_genes:
        raise ValueError(
            "root決定に必要な遺伝子がありません: "
            f"{missing_root_genes}"
        )

    # HVG
    adata_hvg = adata.copy()

    sc.pp.highly_variable_genes(
        adata_hvg,
        n_top_genes=min(
            n_top_genes,
            adata_hvg.n_vars,
        ),
        flavor="seurat",
    )

    if (
        "highly_variable" in adata_hvg.var.columns
        and int(
            adata_hvg.var["highly_variable"].sum()
        ) > 0
    ):
        adata_hvg = adata_hvg[
            :,
            adata_hvg.var["highly_variable"],
        ].copy()
    else:
        print(
            "HVG selection failed: "
            "all genes will be used."
        )

    # PCA / neighbors / UMAP
    sc.pp.scale(
        adata_hvg,
        max_value=10,
    )

    sc.tl.pca(
        adata_hvg,
        svd_solver="arpack",
    )

    usable_pcs = min(
        n_pcs,
        adata_hvg.obsm["X_pca"].shape[1],
    )

    usable_neighbors = min(
        n_neighbors,
        adata_hvg.n_obs - 1,
    )

    sc.pp.neighbors(
        adata_hvg,
        n_neighbors=usable_neighbors,
        n_pcs=usable_pcs,
    )

    sc.tl.umap(
        adata_hvg,
        min_dist=min_dist,
    )

    # root score
    root_score = get_expr(
        adata,
        root_positive_gene,
    ).copy()

    for gene, weight in zip(
        root_negative_genes,
        root_negative_weights,
        strict=True,
    ):
        root_score -= (
            weight
            * get_expr(
                adata,
                gene,
            )
        )

    root_score = root_score.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    if root_score.isna().all():
        raise RuntimeError(
            "root scoreを計算できませんでした。"
        )

    root_obs_name = root_score.idxmax()
    root_index = int(
        adata.obs_names.get_loc(root_obs_name)
    )

    adata_hvg.uns["iroot"] = root_index

    sc.tl.diffmap(adata_hvg)
    sc.tl.dpt(adata_hvg)

    adata.obs[pseudotime_key] = (
        adata_hvg.obs["dpt_pseudotime"]
        .astype(float)
        .values
    )

    adata.obsm[umap_key] = (
        adata_hvg.obsm["X_umap"].copy()
    )

    adata.obs["root_score"] = (
        root_score.reindex(
            adata.obs_names
        ).values
    )

    adata.obs["is_root"] = False
    adata.obs.loc[
        root_obs_name,
        "is_root",
    ] = True

    print("root cell:", root_obs_name)
    print(
        "root cluster:",
        adata.obs.loc[
            root_obs_name,
            cluster_key,
        ],
    )
    print(
        f"{root_positive_gene}:",
        float(
            get_expr(
                adata,
                root_positive_gene,
            ).loc[root_obs_name]
        ),
    )

    for gene in root_negative_genes:
        print(
            f"{gene}:",
            float(
                get_expr(
                    adata,
                    gene,
                ).loc[root_obs_name]
            ),
        )

    print(
        "pseudotime range:",
        (
            float(adata.obs[pseudotime_key].min()),
            float(adata.obs[pseudotime_key].max()),
        ),
    )

    return adata

def plot_pseudotime_umap(
    adata: AnnData,
    *,
    basis: str = "X_umap_traj",
    pseudotime_key: str = "pt",
    cluster_key: str = "leiden_cdc2",
    genes: Sequence[str] = (
        "CCR2",
        "RAMP1",
        "LAMP3",
        "CCR7",
    ),
) -> None:
    """Pseudotime、cluster、遺伝子発現をUMAP上に表示する。"""
    colors: list[str] = []

    if "group" in adata.obs.columns:
        colors.append("group")

    if cluster_key in adata.obs.columns:
        colors.append(cluster_key)

    if pseudotime_key in adata.obs.columns:
        colors.append(pseudotime_key)

    colors.extend(
        gene
        for gene in genes
        if gene in adata.var_names
    )

    sc.pl.embedding(
        adata,
        basis=basis,
        color=colors,
        ncols=4,
        frameon=False,
        vmax="p99",
    )

def plot_pseudotime_trajectory(
    adata: AnnData,
    *,
    pseudotime_key: str = "pt",
    genes: Sequence[str] = (
        "CCR2",
        "RAMP1",
        "LAMP3",
        "CCR7",
    ),
    n_bins: int = 100,
    smooth_window: int = 10,
    xlim: tuple[float, float] = (0.0, 0.16),
    output_folder: str = "FIG/Fig_E",
    output_name: str = "Pseudotime_trajectory_early",
    normalize_pseudotime: bool = False,
) -> tuple[pd.DataFrame, Path]:
    """Pseudotimeに沿った平均発現量と陽性細胞割合を描画する。

    Parameters
    ----------
    adata
        Pseudotimeを含むAnnData。
    pseudotime_key
        adata.obsに保存されたpseudotime列名。
    genes
        描画する遺伝子。
    n_bins
        Pseudotimeを分割するbin数。
    smooth_window
        移動平均の窓幅。
    xlim
        表示するPseudotime範囲。
    output_folder
        プロジェクトルートから見た出力フォルダ。
    output_name
        拡張子を除いた出力ファイル名。
    normalize_pseudotime
        Pseudotimeを0–1に正規化するか。

    Returns
    -------
    tuple
        解析用DataFrame、PDFパス。
    """
    if pseudotime_key not in adata.obs.columns:
        raise KeyError(
            f"adata.obsに'{pseudotime_key}'がありません。"
        )

    if n_bins < 2:
        raise ValueError("n_binsは2以上にしてください。")

    if smooth_window < 1:
        raise ValueError(
            "smooth_windowは1以上にしてください。"
        )

    available_genes = [
        gene
        for gene in genes
        if gene in adata.var_names
    ]

    if not available_genes:
        raise ValueError(
            "指定した遺伝子がadata.var_namesにありません。"
        )

    missing_genes = [
        gene
        for gene in genes
        if gene not in adata.var_names
    ]

    if missing_genes:
        print("存在しない遺伝子:", missing_genes)

    df = pd.DataFrame(
        index=adata.obs_names
    )

    df["pt_raw"] = pd.to_numeric(
        adata.obs[pseudotime_key],
        errors="coerce",
    )

    df = (
        df.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["pt_raw"])
        .copy()
    )

    if df.empty:
        raise RuntimeError(
            "有効なPseudotimeを持つ細胞がありません。"
        )

    pt_min = df["pt_raw"].min()
    pt_max = df["pt_raw"].max()

    if normalize_pseudotime:
        if pt_max == pt_min:
            raise RuntimeError(
                "Pseudotimeの最大値と最小値が同じです。"
            )

        df["pt"] = (
            (df["pt_raw"] - pt_min)
            / (pt_max - pt_min)
        )
    else:
        df["pt"] = df["pt_raw"]

    for gene in available_genes:
        expression = get_expr(
            adata,
            gene,
        )

        df[gene] = (
            expression
            .reindex(df.index)
            .values
        )

        df[f"{gene}_pos"] = (
            df[gene] > 0
        ).astype(float)

    df = df.sort_values("pt").copy()

    actual_bins = min(
        n_bins,
        df.shape[0],
    )

    df["bin"] = pd.qcut(
        df["pt"],
        q=actual_bins,
        duplicates="drop",
    )

    grouped = df.groupby(
        "bin",
        observed=True,
    )

    bin_centers = grouped["pt"].mean()
    mean_expression = grouped[
        available_genes
    ].mean()

    positive_columns = [
        f"{gene}_pos"
        for gene in available_genes
    ]

    fraction_positive = grouped[
        positive_columns
    ].mean()

    fraction_positive.columns = (
        available_genes
    )

    mean_smooth = mean_expression.rolling(
        window=smooth_window,
        center=True,
        min_periods=1,
    ).mean()

    fraction_smooth = fraction_positive.rolling(
        window=smooth_window,
        center=True,
        min_periods=1,
    ).mean()

    x = bin_centers.to_numpy()

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    fig = plt.figure(
        figsize=(7.0, 2.7)
    )

    grid = fig.add_gridspec(
        1,
        2,
        wspace=0.24,
    )

    ax_mean = fig.add_subplot(
        grid[0, 0]
    )
    ax_fraction = fig.add_subplot(
        grid[0, 1]
    )

    for gene in available_genes:
        gene_label = rf"$\it{{{gene}}}$"

        ax_mean.plot(
            x,
            mean_smooth[gene].to_numpy(),
            linewidth=1.7,
            label=gene_label,
        )

        ax_fraction.plot(
            x,
            fraction_smooth[gene].to_numpy(),
            linewidth=1.7,
            label=gene_label + r"$^{+}$",
        )

    ax_mean.set_title(
        "Mean expression",
        fontweight="bold",
    )
    ax_mean.set_ylabel("Expression")
    ax_mean.set_xlabel("Pseudotime")
    ax_mean.legend(
        frameon=False,
        loc="upper left",
    )
    ax_mean.grid(alpha=0.25)

    ax_fraction.set_title(
        "Proportion",
        fontweight="bold",
    )
    ax_fraction.set_ylabel("Fraction")
    ax_fraction.set_xlabel("Pseudotime")
    ax_fraction.set_ylim(-0.02, 0.6)
    ax_fraction.legend(
        frameon=False,
        loc="upper left",
    )
    ax_fraction.grid(alpha=0.25)

    for axis in [
        ax_mean,
        ax_fraction,
    ]:
        axis.set_xlim(xlim)

    output_dir = (
        get_project_dir()
        / output_folder
    )
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    pdf_path = (
        output_dir
        / f"{output_name}.pdf"
    )

    fig.savefig(
        pdf_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
    plt.close(fig)

    cells_in_range = df[
        df["pt"].between(
            xlim[0],
            xlim[1],
            inclusive="both",
        )
    ].shape[0]

    print("Saved:")
    print(pdf_path)
    print("\nDisplayed pseudotime range:", xlim)
    print("Total cells:", df.shape[0])
    print("Cells in displayed range:", cells_in_range)
    print("Available genes:", available_genes)

    return df, pdf_path
