from __future__ import annotations
from pathlib import Path
from typing import Sequence
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from scipy.stats import mannwhitneyu
from .paths import get_project_dir
from .subsets import get_expr, show_zero_immune_samples
def _pvalue_to_symbol(p_value: float) -> str:
    """p値を有意差記号へ変換する。"""
    if pd.isna(p_value):
        return "NA"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"
def _mann_whitney_pvalue(
    group_1: Sequence[float],
    group_2: Sequence[float],
) -> float:
    """有限値のみを用いてMann–Whitney U検定を行う。"""
    values_1 = np.asarray(group_1, dtype=float)
    values_2 = np.asarray(group_2, dtype=float)

    values_1 = values_1[np.isfinite(values_1)]
    values_2 = values_2[np.isfinite(values_2)]

    if len(values_1) < 2 or len(values_2) < 2:
        return np.nan

    return float(
        mannwhitneyu(
            values_1,
            values_2,
            alternative="two-sided",
        ).pvalue
    )
def summarize_cdc2_proportion(
    adata_all: AnnData,
    adata_cdc2: AnnData,
    *,
    selected_clusters: Sequence[str],
    sample_key: str = "sample",
    group_key: str = "group",
    immune_gene: str = "PTPRC",
    cluster_key: str = "leiden_cdc2",
) -> pd.DataFrame:
    """サンプルごとのcDC2/免疫細胞割合を集計する。"""
    required_obs = {sample_key, group_key}
    missing_obs = required_obs - set(adata_all.obs.columns)

    if missing_obs:
        raise KeyError(
            f"adata_all.obsに必要な列がありません: {sorted(missing_obs)}"
        )

    if cluster_key not in adata_cdc2.obs.columns:
        raise KeyError(
            f"adata_cdc2.obsに'{cluster_key}'がありません。"
        )

    selected_clusters = {
        str(cluster)
        for cluster in selected_clusters
    }

    cell_table = pd.DataFrame(
        index=adata_all.obs_names
    )

    cell_table["sample"] = (
        adata_all.obs[sample_key].astype(str)
    )
    cell_table["group"] = (
        adata_all.obs[group_key].astype(str)
    )

    immune_expression = get_expr(
        adata_all,
        immune_gene,
    )
    cell_table["immune"] = (
        immune_expression.reindex(cell_table.index) > 0
    )

    selected_cdc2 = (
        adata_cdc2.obs[cluster_key]
        .astype(str)
        .isin(selected_clusters)
    )

    cell_table["cdc2"] = False

    common_cells = cell_table.index.intersection(
        adata_cdc2.obs_names
    )

    cell_table.loc[common_cells, "cdc2"] = (
        selected_cdc2
        .reindex(common_cells)
        .fillna(False)
        .astype(bool)
        .values
    )

    summary = (
        cell_table
        .groupby(
            ["sample", "group"],
            observed=True,
        )
        .agg(
            n_all=("immune", "size"),
            n_immune=("immune", "sum"),
            n_cdc2=("cdc2", "sum"),
        )
        .reset_index()
    )

    summary["frac_cdc2_in_immune"] = np.where(
        summary["n_immune"] > 0,
        summary["n_cdc2"] / summary["n_immune"],
        np.nan,
    )

    show_zero_immune_samples(summary)

    return summary
def plot_cdc2_proportion(
    summary: pd.DataFrame,
    *,
    group_order: tuple[str, str] = ("HC", "AD"),
    fraction_col: str = "frac_cdc2_in_immune",
    output_folder: str = "FIG/Fig_B",
    output_name: str = "cDC2_proportion",
) -> tuple[float, float, Path]:
    """cDC2割合の箱ひげ図を作成し、MWU検定を行う。"""
    group_1, group_2 = group_order

    values_1 = summary.loc[
        summary["group"] == group_1,
        fraction_col,
    ].to_numpy()

    values_2 = summary.loc[
        summary["group"] == group_2,
        fraction_col,
    ].to_numpy()

    p_value = _mann_whitney_pvalue(
        values_1,
        values_2,
    )
    symbol = _pvalue_to_symbol(p_value)

    mean_difference = (
        np.nanmean(values_2)
        - np.nanmean(values_1)
    )

    print(f"MWU p = {p_value:.4g}")
    print(f"Significance = {symbol}")
    print(
        f"Mean difference "
        f"({group_2} - {group_1}) = "
        f"{mean_difference:.4f}"
    )

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
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
        }
    )

    figure, axis = plt.subplots(
        figsize=(2.5, 2.9)
    )

    sns.boxplot(
        data=summary,
        x="group",
        y=fraction_col,
        order=list(group_order),
        width=0.5,
        fliersize=0,
        linewidth=1.0,
        ax=axis,
    )

    sns.stripplot(
        data=summary,
        x="group",
        y=fraction_col,
        order=list(group_order),
        color="black",
        size=3.0,
        jitter=0.12,
        edgecolor="white",
        linewidth=0.35,
        ax=axis,
    )

    maximum = summary[fraction_col].max()

    if not np.isfinite(maximum):
        maximum = 0.1

    bracket_height = max(
        maximum * 0.05,
        0.005,
    )
    bracket_y = maximum + bracket_height

    axis.plot(
        [0, 0, 1, 1],
        [
            bracket_y,
            bracket_y + bracket_height,
            bracket_y + bracket_height,
            bracket_y,
        ],
        linewidth=1.0,
        color="black",
    )

    axis.text(
        0.5,
        bracket_y + bracket_height,
        symbol,
        horizontalalignment="center",
        verticalalignment="bottom",
        fontsize=10.5,
    )

    axis.set_xlabel("")
    axis.set_ylabel("Fraction")
    axis.set_title(
        r"$\it{CLEC10A}^{+}$ cDC2 / "
        r"$\it{PTPRC}^{+}$ leukocytes",
        fontsize=12,
    )
    axis.set_ylim(
        0,
        bracket_y + bracket_height * 2.5,
    )

    figure.tight_layout()

    output_dir = (
        get_project_dir()
        / output_folder
    )
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path = (
        output_dir
        / f"{output_name}.pdf"
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
    plt.close(figure)

    print("Saved:", output_path)


    return p_value, mean_difference, output_path
