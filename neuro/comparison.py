from __future__ import annotations
from pathlib import Path
from typing import Sequence
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from .paths import get_project_dir
from .subsets import get_expr
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
def _safe_mann_whitney(
    values_1: Sequence[float],
    values_2: Sequence[float],
) -> float:
    """有限値のみを用いてMann–Whitney U検定を行う。"""
    array_1 = np.asarray(values_1, dtype=float)
    array_2 = np.asarray(values_2, dtype=float)

    array_1 = array_1[np.isfinite(array_1)]
    array_2 = array_2[np.isfinite(array_2)]

    if len(array_1) < 2 or len(array_2) < 2:
        return np.nan

    return float(
        mannwhitneyu(
            array_1,
            array_2,
            alternative="two-sided",
        ).pvalue
    )
def plot_clec10a_cdc2_expression_violin(
    adata_cdc2: AnnData,
    *,
    genes: Sequence[str] = (
        "CCR2",
        "CCR7",
        "RAMP1",
    ),
    selected_clusters: Sequence[str] | None = None,
    group_key: str = "group",
    sample_key: str = "sample",
    cluster_key: str = "leiden_cdc2",
    gate_gene: str = "CLEC10A",
    group_order: tuple[str, str] = ("HC", "AD"),
    output_folder: str = "FIG/Supplement/clec10a_cdc2_expression",
    output_name: str = "clec10a_cdc2_expression_violin",
    max_display_cells: int = 2500,
    point_size: float = 2.4,
    point_alpha: float = 0.45,
    dpi: int = 300,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Path]:
    """CLEC10A陽性cDC2のAD vs HC発現をviolin plotで比較する。"""
    required_obs = {
        group_key,
        sample_key,
    }

    if selected_clusters is not None:
        required_obs.add(cluster_key)

    missing_obs = required_obs - set(adata_cdc2.obs.columns)

    if missing_obs:
        raise KeyError(
            f"adata_cdc2.obsに必要な列がありません: "
            f"{sorted(missing_obs)}"
        )

    adata_plot = adata_cdc2.copy()
    adata_plot.var_names = (
        adata_plot.var_names
        .astype(str)
        .str.upper()
    )
    adata_plot.var_names_make_unique()

    gate_gene = gate_gene.upper()

    if gate_gene not in adata_plot.var_names:
        raise ValueError(
            f"{gate_gene}がadata_cdc2.var_namesにありません。"
        )

    for key in required_obs:
        adata_plot.obs[key] = (
            adata_plot.obs[key].astype(str)
        )

    adata_plot = adata_plot[
        adata_plot.obs[group_key].isin(group_order)
    ].copy()

    if selected_clusters is not None:
        selected_clusters = {
            str(cluster)
            for cluster in selected_clusters
        }

        adata_plot = adata_plot[
            adata_plot.obs[cluster_key].isin(
                selected_clusters
            )
        ].copy()

    gate_expression = get_expr(
        adata_plot,
        gate_gene,
    )

    adata_plot = adata_plot[
        gate_expression.reindex(
            adata_plot.obs_names
        ).to_numpy() > 0
    ].copy()

    if adata_plot.n_obs == 0:
        raise RuntimeError(
            f"{gate_gene}陽性cDC2がありません。"
        )

    available_genes = [
        gene.upper()
        for gene in genes
        if gene.upper() in adata_plot.var_names
    ]

    missing_genes = sorted(
        set(gene.upper() for gene in genes)
        - set(available_genes)
    )

    if missing_genes:
        print(
            "見つからない遺伝子は除外しました:",
            missing_genes,
        )

    if not available_genes:
        raise ValueError(
            "指定した遺伝子が見つかりません。"
        )

    expression_df = pd.DataFrame(
        index=adata_plot.obs_names
    )
    expression_df["sample"] = (
        adata_plot.obs[sample_key].values
    )
    expression_df["group"] = (
        adata_plot.obs[group_key].values
    )

    if cluster_key in adata_plot.obs.columns:
        expression_df["cluster"] = (
            adata_plot.obs[cluster_key]
            .astype(str)
            .values
        )

    for gene in available_genes:
        expression_df[gene] = (
            get_expr(adata_plot, gene)
            .reindex(expression_df.index)
            .values
        )

    patient_mean = (
        expression_df
        .groupby(
            ["sample", "group"],
            observed=True,
        )[available_genes]
        .mean()
        .reset_index()
    )

    group_1, group_2 = group_order
    rows = []

    for gene in available_genes:
        cell_p = _safe_mann_whitney(
            expression_df.loc[
                expression_df["group"] == group_1,
                gene,
            ],
            expression_df.loc[
                expression_df["group"] == group_2,
                gene,
            ],
        )

        patient_p = _safe_mann_whitney(
            patient_mean.loc[
                patient_mean["group"] == group_1,
                gene,
            ],
            patient_mean.loc[
                patient_mean["group"] == group_2,
                gene,
            ],
        )

        rows.append(
            {
                "gene": gene,
                "n_cells_hc": int(
                    (expression_df["group"] == group_1).sum()
                ),
                "n_cells_ad": int(
                    (expression_df["group"] == group_2).sum()
                ),
                "n_samples_hc": int(
                    (patient_mean["group"] == group_1).sum()
                ),
                "n_samples_ad": int(
                    (patient_mean["group"] == group_2).sum()
                ),
                "cell_mannwhitney_p": cell_p,
                "patient_mean_mannwhitney_p": patient_p,
                "mean_hc": float(
                    expression_df.loc[
                        expression_df["group"] == group_1,
                        gene,
                    ].mean()
                ),
                "mean_ad": float(
                    expression_df.loc[
                        expression_df["group"] == group_2,
                        gene,
                    ].mean()
                ),
            }
        )

    statistics = pd.DataFrame(rows)

    for column in [
        "cell_mannwhitney_p",
        "patient_mean_mannwhitney_p",
    ]:
        adjusted = np.full(
            len(statistics),
            np.nan,
        )
        pvalues = statistics[column].to_numpy(
            dtype=float
        )
        valid = np.isfinite(pvalues)

        if valid.any():
            adjusted[valid] = multipletests(
                pvalues[valid],
                method="fdr_bh",
            )[1]

        statistics[
            column.replace("_p", "_fdr")
        ] = adjusted

    print("CLEC10A+ cDC2 cells:", adata_plot.n_obs)
    print(statistics.to_string(index=False))

    plot_df = expression_df.melt(
        id_vars=[
            column
            for column in [
                "sample",
                "group",
                "cluster",
            ]
            if column in expression_df.columns
        ],
        value_vars=available_genes,
        var_name="gene",
        value_name="expression",
    )

    displayed_parts = []

    for _, frame in plot_df.groupby(
        ["gene", "group"],
        observed=True,
    ):
        displayed_parts.append(
            frame.sample(
                n=min(max_display_cells, len(frame)),
                random_state=0,
            )
        )

    displayed = pd.concat(
        displayed_parts,
        axis=0,
    )

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    figure, axes = plt.subplots(
        1,
        len(available_genes),
        figsize=(4.8 * len(available_genes), 5.6),
        squeeze=False,
    )

    palette = {
        group_1: "#4C78A8",
        group_2: "#E17C05",
    }

    for index, gene in enumerate(available_genes):
        axis = axes[0, index]
        sub = plot_df[plot_df["gene"] == gene]
        display_sub = displayed[
            displayed["gene"] == gene
        ]

        sns.violinplot(
            data=sub,
            x="group",
            y="expression",
            order=list(group_order),
            hue="group",
            palette=palette,
            legend=False,
            inner="quartile",
            cut=0,
            linewidth=1,
            ax=axis,
        )

        before_collections = len(axis.collections)

        sns.stripplot(
            data=display_sub,
            x="group",
            y="expression",
            order=list(group_order),
            color="black",
            size=point_size,
            alpha=point_alpha,
            jitter=0.18,
            ax=axis,
        )

        for collection in axis.collections[
            before_collections:
        ]:
            collection.set_rasterized(True)

        stat_row = statistics.loc[
            statistics["gene"] == gene
        ].iloc[0]
        adjusted_p = stat_row[
            "cell_mannwhitney_fdr"
        ]

        maximum = sub["expression"].max()
        bracket_height = max(
            maximum * 0.05,
            0.05,
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
            linewidth=1.2,
            color="black",
        )
        axis.text(
            0.5,
            bracket_y + bracket_height,
            _pvalue_to_symbol(adjusted_p),
            ha="center",
            va="bottom",
            fontsize=16,
        )

        axis.set_title(
            rf"$\it{{{gene}}}$",
            fontsize=17,
        )
        axis.set_xlabel("")
        axis.set_ylabel(
            "Expression"
            if index == 0
            else ""
        )
        axis.tick_params(
            axis="both",
            labelsize=13,
        )
        axis.set_ylim(
            0,
            bracket_y + bracket_height * 2.5,
        )

    figure.suptitle(
        rf"$\it{{{gate_gene}}}^+$ cDC2 cells",
        y=1.03,
        fontsize=18,
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
    statistics_path = (
        output_dir
        / f"{output_name}_statistics.csv"
    )
    expression_path = (
        output_dir
        / f"{output_name}_cell_expression.csv"
    )

    figure.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )
    statistics.to_csv(
        statistics_path,
        index=False,
    )
    expression_df.to_csv(
        expression_path,
        index=True,
        index_label="cell",
    )

    plt.show()
    plt.close(figure)

    print("Saved:", output_path)
    print("Saved:", statistics_path)
    print("Saved:", expression_path)

    return (
        expression_df,
        patient_mean,
        statistics,
        output_path,
    )
