#!/usr/bin/env python3
"""Rebuild the public human scRNA-seq components of manuscript Figure 8.

This script deliberately consumes only public GSE153760 matrices and writes
derived outputs outside Git. It is the canonical non-notebook entry point for
the submission-stage repository.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.metadata
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from neuro import (  # noqa: E402
    build_cdc2_subsets,
    compute_cdc2_pseudotime,
    load_10x_samples,
    plot_cdc2_proportion,
    plot_clec10a_cdc2_expression_violin,
    plot_pseudotime_trajectory,
    run_qc_and_global_umap,
)
from neuro.composition import _mann_whitney_pvalue  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        required=True,
        type=Path,
        help="Directory containing the per-sample GSE153760 10x matrix files.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="New or empty local directory for derived figures and tables.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "config" / "analysis_parameters.yaml",
        help="Versioned YAML analysis configuration.",
    )
    return parser.parse_args()


def read_config(path: Path) -> dict:
    with path.open() as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration is not a mapping: {path}")
    return config


def input_inventory(data_dir: Path) -> list[dict[str, int | str]]:
    records = []
    for path in sorted(data_dir.glob("*")):
        if path.is_file() and not path.name.startswith("._"):
            records.append({"path": path.name, "bytes": path.stat().st_size})
    return records


def package_versions() -> dict[str, str]:
    packages = ["anndata", "matplotlib", "numpy", "pandas", "scanpy", "scipy", "seaborn", "statsmodels"]
    return {package: importlib.metadata.version(package) for package in packages}


def write_manifest(output_dir: Path, args: argparse.Namespace, config: dict) -> None:
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository_root": str(ROOT),
        "data_dir": str(args.data_dir.resolve()),
        "config": config,
        "input_inventory": input_inventory(args.data_dir),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": package_versions(),
        "statistical_note": (
            "Primary HC-versus-AD inference is at sample level. Cell-level "
            "tests are retained only as descriptive secondary analyses."
        ),
    }
    (output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def run(args: argparse.Namespace) -> None:
    if not args.data_dir.is_dir():
        raise NotADirectoryError(args.data_dir)
    config = read_config(args.config)
    output_dir = args.output_dir.resolve()
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    write_manifest(output_dir, args, config)

    qc = config["qc"]
    embedding = config["embedding"]
    gate = config["cdc2_gate"]
    trajectory = config["pseudotime"]
    selected_clusters = tuple(gate["retained_leiden_clusters"])

    raw = load_10x_samples(args.data_dir, ad_pattern=config["dataset"]["ad_sample_pattern"])
    adata = run_qc_and_global_umap(raw, **qc, **embedding)
    immune, _, cdc2 = build_cdc2_subsets(
        adata,
        immune_gene=gate["immune_gene"],
        dc_genes=tuple(gate["dc_enrichment_genes"]),
        main_gate_genes=tuple(gate["required_positive_genes"]),
        **embedding,
    )
    del immune

    expression, patient_mean, expression_stats, _ = plot_clec10a_cdc2_expression_violin(
        cdc2,
        selected_clusters=selected_clusters,
        output_folder=str(figures_dir),
        output_name="figure8_cdc2_expression",
    )
    expression.to_csv(tables_dir / "cell_level_expression_descriptive.csv", index=False)
    patient_mean.to_csv(tables_dir / "sample_level_mean_expression.csv", index=False)
    expression_stats.to_csv(tables_dir / "expression_statistics.csv", index=False)

    # Build the sample-level composition table directly from the retained cDC2
    # subset. A sample is the experimental unit for the primary comparison.
    from neuro import summarize_cdc2_proportion  # imported lazily for a compact CLI surface

    proportion = summarize_cdc2_proportion(adata, cdc2, selected_clusters=selected_clusters)
    proportion.to_csv(tables_dir / "sample_level_cdc2_proportion.csv", index=False)
    fraction_p = _mann_whitney_pvalue(
        proportion.loc[proportion["group"] == "HC", "frac_cdc2_in_immune"],
        proportion.loc[proportion["group"] == "AD", "frac_cdc2_in_immune"],
    )
    pd.DataFrame(
        [{"outcome": "CLEC10A_positive_cDC2_fraction", "unit": "sample", "test": "two-sided Mann-Whitney U", "p_value": fraction_p}]
    ).to_csv(tables_dir / "composition_statistics.csv", index=False)
    plot_cdc2_proportion(
        proportion,
        output_folder=str(figures_dir),
        output_name="figure8_cdc2_proportion",
    )

    trajectory_adata = compute_cdc2_pseudotime(
        cdc2,
        selected_clusters=selected_clusters,
        root_positive_gene=trajectory["root_positive_gene"],
        root_negative_genes=tuple(trajectory["root_negative_genes"]),
        root_negative_weights=tuple(trajectory["root_negative_weights"]),
        n_top_genes=embedding["n_top_genes"],
        n_neighbors=trajectory["n_neighbors"],
        n_pcs=trajectory["n_pcs"],
        min_dist=trajectory["min_dist"],
    )
    trajectory_adata.obs[["pt", "root_score", "is_root", "group", "sample", "leiden_cdc2"]].to_csv(
        tables_dir / "pseudotime_cell_metadata.csv"
    )
    trajectory_table, _ = plot_pseudotime_trajectory(
        trajectory_adata,
        n_bins=trajectory["n_bins"],
        smooth_window=trajectory["smooth_window"],
        xlim=tuple(trajectory["early_trajectory_xlim"]),
        output_folder=str(figures_dir),
        output_name="figure8_pseudotime_early",
    )
    trajectory_table.to_csv(tables_dir / "pseudotime_gene_trajectory.csv", index=False)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "run.log").open("w") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        run(args)
    print(f"Completed Figure 8 rebuild: {args.output_dir}")


if __name__ == "__main__":
    main()
