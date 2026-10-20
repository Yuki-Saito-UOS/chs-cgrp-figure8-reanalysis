from __future__ import annotations

from pathlib import Path

import scanpy as sc
from anndata import AnnData

from .paths import get_data_dir


def load_10x_samples(
    folder_name: str | Path,
    *,
    ad_pattern: str = "_AD",
) -> AnnData:
    """複数の10x形式サンプルを読み込み、1つのAnnDataに結合する。"""
    data_path = Path(folder_name).expanduser()
    data_dir = (
        data_path
        if data_path.is_absolute()
        else get_data_dir(str(data_path))
    )

    if not data_dir.exists():
        raise FileNotFoundError(f"フォルダが見つかりません: {data_dir}")
    if not data_dir.is_dir():
        raise NotADirectoryError(f"フォルダではありません: {data_dir}")

    barcode_files = sorted(
        path
        for path in data_dir.glob("*_barcodes.tsv.gz")
        if not path.name.startswith("._")
    )

    sample_names = [
        path.name.removesuffix("_barcodes.tsv.gz")
        for path in barcode_files
    ]

    if not sample_names:
        raise FileNotFoundError(
            f"barcodesファイルが見つかりません: {data_dir}"
        )

    adatas: list[AnnData] = []
    loaded_sample_names: list[str] = []

    for sample_name in sample_names:
        required_files = [
            data_dir / f"{sample_name}_barcodes.tsv.gz",
            data_dir / f"{sample_name}_features.tsv.gz",
            data_dir / f"{sample_name}_matrix.mtx.gz",
        ]

        missing_files = [
            path.name
            for path in required_files
            if not path.exists()
        ]

        if missing_files:
            print(
                f"Skip {sample_name}: "
                f"missing {', '.join(missing_files)}"
            )
            continue

        adata = sc.read_10x_mtx(
            data_dir,
            prefix=f"{sample_name}_",
            var_names="gene_symbols",
            cache=False,
        )

        adata.obs["sample"] = sample_name
        adata.obs["group"] = (
            "AD"
            if ad_pattern in sample_name
            else "HC"
        )

        adatas.append(adata)
        loaded_sample_names.append(sample_name)

        print(f"{sample_name}: {adata.shape}")

    if not adatas:
        raise RuntimeError(
            "読み込み可能なサンプルがありません。"
        )

    combined = sc.concat(
        adatas,
        join="outer",
        label="batch",
        keys=loaded_sample_names,
        index_unique="-",
    )

    combined.var_names_make_unique()

    print(f"\n結合後: {combined.shape}")
    print(combined.obs["group"].value_counts())
    print(
        "サンプル数:",
        combined.obs["sample"].nunique(),
    )

    return combined
