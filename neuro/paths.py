from __future__ import annotations

from pathlib import Path


def get_project_dir() -> Path:
    """プロジェクトのルートフォルダを返す。"""
    return Path(__file__).resolve().parent.parent


def get_data_dir(folder_name: str) -> Path:
    """プロジェクト内のデータフォルダを返す。"""
    data_dir = get_project_dir() / folder_name

    if not data_dir.exists():
        raise FileNotFoundError(
            f"フォルダが見つかりません: {data_dir}"
        )

    if not data_dir.is_dir():
        raise NotADirectoryError(
            f"フォルダではありません: {data_dir}"
        )

    return data_dir


def list_files(
    folder_name: str,
    limit: int | None = None,
) -> list[str]:
    """指定フォルダ内のファイル名を返す。"""
    if limit is not None and limit < 0:
        raise ValueError("limitは0以上にしてください。")

    data_dir = get_data_dir(folder_name)

    files = sorted(
        path.name
        for path in data_dir.iterdir()
        if path.is_file()
        and not path.name.startswith("._")
    )

    if limit is not None:
        files = files[:limit]

    return files