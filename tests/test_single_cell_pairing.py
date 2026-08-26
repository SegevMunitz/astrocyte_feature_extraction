from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


def _load_measure_script():
    script_path = Path(__file__).resolve().parents[1] / "scripts" / "measure_single_cell_crops_full.py"
    spec = importlib.util.spec_from_file_location("measure_single_cell_crops_full", script_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_pair_crops_and_masks_by_shared_stem(tmp_path: Path) -> None:
    module = _load_measure_script()
    crops_dir = tmp_path / "crops"
    masks_dir = tmp_path / "masks"
    crops_dir.mkdir()
    masks_dir.mkdir()
    (crops_dir / "ctrl_photo1_cell001.tif").write_bytes(b"crop")
    (masks_dir / "ctrl_photo1_cell001_mask.tif").write_bytes(b"mask")

    pairs = module._pair_crops_and_masks(crops_dir, masks_dir)
    assert len(pairs) == 1
    assert pairs[0][2] == "ctrl_photo1_cell001"


def test_pair_crops_from_existing_csv_order(tmp_path: Path) -> None:
    module = _load_measure_script()
    crops_dir = tmp_path / "crops"
    masks_dir = tmp_path / "masks"
    crops_dir.mkdir()
    masks_dir.mkdir()
    for sample in ("ctrl_a", "4h_b"):
        (crops_dir / f"{sample}.tif").write_bytes(b"crop")
        (masks_dir / f"{sample}.tif").write_bytes(b"mask")

    existing = pd.DataFrame({"sample": ["4h_b", "ctrl_a"]})
    pairs = module._pair_crops_and_masks(crops_dir, masks_dir, existing=existing)
    assert [sample for _, _, sample in pairs] == ["4h_b", "ctrl_a"]
