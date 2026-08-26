#!/usr/bin/env python3
"""Morphology (slice) → time → bulk RNA program scores."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from collections import Counter
from pathlib import Path

import numpy as np
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, log_loss
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

TIME_ORDER = ("ctrl", "4h", "24h", "72h", "7d")
TIME_TO_RNA = {
    "ctrl": "ctrl",
    "4h": "4H",
    "24h": "24H",
    "72h": "72H",
    "7d": "7D",
}
CELL_ALLOWLIST = (
    "AreaShape_Area",
    "AreaShape_Perimeter",
    "AreaShape_FormFactor",
    "AreaShape_Solidity",
    "AreaShape_Eccentricity",
    "AreaShape_Compactness",
    "skeleton_length_pixels",
    "nontrunk_branch_points",
    "Intensity_MeanIntensity_Channel_0",
)
SKIP_COLUMNS = {
    "sample",
    "original_cell_id",
    "compact_cell_id",
    "Time",
    "cluster",
    "measurement_backend",
    "crop_path",
    "mask_path",
}
# Crop-frame coordinates, not morphology.
SKIP_SUBSTRINGS = (
    "Center_X",
    "Center_Y",
    "BoundingBoxMinimum",
    "BoundingBoxMaximum",
)


def default_programs_path() -> Path:
    return Path(__file__).resolve().parent / "programs.yaml"


def resolve_programs_path(cli_path: Path | None = None, project_root: Path | None = None) -> Path:
    if cli_path is not None:
        if not cli_path.is_file():
            raise SystemExit(f"programs YAML not found: {cli_path}")
        return cli_path
    candidates = [default_programs_path()]
    if project_root is not None:
        candidates.extend(
            [
                project_root / "scripts" / "programs.yaml",
                project_root / "programs.yaml",
            ]
        )
    for path in candidates:
        if path.is_file():
            return path
    raise SystemExit("programs.yaml not found next to this script or under project-root")


def load_program_defs(path: Path | None = None, project_root: Path | None = None) -> dict[str, dict]:
    yaml_path = resolve_programs_path(path, project_root)
    loaded = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict) or not loaded:
        raise SystemExit(f"No programs in {yaml_path}")
    defs: dict[str, dict] = {}
    for name, spec in loaded.items():
        if not isinstance(spec, dict) or not spec.get("genes"):
            raise SystemExit(f"Program {name!r} needs a genes list in {yaml_path}")
        defs[name] = {
            "description": str(spec.get("description") or ""),
            "source": str(spec.get("source") or ""),
            "include": str(spec.get("include") or "always"),
            "genes": tuple(str(g) for g in spec["genes"]),
            "overlap_with": [str(x) for x in (spec.get("overlap_with") or [])],
            "yaml_path": str(yaml_path),
        }
    return defs


def genes_by_program(defs: dict[str, dict]) -> dict[str, tuple[str, ...]]:
    return {name: tuple(spec["genes"]) for name, spec in defs.items()}


PROGRAMS = genes_by_program(load_program_defs())


def _f(value: str) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return float("nan")
    if math.isnan(x) or math.isinf(x):
        return float("nan")
    return x


def _percentile(values: np.ndarray, q: float) -> float:
    if values.size == 0:
        return float("nan")
    return float(np.nanpercentile(values, q))


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 5:
        return 0.0
    a = a[mask]
    b = b[mask]
    if np.std(a) == 0 or np.std(b) == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def drop_collinear(names: list[str], matrix: np.ndarray, thresh: float = 0.9) -> list[int]:
    keep: list[int] = []
    for i, _name in enumerate(names):
        if any(abs(pearson(matrix[:, i], matrix[:, j])) > thresh for j in keep):
            continue
        keep.append(i)
    return keep


def is_feature_column(name: str) -> bool:
    if name in SKIP_COLUMNS:
        return False
    if any(token in name for token in SKIP_SUBSTRINGS):
        return False
    return name.startswith("AreaShape_") or name.startswith("Intensity_") or name in {
        "skeleton_length_pixels",
        "nontrunk_branch_points",
    }


def load_ranked_features(path: Path, k: int) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    names = [row["feature"] for row in rows if row.get("feature")]
    if not names:
        raise SystemExit(f"No features in {path}")
    return names[:k]


def parse_top_k_names(
    ranked: list[str],
) -> tuple[list[str], list[str], tuple[str, ...]]:
    helper_times = tuple(t for t in TIME_ORDER if f"mean_Pcell_{t}" in ranked)
    morph_pooled = [name for name in ranked if not name.startswith("mean_Pcell_")]
    cell_feats: list[str] = []
    seen: set[str] = set()
    for name in morph_pooled:
        if name == "n_cells" or name.startswith("frac_cluster_"):
            continue
        if name.startswith("median_"):
            cell = name[len("median_") :]
        elif name.startswith("p90_"):
            cell = name[len("p90_") :]
        else:
            continue
        if cell not in seen:
            seen.add(cell)
            cell_feats.append(cell)
    return cell_feats, morph_pooled, helper_times


def read_cells(
    path: Path,
    all_features: bool,
    feature_names: list[str] | None = None,
) -> tuple[list[str], dict[str, np.ndarray], np.ndarray, np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    header = list(rows[0].keys()) if rows else []
    if feature_names is not None:
        features = [name for name in feature_names if name in header]
    elif all_features:
        features = [name for name in header if is_feature_column(name)]
    else:
        features = [name for name in CELL_ALLOWLIST if name in header]
    if not features:
        raise SystemExit(f"No usable feature columns in {path}")

    images = np.array([row["sample"] for row in rows], dtype=object)
    times = np.array([row["Time"] for row in rows], dtype=object)
    clusters = np.array([row.get("cluster", "") for row in rows], dtype=object)
    by_feat = {
        name: np.array([_f(row.get(name, "")) for row in rows], dtype=float) for name in features
    }
    return features, by_feat, images, times, clusters


def image_rows(
    features: list[str],
    by_feat: dict[str, np.ndarray],
    images: np.ndarray,
    times: np.ndarray,
    clusters: np.ndarray,
    collinear_thresh: float | None = 0.9,
    wanted_cols: list[str] | None = None,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, list[str]]:
    unique_images = sorted(set(images.tolist()), key=str)
    cluster_levels = sorted({c for c in clusters.tolist() if str(c) != ""}, key=str)
    wanted = set(wanted_cols) if wanted_cols is not None else None
    col_names: list[str] = []
    for name in features:
        median_name = f"median_{name}"
        p90_name = f"p90_{name}"
        if wanted is None or median_name in wanted:
            col_names.append(median_name)
        if wanted is None or p90_name in wanted:
            col_names.append(p90_name)
    if wanted is None or "n_cells" in wanted:
        col_names.append("n_cells")
    for cluster in cluster_levels:
        frac_name = f"frac_cluster_{cluster}"
        if wanted is None or frac_name in wanted:
            col_names.append(frac_name)
    if wanted is not None:
        col_names = [name for name in wanted_cols or [] if name in col_names]

    x_rows = []
    y_time = []
    image_ids = []
    for image in unique_images:
        mask = images == image
        time_vals = times[mask]
        time = Counter(time_vals.tolist()).most_common(1)[0][0]
        if time not in TIME_ORDER:
            continue
        vals_by_feat = {name: by_feat[name][mask] for name in features}
        n = int(mask.sum())
        cl = clusters[mask]
        row: list[float] = []
        for col in col_names:
            if col == "n_cells":
                row.append(float(n))
            elif col.startswith("frac_cluster_"):
                cluster = col[len("frac_cluster_") :]
                row.append(float(np.mean(cl == cluster)))
            elif col.startswith("median_"):
                row.append(float(np.nanmedian(vals_by_feat[col[len("median_") :]])))
            elif col.startswith("p90_"):
                row.append(_percentile(vals_by_feat[col[len("p90_") :]], 90))
            else:
                row.append(float("nan"))
        x_rows.append(row)
        y_time.append(time)
        image_ids.append(image)

    x = np.asarray(x_rows, dtype=float)
    finite_std = np.array(
        [np.nanstd(x[:, j]) > 1e-12 and np.isfinite(x[:, j]).any() for j in range(x.shape[1])]
        if x.size
        else []
    )
    col_names = [name for name, ok in zip(col_names, finite_std) if ok]
    x = x[:, finite_std] if x.size else x
    if collinear_thresh is not None and x.size:
        keep = drop_collinear(col_names, x, collinear_thresh)
        col_names = [col_names[i] for i in keep]
        x = x[:, keep]
    for j in range(x.shape[1]):
        col = x[:, j]
        med = np.nanmedian(col)
        col[~np.isfinite(col)] = med if np.isfinite(med) else 0.0
        x[:, j] = col
    return col_names, x, np.array(y_time, dtype=object), np.array(image_ids, dtype=object), cluster_levels


def read_vst(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        samples = header[1:]
        genes = []
        rows = []
        for row in reader:
            genes.append(row[0])
            rows.append([_f(v) for v in row[1:]])
    return genes, samples, np.asarray(rows, dtype=float)


def read_col_groups(path: Path, samples: list[str]) -> dict[str, str]:
    groups: dict[str, str] = {}
    if not path.is_file():
        return groups
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        id_col = "sample_id" if "sample_id" in fields else fields[0]
        group_col = "group" if "group" in fields else None
        for row in reader:
            sid = row.get(id_col, "")
            if sid.startswith("X") and sid[1:] in samples:
                sid = sid[1:]
            if group_col:
                groups[sid] = row[group_col]
            elif "group" in "".join(fields).lower():
                groups[sid] = row.get("group", "")
    # R write.csv may put rownames in first column named ""
    if not groups:
        handle_groups: dict[str, str] = {}
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            keys = list(row.keys())
            sid = row.get("sample_id") or row[keys[0]]
            grp = row.get("group") or row.get("group.1") or ""
            handle_groups[sid] = grp
        groups = handle_groups
    return groups


def infer_groups_from_names(samples: list[str]) -> dict[str, str]:
    out = {}
    for sample in samples:
        name = sample
        if name.endswith("_S8") and "72h" in name:
            continue
        if name.startswith("BMP10"):
            continue
        if "BMP50_4h" in name or name.startswith("4H"):
            out[name] = "4H"
        elif "BMP50_24h" in name:
            out[name] = "24H"
        elif "BMP50_72h" in name:
            out[name] = "72H"
        elif "BMP50_7d" in name:
            out[name] = "7D"
        elif name.startswith("c_7d") or name.startswith("NS_"):
            out[name] = "ctrl"
        elif "control" in name.lower() or name.lower().startswith("c_"):
            out[name] = "ctrl"
    return out


def _sample_rna_group(sample: str, groups: dict[str, str]) -> str | None:
    grp = groups.get(sample)
    if grp not in TIME_TO_RNA.values():
        grp = groups.get(sample.lstrip("X"))
    if grp not in TIME_TO_RNA.values():
        return None
    return grp


def raw_program_scores(
    genes: list[str],
    samples: list[str],
    vst: np.ndarray,
    groups: dict[str, str],
    programs: dict[str, tuple[str, ...]],
) -> tuple[list[str], list[str], np.ndarray, dict[str, list[str]]]:
    gene_index = {g: i for i, g in enumerate(genes)}
    used: dict[str, list[str]] = {}
    scores_by_sample = []
    program_names = list(programs)
    keep_samples = []
    keep_groups = []
    for j, sample in enumerate(samples):
        grp = _sample_rna_group(sample, groups)
        if grp is None:
            continue
        row = []
        for pname, gset in programs.items():
            present = [g for g in gset if g in gene_index]
            used.setdefault(pname, present)
            if not present:
                row.append(float("nan"))
                continue
            idx = [gene_index[g] for g in present]
            row.append(float(np.nanmean(vst[idx, j])))
        scores_by_sample.append(row)
        keep_samples.append(sample)
        keep_groups.append(grp)
    raw = np.asarray(scores_by_sample, dtype=float)
    return keep_samples, keep_groups, raw, used


def program_matrix(
    genes: list[str],
    samples: list[str],
    vst: np.ndarray,
    groups: dict[str, str],
    programs: dict[str, tuple[str, ...]] | None = None,
) -> tuple[list[str], np.ndarray, dict[str, list[str]]]:
    programs = programs or PROGRAMS
    program_names = list(programs)
    _keep_samples, keep_groups, raw, used = raw_program_scores(
        genes, samples, vst, groups, programs
    )
    if raw.size == 0:
        empty = np.full((len(TIME_TO_RNA), len(program_names)), np.nan)
        return program_names, empty, used
    mu = np.nanmean(raw, axis=0)
    sd = np.nanstd(raw, axis=0)
    sd[sd == 0] = 1.0
    z = (raw - mu) / sd
    y = []
    for grp in TIME_TO_RNA.values():
        mask = np.array([g == grp for g in keep_groups])
        if not mask.any():
            y.append([float("nan")] * len(program_names))
        else:
            y.append(np.nanmean(z[mask], axis=0).tolist())
    return program_names, np.asarray(y, dtype=float), used


def assess_programs(
    defs: dict[str, dict],
    genes: list[str],
    samples: list[str],
    vst: np.ndarray,
    groups: dict[str, str],
    min_genes: int = 3,
    min_time_std: float = 0.15,
) -> tuple[dict[str, tuple[str, ...]], list[dict[str, object]]]:
    """Drop optional programs that are missing or flat on the VST time course."""
    gene_set = set(genes)
    programs = genes_by_program(defs)
    _keep_samples, keep_groups, raw, used = raw_program_scores(
        genes, samples, vst, groups, programs
    )
    time_std = {}
    if raw.size:
        for j, pname in enumerate(programs):
            means = []
            for grp in TIME_TO_RNA.values():
                mask = np.array([g == grp for g in keep_groups])
                vals = raw[mask, j] if mask.any() else np.array([])
                if vals.size == 0 or not np.isfinite(vals).any():
                    means.append(float("nan"))
                else:
                    means.append(float(np.nanmean(vals)))
            vals = np.asarray(means, dtype=float)
            time_std[pname] = float(np.nanstd(vals)) if np.isfinite(vals).sum() >= 2 else 0.0
    else:
        time_std = {name: 0.0 for name in programs}

    kept: dict[str, tuple[str, ...]] = {}
    rows: list[dict[str, object]] = []
    for name, spec in defs.items():
        present = used.get(name, [g for g in spec["genes"] if g in gene_set])
        missing = [g for g in spec["genes"] if g not in gene_set]
        std = time_std.get(name, 0.0)
        include = spec.get("include", "always")
        if include == "if_expressed":
            keep = len(present) >= min_genes and std >= min_time_std
            reason = "kept" if keep else (
                f"dropped: {len(present)} genes in VST, time_std={std:.3f} "
                f"(need ≥{min_genes} genes and time_std≥{min_time_std})"
            )
        else:
            keep = len(present) >= 1
            reason = "kept" if keep else "dropped: no genes present in VST"
        if keep:
            kept[name] = tuple(present) if present else tuple(spec["genes"])
        rows.append(
            {
                "program": name,
                "include": include,
                "n_genes": len(spec["genes"]),
                "n_present": len(present),
                "present": present,
                "missing": missing,
                "time_std_vst": std,
                "kept": keep,
                "reason": reason,
            }
        )
    if not kept:
        raise SystemExit("No RNA programs passed the viability filter")
    return kept, rows


def make_pipeline(C: float) -> Pipeline:
    return Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    penalty="l2",
                    C=C,
                    solver="lbfgs",
                    max_iter=4000,
                    class_weight="balanced",
                ),
            ),
        ]
    )


def fit_cell_helper(
    x_cell: np.ndarray,
    y_cell: np.ndarray,
    classes: tuple[str, ...],
    C: float = 1.0,
) -> Pipeline:
    y_idx = np.array([classes.index(t) for t in y_cell])
    pipe = make_pipeline(C)
    pipe.fit(x_cell, y_idx)
    return pipe


def pooled_with_helper(
    x_image: np.ndarray,
    image_ids: np.ndarray,
    cell_images: np.ndarray,
    cell_x: np.ndarray,
    helper: Pipeline,
    classes: tuple[str, ...],
    helper_times: tuple[str, ...] | None = None,
) -> np.ndarray:
    keep_times = helper_times if helper_times is not None else classes
    idxs = [classes.index(t) for t in keep_times]
    proba = helper.predict_proba(cell_x)
    extra = []
    for image in image_ids:
        mask = cell_images == image
        if mask.any():
            extra.append(proba[mask].mean(axis=0)[idxs])
        else:
            extra.append(np.full(len(idxs), 1.0 / max(len(idxs), 1)))
    return np.hstack([x_image, np.asarray(extra, dtype=float)])


def rank_by_coef(x: np.ndarray, y_idx: np.ndarray, names: list[str], C: float, k: int) -> tuple[list[str], np.ndarray]:
    pipe = make_pipeline(C)
    pipe.fit(x, y_idx)
    importance = np.mean(np.abs(pipe.named_steps["clf"].coef_), axis=0)
    order = np.argsort(importance)[::-1][: min(k, len(names))]
    return [names[i] for i in order], order


def nested_oof(
    x_img: np.ndarray,
    names: list[str],
    y_idx: np.ndarray,
    y_img: np.ndarray,
    image_ids: np.ndarray,
    cell_images: np.ndarray,
    cell_times: np.ndarray,
    cell_x: np.ndarray,
    classes: tuple[str, ...],
    use_helper: bool,
    top_k: int,
    folds: StratifiedKFold,
) -> tuple[np.ndarray, np.ndarray, list[list[str]], bool]:
    oof_pred = np.zeros(len(y_img), dtype=int)
    oof_proba = np.zeros((len(y_img), len(classes)), dtype=float)
    fold_features: list[list[str]] = []
    helper_used = False
    for fold_i, (tr, te) in enumerate(folds.split(x_img, y_idx), start=1):
        x_tr, x_te = x_img[tr], x_img[te]
        feat_names = list(names)
        if use_helper:
            train_images = set(image_ids[tr].tolist())
            cell_mask = np.array([img in train_images for img in cell_images])
            y_cell = cell_times[cell_mask]
            ok = np.array([t in classes for t in y_cell])
            helper = fit_cell_helper(cell_x[cell_mask][ok], y_cell[ok], classes, C=1.0)
            x_tr = pooled_with_helper(x_tr, image_ids[tr], cell_images, cell_x, helper, classes)
            x_te = pooled_with_helper(x_te, image_ids[te], cell_images, cell_x, helper, classes)
            feat_names = names + [f"mean_Pcell_{t}" for t in classes]
            helper_used = True
        ranked, order = rank_by_coef(x_tr, y_idx[tr], feat_names, C=1.0, k=top_k)
        fold_features.append(ranked)
        try:
            C_fold = choose_C(x_tr[:, order], y_img[tr], classes)
        except ValueError:
            C_fold = 3.0
        pipe = make_pipeline(C_fold)
        pipe.fit(x_tr[:, order], y_idx[tr])
        oof_pred[te] = pipe.predict(x_te[:, order])
        oof_proba[te] = pipe.predict_proba(x_te[:, order])
        print(f"nested fold {fold_i}: C={C_fold} kept {len(ranked)} features", flush=True)
    return oof_pred, oof_proba, fold_features, helper_used


def choose_C(x: np.ndarray, y: np.ndarray, classes: tuple[str, ...]) -> float:
    y_idx = np.array([classes.index(t) for t in y])
    folds = StratifiedKFold(n_splits=min(5, max(2, int(np.min(np.bincount(y_idx))))), shuffle=True, random_state=0)
    best_c, best_acc = 1.0, -1.0
    for C in (0.1, 0.3, 1.0, 3.0, 10.0):
        accs = []
        for tr, te in folds.split(x, y_idx):
            pipe = make_pipeline(C)
            pipe.fit(x[tr], y_idx[tr])
            pred = pipe.predict(x[te])
            accs.append(accuracy_score(y_idx[te], pred))
        mean_acc = float(np.mean(accs))
        if mean_acc > best_acc:
            best_acc, best_c = mean_acc, C
    return best_c


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_confusion(cm: np.ndarray, classes: tuple[str, ...], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks(range(len(classes)))
    ax.set_yticks(range(len(classes)))
    ax.set_xticklabels(classes)
    ax.set_yticklabels(classes)
    ax.set_xlabel("Predicted time")
    ax.set_ylabel("True time")
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, int(cm[i, j]), ha="center", va="center")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_top_coefs(coef: np.ndarray, feat_names: list[str], classes: tuple[str, ...], path: Path, n: int = 25) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    importance = np.mean(np.abs(coef), axis=0)
    order = np.argsort(importance)[::-1][:n]
    fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(order))))
    ax.barh([feat_names[i] for i in order][::-1], importance[order][::-1], color="#6a3d9a")
    ax.set_xlabel("mean |coefficient| across times")
    ax.set_title("Top slice features (no image IDs)")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_programs(actual: np.ndarray, pred: np.ndarray, names: list[str], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(names)
    ncols = n if n <= 5 else 4
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(3.2 * ncols, 2.7 * nrows), sharey=True, squeeze=False
    )
    flat = axes.ravel()
    for i, (name, a, p) in enumerate(zip(names, actual.T, pred.T)):
        ax = flat[i]
        ax.plot(TIME_ORDER, a, "-o", color="#444444", label="RNA (actual)")
        ax.plot(TIME_ORDER, p, "-o", color="#6a3d9a", label="from morphology")
        ax.set_title(name.replace("_", " "), fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=8)
        ax.axhline(0, color="#cccccc", lw=0.8)
    for ax in flat[n:]:
        ax.set_visible(False)
    axes[0, 0].set_ylabel("program z-score")
    flat[min(n - 1, ncols - 1)].legend(fontsize=7, loc="best")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq"),
    )
    parser.add_argument(
        "--morph-csv",
        type=Path,
        default=Path(
            "/ems/elsc-labs/habib-n/segev.munitz/astroseg_data/analysis/results/single_cell_full_features.csv"
        ),
    )
    parser.add_argument("--vst-csv", type=Path, default=None)
    parser.add_argument("--coldata-csv", type=Path, default=None)
    parser.add_argument("--use-cell-helper", action="store_true", default=True)
    parser.add_argument("--no-cell-helper", action="store_true")
    parser.add_argument("--all-features", action="store_true")
    parser.add_argument("--top-k", type=int, default=None)
    parser.add_argument(
        "--nested",
        action="store_true",
        help="Select top-k features inside each training fold (honest CV).",
    )
    parser.add_argument("--feature-importance", type=Path, default=None)
    parser.add_argument("--results-subdir", type=str, default=None)
    parser.add_argument("--programs-yaml", type=Path, default=None)
    args = parser.parse_args()
    root = args.project_root
    nested_mode = bool(args.nested)
    top_k_mode = args.top_k is not None and not nested_mode
    nested_k = args.top_k if args.top_k is not None else 35
    if args.results_subdir:
        subdir = args.results_subdir
    elif nested_mode:
        subdir = "results_nested_top_features"
    elif top_k_mode:
        subdir = "results_top_features"
    elif args.all_features:
        subdir = "results_all_features"
    else:
        subdir = "results"
    results = root / subdir
    results.mkdir(parents=True, exist_ok=True)
    shutil.copy2(resolve_programs_path(args.programs_yaml, project_root=root), results / "programs.yaml")
    use_helper = args.use_cell_helper and not args.no_cell_helper
    helper_times: tuple[str, ...] | None = None
    wanted_cols: list[str] | None = None
    selected_cell: list[str] | None = None

    if top_k_mode:
        importance_path = args.feature_importance or (root / "results_all_features" / "feature_importance.csv")
        ranked = load_ranked_features(importance_path, args.top_k)
        selected_cell, wanted_cols, helper_times = parse_top_k_names(ranked)
        use_helper = use_helper and bool(helper_times)
        if not selected_cell and not wanted_cols:
            raise SystemExit(f"Top-k={args.top_k} mapped to no morph columns in {importance_path}")

    features, by_feat, cell_images, cell_times, cell_clusters = read_cells(
        args.morph_csv,
        all_features=args.all_features or top_k_mode or nested_mode,
        feature_names=selected_cell,
    )
    if not args.all_features and not top_k_mode and not nested_mode:
        cell_keep = drop_collinear(features, np.column_stack([by_feat[n] for n in features]), 0.9)
        features = [features[i] for i in cell_keep]
        by_feat = {n: by_feat[n] for n in features}

    if args.all_features or nested_mode:
        collinear_thresh: float | None = 0.98
    elif top_k_mode:
        collinear_thresh = None
    else:
        collinear_thresh = 0.9
    names, x_img, y_img, image_ids, cluster_levels = image_rows(
        features,
        by_feat,
        cell_images,
        cell_times,
        cell_clusters,
        collinear_thresh=collinear_thresh,
        wanted_cols=wanted_cols,
    )
    classes = TIME_ORDER
    y_idx = np.array([classes.index(t) for t in y_img])

    cell_x = np.column_stack([by_feat[n] for n in features])
    for j in range(cell_x.shape[1]):
        med = np.nanmedian(cell_x[:, j])
        cell_x[~np.isfinite(cell_x[:, j]), j] = med if np.isfinite(med) else 0.0

    vst_csv = args.vst_csv or (root / "data" / "vst_matrix.csv")
    coldata_csv = args.coldata_csv or (root / "data" / "colData.csv")
    genes, samples, vst = read_vst(vst_csv)
    groups = read_col_groups(coldata_csv, samples)
    inferred = infer_groups_from_names(samples)
    for sample in samples:
        if sample not in groups or groups[sample] not in TIME_TO_RNA.values():
            if inferred.get(sample):
                groups[sample] = inferred[sample]
            elif inferred.get(sample.lstrip("X")):
                groups[sample] = inferred[sample.lstrip("X")]
    program_defs = load_program_defs(args.programs_yaml, project_root=root)
    selected_programs, viability = assess_programs(
        program_defs, genes, samples, vst, groups
    )
    program_names, y_program, used_genes = program_matrix(
        genes, samples, vst, groups, programs=selected_programs
    )

    C = 3.0
    n_splits = min(5, int(np.min(np.bincount(y_idx))))
    folds = StratifiedKFold(n_splits=max(2, n_splits), shuffle=True, random_state=0)
    fold_features: list[list[str]] = []

    if nested_mode:
        oof_pred, oof_proba, fold_features, helper_used = nested_oof(
            x_img,
            names,
            y_idx,
            y_img,
            image_ids,
            cell_images,
            cell_times,
            cell_x,
            classes,
            use_helper,
            nested_k,
            folds,
        )
    else:
        C = choose_C(x_img, y_img, classes)
        oof_pred = np.zeros(len(y_img), dtype=int)
        oof_proba = np.zeros((len(y_img), len(classes)), dtype=float)
        helper_used = False
        for tr, te in folds.split(x_img, y_idx):
            x_tr, x_te = x_img[tr], x_img[te]
            if use_helper:
                train_images = set(image_ids[tr].tolist())
                cell_mask = np.array([img in train_images for img in cell_images])
                y_cell = cell_times[cell_mask]
                ok = np.array([t in classes for t in y_cell])
                helper = fit_cell_helper(cell_x[cell_mask][ok], y_cell[ok], classes, C=1.0)
                x_tr = pooled_with_helper(
                    x_tr, image_ids[tr], cell_images, cell_x, helper, classes, helper_times
                )
                x_te = pooled_with_helper(
                    x_te, image_ids[te], cell_images, cell_x, helper, classes, helper_times
                )
                helper_used = True
            pipe = make_pipeline(C)
            pipe.fit(x_tr, y_idx[tr])
            oof_pred[te] = pipe.predict(x_te)
            oof_proba[te] = pipe.predict_proba(x_te)

    acc = float(accuracy_score(y_idx, oof_pred))
    try:
        ll = float(log_loss(y_idx, oof_proba, labels=list(range(len(classes)))))
    except ValueError:
        ll = float("nan")
    cm = confusion_matrix(y_idx, oof_pred, labels=list(range(len(classes))))
    chance = float(np.max(np.bincount(y_idx)) / len(y_idx))

    pred_prog = oof_proba @ y_program
    actual_by_time = y_program
    pred_by_time = np.vstack(
        [
            pred_prog[y_img == t].mean(axis=0) if np.any(y_img == t) else np.full(len(program_names), np.nan)
            for t in classes
        ]
    )

    # final model on all images for coefficients
    x_all = x_img
    feat_names = list(names)
    if helper_used:
        helper_all = fit_cell_helper(
            cell_x[np.array([t in classes for t in cell_times])],
            cell_times[np.array([t in classes for t in cell_times])],
            classes,
            C=1.0,
        )
        x_all = pooled_with_helper(
            x_img, image_ids, cell_images, cell_x, helper_all, classes, helper_times
        )
        keep_times = helper_times if helper_times is not None else classes
        feat_names = names + [f"mean_Pcell_{t}" for t in keep_times]
    final = make_pipeline(C)
    final.fit(x_all, y_idx)
    coef = final.named_steps["clf"].coef_

    write_csv(
        results / "image_features.csv",
        ["slice", "time"] + names,
        [
            {"slice": i + 1, "time": y_img[i], **{names[j]: x_img[i, j] for j in range(len(names))}}
            for i in range(len(y_img))
        ],
    )
    write_csv(
        results / "rna_program_scores.csv",
        ["time"] + program_names,
        [{"time": classes[i], **{program_names[j]: y_program[i, j] for j in range(len(program_names))}} for i in range(len(classes))],
    )
    write_csv(
        results / "cv_predictions.csv",
        ["slice", "true_time", "pred_time"] + [f"P_{t}" for t in classes],
        [
            {
                "slice": i + 1,
                "true_time": y_img[i],
                "pred_time": classes[oof_pred[i]],
                **{f"P_{classes[k]}": oof_proba[i, k] for k in range(len(classes))},
            }
            for i in range(len(y_img))
        ],
    )
    write_csv(
        results / "predicted_programs.csv",
        ["slice", "true_time"] + program_names,
        [
            {
                "slice": i + 1,
                "true_time": y_img[i],
                **{program_names[j]: pred_prog[i, j] for j in range(len(program_names))},
            }
            for i in range(len(y_img))
        ],
    )
    write_csv(
        results / "coefficients.csv",
        ["feature"] + list(classes),
        [
            {"feature": feat_names[j], **{classes[k]: coef[k, j] for k in range(len(classes))}}
            for j in range(len(feat_names))
        ],
    )
    importance = np.mean(np.abs(coef), axis=0)
    order = np.argsort(importance)[::-1]
    write_csv(
        results / "feature_importance.csv",
        ["rank", "feature", "mean_abs_coef"],
        [
            {"rank": rank + 1, "feature": feat_names[j], "mean_abs_coef": importance[j]}
            for rank, j in enumerate(order)
        ],
    )

    metrics = {
        "n_cells": int(len(cell_times)),
        "n_slices": int(len(image_ids)),
        "slices_per_time": {t: int(np.sum(y_img == t)) for t in classes},
        "n_cell_features_input": int(len(features)),
        "n_image_features_kept": int(len(feat_names)),
        "cluster_levels": cluster_levels,
        "cell_helper": helper_used,
        "helper_times": list(helper_times) if helper_times is not None and helper_used else (
            list(classes) if helper_used else []
        ),
        "top_k": nested_k if nested_mode else args.top_k,
        "nested": nested_mode,
        "C": C,
        "cv_accuracy": acc,
        "cv_log_loss": ll,
        "majority_chance": chance,
        "confusion": cm.tolist(),
        "classes": list(classes),
        "program_genes_used": used_genes,
        "n_programs": len(program_names),
        "programs_dropped": [row["program"] for row in viability if not row["kept"]],
        "note": "Predicted RNA is P(time|slice) mixed with bulk VST program scores. Not cell-matched transcriptomes.",
    }
    (results / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    (results / "program_viability.json").write_text(
        json.dumps(viability, indent=2, default=str), encoding="utf-8"
    )
    (results / "cell_features_used.txt").write_text("\n".join(features) + "\n", encoding="utf-8")
    (results / "image_features_used.txt").write_text("\n".join(feat_names) + "\n", encoding="utf-8")

    if nested_mode and fold_features:
        freq: Counter[str] = Counter()
        rows = []
        for fold_i, feats in enumerate(fold_features, start=1):
            freq.update(feats)
            for rank, name in enumerate(feats, start=1):
                rows.append({"fold": fold_i, "rank": rank, "feature": name})
        write_csv(results / "nested_fold_features.csv", ["fold", "rank", "feature"], rows)
        n_folds = len(fold_features)
        write_csv(
            results / "nested_feature_frequency.csv",
            ["feature", "n_folds", "fraction"],
            [
                {"feature": name, "n_folds": n, "fraction": n / n_folds}
                for name, n in freq.most_common()
            ],
        )
        (results / "image_features_used.txt").write_text(
            "\n".join(name for name, _n in freq.most_common()) + "\n",
            encoding="utf-8",
        )

    plot_confusion(cm, classes, results / "plots" / "confusion.png")
    plot_programs(actual_by_time, pred_by_time, program_names, results / "plots" / "programs_actual_vs_pred.png")
    plot_top_coefs(coef, feat_names, classes, results / "plots" / "top_features.png")

    summary = {
        "n_slices": metrics["n_slices"],
        "n_cell_features_input": metrics["n_cell_features_input"],
        "n_image_features_kept": metrics["n_image_features_kept"],
        "cv_accuracy": acc,
        "majority_chance": chance,
        "C": C,
        "cell_helper": helper_used,
        "nested": nested_mode,
        "results": str(results),
    }
    print(json.dumps(summary, indent=2))
    if nested_mode and fold_features:
        freq = Counter(name for feats in fold_features for name in feats)
        print(f"nested top-{nested_k}: features selected in all {len(fold_features)} folds:")
        for name, n in freq.most_common():
            if n == len(fold_features):
                print(f"  {n}/{len(fold_features)}  {name}")
        print("selected in most folds:")
        for name, n in freq.most_common(15):
            print(f"  {n}/{len(fold_features)}  {name}")
    elif top_k_mode:
        print("selected features:", ", ".join(feat_names))
    print("top features by |coef| (fit on all slices, descriptive only):")
    for rank, j in enumerate(order[:15], start=1):
        print(f"  {rank:2d}  {importance[j]:.3f}  {feat_names[j]}")
    print("confusion rows=true ctrl,4h,24h,72h,7d")
    print(cm)

    def _acc(path: Path) -> float | None:
        if not path.is_file():
            return None
        try:
            return float(json.loads(path.read_text(encoding="utf-8"))["cv_accuracy"])
        except (KeyError, ValueError, json.JSONDecodeError):
            return None

    short_acc = _acc(root / "results" / "metrics.json")
    full_acc = _acc(root / "results_all_features" / "metrics.json")
    leaky_acc = _acc(root / "results_top_features" / "metrics.json")
    print(
        "accuracy comparison: "
        f"short={short_acc if short_acc is not None else 'NA'}  "
        f"full={full_acc if full_acc is not None else 'NA'}  "
        f"top35_leaky={leaky_acc if leaky_acc is not None else 'NA'}  "
        f"this={acc:.4f}"
    )
    print("Wrote", results)


if __name__ == "__main__":
    main()
