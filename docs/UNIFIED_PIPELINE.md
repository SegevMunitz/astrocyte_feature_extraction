# Unified morph→RNA pipeline

Self-contained ELSC root (default):

`/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end/`

Images stay in their existing folder (`paths.images_dir`). Weights, bundle, code, envs, logs, and outputs live under ROOT.

**Do not run Google Drive FOVs until you intentionally change `images_dir` / filters.**

## Stages (in order)

| Stage | What it does | Outputs |
|-------|----------------|---------|
| `inventory` | Discover 3-channel TIFFs; skip incompatible if configured | `manifests/images.json`, `manifests/skipped_images.json` |
| `segment` | Cellpose with 3ch transfer checkpoint (inference only) | `masks/*_labels.tif` |
| `crops` | Per-cell crops + masks | `cells/<image_id>/cell_*.tif` |
| `crop_measure` | Flatten crop names → CellProfiler SizeShape+Intensity + skeleton → fill `cluster` via assigner | `tables/single_cell_full_features.csv` |
| `predict_rna` | Frozen morph→RNA: all FOVs (no Time drop); labeled accuracy when Time known | `rna_predictions/predicted_*.csv`, `prediction_summary.json` |

Optional FOV `measure` still exists but is **not** on the RNA path.

## Bootstrap (once on ELSC)

```bash
# Sync this repo to the cluster first, then:
bash tools/bootstrap_root.sh
# or: ROOT=... REPO=... bash tools/bootstrap_root.sh
```

Curated copies only (one Cellpose checkpoint, bundle files, cluster CSV, code). Writes `ROOT/logs/bootstrap_manifest.txt`. Fits `cluster_assigner.joblib` when a labeled morph CSV is available.

## Submit via Slurm (not login node)

```bash
export ROOT=/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end
export CONFIG=$ROOT/configs/elsc_unified.yaml
bash $ROOT/slurm/submit_pipeline.sh
```

Logs: `$ROOT/logs/slurm/%x-%j.{out,err}`

Chain: `segment.sbatch` → `crop_measure.sbatch` → `predict_rna.sbatch` (`afterok`).

## Checking classification

After predict, open `rna_predictions/predicted_time.csv` (`true_time` vs `pred_time`) and `prediction_summary.json` (`accuracy_labeled`).
