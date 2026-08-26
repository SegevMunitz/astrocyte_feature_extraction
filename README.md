# Astrocyte feature extraction

End-to-end pipeline on ELSC:

**images → Cellpose (3-channel transfer, inference only) → single-cell crops → CellProfiler (+ skeleton) → frozen morph→RNA predict**

Optional FOV-level CellProfiler `measure` still exists for spatial tables / legacy `tables/cells.csv`, but the **RNA path uses crop-level features**, not FOV `cells.csv`.

Calcium time-series analysis and Cellpose training are outside this repo.

## Unified path (preferred)

Self-contained cluster root (default):

`/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end/`

Images stay on disk under `paths.images_dir` (not copied into ROOT). Weights, morph→RNA bundle, code, envs, logs, and outputs live under ROOT.

Full stage notes: [`docs/UNIFIED_PIPELINE.md`](docs/UNIFIED_PIPELINE.md).

### Stages

| Stage | Role |
|-------|------|
| `inventory` | Discover 3-channel TIFFs; optionally skip incompatible channel counts |
| `segment` | Cellpose 3ch transfer checkpoint |
| `crops` | Per-cell crops + masks |
| `crop_measure` | Crop CellProfiler SizeShape+Intensity + skeleton; cluster assigner → `tables/single_cell_full_features.csv` |
| `predict_rna` | Frozen morph→RNA bundle; all FOVs (no Time drop); labeled accuracy when Time is known |
| `measure` | Optional FOV CellProfiler → `tables/cells.csv` (not used for RNA) |
| `run` | `inventory → segment → crops → crop_measure → predict_rna` |

### Bootstrap once on ELSC

```bash
# Sync this repo to the cluster, then:
bash tools/bootstrap_root.sh
# Optional: ROOT=... REPO=... bash tools/bootstrap_root.sh
# Refresh morph/bundle assets: bash tools/refresh_latest_assets.sh
```

Curated copies only (one Cellpose checkpoint, bundle, cluster table, code). Logs: `ROOT/logs/bootstrap_manifest.txt`.

### Submit via Slurm (not the login node)

```bash
export ROOT=/ems/elsc-labs/habib-n/segev.munitz/astrocyte_end_to_end
export CONFIG=$ROOT/configs/elsc_unified.yaml
bash $ROOT/slurm/submit_pipeline.sh
```

Chain: `segment.sbatch` → `crop_measure.sbatch` → `predict_rna.sbatch` (`afterok`).  
Logs: `$ROOT/logs/slurm/%x-%j.{out,err}`.

Config: [`configs/elsc_unified.yaml`](configs/elsc_unified.yaml).

After predict, check `rna_predictions/predicted_time.csv` (`true_time` vs `pred_time`) and `prediction_summary.json`.

## Imaging-only install (this repo)

```bash
git clone https://github.com/SegevMunitz/astrocyte_feature_extraction.git
cd astrocyte_feature_extraction
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[segmentation,test]'
```

Checkpoint provenance pins Cellpose `3.1.1.1` with training order `GFAP, GFP, DAPI`.  
ELSC configs use `cellpose.input_channel_indices: [1, 2, 0]` and `cellprofiler.image_channel_index: 1` for the current TIFF export order. Revalidate if acquisition order changes.

CellProfiler measurement (FOV or crop) expects CellProfiler `4.2.8` in a separate `.cpvenv` (or ROOT `envs/.cpvenv`).

## Why CellProfiler remains the measurement backend

The header of `AstroResultsAstrocytes_with_nuclei.csv` defines requested column names, but not all settings needed to recalculate them. This project loads Cellpose integer labels as CellProfiler objects and runs the original measurement modules. Python assembles outputs; it does not substitute approximate `scikit-image` formulas (except crop-level skeleton length / branch points used by morph→RNA).

## Legacy imaging configs and commands

Older absolute-path configs: [`configs/elsc.yaml`](configs/elsc.yaml), [`configs/elsc_one_image.yaml`](configs/elsc_one_image.yaml).

```bash
astrocyte-pipeline inventory --config configs/elsc.yaml
astrocyte-pipeline segment --config configs/elsc.yaml
astrocyte-pipeline crops --config configs/elsc.yaml
astrocyte-pipeline spatial --config configs/elsc.yaml
astrocyte-pipeline measure --config configs/elsc.yaml
# or: astrocyte-pipeline run --config configs/elsc_unified.yaml
```

The repository includes [`cellprofiler/Astrocytes_CellposeLabels.cppipe`](cellprofiler/Astrocytes_CellposeLabels.cppipe) for FOV measurement. The pipeline fails rather than silently changing versions, channel count, or feature schema.

## Cell crops

Native crops retain only pixels belonging to the selected Cellpose label. Optional `*_article64.tif` follows the cited 64×64 z-normalized prep. **RNA morph features are measured on native crops**, not on FOV masks alone.

## Outputs

Under the configured output directory:

- `masks/`: Cellpose instance labels;
- `qc/`: overlays;
- `cells/`: native cell crops and masks;
- `tables/single_cell_full_features.csv`: crop-level morph matrix for RNA predict;
- `rna_predictions/`: `predicted_time.csv`, `predicted_programs.csv`, `predicted_panel_genes.csv`, `prediction_summary.json`;
- `tables/cells.csv` / `.parquet`: FOV CellProfiler metrics (optional `measure`);
- `tables/nodes.*` / `edges.*`: spatial graph (optional);
- `manifests/`: inventory, skipped images, schema, provenance;
- `cellprofiler/`: FOV LoadData / export logs when using `measure`.

## Morph→RNA notes

- Frozen bundle: scaler + time classifier + RNA atlas (`model_bundle.joblib`, `rna_by_time.npz`). Prefer **nested** CV feature selection for production (honest ~92% time accuracy), not the leaky non-nested top-k (~98%).
- Cluster mix (`frac_cluster_*`) comes from a frozen cell→cluster assigner when labels are missing.
- Refitting after dropping a training sample: replace nested results + refit bundle (see `slurm/fit_bundle_ex72h2.slurm` pattern); then re-predict.

Analysis scripts live under [`scripts/analysis/morph_to_rnaseq/`](scripts/analysis/morph_to_rnaseq/).

## Validation

```bash
python -m pytest
```

For FOV numerical parity, compare CellProfiler exports by `(ImageNumber, ObjectNumber)` with `assert_numeric_parity` when masks match.
