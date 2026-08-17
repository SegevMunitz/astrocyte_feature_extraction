# Astrocyte feature extraction

Reproducible pipeline for:

1. segmenting static three-channel microscopy images with a custom Cellpose 3 checkpoint;
2. exporting one black-background crop per segmented astrocyte;
3. measuring the requested CellProfiler object features with the original, version-pinned pipeline;
4. producing one feature vector and global position per cell; and
5. exporting a same-image k-nearest-neighbor graph for future spatial models.

Calcium time-series analysis, classifier training, feature attribution, and final visualization are intentionally outside this first pipeline.

## Why CellProfiler remains the measurement backend

The header of `AstroResultsAstrocytes_with_nuclei.csv` defines the requested column names, but it does not encode all settings needed to recalculate those values. Perimeters, intensity edges, texture quantization, Zernike moments, neighbors, and skeleton measurements are version- and pipeline-sensitive. This project therefore loads Cellpose integer labels as CellProfiler objects and runs the original measurement modules. Python code assembles and validates the outputs; it does not substitute approximate `scikit-image` formulas.

## Installation on ELSC

```bash
git clone https://github.com/SegevMunitz/astrocyte_feature_extraction.git
cd astrocyte_feature_extraction
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[segmentation,test]'
```

The supplied checkpoint provenance pins Cellpose `3.1.1.1` with channel order
`GFAP, GFP, DAPI`. The source CellProfiler pipeline uses format revision 428;
the ELSC configuration pins CellProfiler `4.2.8`.

The current ELSC test TIFFs store GFAP at channel index `1`, not index `0`.
`cellpose.input_channel_indices: [1, 2, 0]` maps their stored planes into the
checkpoint's training order, and `cellprofiler.image_channel_index: 1` measures
GFAP intensity from the original image. Revalidate both settings if the
acquisition export order changes.

## Configuration

Cluster paths and explicit inference settings are in [`configs/elsc.yaml`](configs/elsc.yaml). Before the measurement stage:

1. Run inventory.
2. Inspect `manifests/cellprofiler_schema.json`.
3. Set `cellprofiler.expected_version`.
4. Set `paths.cellprofiler_measurement_pipeline` to a copy of the original pipeline that:
   - preserves its preprocessing, calibration, nucleus association, and measurement settings;
   - includes `LoadData`;
   - loads `Image_ObjectsFileName_Astrocytes` / `Image_ObjectsPathName_Astrocytes` as integer objects; and
   - exports the requested object table.
5. If checkpoint provenance gives an exact Cellpose package version, set `cellpose.expected_version`.

The repository includes
[`cellprofiler/Astrocytes_CellposeLabels.cppipe`](cellprofiler/Astrocytes_CellposeLabels.cppipe),
generated from the source `Astrocytes_BMP.cppipe`. It replaces the source
pipeline's obsolete two-channel segmentation with `LoadData`, while preserving
the exact GFAP intensity and advanced size/shape measurement module settings.
The pipeline measures every Cellpose label; the historical filename mentions
nuclei, but its first-line feature schema contains no nucleus measurement.

The pipeline fails rather than silently changing versions, channel count, channel axis, object count, or feature schema.

## Commands

```bash
astrocyte-pipeline inventory --config configs/elsc.yaml
astrocyte-pipeline segment --config configs/elsc.yaml
astrocyte-pipeline crops --config configs/elsc.yaml
astrocyte-pipeline spatial --config configs/elsc.yaml
astrocyte-pipeline measure --config configs/elsc.yaml
```

After the measurement pipeline and versions are configured, the same stages can be run together:

```bash
astrocyte-pipeline run --config configs/elsc.yaml
```

For Slurm, submit the GPU stage followed by the dependent CPU measurement stage:

```bash
bash slurm/submit_pipeline.sh
```

Use [`configs/elsc_one_image.yaml`](configs/elsc_one_image.yaml) for the initial
`ctrl_20x_Multichannel_20240318_437.tif` validation before processing the full
directory.

Override paths without editing scripts:

```bash
REPO_ROOT="$PWD" CONFIG="$PWD/configs/elsc.yaml" bash slurm/submit_pipeline.sh
```

## Cell crops

Native crops retain only pixels belonging to the selected Cellpose label. All neighboring cells, background, and configurable margin pixels are black. Measurements use the full-resolution labels, not crops.

The optional `*_article64.tif` output follows the cited pre-CNN preparation: convert to 8-bit, divide by 255, force-resize to 64×64 with nearest-neighbor behavior, then z-normalize each channel independently with epsilon `1e-5`. The configurable margin is a project choice; the paper did not define a crop buffer.

## Outputs

Under the configured output directory:

- `masks/`: lossless `uint32` Cellpose instance labels;
- `qc/`: red-boundary overlays;
- `cells/`: native cell-only crops, masks, and optional 64×64 tensors;
- `tables/cells.csv` and `.parquet`: ordered CellProfiler metrics plus positions;
- `tables/nodes.*`: one spatial node per cell;
- `tables/edges.*`: directed same-image k-nearest-neighbor edges;
- `manifests/`: image inventory, schema/module mapping, hashes, settings, and software provenance;
- `cellprofiler/`: generated `LoadData` CSV, exact export, and execution log.

Coordinates use image conventions: X is column, Y is row, and the origin is the upper-left. Physical coordinates are emitted only when OME metadata or `images.pixel_size_um` supplies calibration.
Node rows also flag border contact, disconnected label fragments, and configurable minimum/maximum area violations. These are QC fields; objects are not silently discarded.

## Validation

Run local unit tests:

```bash
python -m pytest
```

For numerical parity, run both the source CellProfiler pipeline and this project's measurement pipeline on the same input images and the same Cellpose integer labels. Compare by `(ImageNumber, ObjectNumber)` with `assert_numeric_parity`. The historical CSV is the schema oracle, not a numerical oracle when it was generated from different masks.