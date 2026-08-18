# AGENTS.md

## Cursor Cloud specific instructions

This repo is a **headless Python CLI data pipeline** (`astrocyte-pipeline`, entry point
`astrocyte_feature_extraction.cli:main`). There is no web app, server, or GUI to run — "running the
application" means invoking the CLI stages. See [`README.md`](README.md) for the full stage list and
`configs/` for example configurations.

### Environment
- Python `3.12` with a virtualenv at `.venv` (the system `python3` needs the `python3.12-venv` apt
  package to create venvs, already handled during setup). The update script installs the package
  editable with the `test` extra (`pip install -e '.[test]'`).
- Activate with `source .venv/bin/activate`, or call tools directly as `.venv/bin/<tool>`.

### Running / testing / lint
- Tests: `.venv/bin/python -m pytest` (9 tests, ~0.5s). No GPU or heavy deps needed — the Cellpose
  model is monkeypatched in the tests.
- Lint: **no linter is configured** (no ruff/flake8/pylint config, `setup.cfg`, or tox). The closest
  sanity check is `.venv/bin/python -m compileall -q src scripts tests`.
- Build/run (dev): the editable install *is* the dev build; run stages via
  `.venv/bin/astrocyte-pipeline <stage> --config <yaml>`.

### Non-obvious gotchas
- The committed configs (`configs/elsc*.yaml`) point at **absolute ELSC cluster paths that do not
  exist in this environment**. To run the CLI here you must supply a config whose `paths.*` point at
  real local files/dirs. The `inventory`, `spatial`, and `crops` stages run with only the base
  dependencies.
- Optional extras are **not installed by default** and are only needed for two stages:
  - `segment` needs the `segmentation` extra (`cellpose==3.1.1.1`, pulls in Torch) and normally a GPU.
  - `measure` / `run` need a separate **CellProfiler 4.2.8** install (it is *not* a pip dependency of
    this project; the direct backend lives in `scripts/run_cellprofiler_direct.py` and expects a
    CellProfiler venv referenced from the config). Do not expect these to work without that backend.
- Pipeline stages are ordered: `crops`, `spatial`, and `measure` require segmentation masks at
  `<output_dir>/masks/<image_id>_labels.tif` first (produced by `segment`, or created manually for
  quick end-to-end checks). `image_id` is a slug + sha1 hash of the image's relative path and is
  recorded in `<output_dir>/manifests/images.json`.
- Pipeline outputs go under the configured `output_dir`; `outputs/`, `data/`, and `*.log` are
  gitignored.
