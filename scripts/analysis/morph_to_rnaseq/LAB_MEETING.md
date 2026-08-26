# Morph → RNA: lab meeting pack (BMP50 astrocytes)

Use this document for a 10–15 minute presentation. Figures are under  
`.../morph_to_rnaseq/results_nested_top_features/plots/`.

---

## One-slide summary (copy to slides)

**Question:** Can astrocyte morphology infer bulk transcriptional programs after BMP50?

**Answer:** Yes, indirectly. Morphology predicts **time** (~**92%** nested CV slice accuracy). That time probability mixes a **frozen bulk RNA atlas** → **14 programs** recover with mean Pearson **r ≈ 0.90**.

**Story arc:** Early BMP → proliferation arrest (24h) → reactive + cytoskeletal shift (24–72h) → OPC/lineage change → late maturation (7d).

**Caveat:** Bulk programs, not single-cell RNA. Purple predicted curves are smoother than black bulk RNA (expected).

---

## Core figure (6 programs)

Use `plots/core_six_programs_actual_vs_pred.png` (not all 14 panels).

| # | Program | r | Role in story |
|---|---------|---|---------------|
| 1 | **BMP_Id** | 0.89 | BMP/SMAD on — early pulse |
| 2 | **cell_cycle** | 0.94 | Proliferation arrest at 24h |
| 3 | **reactive_GFAP** | 0.89 | Reactive structure (GFAP, Vim, Aqp4) |
| 4 | **cytoskeleton** | 0.94 | Process/IF/actin RNA — links to morph |
| 5 | **OPC_pdgfra** | 0.95 | OPC compartment (Pdgfra, Cspg4, Gpr17) |
| 6 | **late_7d** | 0.96 | Late maturation at 7d |

**Lineage note:** `notch_fate` (r = 0.81, weakest) is in supplementary text — Dll1 sustained down, Olig2 early down. Main figure uses **OPC_pdgfra** (cleaner recovery).

**Supplementary row:** `plots/supplementary_four_programs_actual_vs_pred.png` — MHC_class_I, Stat3_reactive, glutamate_ion, ECM_remodeling.

**Bar chart:** `plots/core_six_correlations.png`

---

## Actual bulk RNA trajectories (black lines)

Z-scored program scores averaged by time (`rna_program_scores.csv`):

| Program | Peak | Trough | Pattern |
|---------|------|--------|---------|
| **BMP_Id** | **4h** (+1.39) | 7d (−1.0) | Early BMP pulse, fades |
| **cell_cycle** | 7d (+1.16) | **24h** (−1.67) | Transient arrest at 24h, recovery |
| **reactive_GFAP** | **72h** (+1.64) | 4h (−0.77) | Mid-course reactive peak |
| **cytoskeleton** | **24h** (+1.20) | 4h (−0.75) | Early-mid structural shift |
| **OPC_pdgfra** | **ctrl** (+0.81) | **7d** (−1.04) | OPC markers high at baseline, fall with BMP |
| **late_7d** | **7d** (+1.37) | 4h (−1.19) | Late-only activation |
| notch_fate | ctrl (+1.24) | 24h (−1.18) | High at ctrl, suppressed early |
| MHC_class_I | 72h (+1.09) | ctrl (−0.74) | Rises with mid/late reactive window |
| Stat3_reactive | 72h (+1.85) | 4h (−1.12) | Signaling peaks at 72h |
| complement | 72h (+1.79) | ctrl (−0.20) | Innate axis active mid-late |
| interferon_response | 24h (+1.31) | 4h (−0.77) | Brief IFN-like pulse at 24h |
| glutamate_ion | 72h (+1.53) | ctrl/4h (low) | Homeostatic function up at 72h |
| ECM_remodeling | 24h (+1.24) | ctrl (−0.76) | Matrix genes spike at 24h |
| astrocyte_identity | 72h (+1.45) | 4h (−0.64) | Identity markers rise late |

**Morph agreement:** Purple lines follow the same direction for all six core programs. Smoothing is largest where slice-level `P(time)` is uncertain (e.g. 24h vs 72h).

---

## How the model works (30 seconds)

```
morphology (slice) → P(time | morph) → mix frozen bulk atlas → program scores
```

- Morphology never reads genes.
- RNA atlas built once from DESeq2 VST (`programs.yaml`, 14 gene sets).
- See [PROGRAMS.md](PROGRAMS.md) for biology per program.

---

## Morphology connection (“so what”)

Top stable morph features (nested CV, selected in ≥4/5 folds):

- Cell-helper votes: `mean_Pcell_24h`, `mean_Pcell_ctrl`, `mean_Pcell_4h`
- **GFAP intensity:** `p90_Intensity_IntegratedIntensity_Channel_0`
- **Shape:** Zernike moments, mass displacement
- **Cluster mix:** `frac_cluster_1`, `frac_cluster_2`

**Talking point:** Slices that morphology classifies as 24h–72h co-occur with elevated **reactive_GFAP** and **cytoskeleton** bulk scores — without RNA ever entering the morph model.

---

## Slide notes by program (from PROGRAMS.md)

### BMP_Id
Direct BMP/SMAD targets (Id1, Id3, Smad6, Smad7, Nog). Rises at 4h–24h; Smad6 is sustained across the time course in the DE scan.

### cell_cycle
Mki67, Top2a, Ccna2, etc. Classic **24h proliferation arrest** — dominant transient pattern in the genome-wide DE scan.

### reactive_GFAP
Gfap, Vim, Aqp4, Serpina3n, Lcn2. Structural “reactive” astrocyte state; peaks 24h–72h, not sustained at 7d.

### cytoskeleton
Vim, Nes, actin/MT, Rho–ROCK. RNA counterpart of branching and outline morph features. Overlaps Vim with reactive_GFAP — interpret together, not as independent.

### OPC_pdgfra
Pdgfra, Cspg4, Gpr17. OPC/NG2 compartment; high at ctrl, declines under BMP. Split from notch_fate for cleaner biology.

### late_7d
Dio2, Slco1c1, Thy1, Nnat. Late maturation / thyroid hormone axis — flat early, strong at 7d.

---

## Caveats (say these out loud)

1. **Not cell-matched RNA** — bulk wells and imaging slices share time labels only.
2. **Not per-cell transcriptomes** — program scores are bulk signatures mixed by predicted time.
3. **Purple smoother than black** — morphology blends time probabilities; not a bug.
4. **Overlapping programs** — Vim (reactive + cytoskeleton), Lcn2/Serpina3n (reactive + Stat3). Correlated scores are expected.
5. **Quote nested 92% only** — not the leaky 98% feature-selection run.
6. **notch_fate is softer** (r = 0.81) — mention in Q&A, not main claim.

---

## 10–15 minute outline

1. **Motivation** (1 min) — morphology is rich; RNA is expensive; can we link them?
2. **Design** (2 min) — morph → time → frozen atlas; 52 slices, 21 bulk wells, BMP50 time course.
3. **Accuracy** (1 min) — 92% time, r ≈ 0.90 programs; show bar chart.
4. **Core figure** (5 min) — walk through 6 programs in temporal order (table above).
5. **Morph features** (2 min) — GFAP intensity, Zernike, cluster mix.
6. **Caveats + next steps** (2 min) — bulk not single-cell; predict new FOVs with `bundle/`.

---

## Files to open before the meeting

| File | Purpose |
|------|---------|
| `plots/core_six_programs_actual_vs_pred.png` | Main figure |
| `plots/core_six_correlations.png` | r values for core set |
| `plots/rna_program_correlations.png` | All 14 programs |
| `rna_metrics.json` | Exact numbers |
| `program_trajectories.json` | Peak/trough table (auto-generated) |
| [PROGRAMS.md](PROGRAMS.md) | Full biology reference |

---

## Next use (after the meeting)

### Cluster rule (ELSC)

**Always `sbatch`** morph→RNA / evaluation / prediction jobs. Do **not** run long `python ...` on the login node.

Project root on cluster: `/ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq`  
Venv: `$ROOT/.venv`  
Slurm scripts (repo): `scripts/analysis/morph_to_rnaseq/slurm/` → sync to `$ROOT/scripts/slurm/`.

```bash
cd /ems/elsc-labs/habib-n/segev.munitz/morph_to_rnaseq

# Refresh frozen RNA atlas + refit bundle
sbatch scripts/slurm/refresh_rna_atlas.slurm

# Predict on new morphology only
MORPH_CSV=/path/to/NEW_CELLS.csv \
OUT_DIR=/path/to/pred_out \
sbatch scripts/slurm/predict_rna_from_morph.slurm

# Optional: fit-only
MODE=fit sbatch scripts/slurm/predict_rna_from_morph.slurm
```

Local-style CLI (compute node / Slurm job only):

```bash
python predict_rna_from_morph.py predict \
  --morph-csv NEW_CELLS.csv \
  --bundle-dir .../morph_to_rnaseq/bundle
```

Outputs: `predicted_time.csv`, `predicted_programs.csv`, `predicted_panel_genes.csv`.
