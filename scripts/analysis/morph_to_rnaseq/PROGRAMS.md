# RNA programs (morph → time → frozen bulk atlas)

These are **bulk RNA gene-set scores**, not morphology features. Morphology only predicts `P(time)`. Predicted program scores are:

`predicted = P(time | morph) · Y`

where `Y` is a frozen **5 times × K programs** matrix built from VST:

1. For each well and program, `score = mean(VST of genes in the set)`.
2. Z-score each program across wells.
3. Average wells within each time (`ctrl`, `4h`, `24h`, `72h`, `7d`).

Gene lists live in [`programs.yaml`](programs.yaml). Optional programs (`include: if_expressed`) are dropped if fewer than 3 genes are in the VST matrix or the time course is flat. Overlapping genes are intentional and listed below.

`notch_fate` no longer includes `Pdgfra` — that marker moved to `OPC_pdgfra`.

---

## BMP_Id

**Genes:** Id1, Id3, Smad6, Smad7, Nog

Direct BMP/SMAD targets and negative-feedback genes. In this BMP50 time course they rise early (4h–24h); Smad6 is one of the few genes still up at every timepoint, while Id3 can fade by 7d. This is the “BMP is on” dashboard score.

**Overlap:** none.

## cell_cycle

**Genes:** Mki67, Top2a, Ccna2, Cdk1, Rrm2, Pcna

S–G2–M proliferation. The DE scan is dominated by a **transient 24h arrest** (these genes drop, then largely recover by 7d). Expected trajectory: down at 24h, partial recovery later.

**Overlap:** none.

## reactive_GFAP

**Genes:** Gfap, Vim, Aqp4, Serpina3n, Lcn2

Structural reactive-astrocyte state (intermediate filaments, water channel) plus acute-phase genes. Gfap is transient in this experiment (mid-course, not sustained at 7d). Use this for “looks reactive,” not for Stat3 signaling per se.

**Overlap:** `Vim` with `cytoskeleton`; `Serpina3n`/`Lcn2` with `Stat3_reactive`.

## notch_fate

**Genes:** Dll1, Olig2, Sox10

Notch ligand plus oligodendrocyte-lineage transcription factors. Dll1 is suppressed at every BMP50 timepoint; Olig2 drops early. OPC surface markers were split out so this score is fate/transcription, not NG2-cell abundance.

**Overlap:** conceptually with `OPC_pdgfra` (lineage), but no shared genes.

## late_7d

**Genes:** Dio2, Slco1c1, Thy1, Nnat

Late maturation / thyroid-hormone axis. Flat early, activates mainly at 7d. Distinct from the early BMP pulse.

**Overlap:** none.

## cytoskeleton

**Genes:** Vim, Nes, Actg1, Cfl1, Tubb2a, Map1b, Rhoa, Rock1

Combined intermediate-filament, actin, microtubule, and Rho–ROCK genes. This is the RNA counterpart of morph features (branching, Zernike outline, process shape). Gfap is left in `reactive_GFAP` so this set is remodeling machinery plus nestin/vimentin, not a second GFAP score. Actb was omitted (housekeeping); Actg1/Cfl1 carry actin remodeling.

**Overlap:** `Vim` with `reactive_GFAP`.

## MHC_class_I

**Genes:** H2-D1, H2-K1, B2m, Tap1

Classical MHC class I antigen presentation (already in the contrast gene panel). Tests whether BMP shifts antigen-presentation machinery in this culture. If genes are present but the time course is weak, the score will look flat — that is a biological result, not a model failure.

**Overlap:** none.

## Stat3_reactive

**Genes:** Stat3, Lcn2, Serpina3n, Saa3

IL-6–Stat3 reactive **signaling** hub (A1-like literature), not structural GFAP. Kept separate from `reactive_GFAP` so you can see signaling vs filaments. Saa3 may be missing from VST; remaining genes still score the axis.

**Overlap:** `Lcn2`/`Serpina3n` with `reactive_GFAP`.

## complement (optional)

**Genes:** C3, C4b, C1qa, C1qb, C1qc

Innate complement, often produced by astrocytes in gliosis. Included only if the VST viability pass finds ≥3 expressed genes with a non-flat time course. C1q subunits may be absent in a relatively pure astrocyte culture.

**Overlap:** none.

## interferon_response (optional)

**Genes:** Stat1, Isg15, Ifit1, Irf7, Oasl2

Type I IFN / antiviral program, distinct from Stat3. BMP does not typically induce this; the viability pass will drop it if the time course is flat.

**Overlap:** none.

## astrocyte_identity

**Genes:** Aldh1l1, S100b, Sox9

Lineage / identity markers (“this is an astrocyte”), not functional uptake. Nes was placed in `cytoskeleton` (IF), Slc1a3 in `glutamate_ion`.

**Overlap:** none.

## glutamate_ion

**Genes:** Slc1a2, Slc1a3, Kcnj10, Gja1, Glul

Homeostatic function: glutamate uptake (GLT-1/GLAST), Kir4.1, connexin-43, glutamine synthetase. More likely to move with process territory than identity markers.

**Overlap:** none.

## OPC_pdgfra

**Genes:** Pdgfra, Cspg4, Gpr17

OPC / NG2 compartment. Split from `notch_fate` so Pdgfra does not mix Notch ligand dynamics with OPC surface markers.

**Overlap:** conceptually with `notch_fate`.

## ECM_remodeling

**Genes:** Col11a2, Col4a1, Sparc, Fn1, Timp3, Mmp14

Structural remodeling. Col11a2 was flagged in the DE flip-gene scan; the rest are common astrocyte/basement-membrane ECM genes. Complements cytoskeleton (inside the cell) vs matrix (outside).

**Overlap:** none.

---

## What is not a program

Chemokines, MHC class II, A1 vs A2 sets, Wnt, TGF-β, and cholesterol/Apoe were left out unless a later DE scan shows an independent time course. Individual genes are still available in `predicted_panel_genes.csv` (program genes plus top DE from contrasts).
