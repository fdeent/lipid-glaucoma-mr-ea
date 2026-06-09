# Reproducibility Package — Lipid-Glaucoma MR Study

## Overview

This package provides the analysis code and a runnable demonstration for the Mendelian randomization study of genetically proxied lipid levels and glaucoma risk in East Asian populations.

It is a reproducibility code package, not a self-contained archive of all GWAS summary statistics. The included demo data reproduce the core finding using a small real-data subset. The `output/` directory includes the manuscript-reported v7 result CSVs and figure source data for numerical verification. Full raw-data reproduction requires downloading the external GWAS summary statistics and LD reference panel listed in `data/DATA_SOURCES.md`.

## Key figures

### Cross-validation of lipid-glaucoma associations

![Figure 1: cross-validation forest plot](figures/fig1_forest_cross_validation.png)

### Independent lipid effects and subtype specificity

![Figure 2: MVMR and subtype specificity](figures/fig2_mvmr_subtype.png)

## Quick start

```bash
pip install -r requirements.txt
cd demo
python test_all_modules.py
python run_demo.py
```

Expected test result: `26 passed, 0 failed`

## Requirements

- Python 3.8+
- numpy >= 1.21, pandas >= 1.4, scipy >= 1.7
- For full analysis: PLINK 1.9 + 1000 Genomes Phase 3 EAS reference panel

---

## Testing: Two Modes

### Mode 1: Probe Test (~10 seconds)

Verifies functional correctness of all 6 analysis modules.

```bash
cd demo
python test_all_modules.py
```

**Coverage** (26 checks across 6 modules):
- Module 1: Core MR — IVW (OR direction + significance), WM (direction), Egger (intercept null + direction), PRESSO (execution)
- Module 2: MVMR — MVMR-IVW (effect direction), SW conditional F (>0), mediation (Sobel execution + proportion <10%)
- Module 3: Sensitivity — F-statistic (min>10, n_weak=0), Steiger (correct direction), LOO (all significant)
- Module 4: FDR — BH-FDR (monotonicity + first value), Cochran Q (homogeneity), Interaction z (significance)
- Module 5: Drug targets — cis-MR framework (4 targets defined), single-IV Wald ratio
- Module 6: Protein MR — pQTL IVW (execution + direction), Egger sensitivity (execution)

**Expected output**: `26 passed, 0 failed`

### Mode 2: Core Finding Reproduction (~15 seconds)

Reproduces the primary finding using a real GWAS data subset (KOGES total lipids to BBJ Glaucoma).

```bash
cd demo
python run_demo.py
```

**Real data used**: 275 genome-wide significant SNPs from KoGES total lipids GWAS, 9 pre-clumped IVs (r² < 0.01, 10 Mb window, 1000G EAS reference panel).

**Expected output**: IVW OR ≈ 0.87, p < 0.001; Weighted Median OR ≈ 0.85; MR-Egger OR ≈ 0.87, intercept p > 0.05. Minor OR differences from manuscript Table 2 arise because 8 of 9 IVs successfully harmonize in the demo subset.

---

## Full Reproduction

Full reproduction requires external GWAS summary statistics and the LD reference panel. These large files are not redistributed.

### Step 1: Download GWAS data

See `data/DATA_SOURCES.md` for download URLs and format specifications.

### Step 2: Standardize to 12-column TSV

```
SNP  effect_allele  other_allele  beta  se  pval  eaf  n  trait  source  CHR  POS
```

### Step 3: Configure paths

Edit path constants in the scripts under `code/` to point to your downloaded files.

### Step 4: Run modules

```bash
cd code
python 01_run_mr.py            # Core MR: IVW, WM, MR-Egger, MR-PRESSO, mode-based
python 02_run_mvmr.py          # MVMR-IVW and mediation
python 03_sensitivity.py       # F-statistics, Steiger, leave-one-out
python 04_fdr_consistency.py   # BH-FDR, Cochran Q, interaction z-test
python 05_drug_target_mr.py    # Drug-target cis-MR
python 06_protein_mr.py        # Protein pQTL MR
```

---

## Analysis Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| LD clumping r² | 0.01 | Analysis scripts |
| LD clumping window | 10,000 kb (10 Mb) | Analysis scripts |
| Significance threshold (p1) | 5 x 10^-8 | Analysis scripts |
| Secondary threshold (p2) | 1 x 10^-5 | `01_run_mr.py` |
| Minimum IVs | 3 | `01_run_mr.py` |
| F-statistic threshold | 10 | `03_sensitivity.py` |
| WM bootstrap replicates | 1,000 | `01_run_mr.py` |
| Mode bootstrap replicates | 1,000 | `01_run_mr.py` |
| PRESSO simulations | 1,000 | `01_run_mr.py` |
| IOP-Glaucoma b-path (Khawaja 2018) | beta = 0.18, SE = 0.02 | `02_run_mvmr.py` |
| Drug target cis-window | +/- 500 kb | `05_drug_target_mr.py` |
| LD reference | 1000G Phase 3 EAS | PLINK bim/bed/fam |
| SNP matching | rsID-based | Cross-build harmonization |

---

## Code Modules

| File | Function |
|------|----------|
| `01_run_mr.py` | IVW, Weighted Median, MR-Egger, MR-PRESSO, Mode-based, harmonize |
| `02_run_mvmr.py` | MVMR-IVW, SW conditional F, mediation (product-of-coefficients + Sobel) |
| `03_sensitivity.py` | F-statistics, Steiger directionality, Leave-one-out |
| `04_fdr_consistency.py` | Benjamini-Hochberg FDR, Cochran Q, Interaction z-test |
| `05_drug_target_mr.py` | cis-MR for 4 drug targets (gene coordinates + cis-window extraction) |
| `06_protein_mr.py` | pQTL MR-IVW + WM/Egger sensitivity |

---

## Directory Structure

```
.
├── README.md                 This file
├── requirements.txt          Python dependencies
├── code/                     Core analysis scripts (6 modules)
├── data/                     Data source documentation (GWAS downloaded separately)
├── demo/                     Minimal reproducible example with real data subset
│   ├── run_demo.py           Reproduce core finding (KOGES LIP → BBJ Glaucoma)
│   ├── test_all_modules.py   Probe test (26 checks across all modules)
│   ├── run_mr_core.py        Importable MR function module
│   ├── exposure_koges_LIP.tsv  Real exposure data (275 significant SNPs)
│   └── outcome_bbj_Gla.tsv    Real outcome data (207 matched SNPs)
├── figures/                  Rendered manuscript and supplementary figures
└── output/                   Manuscript-reported result CSVs and figure source data
```

## Citation

If you use this code, please cite our manuscript.
