# Demo: KOGES Total Lipids → BBJ Glaucoma

This demo reproduces the core finding of the manuscript using a minimal real data subset.

## What it does

Runs a two-sample Mendelian randomization analysis using:
- **Exposure**: KOGES total lipids (9 genome-wide significant, LD-independent instruments)
- **Outcome**: BBJ Glaucoma (5,765 cases / 171,586 controls)

## Expected result

| Method | OR | p-value |
|--------|-----|---------|
| IVW | 0.87 | < 0.001 |
| Weighted Median | 0.85 | < 0.001 |
| MR-Egger | 0.87 | 0.01 |
| MR-Egger intercept | — | 0.62 (null) |

## How to run

### Probe test (all modules)

```bash
python test_all_modules.py
```

Runs 26 automated checks across all 6 analysis modules. Expected: `26 passed, 0 failed`.

### Core finding reproduction

```bash
python run_demo.py
```

Reproduces KOGES total lipids → BBJ Glaucoma (manuscript Table 2, row 2).

## Files

| File | Description |
|------|-------------|
| `exposure_koges_LIP.tsv` | 275 genome-wide significant SNPs from KoGES total lipids GWAS |
| `outcome_bbj_Gla.tsv` | Corresponding SNP associations from BBJ Glaucoma GWAS |
| `run_mr_core.py` | Core MR functions (IVW, Weighted Median, MR-Egger, MR-PRESSO, harmonize) |
| `run_demo.py` | Demo script using pre-defined clumped instrument set |
| `test_all_modules.py` | Probe test covering all 6 code modules |

## Note on clumping

The demo uses a pre-defined list of 9 LD-independent SNPs (r² < 0.01, 10 Mb window) obtained from the full analysis. Running PLINK LD clumping on arbitrary data requires the 1000 Genomes EAS reference panel. See `../data/DATA_SOURCES.md` for download instructions.
