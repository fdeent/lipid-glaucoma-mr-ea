# Output Files

This directory contains the manuscript-reported v7 result tables and figure source data used for result verification. These files allow readers to inspect the numerical results supporting the manuscript without rerunning the full GWAS-scale pipeline.

These outputs are not a substitute for raw-data reproduction: full reruns still require the external GWAS summary statistics and LD reference panel described in `../data/DATA_SOURCES.md`.

## Main manuscript result files

| File | Purpose |
|------|---------|
| `v7_table2_with_cis.csv` | Final Table 2 source: primary, replication, cross-validation, and FinnGen POAG results |
| `v7_main_from_v6.csv` | v7-filtered 84-row univariable MR analysis matrix |
| `v7_mvmr_results.csv` | MVMR conditional effects and Sanderson-Windmeijer conditional F-statistics |
| `v7_mediation_results.csv` | TG-IOP-glaucoma mediation analysis results |
| `v7_f_statistics.csv` | Instrument strength summary for primary instruments |
| `v7_steiger_results.csv` | Steiger directionality summary |
| `v7_steiger_per_snp.csv` | Per-SNP Steiger directionality details |
| `v7_leave_one_out_all.csv` | Leave-one-out sensitivity output |
| `v7_reverse_mr.csv` | Reverse MR results |
| `v7_stricter_ld_results.csv` | Sensitivity analysis using stricter LD clumping threshold |
| `v7_positive_control_results.csv` | CAD/MI positive-control results |
| `v7_atc_s01e_results.csv` | Ophthalmic medication-use proxy analysis |
| `v7_consistency_tests.csv` | Cross-outcome and subtype consistency tests |

## Figure source data

| File | Figure |
|------|--------|
| `v7_fig1_source.csv` | Figure 1 |
| `v7_fig2_mvmr_source.csv` | Figure 2A |
| `v7_fig2_subtype_source.csv` | Figure 2B |
| `v7_figS1_source.csv` | Supplementary Figure S1 |
| `v7_figS2_scatter_source.csv` | Supplementary Figure S2 |
| `v7_figS3_leave_one_out_source.csv` | Supplementary Figure S3 |

## Supplementary table source files

| File | Supplementary table |
|------|---------------------|
| `table_s1_full_mr.csv` | Supplementary Table S1: complete univariable MR results |
| `table_s2_mvmr.csv` | Supplementary Table S2: MVMR results |
| `table_s3_drugtarget.csv` | Supplementary Table S3: drug-target cis-MR results |
| `table_s4_protein_mr.csv` | Supplementary Table S4: protein MR results |
| `table_s5_sensitivity.csv` | Supplementary Table S5: sensitivity analysis details |
| `table_s6_supporting_poag.csv` | Supplementary Table S6: supporting GBMI East Asian POAG analyses |

## Audit files

| File | Purpose |
|------|---------|
| `v7_paragraph_table_audit.csv` | Numeric audit linking manuscript claims to output files |
| `v7_reference_audit.csv` | Reference content and BibTeX metadata audit |
| `v7_taiwan_reference_audit.csv` | TPMI/TWB reference and description audit |
| `v7_output_coverage_audit.csv` | Output coverage audit classifying final/supporting vs legacy/intermediate files |

## Reproduction levels

- Directly runnable in this package: `../demo/run_demo.py` and `../demo/test_all_modules.py`.
- Directly inspectable here: manuscript-reported v7 outputs and figure source data.
- Requires external data: full GWAS-scale rerun from raw summary statistics and LD clumping.
