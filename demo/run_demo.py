#!/usr/bin/env python3
"""
Demo: Reproduce KOGES Total Lipids → BBJ Glaucoma (OR≈0.86, p≈2.9e-4)
Uses pre-extracted genome-wide significant SNPs (no PLINK clumping needed).
"""
import sys
sys.path.insert(0, '../code')
import numpy as np
import pandas as pd
from run_mr_core import ivw, weighted_median, mr_egger, mr_presso, harmonize

np.random.seed(42)

CLUMPED_SNPS = [
    'rs651821', 'rs1065853', 'rs2980860', 'rs635634', 'rs2738464',
    'rs629301', 'rs79490663', 'rs2278426', 'rs1260326'
]

def main():
    print("=" * 60)
    print("DEMO: KOGES Total Lipids → BBJ Glaucoma")
    print("=" * 60)

    exp = pd.read_csv('exposure_koges_LIP.tsv', sep='\t')
    out = pd.read_csv('outcome_bbj_Gla.tsv', sep='\t')

    exp_iv = exp[exp['SNP'].isin(CLUMPED_SNPS)]
    merged = harmonize(exp_iv, out)
    n_iv = len(merged)
    print(f"\nInstrumental variables after harmonization: {n_iv}")

    if n_iv < 3:
        print("ERROR: Too few IVs. Check data files.")
        return

    beta_exp = merged['beta_exp'].values
    se_exp = merged['se_exp'].values
    beta_out = merged['beta_out'].values
    se_out = merged['se_out'].values

    # F-statistics
    f_stats = (beta_exp / se_exp) ** 2
    print(f"F-statistics: min={f_stats.min():.1f}, mean={f_stats.mean():.1f}")

    # IVW
    b, se, p, q, qp = ivw(beta_exp, se_exp, beta_out, se_out)
    print(f"\n--- IVW (primary method) ---")
    print(f"  OR = {np.exp(b):.4f} (95% CI: {np.exp(b-1.96*se):.4f}"
          f"–{np.exp(b+1.96*se):.4f})")
    print(f"  p  = {p:.2e}")
    print(f"  Cochran Q p = {qp:.3f}")

    # Weighted Median
    bwm, sewm, pwm = weighted_median(beta_exp, se_exp, beta_out, se_out)
    print(f"\n--- Weighted Median ---")
    print(f"  OR = {np.exp(bwm):.4f}, p = {pwm:.4f}")

    # MR-Egger
    beg, seeg, peg, intc, sei, pi = mr_egger(
        beta_exp, se_exp, beta_out, se_out)
    print(f"\n--- MR-Egger ---")
    print(f"  OR = {np.exp(beg):.4f}, p = {peg:.4f}")
    print(f"  Intercept = {intc:.4f}, p = {pi:.4f}")

    # MR-PRESSO
    gp, outliers, bc, sec, pc = mr_presso(
        beta_exp, se_exp, beta_out, se_out)
    print(f"\n--- MR-PRESSO ---")
    print(f"  Global p = {gp:.3f}, outliers = {len(outliers)}")

    # Expected result
    print(f"\n{'=' * 60}")
    print(f"Expected: OR ≈ 0.86, p ≈ 2.9e-4")
    print(f"Got:      OR = {np.exp(b):.3f}, p = {p:.2e}")
    tol = abs(np.exp(b) - 0.86) < 0.02 and p < 0.001
    print(f"Status:   {'PASS' if tol else 'CHECK'}")


if __name__ == '__main__':
    main()
