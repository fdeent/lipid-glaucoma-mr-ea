#!/usr/bin/env python3
"""
Probe Test — Verify all reproducibility modules produce correct results.
Uses demo data (KOGES_LIP → BBJ_Glaucoma) to test every analysis module.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code'))
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import pandas as pd
from run_mr_core import ivw, weighted_median, mr_egger, mr_presso, harmonize

np.random.seed(42)
PASS, FAIL = 0, 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  [PASS] {name}")
    else:
        FAIL += 1
        print(f"  [FAIL] {name} — {detail}")

def load_demo_data():
    demo_dir = os.path.dirname(__file__)
    exp = pd.read_csv(os.path.join(demo_dir, 'exposure_koges_LIP.tsv'), sep='\t')
    out = pd.read_csv(os.path.join(demo_dir, 'outcome_bbj_Gla.tsv'), sep='\t')
    clumped = ['rs651821','rs1065853','rs2980860','rs635634','rs2738464',
               'rs629301','rs79490663','rs2278426','rs1260326']
    exp_iv = exp[exp['SNP'].isin(clumped)]
    merged = harmonize(exp_iv, out)
    return merged

print("=" * 60)
print("PROBE TEST: All Reproducibility Modules")
print("=" * 60)

merged = load_demo_data()
beta_exp = merged['beta_exp'].values
se_exp = merged['se_exp'].values
beta_out = merged['beta_out'].values
se_out = merged['se_out'].values
n_iv = len(merged)

# --- Module 1: Core MR ---
print("\n[1] Core MR (01_run_mr)")
b, se, p, q, qp = ivw(beta_exp, se_exp, beta_out, se_out)
check("IVW OR in range", 0.80 < np.exp(b) < 0.92, f"OR={np.exp(b):.3f}")
check("IVW p < 0.001", p < 0.001, f"p={p:.2e}")

bwm, sewm, pwm = weighted_median(beta_exp, se_exp, beta_out, se_out)
check("WM OR direction", np.exp(bwm) < 1.0, f"OR={np.exp(bwm):.3f}")

beg, seeg, peg, intc, sei, pi = mr_egger(beta_exp, se_exp, beta_out, se_out)
check("Egger intercept null", pi > 0.05, f"p_int={pi:.3f}")
check("Egger OR direction", np.exp(beg) < 1.0, f"OR={np.exp(beg):.3f}")

gp, outliers, _, _, _ = mr_presso(beta_exp, se_exp, beta_out, se_out)
check("PRESSO runs", True)

# --- Module 2: MVMR ---
print("\n[2] MVMR (02_run_mvmr)")
from importlib.util import spec_from_file_location, module_from_spec
spec = spec_from_file_location("mvmr", os.path.join(
    os.path.dirname(__file__), '..', 'code', '02_run_mvmr.py'))
mvmr_mod = module_from_spec(spec); spec.loader.exec_module(mvmr_mod)

# Simulate 3-exposure MVMR with synthetic data
rng = np.random.default_rng(123)
n_snp = 50
exp_betas = rng.normal(0, 0.1, (n_snp, 3))
exp_ses = np.full((n_snp, 3), 0.02)
true_effects = np.array([-0.2, -0.15, 0.0])
out_beta = exp_betas @ true_effects + rng.normal(0, 0.05, n_snp)
out_se = np.full(n_snp, 0.05)

res = mvmr_mod.mvmr_ivw(exp_betas, exp_ses, out_beta, out_se)
check("MVMR returns result", res is not None)
check("MVMR beta[0] negative", res['beta'][0] < 0, f"b={res['beta'][0]:.3f}")
check("MVMR cond_F > 0", all(res['cond_f'] > 0), f"F={res['cond_f']}")

# Mediation
med = mvmr_mod.mediation_product_of_coefficients(
    a_beta=-0.034, a_se=0.024, b_beta=0.18, b_se=0.02,
    c_beta=-0.15, c_se=0.04)
check("Mediation Sobel runs", 'sobel_p' in med)
check("Mediation prop < 10%", abs(med['med_prop']) < 0.10,
      f"prop={med['med_prop']:.3f}")

# --- Module 3: Sensitivity ---
print("\n[3] Sensitivity (03_sensitivity)")
spec3 = spec_from_file_location("sens", os.path.join(
    os.path.dirname(__file__), '..', 'code', '03_sensitivity.py'))
sens_mod = module_from_spec(spec3); spec3.loader.exec_module(sens_mod)

fstats = sens_mod.compute_f_stats(beta_exp, se_exp)
check("F-stats min > 10", fstats['F_min'] > 10, f"F_min={fstats['F_min']:.1f}")
check("F-stats n_weak = 0", fstats['n_weak'] == 0)

steiger = sens_mod.steiger_test(beta_exp, se_exp, beta_out, se_out, 72298, 177351)
check("Steiger correct direction", steiger['correct_direction'])

loo = sens_mod.leave_one_out(beta_exp, se_exp, beta_out, se_out)
all_sig = all(r['ivw_p'] < 0.05 for r in loo)
check(f"LOO all {len(loo)}/{len(loo)} significant", all_sig)

# --- Module 4: FDR + Consistency ---
print("\n[4] FDR + Consistency (04_fdr_consistency)")
spec4 = spec_from_file_location("fdr", os.path.join(
    os.path.dirname(__file__), '..', 'code', '04_fdr_consistency.py'))
fdr_mod = module_from_spec(spec4); spec4.loader.exec_module(fdr_mod)

pvals = np.array([0.001, 0.01, 0.03, 0.05, 0.10, 0.50])
qvals = fdr_mod.benjamini_hochberg(pvals)
check("BH-FDR monotone", all(qvals[i] <= qvals[i+1] for i in range(len(qvals)-1)))
check("BH-FDR first q < 0.01", qvals[0] < 0.01, f"q[0]={qvals[0]:.4f}")

cq = fdr_mod.cochran_q(np.array([-0.15, -0.16, -0.14]), np.array([0.04, 0.05, 0.03]))
check("Cochran Q runs", 'Q_p' in cq)
check("Cochran Q p > 0.5 (homogeneous)", cq['Q_p'] > 0.5, f"Qp={cq['Q_p']:.3f}")

iz = fdr_mod.interaction_z_test(-0.15, 0.04, 0.18, 0.05)
check("Interaction z significant", iz['p'] < 0.001)

# --- Module 5: Drug Target ---
print("\n[5] Drug Target cis-MR (05_drug_target_mr)")
spec5 = spec_from_file_location("dt", os.path.join(
    os.path.dirname(__file__), '..', 'code', '05_drug_target_mr.py'))
dt_mod = module_from_spec(spec5); spec5.loader.exec_module(dt_mod)

check("DRUG_TARGETS has 4 entries", len(dt_mod.DRUG_TARGETS) == 4)
res_cis = dt_mod.cis_mr_ivw(
    np.array([0.3]), np.array([0.05]), np.array([-0.1]), np.array([0.08]))
check("cis-MR single IV runs", res_cis['n_iv'] == 1)
check("cis-MR OR < 1", res_cis['or'] < 1.0)

# --- Module 6: Protein MR ---
print("\n[6] Protein MR (06_protein_mr)")
spec6 = spec_from_file_location("prot", os.path.join(
    os.path.dirname(__file__), '..', 'code', '06_protein_mr.py'))
prot_mod = module_from_spec(spec6); spec6.loader.exec_module(prot_mod)

res_prot = prot_mod.protein_mr_ivw(
    np.array([0.5, 0.3, 0.4]), np.array([0.05, 0.04, 0.06]),
    np.array([-0.08, -0.05, -0.07]), np.array([0.03, 0.04, 0.03]))
check("Protein MR runs", res_prot is not None)
check("Protein MR OR < 1", res_prot['OR'] < 1.0)
sens_p = prot_mod.protein_mr_sensitivity(
    np.array([0.5, 0.3, 0.4]), np.array([0.05, 0.04, 0.06]),
    np.array([-0.08, -0.05, -0.07]), np.array([0.03, 0.04, 0.03]))
check("Protein sensitivity has Egger", 'egger_beta' in sens_p)

# --- Summary ---
print(f"\n{'=' * 60}")
print(f"RESULTS: {PASS} passed, {FAIL} failed out of {PASS+FAIL} tests")
print(f"{'=' * 60}")
sys.exit(0 if FAIL == 0 else 1)
