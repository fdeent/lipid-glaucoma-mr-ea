#!/usr/bin/env python3
"""
Sensitivity Analyses: F-statistics, Steiger directionality, Leave-one-out, r²<0.001
"""
import numpy as np
import pandas as pd
from scipy import stats


def compute_f_stats(beta_exp, se_exp):
    """Per-SNP F-statistic = (beta/se)²"""
    f_vals = (beta_exp / se_exp) ** 2
    return {
        'F_min': f_vals.min(),
        'F_mean': f_vals.mean(),
        'F_max': f_vals.max(),
        'n_weak': int((f_vals <= 10).sum()),
        'all_F_gt_10': bool((f_vals > 10).all()),
        'n_iv': len(f_vals),
    }


def steiger_test(beta_exp, se_exp, beta_out, se_out, n_exp, n_out):
    """Steiger directionality: confirms exposure→outcome direction."""
    r2_exp = np.sum(beta_exp**2 / (beta_exp**2 + n_exp * se_exp**2))
    r2_out = np.sum(beta_out**2 / (beta_out**2 + n_out * se_out**2))
    correct = r2_exp > r2_out
    diff = r2_exp - r2_out
    se_diff = np.sqrt(4 / n_exp + 4 / n_out) * np.sqrt(r2_exp * r2_out)
    z = diff / se_diff if se_diff > 0 else 0
    p = 2 * stats.norm.sf(abs(z))
    return {'correct_direction': correct, 'r2_exp': r2_exp,
            'r2_out': r2_out, 'steiger_p': p}


def leave_one_out(beta_exp, se_exp, beta_out, se_out):
    """Leave-one-out IVW analysis."""
    n = len(beta_exp)
    results = []
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        weights = 1.0 / se_out[mask]**2
        beta_iv = beta_out[mask] / beta_exp[mask]
        beta_ivw = np.sum(weights * beta_iv) / np.sum(weights)
        se_ivw = np.sqrt(1.0 / np.sum(weights))
        q = np.sum(weights * (beta_iv - beta_ivw)**2)
        q_df = n - 2
        phi = max(1.0, q / q_df)
        se_re = se_ivw * np.sqrt(phi)
        p = 2 * stats.norm.sf(abs(beta_ivw / se_re))
        results.append({'excluded_idx': i, 'ivw_beta': beta_ivw,
                       'ivw_se': se_re, 'ivw_p': p,
                       'ivw_or': np.exp(beta_ivw)})
    return results
