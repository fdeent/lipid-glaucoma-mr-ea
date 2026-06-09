#!/usr/bin/env python3
"""
Core MR Analysis — Reproducibility Package
Implements: IVW, Weighted Median, MR-Egger, MR-PRESSO, Mode-based estimation
Input: standardized exposure/outcome TSV files
Output: MR results CSV
"""
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import subprocess, tempfile, warnings

warnings.filterwarnings('ignore')


def ivw(beta_exp, se_exp, beta_out, se_out):
    weights = 1.0 / se_out**2
    beta_iv = beta_out / beta_exp
    beta_ivw = np.sum(weights * beta_iv) / np.sum(weights)
    se_ivw = np.sqrt(1.0 / np.sum(weights))
    # Multiplicative random effects
    q = np.sum(weights * (beta_iv - beta_ivw)**2)
    q_df = len(beta_exp) - 1
    phi = max(1.0, q / q_df)
    se_ivw_re = se_ivw * np.sqrt(phi)
    p = 2 * stats.norm.sf(abs(beta_ivw / se_ivw_re))
    q_p = 1 - stats.chi2.cdf(q, q_df) if q_df > 0 else 1.0
    return beta_ivw, se_ivw_re, p, q, q_p


def weighted_median(beta_exp, se_exp, beta_out, se_out, n_boot=1000):
    beta_iv = beta_out / beta_exp
    weights = (se_exp**2 / beta_exp**2 + se_out**2 / beta_out**2)**(-0.5)
    weights /= weights.sum()
    order = np.argsort(beta_iv)
    cum_w = np.cumsum(weights[order])
    idx = np.searchsorted(cum_w, 0.5)
    beta_wm = beta_iv[order[min(idx, len(order)-1)]]
    # Bootstrap SE
    n = len(beta_exp)
    betas_boot = np.empty(n_boot)
    for i in range(n_boot):
        idx_b = np.random.randint(0, n, n)
        b_iv = beta_iv[idx_b]
        w = weights[idx_b]
        w /= w.sum()
        o = np.argsort(b_iv)
        cw = np.cumsum(w[o])
        j = np.searchsorted(cw, 0.5)
        betas_boot[i] = b_iv[o[min(j, len(o)-1)]]
    se_wm = np.std(betas_boot)
    p = 2 * stats.norm.sf(abs(beta_wm / se_wm))
    return beta_wm, se_wm, p


def mr_egger(beta_exp, se_exp, beta_out, se_out):
    weights = 1.0 / se_out**2
    x = beta_exp
    y = beta_out
    w = weights
    wx = w * x
    n = len(x)
    sum_w = w.sum()
    sum_wx = wx.sum()
    sum_wxx = (w * x * x).sum()
    sum_wy = (w * y).sum()
    sum_wxy = (w * x * y).sum()
    denom = sum_w * sum_wxx - sum_wx**2
    beta_egger = (sum_w * sum_wxy - sum_wx * sum_wy) / denom
    intercept = (sum_wxx * sum_wy - sum_wx * sum_wxy) / denom
    y_pred = intercept + beta_egger * x
    resid = y - y_pred
    s2 = np.sum(w * resid**2) / (n - 2)
    se_beta = np.sqrt(s2 * sum_w / denom)
    se_int = np.sqrt(s2 * sum_wxx / denom)
    p_beta = 2 * stats.t.sf(abs(beta_egger / se_beta), n - 2)
    p_int = 2 * stats.t.sf(abs(intercept / se_int), n - 2)
    return beta_egger, se_beta, p_beta, intercept, se_int, p_int


def mr_presso(beta_exp, se_exp, beta_out, se_out, n_sim=1000):
    n = len(beta_exp)
    beta_iv = beta_out / beta_exp
    weights = 1.0 / se_out**2
    beta_ivw = np.sum(weights * beta_iv) / np.sum(weights)
    rss_obs = np.sum(weights * (beta_iv - beta_ivw)**2)
    rss_sim = np.zeros(n_sim)
    for s in range(n_sim):
        b_out_sim = np.random.normal(beta_ivw * beta_exp, se_out)
        b_iv_sim = b_out_sim / beta_exp
        b_ivw_sim = np.sum(weights * b_iv_sim) / np.sum(weights)
        rss_sim[s] = np.sum(weights * (b_iv_sim - b_ivw_sim)**2)
    global_p = np.mean(rss_sim >= rss_obs)
    # Outlier detection
    outliers = []
    if global_p < 0.05 and n > 3:
        for i in range(n):
            mask = np.ones(n, dtype=bool)
            mask[i] = False
            b_iv_i = beta_iv[mask]
            w_i = weights[mask]
            b_ivw_i = np.sum(w_i * b_iv_i) / np.sum(w_i)
            rss_i = weights[i] * (beta_iv[i] - b_ivw_i)**2
            p_i = 1 - stats.chi2.cdf(rss_i, 1)
            if p_i < 0.05 / n:
                outliers.append(i)
    # Corrected estimate
    if outliers:
        mask = np.ones(n, dtype=bool)
        mask[outliers] = False
        b_corr, se_corr, p_corr, _, _ = ivw(
            beta_exp[mask], se_exp[mask], beta_out[mask], se_out[mask])
    else:
        b_corr, se_corr, p_corr = beta_ivw, None, None
    return global_p, outliers, b_corr, se_corr, p_corr


def mode_based(beta_exp, se_exp, beta_out, se_out, n_boot=1000):
    beta_iv = beta_out / beta_exp
    se_iv = np.sqrt(se_out**2 / beta_exp**2 + beta_out**2 * se_exp**2 / beta_exp**4)
    weights = 1.0 / se_iv
    from scipy.stats import gaussian_kde
    try:
        kde = gaussian_kde(beta_iv, weights=weights, bw_method='silverman')
        x_grid = np.linspace(beta_iv.min() - 1, beta_iv.max() + 1, 1000)
        beta_mode = x_grid[np.argmax(kde(x_grid))]
    except Exception:
        beta_mode = np.median(beta_iv)
    n = len(beta_exp)
    betas_boot = np.empty(n_boot)
    for i in range(n_boot):
        idx = np.random.randint(0, n, n)
        b_iv_b = beta_iv[idx]
        w_b = weights[idx]
        try:
            kde_b = gaussian_kde(b_iv_b, weights=w_b, bw_method='silverman')
            betas_boot[i] = x_grid[np.argmax(kde_b(x_grid))]
        except Exception:
            betas_boot[i] = np.median(b_iv_b)
    se_mode = np.std(betas_boot)
    p = 2 * stats.norm.sf(abs(beta_mode / se_mode)) if se_mode > 0 else 1.0
    return beta_mode, se_mode, p


def harmonize(exp_df, out_df):
    merged = exp_df.merge(out_df, on='SNP', suffixes=('_exp', '_out'))
    flip = merged['effect_allele_exp'] == merged['other_allele_out']
    merged.loc[flip, 'beta_out'] *= -1
    ea_match = (merged['effect_allele_exp'] == merged['effect_allele_out']) | flip
    merged = merged[ea_match].copy()
    merged = merged.dropna(subset=['beta_exp', 'se_exp', 'beta_out', 'se_out'])
    merged = merged[(merged['se_exp'] > 0) & (merged['se_out'] > 0)]
    return merged


def run_mr_for_pair(merged, exp_name, out_name):
    if len(merged) < 3:
        return None
    beta_exp = merged['beta_exp'].values
    se_exp = merged['se_exp'].values
    beta_out = merged['beta_out'].values
    se_out = merged['se_out'].values
    # IVW
    b_ivw, se_ivw, p_ivw, q, q_p = ivw(beta_exp, se_exp, beta_out, se_out)
    # Weighted Median
    b_wm, se_wm, p_wm = weighted_median(beta_exp, se_exp, beta_out, se_out)
    # MR-Egger
    b_eg, se_eg, p_eg, intc, se_int, p_int = mr_egger(
        beta_exp, se_exp, beta_out, se_out)
    # MR-PRESSO
    gp, outliers, b_pr, se_pr, p_pr = mr_presso(
        beta_exp, se_exp, beta_out, se_out)
    # Mode-based
    b_md, se_md, p_md = mode_based(beta_exp, se_exp, beta_out, se_out)
    return {
        'exposure': exp_name, 'outcome': out_name, 'n_iv': len(merged),
        'ivw_beta': b_ivw, 'ivw_se': se_ivw, 'ivw_p': p_ivw,
        'ivw_or': np.exp(b_ivw),
        'wm_beta': b_wm, 'wm_se': se_wm, 'wm_p': p_wm,
        'wm_or': np.exp(b_wm),
        'egger_beta': b_eg, 'egger_se': se_eg, 'egger_p': p_eg,
        'egger_or': np.exp(b_eg),
        'egger_intercept': intc, 'egger_int_se': se_int,
        'egger_int_p': p_int,
        'presso_global_p': gp, 'presso_n_outlier': len(outliers),
        'presso_beta_corr': b_pr, 'presso_p_corr': p_pr,
        'mode_beta': b_md, 'mode_se': se_md, 'mode_p': p_md,
        'mode_or': np.exp(b_md),
        'q': q, 'q_p': q_p,
    }
