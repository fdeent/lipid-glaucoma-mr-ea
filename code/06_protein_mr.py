#!/usr/bin/env python3
"""
Protein MR (pQTL → Glaucoma)
Implements: cis-pQTL MR using BBJ/deCODE protein GWAS instruments
"""
import numpy as np
import pandas as pd
from scipy import stats


def protein_mr_ivw(beta_exp, se_exp, beta_out, se_out):
    """IVW for protein MR (typically few cis-instruments)."""
    n = len(beta_exp)
    if n == 0:
        return None
    if n == 1:
        beta = beta_out[0] / beta_exp[0]
        se = abs(se_out[0] / beta_exp[0])
        p = 2 * stats.norm.sf(abs(beta / se))
        return {'ivw_beta': beta, 'ivw_se': se, 'ivw_p': p,
                'OR': np.exp(beta), 'n_ivs': 1}
    weights = 1.0 / se_out**2
    beta_iv = beta_out / beta_exp
    beta_ivw = np.sum(weights * beta_iv) / np.sum(weights)
    se_ivw = np.sqrt(1.0 / np.sum(weights))
    q = np.sum(weights * (beta_iv - beta_ivw)**2)
    q_df = n - 1
    phi = max(1.0, q / q_df)
    se_re = se_ivw * np.sqrt(phi)
    p = 2 * stats.norm.sf(abs(beta_ivw / se_re))
    q_p = 1 - stats.chi2.cdf(q, q_df)
    i2 = max(0, (q - q_df) / q) if q > 0 else 0
    return {'ivw_beta': beta_ivw, 'ivw_se': se_re, 'ivw_p': p,
            'OR': np.exp(beta_ivw),
            'OR_lo': np.exp(beta_ivw - 1.96 * se_re),
            'OR_hi': np.exp(beta_ivw + 1.96 * se_re),
            'Q': q, 'Q_df': q_df, 'Q_p': q_p, 'I2': i2, 'n_ivs': n}


def protein_mr_sensitivity(beta_exp, se_exp, beta_out, se_out):
    """Run WM + Egger if n_iv >= 3."""
    n = len(beta_exp)
    result = {}
    if n < 3:
        return result
    # Weighted median
    beta_iv = beta_out / beta_exp
    w = (se_exp**2 / beta_exp**2 + se_out**2 / beta_out**2)**(-0.5)
    w /= w.sum()
    order = np.argsort(beta_iv)
    cw = np.cumsum(w[order])
    idx = np.searchsorted(cw, 0.5)
    result['wm_beta'] = beta_iv[order[min(idx, len(order)-1)]]
    # MR-Egger
    W = np.diag(1.0 / se_out**2)
    x, y = beta_exp, beta_out
    sw = (1/se_out**2).sum(); swx = (x/se_out**2).sum()
    swxx = (x*x/se_out**2).sum(); swy = (y/se_out**2).sum()
    swxy = (x*y/se_out**2).sum()
    denom = sw * swxx - swx**2
    result['egger_beta'] = (sw * swxy - swx * swy) / denom
    intercept = (swxx * swy - swx * swxy) / denom
    resid = y - (intercept + result['egger_beta'] * x)
    s2 = np.sum(resid**2 / se_out**2) / (n - 2)
    result['egger_intercept'] = intercept
    se_int = np.sqrt(s2 * swxx / denom)
    result['egger_intercept_p'] = 2 * stats.t.sf(abs(intercept/se_int), n-2)
    return result
