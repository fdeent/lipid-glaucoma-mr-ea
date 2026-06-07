"""
Core MR functions - importable module for demo and full analysis.
"""
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


def ivw(beta_exp, se_exp, beta_out, se_out):
    weights = 1.0 / se_out**2
    beta_iv = beta_out / beta_exp
    beta_ivw = np.sum(weights * beta_iv) / np.sum(weights)
    se_ivw = np.sqrt(1.0 / np.sum(weights))
    q = np.sum(weights * (beta_iv - beta_ivw)**2)
    q_df = len(beta_exp) - 1
    phi = max(1.0, q / q_df)
    se_ivw_re = se_ivw * np.sqrt(phi)
    p = 2 * stats.norm.sf(abs(beta_ivw / se_ivw_re))
    q_p = 1 - stats.chi2.cdf(q, q_df) if q_df > 0 else 1.0
    return beta_ivw, se_ivw_re, p, q, q_p


def weighted_median(beta_exp, se_exp, beta_out, se_out, n_boot=1000):
    beta_iv = beta_out / beta_exp
    w = (se_exp**2 / beta_exp**2 + se_out**2 / beta_out**2)**(-0.5)
    w /= w.sum()
    order = np.argsort(beta_iv)
    cum_w = np.cumsum(w[order])
    idx = np.searchsorted(cum_w, 0.5)
    beta_wm = beta_iv[order[min(idx, len(order)-1)]]
    n = len(beta_exp)
    betas_boot = np.empty(n_boot)
    for i in range(n_boot):
        idx_b = np.random.randint(0, n, n)
        b_iv = beta_iv[idx_b]; wb = w[idx_b]; wb /= wb.sum()
        o = np.argsort(b_iv); cw = np.cumsum(wb[o])
        j = np.searchsorted(cw, 0.5)
        betas_boot[i] = b_iv[o[min(j, len(o)-1)]]
    se_wm = np.std(betas_boot)
    p = 2 * stats.norm.sf(abs(beta_wm / se_wm))
    return beta_wm, se_wm, p


def mr_egger(beta_exp, se_exp, beta_out, se_out):
    w = 1.0 / se_out**2
    x, y, n = beta_exp, beta_out, len(beta_exp)
    sw = w.sum(); swx = (w*x).sum(); swxx = (w*x*x).sum()
    swy = (w*y).sum(); swxy = (w*x*y).sum()
    denom = sw * swxx - swx**2
    beta_eg = (sw * swxy - swx * swy) / denom
    intercept = (swxx * swy - swx * swxy) / denom
    resid = y - (intercept + beta_eg * x)
    s2 = np.sum(w * resid**2) / (n - 2)
    se_beta = np.sqrt(s2 * sw / denom)
    se_int = np.sqrt(s2 * swxx / denom)
    p_beta = 2 * stats.t.sf(abs(beta_eg / se_beta), n - 2)
    p_int = 2 * stats.t.sf(abs(intercept / se_int), n - 2)
    return beta_eg, se_beta, p_beta, intercept, se_int, p_int


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
    outliers = []
    if global_p < 0.05 and n > 3:
        for i in range(n):
            mask = np.ones(n, dtype=bool); mask[i] = False
            b_ivw_i = np.sum(weights[mask]*beta_iv[mask]) / weights[mask].sum()
            rss_i = weights[i] * (beta_iv[i] - b_ivw_i)**2
            p_i = 1 - stats.chi2.cdf(rss_i, 1)
            if p_i < 0.05 / n:
                outliers.append(i)
    return global_p, outliers, None, None, None


def harmonize(exp_df, out_df):
    merged = exp_df.merge(out_df, on='SNP', suffixes=('_exp', '_out'))
    flip = merged['effect_allele_exp'] == merged['other_allele_out']
    merged.loc[flip, 'beta_out'] *= -1
    ea_match = (merged['effect_allele_exp'] == merged['effect_allele_out']) | flip
    merged = merged[ea_match].copy()
    merged = merged.dropna(subset=['beta_exp','se_exp','beta_out','se_out'])
    merged = merged[(merged['se_exp'] > 0) & (merged['se_out'] > 0)]
    return merged
