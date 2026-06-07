#!/usr/bin/env python3
"""
Drug Target cis-MR Analysis
Implements: cis-MR using instruments within ±500kb of target gene
"""
import numpy as np
import pandas as pd
from scipy import stats


# Drug target gene coordinates (hg19)
DRUG_TARGETS = {
    'PCSK9': {'chr': 1, 'start': 55505149, 'end': 55530526, 'drug': 'PCSK9 inhibitor'},
    'HMGCR': {'chr': 5, 'start': 74632154, 'end': 74657929, 'drug': 'statin'},
    'NPC1L1': {'chr': 7, 'start': 44552134, 'end': 44580928, 'drug': 'ezetimibe'},
    'CETP': {'chr': 16, 'start': 56995835, 'end': 57017757, 'drug': 'CETP inhibitor'},
}
CIS_WINDOW = 500_000  # ±500kb


def extract_cis_instruments(exp_df, target_name, p_thresh=5e-6):
    """Extract cis-acting instruments within window of target gene."""
    gene = DRUG_TARGETS[target_name]
    mask = (
        (exp_df['CHR'] == gene['chr']) &
        (exp_df['POS'] >= gene['start'] - CIS_WINDOW) &
        (exp_df['POS'] <= gene['end'] + CIS_WINDOW) &
        (exp_df['pval'] < p_thresh)
    )
    return exp_df[mask].copy()


def cis_mr_ivw(beta_exp, se_exp, beta_out, se_out):
    """IVW for cis-MR (often few instruments, fixed effects)."""
    if len(beta_exp) == 0:
        return {'beta': np.nan, 'se': np.nan, 'p': np.nan, 'or': np.nan, 'n_iv': 0}
    if len(beta_exp) == 1:
        beta = beta_out[0] / beta_exp[0]
        se = se_out[0] / abs(beta_exp[0])
        p = 2 * stats.norm.sf(abs(beta / se))
        return {'beta': beta, 'se': se, 'p': p, 'or': np.exp(beta), 'n_iv': 1}
    weights = 1.0 / se_out**2
    beta_iv = beta_out / beta_exp
    beta_ivw = np.sum(weights * beta_iv) / np.sum(weights)
    se_ivw = np.sqrt(1.0 / np.sum(weights))
    p = 2 * stats.norm.sf(abs(beta_ivw / se_ivw))
    return {'beta': beta_ivw, 'se': se_ivw, 'p': p,
            'or': np.exp(beta_ivw), 'n_iv': len(beta_exp)}
