#!/usr/bin/env python3
"""
FDR Correction + Consistency Tests
Implements: Benjamini-Hochberg FDR, Cochran Q, Interaction z-test
"""
import numpy as np
import pandas as pd
from scipy import stats


def benjamini_hochberg(pvals):
    """BH-FDR correction. Returns q-values."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    fdr_q = np.empty(n)
    cummin = 1.0
    for i in range(n - 1, -1, -1):
        val = min(1.0, ranked[i] * n / (i + 1))
        cummin = min(cummin, val)
        fdr_q[order[i]] = cummin
    return fdr_q


def cochran_q(betas, ses):
    """Cochran Q test for heterogeneity across K estimates."""
    weights = 1.0 / ses**2
    beta_pooled = np.sum(weights * betas) / np.sum(weights)
    Q = np.sum(weights * (betas - beta_pooled)**2)
    df = len(betas) - 1
    p = 1 - stats.chi2.cdf(Q, df) if df > 0 else 1.0
    return {'Q': Q, 'df': df, 'Q_p': p, 'beta_pooled': beta_pooled}


def interaction_z_test(beta1, se1, beta2, se2):
    """Z-test for difference between two MR estimates."""
    diff = beta1 - beta2
    se_diff = np.sqrt(se1**2 + se2**2)
    z = diff / se_diff
    p = 2 * stats.norm.sf(abs(z))
    return {'z': z, 'p': p, 'diff': diff, 'se_diff': se_diff}
