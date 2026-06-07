#!/usr/bin/env python3
"""
Multivariable MR + IOP Mediation Analysis
Implements: MVMR-IVW, Sanderson-Windmeijer conditional F, product-of-coefficients mediation
"""
import numpy as np
import pandas as pd
from scipy import stats


def mvmr_ivw(exposure_betas, exposure_ses, outcome_beta, outcome_se):
    """
    MVMR-IVW estimation for K exposures.
    exposure_betas: array (n_snp, K)
    exposure_ses: array (n_snp, K)
    outcome_beta: array (n_snp,)
    outcome_se: array (n_snp,)
    Returns: dict with betas, ses, pvals, conditional F per exposure
    """
    n, K = exposure_betas.shape
    W = np.diag(1.0 / outcome_se**2)
    X = exposure_betas
    Y = outcome_beta

    XtWX = X.T @ W @ X
    XtWY = X.T @ W @ Y
    try:
        beta_mvmr = np.linalg.solve(XtWX, XtWY)
    except np.linalg.LinAlgError:
        return None

    resid = Y - X @ beta_mvmr
    sigma2 = (resid.T @ W @ resid) / (n - K)
    sigma2 = max(1.0, sigma2)
    cov_beta = sigma2 * np.linalg.inv(XtWX)
    se_mvmr = np.sqrt(np.diag(cov_beta))
    p_mvmr = 2 * stats.norm.sf(np.abs(beta_mvmr / se_mvmr))

    # Sanderson-Windmeijer conditional F
    cond_f = np.zeros(K)
    for k in range(K):
        other_cols = [j for j in range(K) if j != k]
        X_other = X[:, other_cols]
        Q_other = np.eye(n) - X_other @ np.linalg.pinv(X_other.T @ W @ X_other) @ X_other.T @ W
        X_k_resid = Q_other @ X[:, k]
        f_num = (X_k_resid.T @ W @ X_k_resid) / 1.0
        cond_f[k] = f_num / (n - K)

    return {
        'beta': beta_mvmr, 'se': se_mvmr, 'p': p_mvmr,
        'or': np.exp(beta_mvmr), 'cond_f': cond_f, 'n_iv': n
    }


def mediation_product_of_coefficients(a_beta, a_se, b_beta, b_se, c_beta, c_se):
    """
    Two-step mediation: indirect = a * b, SE by delta method, Sobel test.
    a: exposure → mediator (lipid → IOP)
    b: mediator → outcome (IOP → glaucoma), from external source
    c: exposure → outcome (lipid → glaucoma), total effect
    """
    indirect = a_beta * b_beta
    se_indirect = np.sqrt(a_beta**2 * b_se**2 + b_beta**2 * a_se**2)
    z_sobel = indirect / se_indirect if se_indirect > 0 else 0
    p_sobel = 2 * (1 - stats.norm.cdf(abs(z_sobel)))
    med_prop = indirect / c_beta if c_beta != 0 else 0
    return {
        'indirect': indirect, 'indirect_se': se_indirect,
        'sobel_z': z_sobel, 'sobel_p': p_sobel,
        'med_prop': med_prop,
        'a_beta': a_beta, 'b_beta': b_beta, 'c_beta': c_beta,
    }
