#!/usr/bin/env python3
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import beta

from umierrorcorrect2.core.utils import parse_cons_file


def beta_nll(params: tuple[float, float], data: np.ndarray) -> float:
    """Compute negative log-likelihood for beta distribution.

    Args:
        params: Tuple of (alpha, beta) parameters for the beta distribution.
        data: Array of observed values in the range [0, 1].

    Returns:
        Negative log-likelihood value (inf if parameters are invalid).
    """
    a, b = params
    if a <= 0 or b <= 0:
        return np.inf
    log_pdf = beta.logpdf(data, a, b)
    mask = np.isfinite(log_pdf)
    return -log_pdf[mask].sum()


def estimate_beta_parameters(data: np.ndarray) -> tuple[float, float]:
    """Estimate beta distribution parameters using maximum likelihood.

    Uses method of moments for initial parameter estimates, then optimizes
    using Nelder-Mead to find MLE parameters.

    Args:
        data: Array of observed values in the range [0, 1].

    Returns:
        Tuple of (alpha, beta) parameters for the fitted beta distribution.
    """
    m = np.mean(data)
    v = np.var(data)
    # Method of moments for initial guess
    a0 = m * (m * (1 - m) / v - 1)
    b0 = (1 - m) * (m * (1 - m) / v - 1)

    result = minimize(
        beta_nll,
        x0=[max(a0, 0.1), max(b0, 0.1)],
        args=(data,),
        method="Nelder-Mead",
    )
    return result.x[0], result.x[1]


def run_fit_bgmodel(args):
    """Fit beta distribution to model background sequencing noise.

    Reads consensus data and fits a beta distribution to estimate background
    error rates. Positions with known true mutations should be excluded to
    avoid skewing the background estimate.

    Args:
        args: Namespace with attributes:
            - known_mutations_file: Optional path to file listing positions with
              known true mutations (one per line) to exclude from fitting
            - cons_file: Path to consensus file
            - fsize: Family size threshold
            - out_file: Output path for fitted parameters
    """
    if args.known_mutations_file:
        known_mutation_positions = Path(args.known_mutations_file).read_text().splitlines()
    else:
        known_mutation_positions = []

    if not args.cons_file:
        raise ValueError("cons_file is required")

    fsize = int(args.fsize)
    f1, _n1, _a1, pos, _data = parse_cons_file(args.cons_file, fsize, include_position=True)

    pos = np.array(pos)
    f1 = np.array(f1)

    # Exclude known mutation positions from background model fitting
    is_background = ~np.isin(pos, known_mutation_positions)
    alpha, beta_param = estimate_beta_parameters(f1[is_background])

    with Path(args.out_file).open("w") as g:
        g.write(f"{alpha}\n{beta_param}\n")
