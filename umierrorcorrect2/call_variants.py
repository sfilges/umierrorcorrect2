#!/usr/bin/env python3
from argparse import Namespace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.stats import betabinom

from umierrorcorrect2.core.utils import get_sample_name_from_cons, parse_cons_file

# Default beta-binomial model parameters (alpha, beta)
# These were fitted from background error data
DEFAULT_BETABINOM_ALPHA = 2.168215069116764
DEFAULT_BETABINOM_BETA = 3531.588541594945


def write_vcf(
    vcffile: str | Path,
    rout: NDArray[np.str_],
    qsig: NDArray[np.floating] | list[str],
    reference: str,
) -> None:
    """Write variant calls to VCF format."""
    vcf_header = f"""\
##fileformat=VCFv4.2
##reference={reference}
##source=umierrorcorrectV0.1
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AD,Number=1,Type=Float,Description="Alternative Allele Depth">
##INFO=<ID=AF,Number=A,Type=Integer,Description="Alternative Allele Frequency">
##FILTER=<ID=a5,Description="Alternative Allele Depth below 5">
##FILTER=<ID=q10,Description="Variant quality below 10">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Allele Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""
    with Path(vcffile).open("w") as g:
        g.write(vcf_header)

        for record, qual in zip(rout, qsig):
            qual_str = "." if qual == "NA" else str(qual)
            parts = record.split("\t")

            filters = []
            if qual_str != "." and float(qual_str) < 10:
                filters.append("q10")
            if int(parts[-3]) < 5:
                filters.append("a5")
            filter_str = "PASS" if not filters else ",".join(filters)

            vcf_line = "\t".join(
                [
                    parts[1],  # CHROM
                    parts[2],  # POS
                    ".",  # ID
                    parts[4],  # REF
                    parts[-1],  # ALT
                    qual_str,  # QUAL
                    filter_str,  # FILTER
                    f"DP={parts[-5]};AD={parts[-3]};AF={parts[-2]}",  # INFO
                    "DP",  # FORMAT
                    parts[-3],  # SAMPLE
                ]
            )
            g.write(vcf_line + "\n")


def plot_histogram(hist: NDArray[np.floating], plot_filename: str | Path) -> None:
    """Plot histogram of Q-scores and save to file."""
    fig, ax = plt.subplots()
    ax.hist(hist, bins=100, facecolor="dodgerblue", alpha=0.5)
    ax.set_xlabel("Q-score")
    ax.set_ylabel("Frequency")
    ax.set_title("Histogram of Q-scores")
    ax.set_frame_on(False)
    ax.set_xlim(0, 140)
    fig.savefig(plot_filename)
    plt.close(fig)


def _load_betabinom_params(params_file: Path | None) -> tuple[float, float]:
    """Load beta-binomial parameters from file or return defaults."""
    if params_file:
        with Path(params_file).open() as f:
            params = [float(line.rstrip()) for line in f]
            return params[0], params[1]
    return DEFAULT_BETABINOM_ALPHA, DEFAULT_BETABINOM_BETA


def _calculate_qscores(
    n1: NDArray[np.int_],
    a1: NDArray[np.int_],
    alpha: float,
    beta: float,
) -> NDArray[np.floating]:
    """Calculate Q-scores using beta-binomial model."""
    pvalues = 1 - betabinom.cdf(a1 - 1, n1, alpha, beta)
    # Clamp extreme values to avoid log(0) issues
    pvalues = np.clip(pvalues, 1e-10, None)
    pvalues[np.isnan(pvalues)] = 1e-10
    return -10 * np.log10(pvalues)


def run_call_variants(args: Namespace) -> Path:
    """Run variant calling on consensus sequences.

    Args:
        args: Namespace with required fields:
            - cons_file: Path to consensus file (or None to auto-detect)
            - output_path: Output directory
            - sample_name: Sample name (or None to auto-detect)
            - fsize: Family size cutoff
            - vc_method: 'count' or 'bbmodel'
            - count_cutoff: Minimum alt count for count method
            - qvalue_threshold: Q-score threshold for bbmodel method
            - params_file: Path to beta-binomial params file (optional)
            - reference_file: Path to reference genome
            - plot_qscore_histogram: Whether to generate histogram (optional)

    Returns:
        Path to the output VCF file.
    """
    # Auto-detect cons_file if not provided
    if not args.cons_file:
        cons_files = list(Path(args.output_path).glob("*cons.tsv"))
        if not cons_files:
            raise FileNotFoundError(f"No consensus file found in {args.output_path}")
        args.cons_file = str(cons_files[0])

    if not args.sample_name:
        args.sample_name = get_sample_name_from_cons(args.cons_file)

    # Parse consensus file
    fsize = int(args.fsize)
    _f1, n1, a1, data = parse_cons_file(args.cons_file, fsize)
    n1 = np.array(n1)
    a1 = np.array(a1)
    data = np.array(data)

    method = args.vc_method.lower()
    count_cutoff = float(args.count_cutoff)

    if method == "count":
        # Simple count-based filtering
        mask = a1 >= count_cutoff
        rout = data[mask]
        qsig: NDArray[np.floating] | list[str] = ["NA"] * len(rout)
    else:
        # Beta-binomial model
        alpha, beta = _load_betabinom_params(getattr(args, "params_file", None))
        qscores = _calculate_qscores(n1, a1, alpha, beta)

        # Generate histogram if requested
        if getattr(args, "plot_qscore_histogram", False):
            histogram_path = Path(args.output_path) / f"{args.sample_name}.qscore_histogram.png"
            plot_histogram(qscores, histogram_path)

        mask = qscores >= float(args.qvalue_threshold) if method == "bbmodel" else a1 >= count_cutoff

        rout = data[mask]
        qsig = qscores[mask]

    outfilename = Path(args.output_path) / f"{args.sample_name}.vcf"
    write_vcf(outfilename, rout, qsig, args.reference_file)
    return outfilename
