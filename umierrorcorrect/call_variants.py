#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy import inf
from scipy.stats import betabinom

from umierrorcorrect.core.utils import get_sample_name_from_cons, parse_cons_file


def write_vcf(vcffile, rout, Qsig, reference):
    with Path(vcffile).open("w") as g:
        g.write(
            "##fileformat=VCFv4.2\n##reference="
            + reference
            + '\n##source=umierrorcorrectV0.1\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
            + '##INFO=<ID=AD,Number=1,Type=Float,Description="Alternative Allele Depth">\n'
            + '##INFO=<ID=AF,Number=A,Type=Integer,Description="Alternative Allele Frequency">\n'
            + '##FILTER=<ID=a5,Description="Alternative Allele Depth below 5">\n'
            + '##FILTER=<ID=q10,Description="Variant quality below 10">\n'
            + '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Allele Depth">\n'
        )
        g.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]) + "\n")

        for r, q in zip(rout, Qsig):
            if q == "NA":
                q = "."
            parts = r.split("\t")
            vcffilter = []
            if q != "." and float(q) < 10:
                vcffilter.append("q10")
            if int(parts[-3]) < 5:
                vcffilter.append("a5")
            vcffilter = "PASS" if len(vcffilter) == 0 else ",".join(vcffilter)
            newline = "\t".join(
                [
                    parts[1],
                    parts[2],
                    ".",
                    parts[4],
                    parts[-1],
                    str(q),
                    vcffilter,
                    f"DP={parts[-5]};AD={parts[-3]};AF={parts[-2]}",
                    "DP",
                    parts[-3],
                ]
            )
            g.write(newline + "\n")


def plot_histogram(hist, plot_filename):
    num_bins = 100
    n, bins, patches = plt.hist(hist, num_bins, facecolor="dodgerblue", alpha=0.5)
    plt.xlabel("Q-score")
    plt.ylabel("Frequency")
    plt.title("Histogram of Q-scores")
    plt.box(False)
    plt.xlim(0, 140)
    plt.savefig(plot_filename)


def run_call_variants(args):
    if not args.cons_file:
        args.cons_file = str(list(Path(args.output_path).glob("*cons.tsv"))[0])
    if not args.sample_name:
        args.sample_name = get_sample_name_from_cons(args.cons_file)
    args.fsize = int(args.fsize)
    f1, n1, a1, data = parse_cons_file(args.cons_file, args.fsize)
    f1 = np.array(f1)
    n1 = np.array(n1)
    a1 = np.array(a1)
    data = np.array(data)
    if args.vc_method.lower() == "count":
        rout = data[a1 >= float(args.count_cutoff)]
        Qsig = ["NA"] * len(rout)
    # result=get_beta_parameters(f1[np.isin(pos,spikepositions)!=True])
    params = []
    if args.params_file:
        with Path(args.params_file).open() as f:
            for line in f:
                line = line.rstrip()
                params.append(float(line))
    else:
        params = [2.168215069116764, 3531.588541594945]
    a = 1 - betabinom.cdf(a1 - 1, n1, params[0], params[1])
    # print(params[0])
    # print(params[1])
    a[a == inf] = 1e-10
    a[np.isnan(a)] = 1e-10
    a[a == 0] = 1e-10
    Q = -10 * np.log10(a)
    data = np.array(data)
    # plot_histogram(Q,args.output_path+'/'+args.sample_name+'.histogram.png')
    if args.vc_method.lower() == "bbmodel":
        rout = data[float(args.qvalue_threshold) <= Q]
        Qsig = Q[float(args.qvalue_threshold) <= Q]
    else:
        rout = data[a1 >= float(args.count_cutoff)]
        Qsig = Q[a1 >= float(args.count_cutoff)]
    outfilename = Path(args.output_path) / f"{args.sample_name}.vcf"
    write_vcf(str(outfilename), rout, Qsig, args.reference_file)
