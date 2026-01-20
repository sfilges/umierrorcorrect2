#!/usr/bin/env python3
from pathlib import Path

import numpy as np
from scipy.optimize import fmin
from scipy.stats import beta

from umierrorcorrect.core.utils import parse_cons_file


def betaNLL(params, *args):
    a, b = params
    data = np.array(args[0])
    pdf = beta.pdf(data, a, b, loc=0, scale=1)
    lg = np.log(pdf)
    # lg=np.where(lg==-np.inf,0,lg)
    mask = np.isfinite(lg)
    nll = -lg[mask].sum()
    nll = -1 * np.sum(lg)
    return nll


def get_beta_parameters(data):
    m = np.mean(data)
    v = np.var(data)
    a0 = m * (m * (1 - m) / v - 1)
    b0 = (1 - m) * (m * (1 - m) / v - 1)
    result = fmin(betaNLL, [a0, b0], args=(data,))
    return result


def run_fit_bgmodel(args):
    spikepositions = [178952085, 55599321, 7577558, 7577547, 7577538, 7577120]
    if args.nonbgposfile:
        nonbgpos = []
        with Path(args.nonbgposfile).open() as f:
            for line in f:
                line = line.rstrip()
                nonbgpos.append(line)
    else:
        nonbgpos = spikepositions
    if not args.cons_file:
        args.cons_file = str(list(Path(args.output_path).glob("*cons.tsv"))[0])
    args.fsize = int(args.fsize)
    f1, n1, a1, pos, data = parse_cons_file(args.cons_file, args.fsize, include_position=True)
    f1 = np.array(f1)
    n1 = np.array(n1)
    a1 = np.array(a1)
    pos = np.array(pos)
    data = np.array(data)
    result = get_beta_parameters(f1[np.isin(pos, nonbgpos) is not True])
    # a=prob_bb(n1,a1,result[0],result[1])
    print(pos, nonbgpos, np.isin(pos, nonbgpos))
    with Path(args.out_file).open("w") as g:
        g.write(f"{result[0]}\n")
        g.write(f"{result[1]}\n")
    # a[a==inf]=1e-10
    # a[np.isnan(a)]=1e-10
    # Q = -10*np.log10(a)
    # data=np.array(data)
    # plot_histogram(Q,args.output_path+'/'+args.sample_name+'.histogram.png')
    # if args.vc_method.lower()=='bbmodel':
    #    rout=data[Q >= float(args.qvalue_threshold)]
    #    Qsig=Q[Q >= float(args.qvalue_threshold)]
    # else:
    #    rout=data[a1 >= float(args.count_cutoff)]
    #    Qsig=Q[a1 >= float(args.count_cutoff)]
    # outfilename=args.output_path+'/'+args.sample_name+'2.vcf'
    # write_vcf(outfilename,rout,Qsig,args.reference_file)
