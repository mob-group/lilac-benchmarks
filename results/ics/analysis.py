#!/usr/bin/env python

import argparse
import pandas as pd
import scipy.stats.mstats as mstats

def get_speedups(df, expert=False):
    native_times = (
        df.query('implementation == "native"')
          .groupby(['benchmark','platform', 'name'])
          .time.mean())

    if not expert:
        df = df[~df.implementation.str.contains('expert')]

    all_times = (
        df.groupby(['benchmark', 'platform', 'name', 'implementation'])
          .time.mean())

    for k, v in all_times.iteritems():
        all_times[k] = native_times[k[:-1]] / v

    return (
        all_times.to_frame()
                 .reset_index()
                 .rename(index=str, columns={'time': 'speedup'})
    )

def best_per_benchmark(speedups):
    mean_speedups = (
        speedups.groupby(['benchmark', 'platform', 'implementation'])['speedup']
                .apply(mstats.gmean)
                .reset_index())

    max_mean_speedups = (
        mean_speedups.groupby(['benchmark', 'platform'])['speedup']
                     .idxmax())

    return mean_speedups.loc[max_mean_speedups]

# From a group of mean speedups by benchmark, platform, implementation, get the
# speedup for a particular implementation.
def speedup_for_impl(impl, group):
    for r in group.iterrows():
        if r[1]['implementation'] == impl:
            return r[1]['speedup']
    return float('nan')

# Get the performance ratio
def speedup_pairs(g, targets):
    return pd.Series(
        [speedup_for_impl(t[0], g) / speedup_for_impl(t[1], g) for t in targets],
        index=[t[1] for t in targets]
    )

def compare_speedups(speedups, targets, sparse=False):
    mean_speedups = (
        speedups.groupby(['benchmark', 'platform', 'implementation'])['speedup']
                .apply(mstats.gmean)
                .reset_index())

    return (
        mean_speedups.groupby(['benchmark', 'platform'])
                     .apply(lambda g: speedup_pairs(g, targets))
                     .dropna(how='all' if sparse else 'any')
                     .reset_index()
    )

def best_vs_expert(speedups):
    targets = [
        ('gpu', 'opencl-expert'),
        ('mkl', 'openmp-expert')
    ]

    return compare_speedups(speedups, targets)

def vs_marshalled(speedups):
    targets = [
        ('mkl', 'mkl-slow'),
        ('opencl10', 'opencl10-slow'),
        ('opencl00', 'opencl00-slow'),
        ('gpu', 'gpu-slow'),
        ('sparsex', 'sparsex-slow'),
    ]

    return compare_speedups(speedups, targets, sparse=True)

modes = {
    'baseline' : best_per_benchmark,
    'expert' : best_vs_expert,
    'marshall' : vs_marshalled
}

def print_results(df, out=None):
    if out is None:
        print(df.to_csv(index=False))
    else:
        df.to_csv(out, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Data analysis script for ICS results.')
    parser.add_argument('mode', choices=modes, help='Mode to use')
    parser.add_argument('data', type=str, help='Data file to analyse')
    parser.add_argument('-o', '--output', type=str, help='Output file to write')
    args = parser.parse_args()

    df = pd.read_csv(args.data)
    speeds = get_speedups(df, args.mode=='expert')
    bests = modes[args.mode](speeds)

    print_results(bests, args.output)
