#!/usr/bin/env python

import argparse
import pandas as pd
import scipy.stats.mstats as mstats

def get_speedups(df):
    native_times = (
        df.query('implementation == "native"')
          .groupby(['benchmark','platform', 'name'])
          .time.mean())

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Data analysis script for ICS results.')
    parser.add_argument('data', type=str, help='Data file to analyse')
    parser.add_argument('-o', '--output', type=str, help='Output file to write')
    args = parser.parse_args()

    df = pd.read_csv(args.data)
    speeds = get_speedups(df)
    bests = best_per_benchmark(speeds)

    if args.output is None:
        print(bests.to_csv(index=False))
    else:
        bests.to_csv(args.output, index=False)
