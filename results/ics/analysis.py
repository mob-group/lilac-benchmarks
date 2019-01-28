#!/usr/bin/env python

import argparse
import pandas as pd

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Data tidying script for ICS results.')
    parser.add_argument('data', type=str, help='Data file to tidy')
    args = parser.parse_args()

    df = pd.read_csv(args.data)
    print(get_speedups(df))
