#!/usr/bin/env python

import argparse
import pandas as pd

def tidy(filename):
    df = pd.read_csv(filename, header=None)
    df = pd.melt(df, id_vars=range(0,4), value_vars=range(4,9))
    df = df.drop(axis=1, labels='variable')
    df = df.rename(index=str, columns={
        0: 'platform',
        1: 'benchmark',
        2: 'implementation',
        3: 'name',
        'value': 'time'
    })
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Data tidying script for ICS results.')
    parser.add_argument('data', type=str, help='Data file to tidy')
    parser.add_argument('--output', type=str, help='Output file to write')
    args = parser.parse_args()

    tidied = tidy(args.data)

    if args.output is None:
        print(tidied.to_csv(), index=False)
    else:
        tidied.to_csv(args.output, index=False)
