#!/usr/bin/env python
import argparse

from data import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make plots for ICS')
    parser.add_argument('data', type=str, help='Data file to load', nargs='+')
    parser.add_argument('--dump', action='store_true')
    args = parser.parse_args()

    dataset = read_rows(args.data)

    if args.dump:
        for r in dataset:
            print(r)
