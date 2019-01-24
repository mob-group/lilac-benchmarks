#!/usr/bin/env python
import argparse

from collections import defaultdict
from data import *
from scipy.stats.mstats import gmean

def geom_speedup(dataset, plat, impl, benchmark):
    native = query(dataset, platform=plat, implementation='native', benchmark=benchmark)
    comp = query(dataset, platform=plat, implementation=impl, benchmark=benchmark)
    sizes = list_values(native).names.intersection(list_values(comp).names)
    speeds = []
    for s in sizes:
        native_mean = query(native, name=s)[0].mean()
        comp_mean = query(comp, name=s)[0].mean()
        speeds.append(native_mean / comp_mean)
    if speeds:
        return gmean(speeds)

def best_impls(dataset):
    bests = defaultdict(lambda: defaultdict(dict))

    impls = list_values(dataset).implementations
    benches = list_values(dataset).benchmarks
    platforms = list_values(dataset).platforms
    for p in platforms:
        for b in benches:
            best_sp = 0.0
            best_impl = ''
            for i in impls:
                sp = geom_speedup(dataset, p, i, b)
                if sp is not None:
                    if sp > best_sp:
                        best_sp = sp
                        best_impl = i
            bests[p][b] = (best_impl, best_sp)
    return bests

def vs_baseline(dataset):
    bests = best_impls(dataset)
    platforms = list_values(dataset).platforms

    for p in platforms:
        bs = bests[p]
        for b in bs:
            print("{},{},{},{}".format(p, b, bs[b][0], bs[b][1]))

transform_list = [
    vs_baseline
]

transforms = { 
    f.__name__ : f for f in transform_list
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make derived data for ICS')
    parser.add_argument('data', type=str, help='Data file to load', nargs='+')
    parser.add_argument('--dump', action='store_true')
    parser.add_argument('--transform', choices=transforms)
    args = parser.parse_args()

    dataset = read_rows(args.data)

    if args.dump:
        for r in dataset:
            print(r)

    if args.transform:
        transforms[args.transform](dataset)
