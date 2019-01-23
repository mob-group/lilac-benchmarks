#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

from data import *
from scipy.stats.mstats import gmean

def speedup(natives, others):
    ret = dict()
    names = list_values(natives).names
    for n in names:
        nat_mean = query(natives, name=n)
        o_mean = query(others, name=n)
        if len(nat_mean) > 0 and len(o_mean) > 0:
            ret[n] = nat_mean[0].mean() / o_mean[0].mean()
    return ret

# Question: how much faster is each implementation at bfs on Monaco?
# Answer: geomean across all benches for each implementation
def monaco_bfs(dataset):
    on_monaco = query(dataset, platform='monaco', benchmark='bfs')
    native_result = query(on_monaco, implementation='native')

    vals = list_values(on_monaco)
    for impl in vals.implementations:
        if impl != 'native':
            results = query(on_monaco, implementation=impl)
            speeds = speedup(native_result, results)
            print(impl)
            print(gmean([speeds[k] for k in speeds]))

plot_list = [
    monaco_bfs
]

plots = { 
    f.__name__ : f for f in plot_list
}

if __name__ == '__main__':
    sns.set()

    parser = argparse.ArgumentParser(description='Make plots for ICS')
    parser.add_argument('data', type=str, help='Data file to load', nargs='+')
    parser.add_argument('--dump', action='store_true')
    parser.add_argument('--plot', choices=plots)
    args = parser.parse_args()

    dataset = read_rows(args.data)

    if args.dump:
        for r in dataset:
            print(r)

    if args.plot:
        plots[args.plot](dataset)
