#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def all_matrices(frame):
    return set(frame.matrix)

def speedup(sub_f):
    base_time = sub_f.query('method == "native"').time
    speedups = {}
    for _, row in sub_f.iterrows():
        speedups[row.method] = float(base_time / row.time)
        speedups['nnz'] = row.nnz
    return speedups

if __name__ == "__main__":
    frame = pd.read_csv(sys.argv[1], sep=' ')
    results = {}
    for name in all_matrices(frame):
        sub_f = frame.query('matrix == @name')
        results[name] = speedup(sub_f)
    data = pd.DataFrame(results).T
