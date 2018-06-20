#!/usr/bin/env python3

import os
import matplotlib.pyplot as plt
import numpy as np
import sys

def all_times(path):
    with open(path) as f:
        all_times = {}
        for line in f:
            fields = line.split(' ')
            name = fields[0]
            times = {
                'base' : float(fields[1]),
                'integrated': float(fields[2]),
                'gpu': float(fields[3]),
                'nnz': int(fields[4])
            }
            all_times[name] = times
        return all_times

def speedups(all_times):
    speedups = {}
    for matrix in all_times:
        base_time = all_times[matrix]['base']
        speedups[matrix] = {
            'integrated' : base_time / all_times[matrix]['integrated'],
            'gpu' : base_time / all_times[matrix]['gpu'],
            'nnz' : all_times[matrix]['nnz']
        }
    return speedups

if __name__ == "__main__":
    ats = all_times(sys.argv[1])
    speeds = speedups(ats)
    xs = np.array([speeds[name]['nnz'] for name in speeds])
    ys = np.array([speeds[name]['gpu'] for name in speeds])
    fig, ax = plt.subplots(1, 1)
    ax.scatter(xs, ys)
    # fit = np.polyfit(xs, ys, 1)
    # fit_ys = fit[0] * xs + fit[1]
    # ax.plot(xs, fit_ys)
    ax.axhline(1, color='black')
    plt.show()
