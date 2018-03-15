#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import sys

def plot_data(filename, ax):
    frame = pd.read_csv(filename)
    ax.errorbar(frame['threads'], frame['speedup'], yerr=frame['stdev'],
            capsize=2, linewidth=1, ls=':')

if __name__ == "__main__":
    files = sys.argv[1:]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for f in files:
        plot_data(f, ax)
    ax.set_xlabel("No. of MKL Threads")
    ax.set_ylabel("Speedup")
    plt.show()
