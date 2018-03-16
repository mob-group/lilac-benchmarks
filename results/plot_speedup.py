#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import sys

def label(fn):
    return {
        "npb" : "NAS-CG",
        "pfold" : "Wales-PFold",
    }[fn.split("_")[0]]

def plot_data(filename, ax):
    frame = pd.read_csv(filename)
    ax.errorbar(frame['threads'][1::2], frame['speedup'][1::2],
            yerr=frame['stdev'][1::2], capsize=2, linewidth=1, ls=':')
    ax.set_xticks(frame['threads'][1::2])

if __name__ == "__main__":
    files = sys.argv[1:]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    legend = []
    for f in files:
        plot_data(f, ax)
        legend.append(label(f))
    ax.set_xlabel("No. of MKL Threads")
    ax.set_ylabel("Speedup")
    ax.legend(legend)
    plt.show()
