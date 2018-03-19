#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import sys

def label(fn):
    return {
        "npb" : "NAS-CG",
        "pfold" : "Wales-PFold",
        "sparsebench" : "Netlib-CRS",
        "ngt" : "Wales-NGT"
    }[fn.split("_")[0]]

def plot_data(filename, ax):
    s = slice(0,None,1)
    frame = pd.read_csv(filename)
    ax.errorbar(frame['threads'][s], frame['speedup'][s],
            yerr=frame['stdev'][s], capsize=2, linewidth=1, ls=':')
    ax.set_xticks(frame['threads'][s])

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
