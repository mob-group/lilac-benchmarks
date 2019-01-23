#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import sys

def label(fn):
    return {
        "npb" : "NAS-CG",
        "pfold" : "Wales-PFold",
        "ngt" : "Wales-NGT",
        "netlib" : "Netlib-CRS",
        "cnetlib" : "Netlib-C-CRS",
    }[fn.split("_")[0]]

def plot_data(filename, ax):
    s = slice(1,16,1)
    frame = pd.read_csv(filename)
    # ax.errorbar(frame['threads'][s], frame['speedup'][s],
    #         yerr=frame['stdev'][s], capsize=2, linewidth=1, ls=':')
    ax.plot(frame['threads'][s], frame['speedup'][s], linewidth=1, ls='-')
    ax.set_xticks(frame['threads'][s])
    return "{:.2f}".format(frame['speedup'][7])

if __name__ == "__main__":
    files = sys.argv[1:]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    legend = []
    for f in files:
        speedup = plot_data(f, ax)
        legend.append(label(f) + " (" + speedup + "×)")
    ax.legend(legend, fontsize='small')
    ax.axhline(y=1,lw=1,ls=':')
    ax.axvline(x=8,lw=1,ls=':')
    ax.set_xlabel("No. of MKL Threads")
    ax.set_ylabel("Speedup")
    plt.savefig('speedup.pdf', papertype='a4')
