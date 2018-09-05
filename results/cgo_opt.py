#!/usr/bin/env python

import csv
import matplotlib.pyplot as plt
import sys

def plot(data):
    fig, ax = plt.subplots(figsize=(3,1.5))
    bar_style = {
        'width': 0.8,
        'edgecolor': 'black',
        'align': 'edge'
    }

    ax.set_axisbelow(True)
    ax.yaxis.grid(True, ls=':')
    ax.tick_params(axis=u'both', which=u'both',length=0)

    for i, bench in enumerate(data):
        ax.bar(i, float(bench['speedup']), color='grey', **bar_style)

    ticks, labels = zip(*[(i + 0.75, d['name']) for i, d in enumerate(data)])
    ax.set_xticks(ticks)
    ls = ax.set_xticklabels(labels, rotation=-60, fontsize=6)
    # ax.invert_yaxis()
    ax.axhline(1, color='black', lw=1)
    ax.set_ylabel('Speedup')
    ax.set_yticks([0, 2, 4, 6])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    fig.tight_layout()

if __name__ == "__main__":
    with open(sys.argv[1]) as csvfile:
        reader = csv.DictReader(csvfile)
        data = list(reader)
    plot(data)
    plt.savefig('opt.pdf')
    # plt.show()
