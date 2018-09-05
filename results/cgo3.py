#!/usr/bin/env python

import csv
import matplotlib.pyplot as plt
import sys

# Merge data by benchmark, then for each pair plot a pair of bars with opposing
# styles
# add label below

def merged_data(data):
    ret = {}
    for row in data:
        if row['name'] not in ret:
            ret[row['name']] = {}
        ret[row['name']][row['platform']] = float(row['speedup'])
    return ret

def plot(data):
    fig, ax = plt.subplots(figsize=(7,3))
    bar_style = {
        'width': 0.4,
        'edgecolor': 'black',
        'align': 'edge'
    }

    ax.set_axisbelow(True)
    ax.yaxis.grid(True)

    for i, bench in enumerate(data):
        ax.bar(i + 0.1, data[bench].get('Intel', 0), color='grey', **bar_style)
        ax.bar(i + 0.5, data[bench].get('AMD', 0), color='lightgrey', **bar_style)

    ticks, labels = zip(*[(i + 0.5, bench) for i, bench in enumerate(data)])
    ax.set_xticks(ticks)
    ls = ax.set_xticklabels(labels, rotation=-60, fontsize=6)
    # ax.invert_yaxis()
    ax.legend(('Intel', 'AMD'))
    ax.axhline(1, color='black', ls=':')
    ax.set_ylabel('Speedup')
    ax.set_xlabel('Benchmark')
    fig.tight_layout()

if __name__ == "__main__":
    with open(sys.argv[1]) as csvfile:
        reader = csv.DictReader(csvfile)
        data = merged_data(reader)
    plot(data)
    plt.savefig('bar.pdf')
