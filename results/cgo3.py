#!/usr/bin/env python

import csv
import matplotlib.pyplot as plt
import sys

# Merge data by benchmark, then for each pair plot a pair of bars with opposing
# styles
# add label below

def rgb(*args):
    return tuple(v/255.0 for v in args)

LILAC = rgb(200,162,200)
DARK_LILAC = rgb(116,83,117)

def merged_data(data):
    ret = {}
    for row in data:
        if row['name'] not in ret:
            ret[row['name']] = {}
        ret[row['name']][row['platform']] = float(row['speedup'])
    return ret

def plot(data):
    fig, ax = plt.subplots(figsize=(3,2))
    bar_style = {
        'width': 0.4,
        'edgecolor': 'black',
        'align': 'edge'
    }

    ax.set_axisbelow(True)
    ax.yaxis.grid(True, ls=':')
    ax.tick_params(axis=u'both', which=u'both',length=0)

    for i, bench in enumerate(data):
        ax.bar(i + 0.1, data[bench].get('Intel', 0), color=LILAC, **bar_style)
        ax.bar(i + 0.5, data[bench].get('AMD', 0), color=DARK_LILAC, **bar_style)

    ticks, labels = zip(*[(i + 0.75, bench) for i, bench in enumerate(data)])
    ax.set_xticks(ticks)
    ls = ax.set_xticklabels(labels, rotation=-60, fontsize=6)
    # ax.invert_yaxis()
    ax.legend(('Intel', 'AMD'))
    ax.axhline(1, color='black', lw=1)
    ax.set_ylabel('Speedup')
    ax.set_yticks([0, 1, 2, 3])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    fig.tight_layout()

if __name__ == "__main__":
    with open(sys.argv[1]) as csvfile:
        reader = csv.DictReader(csvfile)
        data = merged_data(reader)
    plot(data)
    if len(sys.argv) > 2:
        plt.savefig(sys.argv[2])
    else:
        plt.show()
