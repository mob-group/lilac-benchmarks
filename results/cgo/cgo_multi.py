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
MID_LILAC = rgb(155,120,156)
DARKER_LILAC = rgb(91,60,93)

def flipped(data):
    ret = {
        'intel': [],
        'amd': [],
    }

    for row in data:
        ret['intel'].append((float(row['mkl']), float(row['gpu'])))
        ret['amd'].append((float(row['external']), float(row['integrated'])))
    for k in ret:
        ret[k] = sorted(ret[k])
    return ret

def plot(data):
    fig, ax = plt.subplots(figsize=(3,2))

    scatter_style = {
        'linewidths': 0.5,
    }

    k = 'intel'
    plt.scatter(range(len(data[k])), [p[0] for p in data[k]], color=DARK_LILAC,
            marker='o', **scatter_style)
    plt.scatter(range(len(data[k])), [p[1] for p in data[k]], color=LILAC,
            marker='^', **scatter_style)

    ax.legend(('MKL', 'GPU'))
    ax.axhline(1, color='black', lw=1)
    ax.set_xticks([])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax.set_axisbelow(True)
    ax.yaxis.grid(True, ls=':')
    ax.tick_params(axis=u'both', which=u'both',length=0)

    ax.set_title('Intel')
    ax.set_ylabel('Speedup')

    fig.tight_layout()

    # bar_style = {
    #     'width': 0.15,
    #     'edgecolor': 'black',
    #     'align': 'edge'
    # }


    # for i, bench in enumerate(data):
    #     ax.bar(i + 0, float(bench['mkl']), color=LILAC, **bar_style)
    #     ax.bar(i + 0.25, float(bench['gpu']), color=DARK_LILAC, **bar_style)
    #     ax.bar(i + 0.5, float(bench['integrated']), color='red', **bar_style)
    #     ax.bar(i + 0.75, float(bench['external']), color='blue', **bar_style)

    # ticks, labels = zip(*[(i + 0.75, bench) for i, bench in enumerate(data)])
    # ax.set_xticks(ticks)
    # ls = ax.set_xticklabels(labels, rotation=-60, fontsize=6)
    # # ax.invert_yaxis()
    # ax.legend(('Intel', 'AMD'))
    # ax.axhline(1, color='black', lw=1)
    # ax.set_yticks([0, 2.5, 5, 7.5, 10])

if __name__ == "__main__":
    with open(sys.argv[1]) as csvfile:
        reader = csv.DictReader(csvfile)
        data = flipped(reader)
    plot(data)
    if len(sys.argv) > 2:
        plt.savefig(sys.argv[2])
    else:
        plt.show()
