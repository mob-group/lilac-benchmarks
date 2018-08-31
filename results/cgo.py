#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pylab
import sys

plt.rcParams['hatch.color'] = '#FFFFFF'
plt.rcParams['hatch.linewidth'] = 1

def hatches(n):
    return ['//', '\\\\', '///'][n]

def colors(n):
    return pylab.get_cmap('Greys')(np.linspace(0, 0.4, n))

def make_bars(ax, labels, times):
    cs = colors(len(labels))
    bars = []
    for i, (lab, t) in enumerate(zip(labels, times)):
        bars.append(ax.bar(i, t, color=cs[i], edgecolor='black',
            hatch=hatches(i)))
    ax.set_xticks([])
    # ax.set_xticklabels(labels)
    # ax.axhline(min(times), color='black', linestyle=':')
    ax.set_title('Woo')
    # ax.legend(bars, ('b1', 'b2', 'b3'))

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print('no')
        sys.exit(1)
    labels = sys.argv[1:4]
    times = [float(t) for t in sys.argv[4:]]
    fig, axes = plt.subplots(3, 6, squeeze=False)
    for bench, row in zip(['NAS-CG', 'Netlib', 'Netlib-C'], axes):
        for i, ax in enumerate(row):
            # if i == 0:
            #     ax.set_axis_off()
            #     ax.text(0.7, 0.5, bench, verticalalignment='center',
            #             transform=ax.transAxes)
            # else:
            make_bars(ax, labels, times)
    fig.tight_layout()
    fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.15)
    axes[2][2].legend(['Native', 'MKL', 'GPU'], loc='upper center',
            bbox_to_anchor=(1, -0.2), ncol=3)
    plt.show()
