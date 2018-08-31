#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab

plt.rcParams['hatch.color'] = '#FFFFFF'
plt.rcParams['hatch.linewidth'] = 1

NAMES = {
    'nas': 'NAS-CG',
    'netlib': 'Netlib',
    'netlib-c': 'Netlib-C'
}

HATCHES = [
    '////', '....',
    '////', '....'
]

def colors(n):
    return pylab.get_cmap('Greys')(np.linspace(0, 0.5, n))

def load_data(fn):
    frame = pd.read_csv(fn)
    frame['mklspeed'] = frame.intelcpu / frame.mkl
    frame['gpuspeed'] = frame.intelcpu / frame.gpu
    frame['clspeed'] = frame.amdcpu / frame.opencl
    frame['clgpuspeed'] = frame.amdcpu / frame.clgpu
    return frame

def title_for(data):
    return "{} {}".format(NAMES[data.benchmark], data['size'].upper())

def make_bars(data, ax):
    cs = colors(4)
    vals = (data.mklspeed, data.gpuspeed, data.clspeed, data.clgpuspeed)
    bars = []
    for i, val in enumerate(vals):
        bars.append(ax.bar(i, val, color=cs[3 - (i < 2)*3], edgecolor='black',
            hatch=HATCHES[i]))
    ax.axhline(1, linestyle=':', color='black')
    ax.set_title(title_for(data))
    ax.set_xticks([])
    return bars

if __name__ == "__main__":
    frame = load_data('cgo_data.csv')
    fig, axes = plt.subplots(3, 6)
    fig.set_figheight(4)
    fig.set_figwidth(7)
    for group, row in zip(frame.groupby('benchmark'), axes):
        for (_, data), ax in zip(group[1].iterrows(), row):
            bars = make_bars(data, ax)
    fig.tight_layout()
    fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.15)
    axes[2][2].legend(bars, ['CUDA', 'MKL', 'OpenCL CPU', 'OpenCL GPU'], loc='upper center',
            bbox_to_anchor=(1.5, -0.2), ncol=4)
    # plt.savefig('draft.pdf', dpi=120)
    plt.show()
