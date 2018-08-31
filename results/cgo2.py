#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd

NAMES = {
    'nas': 'NAS-CG',
    'netlib': 'Netlib',
    'netlib-c': 'Netlib-C'
}

def load_data(fn):
    frame = pd.read_csv(fn)
    frame['mklspeed'] = frame.intelcpu / frame.mkl
    frame['gpuspeed'] = frame.intelcpu / frame.gpu
    frame['clspeed'] = frame.amdcpu / frame.opencl
    frame['clgpuspeed'] = frame.amdcpu / frame.clgpu
    return frame

def title_for(data):
    return "{} {}".format(NAMES[data.benchmark], data['size'].upper())

if __name__ == "__main__":
    frame = load_data('cgo_data.csv')
    fig, axes = plt.subplots(3, 6)
    for group, row in zip(frame.groupby('benchmark'), axes):
        for (_, data), ax in zip(group[1].iterrows(), row):
            ax.bar(0, data.mklspeed)
            ax.bar(1, data.gpuspeed)
            ax.bar(2, data.clspeed)
            ax.bar(3, data.clgpuspeed)
            ax.axhline(1, linestyle=':', color='black')
            ax.set_title(title_for(data))
    fig.tight_layout()
    plt.savefig('draft.pdf', figsize=(5,7), dpi=120)
