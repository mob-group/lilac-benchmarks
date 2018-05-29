#!/usr/bin/env python

import matplotlib.pyplot as plt
import re
import statistics

from enum import Enum, auto
from matplotlib.colors import hsv_to_rgb
from typing import Dict, List

plt.style.use('seaborn')

class ResultSet:
    def __init__(self, name, variants):
        self.name = name
        self.__variants = variants

    def __repr__(self):
        return "ResultSet({}, {})".format(self.name, self.__variants)

    def means(self):
        ret = {}
        for var_name in self.__variants:
            ret[var_name] = statistics.mean(self.__variants[var_name])
        return ret

    def speedup(self, base, improved):
        ms = self.means()
        return ms[base] / ms[improved]

    def chart_data(self):
        return self.name, {
            'CUDA' : self.speedup('Native', 'CUDA'),
            'MKL' : self.speedup('Native', 'MKL'),
            'OpenCL-iGPU' : self.speedup('AMD-Native', 'OpenCL-iGPU'),
            'OpenCL-GPU' : self.speedup('AMD-Native', 'OpenCL-GPU'),
        }

def read_results(results_file):
    results = []
    current = {}

    for line in results_file:
        if not line.startswith(' '): 
            if current != {}:
                results.append(ResultSet(name, current))
                current = {}
            name = line.split(':')[0]
        else:
            line = re.sub(' +', ' ', line).strip()
            bench = line.split(':')[0]
            data = [float(d) for d in line.split(':')[1].split(',')]
            current[bench] = data
    
    results.append(ResultSet(name, current))
    return results

def autolabel(ax, rects):
    for rect in rects:
        height = rect.get_height()
        if height < 1:
            text_height = height + 0.25
        else:
            text_height = 1.007 * height
        ax.text(rect.get_x() + rect.get_width()/2., 1.0*text_height,
                '{:.2f}×'.format(height),
                ha='center', va='bottom', fontsize=9)

def plot(data):
    bar_width = 0.15
    fig, axes = plt.subplots(1, len(data), figsize=(10,2))
    for (ax, d) in zip(axes, data):
        # ax.set_axisbelow(True)
        # ax.yaxis.grid(color='gray', linestyle='dotted')
        # ax.xaxis.grid(color='gray', linestyle='dotted')

        rs1 = ax.bar(-3*bar_width/2, d[1]['CUDA'], bar_width*0.9,
                edgecolor='black', color=hsv_to_rgb((0.6, 1, 0.8)))
        rs2 = ax.bar(-bar_width/2, d[1]['MKL'], bar_width*0.9,
                edgecolor='black', color=hsv_to_rgb((0.6, 0.4, 0.8)))
        rs3 = ax.bar(bar_width/2, d[1]['OpenCL-iGPU'], bar_width*0.9,
                edgecolor='black', color=hsv_to_rgb((0., 1, 0.8)))
        rs4 = ax.bar(3*bar_width/2, d[1]['OpenCL-GPU'], bar_width*0.9,
                edgecolor='black', color=hsv_to_rgb((0., 0.4, 0.8)))
        autolabel(ax, rs1)
        autolabel(ax, rs2)
        autolabel(ax, rs3)
        autolabel(ax, rs4)

        ax.axhline(1, color='black')
        ax.tick_params(axis='x', rotation=-45)
        ax.set_xticks([-3*bar_width/2, -bar_width/2, bar_width/2, 3*bar_width/2])
        ax.set_xticklabels(['CUDA', 'MKL', 'iGPU', 'GPU'])
        ax.set_xlabel(d[0])
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    with open('all.txt', 'r') as r_f:
        results = read_results(r_f)
        data = [d.chart_data() for d in results]
        fig = plot(data)
        fig.savefig('perf.pdf')
